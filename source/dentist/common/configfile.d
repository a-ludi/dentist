/**
    Validate and parse config files.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.configfile;

import darg :
    Argument,
    Help,
    isArgumentHandler,
    isOptionHandler,
    Option,
    OptionFlag;
import dentist.commandline : OptionsFor;
import dentist.common.commands :
    DentistCommand,
    dentistCommands;
import dentist.util.jsonschema :
    jsonSchema,
    schemaRef;
import dyaml :
    YAMLLoader = Loader,
    YAML = Node,
    YAMLType = NodeType;
import std.algorithm :
    all,
    canFind,
    endsWith,
    map,
    startsWith;
import std.conv : to;
import std.exception :
    basicExceptionCtors,
    enforce;
import std.format : format;
import std.math : isNaN;
import std.meta : Alias;
import std.range :
    ElementType,
    only;
import std.range.primitives;
import std.stdio : File, stderr;
import std.string :
    split,
    wrap;
import std.traits :
    arity,
    EnumMembers,
    getSymbolsByUDA,
    getUDAs,
    isArray,
    isAssignable,
    isCallable,
    isDynamicArray,
    isFloatingPoint,
    isIntegral,
    isSomeString,
    isStaticArray,
    isUnsigned,
    Parameters,
    ReturnType;
import std.typecons : No, tuple;
import vibe.data.json :
    Json,
    parseJson;


/// Identifier for the config object with default values.
enum configDefaultKey = "__default__";

/// Arguments with this value are assigned the config value.
enum configEmptyArgument = "-";

/// Keys prefixed with this string are ignored.
enum configCommentPrefix = "//";

/// Maximum size of a valid config file.
enum maxConfigSize = toBytes(256, SizeUnit.MiB);


/// Thrown if an error while handling config file occurs.
class ConfigFileException : Exception
{
    ///
    mixin basicExceptionCtors;
}


/// Decorator
struct ConfigType(alias T)
{
    alias Type = T;
}


/// Retroactively initialize `options` from `config`.
///
/// Note: Since this function has no knowledge about the creation process of
///     `options` it assumes that fields are not modified if they have their
///     default value. This means that a CLI option providing the default
///     value of that option does NOT overrule the config value. Here is a
///     small example:
///     ---
///     // config,yaml:
///     collectPileUps:
///       min-spanning-reads: 42
///
///     // invocation
///     dentist collect-pile-ups --config=config.yaml --min-spanning-reads=3 ...
///
///     // --> effectivate value of --min-spanning-reads is 42
///     ---
Options retroInitFromConfig(Options)(ref Options options, in string configFile)
{
    enum defaultOptions = Options.init;
    Options optionsFromConfig = parseConfig!Options(configFile);

    static foreach (member; __traits(allMembers, Options))
    {{
        alias symbol = Alias!(__traits(getMember, options, member));
        enum isMemberAssignable = __traits(compiles,
            __traits(getMember, options, member) = __traits(getMember, options, member)
        );

        static if (isMemberAssignable)
        {
            alias Member = typeof(__traits(getMember, options, member));
            enum unaryMixin(string template_) = format!template_(member);
            enum binaryMixin(string template_) = format!template_(member, member);
            alias assignConfigValue = () => mixin(binaryMixin!"options.%s = optionsFromConfig.%s");

            static if (getUDAs!(symbol, Argument).length > 0)
            {
                static if (isSomeString!Member)
                {
                    if (mixin(unaryMixin!"options.%s == configEmptyArgument"))
                        assignConfigValue();
                }
                else static if (isArray!Member && isSomeString!(ElementType!Member))
                {
                    if (mixin(unaryMixin!"options.%s.all!(v => v == configEmptyArgument)"))
                        assignConfigValue();
                }
            }
            else
            {
                static if (isStaticArray!Member || is(Member == class))
                {
                    if (mixin(binaryMixin!"options.%s == defaultOptions.%s"))
                        assignConfigValue();
                }
                else static if (isFloatingPoint!Member)
                {
                    if (
                        mixin(binaryMixin!"options.%s == defaultOptions.%s") ||
                        (
                            mixin(unaryMixin!"options.%s.isNaN") &&
                            mixin(unaryMixin!"defaultOptions.%s.isNaN")
                        )
                    )
                        assignConfigValue();
                }
                else
                {
                    if (mixin(binaryMixin!"options.%s is defaultOptions.%s"))
                        assignConfigValue();
                }
            }
        }
    }}

    return options;
}


/// Create `Options` from `configFile`.
Options parseConfig(Options)(in string configFile)
{
    enum dentistCommand = Options.commandName;
    auto configValues = parseConfigFile(configFile);

    validateConfig(configValues);

    // Extract options passed via `__default__` key
    auto defaultValues = configDefaultKey in configValues
        ? configValues[configDefaultKey]
        : Json.emptyObject;
    // Extract options passed via command name
    auto commandValues = dentistCommand in configValues
        ? configValues[dentistCommand]
        : Json.emptyObject;

    Options options;

    foreach (member; __traits(allMembers, Options))
    {
        alias symbol = Alias!(__traits(getMember, options, member));
        enum names = configNamesOf!symbol;

        static if (names.length > 0)
        {
            foreach (name; names)
                if (name in commandValues)
                    options.assignConfigValue!member(commandValues[name]);
                else if (name in defaultValues)
                    options.assignConfigValue!member(defaultValues[name]);
        }
    }

    return options;
}


/// Validates the format and structure of `configFile`. These rules must be
/// fulfilled:
/// $(UL
///     $(LI The document must contain a single object (dict type).)
///     $(LI The root object may have a `__default__` key with an object as
///         value.)
///     $(LI Each key in the `__default__` object must be a valid option
///         (secondary or short option names are allowed) or argument name
///         from ANY DENTIST command. The value type must be convertible
///         to the destination type.)
///     $(LI The root object may have a key for each DENTIST command
///         (hyphenated as used in the CLI) with an object as value.)
///     $(LI Each key in the command objects must be a valid option
///         (secondary or short option names are allowed) or argument name
///         from the named DENTIST command. The value type must be convertible
///         to the destination type.)
///     $(LI Keys starting with two slashes (`//`) are ignored on all levels.)
/// )
///
/// Throws: `ConfigFileException` if config is invalid.
void validateConfigFile(in string configFile)
{
    auto configValues = parseConfigFile(configFile);

    validateConfig(configValues);
}

/// ditto
protected void validateConfig(Json configValues)
{
    enforce!ConfigFileException(
        configValues.type == Json.Type.object,
        "config must contain a single object",
    );

    foreach (configKey, configValue; configValues.byKeyValue)
        enforce!ConfigFileException(
            (
                configKey.startsWith(configCommentPrefix) ^
                (configKey == configDefaultKey) ^
                only(dentistCommands).canFind(configKey)
            ),
            format!"invalid key `%s` in config"(configKey),
        );

    if (configDefaultKey in configValues)
        validateConfigDefault(configValues[configDefaultKey]);

    static foreach (command; EnumMembers!DentistCommand)
    {{
        alias Options = OptionsFor!command;

        if (Options.commandName in configValues)
            validateConfigCommand!Options(configValues[Options.commandName]);
    }}
}


/// Validate the `__default__` config values.
private void validateConfigDefault(Json defaultConfig)
{
    enforce!ConfigFileException(
        defaultConfig.type == Json.Type.object,
        "malformed default config: Got JSON of type " ~
        defaultConfig.type.to!string ~ ", expected object",
    );

    configLoop: foreach (configKey, configValue; defaultConfig.byKeyValue)
    {
        if (configKey.startsWith(configCommentPrefix))
            continue;

        static foreach (command; EnumMembers!DentistCommand)
        {{
            alias Options = OptionsFor!command;

            foreach (member; __traits(allMembers, Options))
            {
                alias symbol = Alias!(__traits(getMember, Options, member));
                enum names = configNamesOf!symbol;

                static if (names.length > 0)
                {
                    alias SymbolType = typeof(__traits(getMember, Options, member));

                    foreach (name; names)
                    {
                        try
                        {
                            if (name == configKey)
                            {
                                cast(void) getConfigValue!SymbolType(configValue);
                                continue configLoop;
                            }
                        }
                        catch (Exception cause)
                        {
                            throw new ConfigFileException(
                                format!"malformed config value `%s.%s`: %s"(
                                    configDefaultKey,
                                    configKey,
                                    cause.msg,
                                ),
                                cause,
                            );
                        }
                    }
                }
            }
        }}

        throw new ConfigFileException(
            format!"invalid config key `%s.%s`"(
                configDefaultKey,
                configKey,
            ),
        );
    }
}


/// Validate config for `Options.commandName`.
void validateConfigCommand(Options)(Json commandConfig)
{
    enum commandName = Options.commandName;

    configLoop: foreach (configKey, configValue; commandConfig.byKeyValue)
    {
        if (configKey.startsWith(configCommentPrefix))
            continue;

        foreach (member; __traits(allMembers, Options))
        {
            alias symbol = Alias!(__traits(getMember, Options, member));
            enum names = configNamesOf!symbol;

            static if (names.length > 0)
            {
                alias SymbolType = typeof(__traits(getMember, Options, member));

                foreach (name; names)
                {
                    try
                    {
                        if (name == configKey)
                        {
                            cast(void) getConfigValue!SymbolType(configValue);
                            continue configLoop;
                        }
                    }
                    catch (Exception cause)
                    {
                        throw new ConfigFileException(
                            format!"malformed config value `%s.%s`: %s"(
                                commandName,
                                configKey,
                                cause.msg,
                            ),
                            cause,
                        );
                    }
                }
            }
        }

        throw new ConfigFileException(
            format!"invalid config key `%s.%s`"(
                commandName,
                configKey,
            ),
        );
    }
}


/// Generate a JSON schema for DENTIST's config.
///
/// See_also: $(LINK https://json-schema.org/)
Json getJsonSchema()
{
    auto schema = Json([
        "$schema": Json("https://json-schema.org/draft/2020-12/schema"),
        "$id": Json("uri://a-ludi/dentist/v3/config.schema.json"),
        "title": Json("DENTIST config"),
        "description": Json("Configuration file content for DENTIST."),
        "$defs": Json([
            "option-list": Json([
                "anyOf": Json([
                    Json(["type": Json("string")]),
                    Json(["type": Json([Json("string")])]),
                ]),
            ]),
        ]),
        "type": Json("object"),
    ]);

    auto properties = Json.emptyObject();
    properties[configDefaultKey] = generateJsonSchemaDefault();

    static foreach (command; EnumMembers!DentistCommand)
        properties[OptionsFor!command.commandName] = generateJsonSchemaFor!command();

    schema["properties"] = properties;

    return schema;
}

unittest
{
    auto schema = getJsonSchema();

    assert(schema["properties"]["__default__"]["properties"]["reference"]["type"] == "string");
    assert(schema["properties"]["__default__"]["properties"]["quiet"]["type"] == "boolean");
    assert(schema["properties"]["output"]["properties"]["skip-gaps"]["type"] == schemaRef!"option-list");
}

protected Json generateJsonSchemaDefault()
{
    auto properties = Json.emptyObject();
    bool[string] conflicts;

    static foreach (command; EnumMembers!DentistCommand)
    {{
        alias Options = OptionsFor!command;

        foreach (member; __traits(allMembers, Options))
        {
            alias symbol = Alias!(__traits(getMember, Options, member));
            alias ConfigTypes = getUDAs!(symbol, ConfigType);

            enum propSchemata = propertySchema!(symbol, ConfigTypes);

            foreach (key, propSchema; propSchemata.byKeyValue())
            {
                if (key in properties && properties[key] != propSchema)
                {
                    conflicts[key] = true;
                    if (properties[key].type == Json.Type.array)
                        properties[key] ~= propSchema;
                    else
                        properties[key] = Json([properties[key], propSchema]);
                }
                else
                {
                    properties[key] = propSchema;
                }
            }
        }
    }}

    foreach (key; conflicts.keys)
    {
        auto schemata = properties[key];
        properties.remove(key);
        version (unittest)
            stderr.writefln!"removing multiple conflicting schemata for `%s`:\n%-(%s\n%)"(
                key,
                schemata,
            );
    }

    return Json([
        "type": Json("object"),
        "properties": properties,
    ]);
}


protected Json generateJsonSchemaFor(DentistCommand command, bool reportConflicts = true)()
{
    auto properties = Json.emptyObject();
    bool[string] conflicts;

    alias Options = OptionsFor!command;

    foreach (member; __traits(allMembers, Options))
    {
        alias symbol = Alias!(__traits(getMember, Options, member));
        alias ConfigTypes = getUDAs!(symbol, ConfigType);

        enum propSchemata = propertySchema!(symbol, ConfigTypes);

        foreach (key, propSchema; propSchemata.byKeyValue())
        {
            if (key in properties && properties[key] != propSchema)
            {
                conflicts[key] = true;
                if (properties[key].type == Json.Type.array)
                    properties[key] ~= propSchema;
                else
                    properties[key] = Json([properties[key], propSchema]);
            }
            else
            {
                properties[key] = propSchema;
            }
        }
    }

    foreach (key; conflicts.keys)
    {
        auto schemata = properties[key];
        properties.remove(key);
        if (reportConflicts)
            stderr.writefln!"removing multiple conflicting schemata for `%s`:\n%-(%s\n%)"(
                key,
                schemata,
            );
    }

    return Json([
        "type": Json("object"),
        "properties": properties,
    ]);
}


protected template propertySchema(alias symbol, CT = void)
{
    static if (is(typeof(symbol)) || !is(CT == void))
    {
        static if (is(CT == void))
            alias T = typeof(symbol);
        else
            alias T = CT.Type;
        enum names = configNamesOf!symbol;

        static if (names.length == 0 || is(CT == ConfigType!void))
        {
            enum propertySchema = Json.emptyObject;
        }
        else static if (isOptionHandler!T)
        {
            pragma(msg, "option " ~ names[0] ~ " has an OptionHandler");
            enum propertySchema = Json.emptyObject;
        }
        else
        {
            //enum arguments = getUDAs!(symbol, Argument);
            //enum options = getUDAs!(symbol, Option);
            enum helps = getUDAs!(symbol, Help);

            static if (isArgumentHandler!T)
                enum propSchema = Json(["type": schemaRef!"option-list"]);
            else
                enum propSchema = jsonSchema!T;

            static if (helps.length > 0)
                enum propertySchema = makeFullPropSchema(propSchema, names, helps[0].help);
            else
                enum propertySchema = makeFullPropSchema(propSchema, names);
        }
    }
    else
    {
        enum propertySchema = Json.emptyObject;
    }
}

private Json makeFullPropSchema(Json propSchema, string[] names, string help = null)
{
    if (help !is null)
        propSchema["description"] = help.wrap(size_t.max)[0 .. $ - 1];

    auto schema = Json.emptyObject;
    foreach (name; names)
        schema[name] = propSchema;

    return schema;
}


/// Get an array of names that can be used in a config file to reference
/// `symbol`. If symbol is an argument it will split the argument name by
/// colon (`:`) and return the last part stripped off its last character
/// (which should be a `>`). If symbol is an option it will return the list
/// of names that were provided. Otherwise, an empty array is returned.
protected template configNamesOf(alias symbol)
{
    alias optUDAs = getUDAs!(symbol, Option);
    alias argUDAs = getUDAs!(symbol, Argument);

    static if (argUDAs.length > 0)
        enum argName = argUDAs[0].name.split(":")[$ - 1][0 .. $ - 1];

    static if (optUDAs.length > 0 && argUDAs.length > 0)
    {
        enum configNamesOf = optUDAs[0].names ~ argName;
    }
    else static if (optUDAs.length > 0)
    {
        enum configNamesOf = optUDAs[0].names;
    }
    else static if (argUDAs.length > 0)
    {
        enum configNamesOf = [argName];
    }
    else
    {
        enum string[] configNamesOf = [];
    }
}


/// Convert and assign `configValue` to `member` in `options`.
///
/// Conversion rules:
/// $(UL
///     $(LI if `member` has a zero-argument function: $(UL
///         $(LI call `member` `configValue` times if `configValue` is an integer)
///         $(LI call `member` once if `configValue` is `true`)
///         $(LI do not call `member` if `configValue` is `false`)
///     ))
///     $(LI if `member` has a single-argument function: $(UL
///         $(LI call `member` with each value of `configValue` if
///             `configValue` is a string array)
///         $(LI call `member` with `configValue` if `configValue` is a string)
///     ))
///     $(LI if `member` can be assigned to `configValue` will be converted
///         to `member`'s type using `getConfigValue`)
///     $(LI if none of the above applies a `ConfigFileException` will be thrown)
/// )
protected void assignConfigValue(string member, Options)(ref Options options, Json configValue)
{
    alias SymbolType = typeof(__traits(getMember, options, member));

    static if (isOptionHandler!SymbolType)
    {
        if (configValue.type == Json.Type.int_)
            foreach (i; 0 .. configValue.get!ulong)
                __traits(getMember, options, member)();
        else if (configValue.type == Json.Type.bool_)
        {
            if (configValue.get!bool)
                __traits(getMember, options, member)();
        }
        else
            throw new ConfigFileException(
                "Got JSON of type " ~ configValue.type.to!string ~
                ", expected int_ or bool_."
            );
    }
    else static if (isArgumentHandler!SymbolType)
    {
        if (configValue.type == Json.Type.array)
            foreach (item; configValue.get!(Json[]))
                __traits(getMember, options, member)(item.get!string);
        else if (configValue.type == Json.Type.string)
            __traits(getMember, options, member)(configValue.get!string);
        else
            throw new ConfigFileException(
                "Got JSON of type " ~ configValue.type.to!string ~
                ", expected array or string_."
            );
    }
    else static if (isAssignable!SymbolType)
    {
        __traits(getMember, options, member) = getConfigValue!SymbolType(configValue);
    }
}


/// Convert `configValue` to `SymbolType`.
///
/// Conversion_Rules:
///
/// The rule (SymbolType ← JsonType) is picked based on `SymbolType` and the
/// provided JSON type(s) are permitted.
///
/// $(UL
///     $(LI `OptionFlag`|`bool` ← `bool`)
///     $(LI `enum` ← `string`)
///     $(LI `float`|`double`|`real` ← `int`|`float`)
///     $(LI unsigned|signed integer ← `int`)
///     $(LI `string` ← `string`)
///     $(LI `T[]`|T[n] ← `array`: apply above rules to each member;
///         for static arrays the lengths must match)
/// )
protected auto getConfigValue(SymbolType)(Json configValue)
{
    static if (is(SymbolType == OptionFlag))
        return configValue.get!bool.to!SymbolType;
    else static if (is(SymbolType == enum))
        return configValue.get!string.to!SymbolType;
    else static if (is(SymbolType == OptionFlag) || is(SymbolType : bool))
        return configValue.get!bool.to!SymbolType;
    else static if (isFloatingPoint!SymbolType)
    {
        if (configValue.type == Json.Type.int_)
            return configValue.get!long.to!SymbolType;
        else if (configValue.type == Json.Type.float_)
            return configValue.get!double.to!SymbolType;
        else
            throw new ConfigFileException(
                "Got JSON of type " ~ configValue.type.to!string ~
                ", expected float_ or int_."
            );
    }
    else static if (isIntegral!SymbolType && isUnsigned!SymbolType)
        return configValue.get!ulong.to!SymbolType;
    else static if (isIntegral!SymbolType && !isUnsigned!SymbolType)
        return configValue.get!long.to!SymbolType;
    else static if (isSomeString!SymbolType)
    {
        if (configValue.type == Json.Type.string)
            return configValue.get!string.to!SymbolType;
        else if (configValue.type == Json.Type.null_)
            return null;
        else
            throw new ConfigFileException(
                "Got JSON of type " ~ configValue.type.to!string ~
                ", expected string or null_."
            );
    }
    else static if (isArray!SymbolType)
    {
        SymbolType value;

        static if (isDynamicArray!SymbolType)
            value.length = configValue.length;
        else
            enforce!ConfigFileException(
                configValue.length == value.length,
                "array must have " ~ value.length ~ " elements",
            );

        foreach (size_t i, configElement; configValue.get!(Json[]))
            value[i] = getConfigValue!(ElementType!SymbolType)(configElement);

        return value;
    }
}


/// Parse `configFileName` into a `JSON` struct.
protected Json parseConfigFile(in string configFileName)
{
    auto configContent = readConfigFile(configFileName);

    if (configFileName.endsWith(".yaml", ".yml"))
    {
        auto loader = YAMLLoader.fromFile(configFileName);

        return loader.load().toJson;
    }
    else if (configFileName.endsWith(".json"))
        return parseJson(
            configContent,
            null,
            configFileName,
        );
    else
        throw new ConfigFileException("unknown file type");
}


/// Read `configFileName` into a string respecting a maximum file size of
/// `maxConfigSize` bytes.
protected string readConfigFile(in string configFileName)
{
    auto configFile = File(configFileName, "r");
    auto configFileSize = configFile.size;

    enforce!ConfigFileException(
        configFileSize <= maxConfigSize,
        format!"config file is too large; must be <= %.2f %s"(fromBytes(maxConfigSize).expand),
    );

    auto configContent = configFile.rawRead(new char[configFileSize]);

    return cast(string) configContent;
}


/// Convert `yaml` into `Json` as good as possible. Since YAML is a strict
/// superset of JSON it not always possible BUT every valid config file
/// can be converted.
Json toJson(YAML yaml)
{
    import std.datetime : SysTime;

    alias value(T) = () => yaml.get!(T, No.stringConversion)();

    final switch (yaml.type)
    {
        case YAMLType.null_:
            return Json(null);
        case YAMLType.boolean:
            return Json(value!bool());
        case YAMLType.integer:
            return Json(value!long());
        case YAMLType.string:
            return Json(value!string());
        case YAMLType.mapping:
            auto map = Json.emptyObject;

            foreach (pair; value!(YAML.Pair[])())
            {
                enforce!ConfigFileException(
                    pair.key.type == YAMLType.string,
                    "only string-keys are allowed",
                );

                map[pair.key.get!string] = pair.value.toJson;
            }

            return map;
        case YAMLType.sequence:
            auto seq = Json.emptyArray;

            foreach (element; value!(YAML[])())
                seq ~= element.toJson;

            return seq;
        case YAMLType.decimal:
            return Json(value!real());
        case YAMLType.binary:
        case YAMLType.timestamp:
        case YAMLType.merge:
        case YAMLType.invalid:
            throw new ConfigFileException(format!"unsupported YAML type: %s"(yaml.type));
    }
}

///
unittest {
    enum complexJson = `{
    "int": 42,
    "float": 3.1415,
    "truth": true,
    "void": null,
    "answer": "The answer is 42.",
    "list": ["Apple", "Banana", "Coconut", {"weird": "thing"}],
}`;

    auto expected = ((json) => parseJson(json))(complexJson);
    auto yamlConverted = YAMLLoader.fromString(complexJson).load().toJson();

    assert(expected == yamlConverted);
}


/// Set of size units for `toBytes` and `fromBytes`.
enum SizeUnit
{
    B,
    KiB,
    MiB,
    GiB,
    TiB,
    PiB,
    EiB,
    ZiB,
    YiB,
}


/// Base for size units.
enum size_t sizeUnitBase = 2^^10;


/// Convert `value` with `unit` to bytes.
auto toBytes(in size_t value, in SizeUnit unit)
{
    return value * sizeUnitBase^^unit;
}

/// Convert `bytes` to the smallest unit such that the value is between 0
/// and 1.
auto fromBytes(in size_t bytes)
{
    alias converToUnit = exp => tuple!("value", "unit")(
        bytes.to!double / (sizeUnitBase^^exp),
        exp,
    );

    foreach (exp; EnumMembers!SizeUnit)
    {
        if (bytes <= sizeUnitBase^^exp)
            return converToUnit(exp);
    }

    return converToUnit(SizeUnit.max);
}
