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
    isArgumentHandler,
    isOptionHandler,
    Option,
    OptionFlag;
import dentist.commandline : OptionsFor;
import dentist.common.commands :
    DentistCommand,
    dentistCommands;
import std.algorithm :
    canFind,
    startsWith;
import std.conv : to;
import std.exception :
    basicExceptionCtors,
    enforce;
import std.format : format;
import std.meta : Alias;
import std.range :
    ElementType,
    only;
import std.stdio : File;
import std.string :
    split;
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
import std.typecons : tuple;
import vibe.data.json :
    Json,
    parseJson;


/// Identifier for the config object with default values.
enum configDefaultKey = "__default__";

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


/// Validate config file and initialize options accordingly.
Options initFromConfigFile(Options)(ref Options options, in string configFile)
{
    enum dentistCommand = Options.commandName;
    auto configContent = readConfigFile(configFile);
    auto configValues = parseJson(
        configContent,
        null,
        configFile,
    );

    validateConfig(configValues);

    auto defaultValues = configDefaultKey in configValues
        ? configValues[configDefaultKey]
        : Json.emptyObject;
    auto commandValues = dentistCommand in configValues
        ? configValues[dentistCommand]
        : Json.emptyObject;

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


void validateConfigFile(in string configFile)
{
    auto configContent = readConfigFile(configFile);
    auto configValues = parseJson(
        configContent,
        null,
        configFile,
    );

    validateConfig(configValues);
}

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

void validateConfigDefault(Json defaultConfig)
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

template configNamesOf(alias symbol)
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

void assignConfigValue(string member, Options)(ref Options options, Json configValue)
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

auto getConfigValue(SymbolType)(Json configValue)
{
    static if (is(SymbolType == OptionFlag))
        return configValue.get!bool.to!SymbolType;
    else static if (is(SymbolType == enum))
        return configValue.get!string.to!SymbolType;
    else static if (is(SymbolType == OptionFlag) || is(SymbolType : bool))
        return configValue.get!bool.to!SymbolType;
    else static if (isFloatingPoint!SymbolType)
        return configValue.get!double.to!SymbolType;
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

string readConfigFile(in string configFileName)
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
enum size_t sizeUnitBase = 2^^10;

auto toBytes(in size_t value, in SizeUnit unit)
{
    return value * sizeUnitBase^^unit;
}

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
