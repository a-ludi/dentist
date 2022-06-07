/+ dub.sdl:
    name "jsonschema"
    dependency "vibe-d:data" version="~>0.9.3"
+/
module dentist.util.jsonschema;

import std.format : format;
import std.traits :
    isCallable,
    isDynamicArray,
    isFloatingPoint,
    isIntegral,
    isSomeString,
    isStaticArray,
    Parameters;
import vibe.data.json : Json, toJson = serializeToJson;


/// Generate a basic JSON schema for type `T`.
///
/// This procedure has some special adaptions for DENTIST:
///
/// - `enum : bool` is interpreted as boolean
///
/// See_also: $(LINK https://json-schema.org/)
template jsonSchema(T)
{
    //("null", "boolean", "object", "array", "number", or "string"), or "integer"
    static if (is(T == bool) || (is(T B == enum) && is(B == bool)))
        enum jsonSchema = SchemaType.boolean;
    else static if (is(T == enum))
        enum jsonSchema = SchemaType.enum_([__traits(allMembers, T)].toJson);
    else static if (isSomeString!T)
        enum jsonSchema = SchemaType.string_;
    else static if (isFloatingPoint!T)
        enum jsonSchema = SchemaType.number;
    else static if (isIntegral!T)
        enum jsonSchema = SchemaType.integer;
    else static if (isStaticArray!T)
        enum jsonSchema = SchemaType.array(jsonSchema!(typeof(T.init[0])), T.length);
    else static if (isDynamicArray!T)
        enum jsonSchema = SchemaType.array(jsonSchema!(typeof(T.init[0])));
    else static if (isCallable!T && Parameters!T.length == 0)
        enum jsonSchema = SchemaType.array(jsonSchema!(typeof(T.init[0])));
    else
        static assert(0, "Cannot derive JSON schema from " ~ T.stringof);
}


unittest
{
    enum EnumA { A, B, C, d, e, f }
    enum Flag { no = false, yes = true }

    assert(jsonSchema!Flag == SchemaType.boolean);
    assert(jsonSchema!bool == SchemaType.boolean);
    assert(jsonSchema!string == SchemaType.string_);
    assert(jsonSchema!(char[]) == SchemaType.string_);
    assert(jsonSchema!float == SchemaType.number);
    assert(jsonSchema!double == SchemaType.number);
    assert(jsonSchema!EnumA == SchemaType.enum_(["A","B","C","d","e","f"].toJson));

    assert(jsonSchema!(Flag[]) == SchemaType.array(SchemaType.boolean));
    assert(jsonSchema!(bool[]) == SchemaType.array(SchemaType.boolean));
    assert(jsonSchema!(string[]) == SchemaType.array(SchemaType.string_));
    assert(jsonSchema!(char[][]) == SchemaType.array(SchemaType.string_));
    assert(jsonSchema!(float[]) == SchemaType.array(SchemaType.number));
    assert(jsonSchema!(double[]) == SchemaType.array(SchemaType.number));
    assert(jsonSchema!(EnumA[]) == SchemaType.array(SchemaType.enum_(["A","B","C","d","e","f"].toJson)));

    assert(jsonSchema!(Flag[42]) == SchemaType.array(SchemaType.boolean, 42));
    assert(jsonSchema!(bool[42]) == SchemaType.array(SchemaType.boolean, 42));
    assert(jsonSchema!(string[42]) == SchemaType.array(SchemaType.string_, 42));
    assert(jsonSchema!(char[][42]) == SchemaType.array(SchemaType.string_, 42));
    assert(jsonSchema!(float[42]) == SchemaType.array(SchemaType.number, 42));
    assert(jsonSchema!(double[42]) == SchemaType.array(SchemaType.number, 42));
    assert(jsonSchema!(EnumA[42]) == SchemaType.array(SchemaType.enum_(["A","B","C","d","e","f"].toJson), 42));
}


protected struct SchemaType
{
    enum null_ = ["type": "null"].toJson;
    enum boolean = ["type": "boolean"].toJson;
    enum number = ["type": "number"].toJson;
    enum string_ = ["type": "string"].toJson;
    enum integer = ["type": "integer"].toJson;

    static auto object(Json props)
    {
        return [
            "type": "object".toJson,
            "properties": props,
        ].toJson;
    }

    static auto array(Json items, size_t staticLength = size_t.max)
    {
        auto schema = [
            "type": "array".toJson,
            "items": items,
        ].toJson;

        if (staticLength < size_t.max)
        {
            schema["minItems"] = staticLength.toJson;
            schema["maxItems"] = staticLength.toJson;
        }

        return schema;
    }

    static auto enum_(Json items)
    {
        return ["enum": items].toJson;
    }
}


/// Get a schema reference object.
///
/// ---
/// writeln(schemaRef!"option-list");
/// // {"$ref":"#/$defs/option-list"}
/// writeln(schemaRef!("utc-date", "https://example.com"));
/// // {"$ref":"https://example.com/$defs/utc-date"}
/// ---
///
/// See_also:
///
/// - $(LINK https://json-schema.org/understanding-json-schema/structuring.html#ref)
/// - $(LINK https://json-schema.org/understanding-json-schema/structuring.html#defs)
enum schemaRef(string id, string uri = "#") = ["$ref": uri ~ "/$defs/" ~ id].toJson;
