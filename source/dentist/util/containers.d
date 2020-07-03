/**
    Useful containers.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.containers;

import std.algorithm : all, count, filter, map;
import std.range.primitives;


struct HashSet(T)
{
    private bool[T] _hash;


    void add(T value)
    {
        _hash[value] = true;
    }


    void remove(T value)
    {
        _hash[value] = false;
    }


    bool has(T value) const pure
    {
        return _hash.get(value, false);
    }

    bool opBinaryRight(string op)(T value) const pure if (op == "in")
    {
        return this.has(value);
    }


    bool empty() const pure
    {
        return _hash.values.all;
    }


    void clear() pure
    {
        _hash.clear();
    }


    @property size_t size() const pure
    {
        return _hash.values.count;
    }


    @property auto elements()
    {
        return _hash
            .byKeyValue
            .filter!"a.value"
            .map!"a.key";
    }
}


static HashSet!(ElementType!R) hashSet(R)(R values) if (isInputRange!R)
{
    typeof(return) set;

    foreach (value; values)
        set.add(value);

    return set;
}

unittest
{
    enum numbers = [19, 6, 29, 39, 98, 86, 27, 11, 23, 59, 20, 34, 84, 38, 20];

    auto set = hashSet(numbers);

    foreach (n; numbers)
        assert(n in set);
}
