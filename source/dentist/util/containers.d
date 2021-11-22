/**
    Useful containers. Currently only `HashSet`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.containers;

import std.algorithm : all, count, filter, map;
import std.range.primitives;


/// A set of `T`s implemented using a hash-table.
///
/// In DENTIST, this is used for situations where the memory footprint of
/// `dentist.util.math.NaturalNumberSet` may get too large. They both
/// implement similar interfaces.
///
/// See_also: `hashSet`
struct HashSet(T)
{
    private bool[T] _hash;


    /// Add `value` to this set regardless whether it was present or not.
    void add(T value)
    {
        _hash[value] = true;
    }


    /// Remove `value` from this set regardless whether it was present or not.
    void remove(T value)
    {
        _hash.remove(value);
    }


    /// Return whether `value` is in this set.
    bool has(T value) const pure
    {
        return _hash.get(value, false);
    }

    /// ditto
    bool opBinaryRight(string op)(T value) const pure if (op == "in")
    {
        return this.has(value);
    }


    /// Returns true if this set is empty.
    bool empty() const pure
    {
        return _hash.values.all;
    }


    /// Remove all elements from this set.
    void clear() pure
    {
        _hash.clear();
    }


    /// Return the number of elements in this set.
    @property size_t size() const pure
    {
        return _hash.values.count;
    }


    /// Return a lazy range of the elements in this set.
    @property auto elements()
    {
        return _hash
            .byKeyValue
            .filter!"a.value"
            .map!"a.key";
    }
}


/// Construct a new `HashSet` from `values`.
static HashSet!(ElementType!R) hashSet(R)(R values) if (isInputRange!R)
{
    typeof(return) set;

    foreach (value; values)
        set.add(value);

    return set;
}

///
unittest
{
    enum numbers = [19, 6, 29, 39, 98, 86, 27, 11, 23, 59, 20, 34, 84, 38, 20];

    auto set = hashSet(numbers);

    foreach (n; numbers)
        assert(n in set);
}
