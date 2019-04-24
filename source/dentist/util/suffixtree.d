/**
    Implements a SuffixTree for rapid string lookups.

    Note: Implementation is based on an [example implementation][1] in C++.

    [1]: https://www.sanfoundry.com/cpp-program-implement-suffix-tree/

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.suffixtree;


// Import primitives to make arrays behave like ranges.
import std.range.primitives;

/// Thrown when an error in SuffixTree occurs.
class SuffixTreeException : Exception
{
    import std.exception : basicExceptionCtors;

    ///
    mixin basicExceptionCtors;
}

/// Thrown when a out-of-alphabet char is detected.
class AlphabetException : SuffixTreeException
{
    import std.exception : basicExceptionCtors;

    ///
    mixin basicExceptionCtors;
}

static import std.exception;
protected static alias enforce = std.exception.enforce!SuffixTreeException;


/// Type for characters in texts.
alias Byte = const(ubyte);

/// Valid character types for texts.
enum isSomeByte(T) = (
    is(T : Byte) ||
    is(T : const(byte)) ||
    is(T : const(char))
);

/// Type for texts.
alias Bytes = Byte[];

/// Valid types for texts.
enum isSomeBytes(T) = (
    is(T : Bytes) ||
    is(T : const(byte)[]) ||
    is(T : const(char)[])
);

/// Convert compatible data into `Bytes`.
inout(Bytes) toBytes(T)(inout(T) text) if (isSomeBytes!T)
{
    return cast(typeof(return)) text;
}

/// Convert compatible data back to `T` from `Bytes`.
inout(T) fromBytes(T)(inout(Bytes) text) if (isSomeBytes!T)
{
    return cast(typeof(return)) text;
}

/// Data structure to quickly find substrings in a text.
struct SuffixTree(alias inputAlphabet) if (isSomeBytes!(typeof(inputAlphabet)))
{
    static enum SpecialByte : Byte
    {
        /// Signifies end of text.
        eof = 0x00u,
        /// Used as `0` to build uniquq text id for multiple texts.
        id0 = 0x01u,
        /// Used as `1` to build uniquq text id for multiple texts.
        id1 = 0x02u,
        /// Separates multiple texts from each other.
        rs = 0x03u,
    }

    import std.algorithm : canFind;
    import std.traits : EnumMembers;
    static foreach (specialByte; [EnumMembers!SpecialByte])
        static assert(
            !inputAlphabet.canFind(specialByte),
            format!"byte 0x%x not allowed in inputAlphabet"(specialByte),
        );

    static bool isSpecialByte(in Byte char_)
    {
        import std.algorithm : among;

        return cast(bool) char_.among(cast(Byte) EnumMembers!SpecialByte);
    }

    static enum alphabet = inputAlphabet.toBytes ~ [cast(Byte) EnumMembers!SpecialByte];

    private static struct Node
    {
        size_t begin;
        size_t end;
        size_t depth;
        Node* parent;
        Node* suffixLink;
        Node*[alphabet.length] children;


        @property bool hasChildren() const
        {
            import std.algorithm : canFind;

            return children[].canFind!"a !is null";
        }

        auto ref inout(Node*) getChild(string variable = "text")(in Byte char_) inout
        {
            return children[indexOf!variable(char_)];
        }

        alias opIndex = getChild;

        const(Bytes) str(in Bytes text) const
        {
            return text[begin .. end];
        }

        @property size_t strLength() const
        {
            return end - begin;
        }

        bool contains(in size_t textIndex) const
        {
            return depth <= textIndex && textIndex < depth + strLength;
        }

        @property bool isLeaf() const
        {
            import std.algorithm : all;

            return children[].all!"a == null";
        }

        string humanStr(in Bytes text) const
        {
            import std.regex : Captures, ctRegex, replaceAll;
            import std.string : tr;

            enum terminatorRegex = ctRegex!"\x00(?:([\x01\x02]+)\x03)?";
            auto nodeString = cast(string) str(text);

            static string humanReadableTerminator(Captures!string match)
            {
                import std.format : format;

                if (match[1].length > 0)
                    return format!"$%d;"(fromTextId(match[1].toBytes));
                else
                    return "$";
            }

            return nodeString
                .replaceAll!humanReadableTerminator(terminatorRegex)
                .tr(
                    [cast(const(char)) EnumMembers!SpecialByte],
                    "$01;",
                );
        }

        string toString(in Bytes text, size_t level = 0) const
        {
            import std.algorithm : copy, joiner, map;
            import std.array : appender;
            import std.conv : to;
            import std.format : formattedWrite;
            import std.range : chain, repeat;

            auto txt = appender!string();
            auto indent = chain(
                "  ".repeat(level).joiner,
                "+ ",
            ).map!"cast(char) a";

            enum nodeFormat = `"%s" (%d)`;
            txt.formattedWrite!nodeFormat(humanStr(text), depth);
            if (suffixLink !is null)
                txt.formattedWrite!(" ~> " ~ nodeFormat)(suffixLink.humanStr(text), suffixLink.depth);

            if (isLeaf)
            {
                txt.put('\n');
                indent.copy(txt);
                "EOF".copy(txt);
            }
            else
            {
                foreach (i, child; children)
                {
                    if (child !is null)
                    {
                        txt.put('\n');
                        indent.copy(txt);

                        if (alphabet[i] == SpecialByte.eof)
                            "EOF".copy(txt);
                        else
                            txt.formattedWrite!"%c: %s"(cast(char) alphabet[i], child.toString(text, level + 1));
                    }
                }
            }

            return txt.data;
        }
    }

    private Bytes text;
    private Node* root;
    private size_t[] textEnds;

    /// Returns text (at `textIndex` if multiple texts were given).
    inout(Bytes) getText(in size_t textIndex = 0) inout
    {
        if (!hasMultipleTexts)
            return text[0 .. $ - SpecialByte.eof.sizeof];

        auto begin = textIndex > 0
            ? textEnds[textIndex - 1] + textSeparatorSize
            : 0;
        auto end = textEnds[textIndex];

        return text[begin .. end];
    }

    /// Returns the number of texts in this suffix tree.
    @property size_t numTexts() const
    {
        return hasMultipleTexts
            ? textEnds.length
            : 1;
    }

    /// Returns true if this suffix tree was created with multiple texts.
    @property bool hasMultipleTexts() const
    {
        return textEnds.length > 0;
    }

    import std.range.primitives : ElementType, isForwardRange, isInputRange;
    /**
        Build a suffix tree from a single `inputText`. The text must only
        contain characters in `alphabet`.

        Throws: AlphabetException if `inputText` contains non-alphabet characters
    */
    static SuffixTree build(S)(S inputText)
        if (
            isSomeBytes!S ||
            (isInputRange!S && isSomeByte!(ElementType!S) && hasLength!S)
        )
    {
        static enum Rule
        {
            init,
            continueNode,
            newSuffix,
            nextChar,
        }

        static struct NodesAllocator
        {
            Node[] nodes;

            this(in size_t textLength)
            {
                this.nodes.length = 2 * textLength;
            }

            Node* create(Params...)(Params params)
            {
                // Take reference to current node and ...
                auto allocNode = &nodes[0];
                // ... shift the node off the buffer
                nodes = nodes[1 .. $];

                *allocNode = Node(params);

                return allocNode;
            }
        }

        import std.algorithm : map;
        import std.array : array;
        import std.range : chain, only, takeExactly;

        Bytes text = chain(inputText, only(SpecialByte.eof))
            .takeExactly(inputText.length + 1)
            .map!(c => cast(Byte) c)
            .array;
        auto nodesAllocator = NodesAllocator(text.length);
        auto suffixTree = SuffixTree(text, nodesAllocator.create());
        suffixTree.root.suffixLink = suffixTree.root;

        Node* currentNode = suffixTree.root;
        assert(currentNode !is null, "currentNode must not be null");
        Node* needsSuffixLink;

        Rule lastRule;
        // Base index for `text`
        size_t j = 0;
        alias indexText = i => indexOf(text[i]);

        foreach (i, char_; text)
        {
            size_t charIndex = indexOf(char_);

            for (; j <= i; j++)
            {
                size_t currentDepth = i - j;

                if (lastRule != Rule.nextChar)
                {
                    currentNode = currentNode.suffixLink is null
                        ? currentNode.parent.suffixLink
                        : currentNode.suffixLink;
                    assert(currentNode !is null, "currentNode must not be null");

                    if (currentDepth > 0)
                        for (
                            auto k = j + currentNode.depth + currentNode.strLength;
                            !currentNode.contains(currentDepth - 1);
                            k += currentNode.strLength
                        )
                            currentNode = currentNode.getChild(text[k]);
                    assert(currentNode !is null, "currentNode must not be null");
                }

                if (!currentNode.contains(currentDepth))
                {
                    if (needsSuffixLink !is null)
                    {
                        needsSuffixLink.suffixLink = currentNode;
                        needsSuffixLink = null;
                    }

                    if (currentNode.getChild(char_) is null)
                    {
                        currentNode.getChild(char_) = nodesAllocator.create(
                            i,
                            text.length,
                            currentDepth,
                            currentNode,
                        );
                        lastRule = Rule.newSuffix;
                    }
                    else
                    {
                        currentNode = currentNode.getChild(char_);
                        assert(currentNode !is null, "currentNode must not be null");

                        lastRule = Rule.nextChar;
                        break;
                    }
                }
                else
                {
                    size_t end = currentNode.begin + currentDepth - currentNode.depth;

                    if (text[end] != char_)
                    {
                        auto newNode = nodesAllocator.create(
                            currentNode.begin,
                            end,
                            currentNode.depth,
                            currentNode.parent,
                        );
                        newNode.getChild(char_) = nodesAllocator.create(
                            i,
                            text.length,
                            currentDepth,
                            newNode
                        );

                        newNode.getChild(text[end]) = currentNode;
                        currentNode.parent.getChild(text[currentNode.begin]) = newNode;

                        if (needsSuffixLink !is null)
                            needsSuffixLink.suffixLink = newNode;

                        currentNode.begin = end;
                        currentNode.depth = currentDepth;
                        currentNode.parent = newNode;
                        currentNode = newNode;
                        assert(currentNode !is null, "currentNode must not be null");

                        needsSuffixLink = newNode;
                        lastRule = Rule.newSuffix;
                    }
                    else if (currentNode.end != text.length || (currentNode.begin - currentNode.depth) < j)
                    {
                        lastRule = Rule.nextChar;
                        break;
                    }
                    else
                    {
                        lastRule = Rule.continueNode;
                    }
                }
            }
        }

        suffixTree.root.suffixLink = null;

        return suffixTree;
    }

    ///
    unittest
    {
        auto suffixTree = SuffixTree!"abn"().build("banana");

        assert(suffixTree.findFirst("banana"));
    }


    /**
        Build a suffix tree from multiple `inputTexts`. The texts must only
        contain characters in `alphabet`. Lookup operations return the index
        of the text where a hit was found according to the order of this
        input.

        Throws: AlphabetException if `inputTexts` contain non-alphabet characters
    */
    static SuffixTree build(SS)(SS inputTexts)
        if (
            isForwardRange!SS &&
            (
                isSomeBytes!(ElementType!SS) ||
                (isInputRange!(ElementType!SS) && isSomeByte!(ElementType!(ElementType!SS)) && hasLength!(ElementType!SS))
            )
        )
    {
        import std.algorithm : cumulativeFold, joiner, map, sum;
        import std.array : array;
        import std.range : chain, enumerate, only, takeExactly;
        import std.range.primitives : walkLength;

        auto numTexts = inputTexts.save.walkLength();
        auto totalTextsLength = inputTexts.save.map!"a.length".sum();
        auto preparedTextLength = totalTextsLength + numTexts * textSeparatorSize;
        auto preparedText = inputTexts
            .save
            .enumerate
            .map!(enumText => chain(enumText.value, textSeparator(enumText.index))
                .map!(c => cast(Byte) c))
            .joiner
            .takeExactly(preparedTextLength);

        auto suffixTree = build(preparedText);

        suffixTree.textEnds = inputTexts
            .save
            .map!(text => text.length + textSeparatorSize)
            .cumulativeFold!"a + b"(0UL)
            .map!(textIndex => textIndex - textSeparatorSize)
            .takeExactly(numTexts)
            .array;

        return suffixTree;
    }

    private static string mixinTextId()
    {
        import std.range : repeat;
        import std.string : join;

        enum numBytes = 8 * size_t.sizeof;

        return `
            import std.typecons : Tuple;

            private alias TextId = Tuple!(` ~
            "    ubyte,\n".repeat(numBytes).join ~
            `);`;
    }
    mixin(mixinTextId);

    private static auto textSeparator(in size_t textIndex)
    {
        import std.range : only;

        return only(SpecialByte.eof, toTextId(textIndex).expand, SpecialByte.rs);
    }

    private static enum textSeparatorSize = SpecialByte.eof.sizeof + TextId.length + SpecialByte.rs.sizeof;

    private static TextId toTextId(in size_t recordIndex)
    {
        TextId textId;

        foreach (i, ref idByte; textId)
            idByte = (recordIndex >> i & 1) == 1
                ? SpecialByte.id1
                : SpecialByte.id0;

        return textId;
    }

    private static size_t fromTextId(in TextId textId)
    {
        size_t recordIndex;

        foreach (i, idByte; textId)
            if (idByte == SpecialByte.id1)
                recordIndex |= 1UL << i;

        return recordIndex;
    }

    private static size_t fromTextId(in Bytes textId)
    {
        size_t recordIndex;

        assert(textId.length == TextId.length, "trying to convert textId with wrong length");
        foreach (i, idByte; textId)
            if (idByte == SpecialByte.id1)
                recordIndex |= 1UL << i;

        return recordIndex;
    }

    unittest
    {
        with (SpecialByte)
        {
            enum recordIndex = 0b01010101;
            enum textId = TextId(
                id1, id0, id1, id0, id1, id0, id1, id0,
                id0, id0, id0, id0, id0, id0, id0, id0,
                id0, id0, id0, id0, id0, id0, id0, id0,
                id0, id0, id0, id0, id0, id0, id0, id0,
                id0, id0, id0, id0, id0, id0, id0, id0,
                id0, id0, id0, id0, id0, id0, id0, id0,
                id0, id0, id0, id0, id0, id0, id0, id0,
                id0, id0, id0, id0, id0, id0, id0, id0,
            );

            assert(toTextId(recordIndex) == textId);
            assert(fromTextId(textId) == recordIndex);
        }
    }

    private size_t getTextIndex(in size_t hitEnd) const
    {
        import std.algorithm : countUntil, map;
        import std.range : chain, only, slide;
        import std.typecons : No, tuple;

        auto textBoundaries = chain(only(-textSeparatorSize), textEnds)
            .slide!(No.withPartial)(2)
            .map!(pair => tuple!("begin", "end")(pair[0] + textSeparatorSize, pair[1]));
        auto recordIndex = textBoundaries.countUntil!(txtBound =>
            txtBound.begin <= hitEnd && hitEnd <= txtBound.end);

        assert(recordIndex >= 0, "invalid hitEnd: not contained in any text");

        return cast(size_t) recordIndex;
    }

    size_t memorySize() const
    {
        return memorySize(text.length);
    }

    static size_t memorySize(in size_t totalTextLength, in size_t numTexts = 1)
    {
        auto preparedTextLength = totalTextLength + numTexts * textSeparatorSize;

        return SuffixTree.sizeof +
               preparedTextLength * Byte.sizeof +
               2 * preparedTextLength * Node.sizeof;
    }

    string toString() const
    {
        return root.toString(text);
    }

    /// Find the first occurrence of `needle` in this suffix tree.
    FindResult findFirst(B)(in B needle) const if (isSomeBytes!B)
    {
        import std.typecons : No;

        auto results = findRecursive!(No.findAll)(root, needle.toBytes(), needle.toBytes());

        if (results.length == 0)
            return FindResult.notFound;

        auto result = results[0];

        finishFindResult(result, needle.toBytes());

        return result;
    }

    /// Find substrings of a single text.
    unittest
    {
        auto suffixTree = SuffixTree!"abn"().build("banana");

        assert(suffixTree.findFirst("banana"));
        assert(suffixTree.findFirst("ana"));
        assert(suffixTree.findFirst("ban"));
        assert(!suffixTree.findFirst("naab"));
    }

    unittest
    {
        import std.exception : assertThrown;

        auto suffixTree = SuffixTree!"abn"().build("banana");

        assert(suffixTree.findFirst("banana"));
        assert(suffixTree.findFirst("anana"));
        assert(suffixTree.findFirst("nana"));
        assert(suffixTree.findFirst("ana"));
        assert(suffixTree.findFirst("na"));
        assert(suffixTree.findFirst("a"));
        assertThrown!SuffixTreeException(suffixTree.findFirst(""));
        assert(suffixTree.findFirst("b"));
        assert(suffixTree.findFirst("ba"));
        assert(suffixTree.findFirst("ban"));
        assert(suffixTree.findFirst("bana"));
        assert(suffixTree.findFirst("banan"));
        assert(!suffixTree.findFirst("naab"));
    }

    /// Find substrings of in multiple texts.
    unittest
    {
        import std.range : only;

        auto words = only(
            "banana",
            "apple",
            "strawberry",
        );
        auto suffixTree = SuffixTree!"abelnprstwy"().build(words);

        // Returns matched textId
        assert(suffixTree.findFirst("banana").textId == 0);
        assert(suffixTree.findFirst("apple").textId == 1);
        assert(suffixTree.findFirst("strawberry").textId == 2);

        // Coordinates are relative to the beginning of the designated text.
        assert(suffixTree.findFirst("berry").begin == 5);

        // Cannot find substrings of the concatenation of different words.
        assert(!suffixTree.findFirst("applestraw"));
        assert(!suffixTree.findFirst("berryapple"));
    }

    unittest
    {
        import std.range : only;

        auto words = only(
            "banana",
            "apple",
            "strawberry",
        );
        auto suffixTree = SuffixTree!"abelnprstwy"().build(words);

        auto bananaResult = suffixTree.findFirst("banana");
        assert(bananaResult);
        assert(bananaResult.textId == 0);
        assert(bananaResult.begin == 0);
        assert(bananaResult.end == "banana".length);
        auto appleResult = suffixTree.findFirst("apple");
        assert(appleResult);
        assert(appleResult.textId == 1);
        assert(appleResult.begin == 0);
        assert(appleResult.end == "apple".length);
        auto strawberryResult = suffixTree.findFirst("strawberry");
        assert(strawberryResult);
        assert(strawberryResult.textId == 2);
        assert(strawberryResult.begin == 0);
        assert(strawberryResult.end == "strawberry".length);
        assert(!suffixTree.findFirst("strawapple"));
        assert(!suffixTree.findFirst("bananaapple"));
    }

    /// Find all non-overlapping occurrences of `needle` in this suffix tree.
    FindResult[] findAll(B)(in B needle) const if (isSomeBytes!B)
    {
        import std.typecons : Yes;

        auto results = findRecursive!(Yes.findAll)(root, needle.toBytes(), needle.toBytes());

        foreach (ref result; results)
            finishFindResult(result, needle.toBytes());

        return results;
    }

    /// Find all substrings of a single text.
    unittest
    {
        auto suffixTree = SuffixTree!"abn"().build("bananana");

        assert(suffixTree.findAll("banana").length == 1);
        assert(suffixTree.findAll("a").length == 4);
        assert(suffixTree.findAll("ana").length == 2);
        assert(suffixTree.findAll("ban").length == 1);
        assert(suffixTree.findAll("naab").length == 0);
    }

    /// Find all substrings of a multiple texts.
    unittest
    {
        import std.range : only;

        auto words = only(
            "banana",
            "apple",
            "strawberry",
            "stranberry",
        );
        auto suffixTree = SuffixTree!"abelnprstwy"().build(words);

        assert(suffixTree.findAll("a").length == 6);
        assert(suffixTree.findAll("e").length == 3);
        assert(suffixTree.findAll("an").length == 3);
    }

    unittest
    {
        import std.exception : assertThrown;

        auto suffixTree = SuffixTree!"a"().build("aaaaaa");

        assert(suffixTree.findAll("aaaaaa").length == 1);
        assert(suffixTree.findAll("aaaaa").length == 1);
        assert(suffixTree.findAll("aaaa").length == 1);
        assert(suffixTree.findAll("aaa").length == 2);
        assert(suffixTree.findAll("aa").length == 3);
        assert(suffixTree.findAll("a").length == 6);
        assertThrown!SuffixTreeException(suffixTree.findAll(""));
    }

    import std.typecons : Flag;

    private FindResult[] findRecursive(Flag!"findAll" findAll)(
        in Node* currentNode,
        in Bytes currentNeedle,
        in Bytes originalNeedle,
    ) const
    {
        import dentist.util.algorithm : sliceUntil;
        import std.algorithm : min, startsWith;
        import std.typecons : No, Yes;

        enforce(originalNeedle.length > 0, "cowardly refusing to search for empty string");
        assert(currentNeedle.length <= originalNeedle.length);

        if (currentNode is null)
            return [];

        auto currentNodeStr = currentNode.str(text).sliceUntil!isSpecialByte;
        auto croppedStr = currentNodeStr[0 .. min(currentNeedle.length, $)];
        if (currentNeedle.startsWith(croppedStr))
        {
            if (croppedStr.length == currentNeedle.length)
            {
                auto results = [FindResult(currentNode.begin + croppedStr.length)];

                static if (findAll)
                    results ~= findRestRecursive(
                        currentNode,
                        currentNode.begin + croppedStr.length,
                        originalNeedle,
                    );

                return results;
            }
            else
            {
                auto nextChar = currentNeedle[currentNodeStr.length];
                auto nextNode = currentNode.getChild!"needle"(nextChar);

                return findRecursive!(Yes.findAll)(nextNode, currentNeedle[currentNodeStr.length .. $], originalNeedle);
            }
        }

        return [];
    }

    private FindResult[] findRestRecursive(
        in Node* currentNode,
        in size_t begin,
        in Bytes originalNeedle,
    ) const
    {
        import std.algorithm : filter,find, min, minElement, startsWith;
        import std.typecons : No, Yes;

        assert(0 <= originalNeedle.length && originalNeedle.length <= originalNeedle.length);

        assert(begin >= currentNode.begin, "should not be called on this node");
        auto currentNodeStr = currentNode.str(text);
        auto offset = min(begin - currentNode.begin, currentNodeStr.length);
        const(Node)* nextNode;

        for (
            Bytes nodeSubstr = currentNodeStr[offset .. $];
            nodeSubstr.length > 0;
            nodeSubstr = nodeSubstr[1 .. $].find(originalNeedle[0])
        )
        {
            auto croppedStr = nodeSubstr[0 .. min(originalNeedle.length, $)];

            if (croppedStr.length > 0 && originalNeedle.startsWith(croppedStr))
            {
                if (croppedStr.length == originalNeedle.length)
                {
                    auto matchEnd = currentNode.begin + currentNodeStr.length - nodeSubstr.length + croppedStr.length;

                    return [FindResult(matchEnd)] ~ findRestRecursive(
                        currentNode,
                        matchEnd + 1,
                        originalNeedle,
                    );
                }
                else
                {
                    auto nextChar = originalNeedle[croppedStr.length];
                    nextNode = currentNode.getChild!"needle"(nextChar);

                    return findRecursive!(Yes.findAll)(
                        nextNode,
                        originalNeedle[croppedStr.length .. $],
                        originalNeedle,
                    );
                }
            }
        }

        if (!currentNode.hasChildren)
            return [];

        const(Node)* nearestChild = currentNode
            .children[]
            .filter!(child => child !is null)
            .minElement!(child => child.begin);
        assert(currentNode.end <= nearestChild.begin);

        return findRestRecursive(nearestChild, nearestChild.begin, originalNeedle);
    }

    private void finishFindResult(ref FindResult result, in Bytes needle) const
    {
        assert(result, "cannot finish invalid result");

        result._length = needle.length;

        if (hasMultipleTexts)
        {
            result._textId = getTextIndex(result.end);
            result._end -= result._textId > 0
                ? textEnds[result._textId - 1] + textSeparatorSize
                : 0;
        }
    }

    private enum alphabetIndex = indexAlphabet(alphabet);
    private static auto indexAlphabet(in Bytes alphabet)
    {
        enum numBytes = 2^^(8 * Byte.sizeof);
        long[numBytes] index;

        index[] = -1;
        foreach (i, char_; alphabet)
            index[char_] = i;

        return index;
    }

    /// Find oeprations return this type.
    struct FindResult
    {
        private static enum notFound = FindResult(0, size_t.max);

        private size_t _end;
        private size_t _length;
        private size_t _textId;

        /// Start index of the hit relative to the beginning of the text.
        @property size_t begin() const
        {
            return _end - _length;
        }

        /// End index of the hit relative to the beginning of the text.
        @property size_t end() const
        {
            return _end;
        }

        /// Length of the hit.
        @property size_t length() const
        {
            return _length;
        }

        /// Index of the text where the hit was found.
        @property size_t textId() const
        {
            return _textId;
        }

        private @property size_t length(size_t length)
        {
            return this._length = length;
        }

        /// Evaluates to true if a hit was found.
        bool opCast(T : bool)() const
        {
            return begin <= end;
        }

        string toString() const
        {
            import std.format;

            return format!(FindResult.stringof ~ "(%d, %d, %d)")(begin, end, textId);
        }
    }

    private static size_t indexOf(string variable = "text")(in Byte char_)
    {
        import std.format : format;

        switch (char_)
        {
            static foreach (index, alphabetChar; alphabet)
                case alphabetChar:
                    return index;
            default:
                std.exception.enforce!AlphabetException(
                    0,
                    format!"%s contains illegal character: 0x%x (%c)"(variable, char_, cast(char) char_),
                );
                assert(0);
        }
    }
}

unittest
{
    /// Trigger unit tests for SuffixTree
    SuffixTree!"" dummy;
}

/// Multiple suffix trees combined for parallel usage.
struct MultiSuffixTree(alias inputAlphabet) if (isSomeBytes!(typeof(inputAlphabet)))
{
    private SuffixTree!inputAlphabet[] suffixTrees;

    /**
        Build `numParallel` suffix trees from multiple `inputTexts`
        distributing texts equally across them.

        Throws: AlphabetException if `inputTexts` contain non-alphabet characters
        See_also: SuffixTree.build
    */
    static MultiSuffixTree build(SS)(SS inputTexts, in size_t numParallel)
        if (
            isForwardRange!SS &&
            (
                isSomeBytes!(ElementType!SS) ||
                (isInputRange!(ElementType!SS) && isSomeByte!(ElementType!(ElementType!SS)) && hasLength!(ElementType!SS))
            )
        )
    {
        import dentist.util.math : ceildiv;
        import std.algorithm : cumulativeFold, group, joiner, map, min, sum;
        import std.array : array;
        import std.parallelism : parallel;
        import std.range : chain, dropExactly, enumerate, only, takeExactly;
        import std.range.primitives : walkLength;
        import std.typecons : tuple;

        alias TextBucket = typeof(inputTexts.dropExactly(0).takeExactly(0));

        auto textsLengths = inputTexts.save.map!"a.length";
        auto totalTextsLength = textsLengths.save.sum;
        auto expectedBlockSize = ceildiv(totalTextsLength, numParallel);

        auto textBuckets = new TextBucket[numParallel];
        alias oneLess = n => n > 0 ? n - 1 : 0;
        auto bucketInfoList = textsLengths
            .save
            .cumulativeFold!"a + b"
            .map!(textAmount => min(oneLess(textAmount) / expectedBlockSize, numParallel))
            .group;

        size_t numTextsInBuckets;
        foreach (bucketInfo; bucketInfoList)
        {
            auto bucketIndex = bucketInfo[0];
            auto bucketSize = bucketInfo[1];

            textBuckets[bucketIndex] = inputTexts
                .dropExactly(numTextsInBuckets)
                .takeExactly(bucketSize);

            numTextsInBuckets += bucketSize;
        }

        MultiSuffixTree multiTree;
        multiTree.suffixTrees.length = numParallel;

        foreach (i, textBucket; parallel(textBuckets[]))
            multiTree.suffixTrees[i] = SuffixTree!inputAlphabet.build(textBucket);

        return multiTree;
    }

    unittest
    {
        import std.range : only;

        auto words = only(
            "banana",
            "apple",
            "strawberry",
        );
        auto mSuffixTree = MultiSuffixTree!(2, "abelnprstwy")().build(words);

        assert(mSuffixTree.suffixTrees[0].numTexts == 2);
        assert(mSuffixTree.suffixTrees[1].numTexts == 1);
    }

    static private alias FindResult = SuffixTree!inputAlphabet.FindResult;

    /// Find the first occurrence of `needle` in this suffix tree.
    FindResult findFirst(B)(in B needle) const if (isSomeBytes!B)
    {
        import std.algorithm : find, map;

        auto hits = suffixTrees[]
            .map!(suffixTree => suffixTree.findFirst(needle))
            .find!"a";

        if (hits.empty)
            return FindResult.notFound;

        return hits.front;
    }

    /// Find all non-overlapping occurrences of `needle` in this suffix tree.
    FindResult[] findAll(B)(in B needle) const if (isSomeBytes!B)
    {
        import std.algorithm : joiner, map;
        import std.array : array;

        return suffixTrees[]
            .map!(suffixTree => suffixTree.findAll(needle))
            .joiner
            .array;
    }

    /// Find all substrings of a single text.
    unittest
    {
        import std.range : only;

        auto words = only(
            "banana",
            "apple",
            "strawberry",
        );
        auto mSuffixTree = MultiSuffixTree!(2, "abelnprstwy")().build(words);

        assert(mSuffixTree.findAll("banana").length == 1);
        assert(mSuffixTree.findAll("a").length == 5);
    }
}

unittest
{
    /// Trigger unit tests for SuffixTree
    MultiSuffixTree!(0, "") dummy;
}

version (benchmark)
{
    import std.datetime.stopwatch : benchmark;
    import std.stdio;

    enum alphabet = "acgt";

    void main()
    {
        benchmarkMultiIndexing();
        benchmarkIndexing();
        //benchmarkFindFirst();
    }

    void benchmarkIndexing()
    {
        foreach (genomeSize; [
                 1_000,
                 2_000,
                 3_000,
                 4_000,
                 5_000,
                 6_000,
                 7_000,
                 8_000,
                 9_000,
                10_000,
                20_000,
                30_000,
                40_000,
                50_000,
                60_000,
                70_000,
                80_000,
                90_000,
               100_000,
               200_000,
               300_000,
               400_000,
               500_000,
               600_000,
               700_000,
               800_000,
               900_000,
             1_000_000,
             5_000_000,
            10_000_000,
        ])
        {
            auto genome = getGenome(genomeSize);
            auto result = benchmark!(() => SuffixTree!alphabet.build(genome))(3);

            writefln!"build(len=%,d) took %s; index size %,d"(
                genomeSize,
                result[0],
                SuffixTree!alphabet.memorySize(genomeSize),
            );
        }
    }

    void benchmarkFindFirst()
    {
        import core.time : Duration;
        import std.algorithm : min;
        import std.array : array;
        import std.random : uniform;

        enum genomeSize = 1_000_000;

        auto genome = getGenome(genomeSize).array.toBytes;
        auto suffixTree = SuffixTree!alphabet.build(genome);

        enum numIterations = 1_000;
        enum expectedSubseqLength = 200_000;

        Duration totalRuntime;
        foreach (_; 0 .. numIterations)
        {
            auto begin = uniform(0, genomeSize);
            auto length = min(genomeSize - begin, randGeometric(1.0/expectedSubseqLength));

            void findSubseq()
            {
                auto findResult = suffixTree.findFirst(genome[begin .. begin + length]);

                assert(findResult.begin == begin);
                assert(findResult.length == length);
            }
            auto benchmarkResult = benchmark!findSubseq(3);

            totalRuntime += benchmarkResult[0];
        }

        writefln!"findFirst(len ~ Geo(1/%d)) took %s"(expectedSubseqLength, totalRuntime / numIterations);
    }

    void benchmarkMultiIndexing()
    {
        import std.algorithm : cumulativeFold, map, sum, until;
        import std.array : array;
        import std.parallelism : taskPool;
        import std.range : generate, tee;
        import std.typecons : No, tuple;

        foreach (genomeSize; [
                50_000,
                60_000,
                70_000,
                80_000,
                90_000,
               100_000,
               200_000,
               300_000,
               400_000,
               500_000,
               600_000,
               700_000,
               800_000,
               900_000,
             1_000_000,
             5_000_000,
            10_000_000,
        ])
        {
            enum expectedContigSize = 50_000;

            auto genome = generate!(() => getGenome(randGeometric(1.0/expectedContigSize)))
                .cumulativeFold!((acc, contig) => tuple(acc[0] + contig.length, contig))(tuple(0UL, getGenome(0)))
                .until!(acc => acc[0] > genomeSize)(No.openRight)
                .map!(acc => acc[1].array)
                .array;
            auto realGenomeSize = genome.map!"a.length".sum;
            auto numContigs = genome.length;
            auto numParallel = taskPool.size + 1;
            size_t totalMemorySize;
            auto result = benchmark!(() => {
                auto index = MultiSuffixTree!alphabet.build(genome, numParallel);

                totalMemorySize = index
                    .suffixTrees[]
                    .map!(tree => tree.memorySize())
                    .sum;
            }())(3);

            writefln!"build(totlen=%,d,numtxts=%,d,parallel=%d) took %s; index size %,d"(
                realGenomeSize,
                numContigs,
                numParallel,
                result[0],
                totalMemorySize,
            );
        }
    }

    auto getGenome(in size_t genomeSize)
    {
        import std.random : choice;
        import std.range : generate, only, takeExactly;

        return generate!(() => choice(only('a', 'c', 'g', 't'))).takeExactly(genomeSize);
    }

    size_t randGeometric(in double p)
    {
        import std.math : floor, log, log1p;
        import std.random : uniform, uniform01;

        return cast(size_t) floor(log(uniform01()) / log1p(-p));
    }
}
else
{
    version (Have_dentist) { } else:
    version (unittest) { } else:

    void main(in string[] args)
    {
        import dentist.util.region;
        import std.algorithm;
        import std.array;
        import std.exception;
        import std.range;
        import std.stdio;

        enforce(args.length == 3, "usage: ./suffixtree <reference> <query>");

        auto referenceFile = File(args[1]);
        auto suffixTree = MultiSuffixTree!(8, "acgt").build(referenceFile.byLine.map!(line => line.filter!(c => c.among('a', 'c', 'g', 't')).map!"cast(char) a".array).array);
        alias MatchedRegion = Region!(size_t, size_t, "textId");

        foreach (i, query; File(args[2]).byLine.enumerate)
            if (suffixTree.findFirst(query))
            {
                MatchedRegion m;

                foreach (hit; suffixTree.findAll(query))
                {
                    if (hit)
                    {
                        auto hitInterval = MatchedRegion.TaggedInterval(
                            hit.textId,
                            hit.begin,
                            hit.end,
                        );

                        if (!empty(m & hitInterval))
                            stderr.writefln!"overlap: %03d: %s"(i, hit);

                        writefln!"%03d: %s"(i, hit);
                        m |= hitInterval;
                    }
                    else
                    {
                        stderr.writefln!"%03d: %s"(i, hit);
                    }
                }
            }
    }
}

