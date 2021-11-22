/**
    Some functions to work with FASTA data.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.fasta;

import std.algorithm :
    count,
    equal,
    find,
    joiner,
    startsWith;
import std.ascii : newline;
import std.array :
    appender, array;
import std.conv : to;
import std.exception : enforce;
import std.format :
    format,
    formattedRead;
import std.range :
    chain,
    chunks,
    drop,
    ElementType,
    empty,
    front,
    isBidirectionalRange,
    only,
    popBack,
    popFront,
    take,
    walkLength;
import std.stdio : File;
import std.string : indexOf, lineSplitter, outdent;
import std.traits : isSomeChar, isSomeString;
import std.typecons : tuple, Tuple;


/**
    Gives access to FASTA data. Does not copy the input sequence.
*/
struct Fasta(T) if (isSomeString!T)
{
    /// FASTA headers start with this character.
    static enum headerIndicator = '>';
    /// FASTA data.
    const T data;
    private size_t[] recordIndex;

    alias data this;

    /**
        Build an index in order to give fast access to individual records.
        This is called implicitly when accessing individual records using
        `opIndex` or `length`.
    */
    void buildIndex()
    {
        if ((data.length > 0 && recordIndex.length > 0) || (data.length == 0
                && recordIndex.length == 0))
            return;

        recordIndex.reserve(data.count(headerIndicator) + 1);
        long currIdx = data.indexOf(headerIndicator);

        while (currIdx >= 0)
        {
            recordIndex ~= currIdx.to!size_t;
            currIdx = data.indexOf(headerIndicator, currIdx + 1);
        }

        recordIndex ~= data.length;
    }


    /**
        Get the FASTA record at idx (zero-based).

        Returns: `FastaRecord!T` at index idx.
    */
    FastaRecord!T opIndex(size_t idx)
    {
        assert(0 <= idx && idx < length, "index out of bounds");
        buildIndex();
        auto recordBegin = recordIndex[idx];
        auto recordEnd = recordIndex[idx + 1];

        return data[recordBegin .. recordEnd].parseFastaRecord();
    }


    /// Get the number of FASTA records.
    @property size_t length()
    {
        buildIndex();

        return recordIndex.length - 1;
    }


    /// Returns true iff line starts with '>'.
    static bool isHeaderLine(in T line) pure
    {
        return line.startsWith(only(headerIndicator));
    }
}

///
unittest
{
    auto fasta1 = Fasta!string(q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
        >sequence2
        AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTG
        TAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC
EOF".outdent);
    auto fasta1Records = [
        q"EOF
            >sequence1
            CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
            AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent.parseFastaRecord,
        q"EOF
            >sequence2
            AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTG
            TAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC
EOF".outdent.parseFastaRecord,
    ];

    assert(fasta1.length == 2, fasta1.length.to!string);
    assert(fasta1[0] == fasta1Records[0]);
    assert(fasta1[1] == fasta1Records[1]);
}


/// Convenience wrapper around `Fasta!T(T data)`.
Fasta!T parseFasta(T)(T data)
{
    return typeof(return)(data);
}

///
unittest
{
    string fastaData = q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
        >sequence2
        AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTG
        TAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC
EOF".outdent;
    auto fasta = fastaData.parseFasta();
    auto fastaRecords = [
        FastaRecord!string(q"EOF
            >sequence1
            CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
            AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent),
        FastaRecord!string(q"EOF
            >sequence2
            AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTG
            TAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC
EOF".outdent),
    ];

    assert(fasta.length == 2, fasta.length.to!string);
    assert(fasta[0] == fastaRecords[0]);
    assert(fasta[1] == fastaRecords[1]);
}


/**
    Gives access to a single FASTA record. Does not copy the input sequence.
*/
struct FastaRecord(T) if (isSomeString!T)
{
    /// Unix line separator.
    static enum lineSep = "\n";

    private alias Slice = Tuple!(int, int);

    /// FASTA data.
    const T data;

    alias data this;


    /// Get this record in FASTA format.
    auto toFasta(in size_t lineWidth = 50) pure const
    {
        auto formattedBody = this[].chunks(lineWidth).joiner(lineSep);

        return chain(header, lineSep, formattedBody, lineSep);
    }


    /// Get the complete header line including the leading `>`.
    @property auto header() pure const
    {
        return data.lineSplitter.front;
    }


    /// Get the length of the sequence (in characters).
    @property size_t length() pure const
    {
        return this[].walkLength;
    }

    /// ditto
    @property size_t opDollar(size_t dim : 0)()
    {
        return length;
    }


    /// Get the sequence of this FASTA record without newlines.
    auto opIndex() pure const
    {
        return data.lineSplitter.drop(1).joiner;
    }


    /// Get the sequence character at index `i` of this FASTA record.
    auto opIndex(int i) pure const
    {
        i = normalizeIndex(i);
        assert(0 <= i && i < length,
                format!"index out of bounds: %d not in [-%d, %d)"(i, length, length));

        return this[].drop(i).front;
    }


    /// Get sub-sequence from `i` to `j` (exclusive) of this FASTA record.
    auto opIndex(in Slice slice) pure const
    {
        auto i = normalizeIndex(slice[0]);
        auto j = normalizeIndex(slice[1]);
        assert(0 <= i && i <= j && j <= length,
                format!"index out of bounds: [%d, %d) not in [-%d, %d)"(i, j, length, length));

        return this[].drop(i).take(j - i);
    }


    auto opSlice(size_t dim : 0)(int i, int j)
    {
        return tuple(i, j);
    }

    auto opSlice(size_t dim : 0)(size_t i, size_t j)
    {
        return tuple(i.to!int, j.to!int);
    }


    private int normalizeIndex(int i) const
    {
        auto length = this.length;

        while (i < 0)
            i += length;

        return i;
    }
}

///
unittest
{
    auto fastaRecord1 = q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent.parseFastaRecord;

    assert(fastaRecord1.header == ">sequence1");
    assert(fastaRecord1[].equal("CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC"));
    assert(fastaRecord1[0 .. 5].equal("CTAAC"));
    assert(fastaRecord1.toFasta(13).equal(q"EOF
        >sequence1
        CTAACCCTAACCC
        TAACCCTAACCCT
        AACCCTAACCCTA
        ACCCTAACCCTAA
        CCCTAACCCTAAC
        CCTAACCCTAACC
        CTAACAACCCTAA
        CCCTAACCC
EOF".outdent));

    auto fastaRecord2 = q"EOF
        >sequence2
        AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTG
        TAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC
EOF".outdent.parseFastaRecord;

    assert(fastaRecord2.header == ">sequence2");
    assert(fastaRecord2[].equal("AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTGTAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC"));
    assert(fastaRecord2[5 .. 10].equal("GAGCA"));
    assert(fastaRecord2.toFasta(45).equal(q"EOF
        >sequence2
        AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTG
        TATTGTAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAAC
        AACCCTAACC
EOF".outdent));
}


/// Convenience wrapper around `FastaRecord!T(T data)`.
FastaRecord!T parseFastaRecord(T)(T data)
{
    return typeof(return)(data);
}

///
unittest
{
    string fastaRecordData = q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent;
    auto fastaRecord = fastaRecordData.parseFastaRecord();

    assert(fastaRecord.header == ">sequence1");
    assert(fastaRecord[].equal("CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC"));
    assert(fastaRecord[0 .. 5].equal("CTAAC"));
}


/**
    Calculate the sequence length of the first record in `fastaFile`. Returns
    the length of the next record in `fastaFile` if it is a File object.
*/
size_t getFastaLength(in string fastaFile)
{
    enum headerIndicator = Fasta!string.headerIndicator;
    alias isHeaderLine = (line) => line.length > 0 && line[0] == headerIndicator;

    return File(fastaFile).getFastaLength();
}

/// ditto
size_t getFastaLength(File fastaFile)
{
    enum headerIndicator = Fasta!string.headerIndicator;
    alias isHeaderLine = (line) => line.length > 0 && line[0] == headerIndicator;

    auto fastaLines = fastaFile
        .byLine
        .find!isHeaderLine;

    enforce(!fastaLines.empty, "cannot determine FASTA length: file has no records");

    // ignore header
    fastaLines.popFront();

    char peek()
    {
        import core.stdc.stdio : getc, ungetc;
        import std.exception : errnoEnforce;

        auto c = getc(fastaFile.getFP());
        ungetc(c, fastaFile.getFP());

        errnoEnforce(!fastaFile.error);

        return cast(char) c;
    }

    // sum length of all sequence lines up to next record
    size_t length;

    if (peek() == headerIndicator)
        return 0;
    foreach (line; fastaLines)
    {
        length += line.length;

        if (peek() == headerIndicator)
            break;
    }

    return length;
}

///
unittest
{
    import std.process : pipe;

    string fastaRecordData = q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent;
    auto fastaFile = pipe();
    fastaFile.writeEnd.write(fastaRecordData);
    fastaFile.writeEnd.close();

    auto fastaLength = getFastaLength(fastaFile.readEnd);

    assert(fastaLength == 100);
}

///
unittest
{
    import std.process : pipe;

    string fastaRecordData = q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
        >sequence2
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent;
    auto fastaFile = pipe();
    fastaFile.writeEnd.write(fastaRecordData);
    fastaFile.writeEnd.close();

    auto fastaLength1 = getFastaLength(fastaFile.readEnd);
    auto fastaLength2 = getFastaLength(fastaFile.readEnd);

    assert(fastaLength1 == 100);
    assert(fastaLength2 == 50);
}

unittest
{
    import std.exception : assertThrown;
    import std.process : pipe;

    auto fastaFile = pipe();
    fastaFile.writeEnd.close();

    assertThrown(getFastaLength(fastaFile.readEnd));
}


/// Represents standard PacBio header format:
/// `>{smrtId}/{well}/{hqBegin}_{hqEnd} {readQuality}`
struct PacBioHeader(T) if (isSomeString!T)
{
    static private enum headerFormat = ">%s/%d/%d_%d %s";

    /// Name of the SMRT© Cell.
    T name;

    /// Index of the well where the read occurred.
    size_t well;

    /// Begin of the high quality region.
    size_t qualityRegionBegin;

    /// End of the high quality region.
    size_t qualityRegionEnd;

    /// More information, usually `RQ=0.xx`
    string additionalInformation;


    /// Construct a `PacBioHeader!T` from `header`.
    this(T header)
    {
        this.parse(header);
    }


    /// Assign new `header` data.
    void opAssign(T header)
    {
        this.parse(header);
    }


    /// Builds the header string.
    S to(S : T)() const
    {
        return buildHeader();
    }

    private T buildHeader() const
    {
        return format!headerFormat(
            name,
            well,
            qualityRegionBegin,
            qualityRegionEnd,
            additionalInformation,
        );
    }

    private void parse(in T header)
    {
        auto numMatches = header[].formattedRead!headerFormat(
            name,
            well,
            qualityRegionBegin,
            qualityRegionEnd,
            additionalInformation,
        );

        assert(numMatches == 5);
    }
}

///
unittest
{
    string header = ">name/1/0_1337 RQ=0.75";
    auto pbHeader1 = PacBioHeader!string(header);

    assert(pbHeader1.to!string == ">name/1/0_1337 RQ=0.75");
    assert(pbHeader1.name == "name");
    assert(pbHeader1.well == 1);
    assert(pbHeader1.qualityRegionBegin == 0);
    assert(pbHeader1.qualityRegionEnd == 1337);
    assert(pbHeader1.additionalInformation == "RQ=0.75");

    PacBioHeader!string pbHeader2 = header;

    assert(pbHeader2 == pbHeader1);
}


/// Convenience wrapper around `PacBioHeader!T(T header)`.
PacBioHeader!T parsePacBioHeader(T)(T header)
{
    return typeof(return)(header);
}

///
unittest
{
    string header = ">name/1/0_1337 RQ=0.75";
    auto pbHeader1 = header.parsePacBioHeader();

    assert(pbHeader1.to!string == ">name/1/0_1337 RQ=0.75");
    assert(pbHeader1.name == "name");
    assert(pbHeader1.well == 1);
    assert(pbHeader1.qualityRegionBegin == 0);
    assert(pbHeader1.qualityRegionEnd == 1337);
    assert(pbHeader1.additionalInformation == "RQ=0.75");
}


/**
    Get the complement of a DNA base. Only bases A, T, C, G (case-insensitive)
    will be translated; all other characters are left as is. Replacement
    preserves casing of the characters.
*/
C complement(C)(C base) if (isSomeChar!C)
{
    import std.range : zip;

    enum from = `AGTCagtc`;
    enum to = `TCAGtcag`;

    switch (base)
    {
        static foreach (conv; zip(from, to))
        {
            case conv[0]:
                return conv[1];
        }
        default:
            return base;
    }
}


/**
    Compute the reverse complement of a DNA sequence. Only bases A, T, C, G
    (case-insensitive) will be translated; all other characters are left as
    is. Replacement preserves casing of the characters.
*/
auto reverseComplementer(Range)(Range sequence)
        if (isBidirectionalRange!Range && isSomeChar!(ElementType!Range))
{
    import std.algorithm : map;
    import std.range : retro;

    return sequence
        .retro
        .map!complement;
}

/// ditto
T reverseComplement(T)(in T sequence) if (isSomeString!T)
{
    import std.array : array;

    return sequence[].reverseComplementer.array.to!T;
}


/// Return a copy of `fastaRecord` with reverse-complemented sequence.
FastaRecord!T reverseComplement(T)(in FastaRecord!T fastaRecord) if (isSomeString!T)
{
    enum lineSep = FastaRecord!T.lineSep;
    auto header = fastaRecord.header;
    auto sequence = fastaRecord[].array.reverseComplement;
    auto builder = appender!T;

    builder.reserve(header.length + sequence.length + 2 * lineSep.length);

    builder ~= header;
    builder ~= lineSep;
    builder ~= sequence;
    builder ~= lineSep;

    return typeof(return)(builder.data);
}

///
unittest
{
    auto seq = "GGTTGTAAATTGACTGTTGTCTGCT\ngccaatctactggtgggggagagat";
    auto revComp = "atctctcccccaccagtagattggc\nAGCAGACAACAGTCAATTTACAACC";

    assert(seq.reverseComplement == revComp);
    assert(seq.reverseComplementer.equal(revComp));

    auto fastaRecord1 = q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent.parseFastaRecord;
    auto fastaRecord1RevComp = fastaRecord1.reverseComplement;

    assert(fastaRecord1RevComp.header == ">sequence1");
    assert(fastaRecord1RevComp[].equal("GGGTTAGGGTTAGGGTTGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAG"));
    assert(fastaRecord1RevComp[0 .. 5].equal("GGGTT"));
    assert(fastaRecord1RevComp.toFasta(13).equal(q"EOF
        >sequence1
        GGGTTAGGGTTAG
        GGTTGTTAGGGTT
        AGGGTTAGGGTTA
        GGGTTAGGGTTAG
        GGTTAGGGTTAGG
        GTTAGGGTTAGGG
        TTAGGGTTAGGGT
        TAGGGTTAG
EOF".outdent));
}
