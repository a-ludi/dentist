/**
    Defines bindings and utilities to/for the mummer commands.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.mummer;

import dentist.common.alignments :
    AlignmentChain,
    coord_t,
    diff_t,
    id_t;
import dentist.util.log;
import std.algorithm :
    endsWith,
    filter,
    map,
    sort,
    swap;
import std.array : array;
import std.conv : ConvException, to;
import std.format : format, formattedRead;
import std.math : round;
import std.process :
    Config,
    kill,
    pipeProcess,
    ProcessPipes,
    Redirect,
    wait;
import std.range :
    enumerate,
    only;
import std.range.primitives;
import std.string : outdent;
import std.traits : isSomeString;
import std.uni : isNumber;
import vibe.data.json : toJson = serializeToJson;

debug import std.stdio : writeln;

version (unittest)
{
    enum testDump = "
        /home/alu/.project-spaces/project_ludwig_pb_gaps/data/debug/data/reference.fasta /home/alu/.project-spaces/project_ludwig_pb_gaps/data/debug/results/assembly.fasta
        NUCMER

        [S1]\t[E1]\t[S2]\t[E2]\t[LEN 1]\t[LEN 2]\t[% IDY]\t[LEN R]\t[LEN Q]\t[TAGS]
        1\t8300\t1\t8300\t8300\t8300\t100.00\t837500\t405147\treference/1/0_1337\tscaffold-1
        12401\t169892\t12401\t169897\t157492\t157497\t99.99\t837500\t405147\treference/1/0_1337\tscaffold-1
        174901\t235150\t174898\t235147\t60250\t60250\t99.95\t837500\t405147\treference/1/0_1337\tscaffold-1
        238751\t260150\t238748\t260147\t21400\t21400\t100.00\t837500\t405147\treference/1/0_1337\tscaffold-1
        263651\t405150\t263648\t405147\t141500\t141500\t100.00\t837500\t405147\treference/1/0_1337\tscaffold-1
        408901\t837500\t1\t428592\t428600\t428592\t99.98\t837500\t428592\treference/1/0_1337\tscaffold-12
        61318\t62308\t15279679\t15278653\t991\t1027\t83.40\t15477454\t15465771\treference/2/0_1337\tscaffold-13
    ".outdent;

    alias Contig = AlignmentChain.Contig;
    enum emptyFlags = AlignmentChain.emptyFlags;
    alias Flags = AlignmentChain.Flags;
    enum complement = AlignmentChain.Flag.complement;
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias Locus = LocalAlignment.Locus;

    enum expectedAlignments = [
        AlignmentChain(
            0,
            Contig(1, 837500),
            Contig(1, 405147),
            emptyFlags,
            [LocalAlignment(
                Locus(0, 8300),
                Locus(0, 8300),
                0,
            )]
        ),
        AlignmentChain(
            1,
            Contig(1, 837500),
            Contig(1, 405147),
            emptyFlags,
            [LocalAlignment(
                Locus(12400, 169892),
                Locus(12400, 169897),
                16,
            )]
        ),
        AlignmentChain(
            2,
            Contig(1, 837500),
            Contig(1, 405147),
            emptyFlags,
            [LocalAlignment(
                Locus(174900, 235150),
                Locus(174897, 235147),
                30,
            )]
        ),
        AlignmentChain(
            3,
            Contig(1, 837500),
            Contig(1, 405147),
            emptyFlags,
            [LocalAlignment(
                Locus(238750, 260150),
                Locus(238747, 260147),
                0,
            )]
        ),
        AlignmentChain(
            4,
            Contig(1, 837500),
            Contig(1, 405147),
            emptyFlags,
            [LocalAlignment(
                Locus(263650, 405150),
                Locus(263647, 405147),
                0,
            )]
        ),
        AlignmentChain(
            5,
            Contig(1, 837500),
            Contig(2, 428592),
            emptyFlags,
            [LocalAlignment(
                Locus(408900, 837500),
                Locus(0, 428592),
                86,
            )]
        ),
        AlignmentChain(
            6,
            Contig(2, 15477454),
            Contig(3, 15465771),
            Flags(complement),
            [LocalAlignment(
                Locus(61317, 62308),
                Locus(186092, 187119),
                165,
            )]
        ),
    ];
}

/**
 * Thrown on conversion errors.
 */
class MummerException : Exception
{
    import std.exception : basicExceptionCtors;
    ///
    mixin basicExceptionCtors;
}

AlignmentChain[] getAlignments(in string deltaFile)
{
    auto alignmentChains = readCoordsDump(showCoords(deltaFile)).array;
    alignmentChains.sort!"a < b";

    return alignmentChains;
}

auto readCoordsDump(Range)(Range coordsDump)
{
    id_t numAlignmentChains;
    id_t numAScaffolds = 1;
    id_t numBScaffolds = 1;
    id_t[string] aScaffoldId;
    id_t[string] bScaffoldId;

    return coordsDump
        .enumerate(1)
        .filter!(enumLine => enumLine.value.length > 0 && enumLine.value[0].isNumber)
        .map!(enumLine => {
            auto lineNumber = enumLine.index;
            auto line = enumLine.value;

            coord_t contigAAlignmentLength;
            coord_t contigBAlignmentLength;
            float identity;
            string aHeader;
            string bHeader;
            AlignmentChain alignmentChain;
            alignmentChain.id = numAlignmentChains++;
            alignmentChain.localAlignments.length = 1;

            try
            {
                line.formattedRead!"%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%s\t%s"(
                    alignmentChain.localAlignments[0].contigA.begin,
                    alignmentChain.localAlignments[0].contigA.end,
                    alignmentChain.localAlignments[0].contigB.begin,
                    alignmentChain.localAlignments[0].contigB.end,
                    contigAAlignmentLength,
                    contigBAlignmentLength,
                    identity,
                    alignmentChain.contigA.length,
                    alignmentChain.contigB.length,
                    aHeader,
                    bHeader,
                );
            }
            catch (ConvException e)
            {
                debug writeln(line);

                throw new MummerException(
                    format!"format error in coords dump on line %d: %s"(lineNumber, e.msg),
                    e,
                );
            }

            if (aHeader !in aScaffoldId)
                aScaffoldId[aHeader] = numAScaffolds++;
            alignmentChain.contigA.id = aScaffoldId[aHeader];

            if (bHeader !in bScaffoldId)
                bScaffoldId[bHeader] = numBScaffolds++;
            alignmentChain.contigB.id = bScaffoldId[bHeader];

            with (alignmentChain.localAlignments[0].contigB)
            {
                if (begin > end)
                {
                    alignmentChain.flags.complement = true;
                    begin = alignmentChain.contigB.length - begin;
                    end = alignmentChain.contigB.length - (end - 1);
                }
                else
                {
                    alignmentChain.localAlignments[0].contigB.begin -= 1;
                }
            }

            alignmentChain.localAlignments[0].contigA.begin -= 1;
            alignmentChain.localAlignments[0].numDiffs =
                    round((1.0 - identity / 100.0) * contigAAlignmentLength).to!diff_t;

            return alignmentChain;
        })
        .map!"a()";
}

unittest
{
    import std.string : lineSplitter;

    auto alignmentChains = readCoordsDump(testDump.lineSplitter).array;

    assert(alignmentChains == expectedAlignments);
}

/// Opions for `show-coords`
enum ShowCoordsOptions : string
{
    /// Merges overlapping alignments regardless of match dir or frame and
    /// does not display any idenitity information.
    mergeOverlaps = "-b",
    /// Switch output to btab format
    btab = "-B",
    /// Include percent coverage information in the output
    coverage = "-c",
    /// Display the alignment direction in the additional FRM columns
    /// (default for promer)
    direction = "-d",
    /// Do not print the output header
    noHeader = "-H",
    /// Set minimum percent identity to display
    minIdentity = "-I",
    /// Knockout (do not display) alignments that overlap another alignment
    /// in a different frame by more than 50% of their length, AND have a
    /// smaller percent similarity or are less than 75% of the size of the
    /// other alignment (promer only)
    knockoutOverlaps = "-k",
    /// Include the sequence length information in the output
    length = "-l",
    /// Set minimum alignment length to display
    minLength = "-L",
    /// Annotate maximal alignments between two sequences, i.e. overlaps
    /// between reference and query sequences
    maximalOverlaps = "-o",
    /// Sort output lines by query IDs and coordinates
    sortQuery = "-q",
    /// Sort output lines by reference IDs and coordinates
    sortReference = "-r",
    /// Switch output to tab-delimited format
    tabs = "-T",
}

private
{
    auto showCoords(in string deltaFile)
    {
        return executePipe(only(
            "show-coords",
            cast(string) ShowCoordsOptions.tabs,
            cast(string) ShowCoordsOptions.length,
            deltaFile,
        ));
    }

    auto executePipe(Range)(Range command, in string workdir = null)
            if (isInputRange!Range && isSomeString!(ElementType!Range))
    {
        static struct LinesPipe
        {
            static enum lineTerminator = "\n";

            private const string[] command;
            private const string workdir;
            private ProcessPipes process;
            private string currentLine;

            ~this()
            {
                if (!(process.pid is null))
                    releaseProcess();
            }

            void releaseProcess()
            {
                if (!process.stdout.isOpen)
                    return;

                process.stdout.close();

                version (Posix)
                {
                    import core.sys.posix.signal : SIGKILL;

                    process.pid.kill(SIGKILL);
                }
                else
                {
                    static assert(0, "Only intended for use on POSIX compliant OS.");
                }
                process.pid.wait();
            }

            private void ensureInitialized()
            {
                if (!(process.pid is null))
                    return;

                logJsonDiagnostic(
                    "action", "execute",
                    "type", "pipe",
                    "command", command.toJson,
                    "state", "pre",
                );
                process = pipeProcess(command, Redirect.stdout, null, Config.none, workdir);

                if (!empty)
                    popFront();
            }

            void popFront()
            {
                ensureInitialized();
                assert(!empty, "Attempting to popFront an empty LinesPipe");
                currentLine = process.stdout.readln();

                if (currentLine.empty)
                {
                    currentLine = null;
                    releaseProcess();
                }

                if (currentLine.endsWith(lineTerminator))
                    currentLine = currentLine[0 .. $ - lineTerminator.length];
            }

            @property string front()
            {
                ensureInitialized();
                assert(!empty, "Attempting to fetch the front of an empty LinesPipe");

                return currentLine;
            }

            @property bool empty()
            {
                ensureInitialized();

                if (!process.stdout.isOpen || process.stdout.eof)
                {
                    releaseProcess();

                    return true;
                }
                else
                {
                    return false;
                }
            }
        }

        auto sanitizedCommand = command.filter!"a != null".array;

        return new LinesPipe(sanitizedCommand, workdir);
    }

    unittest
    {
        import std.algorithm : equal;
        import std.range : take;

        auto cheers = executePipe(only("yes", "Cheers!"), ".");
        assert(cheers.take(5).equal([
            "Cheers!",
            "Cheers!",
            "Cheers!",
            "Cheers!",
            "Cheers!",
        ]));

        auto helloWorld = executePipe(only("echo", "Hello World!"), ".");
        assert(helloWorld.equal(["Hello World!"]));
    }
}
