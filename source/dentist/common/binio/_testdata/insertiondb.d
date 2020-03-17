/**
    This package contains test data for `dentist.common.binio.insertiondb`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.binio._testdata.insertiondb;

version (unittest)
{
    import dentist.common : ReferencePoint;
    import dentist.common.binio : CompressedSequence;
    import dentist.common.alignments :
        AlignmentChain,
        AlignmentLocationSeed,
        id_t,
        PileUp,
        SeededAlignment;
    import dentist.common.insertions :
        Insertion,
        InsertionInfo;
    import dentist.common.scaffold :
        ContigNode,
        ContigPart;


    enum numInsertions = 12;
    enum numCompressedBaseQuads = 9036;
    enum numOverlaps = 22;
    enum numLocalAlignments = 22;
    enum numTracePoints = 171;
    enum numReadIds = 66;

    Insertion[] getInsertionsTestData()
    {
        with (AlignmentChain) with (Flag) with (LocalAlignment) with (ContigPart)
            return [
                Insertion(
                    ContigNode(23, end),
                    ContigNode(24, begin),
                    InsertionInfo(
                        CompressedSequence.from("aaaagggaaccaaaaacattttttgctaattaatttcccagtgattattggttgtttataatacattatggctggggttatttttcctgcttcctttctagcctaaaattgagaagaagaaatatttcctctttgggaggctctgtggactctccctgtaccatttaaatgttgcggatattcctttcccactggctctgtttaagaaacttctggaccaaaagccatcattggaagatttaaaagaactcagtccccttttgggaaagtaagtaaatgtttttcctggactttatctagtaagtcttagttgtctcggcataaatgatgttccacacagagacacagacacacaaaaggccttcagtttcaataaatgtatgctctttagatgaatttgaattttagcatgacttatagtcaagtaaaccttcaaaaatggaaaaggccaaatgtttctctaattattgaccagtacaataaacaatgcctgatgtcttcagtttatccggatggaaacagttgctttagacagtatacagtattcttttttttttttttaatatattttttattgattttttacagagaggaagggagagagatagagagtcagaaacatcgatcggctgcctcctgcacaccccccactggggacgcgcccacgaccaaggcacatgcccctgaccggaaccgaacccgggacccctcagtccgcaggccgacgctccagccactgagccaaaccggcttcggcaatatacagtattcttactgcgaatcatacctactatatttattttgaataaattgtatttaatcatatctactttatttaaataaatcagcatgttatcaactttatttttaagggttgaaaaatgcttgttctaagttaaaaggcaaactatatcctttgatgaaaatactggcagtatatgatagacaaaggattagtacagtatctctgatgtataaatatgcagagatgggcaaaagtaggttaaattaccagccatcattatcttaaagataaatcgtaagtaaggatacaaaataagaacagtaaggcagccgaggccggtttggctcagtggatggagcgtcggcctgcggactgaaaggtcctaggttcgattccggtcgagggcatgtgcctgggttgcgggcacatccccagtaggagatgtgcaggaggcggctgatcgatgtttctctctcatcgatgtttctaactctctatctctctccctttctctctgtaaaaaatcaataaaatatatttaaaaaaaaaaaattaaaaaaaaaaaacagtaaggcaataattaataaataaactccatatttcatgtacttacaactgtaaacctacttttgccaaccccggtactgtatacgtgtcatcacacatatgtatatatatggttaatcaaaaaatgaaaatatgctcaaccttactagtaatttaggaaatacaaatgaaatgtcattttttgctttatagattgacaactatcaaaagatttttcataaaatacagagttgacagccagattggggaaaaggattttcccctcatgtactgttagtggaagaacacactggctagatttatactctatgccctgtattgtcaatgtagaggaggggttgagaaatgcctggcatttttgtctggataatcctttattgtggaaggctgtcctgggtattgtaggatgttgagctgcatatgggcttacttaacccactagatgccagtagcaaccctcctctaattgcttttttgttttgttttgtaaccctcctctaattatgacaaccaaagatgtctcccgattttttcaaatatccactgggggttggggggcaaaatcaccctggttgagaatcactggtctagaggaagggactttttattgtaatttctagcaagacacagcagtgtctgtgttagcagtggttcaagcaatcaaaactaaatagaggttaactctagaaagagagaatgggatggtgcaggagaggaaatttacatggcattttatacatttctgtagtgcccccaattattccattgaacatttattacttttgaaatgctttaagaccaataaataaaaattttaaaccaaaattttcactaggaatttgcaagaagttctaaattatgaggctgatgacattggagaagcgctttgcatatatttttctgtgagtactgtcaacccatggttattggcaccaattcagtgccaacaggaatttctgatttggtgtgcaatgaaggaaactattcagtgcaaatcagctcaattgatacagcagagctgtgaaaacagcatggctttgaggcagccctggggccaggcccctcctgctagacagctctgctctgctctcctctcctcctcactcctgtgctctgctctggtctgctccttgcgggtctgctctgggctgggctactctgctttgctccactctggtctgctccatctgagcagtctgcttttgtctgtgctggtctgctgcgctccgccccagtgagacctctttggtgttcagctcggtctgggaaacatagtcctgtcctcagtggaaatggaaagtgcactccccagtcagagaggagctggcttatatagacagaagtccccaaccctggtccctgattgatctatcctcatgcaaatgaggactccatatgcttacagtttgcttggtccttaaggcattgccttgattggtccgtgtagatgttgatggggatatagttgtgcagctcctattggatggagaaagcctcaatccttttggataaaataggattccaggaactcctttagaagggctggctcacatggcagaaacacagtgtaagcaggtagttcagtgcaagcttctcccttacatgcagtttgcatgacaggtccctgtgtagaaacagctgctaggctgttgttgttttttaatttagccccaaggagcccttcttagtaggcatcttctttctcaacgttcacagtactagtgaacaacaaattatgttttagctttagtcttttaaaaaaataatcctacctaataaagagggaatatgctaattgaccctcaggccatcacaaagatggcgatgcccacacccaataaggagggaatatgctaattgactgccccgccctcgaatgtggtgatgcccacagccaataaggcgggaatatgctaattgactatcacaccctcaatgatgacagcccccacaacgaataaggagggaatatgctgattgactgtcacacactcaaagatggtgttgcccacagccaataaggagggaatatgctaattgacttacacaccctcaaaaatggaggcacccacatccactagatggcagcacccagtcccctcagcccccctggggcgcctacctccagagttccctcatcccctcagccccccagccattcagggccgcccgaggctcaggtaaccacagccggctgaggcttgcgctttcggcactggcagcagcagaagtgtgatgtgggcatcgccttccccttgattgccgggtcaccttgcgcccctttgggctccaggactgtgagaggagacaggccaggctgaggaaaccccctccagtacatgaattttcatgcaccaggcctctagttttattatatctcaaatctgggatgccttttttgaatggtactttacttctgtactgtaagactttctctctacctgctcctcattttttaccttatgtctgtgtacccacaaataaagtgcaaatacaagtgattattcaatctatgtaaataactaaatggatcaatgtatgcttaaaaggaatctacaacatatactagggtaggttttaaaatgttttatgtattcaaaaggctgaactgcataaataattatttacgtccagttatttaaattttgttagggaagaacatagagtgtcccaaaacatgtacacacactttgaataattattcattccaatggggtaaatctgaaaagaaagaaacatcaatttgagctatcagctgttaaagtgtgtatgcattttttggatatcctgtaaattggcaataatggtaatatggttgattattgagactttctttcatgtttgctacctgaatatattacagctacgatgggaccaacataatgttgatttaattccaaatgggatctctatacttgtggaccaaaccaacaagtaagttttgagacctagaatatatgcattgaatacacataaaatgtttatgtcactttactttataaatggattttagcttctcagatgatcattttagttcctcatcttttagaaagataaaacagttatagggtaag"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(23, 2189),
                                    Contig(1, 4386),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(1000, 2189),
                                            Locus(0, 1189),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 89),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.back,
                            ),
                            SeededAlignment(
                                AlignmentChain(
                                    2,
                                    Contig(24, 2598),
                                    Contig(1, 4386),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(0, 500),
                                            Locus(3886, 4386),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.front,
                            ),
                        ],
                        [],
                    ),
                ),
                Insertion(
                    ContigNode(24, end),
                    ContigNode(25, begin),
                    InsertionInfo(
                        CompressedSequence.from("tgcttgctagttggatgtccatagcttatttaatctctctgaacctgtttatcgatctgtaaaatggggacaaactatctacctcacggaagttgatgtaaggggtcagcagtacactattagttatttttggagagtcaaaagttaaacacagactttccattgtgtgtggaagaaagggatcatcacccctaaaccccaggttgttcaagggtcaacagtacatctcttgcttgcatctttttgatctattttgtcccgatccttcaacatcaacctcaggtctccacttttccattacactttctccaattcattttcccttgataactttctctctaataaaatttgatactgccctaaagtttttctggaacgttcacatttggatatatgccttgcttgttttattctcacaacagtcctgtaaagcaggcagagcagaaattatgacagccattttatatgtgaggaaaccgaaaactcagaagccaaaagactcaaccaccatcactcgagattacaactctgggcttccaattctacatcagtgtatccttccctatgcttccgtgctcaacaactattacacataaagtcagcaatccaaaatttagaaacaattgctaaagttcatggtgttttagttcttatcactctctctacaaatacaaaaagaggaaaagaattgtgtcctatactaatttgaatcactcacaatatttaatacaatgtggagtctcacgaaatagtcaataaatgcttgttgtttgacccataaaaatcctaccatttgtacaccttagccatttttattttctagaaattcattctttctcaaaaagatacagatcacttctcagtatttttgcacacccactcatctttctcttcctttccctgtcccctgggtcccaaacttctactgtaaaaagaacacttggaaacttatcagaaaggggatccataaattgaggctcaacatcattcttttataattattttctccgtgtttctaatagagaaacaagttttttattattttaaggtatgtttacattttgaggatccttgtaaaagtataatcattgcagaaatgaactgccttattagaattagcaattatctctaggattatagtaataatttctcttcttatagaaacattttagattactgcatctcaaaatctaggaaaccacacaaatatccatcagaagtgactagataagcaacatgtagtacatccatacaatggaatactactcaaccaaagggaatgagctttggttacattcaacaacatgcatgactctcaacataattatgctgactgaaagaggacagacaccaagactacttatgcttccatttctgtaaaattttaggaaattcaaagtactctatagtgacagttgttgcctgggcatgggatagggggcaccagaggcgtagcagcaggggaggagtagggttagaaagggatatgagaaaactttgggggatgatgtgtaggttcattgtctttttaaaaaatatatttttattgatttcagagaggaagggagagggagagagagatagaaacatcaatgaggagagagaatcattgatctgctgcctcctgcacgccccacactggggatcgagtttgcatcctggacatgtgccctgactgggaattgaactgtgacctcctggttcataggtcaacactcaaccactgacccatgttggccaggctattttcattatcttgatgatggtgagagtttcactggtatatatgcatatcaaaacttatcaaattgaaaaaatctatgtgctgtttactgtatctcaagtgtacttcattaaagctgtttaaaaatataaagtgcgtatgaattacttgacaatcttatgaaatacacaggcctgagattctgcttttcacaggtgttgctgaagcagctgctctgtggaccatattttgagtagcaaggctctatatgacctaattttacattatttatttacttaaaatagtgtcattttaaaaattgattttagagagggagataaggagggagagacatcaatgtgaaagagataacatgcattgggttgcctcccatacatgcctgaccgtgcctcaacttggatggcactcacaacctaggtatggcccctgactgggaatcaaatccatgaccttccagtacatggacgatgttccaactaactcagccccactgtccatggcctaattttacatgattatgtatggttctcatgttgtaatagatggtatactttttctgttttgaaaaaaaaataagaaaataaataaaagacttacactgataggaaaccacatgtaggaaccctcttcaggatacatgaacattccatattctgtcttggtcatctcttcaaacatgaagtagaagaactcggactttaccccataacccacagaacgaatttcgttaataaattcaatctgaggaaagagaacactcccagagttaaaatataatgccttgaccttttgtggaagagtttcaaacctcctgcagtcaccagaattcccagaatcccacatctcttacttcactgcaatctgttctccctcaatttgctactgagcagcagctgccatgaacctatttattgtccgattattgtatactgagtgcatggggcttttacaccgccatgaaccaagtgaaggggatgctagcaagacttagtattcaccaattatgttgacccacgtcttaaaccgaattttgatgtgctttaatcctctgggccactaaatcccctgattcttctatcttgttatctgtcttcattgcatcattcagactttacagctttggcccaatcaacaactatttaatttttcctgcctctcccctctttaggacttttttctttcttttttttaagaaagcgaagtaatttattactgttttgccttttatttagaagcgtctctgggctctctgctgggattgctccaatcattggctccaggtgtggtcttgtatcggcgccatgctgacttattgaatatgggaagcctctcgaatgattcacttgaaagagccttcaaattaattacagcatccgtacccactccttccatagaatagagttttagatctccttggaaatatctagcatacagacgagaaattggcaagccataaccaaatccagccagaggggcagctcttgaaggctccaggcttggtctaggagcagtggaatacatgtagttaaaaagacgatctatttttcgaagg"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(24, 2598),
                                    Contig(1, 3296),
                                    emptyFlags,
                                    [
                                        LocalAlignment(
                                            Locus(1600, 2598),
                                            Locus(0, 999),
                                            1,
                                            [
                                                TracePoint(1, 101),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 98),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.back,
                            ),
                            SeededAlignment(
                                AlignmentChain(
                                    1,
                                    Contig(25, 3198),
                                    Contig(1, 3296),
                                    emptyFlags,
                                    [
                                        LocalAlignment(
                                            Locus(0, 700),
                                            Locus(2596, 3296),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.front,
                            ),
                        ],
                        [
                            id_t(1),
                        ],
                    ),
                ),
                Insertion(
                    ContigNode(25, end),
                    ContigNode(26, begin),
                    InsertionInfo(
                        CompressedSequence.from("aaacttcatggcccatgcgcagaagtcggtattttgttcctagccaactcccaaacatgtcaacaggtgcatacagaatctgcaaagccacccagataagttcaccctaaatcaggagatcttctactgaactgcagataaaataaagtattgttattttcagtcattaagctttgggatgggtcattatagccatgggtaactaatatataatcccaatgttcaattaaaacccccattgcttttctctgaagggtaagcaaggaattctggaaatcataaggctgataactataattgcaaatgcaggttcttcctttgctacattaaataattacttctctcccacattagtaaaaatgtccaatgaggctttttaccgtaccaccaattctttccgaaggtcagtgacttcagcttggtctaatcgacgcagagcatcatcaactaggtgacttcgtctgacttgaagtaaaaatttgggtgaaatcctatgcgtcatttgtgactcctagaaaaagaaatattttaggaatcttcactggaatattaggaaatatttaacatttcagataaacccattttgtctaacatacaattatgcctaccactgctacatatgctcatattgtaattacagtatttaaggaaaagagcatcatatgtacacttggaagataaaaaggtggtagaatctgtgagactgctatccaggaattttgtacagaggattaatcatcgtggatgtacagctaaataattaaaagctatctctgatcacttattttagtaaacatgatccttaaaagcaacagtggtgtctatgcagagacattctcccaccagtttttaagaagagtgttagcacaataatcctcaaaaatatatgggcagtttggaagccgagagataatgtgagatgtaggaaacagactagaggaaagacttgctcagagaaaagtgttctctgtcaccctacatagtacctccattaagctcttctaaattctcaacatgccctcaccatttcagctccctgctcttgtttctgctgaactccctgggattagaaatctagtgttggaattattggccttgctttgacctttaaccttcctgtaggacttctgctatatcacaaggacaataaaagaatatttaggattcataaaccaatgccactattgagacagcaagacaggtggataataactatcctagcccaggtttccagacagtgtggttgcccatattagcactctataggcatctcaggtaagtgcatgcatatggctctgggcttgggaaagacctgtcccagatctctgctattaccatatactgagaaattatatgtgattccactagattcttctggtttttccctaaggtctttgttgtaaagaagctccacagccctagctggcttggctcagtgactagagcgtcggcctgcggactgaaaggtcccgggttcaattcccgtcaagggcacatgcctgggttgtgggctcgatccccagtagggggtgtacaggaagcagccatcatcgatgtttctatctctctctccctctcccttcttctctgaaatcaataaaaatatattaaaaaatatgtatttaaaatgaaacaaaacaagaagctccatatttgtagggagatggaaaaatgccatttctagtcttttctggtggtgggtcatctcctcctgtggccttagaatctgcagacacttcaggaatactgactgctaagctacatgtcctgtggctctgagtgacagggagtggccagaaatatcctgcactgactatagcagggaaaggtgtatgcgtcctgctgtggcccaggcttctctcgccttcccatgaattggactggaaacctgtatatcctattttccggtctctctactcccgttctccatgtactggtttcccatttctacatttgtgcatttggattgataatctagtttggagcttacaattttgactaagtccatcatctctcctctggaggacaggattgtggaacctgctttgtgacttcacccatgacaactacaactaagaactttcttgtcctcaaacatttcccaatccttcatagtatgtaacttgcctgggtccgctatgcactcaaccctcactaacagagtaaacaggttagacagaaaatttctatagaagaaaaataatgtagagggaaagatctgcttaaagggaacttaagggaaagttaggtttgggaaggaaagatttaatctattttctcccattacccagaatgtcaggaaaattgaaaaaaaaaaggtggctgggaaacagaagagatgagcaagtcaatctgatgacagacactaataacagtatctaaagggagagtccttaactatggaatgcttttagcaggatgtactttgtccttctaggttgattactgtaaaaactacaggcccggtgcatgaaaattcatgcacttgggggggggggtgtccctcagcccagcctgcgccctcagtccgggagcccttgggggatgtctgactgatggcttagtctcgctcctcagggggagagggccaaagccgcagtctggcctcc"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(25, 3198),
                                    Contig(1, 2609),
                                    emptyFlags,
                                    [
                                        LocalAlignment(
                                            Locus(2700, 3198),
                                            Locus(0, 498),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 98),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.back,
                            ),
                            SeededAlignment(
                                AlignmentChain(
                                    1,
                                    Contig(26, 1128),
                                    Contig(1, 2609),
                                    emptyFlags,
                                    [
                                        LocalAlignment(
                                            Locus(0, 500),
                                            Locus(2109, 2609),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.front,
                            ),
                        ],
                        [
                            id_t(1),
                            id_t(2),
                        ],
                    ),
                ),
                Insertion(
                    ContigNode(26, end),
                    ContigNode(27, begin),
                    InsertionInfo(
                        CompressedSequence.from("ccagcctgcgccctcagtccgggagcccttgggggatgtctgactgatggcttagtctcgctcctcagggggagagggccaaagccgcagtctggcctccctctgtgggaggtgaccggctgggcgatcaggggacagtgctggctggcgagaggccagccctgtccctccgatcactgtcccctccctctgcccgctggcacccaccttggctggcccggcagcgcccactcactggctccaccacccgccgaccagtcattcagccatttggtcaatttgcatattaggcttttattatataagatgtggtgtttatggatttacaggcatggcttattggaaaatgtcatataagctctatgattagtagggaatgaaaaggtaatactttattactgaaatgactccttttgtgcatgtgatgccccaacttttcacttagatggcccatctcctttgaattattagcaaaagtaataggtatgtttctgattcctgtaagtctttaaagtttctcttatagtctttctttatttcttttgttaatcctcatctgaggatatttttccatcgatttttagggagagtagaagagagaaggaaagacagagagaaacattgatgtgagagaaacacatcaattggttgcctcccaggctggggaggagcctgcaaccaagatacatgcccttgctgaaaccggtttggctcagtggatagagcgtcggcctgtggactgaagggtcccaggttcgattccggtcaagggcatgtgcctgggttgtgggcacatccccagtagatgtgcaggaggcagctgatcgatgtttctttctcatcgatgtttctaactctctatccctctcccttcctctctgtaaaaaatcaataaaatatatttaggaaaaaaaaaaaaaaagatacatgcccttgaccagaatcaaacctgggagtctttggtccacaggccaacgctctatccactgagccaaactggctagggcccttgtagtctttctataaacatcgtagaaaatagatacgaacctttcctgaggaatggcctagaaaaccaacttcaatgatcaaaaaatgtttggaaaatctctagcataatgtctgggcatgtggtgagtataactgagagtataaggctgatggtgacaagagagacaaaataccatgtaagtgttttgtgggaataatacaactcccagaaaaagagtggaaggaagagaataaatggaaaacagaataaattcactggaatcaaaataatatcagttggagcttacaaatgaagctccaagactttttaaagaaatcattatttattgattttagagagataggatggggggtgggggggtagagagaaacattgatttgttgttccacttattcattcactgactgctccttgtatatgccatgactaggaatgaacccacaatcttggactattgagatgatgctctaatcaactgagctacccagccaaggctcaagctccaagactcttaaacagatatgaaaatccaaagagtaacaactaagaacttagtatagtcactaaaactgggtagaggagagaaaaaaaagcatttacttgatcatttcagtattgctcatagtaggaaaacaatggatattagctaaagaaagtgaggcctcggccattatataataactagaggcctggtgcatgaattcgtgcatgggtggggtccctaggtctggctggcgatcaaggccgattgtggccatctcaccgagtcccaatcggggccgggccgatcagggccaggcctattggggctagccagccagggtaagcaacagggggagggactgtgggaggttggctgtgggaacacactgaccaccagggggcagctcctgtgctgagcttctgccgcttggtggt"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(26, 1128),
                                    Contig(1, 1930),
                                    emptyFlags,
                                    [
                                        LocalAlignment(
                                            Locus(400, 1128),
                                            Locus(0, 728),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 28),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.back,
                            ),
                            SeededAlignment(
                                AlignmentChain(
                                    1,
                                    Contig(27, 2317),
                                    Contig(1, 1930),
                                    emptyFlags,
                                    [
                                        LocalAlignment(
                                            Locus(0, 1000),
                                            Locus(930, 1930),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.front,
                            ),
                        ],
                        [
                            id_t(1),
                            id_t(2),
                            id_t(3),
                        ],
                    ),
                ),
                Insertion(
                    ContigNode(27, end),
                    ContigNode(28, begin),
                    InsertionInfo(
                        CompressedSequence.from("gaaaataccaagagggctctaaatatttccgtactttaaagggggacacattaatagctatgtaggtatgtaggtatcacattgacaaatcaacatttgaaaatgtccttactggaatgctgaatatgctaaatttgatgatcattgaaaacttaaacatctttcactcacaaaacaaatgcattttcttagaaatttctttctaagtaatttggattcttaaaaacaggttctagaaaaatatagagggaaataatttaaggtgatgaagtaaaggaatgacgttttttctcttcttttaattatgaaaactgttacctgtttgtttcttagataactgcaggaaacttcagtcccattattttcagtgaatttccatttatctttaatttgttatccaaaattaaattatggcaagctgattccattatgaaaatactggtatgtatgccttacatttttattatgaaaaaatttcaaacatatacaaaatagaaataatggtgtaatgaatcccatcacaaagttttagcagtcatcaacattttgctacacttttatttattactgttgttgtttttgtagttgttagtatttttgctgaattattttgaagcaaatcccaggcccactttcttttcacctctccatacttcagcacattatttctggttgtttcactcttatttatgctaacattcatcagtgggttcagctgttgaccacatgatatttccattgtcatgtttcccattaaccttttatttaattactagaggcctggtgcacaaaaatttgtgcactcgggggaaggggggtccctcagcccggcctgtgccctctcacagtctgggacccctcaggagataatgacctgctggcttaggcctgctcccgggtagcagagggcaggcccaatccctaggtgcagcccctggtcgggctcagagcagggccgattggggagttggggcgccgcccctgtcatgcacagagcagggcggattgggaggttgcgatgccaccctcagtcactctcagggtagggccaattggggggttggggcaccgccccctgtcacactcaaggcagggttgatggggaggttgcggcgccaccccctgtcacgcacagagcagggccgatcagggggttggggtgccacaccctgtcacactcagggcagggccaatggggaggttatggctctaccctgtcacacacagagcagggcccgtgggggggggggggggtttggggtgccgcaccctgtcacacacacagccgccgggcgatcagggggttggggcgccttcccctgtcatgaacagagcagggccgatagggaggttgtggccccaccccgtcacacacagagctgcagggtgatcagggggtttgggcactgccccctgtcacgctgatcccggtgccgggaggcctcacggctccgctgatcccggtgctgggaggcatattgcccataggatagaggcctggtgcacgtgtgggggccggctggtttgccctgaagggtgtcctggatcagggtgggggtccccactggggtgcctggccagcctgggtgaggggatggtggctgtttgcagctggcacacactcttcagggtgggggtccccactggggtgcctggccagtctaggtgaggggctgagggctgttttcaggctggcaggtgacggaagctccctacctctcctttttttctttttttttattctgggccagctttagctctgaggctccagctcttaggcctccactgctgaaagcaggtatctggtttgttacggttctataatcgaaacagtgtataactccagttctgagttcccggctccctgaaagcaggtttctggggttttgtttagcttctatatttgttacaatgtttcttaaactgcaagctcagaggtcggcagcggcaggcggggaacgttggtttcctttgtcactgaagcaagcaagcctcatgttagtttcaagctgcctggctgccggctgccatcttggctggcagttaatttgcatatctcgctgattagccaatgggaagggtagcggtcgtatgctaattaccatgtttctcttttattagataggactagtagccctgcgcatgaatccattactagtagctcactgccacttgcctccggtggctacccacccaccgcatatggtagcactcttctgctcatagctcattgctccagcttcctgtagctctccccactcttagcttgctgcccagccctcctgtagctctctgccacttgtagatcagtcgtgacattacagcatcccagccaatttgcatagtcctctattattatattagattccatccaata"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(27, 2317),
                                    Contig(1, 2377),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(1698, 2317),
                                            Locus(0, 619),
                                            2,
                                            [
                                                TracePoint(1, 2),
                                                TracePoint(1, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 17),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.back,
                            ),
                            SeededAlignment(
                                AlignmentChain(
                                    1,
                                    Contig(28, 1062),
                                    Contig(1, 2377),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(0, 600),
                                            Locus(1777, 2377),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.front,
                            ),
                        ],
                        [
                            id_t(1),
                            id_t(2),
                            id_t(3),
                            id_t(4),
                        ],
                    ),
                ),
                Insertion(
                    ContigNode(28, end),
                    ContigNode(28, post),
                    InsertionInfo(
                        CompressedSequence.from("caaaggggtcgcgacccacaggttgagaaccgctgggctaggggttcttctatctcagttcatgtgtatatcatatctgggaactagcagttcttcccgcttacccatttttctatctgttgttagtttttatttaggaggtttcaaactattatgacagctaatcctttatcaacaatatgagttacaaattttttctttaggttattaaaattgtttttactttgtggtatttttgccatgcatgattatttttaacatgtttatgaagcctaatatattagcattttcctttatggtttctgagtttgggtcagacttcaaatagactatctcactcagaaacaaagaaacaaaaacattcccatgctccaaatatagtttttattgtttcatttttatgtttaaatctttgatccatatggagtttgtcttgattaaattcttttctgggtcggtttccacttgtccccagacatttattgaacgatctgtttcattgccattaatatgcaatgccacctttatcatatctaaactctcatgtgcagctggatccatttagggaattgtgctacattgatctgcttttcatgcatttgggtcatactgttttaatggcattagcttgataatatattgtaatgtatattgataggtttattttctgagaattgtcctgtaaaggtacaataatctagactgtaacctcagtaatattaaaatactgctattgatataattttatataaaaatatatccgtatataccaggtgtgtatgaagtctgaagacatagatattattttttacagattgaaaatatccttgatgacattgtatgtgtatgtgtattcacctatgtatgttctcagactttataaataacctatagaactttcagtttgacaatcatgaaatcaatcaagtaatcaagcataaatatgaactagcaacacatttctttagaatcccgaatttgagctttcaattccaatgatacttaatcagggaaaatattttttctgtaataaatgttgtccaggatgatatccttagaaaacgtattttatgggttttgtaaacagaaatggaaaatatctatttctcacattaaagaacccagtagaaactgaaggaatactttttaggggggcagaggagtatatatttatcatggctgattttgcatccctcagagcagtgttgggcatctttacaagagtcttctctcaacacactagtccagatgcttaaatcagccatcattgctcaaatgatttattgggttgaaacacctcagaataactgcaatcttaaagctcttctagaaatgatgaaagaagtgtataaggtaagtgtgcatctaaaggaatttatataccatcttactaaatacgaatctcgtagttcccacccaattggagtcttgttctgttatgtttttcacaacagaaagagaatattttattcagctttaataattcaaaattttgcccgttttgctaactggatagactttggaggaatatattacataatacatataataaattagatacacaaacactactcaaatatgatagtggtttaacattttcattagagtgcacagcatattctagaagtaaccaaggttaggctactaaccttgttaggggacaggagaggagccagagcattttcaataaaacaagagaattcagttgatatcattccattaaaaattagctaaaataagatgcttttacaaataaattgctaaatggtaatcagagcaaagttcaccataaaattttaatctcaataggagaaaaatattacacacaagggttaatgcttttcagaagtactattatagtactatttctaaagtagtatatttctctttttctcaagctttattgagatataattgacatataacattgtgtgactttaaggtgtgcagtgtgttgatttgatacacgttatatttcacaaaatgattaacaccatagctttagctaacacctccatcacctcatgtagtaattactctttgtttattgtggtgataacacttaagatctacttgtaactttcaagactatagtacagtgcagttgactataatcaccatgctgtacattagatccacagaacttattcattttataattggaagattgtacacttggaccaacatctccccatttccctcacccccaagcccctgacaaccatcattctactttctgtttctattagttcagcttttttagctgccacatataagtgaattcataccgcatttgtctttctctgtctgccttatcccacttagcataatgtcctcaaggtccatcaatattgtcccaaaaagccaaaaccggtttggctcagtggatagagcgtcggcctgcagactgaaaggtcccaggttcgattccggtcaagggcatgtacctgggttgcgggcatatccccagtaggagatgtgcaggaggcagctgatcgatgtttctctctcatcgatgtttctaactatctctctcccttcctctctgtgaaaaatcaataaaatatatttttttttaaaaaatattgtcccaaatgtaaaaattttcttctttctcatggctaaataatatttcattatatatgtgacatgatgaaagaagtgtataagatagatgaatgaatatcacatattcttttttcattcatccatccatggatgcttaggttgtttctttcttagctaccataaacaagagcaataaacatgggggtgcagatagtttttgagagagtgattaaatttcccttggatatatacccggaagtgggattactggatcatatggtagttctgttttccattttttgagaaacgtccatgctgtttccatagtggctgcacccatttacattcccaccaacagagcacaagagttcccttttcttcacaccttcaccaacacttgttatctcttgtcattttgatggtagccatcctaacaggtgtgaagtgatatcaccttgtggtttgacttgaatagaaatatctattcagtccctctgctcattttattttttaaatatatttttattgatttcagagaggaagggagagggagagatagaaacatcaatgataaaagagaatcattgattgctgcctcctgcacgccccctactggggatcgagcccacaactagagtatgtgtccttgtctggaattgaacccaagtgccttcagtccgcaggccgatgcttctgtccactgagccaaaccggtcagggcaaagggctgaatttttatggtgcaaattttattctcaagatctgggaactgtttccagaaaacctgtgccccaagaatttctgtactatctccactgccaattggagaaattagagatttaaaaaatggaagttttatagagttctacttcatcattaataaatgctttgggtttatagagctctattctgattacctatgtaatagcactctctccctctttctcatatatgagagaagattacatcatctgttaatagaaaattatgttcacaggtcatattctctttaaggtaaacaaagctacttgtcaactaccaaagaacactttcatcataaatgaactctccactatgttgaactttattgaagtaagaagaagaatgtcctatagagatgacaacctggtaagagtaacattttcacttctgagtgaaaataccaagagggctctaaatatttccgtactttaaagggggacacattaatagctatgtaggtatgtaggtatcacattgacaaatcaacatttgaaaatgtccttactggaatgctgaatatgctaaatttgatgatcattgaaaacttaaacatctttcactcacaaaacaaatgcattttcttagaaatttctttctaagtaatttggattcttaaaaacaggttctagaaaaatatagagggaaataatttaaggtgatgaagtaaaggaatgacgttttttctcttcttttaattatgaaaactgttacctgtttgtttcttagataactgcaggaaacttcagtcccattattttcagtgaatttccatttatctttaatttgttatccaaaattaaattatggcaagctgattccattatgaaaatactggtatgtatgccttacatttttattatgaaaaaatttcaaacatatacaaaatagaaataa"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(28, 1062),
                                    Contig(1, 4224),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(99, 1062),
                                            Locus(0, 963),
                                            0,
                                            [
                                                TracePoint(0, 1),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 62),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.back,
                            ),
                        ],
                        [
                            id_t(1),
                            id_t(2),
                            id_t(3),
                            id_t(4),
                            id_t(5),
                        ],
                    ),
                ),
                Insertion(
                    ContigNode(30, pre),
                    ContigNode(30, begin),
                    InsertionInfo(
                        CompressedSequence.from("cactgccttacaggaacaattatgattgcagtttttgttgggagttgagtgaggggaacctgatagagagaggggacagtgttcattttgtggttgttgtatttcagatgctgtattaagcgccggctcttacttcctctttgaacacttaaacatacagaagatgacaatgctttaattttctttctttcttgattgattataattttatttaaaacatatttttattgatatttagagagaggaagggagagagaaaaacatcaatcagtcgcctcccacacatgcaccaaaaggggatcaaagtcgcaacttcaggtatgtgccctgaccaggaatggaatctgccacttttggttaacaggttgacgctccaaccgaccaaccagggtaacaatatttttattttcaatatgagaaaaccattacttgagaagtaaataacctctatcatgtttgggaagcaggataatcctggactaagagactgaaaaggctggtctctatttcttgagatgaatctcaaatctagttgaggtcaatctgctactaattcatcccttatcctgagaaaactaattgtaagaacactgatgtagtaataataataataataataataatagttaataataataatatagtactttacagtttattgttgcatatatcatcatacttaacccttacgtcatgtagaggcagacactatagttaagagtatctgagatctattagatagtataataatgctagctgttgtcacaaaacaacttctaaacttcagtggcttaatacaatgaagatttattttccacatcatagttccatgaggattgaaaatatgagtcacctgtagagggaaattttttccccattctatccatcttagtttcattgggtgagaccctacaaaatagactgagaaaagacaggttaataggtgaaaaacaaaccaaagtttattcacgtgtgcatccttcgtatatatactcagggataagaaatttgaaggggtggttagaacttgggcttatatagcatcttagtaaagtacaatcatcttgtatagaagtgaccagacaaaggaaaaggactttgagcttctagaggtgccaaattgtgggaaggcaaatatatggggaaatgaatggtaggaaagggctagttagtaagatttgtttgtgtggactttctcagtgccagctcatctctgctgtataggtggttatcctcttcctggtacaggaagttgggaataagaggacaccttcacaaaaggaatttatattcagttttcaggtaaatagaggggaagccagaaagcttgttctagtctgttttttcttacttgccttcagctcaaaatgatcctcaagcctggccagcatggctcagtggttgagtgtcgatctatgaatcaggaggccacagttcaattcccggccagggcacatgcccaggttacgggctcgatctccagtggcggggtgtgcaggatgcagccgatcaatgattctctctcatcattgatgttcctatctctcgctccctctcccttcctccctgaaatgaattaaatacacacacacacacacacacacatacacacatacacacacacacacacatttaaaatgaacctcatgcaggagtggcatcttttggggtggcatagtctgctacccttccaacctctgcacagttgaattcagcagtggcccaggctccttctttctgggttctggggccatgaaaccctccagtagatcttctgtgtttggttggcagattaggaaagagagtaaaggtcatgtgagaggtttctgggggccatgccttgaagtgtcgtacatcacctgcactcacatggcattagctagaacttggtcacatggccccactaactgcgggggcggctagaaaataccatcttccgtgaacccaggaggaagctgcaacaggactggagaggggaacacatggcattttctctgctgaggcagattttagcccagacaggtctgggtatggatccctgctctgtcactttctaattgaaccacctccacaagtccctcagcacctacgtgcacatgggatccacttttatttcctttatgcaccccaatataaaatcggaaggatggctcgcacgtactacttcaggggcatgtgggatgcaaggaatgatgaacaacaggtttcaggaaaactgaatgatgattctagttgggaagtggaggaatgaaacccctggtcccaaaacttaccccgggcctttcttccttggatcccgtttcaactgcagatagctacaggtctcagggatcatctgctcggagctcttccatgccgttctccacaccaagaagcgttgtcagttttccttctcctcccacaatgtcctgtgatgcacgatcctaggaattcgggaagtctggtggttccctttgcacaggctgtttgtaaaatgagtgaccgatctttacgagtcctgcgtaagtttgattcctttcatgttactttatttgaactctcaatagaatcactgtctgaagtagaagttgtagcaaattttgaatagatttcaaaatacacaaatccttgagtacatgatatgtgaattattggaaagataaacctaaagcaattctaatcattactgttataccgtgcaaggcatgtacaatttcttagacttggacatggctttgacatcagagtttagtaaccgtagaggatccagtcacgataaatccactttgaaacattttgcttgtgctgttttaaaactttagacattttaaccagaggggatgaaaaagaaatgttaagttttgtgttcttgtttattttattcactacaataggatgtctttttacaattttaagcaagtgttcaaggtccaaaacttgaaaggaacaaaggatatatacagtgaaatataagccacctagttcccctctcagaaagcaattactgaatatagtttcctgtgtaaactttgaaagatattacacataaaaacatatagaaatgcattttttagtggcagcataccattcgcatgctttggcaacttgcttttccccatttaacagtatattatggagagcattctgctgttgaatcggtagaactgccttattctatgtaatgtctgaatagattgtattgggaggctatactatgaatacttgtttctaatttttgctatcataacagtactaacatgcatttctttgtacatttgtaaatgtacacagtcagggaaaattccagaaggcaaaattactgggtctagttgtacatggattttaaattttgagagattttgccagatgaccctccacagaagtctgggccaccaggtagcccaacttcaaggggccgccatccacactgtatctaagtaaaatgcacccctggagctgggcaaagtagccgccaggacagaggatgtatgtttacccttcacagcagccaatgaattcatgtttctacagtgttatcaaacttttggtctttcccattctgaaaaaaatactcttttaatttaattttctcattaagaatgaggttgtgctggttttcttctttgtttgttttattaatcctcacctgacgatatttttccattgacttttagacataatggaagggagggtgagagacagaaagagagaaacatcaatgtaatagagacacatcgattggttgcctcccgcaagcgccccaactggggactggcaatcgaatctgcaacagaaatatgtgcccttggaccggaaccaaacccaggacccttcagtccttaggcagacgctcctatccactgagccaaacaggctagggcagtggttctcaacctttctaatgc"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(30, 26781),
                                    Contig(1, 3884),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(0, 1600),
                                            Locus(2254, 3884),
                                            38,
                                            [
                                                TracePoint(34, 126),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(1, 101),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(1, 101),
                                                TracePoint(1, 101),
                                                TracePoint(1, 101),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.front,
                            ),
                        ],
                        [
                            id_t(1),
                            id_t(2),
                            id_t(3),
                            id_t(4),
                            id_t(5),
                            id_t(6),
                        ],
                    ),
                ),
                Insertion(
                    ContigNode(30, end),
                    ContigNode(31, begin),
                    InsertionInfo(
                        CompressedSequence.from("actttcgttatgttgatatatataaacttttactctggttattccacaaagtttctattactgggagctcccatgtgtacatataaataagccatgtcttttccactgttaatctgtctattgtcaattaattgcaggtcctcaaatcactgcacccaaattggaaaaagaaatgtttttctttcaacataacctttttgttataatcactgggtcctaatactttacttagtaaaatggagggggggaaatcatatcctcacaaaaatgttttctaggacaaatgatatgacaagcagtctgaattttgtgttcatgaattacacactttggaaatggtcaagagcatgttttctaatactaggaagttaccttgtcaccccctgtggtcattaaatcatgtgagaaattgctgtggcttgaggaatgatcattgaacttgaaacggaaagacgaaggtctagttctacaacttgaagcatatctgacttctttgagttctaggttcaaccatgtgaaattgcggatagttacgttttgatcctgtagaggccattcaacccaatggttttgtttttctcacagggattaaaatagaatcctataaattctaaagccctgtgctatatgtgtatgaagaaagttattgggggagaaccaagatggcggcataggttaacgccggagtttgctgccttgaacaactacttcaaaagtgaaactaaaagacggaacggacatcacccagagccacaggaacactggctgagtggaagtcctacaactaggaggaaagagaaatgcatacagacactcagaggaggcgcagtgctgaagtcaaattctgaggtgcggagtgtgcggagtgggctggcggcggagggcgcggttggcattttcaatcgggagggagtcgcagactctgagctccagatacgggcgagtctttagggacccagactcaaacgggagaagcgggactgtctggcttgggtcagagcgagtgcagctttctctccgaggtttgcagcagttgctgggactcagagaggcagagcccctggggacaggactgagagccaccataactgctctcttcagcccaccctgttgatcctgtgcgacccgccccgcccaagccctgcacagaggcatttgccggatagcctcaggcaaaggctagattagcacctccctagaggacagaagttctctcactgcagacacagctgattctcacagccacttggcctggaggtcaaaccctccctggtattagctacaacaatcaaggcttaactacaagactgcgaacaaagaccactagggggtgcaccaagaaagcataacaaaatgcggagacaaagaaacaggacaaaattgtcaatggaagatatagagttcagaaccacacttttaaggtctctcaagaacgtttagaagccgccgataaacttaatgagatctacaagaaaactaatgagaccctcgatgttatgatggggaaccaactagaaattaagcatacacggactgaaataatgaatattatacagactcccgacagcagaccagaggagcgcaagaatcaagtcaaagatttgaaatgtgaggaagcaaaaaacacccaaccggaggggcggaatgaaaaaagaatccaaaaatgcgaggatagtgtaaggagcctctgggacagcttcaagcgtaccaacatcagaattataggggtgccagaagatgagagagagcaagatattgaaaacctatttgaagaaataatgacagaaaacttcccccacctggtgaaagaaatggacttacaggtccaagaagcgcggagaaccccaaacaaaaggaatccaaagaggaccacaccaagacacatcagaattaaaatgccaagagcaaaagacaaagagagaatcttaaaagcagcaagagaaagaaactcagttacctacaagggaatacccatatgactgtcagctgatttctcaacagaaactttgcaggccagaagggaatggcaagaaatattcaaagtgatgaataccaagaacctacaaccaagattactttacccagcaaagctatcattcagaattgaaggtcagataaagagcttcacagataaggaaaagctaaaggagttcatcaccaccaaaccaggattatatgaaatgctgaaaggtatcctttaagaagaggaagaggaagaaaaaggtaaagatacaaattatgaacaacaaatatgcatctatcaacaagtgaatctaagaatcaagtgaataaataatctgatgaacagaatgaactggtgattataatagaatcagggacatagaaagggaatggactgactattcttgggggggaaaggggtgtgggagatgcgggaagagactggacaaaaatcgtgcacctatggataaaaaaaaaaaaaaagttatttttcttttacatcaataagtgtagaatgaaagttggtacattaaataatgcattaaaagtccagtaaaagtcacttttgcattatatactgttttgaattctgccctgcacctttgtacatgctgatgtcatgttaaaaacagttcaaagctttgaaaaaaagatactacagaaatacaagaacattaaattgtttttctaaaccaaatgtattacctgggatgttttaaatatctgtagcccacttagtgattctaaagaatgattaatagactttgtatattttttatgtttcatttttaatacaggaagataaaaactctggatggtataaaaataatacaagtttcctgtggagaccaccattccctggcattatcagaaggtaaaaaaggtgttttcaatctcctgttggcttttaatgcattggggttcactttaatgcatttagaagtttttaactctccagatttcagaagagatgttgttttcatatgtttttccaatt"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(30, 26781),
                                    Contig(1, 2959),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(26200, 26781),
                                            Locus(0, 581),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 81),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.back,
                            ),
                            SeededAlignment(
                                AlignmentChain(
                                    2,
                                    Contig(31, 7099),
                                    Contig(1, 2959),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(0, 700),
                                            Locus(2259, 2959),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.front,
                            ),
                        ],
                        [
                            id_t(1),
                            id_t(2),
                            id_t(3),
                            id_t(4),
                            id_t(5),
                            id_t(6),
                            id_t(7),
                        ],
                    ),
                ),
                Insertion(
                    ContigNode(31, end),
                    ContigNode(32, begin),
                    InsertionInfo(
                        CompressedSequence.from("gacacaggtgcaggtgcgggcgcggcgccggcaggtggggcccgggagccctgggtcctctgggcaggagcctgcaaagggcgccccgcgccggctgggtgcaacttcgtcagtttcgctttcccagacgcgcctcgggttagtttcggtttccagttcctcttaaaagaaaaactgaaagcgaaacccgtcagtggaactgaactggaccgctcttccctgaagggtgctgtagcccctcccacttttctttggagaccagtccttgtctgcgcaaagctctccagaacggtgaccttccaaatgtggttttggggggttgatacctgtgagtccctcagcaataccctagggcgtcttttttcttttcttttcttttttttaagtattttttattgttaatagtgttacagataccacccgtcccgccccacccccctttaccccttactcttattcagtcacctttgaactttttttttcctgttactttttttttttattaaatctttattgttccgattattacagttgttcctctttttccccccatagctcccctctacccagttcccacctcaccctctgcctttaccccccacccactgtcctcatccataggtgtacgatttttttccagtctcttcccgcaccccccacacccctttccccccgagaattgtcagtccactccctttctatgcccctgattctattatattcaccagtttactctgttcatcagattttttattcacttgatttttggattcacttgttgatagatgtgtatttgttgttcataatttttatctttacccttttcttcttcctcttcttaaagaatacctttcagcatttcatataatactggtttggtggtgatgaactcctttagctctttcttatctgtgaagctctttatctgaccttcaattctgaatgatagctttactggataaagtaatcttggttgtagggtcttgctattcatcactttgaatatttcttgccactcccttctggcctgcatagtttctgttgagaaatcagctgacaatcatatgggtactcccttgtaggtaactaactgtttttctcttgctgcttttaagattctctctttgtcttttgctcttggcattttaattatgatgtgtcttggtgtggccctctttggattccttttgttttgggttctctgtgcttcctggacttataagtctatttctttcaccaggtaggggaagtttcctgtcattgtttcttcaaataggttttcaatatcttgctctctctcatcttctggcacccctataattcggatgttggtgcgcttgaagctgtcccagaggctccttacactatcttcatatttttggattctttttgctttttgcttttctggttgggtgttttttgcttcttcatatttcaaatcttggacttgattcttgcgatcctcttgtctgctgttggaactctgtatattattctttatttcagtcagtgtatgcttaatttctagttggtcttttttcatatcctcgaggttctcactaaatttatcagtggtttctagaaaattcttgaaaaatcttataaccatgattttgaactctatgtccagtagtttgctttcctccatttctgtcatttgtgacctgtttctttgtctccgcattttttacgcttccctgtgttggtagagtggctttgtgtgctgggtgtcctatagggcccagtggctcagactccccaattacctgagatggacactcttggtgcacccctttgtgggctgtgtgcacagtcttgttgtagttaagccttgattgttgttgttatcactggtaggaattgacctccaggccaattgtctgtgagaatc"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(31, 7099),
                                    Contig(1, 1891),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(5800, 7099),
                                            Locus(0, 1299),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 99),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.back,
                            ),
                            SeededAlignment(
                                AlignmentChain(
                                    1,
                                    Contig(32, 18877),
                                    Contig(1, 1891),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(0, 500),
                                            Locus(1390, 1891),
                                            1,
                                            [
                                                TracePoint(1, 101),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.front,
                            ),
                        ],
                        [
                            id_t(1),
                            id_t(2),
                            id_t(3),
                            id_t(4),
                            id_t(5),
                            id_t(6),
                            id_t(7),
                            id_t(8),
                        ],
                    ),
                ),
                Insertion(
                    ContigNode(32, end),
                    ContigNode(33, begin),
                    InsertionInfo(
                        CompressedSequence.from("ccctggagccaacctcccacagtccctccccagccaccagtagcacctggagcagcaccccatggttcaaagggtgtcttcagagttatcagggtccctccggcaggtggggtccctcggcctggcctgcagggatcgagctgaaactggctctccaacatcccctgaggggtctcagagcgcaagagggcactctgcaaagttgctgtcgtacaggatatggaacacacatgtgtgcagtaaaccagattcggggccaacaatgccccaggaacaacataacccatggctgaaggtggaattgggtggccctgtgaggaatcaggctccctcctctctgattccgggatgtatcaccccagaacctccagtgccaagtcatcacagctcagcagctcctgtgttgagcaactgcccctggtggtcagtgcacatcataggtaccagcaggtcagacggttgcatggttgcttagccttatatatatagattgtttctggtgttccatttcaaatgtagttatttcttccacacggtacataacttcttatcttaaaaaattatataagtaatacattccattgtaaaaagattttttttttaatatattttttattaagatattacatatgtgtccttgtccccacgttaccctctaccccccccccccccgctcatgccctcacccccccgttgtccatgtccattgattaggcctatatgcttgcatataagtcctttggttgatctttcccccttacccccaccctcccctaccttccctccaagtcccgacagtcatcagcttatgttgttcattatattccacaaatgagtgagatcatgtggtatttatcagctctccagttctatccatgctgtggcaaatggtaagagttcctttttcttcttcctcttcttaaagaatacctttcagcatttcatataatgctggtttggtggtgatgaactcctttagctttttcttatccgtgaagctctttatctgaccttctattctgaatgatagctttgctggataaagtaatcttggttgcaggttcatgctattcatcactttgaatatttcttgccactctcttctggcctgcatggtttctgttgagaaatcagctgacagtcgtatgggtacacccttgtaggtaactgactgtctttctcttgctgcttttaagattctctctttatcttttgctcttggcattttaattatggtgtgtcttggtgtggtcctctttggattccttttgtttggggttctctgtgcttcctggacttgtaagttcatttctttcaccaggtatgggaagttttctgtcattatttcttcaaaaaggttttcaatatcttgccctctctctccttctggtacccctataattctgatgttggtacgcttgaagttgtcccagaggttccttacactatcttcatatttttggaatctcttttcttttttctttttcggttgggtattttttgtttctttgaatttcaaatctttgacttgattcttgggatcctcttgtctgctgctggatctctgtaaattattctttatttcagtcagtgtatgcttaatttctagttggtcctttttcatgtcctccgtggtctcactatacttattgagggattcattaaatttatcggcggtttccataaaattcttgaaaaaccttataaacgtggccttgaactctatatccaatcgtttgctttcctccatttctgccatttgtgacctgtttctttgtctccgcattttggcagcttccctgtgttgatagagtggccctgcgcgccaggtgtcctgtagggcccagtagctcagcctccccagttacctggggtggacactcttggtgcactcttggtgccttggttgttgtagaatcactgggaggaattgacctccaggccaattggctgtgagaatcagctgtgtctaaagtgggagaacttctgtgctgaagacacccttctggggcaagacttgcttcagtggggctttggtgctcactgagtctgctccctgagtgtgtcccttatggatctgaggagttgcaatctggatggtccagtggctcctggatctaaggaggtgctaattcagcctctgcctgaggccacacagcaggagccacgaagaggctcagctgtgaagcaatgcaagccgcttcggggccttgggccttctcttggaagtttcgggtctgcctggctgggctgcagttaggtacattcatatgcaaaagcctctgccgcagcctgggtggggcggggtctcagggaatcaaagggcggatcagggagctatggccgattctccgcccggagtcgcgggtctcagtgtcccggtaacggctgcaagcacctctgagggaaagccgcattcaagttcgcccgctgccggacagaccagcttctccccgtatgagacctgggtccccagagacttcccggaaccggagttcagagcagtcgggagtttgtgattcaaaaaggcagctgcgtcctcaggcgccaccccctttcctcgcgcgcgagccgccgtacctctgcacttcacctccgcagcgcctctggctctcagcgtgcttctctttccttctagttgtagaatttacactcagccagcctttctgttgttctggatgatgtttgctcagtcctttgcggtaattttcaattgttgtgggaggcgacaatttcccggtgtttacctatgccgccatcttagtttctccgtaaaaagatttttttaaagtcatgcaataacaagtatctagaggcaaaagtatgttacctctcctttccctggtcccattccattccccacagtcccttaagcaacagcttgattaaatccacaagttttcccaaaagtaaaacaattcttggggttatttattctctctcagaactgtaatctaaataaaaatatagtttctctaaaatattttagtgttttaaaatttgtttttattgatttgagagagagacagagagatggagagagaaacatcagtttgttgttccacttatttatgcatcattggttaattcttgaattcttgaaaagttccctgatggggcatcaagcccacaaccttaagatgatgctctaaccaactgagctaccaggcctgagtttagtggggtgtgttgttgtt"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(32, 18877),
                                    Contig(1, 3197),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(18399, 18877),
                                            Locus(0, 478),
                                            0,
                                            [
                                                TracePoint(0, 1),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 77),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.back,
                            ),
                            SeededAlignment(
                                AlignmentChain(
                                    1,
                                    Contig(33, 1084),
                                    Contig(1, 3197),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(0, 600),
                                            Locus(2597, 3197),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.front,
                            ),
                        ],
                        [
                            id_t(1),
                            id_t(2),
                            id_t(3),
                            id_t(4),
                            id_t(5),
                            id_t(6),
                            id_t(7),
                            id_t(8),
                            id_t(9),
                        ],
                    ),
                ),
                Insertion(
                    ContigNode(33, end),
                    ContigNode(34, begin),
                    InsertionInfo(
                        CompressedSequence.from("aatacttcagctttgataatgttaggtcataatcttttgataacaaatctactatcattacaaagaccccatttactattgagatttctcaatagcatgataattgcacccactttcaattttaacttatgacatggcattctcgaaggagtaatactattaagaaattcgatgggaaaattttccttttcagcatcatctgtggaatcaatggaatcatcactcaaatatgtgtgaaaatctccatcgagtataaccaaattttcttcatttaatttgaacatgctcattttttggacaaagaattgcacgtttagatatatttttaatattttctatagaagtggttctcaaccttggctgcacattagaatcacctgggaatctttttaaaatcctgattcctgggcctcatcctccagaaattgtttctttgttatggagtggggccacaacattagtaacaaagaaacagaactaaaactttattgtggggagccgctacagtgttttggacaggttcaagccttggactgatgagagggcacagtcctctcattgccacgtgcaagtcccagcttggagggcagggccctcccagcacaattatagccaacagctgtagttgtaagcttgaacctatggtcagaacacgtggccccatgtgcctgagtaatgtaggcttgatagagttcaggtgcatgtagttgagtaggattgaaacccaggagaatgttgtaaacatcttggtcatgcccttactttgcctcgtagcatctgctataaaataaaggcatggctgtgggcgtggtcactgtctctcagagaagcagcgtcccaccaagacccagctttcattctcttgtctgtcttttctaagcctttcagccacccccactcaggttcacccctggccacactttagaaccgtggtcggcaaactcattagtcaacagagccaaatatcaacagtacaacgattgaaatttcttttgagagccaaattttttaaacttaaactatataggtaggtacattgttattaacttaattagagtattcctaaggcttaggaagagccacactcaaggggccaaagagccgcatttggctcacaagccacagtttgccaaccactgctttagaaaaatatgcctccacgctcccagggcagattcctgcctcaaaccctgcccttacaggaaacagctcggggaaacaacttgggggagggcgcagccgttgccaagacaacccctgcaggactgacttccttctgcttggtcgaggcagttgtggtcatgaggaacaccacccatggtttttatataagtagactagtggcccagtgcacaaattcatgcacattgaaaggaaattaattagaattcctggggcagcactgcagctggaagggtgtctttggagtgaagcgggggaggggaagaaaactcacgctggggtcttctggctcacagcactggctccggtgcagggactaatgcaggcagtcctgagtggtggctggcagggcgggactgggcaagaagggccagacacg"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(33, 1084),
                                    Contig(1, 1561),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(600, 1084),
                                            Locus(0, 484),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 84),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.back,
                            ),
                            SeededAlignment(
                                AlignmentChain(
                                    1,
                                    Contig(34, 3004),
                                    Contig(1, 1561),
                                    Flags(complement),
                                    [
                                        LocalAlignment(
                                            Locus(0, 500),
                                            Locus(1061, 1561),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.front,
                            ),
                        ],
                        [
                            id_t(1),
                            id_t(2),
                            id_t(3),
                            id_t(4),
                            id_t(5),
                            id_t(6),
                            id_t(7),
                            id_t(8),
                            id_t(9),
                            id_t(10),
                        ],
                    ),
                ),
                Insertion(
                    ContigNode(34, end),
                    ContigNode(35, begin),
                    InsertionInfo(
                        CompressedSequence.from("tttgtctttctctgagggcttatttcacttagcgtaatattctactgaccaacaaaataggtccagaggcacagaggcaaggaaaagactgacaaatctcagaggggcagggagggggattggaagagattaattaaagaatttatatacttatatgcatagcccatggacacagacaatagtgtggtgaagacgtgaagtggggtaggggctggatagaggggggcaaagggtggcagaaatggggacatctgtaatactgtcaacaattaaaaagaaaaatattgcctggctggcatggctcagtggttgagcattgacctatgagccagaaggtcacagtttgattctggtcagggcatatgctcagatttcagactcaatccccagtgggaggcgtgcaggaggcagctgatcaatgattctctctcatcattgatgtttctctctctccctctctgaaatatatatcctatctaataatagacaaacatggtaattaactgtacctccgacacgcttcccatcggctaatcagggtgatatgcaaattaactgccaaccaagatggcggccggcagccaggcagctgaagcgaacaggaggcttgcttgctccagtgaaggaggaagccaaggttccccgcctgcactgccggcctctgagcttgcactctaagcaactatgttgcaattatagaagttaaacaaaacccagaaacctgctttcggctagctgggcttcagccaacaggatcgcaacagtgtttcaattatagaactaacacagatacctgctttcagcagccgaggtctcagagctgcagctgggcctcagagctaaagccagctctcagctccagtgacagccatagaaggtaaataaatcccagaataaaaaaagaaaaaagaaaaaaaggagaagttgggagtttcagtcgctcgccagcctgaaaacggccctcagcccctcacccagactggccaggcaccccaggggggacccccaccctgaggagggtgtgaccagctgcaaacagccatcatcccctcatccaggctggccaggcaccccagtggggacccccaccctgaagcaacaagatggcggctaatttgcatactgaaggcagcggcagagggcaaggctgtccatggtccgggacccagatctgtgacctggcgggcggagggcagtgcttgaggctgtccgcggtctgggacccggatccatgacctggcgggtgggggcagcgcttgaggctgtccacgggaccaagacacggatctgtgaccctgccaggtgggggggcagcgcttgaggctgtccgcggtcccaggacccggatccatgaccctggcgagcagggggcagcgcttgaggctgtccacgggaccaagacacggatctgtgaccctgccaggcggggggcagtgcttgaggctgtccgtggtcctgggacccggatctgtgacccttcctggtggggcggcagctctcacaggggcggatcggccagggttgggcagggcaggatgcctggcttccacccagccccagtcactctgggcaggtgagcagcacagagagcctctggacaggggaggcaggctggcggagggcggcctgtaatccccacatggcggtggtgacagctgtgagcgcgccaagcggagcctgcaggctgagctgaagtagcgcctgcaggccatcagcacccaggagtggctgcgcaagatccgcctcctggcccagaaggtgcaggatcgcagggacagccagactgacgtcctcctggctgcaccacgggggtttgcggatctgcaaagccagaggcctccggtgtgcacatccccgcccctccctctcccccccccccccccccccgcgtgcgaacatcccccagcacactgtcaccctcccgcgcctgcaacagttgcaaaccaaggtatgcacaagtcccccctgcctctccactgcacttccatcaacccagcatgcctaacacccccacaggatgtgcacacatcccctctgatgcactcacgtcctttggagagtgtctgggtgttctgggtgctgagggtgggcgcctgtgagatgacagtgagagggcggagtggtgatgccctgttcccatggaggccagtgcagctatgagatcacctctagggggttaatgcaggtgctgaggcaggagcaggtgcagacagacacacctggtgggcgcgggcggccctcgctgaggtgccttcagtcaacgttcaagatcaagaacccagaggtggagaggaattccgtcctctcagtgaaacagtgagatagtcggggtgactgtagaggaggaatgtgcactttgggtgccagcacccatgagcgaaggctgccctgactgagcctcgcttgcttggggcctagcgctctcctgccaggcagcgggtgacaggacgtgttttccaaattaatgagagaaaacgttgattctcctgctgcccacgaaggtagaaaatgaagtggggagacatgtggggagtgagggaggttctctgtctgggagttccaaatctgctgttggtttctttcaccgctgactagagatatacaggatatccccccagaatgtatagacactttgaatacctactaattccaatggttttaagcataagaaagaaacaacaatggagctgttatctgttaaaagtgtgggcacaaacaggtagttggacattccctgagtggtctcagattggagaggttgaaggccagactgaggggcccctccccccaccccaccccagtgcacgaatttcgtgtaccgggctcctaattatatatatatatatatatatatattatatatatatatatatatatatattaagtaaataaataaaaaggggagggggaaaaaaagaagaacagtgtagcagcaatagccccctacagcagaggagaaagagattagtagccacaatgagtgtaaactttgggaagagagagagagctgtttaaaggggatttatagcagtgatggcgaacctatgacatgcgtgtcagaggtgacacgtgaactcatttttttggttgatttttctttgttaaatggcatttaaatatataaaataaatatcaaaaatataagtctttgttttactctggttgcaaatatcaaaaaatttctatatgtgacatggcaccagagttaagttagggtttttcaaaatactgacatgccgagctcaaaaggttcgccatcactgatttataggaaagtttgcctggtaatgaatgtgctgcacccaccaaccccatcctcagcacagcactccacccttattcctaagcactctaacaagcagtactagtatgcagacagactgccaatttgagatctgactcctccatttactgactagatggctttaggcaagtcagtgtacctctctgttcattgatttttctcatctgtaaaaatgctgaaatggtcatgcgtactttacaggagttgtgaggattaactgaattaattctgctccataattgaaacaaagcctgactatgttctcagagctacacagctatctatgtgtgagttagctattatttttctgggccttcataatccaattgatgctgtttgaatgttacttgagcagattgtggctcaaaaataatcttgtagaatttgagatacttctttggttgtcatcattatttgccctgtcacctaatcactcctttctcccagcgttaagcaatcactaagttctgttgcttcatgctcctcagtatttctcaaatccattca"),
                        0,
                        [
                            SeededAlignment(
                                AlignmentChain(
                                    0,
                                    Contig(34, 3004),
                                    Contig(1, 3811),
                                    emptyFlags,
                                    [
                                        LocalAlignment(
                                            Locus(2200, 3004),
                                            Locus(0, 804),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 4),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.back,
                            ),
                            SeededAlignment(
                                AlignmentChain(
                                    1,
                                    Contig(35, 10490),
                                    Contig(1, 3811),
                                    emptyFlags,
                                    [
                                        LocalAlignment(
                                            Locus(0, 600),
                                            Locus(3211, 3811),
                                            0,
                                            [
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                                TracePoint(0, 100),
                                            ],
                                        ),
                                    ],
                                    100,
                                ),
                                AlignmentLocationSeed.front,
                            ),
                        ],
                        [
                            id_t(1),
                            id_t(2),
                            id_t(3),
                            id_t(4),
                            id_t(5),
                            id_t(6),
                            id_t(7),
                            id_t(8),
                            id_t(9),
                            id_t(10),
                            id_t(11),
                        ],
                    ),
                ),
            ];
    }
}
