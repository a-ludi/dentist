/**
    This package contains methods to handle the proprietary binary data
    containers used to store information between stages of the algorithm.
    Currently, there are two containers:

    $(UL
        $(LI `dentist.common.binio.pileupdb` for pile ups,
            i.e. candidate sets of reads for gap closing)
        $(LI `dentist.common.binio.insertiondb` for insertions,
            i.e. consensus sequence and splicing information generated
            from pile ups)
    )

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.binio;


public import dentist.common.binio.common;
public import dentist.common.binio.insertiondb;
public import dentist.common.binio.pileupdb;
