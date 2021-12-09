BEGIN {
    FS = "\t";
}

function max(a, b) {
    return a < b ? b : a;
}

function min(a, b) {
    return a < b ? a : b;
}

function repeat(str, n,  h) {
    if( n == 0 ) {
        return "";
    } else if( n % 2 == 1 ) {
        h = repeat(str, int((n-1)/2));

        return h h str;
    } else {
        h = repeat(str, int(n/2));

        return h h;
    }
}

function gap(n) {
    return repeat("n", n);
}

function printLineWithGaps(lineBegin, line,  lineEnd, gapBegin, gapEnd, effectiveGapBegin, effectiveGapEnd) {
    ++nCalls;
    lineEnd = lineBegin + length(line);

    if (lineBegin == lineEnd) {
        printf "\n";
        return lineBegin;
    }

    if (lineBegin > lineEnd) {
        print "What!?" > "/dev/stderr";
        exit 2;
    }

    # find next gap
    while (gapEnds[headerId][currentGap] <= lineBegin && currentGap <= numGaps[headerId])
        ++currentGap;

    if (currentGap <= numGaps[headerId]) {
        # introduce gap
        gapBegin = gapBegins[headerId][currentGap];
        gapEnd = gapEnds[headerId][currentGap];

        printf "\rprocessing FASTA record %s ... (gap=%d[%d,%d),at=%d,line=%d,calls=%d,bases=%d)", headerId, currentGap, gapBegin, gapEnd, lineBegin, lineNo, nCalls, lineEnd - lineBegin > "/dev/stderr";
        fflush("/dev/stderr");

        effectiveGapBegin = min(max(lineBegin, gapBegin), lineEnd);
        effectiveGapEnd = min(gapEnd, lineEnd);

        # write unmodified prefix of gap
        printf substr(line, 1, effectiveGapBegin - lineBegin);

        # write gap characters (may be zero)
        printf gap(effectiveGapEnd - effectiveGapBegin);

        # process rest of line (maybe there are more gaps)
        return printLineWithGaps(effectiveGapEnd, substr(line, effectiveGapEnd - lineBegin + 1));
    } else {
        # no more gaps -- just echo
        printf "\rprocessing FASTA record %s ... (gap=ALL,at=%d,line=%d,calls=%d,bases=%d)  ", headerId, lineBegin, lineNo, nCalls, lineEnd - lineBegin > "/dev/stderr";
        fflush("/dev/stderr");

        printf line;
        printf "\n";

        return lineEnd;
    }
}

{ isBED = FILENAME ~ /\.bed$/ }

(isBED) {
    ++totalGaps;
    ++numGaps[$1];
    gapBegins[$1][numGaps[$1]] = $2;
    gapEnds[$1][numGaps[$1]] = $3;

    printf "\rcollected %d gap definitions.", totalGaps > "/dev/stderr";
    fflush("/dev/stderr");
}

(!isBED) { ++lineNo; }

(!isBED && $0 ~ /^>/) {
    if (totalGaps > 0) {
        printf "\n" > "/dev/stderr";
        fflush("/dev/stderr");
        totalGaps = 0;
    }

    if (header) {
        printf "\rprocessing FASTA record %s ... done \n", headerId > "/dev/stderr";
        fflush("/dev/stderr");
    }

    header = $0;
    match(header, /^>([^[:blank:]]+).*$/, matchGroups);
    headerId = matchGroups[1];
    numBases = 0;
    baseCount = 0;
    currentGap = 1;

    printf "processing FASTA record %s ...", headerId > "/dev/stderr";
    fflush("/dev/stderr");

    print header;
}

(!isBED && ($0 ~ /^[acgntACGNT]+$/)) {
    baseCount += length($0);
    nextNumBases = printLineWithGaps(numBases, $0);

    if (nextNumBases != numBases + length($0)) {
        print "What!?" > "/dev/stderr";
        exit 3;
    }

    numBases = nextNumBases;

    if (numBases != baseCount) {
        print "What!?" > "/dev/stderr";
        exit 4;
    }
}

END {
    if (header) {
        status = numBases == baseCount ? "done" : "FAILED";

        printf "\rprocessing FASTA record %s ... %s  \n", headerId, status > "/dev/stderr";
        fflush("/dev/stderr");
    }
}
