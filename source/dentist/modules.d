/**
    This module lists all modules in DENTIST. This is used to derive a list of
    all external dependencies in DENTIST.

    DO NOT EDIT! This file is generated by `update-modules.sh`.

    See_also: `dentist.common.external.externalDependencies`
    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.modules;


import std.meta : AliasSeq;
static import dentist;
static import dentist.commandline;
static import dentist.commands;
static import dentist.commands.bed2mask;
static import dentist.commands.buildPartialAssembly;
static import dentist.commands.chainLocalAlignments;
static import dentist.commands.checkResults;
static import dentist.commands.checkScaffolding;
static import dentist.commands.collectPileUps;
static import dentist.commands.collectPileUps.filter;
static import dentist.commands.collectPileUps.pileups;
static import dentist.commands.filterMask;
static import dentist.commands.findClosableGaps;
static import dentist.commands.generateDazzlerOptions;
static import dentist.commands.maskRepetitiveRegions;
static import dentist.commands.mergeInsertions;
static import dentist.commands.mergeMasks;
static import dentist.commands.output;
static import dentist.commands.processPileUps;
static import dentist.commands.processPileUps.cropper;
static import dentist.commands.propagateMask;
static import dentist.commands.showInsertions;
static import dentist.commands.showMask;
static import dentist.commands.showPileUps;
static import dentist.commands.translateCoords;
static import dentist.commands.validateConfig;
static import dentist.commands.validateRegions;
static import dentist.common;
static import dentist.common.alignments;
static import dentist.common.alignments.base;
static import dentist.common.alignments.chaining;
static import dentist.common.binio;
static import dentist.common.binio._testdata.insertiondb;
static import dentist.common.binio._testdata.pileupdb;
static import dentist.common.binio.common;
static import dentist.common.binio.insertiondb;
static import dentist.common.binio.pileupdb;
static import dentist.common.commands;
static import dentist.common.configfile;
static import dentist.common.external;
static import dentist.common.insertions;
static import dentist.common.scaffold;
static import dentist.dazzler;
static import dentist.modules;
static import dentist.swinfo;
static import dentist.util.algorithm;
static import dentist.util.containers;
static import dentist.util.fasta;
static import dentist.util.graphalgo;
static import dentist.util.jsonschema;
static import dentist.util.log;
static import dentist.util.math;
static import dentist.util.process;
static import dentist.util.range;
static import dentist.util.region;
static import dentist.util.saturationmath;
static import dentist.util.string;
static import dentist.util.tempfile;


alias modules = AliasSeq!(
    dentist,
    dentist.commandline,
    dentist.commands,
    dentist.commands.bed2mask,
    dentist.commands.buildPartialAssembly,
    dentist.commands.chainLocalAlignments,
    dentist.commands.checkResults,
    dentist.commands.checkScaffolding,
    dentist.commands.collectPileUps,
    dentist.commands.collectPileUps.filter,
    dentist.commands.collectPileUps.pileups,
    dentist.commands.filterMask,
    dentist.commands.findClosableGaps,
    dentist.commands.generateDazzlerOptions,
    dentist.commands.maskRepetitiveRegions,
    dentist.commands.mergeInsertions,
    dentist.commands.mergeMasks,
    dentist.commands.output,
    dentist.commands.processPileUps,
    dentist.commands.processPileUps.cropper,
    dentist.commands.propagateMask,
    dentist.commands.showInsertions,
    dentist.commands.showMask,
    dentist.commands.showPileUps,
    dentist.commands.translateCoords,
    dentist.commands.validateConfig,
    dentist.commands.validateRegions,
    dentist.common,
    dentist.common.alignments,
    dentist.common.alignments.base,
    dentist.common.alignments.chaining,
    dentist.common.binio,
    dentist.common.binio._testdata.insertiondb,
    dentist.common.binio._testdata.pileupdb,
    dentist.common.binio.common,
    dentist.common.binio.insertiondb,
    dentist.common.binio.pileupdb,
    dentist.common.commands,
    dentist.common.configfile,
    dentist.common.external,
    dentist.common.insertions,
    dentist.common.scaffold,
    dentist.dazzler,
    dentist.modules,
    dentist.swinfo,
    dentist.util.algorithm,
    dentist.util.containers,
    dentist.util.fasta,
    dentist.util.graphalgo,
    dentist.util.jsonschema,
    dentist.util.log,
    dentist.util.math,
    dentist.util.process,
    dentist.util.range,
    dentist.util.region,
    dentist.util.saturationmath,
    dentist.util.string,
    dentist.util.tempfile,
);
