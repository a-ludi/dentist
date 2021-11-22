/**
    Exposes the `execute` function from each command. The functions are
    accessed by the command name, for example,
    `dentist.commands.collectPileUps.execute` is available as `collectPileUps`.


    $(UL
        $(LI `dentist.commands.bed2mask`)
        $(LI `dentist.commands.chainLocalAlignments`)
        $(LI `dentist.commands.collectPileUps`)
        $(LI `dentist.commands.filterMask`)
        $(LI `dentist.commands.generateDazzlerOptions`)
        $(LI `dentist.commands.maskRepetitiveRegions`)
        $(LI `dentist.commands.mergeInsertions`)
        $(LI `dentist.commands.mergeMasks`)
        $(LI `dentist.commands.output`)
        $(LI `dentist.commands.processPileUps`)
        $(LI `dentist.commands.propagateMask`)
        $(LI `dentist.commands.showInsertions`)
        $(LI `dentist.commands.showMask`)
        $(LI `dentist.commands.showPileUps`)
        $(LI `dentist.commands.translateCoords`)
        $(LI `dentist.commands.validateConfig`)
        $(LI `dentist.commands.validateRegions`)
    )

    If DENTIST is compiled with `--config=testing` the following additional
    commands will be available:

    $(UL
        $(LI `dentist.commands.buildPartialAssembly`)
        $(LI `dentist.commands.checkResults`)
        $(LI `dentist.commands.checkScaffolding`)
        $(LI `dentist.commands.findClosableGaps`)
        $(LI `dentist.commands.translocateGaps`)
    )

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands;

public import dentist.commands.bed2mask : bed2mask = execute;
public import dentist.commands.buildPartialAssembly : buildPartialAssembly = execute;
public import dentist.commands.chainLocalAlignments : chainLocalAlignments = execute;
public import dentist.commands.checkResults : checkResults = execute;
public import dentist.commands.checkScaffolding : checkScaffolding = execute;
public import dentist.commands.collectPileUps : collectPileUps = execute;
public import dentist.commands.filterMask : filterMask = execute;
public import dentist.commands.findClosableGaps : findClosableGaps = execute;
public import dentist.commands.generateDazzlerOptions : generateDazzlerOptions = execute;
public import dentist.commands.maskRepetitiveRegions : maskRepetitiveRegions = execute;
public import dentist.commands.mergeInsertions : mergeInsertions = execute;
public import dentist.commands.mergeMasks : mergeMasks = execute;
public import dentist.commands.output : output = execute;
public import dentist.commands.processPileUps : processPileUps = execute;
public import dentist.commands.propagateMask : propagateMask = execute;
public import dentist.commands.showInsertions : showInsertions = execute;
public import dentist.commands.showMask : showMask = execute;
public import dentist.commands.showPileUps : showPileUps = execute;
public import dentist.commands.translateCoords : translateCoords = execute;
public import dentist.commands.translocateGaps : translocateGaps = execute;
public import dentist.commands.validateConfig : validateConfig = execute;
public import dentist.commands.validateRegions : validateRegions = execute;
