{
    if ($0 ~ /full_validation/)
        # turn on full validation
        print "full_validation: true"
    else if ($0 ~ /dentist_container/ && typeof(dentist_container) == "string")
        # adjust Singularity container
        print "dentist_container: " dentist_container
    else if ($0 ~ /dentist_env/ && typeof(dentist_env) == "string")
        # adjust Conda env
        print "dentist_env: " dentist_env
    else if (match($0, /(^[[:blank:]]+- )-s[0-9]+/, matchGroups))
        # reduce DB block size to force multiple blocks
        print matchGroups[1] "-s50"
    else if ($0 ~ /propagate_batch_size/)
        # reduce batch size to force multiple batches
        print "propagate_batch_size: 14"
    else if ($0 ~ /validation_blocks/)
        # reduce block count to avoid empty blocks
        print "validation_blocks: 2"
    else if ($0 ~ /^[[:blank:]]*#/ || length($0) == 0)
        # exclude comments
        next;
    else
        # copy line
        print
}