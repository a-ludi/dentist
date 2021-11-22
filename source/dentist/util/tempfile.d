/**
    Create temporary files and folders. These are wrappers around the
    corresponding functions in `core.sys.posix.stdlib`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.tempfile;

version (Posix)
{
    import std.algorithm : endsWith;
    import std.conv : to;
    import std.exception : errnoEnforce;
    import std.stdio : File;
    import std.string : fromStringz;
    import std.typecons : Tuple, tuple;

    /**
        Generates a uniquely named temporary directory from `templateString`.

        The last six characters of template must be XXXXXX and these are
        replaced with a string that makes the directory name unique. The
        directory is then created with permissions 0700.

        Returns: The generated directory name.
    */
    string mkdtemp(in string templateString) @trusted
    {
        import core.sys.posix.stdlib : c_mkdtemp = mkdtemp;

        char[255] dirnameBuffer;
        auto len = templateString.length;
        assert(len < dirnameBuffer.length);
        assert(templateString.endsWith("XXXXXX"));

        dirnameBuffer[0 .. len] = templateString[];
        dirnameBuffer[len] = 0;

        errnoEnforce(null != c_mkdtemp(dirnameBuffer.ptr), "cannot create temporary directory");

        return to!string(fromStringz(dirnameBuffer.ptr));
    }

    ///
    unittest
    {
        import std.algorithm : startsWith;
        import std.file : isDir, rmdir;

        string tempDir = mkdtemp(".unittest-XXXXXX");

        try
        {
            assert(isDir(tempDir));
            assert(tempDir.startsWith(".unittest-"));
        }
        finally
        {
            rmdir(tempDir);
        }
    }


    // Glibc function since 2.19
    private extern (C) int mkstemps(char*, int);


    /**
        Generates a unique temporary filename from `templateString`, creates
        and opens the file, and returns the open file and generated name.

        The last six characters of template must be "XXXXXX" and these are
        replaced with a string that makes the filename unique.

        The optional `templateSuffix` will be appended to the file name.

        Returns: The open file and generated name.
    */
    Tuple!(File, "file", string, "name") mkstemp(in string templateString) @trusted
    {
        import core.sys.posix.stdlib : mkstemp;

        char[255] dirnameBuffer;
        auto len = templateString.length;
        assert(len < dirnameBuffer.length);
        assert(templateString.endsWith("XXXXXX"));

        dirnameBuffer[0 .. len] = templateString[];
        dirnameBuffer[len] = 0;

        auto fd = mkstemp(dirnameBuffer.ptr);

        errnoEnforce(fd != -1, "cannot create temporary file");

        File tempFile;
        tempFile.fdopen(fd, "r+");

        return tuple!("file", "name")(tempFile, dirnameBuffer.ptr.fromStringz.to!string);
    }

    ///
    unittest
    {
        import std.algorithm : startsWith;
        import std.file : remove;

        auto tempFile = mkstemp(".unittest-XXXXXX");
        scope (exit)
            remove(tempFile.name);

        assert(tempFile.name.startsWith(".unittest-"));
        assert(tempFile.file.isOpen);
        assert(!tempFile.file.error);
        tempFile.file.writeln("foobar");
        tempFile.file.flush();
        tempFile.file.rewind();
        assert(tempFile.file.readln() == "foobar\n");
    }

    /// ditto
    Tuple!(File, "file", string, "name") mkstemp(in string templateString, in string templateSuffix) @trusted
    {
        char[255] dirnameBuffer;
        auto len = templateString.length + templateSuffix.length;
        assert(len < dirnameBuffer.length);
        assert(templateString.endsWith("XXXXXX"));

        dirnameBuffer[0 .. len] = templateString[] ~ templateSuffix[];
        dirnameBuffer[len] = 0;

        auto fd = mkstemps(dirnameBuffer.ptr, templateSuffix.length.to!int);

        errnoEnforce(fd != -1, "cannot create temporary file");

        File tempFile;
        tempFile.fdopen(fd, "r+");

        return tuple!("file", "name")(tempFile, dirnameBuffer.ptr.fromStringz.to!string);
    }

    ///
    unittest
    {
        import std.algorithm : endsWith, startsWith;
        import std.file : remove;

        auto tempFile = mkstemp(".unittest-XXXXXX", ".ext");
        scope (exit)
            remove(tempFile.name);

        assert(tempFile.name.startsWith(".unittest-"));
        assert(tempFile.name.endsWith(".ext"));
        assert(tempFile.file.isOpen);
        assert(!tempFile.file.error);
        tempFile.file.writeln("foobar");
        tempFile.file.flush();
        tempFile.file.rewind();
        assert(tempFile.file.readln() == "foobar\n");
    }
}
else
{
    static assert(0, "Only intended for use on POSIX compliant OS.");
}
