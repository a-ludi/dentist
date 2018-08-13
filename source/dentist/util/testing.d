/**
    This is a collection of helpers for unit testing.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.testing;
import std.traits : isCallable, Parameters, ReturnType;

template MockCallable(T...)
{
    static if (T.length == 1 && isCallable!T)
    {
        alias R = ReturnType!T;
        alias Args = Parameters!T;
    }
    else
    {
        alias R = T[0];
        alias Args = T[1 .. $];
    }

    struct MockCallable
    {

        static struct CallInfo
        {
            R returnValue;
            Args args;
        }

        R returnValue;
        CallInfo[] calls;

        R opCall(Args args)
        {
            calls ~= CallInfo(returnValue, args);

            return returnValue;
        }

        @property bool wasCalled()
        {
            return calls.length > 0;
        }

        void reset()
        {
            R returnValue;

            this.returnValue = returnValue;
            this.calls.length = 0;
        }
    }
}

unittest
{
    int foo(char, string)
    {
        return 0;
    }

    MockCallable!foo fooMock;

    static assert(is(typeof(fooMock(' ', "")) == int));
    fooMock.returnValue = 42;
    assert(fooMock('1', "1") == 42);
    assert(fooMock.calls[0].returnValue == 42);
    assert(fooMock.calls[0].args[0] == '1');
    assert(fooMock.calls[0].args[1] == "1");
    fooMock.returnValue = 1337;
    assert(fooMock('2', "1") == 1337);
    assert(fooMock.calls.length == 2);

    fooMock.reset();
    assert(fooMock.calls.length == 0);
    assert(fooMock.returnValue == 0);
}
