u"""Test utils
"""
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict


def trigger_exception_test(call, args):
    try:
        if len(args) > 0 and type(args) != PKDict and type(args) != str:
            l = call(*args)
        else:
            l = call(args)
    except Exception as e:
        pkdlog('EXCEPTION:{}, with message "{}" triggered by call: {} and args {}', type(e), e, call, args)
        return e
    assert False