u"""Tests for Crystal and CrystalSlice
"""
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
import pykern.pkunit
import pytest
from rslaser.optics import element
from rslaser.pulse import pulse


def test_instantiation01():
    element.Crystal()
    with pykern.pkunit.pkexcept(element.ElementException):
        element.Crystal('fail')
    with pykern.pkunit.pkexcept(element.ElementException):
        element.Crystal(PKDict(slice_params=PKDict()))
    c = element.Crystal()
    for s in c.slice:
        if s.length != c.length/c.nslice:
            pykern.pkunit.pkfail('CrystalSlice had length not equal to Crystal wrapper length/nslice')


def crystal_slice_prop_test(prop_type):
    c = element.CrystalSlice()
    p = pulse.LaserPulse()
    p = c.propagate(p, prop_type)
    if type(p) != pulse.LaserPulse:
        pykern.pkunit.pkfail('Crystal slice propagaition failed to return LaserPulse type')


def test_instantiation02():
    c = element.CrystalSlice()


def test_propagation():
    with pykern.pkunit.pkexcept(element.ElementException):
        crystal_slice_prop_test('default')
    with pykern.pkunit.pkexcept(NotImplementedError):
        crystal_slice_prop_test('attenuate')
    with pykern.pkunit.pkexcept(NotImplementedError):
        crystal_slice_prop_test('placeholder')
    # TODO (gurhar1133): propagation is a work in progress.
    # crystal_slice_prop_test('abcd_lct')
    # crystal_slice_prop_test('n0n2')
    c = element.CrystalSlice()
    with pykern.pkunit.pkexcept(KeyError):
        c.propagate(pulse.LaserPulse(), 'should raise')


def test_instantiation03():
    element.Drift(0.01)


def test_propagation05():
    d = element.Drift(0.01)
    p = pulse.LaserPulse()
    d.propagate(p)
    trigger_prop_fail(
        element.Drift(0.01).propagate,
        pulse.LaserPulse()
        )


def test_instantiation04():
    element.Lens(0.2)


def test_propagation06():
    l = element.Lens(0.2)
    l.propagate(pulse.LaserPulse())
    trigger_prop_fail(
        element.Lens(0.01).propagate,
        pulse.LaserPulse()
        )

def trigger_prop_fail(prop_func, pulse):
    with pykern.pkunit.pkexcept(
        element.ElementException,
        'Invalid element="should raise" should have raised'
        ):
        prop_func(pulse, 'should raise')
