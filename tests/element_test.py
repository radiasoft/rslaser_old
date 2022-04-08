u"""Tests for Crystal and CrystalSlice
"""
from ast import Assert
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
import pytest
from rslaser.optics import element
from rslaser.pulse import pulse
import test_utils


def test_crystal_instantiation():
    c = element.Crystal(
        PKDict(
            n0=.3,
        )
    )


def test_crystal_instantiation2():
    c = element.Crystal()


def test_crystal_instantiation3():
    test_utils.trigger_exception_test(element.Crystal, 'fail')


def test_crystal_instantiation4():
    test_utils.trigger_exception_test(element.Crystal, PKDict(slice_params=PKDict()))


def test_crystal_instantiation5():
    c = element.Crystal()
    for s in c.slice:
        if s.length != c.length/c.nslice:
            raise AssertionError('CrystalSlice had length not equal to Crystal wrapper length/nslice')


def test_crystal_propagation():
    p = pulse.LaserPulse()
    c = element.Crystal()
    p = c.propagate(p, 'abcd')
    assert type(p) == pulse.LaserPulse


def crystal_slice_prop_test(prop_type):
    c = element.CrystalSlice()
    p = pulse.LaserPulse()
    p = c.propagate(p, prop_type)
    assert type(p) == pulse.LaserPulse


def test_crystal_slice_instantiation():
    c = element.CrystalSlice()


def test_crystal_slice_propagate_default():
    test_utils.trigger_exception_test(crystal_slice_prop_test, 'default')


def test_crystal_slice_propagate_attenuate():
    test_utils.trigger_exception_test(crystal_slice_prop_test, 'attenuate')


def test_crystal_slice_propagate_placeholder():
    test_utils.trigger_exception_test(crystal_slice_prop_test, 'placeholder')


def test_crystal_slice_propagate_abcd():
    crystal_slice_prop_test('abcd')


def test_crystal_slice_propagate_exception():
    c = element.CrystalSlice()
    test_utils.trigger_exception_test(c.propagate, [pulse.LaserPulse(), 'should fail'])


def test_drift_instantiate():
    element.Drift(0.01)


def test_drift_propagate():
    d = element.Drift(0.01)
    p = pulse.LaserPulse()
    d.propagate(p)


def test_drift_propagate_fail():
    trigger_prop_fail(
        element.Drift(0.01).propagate,
        pulse.LaserPulse()
        )


def test_lense_instantiate():
    element.Lens(0.2)


def test_lens_propagate():
    l = element.Lens(0.2)
    l.propagate(pulse.LaserPulse())


def test_lens_propagate_fail():
    trigger_prop_fail(
        element.Lens(0.01).propagate,
        pulse.LaserPulse()
        )

def trigger_prop_fail(prop_func, pulse):
    test_utils.trigger_exception_test(
        prop_func,
        [pulse, 'should fail']
    )
