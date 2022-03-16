u"""Tests for Crystal and CrystalSlice
"""
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
import pytest
from rslaser.optics import element
from rslaser.pulse import pulse
import test_utils

# TODO: (gurhar1133): may want to reorganize the test files
# ie. move some stuff out of instantiation_test and create a
# pulse_test.py file

def crystal_slice_prop_test(prop_type):
    c = element.CrystalSlice('test', 0.01)
    p = pulse.LaserPulse()
    p = c.propagate(p, prop_type)
    assert type(p) == pulse.LaserPulse


def test_crystal_slice_instantiation():
    c = element.CrystalSlice('test', 0.01)


def test_crystal_slice_propagate_default():
    crystal_slice_prop_test('default')


def test_crystal_slice_propagate_attenuate():
    crystal_slice_prop_test('attenuate')


def test_crystal_slice_propagate_placeholder():
    crystal_slice_prop_test('placeholder')


def test_crystal_slice_propagate_abcd():
    crystal_slice_prop_test('abcd')


def test_crystal_slice_propagate_exception():
    c = element.CrystalSlice('test', 0.01)
    test_utils.trigger_exception_test(c.propagate, 'should fail')


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
