u"""Tests for Crystal and CrystalSlice
"""
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
import pykern.pkunit
import pytest
from rslaser.optics import element, lens, drift, crystal
from rslaser.pulse import pulse


def test_instantiation01():
    crystal.Crystal()
    with pykern.pkunit.pkexcept(element.ElementException):
        crystal.Crystal('fail')
    with pykern.pkunit.pkexcept(element.ElementException):
        crystal.Crystal(PKDict(slice_params=PKDict()))
    c = crystal.Crystal()
    for s in c.slice:
        if s.length != c.length/c.nslice:
            pykern.pkunit.pkfail('CrystalSlice had length not equal to Crystal wrapper length/nslice')


def test_crystal_nslice():
    crystal.Crystal(PKDict(nslice=51))
    crystal.Crystal(PKDict(n0=[1], n2=[1]))
    with pykern.pkunit.pkexcept(element.ElementException, "you've specified"):
        crystal.Crystal(PKDict(nslice=51, n0=[1], n2=[1]))


def crystal_slice_prop_test(prop_type):
    c = crystal.CrystalSlice()
    p = pulse.LaserPulse()
    p = c.propagate(p, prop_type)
    if type(p) != pulse.LaserPulse:
        pykern.pkunit.pkfail('Crystal slice propagaition failed to return LaserPulse type')


def test_instantiation02():
    c = crystal.CrystalSlice()


def test_propagation():
    with pykern.pkunit.pkexcept(element.ElementException):
        crystal_slice_prop_test('default')
    with pykern.pkunit.pkexcept(NotImplementedError):
        crystal_slice_prop_test('attenuate')
    with pykern.pkunit.pkexcept(NotImplementedError):
        crystal_slice_prop_test('placeholder')
    c = crystal.CrystalSlice()
    with pykern.pkunit.pkexcept(KeyError):
        c.propagate(pulse.LaserPulse(), 'should raise')


def test_prop_with_gain():
    # TODO (gurhar1133): add ndiff data
    from pykern import pkio

    data_dir = pykern.pkunit.data_dir()
    def _prop(prop_type):

        c = crystal.Crystal(
            PKDict(
                n2=[16],
                # TODO (gurhar): change default to 0.01
                l_scale=0.001,
            )
        )
        p = pulse.LaserPulse()
        c.propagate(p, prop_type, calc_gain=True)
        r = ""
        # for s in p.slice:
        #     for k in s.__dict__:
        #         r += " " + k + " "
        for s in p.slice:
            for k in s.wfr.__dict__:
                r += f" {k} "
        # r = [x.wfr.__dict__ for x in p.slice]

        pkio.write_text("TEST_DATA"+prop_type+".ndiff", str(r))
        # print(p.slice[0].__dict__)
    _prop("n0n2_srw")
    # _prop("n0n2_lct")
    # _prop("gain_calc")


def test_instantiation03():
    drift.Drift(0.01)


def test_propagation05():
    d = drift.Drift(0.01)
    p = pulse.LaserPulse()
    d.propagate(p)
    trigger_prop_fail(
        drift.Drift(0.01).propagate,
        pulse.LaserPulse()
        )


def test_instantiation04():
    lens.Lens(0.2)


def test_propagation06():
    l = lens.Lens(0.2)
    l.propagate(pulse.LaserPulse())
    trigger_prop_fail(
        lens.Lens(0.01).propagate,
        pulse.LaserPulse()
        )

def trigger_prop_fail(prop_func, pulse):
    with pykern.pkunit.pkexcept(
        element.ElementException,
        'Invalid element="should raise" should have raised'
        ):
        prop_func(pulse, 'should raise')
