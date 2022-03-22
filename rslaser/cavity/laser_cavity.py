from rslaser.optics import element
from rslaser.pulse import pulse
from array import array
from pykern.pkcollections import PKDict
from rslaser.utils.validator import ValidatorBase

class InvalidLaserCavityInputError(Exception):
    pass


class LaserCavity(ValidatorBase):
    """
    create laser cavity

    Args:
        params (PKDict):
            required fields:
                drift_right_length (float): length between crystal and right mirror [m]
                drift_left_length (float): length between crystal and left mirror [m]
                lens_left_focal_length (float): focal length of mirror on left (represented as lens) [m]
                lens_right_focal_length (float): focal length of mirror on right (represented as lens) [m]
                n0 (float): crystal index of refraction on axis
                n2 (float): crystal index of refraction quadratic dependence n(r) = n0 - 1/2 n2 r^2 [m^-2]
                L_half_cryst (float): half length of crystal [m]
                crystal_slices: to be passed as nslice to Crystal
                pulse_params (PKDict): see LaserPulse docs
    """
    _DEFAULTS = PKDict(
        drift_right_length=0.5,
        drift_left_length=0.5,
        lens_left_focal_length=0.2,
        lens_right_focal_length=0.2,
        n0 = 1.75,
        n2 = 0.001,
        L_half_cryst=0.2,
        crystal_slices=3,
        pulse_params=PKDict()
    )
    _INPUT_ERROR = InvalidLaserCavityInputError

    def __init__(self, params=None):
        params = self._get_params(params)
        self._validate_params(params)
        crystal_params = PKDict(
            n0=params.n0,
            n2=params.n2,
            length=params.L_half_cryst,
            nslice=params.crystal_slices,
        )
        self.laser_pulse = pulse.LaserPulse(params.pulse_params)
        self.crystal_right = element.Crystal(crystal_params)
        self.crystal_left = element.Crystal(crystal_params)
        self.drift_right = element.Drift(params.drift_right_length)
        self.drift_left = element.Drift(params.drift_left_length)
        self.lens_right = element.Lens(params.lens_right_focal_length)
        self.lens_left  = element.Lens(params.lens_left_focal_length)

    def propagate(self, num_cycles, callback=None):
        l = self.laser_pulse
        l._sxvals = []
        l._syvals = []
        current_position = 0

        # initial wavefront rms values
        vals = l.compute_middle_slice_intensity()
        positions = [current_position]
        if callback:
            callback(current_position, vals)

        for n in range(num_cycles):
            for e in (
                self.crystal_right,
                self.drift_right,
                self.lens_right,
                self.drift_right,
                self.crystal_right,
                self.crystal_left,
                self.drift_left,
                self.lens_left,
                self.drift_left,
                self.crystal_left,
            ):
                if type(e) == element.Crystal:
                    e.propagate(l, 'abcd')
                else:
                    e.propagate(l)
                current_position += e.length
                vals = l.compute_middle_slice_intensity()
                positions.append(current_position)
                if callback:
                    callback(current_position, vals)

        return(positions, l._sxvals, l._syvals)
