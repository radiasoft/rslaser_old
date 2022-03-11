from rslaser.optics import element
from rslaser.pulse import pulse
from array import array
from pykern.pkcollections import PKDict


class InvalidLaserCavityInputError(Exception):
    pass


class LaserCavity(pulse.LaserBase):
    """
    create laser cavity

    Args:
        params (PKDict):
            required fields:
                drift_right_length
                drift_left_length
                lens_left_focal_length
                lens_right_focal_length
                n0
                n2
                L_half_cryst
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
        pulse_params=PKDict()
    )
    _INPUT_ERROR = InvalidLaserCavityInputError

    def __init__(self, params=None):
        params = self._get_params(params)
        self._validate_params(params)
        self.laser_pulse = pulse.LaserPulse(params.pulse_params)
        self.crystal_right = element.Crystal(params.n0,params.n2,params.L_half_cryst)
        self.crystal_left = element.Crystal(params.n0,params.n2,params.L_half_cryst)
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
            for element in (
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
                element.propagate(l)
                current_position += element.length
                vals = l.compute_middle_slice_intensity()
                positions.append(current_position)
                if callback:
                    callback(current_position, vals)

        return(positions, l._sxvals, l._syvals)
