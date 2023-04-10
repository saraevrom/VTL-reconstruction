import unittest
import numpy.random as rng
import numpy as np
from stellar_math import radec_to_eci, eci_to_ecef, eci_to_ocef


def generate_random_direction():
    dec = (rng.rand()*2-1)*np.pi/2
    ra = (rng.rand()*2-1)*np.pi
    return radec_to_eci(ra, dec)

class StarTest(unittest.TestCase):
    def test_radec_length(self):
        x,y,z = generate_random_direction()
        self.assertAlmostEqual(x**2+y**2+z**2, 1.0)

    def test_eci_negation(self):
        dec = (rng.rand() * 2 - 1) * np.pi / 2
        ra = (rng.rand() * 2 - 1) * np.pi
        x,y,z = radec_to_eci(ra, dec)
        _,y1,_ = eci_to_ecef(x, y, z, ra)
        self.assertAlmostEqual(y1, 0.0)

    def test_eci_ecef_z_consistency(self):
        x_eci, y_eci, z_eci = generate_random_direction()
        era = rng.rand()*2*np.pi
        x_ecef, y_ecef, z_ecef = eci_to_ecef(x_eci, y_eci, z_eci, era)
        self.assertAlmostEqual(z_ecef, z_eci)

    def test_eci_ecef_rotation_consistency(self):
        x_eci, y_eci, z_eci = generate_random_direction()
        era = rng.rand() * 2 * np.pi
        x_ecef, y_ecef, z_ecef = eci_to_ecef(x_eci, y_eci, z_eci, era)
        x_eci1, y_eci1, z_eci1 = eci_to_ecef(x_ecef, y_ecef, z_ecef, -era)
        self.assertAlmostEqual(x_eci, x_eci1)
        self.assertAlmostEqual(y_eci, y_eci1)
        self.assertAlmostEqual(z_eci, z_eci1)

    def test_zenith(self):
        dec = (rng.rand() * 2 - 1) * np.pi / 2
        ra = (rng.rand() * 2 - 1) * np.pi
        era = 0.0
        x_eci, y_eci, z_eci = radec_to_eci(ra, dec)
        x_local, y_local, z_local = eci_to_ocef(x_eci, y_eci, z_eci, era, lat=dec, lon=ra)
        self.assertAlmostEqual(x_local, 1.0)

    def test_north(self):
        dec = (rng.rand() * (0.98+1) - 1) * np.pi / 2
        ra = (rng.rand() * 2 - 1) * np.pi
        era = 0.0
        x_eci, y_eci, z_eci = radec_to_eci(ra, dec+0.01)
        x_local, y_local, z_local = eci_to_ocef(x_eci, y_eci, z_eci, era, lat=dec, lon=ra)
        self.assertGreater(z_local, 0.0)
        self.assertAlmostEqual(y_local, 0.0)

    def test_east(self):
        dec = (rng.rand() * 2 - 1) * np.pi / 2
        ra = (rng.rand() * 2 - 1) * np.pi
        era = 0.0
        x_eci, y_eci, z_eci = radec_to_eci(ra+0.09, dec)
        x_local, y_local, z_local = eci_to_ocef(x_eci, y_eci, z_eci, era, lat=dec, lon=ra)
        self.assertGreater(y_local, 0.0)
        # self.assertAlmostEqual(z_local, 0.0)
