import unittest
import numpy as np
from scipy import interpolate
from scipy import linalg as LA

class TestFitpack(unittest.TestCase):

    def test_fitpack(self):
        # load the tck tuple written to disk
        t_vec = np.loadtxt('./data/knots.dat')
        c_vec = np.loadtxt('./data/coefs.dat')
        x_vec = np.loadtxt('./data/x_vec.dat')

        # create the tck tuple
        tck = (t_vec, c_vec, 3)

        # load the spline evaluated in fortran
        y_for = np.loadtxt('./data/y_vec.dat')

        # evaluate the spline in python
        y_sci = interpolate.splev(x_vec, tck)

        # make sure the fortran and python implementations
        # are within some desired tolerance
        self.assertTrue(LA.norm(y_sci - y_for, np.inf) < 1e-3)

if __name__ == '__main__':
    unittest.main()
