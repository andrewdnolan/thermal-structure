import numpy as np 
from scipy import linalg as LA
from thermal.gradient import calc_gradient


def test_constant_field(): 
    """ gradient of a constant scalar field should be ==0 
        for all dimensions. 

        Still, test calc_gradient function: 
            - by irregularly spaced x and y coordinates
            - a non-square matrix
    """

    # dimension length:
    # make non square matrices to make things challenging 
    nx, ny = 500, 300

    # have coordinates be irregularly spaced
    x = np.linspace(0, 1, nx)**2
    y = np.linspace(0, 1, ny)**4

    # mesh the dimensions into (nx,ny) coordinate array
    xx, yy = np.meshgrid(x,y)

    # input (nx, ny) field
    F = 555 * np.zeros_like(xx)

    # in house function for calculating gradient 
    dFdy, dFdx = calc_gradient(F, yy, xx)

    # gradient should be zero, so just check norm of gradient 
    # instead of norm of difference
    assert LA.norm(dFdx) <= 1e-6
    assert LA.norm(dFdy) <= 1e-6

def test_variable_field(): 
    nx, ny = 100, 1000

    # regularly spaced dimensions of the 2-D data 
    x = np.linspace(0, np.pi, nx)
    y = np.linspace(0, np.pi, ny)


    # mesh the dimensions into (nx,ny) coordinate array
    xx, yy = np.meshgrid(x,y)

    # input (nx, ny) field
    F =  np.sin(xx*yy) + np.cos(xx*yy)

    # we know the analytical gradient 
    dFdx_anl = yy*(np.cos(xx*yy) - np.sin(xx*yy))
    dFdy_anl = xx*(np.cos(xx*yy) - np.sin(xx*yy))

    # create empty array to store numpy grad
    dFdx_np = np.zeros_like(dFdx_anl)
    dFdy_np = np.zeros_like(dFdy_anl)

    for i in range(ny): 
        dFdx_np[i,:] = np.gradient(F[i,:], xx[i,:])

    for j in range(nx): 
        dFdy_np[:,j] = np.gradient(F[:,j], yy[:,j])

    # in house function for calculating gradient 
    dFdy, dFdx = calc_gradient(F, yy, xx)

    # print(LA.norm(dFdx_anl) - LA.norm(dFdx_np))
    # print(LA.norm(dFdy_anl) - LA.norm(dFdy_np))

    # make sure the in house gradient returns the same result as numpy
    assert LA.norm(dFdx - dFdx_np) <= 1e-6
    assert LA.norm(dFdy - dFdy_np) <= 1e-6