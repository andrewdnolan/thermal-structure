import numpy  as np

def __calc_coefs(shape, axis, dist, edge_order=1): 
    """
    input: 
        z (nz, nx) : 
    
    output: 
        c (nz, nz, nx): i.e. a (nz,nz) coefficient matrix for every nx point in the domain
        
    input a two dimmesion Z coordinate array
    
    uses the same variable coefficient formula as numpy for interior points:
        https://numpy.org/doc/stable/reference/generated/numpy.gradient.html#

    used to calculate the gridcell spacing, and with those spacing fill in the appropriate finite difference operator matrix 
    """
    # unpack the shape tuple
    nz, nx = shape
    
    # determine the size of the interior stencil matrix 
    if axis==0: 
        # i.e. taking gradient along vertical (z) coordinate
        int_m = nz-2
        int_n = nz
        # get the vertical gridcell spacing 
        hd = dist[ 1:, :]  # + an element, vertically 
        hs = dist[:-1, :]  # - an element, vertically 
        # summation indexes for doing array broadcasting 
        einsum_path = 'ij,ik->ikj'
    elif axis==1: 
        # i.e. taking gradient along horizontal (x) coordinate
        int_m = nx-2
        int_n = nx
        # get the vertical gridcell spacing 
        hd = dist[:, 1 :]  # + an element, horizontally 
        hs = dist[:, :-1]  # - an element, horizontally
        # summation indexes for doing array broadcasting 
        einsum_path = 'ij,jk->ijk'

    # create the template interior stencil array (nz-2, nz)
    # i.e. it's the interior stencil for a single x coordinate
    sup  =  np.eye(int_m, int_n, 2)
    main =  np.eye(int_m, int_n, 1)
    sub  = -np.eye(int_m, int_n, 0) 
    # ref: https://stackoverflow.com/a/55513250/10221482

    # get the denomaintor 
    denom = (hs*hd*(hs+hd))

    # broadcast template interior stencil arrays (nz-2, nz)  and 
    # gridcellspacing array(nz-2, nx) to be (nz-2, nz, nx) 
    # where: 
    #      i -- nz-2
    #      j -- nx
    #      k -- nz
    interior = np.einsum(einsum_path, hs**2/denom, sup)  + \
               np.einsum(einsum_path, (hd**2 - hs**2)/denom, main) +\
               np.einsum(einsum_path, hd**2/denom, sub)

    # create the full coefs vector (nz, nz, nx),
    # i.e. a (nz, nz) matrix for every x gridpoint 
    # which is empty for now 
    if axis==0: 
        coefs = np.zeros((nz,nz,nx))
        # dump the interior matrix, i.e. centered finite difference 
        coefs[1:-1, :, :] = interior

        # use a first order difference along the boundaries:
        # to add support for three or 5 point differences at boudnary, check: 
        # https://docs.sympy.org/latest/explanation/special_topics/finite_diff_derivatives.html#a-direct-method-using-sympy-matrices
        # for deriving the coefficents using variable spacing
        coefs[0,  0:2, :] = (np.array([[-1], [1]]) * np.ones((1,nx))) / dist[0,  :]
        coefs[-1, -2:, :] = (np.array([[-1], [1]]) * np.ones((1,nx))) / dist[-1, :]
    
    elif axis==1: 
        coefs = np.zeros((nz,nx,nx))
        # dump the interior matrix, i.e. centered finite difference 
        coefs[:,1:-1,:] = interior
        
        # use a first order difference along the boundaries
        coefs[:, 0,  0:2] = (np.array([[-1, 1]]) * np.ones((nz,1))) / dist[:,   0:1]
        coefs[:, -1, -2:] = (np.array([[-1, 1]]) * np.ones((nz,1))) / dist[:, -2:-1]
    
    return coefs

def calc_gradient(F, Z, X): 
    """ Take the gradient of a multidimensional scalar field. 
        
        Note: this function does not support time dependent data. 
              To calculate the gradient as a fucntion of time first group by 'time' 
              then apply this function over:
              ```
              da.groupby('t').apply(graient)
              ```
        where: 

            da (nz, nx) --> is the multidimensional dataarray to take the gradient of 
            X  (nz, nx) --> is the x coordinate array, defined at every node
            Z  (nz, nx) --> is the z coordinate array, defined ar every node
    
    """
    assert X.shape == Z.shape
    
    nz, nx = Z.shape
    
    # vertical gridceel spacing (nz-1, nx)
    dz = np.diff(Z, axis=0)
    # finite difference coefficent matrix (nz, nz, nx)
    coefs = __calc_coefs((nz,nx), 0, dz)
    # compute the gradient!!
    grad_z = np.einsum('ijk, jk -> ik', coefs, F)

    # horizontal gridcell spacing (nz, nx-1)
    dx = np.diff(X, axis=1)
    # finite difference coefficent matrix (nz, nx, nx)
    coefs = __calc_coefs((nz,nx), 1, dx)
    # compute the gradient!!
    grad_x = np.einsum('ijk, ik -> ij', coefs, F)
    
    
    return grad_z, grad_x