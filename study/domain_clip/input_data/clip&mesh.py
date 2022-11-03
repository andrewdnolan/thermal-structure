#!/usr/bin/env python3


import sys
sys.path.append('../../../src')

import os
import glob
import subprocess
import numpy as np
import thermal.mesh as meshing
import matplotlib.pyplot as plt

# plt.rcParams['text.usetex'] = True

# Deal with relative path problems when using a module
thermal_fp = os.path.dirname(__file__)
git_top    = os.path.join(thermal_fp, '../../../')

class domain:
    def __init__(self, xx, zs, zb):

        # store the arrays
        self.xx = xx
        self.zs = zs
        self.zb = zb

    def clip(self, zs_max, xx_conect, plot=False, out_fp=None):

        # mask the arrays with the passed parameters
        mask1 = self.zs <= zs_max
        mask2 = self.xx >= xx_conect

        # get the two points to connect
        x1, x2 = (self.xx[mask1 & mask2][0], self.xx[mask1 & mask2][-1])
        y1, y2 = (self.zb[mask1 & mask2][0], self.zs[mask1 & mask2][-1])


        # https://math.stackexchange.com/questions/4014286/exponential-between-2-points
        k = 1.0/200
        a = (y2 - y1) / (np.exp(k * x2) - np.exp(k * x1))
        b = y1 - a * np.exp(k * x1)

        xx_new = self.xx[mask1]
        zs_new = self.zs[mask1]

        zb_new = np.concatenate((self.zb[mask1 & (~mask2)],
                                a * np.exp(k * self.xx[mask1 & mask2]) + b))

        if plot:
            fig, ax = plt.subplots()

            ax.plot(self.xx, self.zs, lw=1.0)
            ax.plot(self.xx, self.zb, lw=1.0)
            ax.scatter([x1, x2], [y1, y2], s=15.0, c='tab:red')

            ax.plot(xx_new, zb_new, lw=1.0)

            plt.show()

        if out_fp:

            np.savetxt(os.path.join(out_fp, f'crmpt12_{int(zs_max)}_clip_bed.dat'),
                       np.hstack([xx_new[:,None], zb_new[:,None]]))

            np.savetxt(os.path.join(out_fp, f'crmpt12_{int(zs_max)}_clip_surf.dat'),
                       np.hstack([xx_new[:,None], zs_new[:,None]]))

        return xx_new, zs_new, zb_new

def make_mesh(key, dx, out_fp, Nz=15, force=False):
    """Make the mesh directory. Currently only suports first-order
       quadrilateral elments (Elmer 404 nodes)

    Inputs:
        key    (str)  --> glacier identifier (see input_data repo)
        dx     (int)  --> horzontal mesh resolution
        out_fp (str)  --> file path to directory to place the mesh.* files
        force  (bool) --> whether to overwrite mesh if it already exists

    return:
        result (...)  --> instance of subprocess.CompletedProcess class
                          has attributes stderr and stdout which record the
                          stderr and stdout respectively from the mesh generation
                          call. If you mesh isn't made correctly check these attrs
    """

    # Really only need to check about one directory up.....
    # ElmerGrid can create the right most directory (i.e. .../mesh_dx{dx})
    git_top_fp = '/'.join(out_fp.split('/')[:-1])
    if not os.path.exists(os.path.join(git_top, git_top_fp)):
        os.makedirs(os.path.join(git_top, git_top_fp))

    # check if the mesh files alrady exists.  no need to create the mesh if it
    # already exists To overwrite exiting mesh set force == True
    if ( os.path.exists(os.path.join(git_top, out_fp)) ) and \
       ( len(glob.glob( os.path.join(git_top, out_fp, 'mesh.*'))) == 4 ) and \
       ( not force ):
       return

    # find the x coordinate of the approraite bed file
    x = np.loadtxt(f'./topography/{key}_bed.dat')[:,0]
    L  = x.max() - x.min()
    Nx = int(round(L/dx))

    # create the glacier specific mesh file
    ElmerGird = meshing.update_template(L, Nx, Nz)

    # write the ElmerGrid string to a temporary file
    with open(os.path.join(git_top,"tmp.grd"), 'w+') as f:
        f.write(ElmerGird)

    pwd = os.getcwd()
    os.chdir(git_top)

    cmd   = "ElmerGrid 1 2 tmp.grd -out {} -autoclean".format(out_fp)
    cmd  += " && rm tmp.grd"
    result = subprocess.run([cmd],
                            shell=True,
                            capture_output=True,
                            text=True
                            )
    os.chdir(pwd)

    if "ElmerGrid: command not found" in result.stderr:
        raise OSError('Elmer/Ice executables are not available. Make sure you are '\
                      'in the correct computing environment')

    return result

# if __name__ == '__main__':

# read in reference topography
ref_zs = np.loadtxt('topography/crmpt12_surf.dat')[:,1]
ref_zb = np.loadtxt('topography/crmpt12_bed.dat' )[:,1]
ref_xx = np.loadtxt('topography/crmpt12_bed.dat' )[:,0]

# initialize our simple class for clipping
Crompton_12 = domain(ref_xx, ref_zs, ref_zb)


fig, ax = plt.subplots(1,1, figsize=(4,3), sharex=True, sharey=True, constrained_layout=True)

# plot the reference
ax.plot(Crompton_12.xx, Crompton_12.zs, lw=1.0)
ax.plot(Crompton_12.xx, Crompton_12.zb, lw=1.0, label = 'reference $z_{\\rm b}$')

for i, params in enumerate([(2600.0, 9800.0), (2500.0, 9000.0)]):
    zs_max, xx_conect = params

    xx_new, zs_new, zb_new = Crompton_12.clip(zs_max, xx_conect, out_fp='./topography')

    label = '$z_{\\rm b}$ w/ $z_{\\rm s_{\\rm max}}$ = ' + f'{int(zs_max)}'
    ax.plot(xx_new, zb_new, lw=1.0, label=label)


ax.set_xlim(3e3, 11.1e3)
ax.set_ylim(1850, 2800)
ax.invert_xaxis()

ax.set_ylabel('Elevation [m a.s.l.]')
ax.set_xlabel('Distance  [m]')
plt.legend()

fig.savefig('clipped_topo.png', dpi=400, bbox_inches='tight')
# plt.show()
plt.close()

for i, params in enumerate([(2600.0, 9800.0), (2500.0, 9000.0)]):

    zs_max, xx_conect = params

    key = f'crmpt12_{int(zs_max)}_clip'
    dx  = 50
    out_fp = f'study/domain_clip/result/{key}/mesh_dx50/'

    print(out_fp)

    # make the mesh for the parameters
    make_mesh(key=key, dx=dx, out_fp=out_fp)
