import numpy as np 
from os import path 
import matplotlib.pyplot as plt 
import matplotlib.colors as colors

def parse_params(run_name): 
    info = run_name.split('1aTST_')[-1]
    
    param = '_'.join(info.split('_')[:-1])
    value = info.split('_')[-1]
    
    return param, value

def preprocess(src): 
    # filename 
    fn = src.encoding['source']
    # remove the file extension 
    bn = path.splitext(fn)[0]
    
    # parse the parameter values
    param, value = parse_params(bn)
    
    # expand the dimmension along parameter values to concatentate along
    src = src.expand_dims(param).assign_coords({param: (param, [float(value)])})
    
    return src


class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        # Note also that we must extrapolate beyond vmin/vmax
        center = np.abs(self.vcenter / (self.vmax-self.vmin))
        x, y = [self.vmin, self.vcenter, self.vmax], [0, center, 1.]
        return np.ma.masked_array(np.interp(value, x, y,
                                            left=-np.inf, right=np.inf))

    def inverse(self, value):
        y, x = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.interp(value, x, y, left=-np.inf, right=np.inf)

def make_colorbar(array, ref, IC=False): 
    
    pmin = array.min()
    pmax = array.max()
    pmid = ref
    pN   = array.size
    
    if not IC: 
        norm  = colors.TwoSlopeNorm(vmin=pmin, vcenter=pmid, vmax=pmax)
        cmap  = plt.matplotlib.cm.RdBu_r

    else: 
        norm  = colors.Normalize(vmin=pmin, vmax=pmax)
        cmap  = plt.matplotlib.cm.Blues_r
        Blues_r = cmap.resampled(256)
        updated = Blues_r(np.linspace(0,1,256))
        black   = np.array([0, 0, 0, 1])
        updated[-23:, :] = black
        cmap = colors.ListedColormap(updated)
        
    s_map = plt.matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    s_map.set_array(array)
    
    dp    = np.diff(array)
    left  = np.round(array[:-1] - dp/2, 2) 
    right = np.round(array[1:]  + dp/2, 2)

    bounds = np.union1d(left, right)

    return cmap, norm, s_map, bounds
