#!/usr/bin/env python3

import json
# heuristically deterimned mass balance offset for coarse gridsearch
offsets = {

    'crmpt12'   : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"
                   },

    'crmpt18-a' : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},

    'crmpt18-b' : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},

    'glc1-a'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},

    'glc1-b'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},

    'lilk-a'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},

    'lilk-b'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},

    'klun-a'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},

    'klun-b'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},

    'sprg'      : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},

    'fish'      : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},
    'klut-a'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},
    'klut-b'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},

    'twds-a'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"},

    'twds-b'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0,
                   "TT": 1000,
                   "fit": "cubic_spline"}
}


for key in offsets:

    params = offsets[key]

    with open(f"params/{key}.json", "w") as f:
        json.dump(params, f, indent=2)
