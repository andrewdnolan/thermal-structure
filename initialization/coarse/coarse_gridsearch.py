#!/usr/bin/env python3

import json
# heuristically deterimned mass balance offset for coarse gridsearch
offsets = {

    'crmpt12'   : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0
                   },

    'crmpt18-a' : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},

    'crmpt18-b' : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},

    'glc1-a'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},

    'glc1-b'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},

    'lilk-a'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},

    'lilk-b'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},

    'klun-a'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},

    'klun-b'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},

    'sprg'      : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},

    'fish'      : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},
    'klut-a'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},
    'klut-b'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},

    'twds-a'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0},

    'twds-b'    : {"MB": { 'range': [-1.5, 0.1, -0.5],
                           'value': -1.0 },
                   "dx": 50,
                   "dt": 0}
}


for key in offsets:

    params = offsets[key]

    with open(f"params/{key}.json", "w") as f:
        json.dump(params, f, indent=2)
