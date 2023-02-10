#!/usr/bin/env python3

import os
import sys
import json
import argparse
import numpy as np
from itertools import product
from dataclasses import dataclass

cmd = "./sensitivity.py -dx {dx} --key \"{key}\" -t_f {t_f} -dt {dt} -Dynamic_int {dyn_int} -off {off} -T_ma {T_ma} "\
      "-C_firn {C_firn:1.4f} -f_dd {f_dd} -w_en {w_en} -w_aq {w_aq} -IC {IC}"

def find_precision(stride):
    """ Find the floating point precision needed
    Input:
        stride (float): parameter increment (i.e., stride)
    Output:
        NDP      (int): (N)umber of (D)ecimal (P)laces
    """
    NDP = len(str(stride).split('.')[-1])
    return NDP

class sensitivity_test:

    def __init__(self, key):
        # Glacier identifer (e.g. crmpt12)
        self.key = key
        # parameter file for the specified glacier (type: dictionary)
        self.json = self.parse_json()

    def parse_json(self):
        """ Read and parse the input parameter files
        """

        # open the json file
        with open(f'params/{self.key}.json') as file:
            j = json.load(file)

        return j

    def fill_template(self):
        """ fill the template file with the approriate parameters
        Inputs:
            key      (str) --> Glacier identifer (e.g. crmpt12)
            M        (str) --> Number of  jobs in SLURM job array
            S        (str) --> Stride for jobs in SLURM job array
            WT       (str) --> Walltime allocation (hh:mm:ss)
            MEM      (str) --> Memeory  allocation (GB)
            run_type (str) --> Run type (gridsearch or linesearch)
        """

        # read the submission template
        with open(f'./run/submission.template', 'r') as f:
            template = f.read()

        # pack the json params into another dict
        params = dict(KEY      = self.key,
                      run_name = self.key,
                      M        = self.json['memory'],
                      S        = self.json['stride'],
                      WT       = self.json['walltime'],
                      MEM      = self.json['memory'])

        # fill the submission template with the set parameters and write to disk
        with open(f'./run/{key}/{key}_{run_type}.sh', 'w') as f:
            f.write(template.format(**params))

    def write_cmdfile(self, search_type):

        # open the commands file used by the job array
        with open(f'./run/{self.key}/{search_type}.commands', 'w') as f:

            # dump the command in the text file
            f.write(self.commands)

    def write_jobscript(self, run_type, stride=None):
        # read the submission template
        with open(f'./run/submission.template', 'r') as f:
            template = f.read()
        
        # number of model runs
        M = len(self.commands.split('\n'))-1

        # if the stride is not set, run all the model runs at the same time
        if not stride: 
            stride=M

        # pack the json params into another dict
        params = dict(KEY         = self.key,
                      run_name    = self.key,
                      M           = M,
                      S           = stride,
                      WT          = self.json['walltime'],
                      MEM         = self.json['memory'],
                      search_type = run_type)

        # fill the submission template with the set parameters and write to disk
        with open(f'./run/{self.key}/{run_type}_submit.sh', 'w') as f:
            f.write(template.format(**params))

    def prep_gridsearch(self):
        """
        """

        # fixed parameters during the gridsearch
        fixed = dict(key     = self.key,
                     t_f     = self.json['TT'],
                     dx      = self.json['dx'],
                     dt      = self.json['dt'],
                     dyn_int = self.json['Dynamic_int'],
                     off     = self.json['MB'],
                     T_ma    = self.json['T_ma'])
        
        # dictionary of reference parameter values
        default = {key: self.json['params'][key]['reference'] for key in self.json['params']}

        # create an empty string to start concatenating too
        commands = ""
        
        # function to unpack the parameter range info contained in the json
        unpack_params = lambda x: getattr(np, x['spacing'])(x['start'], x['stop'], x['samples'])

        # loop over the parameters we are testing 
        for key in self.json['params']: 
            # unpack th parameter values, and loop over them.
            for val in unpack_params(self.json['params'][key]): 
                
                # make a copy of the default dictionary
                params = default.copy()
                # update with the current value of the current 
                params[key] = val

                # dump the default and varied params into the string, with end of line character
                commands += cmd.format(**fixed, **params) + '\n'

        # dump the commands str as an attribute
        self.commands = commands


def parametric_sensitivity(args):

    # initialize our class
    sens_test = sensitivity_test(args.key)
    # from the passed parameter get the run commands
    sens_test.prep_gridsearch()
    # dump the run commands into a command file
    sens_test.write_cmdfile('parametric_sensitivity')
    # update the job submission script
    sens_test.write_jobscript('parametric_sensitivity', stride=args.job_s)


def main(argv):

    # prepare to parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--key', help='Glacier identifier (e.g. crmpt12)', default=None)
    parser.add_argument('-job_s', type=int, help='slurm job array stride', default=None)

    # parse the arugments from the command line
    args, _ = parser.parse_known_args(argv)

    # make sure a search type was specified
    if args.key is None:
        raise ValueError('No glacier identifier set, must use --key flag')

    # run the function to generate input file for parametric sens. tests
    parametric_sensitivity(args)

if __name__ == '__main__':
    main(sys.argv[1:])