#!/usr/bin/env python3

import os
import sys
import json
import argparse
import numpy as np
from itertools import product
from dataclasses import dataclass

# base command used to execute a single simulation
cmd = "./initialize.py -dx {dx} --key \"{key}\" -t_f {t_f} -dt {dt} -Dynamic_int {dyn_int} -off {off} -T_ma {T_ma} "

def find_precision(stride):
    """ Find the floating point precision needed

    Input:
        stride (float): parameter increment (i.e., stride)
    Output:
        NDP      (int): (N)umber of (D)ecimal (P)laces
    """
    NDP = len(str(stride).split('.')[-1])
    return NDP

class initialization:

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

    def write_jobscript(self, run_type):
        # read the submission template
        with open(f'./run/submission.template', 'r') as f:
            template = f.read()

        # pack the json params into another dict
        params = dict(KEY         = self.key,
                      run_name    = self.key,
                      M           = len(self.commands.split('\n'))-1,
                      S           = self.json['stride'],
                      WT          = self.json['walltime'],
                      MEM         = self.json['memory'],
                      search_type = run_type)

        # fill the submission template with the set parameters and write to disk
        with open(f'./run/{self.key}/{run_type}_submit.sh', 'w') as f:
            f.write(template.format(**params))

    def prep_gridsearch(self):
        """
        """

        MB_0, MB_s, MB_f = self.json['MB']['range']
        T_ma_0, T_ma_s, T_ma_f = self.json['T_ma']['range']

        # find the precision for mass balance and air temp
        MB_fmt   = find_precision(MB_s)
        T_ma_fmt = find_precision(T_ma_s)

        # create an empty string to start concatenating too
        commands = ""

        # cartesian product of the air temp and the mass balance arrays
        for T_ma, MB in product(np.arange(T_ma_0, T_ma_f+T_ma_s, T_ma_s),
                                np.arange(MB_0,    MB_f+MB_s,    MB_s)):

            # format the parameter values based on precision determined above
            T_ma = f'{T_ma:3.{T_ma_fmt}f}'
            MB   = f'{MB:3.{MB_fmt}f}'

            # pack all the parameters into a dictionary
            params = dict(key     = self.key,
                          t_f     = self.json['TT'],
                          dx      = self.json['dx'],
                          dt      = self.json['dt'],
                          dyn_int = self.json['Dynamic_int'],
                          off     = MB,
                          T_ma    = T_ma)

            # dump the parameters into the string, with end of line character
            commands += cmd.format(**params) + '\n'

        # dump the commands str as an attribute
        self.commands = commands

    def prep_linesearch(self, MB_s, job_s):

        MB_0, MB_s, MB_f       = self.json['MB']['range']

        # number of mass balance gridpoints
        MB_N = N_points(MB_0, MB_f, MB_s)

        # dump the commands str as an attribute
        self.commands = commands

def gridsearch(args):

    # initialize our class
    init = initialization(args.key)
    # from the passed parameter get the run commands
    init.prep_gridsearch()
    # dump the run commands into a command file
    init.write_cmdfile('gridsearch')
    # update the job submission script
    init.write_jobscript('gridsearch')


def linesearch(args):

    init = initialization(args.key)

    # allow for job and mass balance stride to be overwritten from .json values for linesearch

    if args.MB_s == None:
        args.MB_s = init.json['MB']['range'][1]
    if args.job_s == None:
        args.job_s = init.json['stride']

    init.prep_linesearch(MB_s=args.MB_s, job_s=args.job_s)

    init.write_cmdfile('linesearch')


def make_parsers():
    # create the top-level parser
    parser = argparse.ArgumentParser(prog='make_submisssion')
    subparsers = parser.add_subparsers(dest='search_type')

    # create linesearch subparser
    gs_parser = subparsers.add_parser('gridsearch', help='a help')
    # set execution function for gridsearch
    gs_parser.set_defaults(func=gridsearch)
    gs_parser.add_argument('--key', help='Glacier identifer (e.g. crmpt12)')

    # create linesearch subparser
    ls_parser = subparsers.add_parser('linesearch', help='a help')
    ls_parser.add_argument('--key', help='Glacier identifer (e.g. crmpt12)')
    ls_parser.add_argument('-MB_s',  type=int, help='mass balance stride')
    ls_parser.add_argument('-job_s', type=int, help='slurm job array stride')
    # set execution function for linesearch
    ls_parser.set_defaults(func=linesearch)

    return parser

def main(argv):

    # prepare to parse command line arguments
    parser = make_parsers()
    # parse the arugments from the command line
    args, _ = parser.parse_known_args(argv)

    # make sure a search type was specified
    if args.search_type is None:
        raise ValueError('Specify "gridsearch" or "linesearch"')

    # based on the search type passed, run the wrapper functions
    args.func(args)

if __name__ == '__main__':
    main(sys.argv[1:])
