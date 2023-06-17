#!/usr/bin/env python3

"""
periodic_surge.py:

Using the bash functions in the periodic_surge.sh script, run a periodic surge
sequence for a given parameter combination. initialization sequence includes
    - diagnostic simulation
    - prognostic simulation
"""


import os
import sys
import json
import argparse
import subprocess
from argparse import RawTextHelpFormatter

def main(argv):

    # list of tuples where first entry in tuple is the
    # bash (environmental) variables, and the second entry
    # is a more descriptive name, which makes the help
    # docstring be more readable
    varibales = [("KEY",            "key"              ),
                 ("dx",             "horzontal_spacing"),
                 ("TT",             "time_final"       ),
                 ("T0",             "restart_time"     ),
                 ("SP",             "surge_period"     ),
                 ("QP",             "quies_period"     ),
                 ("ST_dt",          "ST_dt"            ),
                 ("SD_dt",          "SD_dt"            ),
                 ("QT_dt",          "QT_dt"            ),
                 ("QD_dt",          "QD_dt"            ),
                 ("beta",           "slip_coef"        ),
                 ("z_lim",          "z_limit"          ),
                 ("offset",         "offset"           ),
                 ("T_ma",           "Temp_mean_annual" ),
                 ("RESTART",        "RESTART"          ),
                 ("SS_itters",      "SS_iterations"    )]

    # create a dictionary of dictionaries, to store passed values
    # and map between flag name for variables and the bash environmental
    # variable names
    vars = {bash_var : dict(flag_var=flag_var, val=None) \
            for (bash_var, flag_var) in varibales}

    #---------------------------------------------------------------------------
    # Specify command line arguments
    #---------------------------------------------------------------------------
    # parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)

    # NOTE: for all the arguments, the dictionary of variable mappings
    #       is used to access the "long name" of the variable. This way
    #       only needs to be done once during the list definitions above.
    #
    # NOTE: if the bash variable names change, then the corresponding key
    #       needs to change below!!

    parser.add_argument('-k',f'--{vars["KEY"]["flag_var"]}', type=str,
                        help = "glacier identifier")

    parser.add_argument('-dx',f'--{vars["dx"]["flag_var"]}', type=int,
                        help = "Horzontal gridcell spacing [m]. \n"\
                               "Mesh directory for $dx should already exist.")

    parser.add_argument('-TT', f'--{vars["TT"]["flag_var"]}', type=float,
                        help = "Final time [a] for the transient simulation. \n")

    parser.add_argument('-SP',f'--{vars["SP"]["flag_var"]}', type=str,
                        help = "Scalar offset to the mass balance curve. [m a-1] \n")

    parser.add_argument('-QP',f'--{vars["QP"]["flag_var"]}', type=str,
                        help = "Scalar offset to the mass balance curve. [m a-1] \n")

    parser.add_argument('-z_lim',f'--{vars["z_lim"]["flag_var"]}', type=str,
                        help = "Vertical limit of where sliding can occur [m a.s.l.] \n"\
                               "Needed to prevent the accumulation zone from"\
                               " 'sliding away', when headwall is fully temperate")

    parser.add_argument('-beta',f'--{vars["beta"]["flag_var"]}', type=str,
                        help = "Scalar offset to the mass balance curve. [m a-1] \n")

    parser.add_argument(f'-{vars["ST_dt"]["flag_var"]}', type=float,
                        help = "Thermal timestep length [a] during surging \n")

    parser.add_argument(f'-{vars["SD_dt"]["flag_var"]}', type=float,
                        help = "Dynamic timestep length [a] during surging \n")

    parser.add_argument(f'-{vars["QT_dt"]["flag_var"]}', type=float,
                        help = "Thermal timestep length [a] during quiesience \n")

    parser.add_argument(f'-{vars["QD_dt"]["flag_var"]}', type=float,
                        help = "Dynamic timestep length [a] during quiesience \n")

    parser.add_argument('-off',f'--{vars["offset"]["flag_var"]}', type=str,
                        help = "Scalar offset to the mass balance curve. [m a-1] \n")

    parser.add_argument('-T_ma',f'--{vars["T_ma"]["flag_var"]}', type=str,
                        help = "Mean annual air temp [C] at 'z_ref' found in the"\
                               " `params/ref_params.sif` file. \nReference value for"\
                               " 'z_ref' = 2193 [m a.s.l.] (mean elev. of Kaskawulsh)")

    parser.add_argument('-RESTART',f'--{vars["RESTART"]["flag_var"]}', type=str,
                        help = "Name of the restart (.result) file to use for the simulation"\
                               " should be in the results/$KEY/mesh_dx$dx folder")
    
    parser.add_argument('-T0',f'--{vars["T0"]["flag_var"]}', type=str,
                        help = "Restart time [a] from the previous simulation")
    
    parser.add_argument('-itters',f'--{vars["SS_itters"]["flag_var"]}', type=int,
                        help = "Number of (S)teady (S)tate itterations",
                        default = 10)


    args, _ = parser.parse_known_args(argv)


    # get the glacier identifer
    key = args.__getattribute__(vars['KEY']['flag_var'])

    # raise error if no key is given
    if key is None:
        raise EnvironmentError("${KEY} not set")

    # check if there is a param dictionary for the key
    elif os.path.exists(f'./params/{key}.json'):
        # if so, open and store it 
        with open(f'./params/{key}.json', "r") as f:
            json_params = json.load(f)

    # set the environmental varable, values have to be string
    os.environ['KEY'] = str(key)

    # after we've already access, remove KEY from dictionary
    vars.pop('KEY', None)


    for var in vars:
        # using the variale mapping dictionary, get the argument
        # value passed over the command line
        value = args.__getattribute__(vars[var]['flag_var'])

        # make sure all the varibales without default values are set
        if value is None:
            # check the json params dictionary for a default values
            if var in json_params:
                value = json_params[var]

        #
        if value is None:
            raise EnvironmentError(f"\"{vars[var]['flag_var']}\" not set")

        # set the environmental varable, values have to be string
        os.environ[var] = str(value)


    # load the functions from the initialization script,
    # then w/ the environmental variables set,
    # run a single initialization simulation
    cmd   = "source periodic_surge.sh; periodic_simulation"

    # Actually run Elmer/Ice!
    result = subprocess.run([cmd],
                            shell=True,
                            executable='/bin/bash',
                            capture_output=False,
                            text=True
                            )

if __name__ == '__main__':
    main(sys.argv[1:])
