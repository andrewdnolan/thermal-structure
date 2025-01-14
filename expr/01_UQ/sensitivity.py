#!/usr/bin/env python3

"""
initialize.py:

Using the bash functions in the initialize.sh script, run the initialization
sequence for a given parameter combination. initialization sequence includes
    - diagnostic simulation
    - prognostic simulation
"""


import os
import sys
import argparse
import subprocess
from argparse import RawTextHelpFormatter

def main(argv):

    # list of tuples where first entry in tuple is the
    # bash (environmental) variables, and the second entry
    # is a more descriptive name, which makes the help
    # docstring be more readable
    varibales = [("SS_itters",      "prognostic_SS_iterations"),
                 ("KEY",            "key"                     ),
                 ("dx",             "horzontal_spacing"       ),
                 ("offset",         "offset"                  ),
                 ("T_ma",           "Temp_mean_annual"        ),
                 ("dt",             "time_step"               ),
                 ("t_f",            "time_final"              ),
                 ("Dynamic_int",    "Dynamic_exec_interval"   ),
                 ("C_firn",         "C_firn"                  ), 
                 ("f_dd",           "f_dd"                    ),
                 ("w_en",           "w_en"                    ),
                 ("w_aq",           "w_aq"                    ),
                 ("IC",             "IC"                      ),
                 ]

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

    parser.add_argument('-off',f'--{vars["offset"]["flag_var"]}', type=str,
                        help = "Scalar offset to the mass balance curve. [m a-1] \n")

    parser.add_argument('-T_ma',f'--{vars["T_ma"]["flag_var"]}', type=str,
                        help = "Mean annual air temp [C] at 'z_ref' found in the"\
                               " `params/ref_params.sif` file. \nReference value for"\
                               " 'z_ref' = 2193 [m a.s.l.] (mean elev. of Kaskawulsh)")

    parser.add_argument('-dt',f'--{vars["dt"]["flag_var"]}', type=float,
                        help = "Timestep length [a] for the transient simulation",
                        default = 1.0)

    parser.add_argument('-t_f',f'--{vars["t_f"]["flag_var"]}', type=float,
                        help = "Final time [a] for the transient simulation. \n"\
                               "The number of time steps ($NT) is calculated as \n"\
                               "        `NT = t_f // dt`")

    parser.add_argument('-Dynamic_int','--Dynamic_exec_interval', type=int,
                        help = "Interval to execute the dynamic solver, relative"\
                               "to the thermodynamic solver",
                        default = 10)

    parser.add_argument('-it_prog','--prognostic_SS_iterations', type=int,
                        help = "Number of (S)teady (S)tate itterations for the \n"\
                               "prognostic simulation",
                        default = 1)

    parser.add_argument('-C_firn', type=float,
                        help = "firn densification constant",
                        default = None)

    parser.add_argument('-f_dd', type=float,
                        help = "Degree day factor for snow",
                        default = None)

    parser.add_argument('-w_en', type=float,
                        help = "Maximum englacial water content [-]. (i.e. threshold"\
                               "for instantaneous runoff within glacier ice",
                        default = None)

    parser.add_argument('-w_aq', type=float,
                        help = "Maximum water content [-]. within firn (i.e. threshold"\
                               "for instantaneous runoff within parameterized firn column",
                        default = None)

    parser.add_argument('-IC', type=float,
                        help = "Temperature [Deg. C] used the isothermal IC.",
                        default = None)

    args, _ = parser.parse_known_args(argv)

    for var in vars:
        # using the variale mapping dictionary, get the argument
        # value passed over the command line
        value = args.__getattribute__(vars[var]['flag_var'])

        # make sure all the varibales without default values are set
        if value is None:
            raise EnvironmentError(f"\"{vars[var]['flag_var']}\" not set")

        # set the environmental varable, values have to be string
        os.environ[var] = str(value)


    # load the functions from the initialization script,
    # then w/ the environmental variables set,
    # run a single initialization simulation
    cmd   = "source sensitivity.sh; test_sensitivity"

    # Actually run Elmer/Ice!
    result = subprocess.run([cmd],
                            shell=True,
                            executable='/bin/bash',
                            capture_output=False,
                            text=True
                            )

if __name__ == '__main__':
    main(sys.argv[1:])
