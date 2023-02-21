#!/usr/bin/env python3

import click
from thermal import meshing 
from thermal.meshing import make_mesh as mesh_maker

@click.command()
@click.option('-k','--key', type=click.STRING,
              help = "glacier identifier")
@click.option('-dx','--horzontal_spacing', type=click.INT,
              help = "Horzontal gridcell spacing [m]")
@click.option('-Nz','--vertical_layers', type=click.INT,
              help = "Number of vertical layers  [-]")
@click.option('-o','--out_fp', type=click.Path(),
              help = "path folder where mesh.* files will be written."\
                     " path should be relative to git trunk.")
def make_mesh(key, horzontal_spacing, vertical_layers, out_fp): 

    meshing.make_mesh(key=key, 
              dx=horzontal_spacing, 
              Nz=vertical_layers, 
              out_fp=out_fp)
