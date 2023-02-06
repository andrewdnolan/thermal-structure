# Copyright (C) 2022 by Andrew Nolan <anolan@sfu.ca>

from setuptools import setup, find_packages

setup(
    name="thermal",
    version="0.0.0",
    license="GPL v3",
    description="Pre and post processing Elmer/Ice results",
    author="Andrew Nolan",
    packages=find_packages(where='thermal'),
    install_requires=[
        "click",
        "numpy",
        "scipy",
        "matplotlib",
        "xarray",
        "dask",
        "zarr",
        "distributed"
    ],
    entry_points={
        'console_scripts' : [
        'downsample.py = thermal.scripts.downsample:downsample',
        'grid_data.py  = thermal.scripts.grid_data:grid_data'
        ],
    }
)
