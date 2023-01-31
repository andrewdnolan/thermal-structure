# Copyright (C) 2022 by Andrew Nolan <anolan@sfu.ca>

from setuptools import setup, find_packages

setup(
    name="thermal",
    version="0.0.0",
    license="GPL v3",
    description="Pre and post processing Elmer/Ice results",
    author="Andrew Nolan",
    # url="https://github.com/icepack/icepack",
    packages=find_packages(where='thermal', exclude=["add_attr.py","grid_data.py", "mesh.py"]),
    install_requires=[
        "click",
        "numpy",
        "scipy",
        "matplotlib",
        "netCDF4",
        "xarray",
        "dask",
        "dask_jobqueue",
        "bokeh",
        "seaborn",
        "zarr",
        "distributed"
    ],
    entry_points="""
        [console_scripts]
        downsample=thermal.scripts.downsample:downsample
    """,
)
