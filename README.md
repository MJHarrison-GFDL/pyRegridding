# pyRegridding : Python API to MOM6 regridding algorithms

SUBMODULES:
==========

- pkg/pyRemapping: python interfaces to remapping code
- pkg/MOM6: MOM6 src code



REQUIREMENTS:
============

- gcc/gfortran compiler environment, e.g.
$ conda create --name pyRegridding-dev python=3.10 numpy=2.0 gcc gfortran matplotlib xarray netcdf4



INSTRUCTIONS:
============

Regular install:


$ (git submodule init; git submodule update)
$ (cd pkg/pyRemapping;make)
$ (make;make test)
