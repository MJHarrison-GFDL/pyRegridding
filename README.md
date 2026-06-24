# pyRegridding : Python API to MOM6 regridding algorithms

SUBMODULES:
==========

- pkg/pyRemapping: python interfaces to remapping code
- pkg/MOM6: MOM6 src code



REQUIREMENTS:
============

- gcc/gfortran compiler environment, e.g.
$  conda create --name pyRegridding-dev python=3.9 numpy=1.19 gcc gfortran matplotlib



INSTRUCTIONS:
============

Regular install:


$ (git submodule init; git submodule update)
$ (cd pkg/pyRemapping;make)
$ make test
