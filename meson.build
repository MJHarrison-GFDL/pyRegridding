project('pyRegridding', ['c', 'fortran'],
  version: '1.0.0',
  default_options: [
    'buildtype=release',
    'warning_level=0',
    'fortran_std=legacy',
  ],
)

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

# Locate numpy and f2py include directories at configure time
incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true,
).stdout().strip()

incdir_f2py = run_command(py,
  ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
  check: true,
).stdout().strip()

numpy_inc = include_directories(incdir_numpy)
f2py_inc  = include_directories(incdir_f2py)

# Fortran flags matching the original Makefile
add_project_arguments(
  [
    '-fdefault-real-8',
    '-fno-second-underscore',
    '-fbacktrace',
    '-fno-align-commons',
    '-fbounds-check',
  ],
  language: 'fortran',
)

# Header/module include paths required by MOM6 sources
fortran_inc = include_directories(
  'pkg/MOM6/config_src/memory/dynamic_symmetric',
  'pkg/MOM6/src/framework',
)

# All supporting Fortran sources (deduplicated from Makefile SRCS)
regrid_srcs = files(
  'MOM_error_handler.F90',
  'MOM_document.F90',
  'pkg/MOM6/src/framework/MOM_string_functions.F90',
  'MOM_file_parser.F90',
  'pkg/MOM6/src/framework/MOM_unit_scaling.F90',
  'MOM_io.F90',
  'pkg/MOM6/src/core/MOM_verticalGrid.F90',
  'MOM_get_input.F90',
  'pkg/pyRemapping/remap.f90',
  'pkg/MOM6/src/ALE/polynomial_functions.F90',
  'pkg/MOM6/src/ALE/regrid_solvers.F90',
  'pkg/MOM6/src/ALE/regrid_edge_values.F90',
  'pkg/MOM6/src/ALE/regrid_consts.F90',
  'pkg/MOM6/src/ALE/P1M_functions.F90',
  'pkg/MOM6/src/ALE/PCM_functions.F90',
  'pkg/MOM6/src/ALE/P3M_functions.F90',
  'pkg/MOM6/src/ALE/PLM_functions.F90',
  'pkg/MOM6/src/ALE/PPM_functions.F90',
  'pkg/MOM6/src/ALE/PQM_functions.F90',
  'pkg/MOM6/src/ALE/regrid_interp.F90',
  'pkg/MOM6/src/ALE/MOM_remapping.F90',
  'pkg/MOM6/src/ALE/coord_zlike.F90',
  'pkg/MOM6/src/ALE/coord_sigma.F90',
  'MOM_EOS_linear.F90',
  'MOM_EOS_Wright.F90',
  'MOM_EOS.F90',
  'pkg/MOM6/src/ALE/coord_hycom.F90',
  'pkg/MOM6/src/ALE/coord_rho.F90',
  'MOM_grid.F90',
  'MOM_regridding.F90',
)

# Generate the f2py C and Fortran wrapper sources from the interface file.
# --lower ensures identifiers are lowercased (required for correct linking with
# Fortran 90 modules under gfortran).
# --build-dir writes the generated sources to the meson output directory without
# invoking the compiler, so meson can compile them as part of the extension.
regrid_f2py = custom_target('regrid_f2py',
  input: 'regrid.f90',
  output: ['regridmodule.c', 'regrid-f2pywrappers2.f90'],
  command: [
    py, '-m', 'numpy.f2py',
    '@INPUT@',
    '-m', 'regrid',
    '--lower',
    '--build-dir', '@OUTDIR@',
  ],
)

# Build the Python extension module.  All Fortran sources are compiled together
# so meson can automatically resolve Fortran USE/module dependencies.
# fortranobject.c provides _PyFortran_Type and other runtime support required
# by every f2py-generated extension.
py.extension_module('regrid',
  sources: [regrid_f2py, incdir_f2py / 'fortranobject.c', 'regrid.f90'] + regrid_srcs,
  include_directories: [numpy_inc, f2py_inc, fortran_inc],
  dependencies: [py_dep],
  install: true,
)
