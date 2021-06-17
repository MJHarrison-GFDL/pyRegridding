!> Provides the ocean grid type
module MOM_grid

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

#include <MOM_memory.h>

public MOM_grid_init
! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Ocean grid type. See mom_grid for details.
type, public :: ocean_grid_type

  integer :: isc !< The start i-index of cell centers within the computational domain
  integer :: iec !< The end i-index of cell centers within the computational domain
  integer :: jsc !< The start j-index of cell centers within the computational domain
  integer :: jec !< The end j-index of cell centers within the computational domain

  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain

  integer :: isg !< The start i-index of cell centers within the global domain
  integer :: ieg !< The end i-index of cell centers within the global domain
  integer :: jsg !< The start j-index of cell centers within the global domain
  integer :: jeg !< The end j-index of cell centers within the global domain

  integer :: IscB !< The start i-index of cell vertices within the computational domain
  integer :: IecB !< The end i-index of cell vertices within the computational domain
  integer :: JscB !< The start j-index of cell vertices within the computational domain
  integer :: JecB !< The end j-index of cell vertices within the computational domain

  integer :: IsdB !< The start i-index of cell vertices within the data domain
  integer :: IedB !< The end i-index of cell vertices within the data domain
  integer :: JsdB !< The start j-index of cell vertices within the data domain
  integer :: JedB !< The end j-index of cell vertices within the data domain

  integer :: IsgB !< The start i-index of cell vertices within the global domain
  integer :: IegB !< The end i-index of cell vertices within the global domain
  integer :: JsgB !< The start j-index of cell vertices within the global domain
  integer :: JegB !< The end j-index of cell vertices within the global domain

  integer :: isd_global !< The value of isd in the global index space (decompoistion invariant).
  integer :: jsd_global !< The value of isd in the global index space (decompoistion invariant).
  integer :: idg_offset !< The offset between the corresponding global and local i-indices.
  integer :: jdg_offset !< The offset between the corresponding global and local j-indices.
  integer :: ke         !< The number of layers in the vertical.
  logical :: symmetric  !< True if symmetric memory is used.
  logical :: nonblocking_updates  !< If true, non-blocking halo updates are
                                  !! allowed.  The default is .false. (for now).
  integer :: first_direction !< An integer that indicates which direction is
                             !! to be updated first in directionally split
                             !! parts of the calculation.  This can be altered
                             !! during the course of the run via calls to
                             !! set_first_direction.

  real :: Rad_Earth = 6.378e6 !< The radius of the planet [m].
  real :: max_depth     !< The maximum depth of the ocean in depth units [Z ~> m].

  real, dimension(:,:), allocatable :: bathyT
  real, dimension(:,:), allocatable :: mask2dT
end type ocean_grid_type


type :: hor_index_type

  integer :: isc !< The start i-index of cell centers within the computational domain
  integer :: iec !< The end i-index of cell centers within the computational domain
  integer :: jsc !< The start j-index of cell centers within the computational domain
  integer :: jec !< The end j-index of cell centers within the computational domain
end type hor_index_type

contains

!> MOM_grid_init initializes the ocean grid array sizes and grid memory.
subroutine MOM_grid_init(G, param_file, US, HI, global_indexing, bathymetry_at_vel)
  type(ocean_grid_type), intent(inout) :: G          !< The horizontal grid type
  type(param_file_type), intent(in)    :: param_file !< Parameter file handle
  type(unit_scale_type), optional, pointer :: US !< A dimensional unit scaling type
  type(hor_index_type), &
                  optional, intent(in) :: HI !< A hor_index_type for array extents
  logical,        optional, intent(in) :: global_indexing !< If true use global index
                             !! values instead of having the data domain on each
                             !! processor start at 1.
  logical,        optional, intent(in) :: bathymetry_at_vel !< If true, there are
                             !! separate values for the ocean bottom depths at
                             !! velocity points.  Otherwise the effects of topography
                             !! are entirely determined from thickness points.


end subroutine MOM_grid_init


end module MOM_grid
