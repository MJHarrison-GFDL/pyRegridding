!> Provides subroutines for quantities specific to the equation of state
module MOM_EOS

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_string_functions, only : uppercase
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private





!> A control structure for the equation of state
type, public :: EOS_type ; private
  integer :: form_of_EOS = 0 !< The equation of state to use.
  integer :: form_of_TFreeze = 0 !< The expression for the potential temperature
                             !! of the freezing point.
  logical :: EOS_quadrature  !< If true, always use the generic (quadrature)
                             !! code for the integrals of density.
  logical :: Compressible = .true. !< If true, in situ density is a function of pressure.
! The following parameters are used with the linear equation of state only.
  real :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3]
  real :: dRho_dT   !< The partial derivative of density with temperature [kg m-3 degC-1]
  real :: dRho_dS   !< The partial derivative of density with salinity [kg m-3 ppt-1]
! The following parameters are use with the linear expression for the freezing
! point only.
  real :: TFr_S0_P0 !< The freezing potential temperature at S=0, P=0 [degC]
  real :: dTFr_dS   !< The derivative of freezing point with salinity [degC ppt-1]
  real :: dTFr_dp   !< The derivative of freezing point with pressure [degC Pa-1]

! Unit conversion factors (normally used for dimensional testing but could also allow for
! change of units of arguments to functions)
  real :: m_to_Z = 1.      !< A constant that translates distances in meters to the units of depth.
  real :: kg_m3_to_R = 1.  !< A constant that translates kilograms per meter cubed to the units of density.
  real :: R_to_kg_m3 = 1.  !< A constant that translates the units of density to kilograms per meter cubed.
  real :: RL2_T2_to_Pa = 1.!< Convert pressures from R L2 T-2 to Pa.
  real :: L_T_to_m_s = 1.  !< Convert lateral velocities from L T-1 to m s-1.

!  logical :: test_EOS = .true. ! If true, test the equation of state
end type EOS_type

! The named integers that might be stored in eqn_of_state_type%form_of_EOS.
integer, parameter, public :: EOS_LINEAR = 1 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_UNESCO = 2 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_WRIGHT = 3 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_TEOS10 = 4 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_NEMO   = 5 !< A named integer specifying an equation of state

character*(10), parameter :: EOS_LINEAR_STRING = "LINEAR" !< A string for specifying the equation of state
character*(10), parameter :: EOS_UNESCO_STRING = "UNESCO" !< A string for specifying the equation of state
character*(10), parameter :: EOS_WRIGHT_STRING = "WRIGHT" !< A string for specifying the equation of state
character*(10), parameter :: EOS_TEOS10_STRING = "TEOS10" !< A string for specifying the equation of state
character*(10), parameter :: EOS_NEMO_STRING   = "NEMO"   !< A string for specifying the equation of state
character*(10), parameter :: EOS_DEFAULT = EOS_WRIGHT_STRING !< The default equation of state

integer, parameter :: TFREEZE_LINEAR = 1  !< A named integer specifying a freezing point expression
integer, parameter :: TFREEZE_MILLERO = 2 !< A named integer specifying a freezing point expression
integer, parameter :: TFREEZE_TEOS10 = 3  !< A named integer specifying a freezing point expression
character*(10), parameter :: TFREEZE_LINEAR_STRING = "LINEAR" !< A string for specifying the freezing point expression
character*(10), parameter :: TFREEZE_MILLERO_STRING = "MILLERO_78" !< A string for specifying
                                                              !! freezing point expression
character*(10), parameter :: TFREEZE_TEOS10_STRING = "TEOS10" !< A string for specifying the freezing point expression
character*(10), parameter :: TFREEZE_DEFAULT = TFREEZE_LINEAR_STRING !< The default freezing point expression

contains

subroutine extract_member_EOS(EOS, form_of_EOS, form_of_TFreeze, EOS_quadrature, Compressible, &
                              Rho_T0_S0, drho_dT, dRho_dS, TFr_S0_P0, dTFr_dS, dTFr_dp)
  type(EOS_type),    pointer     :: EOS !< Equation of state structure
  integer, optional, intent(out) :: form_of_EOS !< A coded integer indicating the equation of state to use.
  integer, optional, intent(out) :: form_of_TFreeze !< A coded integer indicating the expression for
                                       !! the potential temperature of the freezing point.
  logical, optional, intent(out) :: EOS_quadrature !< If true, always use the generic (quadrature)
                                       !! code for the integrals of density.
  logical, optional, intent(out) :: Compressible !< If true, in situ density is a function of pressure.
  real   , optional, intent(out) :: Rho_T0_S0 !< Density at T=0 degC and S=0 ppt [kg m-3]
  real   , optional, intent(out) :: drho_dT   !< Partial derivative of density with temperature
                                              !! in [kg m-3 degC-1]
  real   , optional, intent(out) :: dRho_dS   !< Partial derivative of density with salinity
                                              !! in [kg m-3 ppt-1]
  real   , optional, intent(out) :: TFr_S0_P0 !< The freezing potential temperature at S=0, P=0 [degC]
  real   , optional, intent(out) :: dTFr_dS   !< The derivative of freezing point with salinity
                                              !! [degC PSU-1]
  real   , optional, intent(out) :: dTFr_dp   !< The derivative of freezing point with pressure
                                              !! [degC Pa-1]

  if (present(form_of_EOS    ))  form_of_EOS     = EOS%form_of_EOS
  if (present(form_of_TFreeze))  form_of_TFreeze = EOS%form_of_TFreeze
  if (present(EOS_quadrature ))  EOS_quadrature  = EOS%EOS_quadrature
  if (present(Compressible   ))  Compressible    = EOS%Compressible
  if (present(Rho_T0_S0      ))  Rho_T0_S0       = EOS%Rho_T0_S0
  if (present(drho_dT        ))  drho_dT         = EOS%drho_dT
  if (present(dRho_dS        ))  dRho_dS         = EOS%dRho_dS
  if (present(TFr_S0_P0      ))  TFr_S0_P0       = EOS%TFr_S0_P0
  if (present(dTFr_dS        ))  dTFr_dS         = EOS%dTFr_dS
  if (present(dTFr_dp        ))  dTFr_dp         = EOS%dTFr_dp

end subroutine extract_member_EOS

end module MOM_EOS

!> \namespace mom_eos
!!
!! The MOM_EOS module is a wrapper for various equations of state (e.g. Linear,
!! Wright, UNESCO) and provides a uniform interface to the rest of the model
!! independent of which equation of state is being used.
