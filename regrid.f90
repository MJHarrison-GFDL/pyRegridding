module regrid_mod

  use MOM_verticalGrid, only : verticalGrid_type
  use MOM_regridding, only : regridding_CS
  use MOM_regridding, only : initialize_regridding, end_regridding, regridding_main
  use MOM_regridding, only : setCoordinateResolution, set_regrid_params, get_regrid_size
  use MOM_regridding, only : set_regrid_max_depths, set_regrid_max_thickness, getCoordinateResolution
  use MOM_regridding, only : getCoordinateInterfaces, get_zlike_CS, get_sigma_CS!, get_rho_CS
  use MOM_regridding, only : uniformResolution, DEFAULT_COORDINATE_MODE
  use MOM_regridding, only : thermo_var_ptrs
  use MOM_remapping, only : initialize_remapping, remapping_CS
  use MOM_eos, only : EOS_type, EOS_init
  use MOM_grid, only : ocean_grid_type
  use MOM_error_handler, only : MOM_error, FATAL
  use MOM_unit_scaling,     only : unit_scale_type, unit_scaling_init
  use MOM_file_parser, only : param_file_type, get_param
  use MOM_get_input, only : directories, get_MOM_input


  implicit none

  private

  public :: update_grid


contains

  function update_grid(zi,T,S,zbot,ps,frac_shelf_h,coord_mode,coord_resolution,remapping_scheme)
    real(kind=8), dimension(:,:,:), intent(in) :: zi !< The original interface positions
    real(kind=8), target, dimension(:,:,:), intent(in) :: T  !< Temperature on original grid (degC)
    real(kind=8), target, dimension(:,:,:), intent(in) :: S  !< Salinity on original grid (g kg-1)
    real(kind=8), dimension(:,:), intent(in) :: zbot !< The original interface positions
    real(kind=8), target, dimension(:,:), intent(in) :: ps  !< Surface pressure force (kg m s-2 )
    real(kind=8), dimension(:,:), intent(in) :: frac_shelf_h  !< Fractional ice shelf area (nondim)
    character(len=*), intent(in) :: coord_mode !< The coordinate mode ('REGRIDDING_ZSTAR', etc.')
    real(kind=8), dimension(:) :: coord_resolution !< The nominal coordinate resolution in units defined by coord_mode
    character(len=*), intent(in) :: remapping_scheme !< The remapping scheme to use for coordinate construction ('PLM','PPM_IH4',etc)
    real(kind=8), dimension(size(zi,1),size(zi,2),size(zi,3)) :: update_grid !< The resulting grid interface positions

    real, parameter :: epsln=1.e-10
    real, parameter :: min_thickness=1.e-9
    real :: max_depth
    type(regridding_CS) :: CS
    type(verticalGrid_type) :: GV
    type(ocean_grid_type) :: G
    type(unit_scale_type), pointer :: US=>NULL()

    type(param_file_type) :: PF
    type(directories) :: dirs

    type(thermo_var_ptrs) :: tv
    type(remapping_CS) :: remapCS

    real(kind=8), dimension(size(zi,1),size(zi,2),size(zi,3)-1) :: h0
    real(kind=8), dimension(size(zi,1),size(zi,2),size(zi,3)-1) :: h_new
    real(kind=8), dimension(size(zi,1),size(zi,2),size(zi,3)) :: dzInterface


    integer :: nk,ni,nj,i,j,k,imethod, degree,nk2, n1, n2
    integer :: ierr
    logical :: conv_adjust = .false.

    ni=size(zi,1);nj=size(zi,2);nk=size(zi,3)-1

    if (ni /= size(T,1) .or. nj .ne. size(T,2) .or. nk /= size(T,3)) call MOM_error(FATAL,'size mismatch zi/T')
    if (ni /= size(S,1) .or. nj .ne. size(S,2) .or. nk /= size(S,3)) call MOM_error(FATAL,'size mismatch zi/S')
    if (ni /= size(ps,1) .or. nj .ne. size(ps,2)) call MOM_error(FATAL,'size mismatch zi/ps')

    G%isc=1;G%iec=ni;G%jsc=1;G%jec=nj
    G%isd=1;G%ied=ni;G%jsd=1;G%jed=nj
    allocate(G%bathyT(G%isc:G%iec,G%jsc:G%jec))
    allocate(G%mask2dT(G%isc:G%iec,G%jsc:G%jec))
    G%bathyT(:,:)=zbot(:,:)

    G%mask2dT(:,:)=0.0
    do j=1,nj; do i=1,ni
      if (G%bathyT(i,j)>0.) G%mask2dT(i,j)=1
    enddo; enddo

    max_depth=-minval(zi)
    print *,'max_depth=',max_depth
    call get_MOM_input(PF, dirs)
    call unit_scaling_init(PF,US)
    GV%g_Earth=9.8
    GV%Boussinesq = .true.
    GV%Angstrom_m = 1.e-12
    GV%H_to_m = 1.0
!    GV%mks_g_Earth = US%L_T_to_m_s**2*US%m_to_Z * GV%g_Earth
    GV%ke = nk

    if (GV%Boussinesq) then
       GV%H_to_kg_m2 = GV%H_to_m
       GV%kg_m2_to_H = 1.0 / GV%H_to_kg_m2
       GV%m_to_H = 1.0 / GV%H_to_m
       GV%Angstrom_H = GV%m_to_H * GV%Angstrom_m
       GV%H_to_MKS = GV%H_to_m
    else
       GV%kg_m2_to_H = 1.0 / GV%H_to_kg_m2
       GV%m_to_H = GV%kg_m2_to_H
       GV%H_to_m = GV%H_to_kg_m2 / (GV%Rho0)
       GV%Angstrom_H = GV%Angstrom_m*1000.0*GV%kg_m2_to_H
       GV%H_to_MKS = GV%H_to_kg_m2
    endif
    GV%H_subroundoff = 1e-20 * max(GV%Angstrom_H,GV%m_to_H*1e-17)
    GV%H_to_Pa = GV%g_Earth * GV%H_to_kg_m2

    GV%H_to_Z = GV%H_to_m
    GV%Z_to_H = GV%m_to_H
    GV%Angstrom_Z = GV%Angstrom_m
    GV%H_to_RZ = GV%H_to_kg_m2
    GV%RZ_to_H = GV%kg_m2_to_H

    allocate( GV%sInterface(nk+1) )
    allocate( GV%sLayer(nk) )
    allocate( GV%g_prime(nk+1) ) ; GV%g_prime(:) = 0.0
    allocate( GV%Rlay(nk) )      ; GV%Rlay(:) = 0.0

    !call MOM_domains_init(G%domain, PF, symmetric=.true.,domain_name='MOM')
    call EOS_init(PF, tv%eqn_of_state,US)
    call set_regrid_params(CS)
    call initialize_remapping(remapCS,remapping_scheme)

    call initialize_regridding(CS, GV, US, max_depth, PF, 'MOM', coord_mode )
!    CS%coordinateResolution(:)=coord_resolution(:)


    tv%T => T
    tv%S => S

    do j=1,nj
      do i=1,ni
        do k=1,nk
          h0(i,j,k)=zi(i,j,k)-zi(i,j,k+1)
        enddo
      enddo
    enddo
    call regridding_main(remapCS, CS, G, GV, h0, tv, h_new, dzInterface, frac_shelf_h, conv_adjust)


    update_grid=zi+dzInterface



  end function update_grid



end module regrid_mod
