!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.

!*******************************************************************************
module time_average
!*******************************************************************************
use types, only : rprec
use param, only : nx, ny, nz, lbz
#ifdef PPCGNS
use cgns
#endif

private
public :: tavg_t

real(rprec), allocatable, dimension(:,:,:) :: w_uv, u_w, v_w
real(rprec), allocatable, dimension(:,:,:) :: vortx, vorty, vortz
real(rprec), allocatable, dimension(:,:,:) :: pres_real
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
real(rprec), allocatable, dimension(:,:,:) :: fza_uv
#endif
#ifdef PPSCALARS
real(rprec), allocatable, dimension(:,:,:)::theta_w
#endif
#ifdef PPBUDGET
real(rprec), allocatable, dimension(:,:,:)::uw_uv, vw_uv, w2_uv,txz_uv,tyz_uv
real(rprec), allocatable, dimension(:,:,:)::dudz_uv, dvdz_uv, dwdx_uv, dwdy_uv
#endif
!  Sums performed over time
type tavg_t
    real(rprec), dimension(:,:,:), allocatable :: u, v, w, u_w, v_w, w_uv
    real(rprec), dimension(:,:,:), allocatable :: u2, v2, w2, uv, uw, vw
    real(rprec), dimension(:,:,:), allocatable :: txx, tyy, tzz, txy, txz, tyz
    real(rprec), dimension(:,:,:), allocatable :: p, fx, fy, fz, cs_opt2
    real(rprec), dimension(:,:,:), allocatable :: vortx, vorty, vortz
    real(rprec), dimension(:,:,:), allocatable :: p2,pu,pv,pw                       !add by Mingwei
    real(rprec) :: total_time
    ! Time between calls of tavg_compute, built by summing dt
    real(rprec) :: dt
    ! Switch for determining if time averaging has been initialized
    logical :: initialized = .false.
#ifdef PPSCALARS
    real(rprec), dimension(:,:,:), allocatable :: theta, theta_w                        !add by Mingwei
    real(rprec), dimension(:,:,:), allocatable :: theta2,thetau,thetav,thetaw,thetap    !add by Mingwei
    real(rprec), dimension(:,:,:), allocatable :: sgs_t3                                !add by Mingwei 
#endif
#ifdef PPBUDGET
    real(rprec), dimension(:,:,:), allocatable :: uw_uv,vw_uv,w2_uv,txz_uv,tyz_uv!uv-grid statistics
    real(rprec), dimension(:,:,:), allocatable :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz!for budget ppSijp,taupduipdxjp
    real(rprec), dimension(:,:,:), allocatable :: u3,v3,w3,u2v,u2w,v2u,v2w,w2u,w2v,uvw!for budget uipujpukp
    real(rprec), dimension(:,:,:), allocatable :: pS11, pS12, pS13, pS22, pS23, pS33!for budget ppSijp
    real(rprec), dimension(:,:,:), allocatable :: utau11, utau12, utau13, utau22, utau23, utau33, &
                                                  vtau11, vtau12, vtau13, vtau22, vtau23, vtau33, &
                                                  wtau11, wtau12, wtau13, wtau22, wtau23, wtau33  !for budget uptaup
    real(rprec), dimension(:,:,:), allocatable :: tau11dudx, tau11dvdx, tau11dwdx,&
                                                  tau12dudy, tau12dvdy, tau12dwdy,&
                                                  tau13dudz, tau13dvdz, tau13dwdz,&
                                                  tau21dudx, tau21dvdx, tau21dwdx,&
                                                  tau22dudy, tau22dvdy, tau22dwdy,&
                                                  tau23dudz, tau23dvdz, tau23dwdz,&
                                                  tau31dudx, tau31dvdx, tau31dwdx,&
                                                  tau32dudy, tau32dvdy, tau32dwdy,&
                                                  tau33dudz, tau33dvdz, tau33dwdz
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
    real(rprec), dimension(:,:,:), allocatable :: fxu,fyu,fzu,fxv,fyv,fzv,fxw,fyw,fzw
#endif
#endif
contains
    procedure, public :: init
    procedure, public :: compute
    procedure, public :: finalize
    procedure, public :: checkpoint
end type tavg_t

character(:), allocatable :: checkpoint_tavg_file

contains

!*******************************************************************************
subroutine init(this)
!*******************************************************************************
use messages
use string_util
use param, only : read_endian, coord, path

implicit none

class(tavg_t), intent(inout) :: this

character(:), allocatable :: ftavg_in
#ifdef PPMPI
character (*), parameter :: MPI_suffix = '.c'
#endif
character (128) :: fname
! integer :: i, j, k

logical :: exst

! Create file name
allocate(ftavg_in, source = path // 'tavg.out')
allocate(checkpoint_tavg_file, source=path // 'tavg.out')

! Allocate and initialize
allocate( this%u(nx,ny,lbz:nz) ); this%u(:,:,:) = 0._rprec
allocate( this%v(nx,ny,lbz:nz) ); this%v(:,:,:) = 0._rprec
allocate( this%w(nx,ny,lbz:nz) ); this%w(:,:,:) = 0._rprec
allocate( this%u_w(nx,ny,lbz:nz) ); this%u_w(:,:,:) = 0._rprec
allocate( this%v_w(nx,ny,lbz:nz) ); this%v_w(:,:,:) = 0._rprec
allocate( this%w_uv(nx,ny,lbz:nz) ); this%w_uv(:,:,:) = 0._rprec
allocate( this%u2(nx,ny,lbz:nz) ); this%u2(:,:,:) = 0._rprec
allocate( this%v2(nx,ny,lbz:nz) ); this%v2(:,:,:) = 0._rprec
allocate( this%w2(nx,ny,lbz:nz) ); this%w2(:,:,:) = 0._rprec
allocate( this%uv(nx,ny,lbz:nz) ); this%uv(:,:,:) = 0._rprec
allocate( this%uw(nx,ny,lbz:nz) ); this%uw(:,:,:) = 0._rprec
allocate( this%vw(nx,ny,lbz:nz) ); this%vw(:,:,:) = 0._rprec
allocate( this%txx(nx,ny,lbz:nz) ); this%txx(:,:,:) = 0._rprec
allocate( this%tyy(nx,ny,lbz:nz) ); this%tyy(:,:,:) = 0._rprec
allocate( this%tzz(nx,ny,lbz:nz) ); this%tzz(:,:,:) = 0._rprec
allocate( this%txy(nx,ny,lbz:nz) ); this%txy(:,:,:) = 0._rprec
allocate( this%txz(nx,ny,lbz:nz) ); this%txz(:,:,:) = 0._rprec
allocate( this%tyz(nx,ny,lbz:nz) ); this%tyz(:,:,:) = 0._rprec
allocate( this%p(nx,ny,lbz:nz) ); this%p(:,:,:) = 0._rprec
allocate( this%fx(nx,ny,lbz:nz) ); this%fx(:,:,:) = 0._rprec
allocate( this%fy(nx,ny,lbz:nz) ); this%fy(:,:,:) = 0._rprec
allocate( this%fz(nx,ny,lbz:nz) ); this%fz(:,:,:) = 0._rprec
allocate( this%cs_opt2(nx,ny,lbz:nz) ); this%cs_opt2(:,:,:) = 0._rprec
allocate( this%vortx(nx,ny,lbz:nz) ); this%vortx(:,:,:) = 0._rprec
allocate( this%vorty(nx,ny,lbz:nz) ); this%vorty(:,:,:) = 0._rprec
allocate( this%vortz(nx,ny,lbz:nz) ); this%vortz(:,:,:) = 0._rprec
allocate( this%p2(nx,ny,lbz:nz) ); this%p2(:,:,:) = 0._rprec
allocate( this%pu(nx,ny,lbz:nz) ); this%pu(:,:,:) = 0._rprec
allocate( this%pv(nx,ny,lbz:nz) ); this%pv(:,:,:) = 0._rprec
allocate( this%pw(nx,ny,lbz:nz) ); this%pw(:,:,:) = 0._rprec
#ifdef PPSCALARS
    allocate( this%theta(nx,ny,lbz:nz) ); this%theta(:,:,:) = 0._rprec
    allocate( this%theta_w(nx,ny,lbz:nz) ); this%theta_w(:,:,:) = 0._rprec
    allocate( this%theta2(nx,ny,lbz:nz) ); this%theta2(:,:,:) = 0._rprec
    allocate( this%thetau(nx,ny,lbz:nz) ); this%thetau(:,:,:) = 0._rprec
    allocate( this%thetav(nx,ny,lbz:nz) ); this%thetav(:,:,:) = 0._rprec
    allocate( this%thetaw(nx,ny,lbz:nz) ); this%thetaw(:,:,:) = 0._rprec
    allocate( this%thetap(nx,ny,lbz:nz) ); this%thetap(:,:,:) = 0._rprec
    allocate( this%sgs_t3(nx,ny,lbz:nz) ); this%sgs_t3(:,:,:) = 0._rprec
#endif
#ifdef PPBUDGET
    allocate( this%uw_uv(nx,ny,lbz:nz) ); this%uw_uv(:,:,:) = 0._rprec
    allocate( this%vw_uv(nx,ny,lbz:nz) ); this%vw_uv(:,:,:) = 0._rprec
    allocate( this%w2_uv(nx,ny,lbz:nz) ); this%w2_uv(:,:,:) = 0._rprec
    allocate( this%txz_uv(nx,ny,lbz:nz) ); this%txz_uv(:,:,:) = 0._rprec
    allocate( this%tyz_uv(nx,ny,lbz:nz) ); this%tyz_uv(:,:,:) = 0._rprec

    allocate( this%dudx(nx,ny,lbz:nz) ); this%dudx(:,:,:) = 0._rprec
    allocate( this%dudy(nx,ny,lbz:nz) ); this%dudy(:,:,:) = 0._rprec
    allocate( this%dudz(nx,ny,lbz:nz) ); this%dudz(:,:,:) = 0._rprec
    allocate( this%dvdx(nx,ny,lbz:nz) ); this%dvdx(:,:,:) = 0._rprec
    allocate( this%dvdy(nx,ny,lbz:nz) ); this%dvdy(:,:,:) = 0._rprec
    allocate( this%dvdz(nx,ny,lbz:nz) ); this%dvdz(:,:,:) = 0._rprec
    allocate( this%dwdx(nx,ny,lbz:nz) ); this%dwdx(:,:,:) = 0._rprec
    allocate( this%dwdy(nx,ny,lbz:nz) ); this%dwdy(:,:,:) = 0._rprec
    allocate( this%dwdz(nx,ny,lbz:nz) ); this%dwdz(:,:,:) = 0._rprec

    allocate( this%u3(nx,ny,lbz:nz) ); this%u3(:,:,:) = 0._rprec
    allocate( this%v3(nx,ny,lbz:nz) ); this%v3(:,:,:) = 0._rprec
    allocate( this%w3(nx,ny,lbz:nz) ); this%w3(:,:,:) = 0._rprec
    allocate( this%u2v(nx,ny,lbz:nz) ); this%u2v(:,:,:) = 0._rprec
    allocate( this%u2w(nx,ny,lbz:nz) ); this%u2w(:,:,:) = 0._rprec
    allocate( this%v2u(nx,ny,lbz:nz) ); this%v2u(:,:,:) = 0._rprec
    allocate( this%v2w(nx,ny,lbz:nz) ); this%v2w(:,:,:) = 0._rprec
    allocate( this%w2u(nx,ny,lbz:nz) ); this%w2u(:,:,:) = 0._rprec
    allocate( this%w2v(nx,ny,lbz:nz) ); this%w2v(:,:,:) = 0._rprec
    allocate( this%uvw(nx,ny,lbz:nz) ); this%uvw(:,:,:) = 0._rprec

    allocate( this%pS11(nx,ny,lbz:nz) ); this%pS11(:,:,:) = 0._rprec
    allocate( this%pS12(nx,ny,lbz:nz) ); this%pS12(:,:,:) = 0._rprec
    allocate( this%pS13(nx,ny,lbz:nz) ); this%pS13(:,:,:) = 0._rprec
    allocate( this%pS22(nx,ny,lbz:nz) ); this%pS22(:,:,:) = 0._rprec
    allocate( this%pS23(nx,ny,lbz:nz) ); this%pS23(:,:,:) = 0._rprec
    allocate( this%pS33(nx,ny,lbz:nz) ); this%pS33(:,:,:) = 0._rprec

    allocate( this%utau11(nx,ny,lbz:nz) ); this%utau11(:,:,:) = 0._rprec
    allocate( this%utau12(nx,ny,lbz:nz) ); this%utau12(:,:,:) = 0._rprec
    allocate( this%utau13(nx,ny,lbz:nz) ); this%utau13(:,:,:) = 0._rprec
    allocate( this%utau22(nx,ny,lbz:nz) ); this%utau22(:,:,:) = 0._rprec
    allocate( this%utau23(nx,ny,lbz:nz) ); this%utau23(:,:,:) = 0._rprec
    allocate( this%utau33(nx,ny,lbz:nz) ); this%utau33(:,:,:) = 0._rprec
    allocate( this%vtau11(nx,ny,lbz:nz) ); this%vtau11(:,:,:) = 0._rprec
    allocate( this%vtau12(nx,ny,lbz:nz) ); this%vtau12(:,:,:) = 0._rprec
    allocate( this%vtau13(nx,ny,lbz:nz) ); this%vtau13(:,:,:) = 0._rprec
    allocate( this%vtau22(nx,ny,lbz:nz) ); this%vtau22(:,:,:) = 0._rprec
    allocate( this%vtau23(nx,ny,lbz:nz) ); this%vtau23(:,:,:) = 0._rprec
    allocate( this%vtau33(nx,ny,lbz:nz) ); this%vtau33(:,:,:) = 0._rprec
    allocate( this%wtau11(nx,ny,lbz:nz) ); this%wtau11(:,:,:) = 0._rprec
    allocate( this%wtau12(nx,ny,lbz:nz) ); this%wtau12(:,:,:) = 0._rprec
    allocate( this%wtau13(nx,ny,lbz:nz) ); this%wtau13(:,:,:) = 0._rprec
    allocate( this%wtau22(nx,ny,lbz:nz) ); this%wtau22(:,:,:) = 0._rprec
    allocate( this%wtau23(nx,ny,lbz:nz) ); this%wtau23(:,:,:) = 0._rprec
    allocate( this%wtau33(nx,ny,lbz:nz) ); this%wtau33(:,:,:) = 0._rprec

    allocate( this%tau11dudx(nx,ny,lbz:nz) ); this%tau11dudx(:,:,:) = 0._rprec
    allocate( this%tau11dvdx(nx,ny,lbz:nz) ); this%tau11dvdx(:,:,:) = 0._rprec
    allocate( this%tau11dwdx(nx,ny,lbz:nz) ); this%tau11dwdx(:,:,:) = 0._rprec
    allocate( this%tau12dudy(nx,ny,lbz:nz) ); this%tau12dudy(:,:,:) = 0._rprec
    allocate( this%tau12dvdy(nx,ny,lbz:nz) ); this%tau12dvdy(:,:,:) = 0._rprec
    allocate( this%tau12dwdy(nx,ny,lbz:nz) ); this%tau12dwdy(:,:,:) = 0._rprec
    allocate( this%tau13dudz(nx,ny,lbz:nz) ); this%tau13dudz(:,:,:) = 0._rprec
    allocate( this%tau13dvdz(nx,ny,lbz:nz) ); this%tau13dvdz(:,:,:) = 0._rprec
    allocate( this%tau13dwdz(nx,ny,lbz:nz) ); this%tau13dwdz(:,:,:) = 0._rprec

    allocate( this%tau21dudx(nx,ny,lbz:nz) ); this%tau21dudx(:,:,:) = 0._rprec
    allocate( this%tau21dvdx(nx,ny,lbz:nz) ); this%tau21dvdx(:,:,:) = 0._rprec
    allocate( this%tau21dwdx(nx,ny,lbz:nz) ); this%tau21dwdx(:,:,:) = 0._rprec
    allocate( this%tau22dudy(nx,ny,lbz:nz) ); this%tau22dudy(:,:,:) = 0._rprec
    allocate( this%tau22dvdy(nx,ny,lbz:nz) ); this%tau22dvdy(:,:,:) = 0._rprec
    allocate( this%tau22dwdy(nx,ny,lbz:nz) ); this%tau22dwdy(:,:,:) = 0._rprec
    allocate( this%tau23dudz(nx,ny,lbz:nz) ); this%tau23dudz(:,:,:) = 0._rprec
    allocate( this%tau23dvdz(nx,ny,lbz:nz) ); this%tau23dvdz(:,:,:) = 0._rprec
    allocate( this%tau23dwdz(nx,ny,lbz:nz) ); this%tau23dwdz(:,:,:) = 0._rprec

    allocate( this%tau31dudx(nx,ny,lbz:nz) ); this%tau31dudx(:,:,:) = 0._rprec
    allocate( this%tau31dvdx(nx,ny,lbz:nz) ); this%tau31dvdx(:,:,:) = 0._rprec
    allocate( this%tau31dwdx(nx,ny,lbz:nz) ); this%tau31dwdx(:,:,:) = 0._rprec
    allocate( this%tau32dudy(nx,ny,lbz:nz) ); this%tau32dudy(:,:,:) = 0._rprec
    allocate( this%tau32dvdy(nx,ny,lbz:nz) ); this%tau32dvdy(:,:,:) = 0._rprec
    allocate( this%tau32dwdy(nx,ny,lbz:nz) ); this%tau32dwdy(:,:,:) = 0._rprec
    allocate( this%tau33dudz(nx,ny,lbz:nz) ); this%tau33dudz(:,:,:) = 0._rprec
    allocate( this%tau33dvdz(nx,ny,lbz:nz) ); this%tau33dvdz(:,:,:) = 0._rprec
    allocate( this%tau33dwdz(nx,ny,lbz:nz) ); this%tau33dwdz(:,:,:) = 0._rprec

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
    allocate( this%fxu(nx,ny,lbz:nz) ); this%fxu(:,:,:) = 0._rprec
    allocate( this%fyu(nx,ny,lbz:nz) ); this%fyu(:,:,:) = 0._rprec
    allocate( this%fzu(nx,ny,lbz:nz) ); this%fzu(:,:,:) = 0._rprec
    allocate( this%fxv(nx,ny,lbz:nz) ); this%fxv(:,:,:) = 0._rprec
    allocate( this%fyv(nx,ny,lbz:nz) ); this%fyv(:,:,:) = 0._rprec
    allocate( this%fzv(nx,ny,lbz:nz) ); this%fzv(:,:,:) = 0._rprec
    allocate( this%fxw(nx,ny,lbz:nz) ); this%fxw(:,:,:) = 0._rprec
    allocate( this%fyw(nx,ny,lbz:nz) ); this%fyw(:,:,:) = 0._rprec
    allocate( this%fzw(nx,ny,lbz:nz) ); this%fzw(:,:,:) = 0._rprec
#endif
#endif
fname = ftavg_in
#ifdef PPMPI
call string_concat(fname, MPI_suffix, coord)
#endif

inquire (file=fname, exist=exst)
if (.not. exst) then
    !  Nothing to read in
    if (coord == 0) then
        write(*,*) ' '
        write(*,*)'No previous time averaged data - starting from scratch.'
    end if
    this%total_time = 0._rprec
else
    open(1, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(1) this%total_time
    read(1) this%u
    read(1) this%v
    read(1) this%w
    read(1) this%u_w
    read(1) this%v_w
    read(1) this%w_uv
    read(1) this%u2
    read(1) this%v2
    read(1) this%w2
    read(1) this%uv
    read(1) this%uw
    read(1) this%vw
    read(1) this%txx
    read(1) this%tyy
    read(1) this%tzz
    read(1) this%txy
    read(1) this%txz
    read(1) this%tyz
    read(1) this%fx
    read(1) this%fy
    read(1) this%fz
    read(1) this%cs_opt2
    read(1) this%vortx
    read(1) this%vorty
    read(1) this%vortz
    read(1) this%p2
    read(1) this%pu
    read(1) this%pv
    read(1) this%pw

#ifdef PPSCALARS
    read(1) this%theta
    read(1) this%theta_w
    read(1) this%theta2
    read(1) this%thetau
    read(1) this%thetav
    read(1) this%thetaw
    read(1) this%thetap
    read(1) this%sgs_t3
#endif
#ifdef PPBUDGET
    read(1) this%uw_uv
    read(1) this%vw_uv
    read(1) this%w2_uv
    read(1) this%txz_uv
    read(1) this%tyz_uv

    read(1) this%dudx
    read(1) this%dudy
    read(1) this%dudz
    read(1) this%dvdx
    read(1) this%dvdy
    read(1) this%dvdz
    read(1) this%dwdx
    read(1) this%dwdy
    read(1) this%dwdz

    read(1) this%u3
    read(1) this%v3
    read(1) this%w3
    read(1) this%u2v
    read(1) this%u2w
    read(1) this%v2u
    read(1) this%v2w
    read(1) this%w2u
    read(1) this%w2v
    read(1) this%uvw

    read(1) this%pS11
    read(1) this%pS12
    read(1) this%pS13
    read(1) this%pS22
    read(1) this%pS23
    read(1) this%pS33

    read(1) this%utau11
    read(1) this%utau12
    read(1) this%utau13
    read(1) this%utau22
    read(1) this%utau23
    read(1) this%utau33
    read(1) this%vtau11
    read(1) this%vtau12
    read(1) this%vtau13
    read(1) this%vtau22
    read(1) this%vtau23
    read(1) this%vtau33
    read(1) this%wtau11
    read(1) this%wtau12
    read(1) this%wtau13
    read(1) this%wtau22
    read(1) this%wtau23
    read(1) this%wtau33

    read(1) this%tau11dudx
    read(1) this%tau11dvdx
    read(1) this%tau11dwdx
    read(1) this%tau12dudy
    read(1) this%tau12dvdy
    read(1) this%tau12dwdy
    read(1) this%tau13dudz
    read(1) this%tau13dvdz
    read(1) this%tau13dwdz

    read(1) this%tau21dudx
    read(1) this%tau21dvdx
    read(1) this%tau21dwdx
    read(1) this%tau22dudy
    read(1) this%tau22dvdy
    read(1) this%tau22dwdy
    read(1) this%tau23dudz
    read(1) this%tau23dvdz
    read(1) this%tau23dwdz

    read(1) this%tau31dudx
    read(1) this%tau31dvdx
    read(1) this%tau31dwdx
    read(1) this%tau32dudy
    read(1) this%tau32dvdy
    read(1) this%tau32dwdy
    read(1) this%tau33dudz
    read(1) this%tau33dvdz
    read(1) this%tau33dwdz

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
    read(1) this%fxu
    read(1) this%fyu
    read(1) this%fzu
    read(1) this%fxv
    read(1) this%fyv
    read(1) this%fzv
    read(1) this%fxw
    read(1) this%fyw
    read(1) this%fzw
#endif
#endif
    close(1)
end if

! Allocate computation arrays
allocate(w_uv(nx,ny,lbz:nz), u_w(nx,ny,lbz:nz), v_w(nx,ny,lbz:nz))
allocate(vortx(nx,ny,lbz:nz), vorty(nx,ny,lbz:nz), vortz(nx,ny,lbz:nz))
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
allocate(fza_uv(nx,ny,lbz:nz))
#endif
#ifdef PPSCALARS
    allocate(theta_w(nx,ny,lbz:nz))
#endif
#ifdef PPBUDGET
    allocate(uw_uv(nx,ny,lbz:nz), vw_uv(nx,ny,lbz:nz), w2_uv(nx,ny,lbz:nz),txz_uv(nx,ny,lbz:nz),tyz_uv(nx,ny,lbz:nz))
    allocate(dudz_uv(nx,ny,lbz:nz),dvdz_uv(nx,ny,lbz:nz),dwdx_uv(nx,ny,lbz:nz),dwdy_uv(nx,ny,lbz:nz))
#endif



! Initialize dt
this%dt = 0._rprec

! Set global switch that tavg as been initialized
this%initialized = .true.

end subroutine init

!*******************************************************************************
subroutine compute(this)
!*******************************************************************************
!
!  This subroutine collects the stats for each flow
!  variable quantity
!
use param, only : ubc_mom, lbc_mom, coord, nproc
use sgs_param, only : Cs_opt2
use sim_param, only : u, v, w, p, u_f, v_f, w_f
use sim_param, only : txx, txy, tyy, txz, tyz, tzz
use sim_param, only : dudy, dudz, dvdx, dvdz, dwdx, dwdy
#ifdef PPBUDGET
use sim_param, only : dudx, dvdy, dwdz
#endif
#ifdef PPSCALARS
use scalars, only : theta, pi_z
#endif

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
use sim_param, only : fxa, fya, fza
#endif
use functions, only : interp_to_uv_grid, interp_to_w_grid
implicit none

class(tavg_t), intent(inout) :: this

! integer :: i, j, k
! real(rprec) :: u_p, u_p2, v_p, v_p2, w_p, w_p2

! Interpolation onto other grids
w_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz )
u_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(u(1:nx,1:ny,lbz:nz), lbz )
v_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(v(1:nx,1:ny,lbz:nz), lbz )
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
fza_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(fza(1:nx,1:ny,lbz:nz), lbz )
#endif
#ifdef PPSCALARS
    theta_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(theta(1:nx,1:ny,lbz:nz), lbz )
#endif
#ifdef PPBUDGET
    uw_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(u_w(1:nx,1:ny,lbz:nz)*w(1:nx,1:ny,lbz:nz), lbz )
    vw_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(v_w(1:nx,1:ny,lbz:nz)*w(1:nx,1:ny,lbz:nz), lbz )
    w2_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz)*w(1:nx,1:ny,lbz:nz), lbz )
    txz_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(txz(1:nx,1:ny,lbz:nz), lbz )
    tyz_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(tyz(1:nx,1:ny,lbz:nz), lbz )
    dudz_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(dudz(1:nx,1:ny,lbz:nz), lbz )
    dvdz_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(dvdz(1:nx,1:ny,lbz:nz), lbz )
    dwdx_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(dwdx(1:nx,1:ny,lbz:nz), lbz )
    dwdy_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(dwdy(1:nx,1:ny,lbz:nz), lbz )
#endif

! Compute vorticity
! Use vortx as an intermediate step for performing uv-w interpolation
! Vorticity is written in w grid
vortz(1:nx,1:ny,lbz:nz) = dvdx(1:nx,1:ny,lbz:nz) - dudy(1:nx,1:ny,lbz:nz) !! w grid
vortz(1:nx,1:ny,lbz:nz) = interp_to_w_grid( vortx(1:nx,1:ny,lbz:nz), lbz) !! w grid
vortx(1:nx,1:ny,lbz:nz) = dwdy(1:nx,1:ny,lbz:nz) - dvdz(1:nx,1:ny,lbz:nz) !! w grid
vorty(1:nx,1:ny,lbz:nz) = dudz(1:nx,1:ny,lbz:nz) - dwdx(1:nx,1:ny,lbz:nz) !! w grid

if (coord == 0) then
    vortz(1:nx,1:ny, 1) = 0._rprec
end if

! note: u_w not necessarily zero on walls, but only mult by w=0 vu u'w', so OK
! can zero u_w at BC anyway:
if(coord==0       .and. lbc_mom>0) u_w(:,:,1)  = 0._rprec
if(coord==nproc-1 .and. ubc_mom>0) u_w(:,:,nz) = 0._rprec
if(coord==0       .and. lbc_mom>0) v_w(:,:,1)  = 0._rprec
if(coord==nproc-1 .and. ubc_mom>0) v_w(:,:,nz) = 0._rprec

this%u(:,:,:) = this%u(:,:,:) + u(1:nx,1:ny,:)*this%dt          !! uv grid
this%v(:,:,:) = this%v(:,:,:) + v(1:nx,1:ny,:)*this%dt          !! uv grid
this%w(:,:,:) = this%w(:,:,:) + w(1:nx,1:ny,:)*this%dt          !! w grid
this%w_uv(:,:,:) = this%w(:,:,:) + w_uv(1:nx,1:ny,:)*this%dt    !! uv grid

this%u2(:,:,:) = this%u2(:,:,:) + u(1:nx,1:ny,:)*u(1:nx,1:ny,:)*this%dt     !! uv grid
this%v2(:,:,:) = this%v2(:,:,:) + v(1:nx,1:ny,:)*v(1:nx,1:ny,:)*this%dt     !! uv grid
this%w2(:,:,:) = this%w2(:,:,:) + w(1:nx,1:ny,:)*w(1:nx,1:ny,:)*this%dt     !! w grid
this%uv(:,:,:) = this%uv(:,:,:) + u(1:nx,1:ny,:)*v(1:nx,1:ny,:)*this%dt     !! uv grid
this%uw(:,:,:) = this%uw(:,:,:) + u_w(1:nx,1:ny,:)*w(1:nx,1:ny,:)*this%dt   !! w grid
this%vw(:,:,:) = this%vw(:,:,:) + v_w(1:nx,1:ny,:)*w(1:nx,1:ny,:)*this%dt   !! w grid

this%txx(:,:,:) = this%txx(:,:,:) + txx(1:nx,1:ny,:)*this%dt    !! uv grid
this%tyy(:,:,:) = this%tyy(:,:,:) + tyy(1:nx,1:ny,:)*this%dt    !! uv grid
this%tzz(:,:,:) = this%tzz(:,:,:) + tzz(1:nx,1:ny,:)*this%dt    !! uv grid
this%txy(:,:,:) = this%txy(:,:,:) + txy(1:nx,1:ny,:)*this%dt    !! uv grid
this%txz(:,:,:) = this%txz(:,:,:) + txz(1:nx,1:ny,:)*this%dt    !! w grid
this%tyz(:,:,:) = this%tyz(:,:,:) + tyz(1:nx,1:ny,:)*this%dt    !! w grid

this%p(:,:,:) = this%p(:,:,:) + p(1:nx,1:ny,:)*this%dt

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
this%fx(:,:,1:) = this%fx(:,:,1:) + fxa(1:nx,1:ny,1:)*this%dt !!uv grid
this%fy(:,:,1:) = this%fy(:,:,1:) + fya(1:nx,1:ny,1:)*this%dt !!uv grid
this%fz(:,:,1:) = this%fz(:,:,1:) + fza_uv(1:nx,1:ny,1:)*this%dt !!uv grid
#endif

this%cs_opt2(:,:,1:) = this%cs_opt2(:,:,1:) + Cs_opt2(1:nx,1:ny,1:)*this%dt

this%vortx(:,:,:) = this%vortx(:,:,:) + vortx(1:nx,1:ny,:) * this%dt !! w grid
this%vorty(:,:,:) = this%vorty(:,:,:) + vorty(1:nx,1:ny,:) * this%dt !! w grid
this%vortz(:,:,:) = this%vortz(:,:,:) + vortz(1:nx,1:ny,:) * this%dt !! w grid

this%p2(:,:,:) = this%p2(:,:,:) + p(1:nx,1:ny,:) * p(1:nx,1:ny,:) * this%dt !!uv grid
this%pu(:,:,:) = this%pu(:,:,:) + p(1:nx,1:ny,:) * u_f(1:nx,1:ny,:) * this%dt !!uv grid
this%pv(:,:,:) = this%pv(:,:,:) + p(1:nx,1:ny,:) * v_f(1:nx,1:ny,:) * this%dt !!uv grid
this%pw(:,:,:) = this%pw(:,:,:) + p(1:nx,1:ny,:) * interp_to_uv_grid(w_f(1:nx,1:ny,lbz:nz), lbz ) * this%dt !!uv grid

#ifdef PPSCALARS
    this%theta(:,:,:) = this%theta(:,:,:) + theta(1:nx,1:ny,:)*this%dt  !! uv grid
    this%theta_w(:,:,:) = this%theta_w(:,:,:) + theta_w(1:nx,1:ny,:)*this%dt  !! w grid
    this%theta2(:,:,:) = this%theta2(:,:,:) + theta(1:nx,1:ny,:)*theta(1:nx,1:ny,:)*this%dt !! uv grid
    this%thetau(:,:,:) = this%thetau(:,:,:) + theta(1:nx,1:ny,:)*u(1:nx,1:ny,:)*this%dt  !! uv grid
    this%thetav(:,:,:) = this%thetav(:,:,:) + theta(1:nx,1:ny,:)*v(1:nx,1:ny,:)*this%dt  !! uv grid
    this%thetaw(:,:,:) = this%thetaw(:,:,:) + theta_w(1:nx,1:ny,:)*w(1:nx,1:ny,:)*this%dt !! w grid
    this%thetap(:,:,:) = this%thetap(:,:,:) + theta(1:nx,1:ny,:)*p(1:nx,1:ny,:)*this%dt !! uv grid
    this%sgs_t3(:,:,:) = this%sgs_t3(:,:,:) + pi_z(1:nx,1:ny,:)*this%dt !! uv grid
#endif

#ifdef PPBUDGET
    this%uw_uv(:,:,:) = this%uw_uv(:,:,:) + u(1:nx,1:ny,:)*w_uv(1:nx,1:ny,:)*this%dt   !! uv grid
    this%vw_uv(:,:,:) = this%vw_uv(:,:,:) + v(1:nx,1:ny,:)*w_uv(1:nx,1:ny,:)*this%dt   !! uv grid
    this%w2_uv(:,:,:) = this%w2_uv(:,:,:) + w_uv(1:nx,1:ny,:)*w_uv(1:nx,1:ny,:)*this%dt     !! uv grid
    this%txz_uv(:,:,:) = this%txz_uv(:,:,:) + txz_uv(1:nx,1:ny,:)*this%dt    !! uv grid
    this%tyz_uv(:,:,:) = this%tyz_uv(:,:,:) + tyz_uv(1:nx,1:ny,:)*this%dt    !! uv grid

    this%dudx(:,:,:) = this%dudx(:,:,:) + dudx(1:nx,1:ny,:)*this%dt   !! uv grid
    this%dudy(:,:,:) = this%dudy(:,:,:) + dudy(1:nx,1:ny,:)*this%dt   !! uv grid
    this%dudz(:,:,:) = this%dudz(:,:,:) + dudz_uv(1:nx,1:ny,:)*this%dt   !! uv grid
    this%dvdx(:,:,:) = this%dvdx(:,:,:) + dvdx(1:nx,1:ny,:)*this%dt   !! uv grid
    this%dvdy(:,:,:) = this%dvdy(:,:,:) + dvdy(1:nx,1:ny,:)*this%dt   !! uv grid
    this%dvdz(:,:,:) = this%dvdz(:,:,:) + dvdz_uv(1:nx,1:ny,:)*this%dt   !! uv grid
    this%dwdx(:,:,:) = this%dwdx(:,:,:) + dwdx_uv(1:nx,1:ny,:)*this%dt   !! uv grid
    this%dwdy(:,:,:) = this%dwdy(:,:,:) + dwdy_uv(1:nx,1:ny,:)*this%dt   !! uv grid
    this%dwdz(:,:,:) = this%dwdz(:,:,:) + dwdz(1:nx,1:ny,:)*this%dt   !! uv grid

    this%u3(:,:,:) = this%u3(:,:,:) + u(1:nx,1:ny,:)*u(1:nx,1:ny,:)*u(1:nx,1:ny,:)*this%dt   !! uv grid
    this%v3(:,:,:) = this%v3(:,:,:) + v(1:nx,1:ny,:)*v(1:nx,1:ny,:)*v(1:nx,1:ny,:)*this%dt   !! uv grid
    this%w3(:,:,:) = this%w3(:,:,:) + w_uv(1:nx,1:ny,:)*w_uv(1:nx,1:ny,:)*w_uv(1:nx,1:ny,:)*this%dt   !! uv grid
    this%u2v(:,:,:) = this%u2v(:,:,:) + u(1:nx,1:ny,:)*u(1:nx,1:ny,:)*v(1:nx,1:ny,:)*this%dt   !! uv grid
    this%u2w(:,:,:) = this%u2w(:,:,:) + u(1:nx,1:ny,:)*u(1:nx,1:ny,:)*w_uv(1:nx,1:ny,:)*this%dt   !! uv grid
    this%v2u(:,:,:) = this%v2u(:,:,:) + v(1:nx,1:ny,:)*v(1:nx,1:ny,:)*u(1:nx,1:ny,:)*this%dt   !! uv grid
    this%v2w(:,:,:) = this%v2w(:,:,:) + v(1:nx,1:ny,:)*v(1:nx,1:ny,:)*w_uv(1:nx,1:ny,:)*this%dt   !! uv grid
    this%w2u(:,:,:) = this%w2u(:,:,:) + w_uv(1:nx,1:ny,:)*w_uv(1:nx,1:ny,:)*u(1:nx,1:ny,:)*this%dt   !! uv grid
    this%w2v(:,:,:) = this%w2v(:,:,:) + w_uv(1:nx,1:ny,:)*w_uv(1:nx,1:ny,:)*v(1:nx,1:ny,:)*this%dt   !! uv grid
    this%uvw(:,:,:) = this%uvw(:,:,:) + u(1:nx,1:ny,:)*v(1:nx,1:ny,:)*w_uv(1:nx,1:ny,:)*this%dt   !! uv grid

    this%pS11(:,:,:) = this%pS11(:,:,:) + 0.5_rprec*(dudx(1:nx,1:ny,:)+dudx(1:nx,1:ny,:))*p(1:nx,1:ny,:) * this%dt !! uv grid
    this%pS12(:,:,:) = this%pS12(:,:,:) + 0.5_rprec*(dudy(1:nx,1:ny,:)+dvdx(1:nx,1:ny,:))*p(1:nx,1:ny,:) * this%dt !! uv grid
    this%pS13(:,:,:) = this%pS13(:,:,:) + 0.5_rprec*(dudz_uv(1:nx,1:ny,:)+dwdx_uv(1:nx,1:ny,:))*p(1:nx,1:ny,:) * this%dt !! uv grid
    this%pS22(:,:,:) = this%pS22(:,:,:) + 0.5_rprec*(dvdy(1:nx,1:ny,:)+dvdy(1:nx,1:ny,:))*p(1:nx,1:ny,:) * this%dt !!uv grid
    this%pS23(:,:,:) = this%pS23(:,:,:) + 0.5_rprec*(dvdz_uv(1:nx,1:ny,:)+dwdy_uv(1:nx,1:ny,:))*p(1:nx,1:ny,:) * this%dt !! uv grid
    this%pS33(:,:,:) = this%pS33(:,:,:) + 0.5_rprec*(dwdz(1:nx,1:ny,:)+dwdz(1:nx,1:ny,:))*p(1:nx,1:ny,:) * this%dt !! uv grid

    this%utau11(:,:,:) = this%utau11(:,:,:) + u(1:nx,1:ny,:)*txx(1:nx,1:ny,:) * this%dt !! uv grid
    this%utau12(:,:,:) = this%utau12(:,:,:) + u(1:nx,1:ny,:)*txy(1:nx,1:ny,:) * this%dt !! uv grid
    this%utau13(:,:,:) = this%utau13(:,:,:) + u(1:nx,1:ny,:)*txz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%utau22(:,:,:) = this%utau22(:,:,:) + u(1:nx,1:ny,:)*tyy(1:nx,1:ny,:) * this%dt !! uv grid
    this%utau23(:,:,:) = this%utau23(:,:,:) + u(1:nx,1:ny,:)*tyz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%utau33(:,:,:) = this%utau33(:,:,:) + u(1:nx,1:ny,:)*tzz(1:nx,1:ny,:) * this%dt !! uv grid

    this%vtau11(:,:,:) = this%vtau11(:,:,:) + v(1:nx,1:ny,:)*txx(1:nx,1:ny,:) * this%dt !! uv grid
    this%vtau12(:,:,:) = this%vtau12(:,:,:) + v(1:nx,1:ny,:)*txy(1:nx,1:ny,:) * this%dt !! uv grid
    this%vtau13(:,:,:) = this%vtau13(:,:,:) + v(1:nx,1:ny,:)*txz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%vtau22(:,:,:) = this%vtau22(:,:,:) + v(1:nx,1:ny,:)*tyy(1:nx,1:ny,:) * this%dt !! uv grid
    this%vtau23(:,:,:) = this%vtau23(:,:,:) + v(1:nx,1:ny,:)*tyz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%vtau33(:,:,:) = this%vtau33(:,:,:) + v(1:nx,1:ny,:)*tzz(1:nx,1:ny,:) * this%dt !! uv grid

    this%wtau11(:,:,:) = this%wtau11(:,:,:) + w_uv(1:nx,1:ny,:)*txx(1:nx,1:ny,:) * this%dt !! uv grid
    this%wtau12(:,:,:) = this%wtau12(:,:,:) + w_uv(1:nx,1:ny,:)*txy(1:nx,1:ny,:) * this%dt !! uv grid
    this%wtau13(:,:,:) = this%wtau13(:,:,:) + w_uv(1:nx,1:ny,:)*txz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%wtau22(:,:,:) = this%wtau22(:,:,:) + w_uv(1:nx,1:ny,:)*tyy(1:nx,1:ny,:) * this%dt !! uv grid
    this%wtau23(:,:,:) = this%wtau23(:,:,:) + w_uv(1:nx,1:ny,:)*tyz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%wtau33(:,:,:) = this%wtau33(:,:,:) + w_uv(1:nx,1:ny,:)*tzz(1:nx,1:ny,:) * this%dt !! uv grid

    this%tau11dudx(:,:,:) = this%tau11dudx(:,:,:) + dudx(1:nx,1:ny,:)*txx(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau11dvdx(:,:,:) = this%tau11dvdx(:,:,:) + dvdx(1:nx,1:ny,:)*txx(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau11dwdx(:,:,:) = this%tau11dwdx(:,:,:) + dwdx_uv(1:nx,1:ny,:)*txx(1:nx,1:ny,:) * this%dt !! uv grid

    this%tau12dudy(:,:,:) = this%tau12dudy(:,:,:) + dudy(1:nx,1:ny,:)*txy(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau12dvdy(:,:,:) = this%tau12dvdy(:,:,:) + dvdy(1:nx,1:ny,:)*txy(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau12dwdy(:,:,:) = this%tau12dwdy(:,:,:) + dwdy_uv(1:nx,1:ny,:)*txy(1:nx,1:ny,:) * this%dt !! uv grid

    this%tau13dudz(:,:,:) = this%tau13dudz(:,:,:) + dudz_uv(1:nx,1:ny,:)*txz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau13dvdz(:,:,:) = this%tau13dvdz(:,:,:) + dvdz_uv(1:nx,1:ny,:)*txz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau13dwdz(:,:,:) = this%tau13dwdz(:,:,:) + dwdz(1:nx,1:ny,:)*txz_uv(1:nx,1:ny,:) * this%dt !! uv grid

    this%tau21dudx(:,:,:) = this%tau21dudx(:,:,:) + dudx(1:nx,1:ny,:)*txy(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau21dvdx(:,:,:) = this%tau21dvdx(:,:,:) + dvdx(1:nx,1:ny,:)*txy(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau21dwdx(:,:,:) = this%tau21dwdx(:,:,:) + dwdx_uv(1:nx,1:ny,:)*txy(1:nx,1:ny,:) * this%dt !! uv grid

    this%tau22dudy(:,:,:) = this%tau22dudy(:,:,:) + dudy(1:nx,1:ny,:)*tyy(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau22dvdy(:,:,:) = this%tau22dvdy(:,:,:) + dvdy(1:nx,1:ny,:)*tyy(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau22dwdy(:,:,:) = this%tau22dwdy(:,:,:) + dwdy_uv(1:nx,1:ny,:)*tyy(1:nx,1:ny,:) * this%dt !! uv grid

    this%tau23dudz(:,:,:) = this%tau23dudz(:,:,:) + dudz_uv(1:nx,1:ny,:)*tyz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau23dvdz(:,:,:) = this%tau23dvdz(:,:,:) + dvdz_uv(1:nx,1:ny,:)*tyz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau23dwdz(:,:,:) = this%tau23dwdz(:,:,:) + dwdz(1:nx,1:ny,:)*tyz_uv(1:nx,1:ny,:) * this%dt !! uv grid

    this%tau31dudx(:,:,:) = this%tau31dudx(:,:,:) + dudx(1:nx,1:ny,:)*txz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau31dvdx(:,:,:) = this%tau31dvdx(:,:,:) + dvdx(1:nx,1:ny,:)*txz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau31dwdx(:,:,:) = this%tau31dwdx(:,:,:) + dwdx_uv(1:nx,1:ny,:)*txz_uv(1:nx,1:ny,:) * this%dt !! uv grid

    this%tau32dudy(:,:,:) = this%tau32dudy(:,:,:) + dudy(1:nx,1:ny,:)*tyz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau32dvdy(:,:,:) = this%tau32dvdy(:,:,:) + dvdy(1:nx,1:ny,:)*tyz_uv(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau32dwdy(:,:,:) = this%tau32dwdy(:,:,:) + dwdy_uv(1:nx,1:ny,:)*tyz_uv(1:nx,1:ny,:) * this%dt !! uv grid

    this%tau33dudz(:,:,:) = this%tau33dudz(:,:,:) + dudz_uv(1:nx,1:ny,:)*tzz(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau33dvdz(:,:,:) = this%tau33dvdz(:,:,:) + dvdz_uv(1:nx,1:ny,:)*tzz(1:nx,1:ny,:) * this%dt !! uv grid
    this%tau33dwdz(:,:,:) = this%tau33dwdz(:,:,:) + dwdz(1:nx,1:ny,:)*tzz(1:nx,1:ny,:) * this%dt !! uv grid

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
    this%fxu(:,:,:) = this%fxu(:,:,:) + fxa(1:nx,1:ny,:) * u(1:nx,1:ny,:) * this%dt
    this%fyu(:,:,:) = this%fyu(:,:,:) + fya(1:nx,1:ny,:) * u(1:nx,1:ny,:) * this%dt
    this%fzu(:,:,:) = this%fzu(:,:,:) + fza_uv(1:nx,1:ny,:) * u(1:nx,1:ny,:) * this%dt

    this%fxv(:,:,:) = this%fxv(:,:,:) + fxa(1:nx,1:ny,:) * v(1:nx,1:ny,:) * this%dt
    this%fyv(:,:,:) = this%fyv(:,:,:) + fya(1:nx,1:ny,:) * v(1:nx,1:ny,:)* this%dt
    this%fzv(:,:,:) = this%fzv(:,:,:) + fza_uv(1:nx,1:ny,:) * v(1:nx,1:ny,:) * this%dt

    this%fxw(:,:,:) = this%fxv(:,:,:) + fxa(1:nx,1:ny,:) * w_uv(1:nx,1:ny,:) * this%dt
    this%fyw(:,:,:) = this%fyw(:,:,:) + fya(1:nx,1:ny,:) * w_uv(1:nx,1:ny,:) * this%dt
    this%fzw(:,:,:) = this%fzw(:,:,:) + fza_uv(1:nx,1:ny,:) * w_uv(1:nx,1:ny,:) * this%dt
#endif

#endif
!
! do k = lbz, jzmax     ! lbz = 0 for mpi runs, otherwise lbz = 1
! do j = 1, ny
! do i = 1, nx
!     u_p = u(i,j,k)       !! uv grid
!     u_p2= u_w(i,j,k)     !! w grid
!     v_p = v(i,j,k)       !! uv grid
!     v_p2= v_w(i,j,k)     !! w grid
!     w_p = w(i,j,k)       !! w grid
!     w_p2= w_uv(i,j,k)    !! uv grid

    ! this%u(i,j,k) = this%u(i,j,k) + u_p * this%dt !! uv grid
    ! this%v(i,j,k) = this%v(i,j,k) + v_p * this%dt !! uv grid
    ! this%w(i,j,k) = this%w(i,j,k) + w_p * this%dt !! w grid
    ! this%w_uv(i,j,k) = this%w_uv(i,j,k) + w_p2 * this%dt !! uv grid

    ! Note: compute u'w' on w-grid because stresses on w-grid --pj
    ! this%u2(i,j,k) = this%u2(i,j,k) + u_p * u_p * this%dt !! uv grid
    ! this%v2(i,j,k) = this%v2(i,j,k) + v_p * v_p * this%dt !! uv grid
    ! this%w2(i,j,k) = this%w2(i,j,k) + w_p * w_p * this%dt !! w grid
    ! this%uv(i,j,k) = this%uv(i,j,k) + u_p * v_p * this%dt !! uv grid
    ! this%uw(i,j,k) = this%uw(i,j,k) + u_p2 * w_p * this%dt !! w grid
    ! this%vw(i,j,k) = this%vw(i,j,k) + v_p2 * w_p * this%dt !! w grid

    ! this%txx(i,j,k) = this%txx(i,j,k) + txx(i,j,k) * this%dt !! uv grid
    ! this%tyy(i,j,k) = this%tyy(i,j,k) + tyy(i,j,k) * this%dt !! uv grid
    ! this%tzz(i,j,k) = this%tzz(i,j,k) + tzz(i,j,k) * this%dt !! uv grid
    ! this%txy(i,j,k) = this%txy(i,j,k) + txy(i,j,k) * this%dt !! uv grid
    ! this%txz(i,j,k) = this%txz(i,j,k) + txz(i,j,k) * this%dt !! w grid
    ! this%tyz(i,j,k) = this%tyz(i,j,k) + tyz(i,j,k) * this%dt !! w grid

    ! this%p(i,j,k) = this%p(i,j,k) + pres_real(i,j,k) * this%dt

! #if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
!     this%fx(i,j,k) = this%fx(i,j,k) + fxa(i,j,k) * this%dt
!     this%fy(i,j,k) = this%fy(i,j,k) + fya(i,j,k) * this%dt
!     this%fz(i,j,k) = this%fz(i,j,k) + fza_uv(i,j,k) * this%dt
! #endif
!
!     this%cs_opt2(i,j,k) = this%cs_opt2(i,j,k) + Cs_opt2(i,j,k) * this%dt

    ! this%vortx(i,j,k) = this%vortx(i,j,k) + vortx(i,j,k) * this%dt
    ! this%vorty(i,j,k) = this%vorty(i,j,k) + vorty(i,j,k) * this%dt
    ! this%vortz(i,j,k) = this%vortz(i,j,k) + vortz(i,j,k) * this%dt
! end do
! end do
! end do



! deallocate(w_uv, u_w, v_w)
! deallocate(vortx, vorty, vortz)
! #ifdef PPSCALARS
! deallocate(theta_w)
! #endif
! #ifdef defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
! deallocate(fza_uv)
! #endif
! #ifdef PPBUDGET
! deallocate(dudz_uv, dvdz_uv, dwdx_uv, dwdy_uv)
! deallocate(uw_uv, vw_uv, w2_uv, txz_uv,tyz_uv)
! #endif
! Update this%total_time for variable time stepping
this%total_time = this%total_time + this%dt

! Set this%dt back to zero for next increment
this%dt = 0._rprec

end subroutine compute

!*******************************************************************************
subroutine finalize(this)
!*******************************************************************************
use grid_m
use param, only : write_endian, lbz, path, coord, nproc
use string_util
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array,MPI_SYNC_DOWNUP
use param, only : ierr,comm
#endif
implicit none

class(tavg_t), intent(inout) :: this

#ifndef PPCGNS
character(64) :: bin_ext
#endif

character(64) :: fname_vel, fname_velw, fname_vel2, fname_tau, fname_pres
character(64) :: fname_f, fname_rs, fname_cs, fname_vort
!add by Mingwei for pressure term pp,pu,pv,pw all on uv grid
!ps2 for the pu,pv,pw,pp
!ps for the p'u',p'v',p'w',p'p'
character(64) :: fname_ps2, fname_ps
#ifdef PPSCALARS
character(64) :: fname_sca, fname_sca2, fname_sgs_t3
character(64) :: fname_rs_sca
#endif
#ifdef PPBUDGET
character(64) :: fname_ppSijp, fname_uipujpukp, fname_uptaup, fname_taupduidxjp 
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
character(64) :: fname_fipujp
#endif
#endif
! integer :: i,j,k
! Where to end with nz index.
integer :: nz_end

real(rprec), pointer, dimension(:) :: x, y, z, zw

real(rprec), allocatable, dimension(:,:,:) :: up2, vp2, wp2, upvp, upwp, vpwp
real(rprec), allocatable, dimension(:,:,:) :: pp2, ppup, ppvp, ppwp
#ifdef PPSCALARS
real(rprec), allocatable, dimension(:,:,:) :: thetap2,thetapup,thetapvp,thetapwp,thetappp
#endif
#ifdef PPBUDGET
real(rprec), allocatable, dimension(:,:,:) :: ppS11p, ppS12p, ppS13p, ppS22p, ppS23p, ppS33p
real(rprec), allocatable, dimension(:,:,:) :: upupup, vpvpvp, wpwpwp, upupvp, upupwp, vpvpup, vpvpwp, wpwpup, wpwpvp, upvpwp
real(rprec), allocatable, dimension(:,:,:) :: uptau11p, uptau12p, uptau13p, uptau22p, uptau23p, uptau33p,&
                   vptau11p, vptau12p, vptau13p, vptau22p, vptau23p, vptau33p,&
                   wptau11p, wptau12p, wptau13p, wptau22p, wptau23p, wptau33p
real(rprec), allocatable, dimension(:,:,:) :: tau11pdudxp, tau11pdvdxp, tau11pdwdxp,&
                   tau12pdudyp, tau12pdvdyp, tau12pdwdyp,&
                   tau13pdudzp, tau13pdvdzp, tau13pdwdzp,&
                   tau21pdudxp, tau21pdvdxp, tau21pdwdxp,&
                   tau22pdudyp, tau22pdvdyp, tau22pdwdyp,&
                   tau23pdudzp, tau23pdvdzp, tau23pdwdzp,&
                   tau31pdudxp, tau31pdvdxp, tau31pdwdxp,&
                   tau32pdudyp, tau32pdvdyp, tau32pdwdyp,&
                   tau33pdudzp, tau33pdvdzp, tau33pdwdzp
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
    real(rprec), allocatable, dimension(:,:,:) :: fxpup,fypup,fzpup,fxpvp,fypvp,fzpvp,fxpwp,fypwp,fzpwp
#endif
    
#endif
    !end for budget
nullify(x,y,z,zw)

x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

#ifdef PPMPI
! This adds one more element to the last processor (which contains an extra one)
! Processor nproc-1 has data from 1:nz
! Rest of processors have data from 1:nz-1
if ( coord == nproc-1 ) then
    nz_end = 0
else
    nz_end = 1
end if
#else
nz_end = 0
#endif

! Common file name
fname_vel = path // 'output/veluv_avg'
fname_velw = path // 'output/velw_avg'
fname_vel2 = path // 'output/vel2_avg'
fname_tau = path // 'output/tau_avg'
fname_f = path // 'output/force_avg'
fname_pres = path // 'output/pre_uv_avg'
fname_rs = path // 'output/rs'
fname_cs = path // 'output/cs_opt2'
fname_vort = path // 'output/vort_avg'
fname_ps2 = path // 'output/ps2_uv_avg'
fname_ps = path // 'output/ps_uv_avg'
#ifdef PPSCALARS
fname_sca = path // 'output/scalar_avg'
fname_sca2 = path // 'output/scalar2_avg'
fname_sgs_t3 = path // 'output/scalar_sgs'
fname_rs_sca = path // 'output/scalar_rs'
#endif
#ifdef PPBUDGET
fname_ppSijp = path // 'output/ppSijp'
fname_uipujpukp = path // 'output/uipujpukp'
fname_uptaup = path // 'output/uptaup'
fname_taupduidxjp = path // 'output/taupduidxjp'
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
fname_fipujp = path // 'output/fipujp'
#endif
#endif
! CGNS
#ifdef PPCGNS
call string_concat(fname_vel, '.cgns')
call string_concat(fname_velw, '.cgns')
call string_concat(fname_vel2, '.cgns')
call string_concat(fname_tau, '.cgns')
call string_concat(fname_pres, '.cgns')
call string_concat(fname_f, '.cgns')
call string_concat(fname_rs, '.cgns')
call string_concat(fname_cs, '.cgns')
call string_concat(fname_vort, '.cgns')

! Binary
#else
#ifdef PPMPI
call string_splice(bin_ext, '.c', coord, '.bin')
#else
bin_ext = '.bin'
#endif
call string_concat(fname_vel, bin_ext)
call string_concat(fname_velw, bin_ext)
call string_concat(fname_vel2, bin_ext)
call string_concat(fname_tau, bin_ext)
call string_concat(fname_pres, bin_ext)
call string_concat(fname_f, bin_ext)
call string_concat(fname_rs, bin_ext)
call string_concat(fname_cs, bin_ext)
call string_concat(fname_vort, bin_ext)
call string_concat(fname_ps2, bin_ext)
call string_concat(fname_ps, bin_ext)
#ifdef PPSCALARS
call string_concat(fname_sca, bin_ext)
call string_concat(fname_sca2, bin_ext)
call string_concat(fname_sgs_t3, bin_ext)
call string_concat(fname_rs_sca, bin_ext)
#endif
#ifdef PPBUDGET
call string_concat(fname_ppSijp, bin_ext)
call string_concat(fname_uipujpukp, bin_ext)
call string_concat(fname_uptaup, bin_ext)
call string_concat(fname_taupduidxjp, bin_ext)
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
call string_concat(fname_fipujp, bin_ext)
#endif
#endif
#endif

! Final checkpoint all restart data
call this%checkpoint()

#ifdef PPMPI
call mpi_barrier(comm, ierr)
#endif

!  Perform time averaging operation
this%u(:,:,:) = this%u(:,:,:) /  this%total_time
this%v(:,:,:) = this%v(:,:,:) /  this%total_time
this%w(:,:,:) = this%w(:,:,:) /  this%total_time
this%u_w(:,:,:)  = this%u_w(:,:,:)  /  this%total_time
this%v_w(:,:,:)  = this%v_w(:,:,:)  /  this%total_time
this%w_uv(:,:,:) = this%w_uv(:,:,:) /  this%total_time
this%u2(:,:,:) = this%u2(:,:,:) /  this%total_time
this%v2(:,:,:) = this%v2(:,:,:) /  this%total_time
this%w2(:,:,:) = this%w2(:,:,:) /  this%total_time
this%uv(:,:,:) = this%uv(:,:,:) /  this%total_time
this%uw(:,:,:) = this%uw(:,:,:) /  this%total_time
this%vw(:,:,:) = this%vw(:,:,:) /  this%total_time
this%txx(:,:,:) = this%txx(:,:,:) /  this%total_time
this%tyy(:,:,:) = this%tyy(:,:,:) /  this%total_time
this%tzz(:,:,:) = this%tzz(:,:,:) /  this%total_time
this%txy(:,:,:) = this%txy(:,:,:) /  this%total_time
this%txz(:,:,:) = this%txz(:,:,:) /  this%total_time
this%tyz(:,:,:) = this%tyz(:,:,:) /  this%total_time
this%p(:,:,:) = this%p(:,:,:) /  this%total_time
this%fx(:,:,:) = this%fx(:,:,:) /  this%total_time
this%fy(:,:,:) = this%fy(:,:,:) /  this%total_time
this%fz(:,:,:) = this%fz(:,:,:) /  this%total_time
this%cs_opt2(:,:,:) = this%cs_opt2(:,:,:) /  this%total_time
this%vortx(:,:,:) = this%vortx(:,:,:) /  this%total_time
this%vorty(:,:,:) = this%vorty(:,:,:) /  this%total_time
this%vortz(:,:,:) = this%vortz(:,:,:) /  this%total_time

this%p2(:,:,:) = this%p2(:,:,:) /  this%total_time
this%pu(:,:,:) = this%pu(:,:,:) /  this%total_time
this%pv(:,:,:) = this%pv(:,:,:) /  this%total_time
this%pw(:,:,:) = this%pw(:,:,:) /  this%total_time

#ifdef PPSCALARS
    this%theta(:,:,:) = this%theta(:,:,:) /  this%total_time
    this%theta_w(:,:,:) = this%theta_w(:,:,:) /  this%total_time
    this%theta2(:,:,:) = this%theta2(:,:,:) /  this%total_time
    this%thetau(:,:,:) = this%thetau(:,:,:) /  this%total_time
    this%thetav(:,:,:) = this%thetav(:,:,:) /  this%total_time
    this%thetaw(:,:,:) = this%thetaw(:,:,:) /  this%total_time
    this%thetap(:,:,:) = this%thetap(:,:,:) /  this%total_time
    this%sgs_t3(:,:,:) = this%sgs_t3(:,:,:) /  this%total_time
#endif

#ifdef PPBUDGET
    this%uw_uv(:,:,:) = this%uw_uv(:,:,:) /  this%total_time
    this%vw_uv(:,:,:) = this%vw_uv(:,:,:) /  this%total_time
    this%w2_uv(:,:,:) = this%w2_uv(:,:,:) /  this%total_time
    this%txz_uv(:,:,:) = this%txz_uv(:,:,:) /  this%total_time
    this%tyz_uv(:,:,:) = this%tyz_uv(:,:,:) /  this%total_time

    this%dudx(:,:,:) = this%dudx(:,:,:) /  this%total_time
    this%dudy(:,:,:) = this%dudy(:,:,:) /  this%total_time
    this%dudz(:,:,:) = this%dudz(:,:,:) /  this%total_time
    this%dvdx(:,:,:) = this%dvdx(:,:,:) /  this%total_time
    this%dvdy(:,:,:) = this%dvdy(:,:,:) /  this%total_time
    this%dvdz(:,:,:) = this%dvdz(:,:,:) /  this%total_time
    this%dwdx(:,:,:) = this%dwdx(:,:,:) /  this%total_time
    this%dwdy(:,:,:) = this%dwdy(:,:,:) /  this%total_time
    this%dwdz(:,:,:) = this%dwdz(:,:,:) /  this%total_time

    this%u3(:,:,:) = this%u3(:,:,:) /  this%total_time
    this%v3(:,:,:) = this%v3(:,:,:) /  this%total_time
    this%w3(:,:,:) = this%w3(:,:,:) /  this%total_time
    this%u2v(:,:,:) = this%u2v(:,:,:) /  this%total_time
    this%u2w(:,:,:) = this%u2w(:,:,:) /  this%total_time
    this%v2u(:,:,:) = this%v2u(:,:,:) /  this%total_time
    this%v2w(:,:,:) = this%v2w(:,:,:) /  this%total_time
    this%w2u(:,:,:) = this%w2u(:,:,:) /  this%total_time
    this%w2v(:,:,:) = this%w2v(:,:,:) /  this%total_time
    this%uvw(:,:,:) = this%uvw(:,:,:) /  this%total_time

    this%pS11(:,:,:) = this%pS11(:,:,:) /  this%total_time
    this%pS12(:,:,:) = this%pS12(:,:,:) /  this%total_time
    this%pS13(:,:,:) = this%pS13(:,:,:) /  this%total_time
    this%pS22(:,:,:) = this%pS22(:,:,:) /  this%total_time
    this%pS23(:,:,:) = this%pS23(:,:,:) /  this%total_time
    this%pS33(:,:,:) = this%pS33(:,:,:) /  this%total_time

    this%utau11(:,:,:) = this%utau11(:,:,:) /  this%total_time
    this%utau12(:,:,:) = this%utau12(:,:,:) /  this%total_time
    this%utau13(:,:,:) = this%utau13(:,:,:) /  this%total_time
    this%utau22(:,:,:) = this%utau22(:,:,:) /  this%total_time
    this%utau23(:,:,:) = this%utau23(:,:,:) /  this%total_time
    this%utau33(:,:,:) = this%utau33(:,:,:) /  this%total_time

    this%vtau11(:,:,:) = this%vtau11(:,:,:) /  this%total_time
    this%vtau12(:,:,:) = this%vtau12(:,:,:) /  this%total_time
    this%vtau13(:,:,:) = this%vtau13(:,:,:) /  this%total_time
    this%vtau22(:,:,:) = this%vtau22(:,:,:) /  this%total_time
    this%vtau23(:,:,:) = this%vtau23(:,:,:) /  this%total_time
    this%vtau33(:,:,:) = this%vtau33(:,:,:) /  this%total_time

    this%wtau11(:,:,:) = this%wtau11(:,:,:) /  this%total_time
    this%wtau12(:,:,:) = this%wtau12(:,:,:) /  this%total_time
    this%wtau13(:,:,:) = this%wtau13(:,:,:) /  this%total_time
    this%wtau22(:,:,:) = this%wtau22(:,:,:) /  this%total_time
    this%wtau23(:,:,:) = this%wtau23(:,:,:) /  this%total_time
    this%wtau33(:,:,:) = this%wtau33(:,:,:) /  this%total_time

    this%tau11dudx(:,:,:) = this%tau11dudx(:,:,:) /  this%total_time
    this%tau11dvdx(:,:,:) = this%tau11dvdx(:,:,:) /  this%total_time
    this%tau11dwdx(:,:,:) = this%tau11dwdx(:,:,:) /  this%total_time

    this%tau12dudy(:,:,:) = this%tau12dudy(:,:,:) /  this%total_time
    this%tau12dvdy(:,:,:) = this%tau12dvdy(:,:,:) /  this%total_time
    this%tau12dwdy(:,:,:) = this%tau12dwdy(:,:,:) /  this%total_time

    this%tau13dudz(:,:,:) = this%tau13dudz(:,:,:) /  this%total_time
    this%tau13dvdz(:,:,:) = this%tau13dvdz(:,:,:) /  this%total_time
    this%tau13dwdz(:,:,:) = this%tau13dwdz(:,:,:) /  this%total_time

    this%tau21dudx(:,:,:) = this%tau21dudx(:,:,:) /  this%total_time
    this%tau21dvdx(:,:,:) = this%tau21dvdx(:,:,:) /  this%total_time
    this%tau21dwdx(:,:,:) = this%tau21dwdx(:,:,:) /  this%total_time

    this%tau22dudy(:,:,:) = this%tau22dudy(:,:,:) /  this%total_time
    this%tau22dvdy(:,:,:) = this%tau22dvdy(:,:,:) /  this%total_time
    this%tau22dwdy(:,:,:) = this%tau22dwdy(:,:,:) /  this%total_time

    this%tau23dudz(:,:,:) = this%tau23dudz(:,:,:) /  this%total_time
    this%tau23dvdz(:,:,:) = this%tau23dvdz(:,:,:) /  this%total_time
    this%tau23dwdz(:,:,:) = this%tau23dwdz(:,:,:) /  this%total_time

    this%tau31dudx(:,:,:) = this%tau31dudx(:,:,:) /  this%total_time
    this%tau31dvdx(:,:,:) = this%tau31dvdx(:,:,:) /  this%total_time
    this%tau31dwdx(:,:,:) = this%tau31dwdx(:,:,:) /  this%total_time

    this%tau32dudy(:,:,:) = this%tau32dudy(:,:,:) /  this%total_time
    this%tau32dvdy(:,:,:) = this%tau32dvdy(:,:,:) /  this%total_time
    this%tau32dwdy(:,:,:) = this%tau32dwdy(:,:,:) /  this%total_time

    this%tau33dudz(:,:,:) = this%tau33dudz(:,:,:) /  this%total_time
    this%tau33dvdz(:,:,:) = this%tau33dvdz(:,:,:) /  this%total_time
    this%tau33dwdz(:,:,:) = this%tau33dwdz(:,:,:) /  this%total_time

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
    this%fxu(:,:,:) = this%fxu(:,:,:) /  this%total_time
    this%fyu(:,:,:) = this%fyu(:,:,:) /  this%total_time
    this%fzu(:,:,:) = this%fzu(:,:,:) /  this%total_time

    this%fxv(:,:,:) = this%fxv(:,:,:) /  this%total_time
    this%fyv(:,:,:) = this%fyv(:,:,:) /  this%total_time
    this%fzv(:,:,:) = this%fzv(:,:,:) /  this%total_time

    this%fxw(:,:,:) = this%fxv(:,:,:) /  this%total_time
    this%fyw(:,:,:) = this%fyw(:,:,:) /  this%total_time
    this%fzw(:,:,:) = this%fzw(:,:,:) /  this%total_time
#endif

#endif

#ifdef PPMPI
call mpi_barrier( comm, ierr )
#endif

!  Sync entire tavg structure
#ifdef PPMPI
call mpi_sync_real_array( this%u(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%v(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%w(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%u2(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%v2(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%w2(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%uw(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vw(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%uv(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%p(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%fx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%cs_opt2(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vortx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vorty(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vortz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
!add by Mingwei
call mpi_sync_real_array( this%p2(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%pu(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%pv(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%pw(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )

#ifdef PPSCALARS
call mpi_sync_real_array( this%theta(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%theta_w(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%theta2(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%thetau(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%thetav(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%thetaw(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%thetap(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%sgs_t3(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
#endif
#ifdef PPBUDGET
call mpi_sync_real_array( this%uw_uv(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vw_uv(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%w2_uv(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%txz_uv(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tyz_uv(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( this%u3(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%v3(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%w3(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%u2v(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%u2w(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%v2w(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%v2u(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%w2u(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%w2v(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%uvw(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( this%dudx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%dudy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%dudz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%dvdx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%dvdy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%dvdz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%dwdx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%dwdy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%dwdz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( this%pS11(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%pS12(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%pS13(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%pS22(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%pS23(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%pS33(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( this%utau11(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%utau12(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%utau13(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%utau22(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%utau23(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%utau33(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( this%vtau11(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vtau12(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vtau13(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vtau22(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vtau23(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vtau33(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( this%wtau11(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%wtau12(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%wtau13(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%wtau22(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%wtau23(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%wtau33(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( this%tau11dudx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau11dvdx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau11dwdx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau12dudy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau12dvdy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau12dwdy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau13dudz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau13dvdz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau13dwdz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( this%tau21dudx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau21dvdx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau21dwdx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau22dudy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau22dvdy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau22dwdy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau23dudz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau23dvdz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau23dwdz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )

call mpi_sync_real_array( this%tau31dudx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau31dvdx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau31dwdx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau32dudy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau32dvdy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau32dwdy(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau33dudz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau33dvdz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%tau33dwdz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
call mpi_sync_real_array( this%fxu(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%fyu(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%fzu(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%fxv(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%fyv(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%fzv(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%fxw(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%fyw(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%fzw(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
#endif
#endif
#endif

! Write all the 3D data
#ifdef PPCGNS
! Write CGNS Data
call write_parallel_cgns (fname_vel ,nx, ny, nz - nz_end, nz_tot,              &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 3,                                  &
    (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                               &
    (/ this%u(1:nx,1:ny,1:nz-nz_end),                                          &
       this%v(1:nx,1:ny,1:nz-nz_end),                                          &
       this%w_uv(1:nx,1:ny,1:nz-nz_end) /) )

call write_parallel_cgns (fname_velw ,nx, ny, nz - nz_end, nz_tot,             &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ),                                    &
    1, (/ 'VelocityZ' /), (/ this%w(1:nx,1:ny,1:nz-nz_end) /) )

call write_parallel_cgns(fname_vel2,nx,ny,nz- nz_end,nz_tot,                   &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 6,                                 &
    (/ 'Mean--uu', 'Mean--vv', 'Mean--ww','Mean--uw','Mean--vw','Mean--uv'/),  &
    (/ this%u2(1:nx,1:ny,1:nz-nz_end),                                         &
       this%v2(1:nx,1:ny,1:nz-nz_end),                                         &
       this%w2(1:nx,1:ny,1:nz-nz_end),                                         &
       this%uw(1:nx,1:ny,1:nz-nz_end),                                         &
       this%vw(1:nx,1:ny,1:nz-nz_end),                                         &
       this%uv(1:nx,1:ny,1:nz-nz_end) /) )

call write_parallel_cgns(fname_tau,nx,ny,nz- nz_end,nz_tot,                    &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 6,                                 &
    (/ 'Tau--txx', 'Tau--txy', 'Tau--tyy','Tau--txz','Tau--tyz','Tau--tzz'/),  &
    (/ this%txx(1:nx,1:ny,1:nz-nz_end),                                        &
       this%txy(1:nx,1:ny,1:nz-nz_end),                                        &
       this%tyy(1:nx,1:ny,1:nz-nz_end),                                        &
       this%txz(1:nx,1:ny,1:nz-nz_end),                                        &
       this%tyz(1:nx,1:ny,1:nz-nz_end),                                        &
       this%tzz(1:nx,1:ny,1:nz-nz_end) /) )

call write_parallel_cgns(fname_pres,nx,ny,nz- nz_end,nz_tot,                   &
   (/ 1, 1,   (nz-1)*coord + 1 /),                                             &
   (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                                &
   x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 1,                                  &
   (/ 'pressure' /),                                                           &
   (/ this%p(1:nx,1:ny,1:nz-nz_end) /) )

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
call write_parallel_cgns(fname_f,nx,ny,nz- nz_end,nz_tot,                      &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 3,                                 &
    (/ 'bodyForX', 'bodyForY', 'bodyForZ' /),                                  &
    (/ this%fx(1:nx,1:ny,1:nz-nz_end),                                         &
       this%fy(1:nx,1:ny,1:nz-nz_end),                                         &
       this%fz(1:nx,1:ny,1:nz-nz_end) /) )
#endif

call write_parallel_cgns(fname_cs,nx,ny,nz- nz_end,nz_tot,                     &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 1,                                 &
    (/ 'Cs_Coeff'/),  (/ this%cs_opt2(1:nx,1:ny,1:nz- nz_end) /) )

call write_parallel_cgns(fname_vort,nx,ny,nz- nz_end,nz_tot,                   &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 3,                                 &
    (/ 'VorticityX', 'VorticityY', 'VorticityZ' /),                            &
    (/ this%vortx(1:nx,1:ny,1:nz-nz_end),                                      &
       this%vorty(1:nx,1:ny,1:nz-nz_end),                                      &
       this%vortz(1:nx,1:ny,1:nz-nz_end) /)

#else
! Write binary data
open(unit=13, file=fname_vel, form='unformatted', convert=write_endian,        &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%u(:nx,:ny,1:nz)
write(13,rec=2) this%v(:nx,:ny,1:nz)
write(13,rec=3) this%w_uv(:nx,:ny,1:nz)
close(13)

! Write binary data
open(unit=13, file=fname_velw, form='unformatted', convert=write_endian,       &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%w(:nx,:ny,1:nz)
close(13)

open(unit=13, file=fname_vel2, form='unformatted', convert=write_endian,       &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%u2(:nx,:ny,1:nz)
write(13,rec=2) this%v2(:nx,:ny,1:nz)
write(13,rec=3) this%w2(:nx,:ny,1:nz)
write(13,rec=4) this%uw(:nx,:ny,1:nz)
write(13,rec=5) this%vw(:nx,:ny,1:nz)
write(13,rec=6) this%uv(:nx,:ny,1:nz)
close(13)

open(unit=13, file=fname_tau, form='unformatted', convert=write_endian,        &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%txx(:nx,:ny,1:nz)
write(13,rec=2) this%txy(:nx,:ny,1:nz)
write(13,rec=3) this%tyy(:nx,:ny,1:nz)
write(13,rec=4) this%txz(:nx,:ny,1:nz)
write(13,rec=5) this%tyz(:nx,:ny,1:nz)
write(13,rec=6) this%tzz(:nx,:ny,1:nz)
close(13)

open(unit=13, file=fname_pres, form='unformatted', convert=write_endian,       &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%p(:nx,:ny,1:nz)
close(13)

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
open(unit=13, file=fname_f, form='unformatted', convert=write_endian,          &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%fx(:nx,:ny,1:nz)
write(13,rec=2) this%fy(:nx,:ny,1:nz)
write(13,rec=3) this%fz(:nx,:ny,1:nz)
close(13)
#endif

open(unit=13, file=fname_cs, form='unformatted', convert=write_endian,         &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%cs_opt2(:nx,:ny,1:nz)
close(13)

open(unit=13, file=fname_vort, form='unformatted', convert=write_endian,       &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%vortx(:nx,:ny,1:nz)
write(13,rec=2) this%vorty(:nx,:ny,1:nz)
write(13,rec=3) this%vortz(:nx,:ny,1:nz)
close(13)

open(unit=13, file=fname_ps2, form='unformatted', convert=write_endian,          &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%p2(:nx,:ny,1:nz)
write(13,rec=2) this%pu(:nx,:ny,1:nz)
write(13,rec=3) this%pv(:nx,:ny,1:nz)
write(13,rec=4) this%pw(:nx,:ny,1:nz)
close(13)

#ifdef PPSCALARS
open(unit=13, file=fname_sca, form='unformatted', convert=write_endian,         &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%theta(:nx,:ny,1:nz)
write(13,rec=2) this%theta_w(:nx,:ny,1:nz)
close(13)

open(unit=13, file=fname_sca2, form='unformatted', convert=write_endian,          &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%theta2(:nx,:ny,1:nz)
write(13,rec=2) this%thetau(:nx,:ny,1:nz)
write(13,rec=3) this%thetav(:nx,:ny,1:nz)
write(13,rec=4) this%thetaw(:nx,:ny,1:nz)
write(13,rec=5) this%thetap(:nx,:ny,1:nz)
close(13)

open(unit=13, file=fname_sgs_t3, form='unformatted', convert=write_endian,         &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%sgs_t3(:nx,:ny,1:nz)
close(13)

#endif

#endif

#ifdef PPMPI
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
#endif

! Do the Reynolds stress calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
allocate( up2(nx,ny,lbz:nz) )
allocate( vp2(nx,ny,lbz:nz) )
allocate( wp2(nx,ny,lbz:nz) )
allocate( upvp(nx,ny,lbz:nz) )
allocate( upwp(nx,ny,lbz:nz) )
allocate( vpwp(nx,ny,lbz:nz) )
up2 = this%u2 - this%u * this%u
vp2 = this%v2 - this%v * this%v
wp2 = this%w2 - this%w * this%w
upvp = this%uv - this%u * this%v
!! using u_w and v_w below instead of u and v ensures that the Reynolds
!! stresses are on the same grid as the squared velocities (i.e., w-grid)
upwp = this%uw - this%u_w * this%w
vpwp = this%vw - this%v_w * this%w

#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_rs,nx,ny,nz- nz_end,nz_tot,                     &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                                  &
    (/ 'Meanupup', 'Meanvpvp', 'Meanwpwp','Meanupwp','Meanvpwp','Meanupvp'/),  &
    (/ up2(1:nx,1:ny,1:nz- nz_end) ,                                         &
    vp2(1:nx,1:ny,1:nz- nz_end) ,                                              &
    wp2(1:nx,1:ny,1:nz- nz_end) ,                                              &
    upwp(1:nx,1:ny,1:nz- nz_end) ,                                             &
    vpwp(1:nx,1:ny,1:nz- nz_end) ,                                             &
    upvp(1:nx,1:ny,1:nz- nz_end)  /) )
#else
! Write binary data
open(unit=13, file=fname_rs, form='unformatted', convert=write_endian,         &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) up2(:nx,:ny,1:nz)
write(13,rec=2) vp2(:nx,:ny,1:nz)
write(13,rec=3) wp2(:nx,:ny,1:nz)
write(13,rec=4) upwp(:nx,:ny,1:nz)
write(13,rec=5) vpwp(:nx,:ny,1:nz)
write(13,rec=6) upvp(:nx,:ny,1:nz)
close(13)
#endif
! Do the pressure-velocity correlation calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
allocate( pp2(nx,ny,lbz:nz) )
allocate( ppup(nx,ny,lbz:nz) )
allocate( ppvp(nx,ny,lbz:nz) )
allocate( ppwp(nx,ny,lbz:nz) )
pp2 = this % p2 - this % p * this % p
ppup = this % pu - this % p * this % u
ppvp = this % pv - this % p * this % v
ppwp = this % pw - this % p * this % w_uv  !the average value w_uv_f is the same

#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_ps,nx,ny,nz- nz_end,nz_tot,                     &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                                  &
    (/ 'Meanupup', 'Meanvpvp', 'Meanwpwp','Meanupwp','Meanvpwp','Meanupvp'/),  &
    (/ pp2(1:nx,1:ny,1:nz- nz_end) ,                                         &
    ppup(1:nx,1:ny,1:nz- nz_end) ,                                              &
    ppvp(1:nx,1:ny,1:nz- nz_end) ,                                              &
    ppwp(1:nx,1:ny,1:nz- nz_end)
     /) )
#else
! Write binary data
open(unit=13, file=fname_ps, form='unformatted', convert=write_endian,         &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) pp2(:nx,:ny,1:nz)
write(13,rec=2) ppup(:nx,:ny,1:nz)
write(13,rec=3) ppvp(:nx,:ny,1:nz)
write(13,rec=4) ppwp(:nx,:ny,1:nz)
close(13)
#endif
#ifdef PPSCALARS
! Do the theta-velocity correlation calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
allocate( thetap2(nx,ny,lbz:nz) )
allocate( thetapup(nx,ny,lbz:nz) )
allocate( thetapvp(nx,ny,lbz:nz) )
allocate( thetapwp(nx,ny,lbz:nz) )
allocate( thetappp(nx,ny,lbz:nz) )

thetap2 = this % theta2 - this % theta * this % theta
thetapup = this % thetau - this % u * this % theta
thetapvp = this % thetav - this % v * this % theta
thetapwp = this % thetaw - this % w * this % theta_w !! w grid
thetappp = this % thetap - this % p * this % theta

#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_rs_sca,nx,ny,nz- nz_end,nz_tot,                     &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                                  &
    (/ 'Meanupup', 'Meanvpvp', 'Meanwpwp','Meanupwp','Meanvpwp','Meanupvp'/),  &
    (/ thetap2(1:nx,1:ny,1:nz- nz_end) ,                                         &
    thetapup(1:nx,1:ny,1:nz- nz_end) ,                                              &
    thetapvp(1:nx,1:ny,1:nz- nz_end) ,                                              &
    thetapwp(1:nx,1:ny,1:nz- nz_end) ,                                             &
    thetappp(1:nx,1:ny,1:nz- nz_end) 
     /) )
#else
! Write binary data
open(unit=13, file=fname_rs_sca, form='unformatted', convert=write_endian,         &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) thetap2(:nx,:ny,1:nz)
write(13,rec=2) thetapup(:nx,:ny,1:nz)
write(13,rec=3) thetapvp(:nx,:ny,1:nz)
write(13,rec=4) thetapwp(:nx,:ny,1:nz)
write(13,rec=5) thetappp(:nx,:ny,1:nz)
close(13)
#endif
#endif

#ifdef PPBUDGET
! Do the p strain rate tensor correlation calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
allocate( ppS11p(nx,ny,lbz:nz) )
allocate( ppS12p(nx,ny,lbz:nz) )
allocate( ppS13p(nx,ny,lbz:nz) )
allocate( ppS22p(nx,ny,lbz:nz) )
allocate( ppS23p(nx,ny,lbz:nz) )
allocate( ppS33p(nx,ny,lbz:nz) )

ppS11p = this % pS11 - this % p * 0.5_rprec*(this % dudx + this % dudx)
ppS12p = this % pS12 - this % p * 0.5_rprec*(this % dudy + this % dvdx)
ppS13p = this % pS13 - this % p * 0.5_rprec*(this % dudz + this % dwdx)
ppS22p = this % pS22 - this % p * 0.5_rprec*(this % dvdy + this % dvdy)
ppS23p = this % pS23 - this % p * 0.5_rprec*(this % dvdz + this % dwdy)  
ppS33p = this % pS33 - this % p * 0.5_rprec*(this % dwdz + this % dwdz)

#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_ppSijp,nx,ny,nz- nz_end,nz_tot,                     &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                                  &
    (/ 'Meanupup', 'Meanvpvp', 'Meanwpwp','Meanupwp','Meanvpwp','Meanupvp'/),  &
    (/ ppS11p(1:nx,1:ny,1:nz- nz_end) ,                                         &
    ppS12p(1:nx,1:ny,1:nz- nz_end) ,                                              &
    ppS13p(1:nx,1:ny,1:nz- nz_end) ,                                              &
    ppS22p(1:nx,1:ny,1:nz- nz_end) ,                                             &
    ppS23p(1:nx,1:ny,1:nz- nz_end) ,                                             &
    ppS33p(1:nx,1:ny,1:nz- nz_end)  /) )
#else
! Write binary data
open(unit=13, file=fname_ppSijp, form='unformatted', convert=write_endian,         &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) ppS11p(:nx,:ny,1:nz)
write(13,rec=2) ppS12p(:nx,:ny,1:nz)
write(13,rec=3) ppS13p(:nx,:ny,1:nz)
write(13,rec=4) ppS22p(:nx,:ny,1:nz)
write(13,rec=5) ppS23p(:nx,:ny,1:nz)
write(13,rec=6) ppS33p(:nx,:ny,1:nz)
close(13)
#endif
! Do the triple velocity correlation calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
allocate( upupup(nx,ny,lbz:nz) )
allocate( vpvpvp(nx,ny,lbz:nz) )
allocate( wpwpwp(nx,ny,lbz:nz) )
allocate( upupvp(nx,ny,lbz:nz) )
allocate( upupwp(nx,ny,lbz:nz) )
allocate( vpvpup(nx,ny,lbz:nz) )
allocate( vpvpwp(nx,ny,lbz:nz) )
allocate( wpwpup(nx,ny,lbz:nz) )
allocate( wpwpvp(nx,ny,lbz:nz) )
allocate( upvpwp(nx,ny,lbz:nz) )

upupup = this % u3- 3 * this % u2 * this % u + 2 * (this % u)**3
vpvpvp = this % v3- 3 * this % v2 * this % v + 2 * (this % v)**3
wpwpwp = this % w3- 3 * this % w2_uv * this % w_uv + 2 * (this % w_uv)**3

upupvp = this % u2v - this % u2 * this % v - 2 * this % uv * this % u + &
2 * (this % u)**2 * (this % v)
upupwp = this % u2w - this % u2 * this % w_uv - 2 * this % uw_uv * this % u +&
 2 * (this % u)**2 * (this % w_uv)
vpvpup = this % v2u - this % v2 * this % u - 2 * this % uv * this % v + &
2 * (this % v)**2 * (this % u)
vpvpwp = this % v2w - this % v2 * this % w_uv - 2 * this % vw_uv * this % v + &
2 * (this % v)**2 * (this % w_uv)
wpwpup = this % w2u - this % w2_uv * this % u - 2 * this % uw_uv * this % w_uv +&
 2 * (this % w_uv)**2 * (this % u)
wpwpvp = this % w2v - this % w2_uv * this % v - 2 * this % vw_uv * this % w_uv +&
 2 * (this % w_uv)**2 * (this % v)

upvpwp = this % uvw - this % u * this % vw_uv - this % v * this % uw_uv -&
 this % w_uv * this % uv + 2 * this % u * this % v * this % w_uv

#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_rs,nx,ny,nz- nz_end,nz_tot,                     &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                                  &
    (/ 'Meanupup', 'Meanvpvp', 'Meanwpwp','Meanupwp','Meanvpwp','Meanupvp'/),  &
    (/ up2(1:nx,1:ny,1:nz- nz_end) ,                                         &
    vp2(1:nx,1:ny,1:nz- nz_end) ,                                              &
    wp2(1:nx,1:ny,1:nz- nz_end) ,                                              &
    upwp(1:nx,1:ny,1:nz- nz_end) ,                                             &
    vpwp(1:nx,1:ny,1:nz- nz_end) ,                                             &
    upvp(1:nx,1:ny,1:nz- nz_end)  /) )
#else
! Write binary data
open(unit=13, file=fname_uipujpukp, form='unformatted', convert=write_endian,         &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) upupup(:nx,:ny,1:nz)
write(13,rec=2) vpvpvp(:nx,:ny,1:nz)
write(13,rec=3) wpwpwp(:nx,:ny,1:nz)
write(13,rec=4) upupvp(:nx,:ny,1:nz)
write(13,rec=5) upupwp(:nx,:ny,1:nz)
write(13,rec=6) vpvpup(:nx,:ny,1:nz)
write(13,rec=7) vpvpwp(:nx,:ny,1:nz)
write(13,rec=8) wpwpup(:nx,:ny,1:nz)
write(13,rec=9) wpwpvp(:nx,:ny,1:nz)
write(13,rec=10) upvpwp(:nx,:ny,1:nz)
close(13)
#endif
! Do the triple velocity correlation calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
allocate( uptau11p(nx,ny,lbz:nz) )
allocate( uptau12p(nx,ny,lbz:nz) )
allocate( uptau13p(nx,ny,lbz:nz) )
allocate( uptau22p(nx,ny,lbz:nz) )
allocate( uptau23p(nx,ny,lbz:nz) )
allocate( uptau33p(nx,ny,lbz:nz) )

allocate( vptau11p(nx,ny,lbz:nz) )
allocate( vptau12p(nx,ny,lbz:nz) )
allocate( vptau13p(nx,ny,lbz:nz) )
allocate( vptau22p(nx,ny,lbz:nz) )
allocate( vptau23p(nx,ny,lbz:nz) )
allocate( vptau33p(nx,ny,lbz:nz) )

allocate( wptau11p(nx,ny,lbz:nz) )
allocate( wptau12p(nx,ny,lbz:nz) )
allocate( wptau13p(nx,ny,lbz:nz) )
allocate( wptau22p(nx,ny,lbz:nz) )
allocate( wptau23p(nx,ny,lbz:nz) )
allocate( wptau33p(nx,ny,lbz:nz) )

uptau11p = this % utau11 - this % u * this % txx
uptau12p = this % utau12 - this % u * this % txy
uptau13p = this % utau13 - this % u * this % txz_uv
uptau22p = this % utau22 - this % u * this % tyy
uptau23p = this % utau23 - this % u * this % tyz_uv 
uptau33p = this % utau33 - this % u * this % tzz

vptau11p = this % vtau11 - this % v * this % txx
vptau12p = this % vtau12 - this % v * this % txy
vptau13p = this % vtau13 - this % v * this % txz_uv
vptau22p = this % vtau22 - this % v * this % tyy
vptau23p = this % vtau23 - this % v * this % tyz_uv 
vptau33p = this % vtau33 - this % v * this % tzz

wptau11p = this % wtau11 - this % w_uv * this % txx
wptau12p = this % wtau12 - this % w_uv * this % txy
wptau13p = this % wtau13 - this % w_uv * this % txz_uv
wptau22p = this % wtau22 - this % w_uv * this % tyy
wptau23p = this % wtau23 - this % w_uv * this % tyz_uv
wptau33p = this % wtau33 - this % w_uv * this % tzz
#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_rs,nx,ny,nz- nz_end,nz_tot,                     &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                                  &
    (/ 'Meanupup', 'Meanvpvp', 'Meanwpwp','Meanupwp','Meanvpwp','Meanupvp'/),  &
    (/ up2(1:nx,1:ny,1:nz- nz_end) ,                                         &
    vp2(1:nx,1:ny,1:nz- nz_end) ,                                              &
    wp2(1:nx,1:ny,1:nz- nz_end) ,                                              &
    upwp(1:nx,1:ny,1:nz- nz_end) ,                                             &
    vpwp(1:nx,1:ny,1:nz- nz_end) ,                                             &
    upvp(1:nx,1:ny,1:nz- nz_end)  /) )
#else
! Write binary data
    open(unit=13, file=fname_uptaup, form='unformatted', convert=write_endian,         &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) uptau11p(:nx,:ny,1:nz)
write(13,rec=2) uptau12p(:nx,:ny,1:nz)
write(13,rec=3) uptau13p(:nx,:ny,1:nz)
write(13,rec=4) uptau22p(:nx,:ny,1:nz)
write(13,rec=5) uptau23p(:nx,:ny,1:nz)
write(13,rec=6) uptau33p(:nx,:ny,1:nz)

write(13,rec=7) vptau11p(:nx,:ny,1:nz)
write(13,rec=8) vptau12p(:nx,:ny,1:nz)
write(13,rec=9) vptau13p(:nx,:ny,1:nz)
write(13,rec=10) vptau22p(:nx,:ny,1:nz)
write(13,rec=11) vptau23p(:nx,:ny,1:nz)
write(13,rec=12) vptau33p(:nx,:ny,1:nz)

write(13,rec=13) wptau11p(:nx,:ny,1:nz)
write(13,rec=14) wptau12p(:nx,:ny,1:nz)
write(13,rec=15) wptau13p(:nx,:ny,1:nz)
write(13,rec=16) wptau22p(:nx,:ny,1:nz)
write(13,rec=17) wptau23p(:nx,:ny,1:nz)
write(13,rec=18) wptau33p(:nx,:ny,1:nz)
close(13)
#endif
! Do the triple velocity correlation calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
allocate( tau11pdudxp(nx,ny,lbz:nz) )
allocate( tau11pdvdxp(nx,ny,lbz:nz) )
allocate( tau11pdwdxp(nx,ny,lbz:nz) )
allocate( tau12pdudyp(nx,ny,lbz:nz) )
allocate( tau12pdvdyp(nx,ny,lbz:nz) )
allocate( tau12pdwdyp(nx,ny,lbz:nz) )
allocate( tau13pdudzp(nx,ny,lbz:nz) )
allocate( tau13pdvdzp(nx,ny,lbz:nz) )
allocate( tau13pdwdzp(nx,ny,lbz:nz) )

allocate( tau21pdudxp(nx,ny,lbz:nz) )
allocate( tau21pdvdxp(nx,ny,lbz:nz) )
allocate( tau21pdwdxp(nx,ny,lbz:nz) )
allocate( tau22pdudyp(nx,ny,lbz:nz) )
allocate( tau22pdvdyp(nx,ny,lbz:nz) )
allocate( tau22pdwdyp(nx,ny,lbz:nz) )
allocate( tau23pdudzp(nx,ny,lbz:nz) )
allocate( tau23pdvdzp(nx,ny,lbz:nz) )
allocate( tau23pdwdzp(nx,ny,lbz:nz) )

allocate( tau31pdudxp(nx,ny,lbz:nz) )
allocate( tau31pdvdxp(nx,ny,lbz:nz) )
allocate( tau31pdwdxp(nx,ny,lbz:nz) )
allocate( tau32pdudyp(nx,ny,lbz:nz) )
allocate( tau32pdvdyp(nx,ny,lbz:nz) )
allocate( tau32pdwdyp(nx,ny,lbz:nz) )
allocate( tau33pdudzp(nx,ny,lbz:nz) )
allocate( tau33pdvdzp(nx,ny,lbz:nz) )
allocate( tau33pdwdzp(nx,ny,lbz:nz) )

tau11pdudxp = this %tau11dudx - this %dudx * this %txx
tau11pdvdxp = this %tau11dvdx - this %dvdx * this %txx
tau11pdwdxp = this %tau11dwdx - this %dwdx * this %txx

tau12pdudyp = this %tau12dudy - this %dudy * this %txy
tau12pdvdyp = this %tau12dvdy - this %dvdy * this %txy
tau12pdwdyp = this %tau12dwdy - this %dwdy * this %txy

tau13pdudzp = this %tau13dudz - this %dudz * this %txz_uv
tau13pdvdzp = this %tau13dvdz - this %dvdz * this %txz_uv
tau13pdwdzp = this %tau13dwdz - this %dwdz * this %txz_uv

tau21pdudxp = this %tau21dudx - this %dudx * this %txy
tau21pdvdxp = this %tau21dvdx - this %dvdx * this %txy
tau21pdwdxp = this %tau21dwdx - this %dwdx * this %txy

tau22pdudyp = this %tau22dudy - this %dudy * this %tyy
tau22pdvdyp = this %tau22dvdy - this %dvdy * this %tyy
tau22pdwdyp = this %tau22dwdy - this %dwdy * this %tyy

tau23pdudzp = this %tau23dudz - this %dudz * this %tyz_uv
tau23pdvdzp = this %tau23dvdz - this %dvdz * this %tyz_uv
tau23pdwdzp = this %tau23dwdz - this %dwdz * this %tyz_uv

tau31pdudxp = this %tau31dudx - this %dudx * this %txz_uv
tau31pdvdxp = this %tau31dvdx - this %dvdx * this %txz_uv
tau31pdwdxp = this %tau31dwdx - this %dwdx * this %txz_uv

tau32pdudyp = this %tau32dudy - this %dudy * this %tyz_uv
tau32pdvdyp = this %tau32dvdy - this %dvdy * this %tyz_uv
tau32pdwdyp = this %tau32dwdy - this %dwdy * this %tyz_uv

tau33pdudzp = this %tau33dudz - this %dudz * this %tzz
tau33pdvdzp = this %tau33dvdz - this %dvdz * this %tzz
tau33pdwdzp = this %tau33dwdz - this %dwdz * this %tzz
#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_rs,nx,ny,nz- nz_end,nz_tot,                     &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                                  &
    (/ 'Meanupup', 'Meanvpvp', 'Meanwpwp','Meanupwp','Meanvpwp','Meanupvp'/),  &
    (/ up2(1:nx,1:ny,1:nz- nz_end) ,                                         &
    vp2(1:nx,1:ny,1:nz- nz_end) ,                                              &
    wp2(1:nx,1:ny,1:nz- nz_end) ,                                              &
    upwp(1:nx,1:ny,1:nz- nz_end) ,                                             &
    vpwp(1:nx,1:ny,1:nz- nz_end) ,                                             &
    upvp(1:nx,1:ny,1:nz- nz_end)  /) )
#else
! Write binary data
open(unit=13, file=fname_taupduidxjp, form='unformatted', convert=write_endian,         &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) tau11pdudxp(:nx,:ny,1:nz)
write(13,rec=2) tau11pdvdxp(:nx,:ny,1:nz)
write(13,rec=3) tau11pdwdxp(:nx,:ny,1:nz)

write(13,rec=4) tau12pdudyp(:nx,:ny,1:nz)
write(13,rec=5) tau12pdvdyp(:nx,:ny,1:nz)
write(13,rec=6) tau12pdwdyp(:nx,:ny,1:nz)

write(13,rec=7) tau13pdudzp(:nx,:ny,1:nz)
write(13,rec=8) tau13pdvdzp(:nx,:ny,1:nz)
write(13,rec=9) tau13pdwdzp(:nx,:ny,1:nz)

write(13,rec=10) tau21pdudxp(:nx,:ny,1:nz)
write(13,rec=11) tau21pdvdxp(:nx,:ny,1:nz)
write(13,rec=12) tau21pdwdxp(:nx,:ny,1:nz)

write(13,rec=13) tau22pdudyp(:nx,:ny,1:nz)
write(13,rec=14) tau22pdvdyp(:nx,:ny,1:nz)
write(13,rec=15) tau22pdwdyp(:nx,:ny,1:nz)

write(13,rec=16) tau23pdudzp(:nx,:ny,1:nz)
write(13,rec=17) tau23pdvdzp(:nx,:ny,1:nz)
write(13,rec=18) tau23pdwdzp(:nx,:ny,1:nz)

write(13,rec=19) tau31pdudxp(:nx,:ny,1:nz)
write(13,rec=20) tau31pdvdxp(:nx,:ny,1:nz)
write(13,rec=21) tau31pdwdxp(:nx,:ny,1:nz)

write(13,rec=22) tau32pdudyp(:nx,:ny,1:nz)
write(13,rec=23) tau32pdvdyp(:nx,:ny,1:nz)
write(13,rec=24) tau32pdwdyp(:nx,:ny,1:nz)

write(13,rec=25) tau33pdudzp(:nx,:ny,1:nz)
write(13,rec=26) tau33pdvdzp(:nx,:ny,1:nz)
write(13,rec=27) tau33pdwdzp(:nx,:ny,1:nz)
close(13)
#endif
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
! Do the p strain rate tensor correlation calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
allocate( fxpup(nx,ny,lbz:nz) )
allocate( fypup(nx,ny,lbz:nz) )
allocate( fzpup(nx,ny,lbz:nz) )
allocate( fxpvp(nx,ny,lbz:nz) )
allocate( fypvp(nx,ny,lbz:nz) )
allocate( fzpvp(nx,ny,lbz:nz) )
allocate( fxpwp(nx,ny,lbz:nz) )
allocate( fypwp(nx,ny,lbz:nz) )
allocate( fzpwp(nx,ny,lbz:nz) )

fxpup = this %fxu - this %fx * this %u
fypup = this %fyu - this %fy * this %u
fzpup = this %fzu - this %fz * this %u

fxpvp = this %fxv - this %fx * this %v
fypvp = this %fyv - this %fy * this %v
fzpvp = this %fzv - this %fz * this %v

fxpwp = this %fxw - this %fx * this %w
fypwp = this %fyw - this %fy * this %w
fzpwp = this %fzw - this %fz * this %w

#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_ppSijp,nx,ny,nz- nz_end,nz_tot,                     &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                                  &
    (/ 'Meanupup', 'Meanvpvp', 'Meanwpwp','Meanupwp','Meanvpwp','Meanupvp'/),  &
    (/ ppS11p(1:nx,1:ny,1:nz- nz_end) ,                                         &
    ppS12p(1:nx,1:ny,1:nz- nz_end) ,                                              &
    ppS13p(1:nx,1:ny,1:nz- nz_end) ,                                              &
    ppS22p(1:nx,1:ny,1:nz- nz_end) ,                                             &
    ppS23p(1:nx,1:ny,1:nz- nz_end) ,                                             &
    ppS33p(1:nx,1:ny,1:nz- nz_end)  /) )
#else
! Write binary data
    open(unit=13, file=fname_fipujp, form='unformatted', convert=write_endian,         &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) fxpup(:nx,:ny,1:nz)
write(13,rec=2) fypup(:nx,:ny,1:nz)
write(13,rec=3) fzpup(:nx,:ny,1:nz)
write(13,rec=4) fxpvp(:nx,:ny,1:nz)
write(13,rec=5) fypvp(:nx,:ny,1:nz)
write(13,rec=6) fzpvp(:nx,:ny,1:nz)
write(13,rec=7) fxpwp(:nx,:ny,1:nz)
write(13,rec=8) fypwp(:nx,:ny,1:nz)
write(13,rec=9) fzpwp(:nx,:ny,1:nz)
close(13)
#endif
#endif
#endif

#ifdef PPMPI
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
#endif

end subroutine finalize

!*******************************************************************************
subroutine checkpoint(this)
!*******************************************************************************
!
! This subroutine writes the restart data and is to be called by 'checkpoint'
! for intermediate checkpoints and by 'tavg_finalize' at the end of the
! simulation.
!
use param, only : write_endian, coord
use string_util
implicit none

class(tavg_t), intent(inout) :: this

character(64) :: fname

fname = checkpoint_tavg_file
#ifdef PPMPI
call string_concat( fname, '.c', coord)
#endif

!  Write data to tavg.out
open(1, file=fname, action='write', position='rewind',form='unformatted',      &
    convert=write_endian)
write(1) this%total_time
write(1) this%u
write(1) this%v
write(1) this%w
write(1) this%u_w
write(1) this%v_w
write(1) this%w_uv
write(1) this%u2
write(1) this%v2
write(1) this%w2
write(1) this%uv
write(1) this%uw
write(1) this%vw
write(1) this%txx
write(1) this%tyy
write(1) this%tzz
write(1) this%txy
write(1) this%txz
write(1) this%tyz
write(1) this%fx
write(1) this%fy
write(1) this%fz
write(1) this%cs_opt2
write(1) this%vortx
write(1) this%vorty
write(1) this%vortz
write(1) this%p2
write(1) this%pu
write(1) this%pv
write(1) this%pw
#ifdef PPSCALARS
    write(1) this%theta
    write(1) this%theta_w
    write(1) this%theta2
    write(1) this%thetau
    write(1) this%thetav
    write(1) this%thetaw
    write(1) this%thetap
    write(1) this%sgs_t3
#endif
#ifdef PPBUDGET
    write(1) this%uw_uv
    write(1) this%vw_uv
    write(1) this%w2_uv
    write(1) this%txz_uv
    write(1) this%tyz_uv

    write(1) this%dudx
    write(1) this%dudy
    write(1) this%dudz
    write(1) this%dvdx
    write(1) this%dvdy
    write(1) this%dvdz
    write(1) this%dwdx
    write(1) this%dwdy
    write(1) this%dwdz

    write(1) this%u3
    write(1) this%v3
    write(1) this%w3
    write(1) this%u2v
    write(1) this%u2w
    write(1) this%v2u
    write(1) this%v2w
    write(1) this%w2u
    write(1) this%w2v
    write(1) this%uvw

    write(1) this%pS11
    write(1) this%pS12
    write(1) this%pS13
    write(1) this%pS22
    write(1) this%pS23
    write(1) this%pS33

    write(1) this%utau11
    write(1) this%utau12
    write(1) this%utau13
    write(1) this%utau22
    write(1) this%utau23
    write(1) this%utau33
    write(1) this%vtau11
    write(1) this%vtau12
    write(1) this%vtau13
    write(1) this%vtau22
    write(1) this%vtau23
    write(1) this%vtau33
    write(1) this%wtau11
    write(1) this%wtau12
    write(1) this%wtau13
    write(1) this%wtau22
    write(1) this%wtau23
    write(1) this%wtau33

    write(1) this%tau11dudx
    write(1) this%tau11dvdx
    write(1) this%tau11dwdx
    write(1) this%tau12dudy
    write(1) this%tau12dvdy
    write(1) this%tau12dwdy
    write(1) this%tau13dudz
    write(1) this%tau13dvdz
    write(1) this%tau13dwdz

    write(1) this%tau21dudx
    write(1) this%tau21dvdx
    write(1) this%tau21dwdx
    write(1) this%tau22dudy
    write(1) this%tau22dvdy
    write(1) this%tau22dwdy
    write(1) this%tau23dudz
    write(1) this%tau23dvdz
    write(1) this%tau23dwdz

    write(1) this%tau31dudx
    write(1) this%tau31dvdx
    write(1) this%tau31dwdx
    write(1) this%tau32dudy
    write(1) this%tau32dvdy
    write(1) this%tau32dwdy
    write(1) this%tau33dudz
    write(1) this%tau33dvdz
    write(1) this%tau33dwdz

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
    write(1) this%fxu
    write(1) this%fyu
    write(1) this%fzu
    write(1) this%fxv
    write(1) this%fyv
    write(1) this%fzv
    write(1) this%fxw
    write(1) this%fyw
    write(1) this%fzw
#endif
#endif
close(1)

end subroutine checkpoint

#ifdef PPCGNS
#ifdef PPMPI
!*******************************************************************************
subroutine write_parallel_cgns (file_name, nx, ny, nz, nz_tot, start_n_in,     &
    end_n_in, xin, yin, zin, num_fields, fieldNames, input )
!*******************************************************************************
use param, only : coord
implicit none

integer, intent(in) :: nx, ny, nz, nz_tot, num_fields
! Name of file to be written
character(*), intent(in) :: file_name
! Name of fields we are writing
character(*), intent(in), dimension(:) :: fieldNames
! Data to be written
real(rprec), intent(in), dimension(:) :: input
! Coordinates to write
real(rprec), intent(in), dimension(:) :: xin, yin, zin
! Where the total node counter starts nodes
integer, intent(in) :: start_n_in(3)
! Where the total node counter ends nodes
integer, intent(in) :: end_n_in(3)

integer :: fn=1        ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base=1      ! base number
integer :: zone=1      ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1      ! solution number
integer :: field       ! section number
integer(cgsize_t) :: sizes(3,3)  ! Sizes

! Convert input to right data type
integer(cgsize_t) :: start_n(3)  ! Where the total node counter starts nodes
integer(cgsize_t) :: end_n(3)  ! Where the total node counter ends nodes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: xyz

!  ! Set the parallel communicator
!  call cgp_mpi_comm_f(cgnsParallelComm, ierr)

! Convert types such that CGNS libraries can handle the input
start_n(1) = int(start_n_in(1), cgsize_t)
start_n(2) = int(start_n_in(2), cgsize_t)
start_n(3) = int(start_n_in(3), cgsize_t)
end_n(1) = int(end_n_in(1), cgsize_t)
end_n(2) = int(end_n_in(2), cgsize_t)
end_n(3) = int(end_n_in(3), cgsize_t)

! The total number of nodes in this processor
nnodes = nx*ny*nz

! Sizes, used to create zone
sizes(:,1) = (/int(nx, cgsize_t),int(ny, cgsize_t),int(nz_tot, cgsize_t)/)
sizes(:,2) = (/int(nx-1, cgsize_t),int(ny-1, cgsize_t),int(nz_tot-1, cgsize_t)/)
sizes(:,3) = (/int(0, cgsize_t) , int(0, cgsize_t), int(0, cgsize_t)/)

! Open CGNS file
call cgp_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
end if

! Create data nodes for coordinates
call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

! This is done for the 3 dimensions x,y and z
! It writes the coordinates
! Create grid points
do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = xin(i)
end do
end do
end do

call cgp_coord_write_data_f(fn, base, zone, 1,                                 &
    start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = yin(j)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 2,   &
    start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = zin(k)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 3,   &
                            start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i=1,num_fields
    call cgp_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),     &
        field, ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

    call cgp_field_write_data_f(fn, base, zone, sol, field, start_n, end_n,    &
        input((i-1)*nnodes+1:(i)*nnodes), ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

end do

! Close the file
call cgp_close_f(fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

end subroutine write_parallel_cgns

!*******************************************************************************
subroutine write_null_cgns (file_name, nx, ny, nz, nz_tot, start_n_in,         &
    end_n_in, xin, yin, zin, num_fields, fieldNames )
!*******************************************************************************
use param, only : coord
implicit none

integer, intent(in) :: nx, ny, nz, nz_tot, num_fields
! Name of file to be written
character(*), intent(in) :: file_name
! Name of fields we are writing
character(*), intent(in), dimension(:) :: fieldNames
! Coordinates to write
real(rprec), intent(in), dimension(:) :: xin, yin, zin
! Where the total node counter starts nodes
integer, intent(in) :: start_n_in(3)
! Where the total node counter ends nodes
integer, intent(in) :: end_n_in(3)

integer :: fn=1        ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base=1      ! base number
integer :: zone=1      ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1      ! solution number
integer :: field       ! section number
integer(cgsize_t) :: sizes(3,3)  ! Sizes

! Convert input to right data type
integer(cgsize_t) :: start_n(3)  ! Where the total node counter starts nodes
integer(cgsize_t) :: end_n(3)  ! Where the total node counter ends nodes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: xyz

!  ! Set the parallel communicator
!  call cgp_mpi_comm_f(cgnsParallelComm, ierr)

! Convert types such that CGNS libraries can handle the input
start_n(1) = int(start_n_in(1), cgsize_t)
start_n(2) = int(start_n_in(2), cgsize_t)
start_n(3) = int(start_n_in(3), cgsize_t)
end_n(1) = int(end_n_in(1), cgsize_t)
end_n(2) = int(end_n_in(2), cgsize_t)
end_n(3) = int(end_n_in(3), cgsize_t)

! The total number of nodes in this processor
nnodes = nx*ny*nz

! Sizes, used to create zone
sizes(:,1) = (/int(nx, cgsize_t),int(ny, cgsize_t),int(nz_tot, cgsize_t)/)
sizes(:,2) = (/int(nx-1, cgsize_t),int(ny-1, cgsize_t),int(nz_tot-1, cgsize_t)/)
sizes(:,3) = (/int(0, cgsize_t) , int(0, cgsize_t), int(0, cgsize_t)/)

! Open CGNS file
call cgp_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
end if

! Create data nodes for coordinates
call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! This is done for the 3 dimensions x,y and z
! It writes the coordinates
! Create grid points
do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = xin(i)
end do
end do
end do

call cgp_coord_write_data_f(fn, base, zone, 1, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = yin(j)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 2, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = zin(k)
end do
end do
end do

call cgp_coord_write_data_f(fn, base, zone, 3, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i = 1, num_fields
    call cgp_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),     &
                           field, ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

    call cgp_field_write_data_f(fn, base, zone, sol, field, start_n, end_n,    &
                                %VAL(0), ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

end do

! Close the file
call cgp_close_f(fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

write(*,*) "end of write_null_cgns"

end subroutine write_null_cgns
#endif
#endif

end module time_average
