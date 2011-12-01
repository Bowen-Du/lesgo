!**********************************************************************
module input_util
!**********************************************************************
implicit none

save 
private

public read_input_conf

character (*), parameter :: input_conf = 'lesgo.conf'
character (*), parameter :: comment = '!'
!character (*), parameter :: ldelim = '('  !--no whitespace allowed
!character (*), parameter :: rdelim = ')'  !--no whitespace allowed
character (*), parameter :: block_entry = '{'
character (*), parameter :: block_exit = '}'
character (*), parameter :: equal = '='
character (*), parameter :: esyntax = 'syntax error at line'

! Delimiters used for reading vectors and points
character(*), parameter :: delim_minor=','
character(*), parameter :: delim_major='//'

! Default buffer length for characters of unknown length
integer, parameter :: BUFF_LEN = 256

interface parse_vector
  module procedure parse_vector_real, parse_vector_point3D
end interface

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_input_conf ()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
use messages
use string_util, only : eat_whitespace, uppercase
implicit none

character (*), parameter :: sub = 'read_input_conf'

integer, parameter :: lun = 1


character (BUFF_LEN) :: buff

integer :: block_entry_pos, block_exit_pos, equal_pos

integer :: ios
integer :: line

logical :: exst

! Check that the configuration file exists
inquire (file=input_conf, exist=exst)

if (exst) then
  open (lun, file=input_conf, action='read')
else
  call error (sub, 'file ' // input_conf // ' does not exist')
end if

line = 0
do

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )

  if (ios /= 0) exit

  if (block_entry_pos == 0) then  !--for now, invalid format if no block entry found
    call error (sub, 'block entry not found on line', line) 
  end if

  ! Find block
  select case (uppercase(buff(1:block_entry_pos-1)))

  case ('DOMAIN')

     call domain_block()

  case ('MODEL')
         
     call model_block()

  case ('TIME')
     
     call time_block()

  case ('FLOW_COND')

     call flow_cond_block()      
      
  case ('OUTPUT')

     call output_block()

  $if($LVLSET)
  case ('LEVEL_SET')
     
     call level_set_block()

    $if($RNS_LS)
    case ('RNS')
       call rns_block()
    $endif

    $if($CYL_SKEW_LS)
    case ('CYL_SKEW')
       call cyl_skew_block()
    $endif

  $endif

  case default

     call mesg( sub, 'Found unused input block: ' // buff(1:block_entry_pos-1) )
     ! Now need to 'fast-forward' untile we reach the end of the block
     do while ( block_exit_pos == 0 )
        call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
             equal_pos, ios )
        if (ios /= 0) exit ! exit if end of file is reached
     enddo

  end select
  
end do

close (lun)

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine domain_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use param
implicit none

character(*), parameter :: block_name = 'DOMAIN'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry() 
     
     select case (uppercase(buff(1:equal_pos-1)))

     case ('NX')
        read (buff(equal_pos+1:), *) Nx
     case ('NY')
        read (buff(equal_pos+1:), *) Ny
     case ('NZ') 
        read (buff(equal_pos+1:), *) Nz_tot
     case ('Z_I')
        read (buff(equal_pos+1:), *) z_i
     case ('LX')
        read (buff(equal_pos+1:), *) L_x
     case ('LY')
        read (buff(equal_pos+1:), *) L_y
     case ('LZ')
        read (buff(equal_pos+1:), *) L_z
     case default
        call mesg( sub, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )
     end select

  elseif ( block_exit_pos == 1 ) then

     ! === Set dependant variables ===
     
     ! Set the processor owned vertical grid spacing
     nz = ceiling ( real( nz_tot, rprec ) / nproc ) + 1
     nz_tot = ( nz - 1 ) * nproc + 1 
     ! Grid size for dealiasing
     nx2 = 3 * nx / 2
     ny2 = 3 * ny / 2
     ! Grid size for FFT's
     lh = nx / 2 + 1
     ld = 2 * lh
     lh_big = nx2 / 2 + 1
     ld_big = 2 * lh_big
     
     ! Grid spacing
     dx = L_x / nx
     dy = L_y / ny
     dz = L_z / ( nz_tot - 1 )

     return

  else

     call error( sub, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif
     
enddo

return
end subroutine domain_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine model_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

character(*), parameter :: block_name = 'MODEL'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block') 

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()  

     select case (uppercase(buff(1:equal_pos-1)))

     case ('SGS_MODEL')
        read (buff(equal_pos+1:), *) sgs_model
     case ('WALL_DAMP_EXP') 
        read (buff(equal_pos+1:), *) wall_damp_exp
     case ('CS_COUNT')
        read (buff(equal_pos+1:), *) cs_count
     case ('DYN_INIT')
        read (buff(equal_pos+1:), *) DYN_init
     case ('CO')
        read (buff(equal_pos+1:), *) Co
     case ('IFILTER')
        read (buff(equal_pos+1:), *) ifilter
     case ('U_STAR')
        read (buff(equal_pos+1:), *) u_star
     case ('VONK')
        read (buff(equal_pos+1:), *) vonk
     case ('CORIOLIS_FORCING')
        read (buff(equal_pos+1:), *) coriolis_forcing
     case ('CORIOL')
        read (buff(equal_pos+1:), *) coriol
     case ('UG')
        read (buff(equal_pos+1:), *) ug
     case ('VG')
        read (buff(equal_pos+1:), *) vg
     case ('NU_MOLEC')
        read (buff(equal_pos+1:), *) nu_molec
     case ('MOLEC')
        read (buff(equal_pos+1:), *) molec
     case ('SGS')
        read (buff(equal_pos+1:), *) sgs
     case ('DNS_BC')
        read (buff(equal_pos+1:), *) dns_bc
     case default

        call mesg( sub, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )

     end select

  elseif( block_exit_pos == 1 ) then

     return
     
  else

     call error( sub, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine model_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine time_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use param
implicit none

character(*), parameter :: block_name = 'TIME'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')


  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('NSTEPS')
        read (buff(equal_pos+1:), *) nsteps

     case ('USE_CFL_DT') 
        read (buff(equal_pos+1:), *) use_cfl_dt
        
     case ('CFL')
        read (buff(equal_pos+1:), *) cfl

     case('DT')
        read (buff(equal_pos+1:), *) dt

     case('CUMULATIVE_TIME')
        read (buff(equal_pos+1:), *) cumulative_time
     case default

        call mesg( sub, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )
     end select

  elseif( block_exit_pos == 1 ) then

     ! Set dependent data
     if( .not. use_cfl_dt ) then
       ! Set dimensional time step
       dt_dim = dt * z_i / u_star
       ! Set AB2 integration coefficients
       tadv1 = 1.5_rprec
       tadv2 = 1.0_rprec - tadv1
     endif

     return

  else

     call error( sub, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  time_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine flow_cond_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

character(*), parameter :: block_name = 'FLOW_COND'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('INITU')
        read (buff(equal_pos+1:), *) initu
     case ('INILAG')
        read (buff(equal_pos+1:), *) inilag
     case ('UBC')
        read (buff(equal_pos+1:), *) ubc
     case ('LBC_MOM')
        Read (buff(equal_pos+1:), *) lbc_mom
     case ('ZO')
        read (buff(equal_pos+1:), *) zo
     case ('INFLOW')
        read (buff(equal_pos+1:), *) inflow
     case ('FRINGE_REGION_END')
        read (buff(equal_pos+1:), *) fringe_region_end
     case ('FRINGE_REGION_LEN')
        read (buff(equal_pos+1:), *) fringe_region_len
     case ('UNIFORM_INFLOW')
        read (buff(equal_pos+1:), *) uniform_inflow
     case ('INFLOW_VELOCITY')
        read (buff(equal_pos+1:), *) inflow_velocity
     case ('FORCE_TOP_BOT')
        read (buff(equal_pos+1:), *) force_top_bot
     case ('USE_MEAN_P_FORCE')
        read (buff(equal_pos+1:), *) use_mean_p_force
     case ('MEAN_P_FORCE')
        read (buff(equal_pos+1:), *) mean_p_force

     case default

        call mesg( sub, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )

     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  flow_cond_block

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param
implicit none

character(*), parameter :: block_name = 'OUTPUT'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('WBASE')
        read (buff(equal_pos+1:), *) wbase

     case ('NENERGY')
        read (buff(equal_pos+1:), *) nenergy

     case ('LAG_CFL_COUNT')
        read (buff(equal_pos+1:), *) lag_cfl_count

     case ('TAVG_CALC')
        read (buff(equal_pos+1:), *) tavg_calc
     case ('TAVG_NSTART')
        read (buff(equal_pos+1:), *) tavg_nstart
     case ('TAVG_NEND')
        read (buff(equal_pos+1:), *) tavg_nend

     case ('POINT_CALC')
        read (buff(equal_pos+1:), *) point_calc
     case ('POINT_NSTART')
        read (buff(equal_pos+1:), *) point_nstart
     case ('POINT_NEND')
        read (buff(equal_pos+1:), *) point_nend
     case ('POINT_NSKIP')
        read (buff(equal_pos+1:), *) point_nskip
     case ('POINT_NLOC')
        read (buff(equal_pos+1:), *) point_nloc
     case ('POINT_LOC')
        allocate( point_loc( point_nloc ) )
        call parse_vector( buff(equal_pos+1:), point_loc )

     case ('DOMAIN_CALC')
        read (buff(equal_pos+1:), *) domain_calc
     case ('DOMAIN_NSTART')
        read (buff(equal_pos+1:), *) domain_nstart
     case ('DOMAIN_NEND')
        read (buff(equal_pos+1:), *) domain_nend
     case ('DOMAIN_NSKIP')
        read (buff(equal_pos+1:), *) domain_nskip

     case ('XPLANE_CALC')
        read (buff(equal_pos+1:), *) xplane_calc
     case ('XPLANE_NSTART')
        read (buff(equal_pos+1:), *) xplane_nstart
     case ('XPLANE_NEND')
        read (buff(equal_pos+1:), *) xplane_nend
     case ('XPLANE_NSKIP')
        read (buff(equal_pos+1:), *) xplane_nskip
     case ('XPLANE_NLOC')
        read (buff(equal_pos+1:), *) xplane_nloc
     case ('XPLANE_LOC')
        allocate( xplane_loc( xplane_nloc ) )
        call parse_vector( buff(equal_pos+1:), xplane_loc )

     case ('YPLANE_CALC')
        read (buff(equal_pos+1:), *) yplane_calc
     case ('YPLANE_NSTART')
        read (buff(equal_pos+1:), *) yplane_nstart
     case ('YPLANE_NEND')
        read (buff(equal_pos+1:), *) yplane_nend
     case ('YPLANE_NSKIP')
        read (buff(equal_pos+1:), *) yplane_nskip
     case ('YPLANE_NLOC')
        read (buff(equal_pos+1:), *) yplane_nloc
     case ('YPLANE_LOC')
        allocate( yplane_loc( yplane_nloc ) )
        call parse_vector( buff(equal_pos+1:), yplane_loc )

     case ('ZPLANE_CALC')
        read (buff(equal_pos+1:), *) zplane_calc
     case ('ZPLANE_NSTART')
        read (buff(equal_pos+1:), *) zplane_nstart
     case ('ZPLANE_NEND')
        read (buff(equal_pos+1:), *) zplane_nend
     case ('ZPLANE_NSKIP')
        read (buff(equal_pos+1:), *) zplane_nskip
     case ('ZPLANE_NLOC')
        read (buff(equal_pos+1:), *) zplane_nloc
     case ('ZPLANE_LOC')
        allocate( zplane_loc( zplane_nloc ) )
        call parse_vector( buff(equal_pos+1:), zplane_loc )

     case ('SPECTRA_CALC')
        read (buff(equal_pos+1:), *) spectra_calc
     case ('SPECTRA_NSTART')
        read (buff(equal_pos+1:), *) spectra_nstart
     case ('SPECTRA_NEND')
        read (buff(equal_pos+1:), *) spectra_nend
     case ('SPECTRA_NLOC')
        read (buff(equal_pos+1:), *) spectra_nloc
     case ('SPECTRA_LOC')
        allocate( spectra_loc( spectra_nloc ) )
        call parse_vector( buff(equal_pos+1:), spectra_loc )

     case default

        call mesg( sub, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )
     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  output_block

$if($LVLSET)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine level_set_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use level_set_base
implicit none

character(*), parameter :: block_name = 'LEVEL_SET'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('GLOBAL_CD_CALC') 
        read (buff(equal_pos+1:), *) global_CD_calc
     case ('LDIR')
        read (buff(equal_pos+1:), *) Ldir
     case ('VEL_BC')
        read (buff(equal_pos+1:), *) vel_bc
     case ('USE_LOG_PROFILE')
        Read (buff(equal_pos+1:), *) use_log_profile
     case ('USE_ENFORCE_UN')
        read (buff(equal_pos+1:), *) use_enforce_un
     case ('PHYSBC')
        read (buff(equal_pos+1:), *) physBC
     case ('USE_SMOOTH_TAU')
        read (buff(equal_pos+1:), *) use_smooth_tau
     case ('USE_EXTRAP_TAU_LOG')
        read (buff(equal_pos+1:), *) use_extrap_tau_log
     case ('USE_EXTRAP_TAU_SIMPLE')
        read (buff(equal_pos+1:), *) use_extrap_tau_simple
     case ('USE_MODIFY_DUTDN')
        read (buff(equal_pos+1:), *) use_modify_dutdn
     case ('LAG_DYN_MODIFY_BETA')
        read (buff(equal_pos+1:), *) lag_dyn_modify_beta
     case ('SMOOTH_MODE')
        read (buff(equal_pos+1:), *) smooth_mode
     case ('ZO_LEVEL_SET')
        read (buff(equal_pos+1:), *) zo_level_set
     
     $if($MPI)
     case ('NPHITOP')
        read (buff(equal_pos+1:), *) nphitop
     case ('NPHIBOT')
        read (buff(equal_pos+1:), *) nphibot
     case ('NVELTOP')
        read (buff(equal_pos+1:), *) nveltop
     case ('NVELBOT')
        read (buff(equal_pos+1:), *) nvelbot
     case ('NTAUTOP')
        read (buff(equal_pos+1:), *) ntautop
     case ('NTAUBOT')
        read (buff(equal_pos+1:), *) ntaubot
     case ('NFMMTOP')
        read (buff(equal_pos+1:), *) nFMMtop
     case ('NFMMBOT')
        read (buff(equal_pos+1:), *) nFMMbot
     $endif

     case default

        call mesg( sub, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )

     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  level_set_block

$if($RNS_LS)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine rns_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use rns_base_ls
implicit none

character(*), parameter :: block_name = 'RNS'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )
  if (ios /= 0) call error( sub, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('RNS_NTREE') 
        read (buff(equal_pos+1:), *) rns_ntree
     case ('RNS_TREE_LAYOUT')
        read (buff(equal_pos+1:), *) rns_tree_layout
     case ('TEMPORAL_WEIGHT')
        read (buff(equal_pos+1:), *) temporal_weight
     case ('TCONST')
        Read (buff(equal_pos+1:), *) Tconst
     case ('WEIGHT_NSTART')
        read (buff(equal_pos+1:), *) weight_nstart
     case ('TEMPORAL_MODEL')
        read (buff(equal_pos+1:), *) temporal_model
     case ('SPATIAL_MODEL')
        read (buff(equal_pos+1:), *) spatial_model
     case ('OUTPUT_NSKIP')
        read (buff(equal_pos+1:), *) output_nskip
     case ('CD_RAMP_NSTEP')
        read (buff(equal_pos+1:), *) CD_ramp_nstep
     case ('ALPHA_WIDTH')
        read (buff(equal_pos+1:), *) alpha_width
     case ('ALPHA_DIST')
        read (buff(equal_pos+1:), *) alpha_dist
     case ('CHI_CUTOFF')
        read (buff(equal_pos+1:), *) chi_cutoff

     case default

        call mesg( sub, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )

     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  rns_block
$endif

$if($CYL_SKEW_LS)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine cyl_skew_block()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use param, only : pi
use cyl_skew_base_ls, only : zrot_angle, skew_angle, use_bottom_surf, &
                             z_bottom_surf, ntree, tree_location, &
                             ngen, ngen_reslv, nbranch, d, l, offset, &
                             scale_fact, filter_chi, filt_width
implicit none

character(*), parameter :: block_name = 'CYL_SKEW'

do 

  call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                 equal_pos, ios )

  if (ios /= 0) call error( sub, 'Bad read in block')

  if( block_exit_pos == 0 ) then

     ! Check that the data entry conforms to correct format
     call checkentry()

     select case (uppercase(buff(1:equal_pos-1)))

     case ('ZROT_ANGLE') 
        read (buff(equal_pos+1:), *) zrot_angle
        ! Convert to radians
        zrot_angle = pi * zrot_angle / 180.0_rprec
     case ('SKEW_ANGLE')
        read (buff(equal_pos+1:), *) skew_angle
        ! Convert to radians
        skew_angle = pi * skew_angle / 180.0_rprec
     case ('USE_BOTTOM_SURF')
        read (buff(equal_pos+1:), *) use_bottom_surf
     case ('Z_BOTTOM_SURF')
        Read (buff(equal_pos+1:), *) z_bottom_surf
     case ('NTREE')
        read (buff(equal_pos+1:), *) ntree
     case ('TREE_LOCATION')
        allocate( tree_location( ntree ) )
        call parse_vector(buff(equal_pos+1:), tree_location )
     case ('NGEN')
        read (buff(equal_pos+1:), *) ngen
     case ('NGEN_RESLV')
        read (buff(equal_pos+1:), *) ngen_reslv
     case ('NBRANCH')
        read (buff(equal_pos+1:), *) nbranch
     case ('D')
        read (buff(equal_pos+1:), *) d
     case ('L')
        read (buff(equal_pos+1:), *) l
     case ('OFFSET')
        read (buff(equal_pos+1:), *) offset
     case ('SCALE_FACT')
        read (buff(equal_pos+1:), *) scale_fact
     case ('FILTER_CHI')
        read (buff(equal_pos+1:), *) filter_chi
     case ('FILT_WIDTH')
        read (buff(equal_pos+1:), *) filt_width

     case default

        call mesg( sub, 'Found unused data value in ' // block_name // ' block: ' // buff(1:equal_pos-1) )

     end select

  elseif( block_exit_pos == 1 ) then

     return

  else

     call error( sub, block_name // ' data block not formatted correctly: ' // buff(1:equal_pos-1) )

  endif

enddo

return
end subroutine  cyl_skew_block
$endif

$endif

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine checkentry()
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

if( equal_pos == 0 ) call error( sub, 'Bad read in block at line', line, ': ' // trim(adjustl(buff)))
!--invalid if nothing after equals
if (len_trim (buff) == equal_pos) call error (sub, 'nothing after equals sign in line', line) 

return
end subroutine checkentry  

end subroutine read_input_conf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine readline(lun, line, buff, block_entry_pos, &
                    block_exit_pos, equal_pos, ios )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! This subroutine reads the specified line and determines the attributes
! of the contents of the line.
!
use string_util, only : eat_whitespace
implicit none

integer, intent(in) :: lun
integer, intent(inout) :: line

character(*), intent(inout) :: buff
integer, intent(out) :: block_entry_pos, &
                        block_exit_pos, &
                        equal_pos, ios

block_entry_pos = 0
block_exit_pos = 0
equal_pos = 0
ios = -1

do     

  line = line + 1
  read (lun, '(a)', iostat=ios) buff
  if (ios /= 0) exit

  call eat_whitespace (buff)

  if (verify (buff, ' ') == 0) cycle  !--drop blank lines
  
  if (buff (1:len (comment)) == comment) cycle  !--drop comment lines

  block_entry_pos = index( buff, block_entry )
  block_exit_pos  = index( buff, block_exit )
  equal_pos       = index( buff, equal )

  exit

enddo 
return
end subroutine readline

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine parse_vector_real( string, vector )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec
use string_util, only : split_string

implicit none

character(*), intent(in) :: string
real(rprec), dimension(:), intent(inout) :: vector
character(BUFF_LEN), dimension(:), allocatable :: svector

integer :: nelem

! Get the number of elements in the vector
nelem = size(vector,1)
allocate( svector( nelem ) )

call split_string( string, delim_minor, nelem, svector )
read( svector, * ) vector

deallocate(svector)

return
end subroutine parse_vector_real

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine parse_vector_point3D( string, vector )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use types, only : rprec, point3D
use string_util, only : split_string

implicit none

character(*), intent(in) :: string
type(point3D), dimension(:), intent(inout) :: vector
character(BUFF_LEN), allocatable, dimension(:) :: svector

integer :: n, nelem
real(rprec), dimension(3) :: vector_minor

! Get the number of elements in the vector
nelem = size(vector,1)

allocate( svector( nelem ) )

! Split based on major delimiter
call split_string( string, delim_major, nelem, svector )
! Now parse result string 
do n=1, nelem
   call parse_vector_real( svector(n), vector_minor )
   vector(n) = point3D( (/ vector_minor(1), vector_minor(2), vector_minor(3) /) )
enddo

deallocate(svector)

return
end subroutine parse_vector_point3D

end module input_util
