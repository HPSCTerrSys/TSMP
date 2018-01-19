MODULE oas_icon_define

USE mo_kind,  ONLY : wp
USE netcdf
USE mod_oasis

IMPLICIT NONE

INTEGER               :: oas_comp_id
CHARACTER(len=2)      :: oas_comp_name="icon"
INTEGER               :: kl_comm              ! local communicator 
INTEGER               :: oas_error            ! return error code
INTEGER               :: oas_nlat = 2949120
INTEGER               :: oas_nlon = 1
INTEGER               :: oas_var_nodims(2)
INTEGER               :: oas_vshape(2)
INTEGER               :: oas_part_id
INTEGER, ALLOCATABLE  :: oas_part(:)
TYPE :: t_oas_field
  CHARACTER(len = 8)  :: clpname
  INTEGER             :: vid
END TYPE t_oas_field
TYPE(t_oas_field), DIMENSION(7)  :: oas_snd_fields
TYPE(t_oas_field), DIMENSION(5)  :: oas_rcv_fields
REAL(wp), POINTER     :: &
  oas_sw_snd(:), &
  oas_lw_snd(:), &
  oas_rain_snd(:), &
  oas_t_snd(:), &
  oas_p_snd(:), &
  oas_qv_snd(:), &
  oas_u_snd(:), &
  oas_sh_rcv(:), &
  oas_lh_rcv(:), &
  oas_tau_rcv(:), &
  oas_lw_rcv(:), &
  oas_alb_rcv(:)
  

PUBLIC :: &
 oas_comp_id, oas_comp_name, kl_comm, oas_error, oas_nlat, &
 oas_nlon, oas_var_nodims, oas_vshape, oas_part, &
 oas_snd_fields, oas_rcv_fields

END MODULE oas_icon_define
