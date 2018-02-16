MODULE oas_icon_define

!USE mo_kind,  ONLY : wp
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
TYPE(t_oas_field), DIMENSION(11)  :: oas_snd_meta
TYPE(t_oas_field), DIMENSION(7)  :: oas_rcv_meta
REAL, ALLOCATABLE     :: oas_snd_field(:,:), oas_rcv_field(:,:)

PUBLIC :: &
 oas_comp_id, oas_comp_name, kl_comm, oas_error, oas_nlat, &
 oas_nlon, oas_var_nodims, oas_vshape, oas_part, &
 oas_snd_field, oas_rcv_field, oas_snd_meta, oas_rcv_meta

END MODULE oas_icon_define
