MODULE oas_clm_define
#ifdef COUP_OAS_CLM

USE mo_kind,  ONLY : wp
USE netcdf

USE mod_oasis

IMPLICIT NONE

INTEGER               :: oas_comp_id
CHARACTER(len=4)      :: oas_comp_name="clm"
INTEGER               :: kl_comm              ! local communicator 
INTEGER               :: oas_error            ! return error code
INTEGER               :: oas_nlat
INTEGER               :: oas_nlon
INTEGER               :: oas_var_nodims(2)
INTEGER               :: oas_vshape(2)
INTEGER               :: oas_part_id
INTEGER, ALLOCATABLE  :: oas_part(:)
TYPE :: t_oas_field
  CHARACTER(len = 8)  :: clpname
  INTEGER             :: vid
END TYPE t_oas_field
TYPE(t_oas_field), DIMENSION(11)  :: oas_snd_meta
TYPE(t_oas_field), DIMENSION(8)  :: oas_rcv_meta
REAL, ALLOCATABLE     :: oas_snd_field(:,:), oas_rcv_field(:,:)
REAL(wp), ALLOCATABLE :: oas_rcv_field_icon(:,:,:)

PUBLIC :: &
 oas_comp_id, oas_comp_name, kl_comm, oas_error, oas_nlat, &
 oas_nlon, oas_var_nodims, oas_vshape, oas_part, &
 oas_snd_field, oas_rcv_field, oas_snd_meta, oas_rcv_meta
#endif

END MODULE oas_clm_define
