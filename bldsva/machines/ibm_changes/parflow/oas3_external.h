extern int oas_pfl_vardef_mp_comp_id_;
extern int oas_pfl_vardef_mp_ierror_;
extern int localcomm;

#define oas_pfl_init oas_pfl_init
#define CALL_oas_pfl_init(arg1) oas_pfl_init(arg1);
void oas_pfl_init(int *arg1);
#define oas_pfl_finalize oas_pfl_finalize
#define CALL_oas_pfl_finalize(arg2) oas_pfl_finalize(arg2);
void oas_pfl_finalize(int *arg2);
