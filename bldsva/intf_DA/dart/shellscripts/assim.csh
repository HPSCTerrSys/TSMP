./quickbuild.csh -mpi
rm input.nml.*
./create_obs_sequence
cat set_def.out
/create_fixed_network_seq
cat obs_seq.in
./cosmo_to_dart
./perfect_model_obs
# creates 3 files including obs_seq.perfect 

