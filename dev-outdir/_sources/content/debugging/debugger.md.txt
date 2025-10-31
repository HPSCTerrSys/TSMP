# Debugging with TotalView

It is possible to debug TSMP with multiple coupled components with the mean of the TotalView debugger. 

Please contact the SDLTS for support for the initial set up TotalView on JSC machines. 

1) Start interactive session. E.g. on JUWELS machine 
```sh
salloc --partition=devel --nodes=5 --account=PROJECT --time=01:30:00
```

2) Setup and modifiy the submission script of TSMP
```sh
totalview -args srun --multi-prog slm_multiprog_mapping.conf
```
