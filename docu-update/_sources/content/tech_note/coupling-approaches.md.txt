# Coupling approaches

## Coupling between COSMO and CLM

### Update of transfer coefficients in COSMO based on CLM surface fluxes.

1\. Flux Inversion Scheme (INV)

```text
    -> CLM surface fluxes are inverted to obtain the transfer coefficients.
    -> Use --cplscheme=false when building the source code
```

2\. Transfer Coefficient Scheme (TRN)

```text
    -> CLM transfer coefficients are directly used to update the COSMO transfer coefficients.
    -> Use --cplscheme=true --readclm=true when building the source code
```

### Coupling Time Step

```text
The coupling frequency should be always less than the radiation update frequency in COSMO controlled by the namelist
(hincrad) in COSMO, which is a time step in hours.

Within the uncertainty of the simulated surface fluxes and the vertical transport of heat/moisture into the boundary layer,
a coupling time step of 60 s appears to be quite appropriate in terms of rmsd and computational efficiency. 
The INV scheme is less sensitive to larger coupling intervals than the TRN scheme while both perform equally well at 
coupling time steps <= 60 s. The findings are based from an unpublished study of idealized simulations (idealRTD) 
with the above two schemes using different coupling frequency. 
```

### Use of same grid between COSMO and CLM

COSMO and CLM uses rotated and regular geographic coordinates respectively. So, there is always a need of spatial interpolation for the coupling between the component models. However, there exists an additional flag in the TerrSysMP interface, which allows coupling without explicit interpolation between the two component models.

```text
   -> Use --readclm='True' when building the source code. In this setup, OASIS coupler will use the same grid/mask for 
        both component models. This grid/mask is based on CLM grid data.
   -> The CLM domain size for this setup should be 2*nboundlines smaller than the COSMO domain in both X and Y extents. 
        This is basically cropping the outer halo regions of the COSMO domain and 'nboundlines' is a parameter in the namelist 
        of COSMO, which specifies the number of halos for each subdomains. A minimum of 3 boundary lines are required for 
        the Runge-Kutta dynamical core in COSMO.
```

### Interpolation between COSMO and CLM grids

Bilinear Interpolation — (COSMO to CLM grids)  
Distance Weighted Interpolation — (CLM to COSMO grids)

```text
    Remapping files are generated at runtime, if they are not existing in the run directory. The distance weighted interpolation scheme uses a small number 'epsilon' to prevent division by 'zero', so it will never generate a remapping weight of '1' even if COSMO and CLM uses same grid. So, when COSMO and CLM are using the same grids, the weights in the remapping file has to be forced to '1'.  A pre-processing script is available in tps for fixing the weights in the remapping file. 
```

## Coupling between CLM and ParFlow

1\. The time step of the two component models and the coupling time step should always be the same.

2\. CLM and ParFlow are only coupled for the top 10 layers. The vertical discretization of this top 10 layers are specified ‘iniTimeConst.F90’. This has to always match with the discretization of top 10 layers specified in the TCL script of ParFlow. If the user wants to change the vertical discretization, it needs to be changed accordingly.

3\. The porosity in CLM and ParFlow are specified via the input surface data and the TCL script respectively. In CLM, pedotransfer functions are used to convert the sand percentage in surface data to porosity. This porosity has to match with the porosity specified in TCL script for the first 10 layers.
