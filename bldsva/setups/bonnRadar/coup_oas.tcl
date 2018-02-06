# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*
#

#--------------------------------------------------------
pfset FileVersion 4
#-------------------------------------------------------
pfset Process.Topology.P __nprocx_pfl_bldsva__ 
pfset Process.Topology.Q __nprocy_pfl_bldsva__
pfset Process.Topology.R 1

# THE COMPUTATIONAL GRID IS A (BOX) THT CONTAINS THE MAIN PROBLEM. THIS CAN EITHER BE EXACTLY THE SIZE
# OF THE PROBLEM OR LARGER. A BOX GEOMETRY IN PARFLOW CAN BE ASIGNED BY EITHER SPECIFYING COORDINATES FOR
# TWO CORNERS OF THE BOX OR GRID SIZE AND NUMBER OF CELLS IN X,Y, AND Z.
#------------------------------------------------------------------------
# Computational Grid: It Defines The Grid Resolutions within The Domain
#------------------------------------------------------------------------
pfset ComputationalGrid.Lower.X			 0.0
pfset ComputationalGrid.Lower.Y			 0.0
pfset ComputationalGrid.Lower.Z			 -30.

pfset ComputationalGrid.DX			 1132.
pfset ComputationalGrid.DY		         1132. 
pfset ComputationalGrid.DZ			 1.00

pfset ComputationalGrid.NX			 __ngpflx_bldsva__ 
pfset ComputationalGrid.NY			 __ngpfly_bldsva__ 
pfset ComputationalGrid.NZ			 30 

# DOMAIN GEOMETRY IS THE (EXACTLY) OUTER DOMAIN OR BOUNDARY OF THE MODEL PROBLEM. IT HAS TO BE CONTAINED WITHIN THE COMPUTATIONAL GRID (i.e.
# OR HAS THE SAME SIZE OF IT). THE DOMAIN GEOMETRY COULD BE A BOX OR IT COULD BE A SOLID-FILE.
# BOUNDARY CONDITIONS ARE ASSIGNED TO THE DOMAIN SIDES WITH SOMETHING CALLED (PATCHES) IN THE TCL-SCRIPT.
# A BOX HAS SIX (6) SIDES AND (6) PATCHES WHILE A SOLID-FILE CAN HAVE ANY NUMBER OF PATCHES.
#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------
pfset Domain.GeomName                            domain

#---------------------------------------------------------
# Domain Geometry Input 
#---------------------------------------------------------
pfset GeomInput.Names	                     "domaininput indinput"
pfset GeomInput.domaininput.InputType         Box
pfset GeomInput.domaininput.GeomName          domain

pfset Geom.domain.Lower.X                     [pfget ComputationalGrid.Lower.X] 
pfset Geom.domain.Lower.Y                     [pfget ComputationalGrid.Lower.Y] 
pfset Geom.domain.Lower.Z                     [pfget ComputationalGrid.Lower.Z] 

set DX                                        [pfget ComputationalGrid.DX]
set DY                                        [pfget ComputationalGrid.DY]
set NX                                        [pfget ComputationalGrid.NX]
set NY                                        [pfget ComputationalGrid.NY]

pfset Geom.domain.Upper.X                     [expr $NX*$DX]
pfset Geom.domain.Upper.Y                     [expr $NY*$DY] 
pfset Geom.domain.Upper.Z                     0.0
pfset Geom.domain.Patches	             "x-lower x-upper y-lower y-upper z-lower z-upper"

pfset GeomInput.indinput.InputType  IndicatorField
pfset GeomInput.indinput.GeomNames  "sand sloam loam cloam clay ac1 ac2 ac3 ac4 ac5 ac6"
pfset Geom.indinput.FileName  "__forcingdir__/bonnRadarSoilInd.pfb"

pfset GeomInput.sand.Value               3
pfset GeomInput.sloam.Value              4
pfset GeomInput.loam.Value               5
pfset GeomInput.cloam.Value              6
pfset GeomInput.clay.Value               7
pfset GeomInput.ac1.Value               11
pfset GeomInput.ac2.Value               12
pfset GeomInput.ac3.Value               13
pfset GeomInput.ac4.Value               14
pfset GeomInput.ac5.Value               15
pfset GeomInput.ac6.Value               16
#-----------------------------------------------------------------------------
# VARIABLE dz ASSIGNMENTS
#-----------------------------------------------------------------------------
pfset Solver.Nonlinear.VariableDz            True
pfset dzScale.GeomNames                      domain
pfset dzScale.Type                           nzList
pfset dzScale.nzListNumber                   [pfget ComputationalGrid.NZ]
pfset Cell.0.dzScale.Value                   1.35
pfset Cell.1.dzScale.Value                   1.35
pfset Cell.2.dzScale.Value                   1.35
pfset Cell.3.dzScale.Value                   1.35
pfset Cell.4.dzScale.Value                   1.35
pfset Cell.5.dzScale.Value                   1.35
pfset Cell.6.dzScale.Value                   1.35
pfset Cell.7.dzScale.Value                   1.35
pfset Cell.8.dzScale.Value                   1.35
pfset Cell.9.dzScale.Value                   1.35
pfset Cell.10.dzScale.Value                   1.35
pfset Cell.11.dzScale.Value                   1.35
pfset Cell.12.dzScale.Value                   1.35
pfset Cell.13.dzScale.Value                   1.35
pfset Cell.14.dzScale.Value                   1.35
pfset Cell.15.dzScale.Value                   1.35
pfset Cell.16.dzScale.Value                   1.35
pfset Cell.17.dzScale.Value                   1.35
pfset Cell.18.dzScale.Value                   1.35
pfset Cell.19.dzScale.Value                   1.35
pfset Cell.20.dzScale.Value                   1.00
pfset Cell.21.dzScale.Value                   0.70
pfset Cell.22.dzScale.Value                   0.50
pfset Cell.23.dzScale.Value                   0.30
pfset Cell.24.dzScale.Value                   0.20
pfset Cell.25.dzScale.Value                   0.13
pfset Cell.26.dzScale.Value                   0.07
pfset Cell.27.dzScale.Value                   0.05
pfset Cell.28.dzScale.Value                   0.03
pfset Cell.29.dzScale.Value                   0.02

# TIME INFORMATION &
# TIME SETUP
#-----------------------------------------------------------------------------
# Setup timing info
#-----------------------------------------------------------------------------
pfset TimingInfo.BaseUnit		 __base_pfl__
pfset TimingInfo.StartCount		 __start_cnt_pfl__
pfset TimingInfo.StartTime		 0.0
pfset TimingInfo.StopTime		 __stop_pfl_bldsva__ 
pfset TimeStep.Type			 Constant
pfset TimeStep.Value			 __dt_pfl_bldsva__ 
pfset TimingInfo.DumpInterval		 __dump_pfl_interval__

# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names "constant"
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat             -1

#  HYDROLOGICAL PARAMETERS
#Schaap and Leiz (1998), Soil Science
#  SETUP AND VALUES
#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------
pfset Geom.Perm.Names			 "sand sloam loam cloam clay ac1 ac2 ac3 ac4 ac5 ac6"

pfset Geom.sand.Perm.Type		 Constant
pfset Geom.sand.Perm.Value		 0.269

pfset Geom.sloam.Perm.Type               Constant
pfset Geom.sloam.Perm.Value              0.0158

pfset Geom.loam.Perm.Type                Constant
pfset Geom.loam.Perm.Value               0.0050

pfset Geom.cloam.Perm.Type               Constant
pfset Geom.cloam.Perm.Value              0.0034

pfset Geom.clay.Perm.Type                Constant
pfset Geom.clay.Perm.Value               0.0062

pfset Geom.ac1.Perm.Type              Constant
pfset Geom.ac1.Perm.Value             0.02

pfset Geom.ac2.Perm.Type             Constant
pfset Geom.ac2.Perm.Value            0.01

pfset Geom.ac3.Perm.Type             Constant
pfset Geom.ac3.Perm.Value            0.02

pfset Geom.ac4.Perm.Type             Constant
pfset Geom.ac4.Perm.Value            0.01

pfset Geom.ac5.Perm.Type             Constant
pfset Geom.ac5.Perm.Value            0.0005

pfset Geom.ac6.Perm.Type             Constant
pfset Geom.ac6.Perm.Value            0.0001

pfset Perm.TensorType			 TensorByGeom
pfset Geom.Perm.TensorByGeom.Names	 "domain"
pfset Geom.domain.Perm.TensorValX	 10.0
pfset Geom.domain.Perm.TensorValY	 10.0
pfset Geom.domain.Perm.TensorValZ	 1.0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------
pfset SpecificStorage.Type			 Constant
pfset SpecificStorage.GeomNames			 "domain"
pfset Geom.domain.SpecificStorage.Value		 1.0e-4

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------
pfset Phase.Names			 "water"
pfset Phase.water.Density.Type		 Constant
pfset Phase.water.Density.Value		 1.0
pfset Phase.water.Viscosity.Type	 Constant
pfset Phase.water.Viscosity.Value	 1.0

#-----------------------------------------------------------------------------
# Gravity
#-----------------------------------------------------------------------------
pfset Gravity				 1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------
pfset Contaminants.Names		 ""

#-----------------------------------------------------------------------------
# Retardation
#-----------------------------------------------------------------------------
pfset Geom.Retardation.GeomNames	 ""

#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------
pfset Geom.Porosity.GeomNames           "sand sloam loam cloam clay ac1 ac2 ac3 ac4 ac5 ac6"

pfset Geom.sand.Porosity.Type          Constant
pfset Geom.sand.Porosity.Value         0.3756

pfset Geom.sloam.Porosity.Type          Constant
pfset Geom.sloam.Porosity.Value        0.4071

pfset Geom.loam.Porosity.Type          Constant
pfset Geom.loam.Porosity.Value         0.4386

pfset Geom.cloam.Porosity.Type          Constant
pfset Geom.cloam.Porosity.Value        0.4449

pfset Geom.clay.Porosity.Type          Constant
pfset Geom.clay.Porosity.Value         0.4701

pfset Geom.ac1.Porosity.Type          Constant
pfset Geom.ac1.Porosity.Value         0.33

pfset Geom.ac2.Porosity.Type          Constant
pfset Geom.ac2.Porosity.Value         0.33

pfset Geom.ac3.Porosity.Type          Constant
pfset Geom.ac3.Porosity.Value         0.33

pfset Geom.ac4.Porosity.Type          Constant
pfset Geom.ac4.Porosity.Value         0.33

pfset Geom.ac5.Porosity.Type          Constant
pfset Geom.ac5.Porosity.Value         0.33

pfset Geom.ac6.Porosity.Type          Constant
pfset Geom.ac6.Porosity.Value         0.33

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------
pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          "sand sloam loam cloam clay ac1 ac2 ac3 ac4 ac5 ac6"

pfset Geom.sand.RelPerm.Alpha         3.55 
pfset Geom.sand.RelPerm.N             4.16

pfset Geom.sloam.RelPerm.Alpha        2.69
pfset Geom.sloam.RelPerm.N            2.45

pfset Geom.loam.RelPerm.Alpha         1.12
pfset Geom.loam.RelPerm.N             2.48

pfset Geom.cloam.RelPerm.Alpha        1.58
pfset Geom.cloam.RelPerm.N            2.41

pfset Geom.clay.RelPerm.Alpha         1.51
pfset Geom.clay.RelPerm.N             2.26

pfset Geom.ac1.RelPerm.Alpha        1.0
pfset Geom.ac1.RelPerm.N            3.0

pfset Geom.ac2.RelPerm.Alpha        1.0
pfset Geom.ac2.RelPerm.N            3.0

pfset Geom.ac3.RelPerm.Alpha        1.0
pfset Geom.ac3.RelPerm.N            3.0

pfset Geom.ac4.RelPerm.Alpha        1.0
pfset Geom.ac4.RelPerm.N            3.0

pfset Geom.ac5.RelPerm.Alpha        1.0
pfset Geom.ac5.RelPerm.N            3.0

pfset Geom.ac6.RelPerm.Alpha        1.0
pfset Geom.ac6.RelPerm.N            3.0

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------
pfset Phase.Saturation.Type              VanGenuchten
pfset Phase.Saturation.GeomNames         "sand sloam loam cloam clay ac1 ac2 ac3 ac4 ac5 ac6"

pfset Geom.sand.Saturation.Alpha        3.55
pfset Geom.sand.Saturation.N            4.16
pfset Geom.sand.Saturation.SRes         0.14
pfset Geom.sand.Saturation.SSat         1.0

pfset Geom.sloam.Saturation.Alpha       2.69
pfset Geom.sloam.Saturation.N           2.45
pfset Geom.sloam.Saturation.SRes        0.1
pfset Geom.sloam.Saturation.SSat        1.0

pfset Geom.loam.Saturation.Alpha        1.12
pfset Geom.loam.Saturation.N            2.48
pfset Geom.loam.Saturation.SRes         0.15
pfset Geom.loam.Saturation.SSat         1.0

pfset Geom.cloam.Saturation.Alpha       1.58
pfset Geom.cloam.Saturation.N           2.41
pfset Geom.cloam.Saturation.SRes        0.18
pfset Geom.cloam.Saturation.SSat        1.0

pfset Geom.clay.Saturation.Alpha        1.51
pfset Geom.clay.Saturation.N            2.26
pfset Geom.clay.Saturation.SRes         0.21
pfset Geom.clay.Saturation.SSat         1.0

pfset Geom.ac1.Saturation.Alpha       1.0
pfset Geom.ac1.Saturation.N           3.0
pfset Geom.ac1.Saturation.SRes        0.001
pfset Geom.ac1.Saturation.SSat        1.0

pfset Geom.ac2.Saturation.Alpha       1.0
pfset Geom.ac2.Saturation.N           3.0
pfset Geom.ac2.Saturation.SRes        0.001
pfset Geom.ac2.Saturation.SSat        1.0

pfset Geom.ac3.Saturation.Alpha       1.0
pfset Geom.ac3.Saturation.N           3.0
pfset Geom.ac3.Saturation.SRes        0.001
pfset Geom.ac3.Saturation.SSat        1.0

pfset Geom.ac4.Saturation.Alpha       1.0
pfset Geom.ac4.Saturation.N           3.0
pfset Geom.ac4.Saturation.SRes        0.001
pfset Geom.ac4.Saturation.SSat        1.0

pfset Geom.ac5.Saturation.Alpha       1.0
pfset Geom.ac5.Saturation.N           3.0
pfset Geom.ac5.Saturation.SRes        0.001
pfset Geom.ac5.Saturation.SSat        1.0

pfset Geom.ac6.Saturation.Alpha       1.0
pfset Geom.ac6.Saturation.N           3.0
pfset Geom.ac6.Saturation.SRes        0.001
pfset Geom.ac6.Saturation.SSat        1.0

#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names				 ""

#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------

pfset Geom.domain.Patches             "x-lower x-upper y-lower y-upper z-lower z-upper"
pfset BCPressure.PatchNames [pfget Geom.domain.Patches]

pfset Patch.x-lower.BCPressure.Type                   FluxConst
pfset Patch.x-lower.BCPressure.Cycle                  "constant"
pfset Patch.x-lower.BCPressure.alltime.Value          0.0

pfset Patch.y-lower.BCPressure.Type                   FluxConst
pfset Patch.y-lower.BCPressure.Cycle                  "constant"
pfset Patch.y-lower.BCPressure.alltime.Value          0.0

pfset Patch.z-lower.BCPressure.Type                   FluxConst
pfset Patch.z-lower.BCPressure.Cycle                  "constant"
pfset Patch.z-lower.BCPressure.alltime.Value          0.0

pfset Patch.x-upper.BCPressure.Type                   FluxConst
pfset Patch.x-upper.BCPressure.Cycle                  "constant"
pfset Patch.x-upper.BCPressure.alltime.Value          0.0

pfset Patch.y-upper.BCPressure.Type                   FluxConst
pfset Patch.y-upper.BCPressure.Cycle                  "constant"
pfset Patch.y-upper.BCPressure.alltime.Value          0.0

pfset Patch.z-upper.BCPressure.Type                   OverlandFlow
pfset Patch.z-upper.BCPressure.Cycle                  "constant"
pfset Patch.z-upper.BCPressure.alltime.Value           0.0

#  TOPOGRAPHY & SLOPES IN
#  BOTH X- & Y- DIRECTIONS
#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------
pfset TopoSlopesX.Type			 "PFBFile"
pfset TopoSlopesX.GeomNames		 "domain"
pfset TopoSlopesX.FileName		 "__forcingdir__/xslope.pfb"

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------
pfset TopoSlopesY.Type			 "PFBFile"
pfset TopoSlopesY.GeomNames		 "domain"
pfset TopoSlopesY.FileName		 "__forcingdir__/yslope.pfb"

#---------------------------------------------------------
# Mannings coefficient
#---------------------------------------------------------
pfset Mannings.Type			 "Constant"
pfset Mannings.GeomNames		 "domain"
pfset Mannings.Geom.domain.Value	 5.52e-6

#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------
#
pfset ICPressure.Type                    __pfl_ICPpressureType__
pfset ICPressure.GeomNames               domain
pfset Geom.domain.ICPressure.Value       __pfl_ICPpressureValue__
pfset Geom.domain.ICPressure.FileName    "__pfl_ICPpressureFileName__"
pfset Geom.domain.ICPressure.RefGeom     domain
pfset Geom.domain.ICPressure.RefPatch    z-upper

#pfset ICPressure.Type                   "PFBFile"
#pfset ICPressure.GeomNames              "domain"
#pfset Geom.domain.ICPressure.FileName   "./bonnRadar_ic_press.pfb"

#
#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------
pfset PhaseSources.water.Type			 Constant
pfset PhaseSources.water.GeomNames		 domain
pfset PhaseSources.water.Geom.domain.Value	 0.0

#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------
pfset KnownSolution				 NoKnownSolution

##---------------------------------------------------------
## Spinup Flags
###---------------------------------------------------------
#pfset OverlandFlowSpinUp           1
#pfset OverlandSpinupDampP1  10.0
#pfset OverlandSpinupDampP2  0.1
#

# Set solver parameters
#-----------------------------------------------------------------------------
pfset Solver					 Richards
pfset Solver.MaxIter				 10000

pfset Solver.TerrainFollowingGrid                True

pfset Solver.Nonlinear.MaxIter			 500
pfset Solver.Nonlinear.ResidualTol		 1e-5
pfset Solver.Nonlinear.EtaChoice		 Walker1
pfset Solver.Nonlinear.EtaChoice		 EtaConstant
pfset Solver.Nonlinear.EtaValue			 0.001
pfset Solver.Nonlinear.UseJacobian		 False
pfset Solver.Nonlinear.DerivativeEpsilon	 1e-16
pfset Solver.Nonlinear.StepTol			 1e-12
pfset Solver.Nonlinear.Globalization		 LineSearch
pfset Solver.Linear.KrylovDimension		 70
pfset Solver.Linear.MaxRestart			 7

pfset Solver.Linear.Preconditioner                       PFMGOctree
pfset Solver.PrintSubsurf				 False
pfset Solver.Drop					 1E-20
pfset Solver.AbsTol					 1E-12

pfset Solver.PrintSaturation                            True
pfset Solver.PrintSubsurf                               True
pfset Solver.PrintPressure                              True 
pfset Solver.Nonlinear.PrintFlag                        LowVerbosity

pfset Solver.WriteSiloSubsurfData		        False	
pfset Solver.WriteSiloPressure				False
pfset Solver.WriteSiloSaturation		        False	
pfset Solver.WriteSiloMask			        False	
pfset Solver.WriteCLMBinary			        False	

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------
pfwritedb rurlaf 
