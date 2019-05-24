# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*
#
#For normal soils, need to fix near saturation CPS
set VG_points 20000
set VG_pmin -50.0
#set VG_interpolation_method "Spline"

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

pfset ComputationalGrid.DX			 500.
pfset ComputationalGrid.DY		         500. 
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
pfset GeomInput.indinput.GeomNames  "v1 v2 v3 v4 v5 v6 v7 v8 v9 a1 a2 a3 a4 a5 a6"
pfset Geom.indinput.FileName  "__forcingdir__/rurSoilInd.pfb"

pfset GeomInput.v1.Value              1
pfset GeomInput.v2.Value              2
pfset GeomInput.v3.Value              3 
pfset GeomInput.v4.Value              4 
pfset GeomInput.v5.Value              5
pfset GeomInput.v6.Value              6
pfset GeomInput.v7.Value              7
pfset GeomInput.a1.Value              11
pfset GeomInput.a2.Value              12
pfset GeomInput.a3.Value              13
pfset GeomInput.a4.Value              14
pfset GeomInput.a5.Value              15
pfset GeomInput.a6.Value              16
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
pfset Cycle.Names                       "constant"
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat             -1

#  HYDROLOGICAL PARAMETERS
#Schaap and Leiz (1998), Soil Science
#  SETUP AND VALUES
#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------
pfset Geom.Perm.Names			"v1 v2 v3 v4 v5 v6 v7 v8 v9 a1 a2 a3 a4 a5 a6"

pfset Geom.v1.Perm.Type              Constant
pfset Geom.v1.Perm.Value             0.0328

pfset Geom.v2.Perm.Type              Constant
pfset Geom.v2.Perm.Value             0.0328

pfset Geom.v3.Perm.Type              Constant
pfset Geom.v3.Perm.Value             0.0328

pfset Geom.v4.Perm.Type              Constant
pfset Geom.v4.Perm.Value             0.0328

pfset Geom.v5.Perm.Type              Constant
pfset Geom.v5.Perm.Value             0.0328

pfset Geom.v6.Perm.Type              Constant
pfset Geom.v6.Perm.Value             0.0328

pfset Geom.v7.Perm.Type              Constant
pfset Geom.v7.Perm.Value             0.0328

pfset Geom.v8.Perm.Type              Constant
pfset Geom.v8.Perm.Value             0.0328

pfset Geom.v9.Perm.Type              Constant
pfset Geom.v9.Perm.Value             0.0328

pfset Geom.a1.Perm.Type              Constant
pfset Geom.a1.Perm.Value             0.0045

pfset Geom.a2.Perm.Type              Constant
pfset Geom.a2.Perm.Value             0.0045

pfset Geom.a3.Perm.Type              Constant
pfset Geom.a3.Perm.Value             0.0045

pfset Geom.a4.Perm.Type              Constant
pfset Geom.a4.Perm.Value             0.0045

pfset Geom.a5.Perm.Type              Constant
pfset Geom.a5.Perm.Value             0.0045

pfset Geom.a6.Perm.Type              Constant
pfset Geom.a6.Perm.Value             0.0045

pfset Perm.TensorType			 TensorByGeom
pfset Geom.Perm.TensorByGeom.Names	 "domain"
pfset Geom.domain.Perm.TensorValX	 1.0
pfset Geom.domain.Perm.TensorValY	 1.0
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
pfset Geom.Porosity.GeomNames         "v1 v2 v3 v4 v5 v6 v7 v8 v9 a1 a2 a3 a4 a5 a6" 

pfset Geom.v1.Porosity.Type          Constant
pfset Geom.v1.Porosity.Value         0.4386

pfset Geom.v2.Porosity.Type          Constant
pfset Geom.v2.Porosity.Value         0.4386

pfset Geom.v3.Porosity.Type          Constant
pfset Geom.v3.Porosity.Value         0.4386

pfset Geom.v4.Porosity.Type          Constant
pfset Geom.v4.Porosity.Value         0.4386

pfset Geom.v5.Porosity.Type          Constant
pfset Geom.v5.Porosity.Value         0.4386

pfset Geom.v6.Porosity.Type          Constant
pfset Geom.v6.Porosity.Value         0.4386

pfset Geom.v7.Porosity.Type          Constant
pfset Geom.v7.Porosity.Value         0.4386

pfset Geom.v8.Porosity.Type          Constant
pfset Geom.v8.Porosity.Value         0.4386

pfset Geom.v9.Porosity.Type          Constant
pfset Geom.v9.Porosity.Value         0.4386

pfset Geom.a1.Porosity.Type          Constant
pfset Geom.a1.Porosity.Value         0.4071

pfset Geom.a2.Porosity.Type          Constant
pfset Geom.a2.Porosity.Value         0.4071

pfset Geom.a3.Porosity.Type          Constant
pfset Geom.a3.Porosity.Value         0.4071

pfset Geom.a4.Porosity.Type          Constant
pfset Geom.a4.Porosity.Value         0.4071

pfset Geom.a5.Porosity.Type          Constant
pfset Geom.a5.Porosity.Value         0.4071

pfset Geom.a6.Porosity.Type          Constant
pfset Geom.a6.Porosity.Value         0.4071

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------
pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          "v1 v2 v3 v4 v5 v6 v7 v8 v9 a1 a2 a3 a4 a5 a6" 

pfset Geom.v1.RelPerm.Alpha         0.4597
pfset Geom.v1.RelPerm.N             1.4652
pfset Geom.v1.RelPerm.NumSamplePoints   $VG_points
pfset Geom.v1.RelPerm.MinPressureHead   $VG_pmin 

pfset Geom.v2.RelPerm.Alpha         0.4597
pfset Geom.v2.RelPerm.N             1.4652
pfset Geom.v2.RelPerm.NumSamplePoints   $VG_points
pfset Geom.v2.RelPerm.MinPressureHead   $VG_pmin

pfset Geom.v3.RelPerm.Alpha         0.4597
pfset Geom.v3.RelPerm.N             1.4652
pfset Geom.v3.RelPerm.NumSamplePoints   $VG_points
pfset Geom.v3.RelPerm.MinPressureHead   $VG_pmin

pfset Geom.v4.RelPerm.Alpha         0.4597
pfset Geom.v4.RelPerm.N             1.4652
pfset Geom.v4.RelPerm.NumSamplePoints   $VG_points
pfset Geom.v4.RelPerm.MinPressureHead   $VG_pmin

pfset Geom.v5.RelPerm.Alpha         0.4597
pfset Geom.v5.RelPerm.N             1.4652
pfset Geom.v5.RelPerm.NumSamplePoints   $VG_points
pfset Geom.v5.RelPerm.MinPressureHead   $VG_pmin

pfset Geom.v6.RelPerm.Alpha         0.4597
pfset Geom.v6.RelPerm.N             1.4652
pfset Geom.v6.RelPerm.NumSamplePoints   $VG_points
pfset Geom.v6.RelPerm.MinPressureHead   $VG_pmin

pfset Geom.v7.RelPerm.Alpha         0.4597
pfset Geom.v7.RelPerm.N             1.4652
pfset Geom.v7.RelPerm.NumSamplePoints   $VG_points
pfset Geom.v7.RelPerm.MinPressureHead   $VG_pmin

pfset Geom.v8.RelPerm.Alpha         0.4597
pfset Geom.v8.RelPerm.N             1.4652
pfset Geom.v8.RelPerm.NumSamplePoints   $VG_points
pfset Geom.v8.RelPerm.MinPressureHead   $VG_pmin

pfset Geom.v9.RelPerm.Alpha         0.4597
pfset Geom.v9.RelPerm.N             1.4652
pfset Geom.v9.RelPerm.NumSamplePoints   $VG_points
pfset Geom.v9.RelPerm.MinPressureHead   $VG_pmin

pfset Geom.a1.RelPerm.Alpha       0.39
pfset Geom.a1.RelPerm.N           1.4
pfset Geom.a1.RelPerm.NumSamplePoints   $VG_points
pfset Geom.a1.RelPerm.MinPressureHead   $VG_pmin 

pfset Geom.a2.RelPerm.Alpha       0.39
pfset Geom.a2.RelPerm.N           1.4
pfset Geom.a2.RelPerm.NumSamplePoints   $VG_points
pfset Geom.a2.RelPerm.MinPressureHead   $VG_pmin 

pfset Geom.a3.RelPerm.Alpha       0.39
pfset Geom.a3.RelPerm.N           1.4
pfset Geom.a3.RelPerm.NumSamplePoints   $VG_points
pfset Geom.a3.RelPerm.MinPressureHead   $VG_pmin 

pfset Geom.a4.RelPerm.Alpha       0.39
pfset Geom.a4.RelPerm.N           1.4
pfset Geom.a4.RelPerm.NumSamplePoints   $VG_points
pfset Geom.a4.RelPerm.MinPressureHead   $VG_pmin 

pfset Geom.a5.RelPerm.Alpha       0.39
pfset Geom.a5.RelPerm.N           1.4
pfset Geom.a5.RelPerm.NumSamplePoints   $VG_points
pfset Geom.a5.RelPerm.MinPressureHead   $VG_pmin 

pfset Geom.a6.RelPerm.Alpha       0.39
pfset Geom.a6.RelPerm.N           1.4
pfset Geom.a6.RelPerm.NumSamplePoints   $VG_points
pfset Geom.a6.RelPerm.MinPressureHead   $VG_pmin 

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------
pfset Phase.Saturation.Type              VanGenuchten
pfset Phase.Saturation.GeomNames         "v1 v2 v3 v4 v5 v6 v7 v8 v9 a1 a2 a3 a4 a5 a6"

pfset Geom.v1.Saturation.Alpha        0.4597
pfset Geom.v1.Saturation.N            1.4652
pfset Geom.v1.Saturation.SRes         0.21
pfset Geom.v1.Saturation.SSat         1.0

pfset Geom.v2.Saturation.Alpha        0.4597
pfset Geom.v2.Saturation.N            1.4652
pfset Geom.v2.Saturation.SRes         0.21
pfset Geom.v2.Saturation.SSat         1.0

pfset Geom.v3.Saturation.Alpha        0.4597
pfset Geom.v3.Saturation.N            1.4652
pfset Geom.v3.Saturation.SRes         0.21
pfset Geom.v3.Saturation.SSat         1.0

pfset Geom.v4.Saturation.Alpha        0.4597
pfset Geom.v4.Saturation.N            1.4652
pfset Geom.v4.Saturation.SRes         0.21
pfset Geom.v4.Saturation.SSat         1.0

pfset Geom.v5.Saturation.Alpha        0.4597
pfset Geom.v5.Saturation.N            1.4652
pfset Geom.v5.Saturation.SRes         0.21
pfset Geom.v5.Saturation.SSat         1.0

pfset Geom.v6.Saturation.Alpha        0.4597
pfset Geom.v6.Saturation.N            1.4652
pfset Geom.v6.Saturation.SRes         0.21
pfset Geom.v6.Saturation.SSat         1.0

pfset Geom.v7.Saturation.Alpha        0.4597
pfset Geom.v7.Saturation.N            1.4652
pfset Geom.v7.Saturation.SRes         0.21
pfset Geom.v7.Saturation.SSat         1.0

pfset Geom.v8.Saturation.Alpha        0.4597
pfset Geom.v8.Saturation.N            1.4652
pfset Geom.v8.Saturation.SRes         0.21
pfset Geom.v8.Saturation.SSat         1.0

pfset Geom.v9.Saturation.Alpha        0.4597
pfset Geom.v9.Saturation.N            1.4652
pfset Geom.v9.Saturation.SRes         0.21
pfset Geom.v9.Saturation.SSat         1.0

pfset Geom.a1.Saturation.Alpha        0.39
pfset Geom.a1.Saturation.N            1.4
pfset Geom.a1.Saturation.SRes         0.1
pfset Geom.a1.Saturation.SSat         1.0

pfset Geom.a2.Saturation.Alpha        0.39
pfset Geom.a2.Saturation.N            1.4
pfset Geom.a2.Saturation.SRes         0.1
pfset Geom.a2.Saturation.SSat         1.0

pfset Geom.a3.Saturation.Alpha        0.39
pfset Geom.a3.Saturation.N            1.4
pfset Geom.a3.Saturation.SRes         0.1
pfset Geom.a3.Saturation.SSat         1.0

pfset Geom.a4.Saturation.Alpha        0.39
pfset Geom.a4.Saturation.N            1.4
pfset Geom.a4.Saturation.SRes         0.1
pfset Geom.a4.Saturation.SSat         1.0

pfset Geom.a5.Saturation.Alpha        0.39
pfset Geom.a5.Saturation.N            1.4
pfset Geom.a5.Saturation.SRes         0.1
pfset Geom.a5.Saturation.SSat         1.0

pfset Geom.a6.Saturation.Alpha        0.39
pfset Geom.a6.Saturation.N            1.4
pfset Geom.a6.Saturation.SRes         0.1
pfset Geom.a6.Saturation.SSat         1.0

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

# Set solver parameters
#-----------------------------------------------------------------------------
pfset Solver					 Richards
pfset Solver.MaxIter				 10000

pfset Solver.TerrainFollowingGrid                True

pfset Solver.Nonlinear.MaxIter			 100
pfset Solver.Nonlinear.ResidualTol		 1e-5
pfset Solver.Nonlinear.EtaChoice		 Walker1
pfset Solver.Nonlinear.EtaChoice		 EtaConstant
pfset Solver.Nonlinear.EtaValue			 0.001
pfset Solver.Nonlinear.UseJacobian		 False
pfset Solver.Nonlinear.DerivativeEpsilon	 1e-16
pfset Solver.Nonlinear.StepTol			 1e-12
pfset Solver.Nonlinear.Globalization		 LineSearch
pfset Solver.Linear.KrylovDimension		 30
pfset Solver.Linear.MaxRestart			 2

pfset Solver.Linear.Preconditioner                       PFMG
#pfset Solver.Linear.Preconditioner			 MGSemi
pfset Solver.Linear.Preconditioner.MGSemi.MaxIter	 1
pfset Solver.Linear.Preconditioner.MGSemi.MaxLevels	 10
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
#pfrun default_single
#pfundist default_single
pfwritedb rurlaf 
