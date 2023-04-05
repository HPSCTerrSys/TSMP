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

pfset ComputationalGrid.DX			 50.
pfset ComputationalGrid.DY		     50.
pfset ComputationalGrid.DZ			 1.00

pfset ComputationalGrid.NX			 90
pfset ComputationalGrid.NY			 90
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
pfset GeomInput.Names	                     "domaininput"
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
pfset TimingInfo.BaseUnit		     0.0025
pfset TimingInfo.StartCount		     0.0
pfset TimingInfo.StartTime		     0.0
pfset TimingInfo.StopTime		     1.0025 
pfset TimeStep.Type			         Constant
pfset TimeStep.Value		     	 0.0125 
pfset TimingInfo.DumpInterval		 1

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
pfset Geom.Perm.Names			       "domain"

#pfset Geom.domain.Perm.Type            Constant
#pfset Geom.domain.Perm.Value           0.0036 

pfset Geom.domain.Perm.Type            "TurnBands"
pfset Geom.domain.Perm.LambdaX         1000.0
pfset Geom.domain.Perm.LambdaY         1000.0
pfset Geom.domain.Perm.LambdaZ         1.5
pfset Geom.domain.Perm.GeomMean        4.86
pfset Geom.domain.Perm.Sigma           1.43
pfset Geom.domain.Perm.NumLines        100
pfset Geom.domain.Perm.RZeta           5.0
pfset Geom.domain.Perm.KMax            100.0
pfset Geom.domain.Perm.DelK            1.0
pfset Geom.domain.Perm.Seed            1 
pfset Geom.domain.Perm.LogNormal       LogTruncated
pfset Geom.domain.Perm.StratType       Bottom
pfset Geom.domain.Perm.LowCutoff       0.0
pfset Geom.domain.Perm.HighCutoff      100.0


pfset Perm.TensorType			       TensorByGeom
pfset Geom.Perm.TensorByGeom.Names	   "domain"
pfset Geom.domain.Perm.TensorValX	   1.0
pfset Geom.domain.Perm.TensorValY	   1.0
pfset Geom.domain.Perm.TensorValZ	   1.0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------
pfset SpecificStorage.Type			       Constant
pfset SpecificStorage.GeomNames			   "domain"
pfset Geom.domain.SpecificStorage.Value	   1.0e-4

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------
pfset Phase.Names			         "water"
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
pfset Geom.Porosity.GeomNames               "domain"

pfset Geom.domain.Porosity.Type             Constant
pfset Geom.domain.Porosity.Value            0.451

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------
pfset Phase.RelPerm.Type                    VanGenuchten
pfset Phase.RelPerm.GeomNames               "domain"

pfset Geom.domain.RelPerm.Alpha             2.0
pfset Geom.domain.RelPerm.N                 2.0
pfset Geom.domain.RelPerm.NumSamplePoints   $VG_points
pfset Geom.domain.RelPerm.MinPressureHead   $VG_pmin 

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------
pfset Phase.Saturation.Type                 VanGenuchten
pfset Phase.Saturation.GeomNames            "domain"

pfset Geom.domain.Saturation.Alpha          2.0
pfset Geom.domain.Saturation.N              2.0
pfset Geom.domain.Saturation.SRes           0.1
pfset Geom.domain.Saturation.SSat           1.0

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

pfset Patch.x-upper.BCPressure.Type                   FluxConst
pfset Patch.x-upper.BCPressure.Cycle                  "constant"
pfset Patch.x-upper.BCPressure.alltime.Value          0.0

pfset Patch.y-lower.BCPressure.Type                   FluxConst
pfset Patch.y-lower.BCPressure.Cycle                  "constant"
pfset Patch.y-lower.BCPressure.alltime.Value          0.0

pfset Patch.y-upper.BCPressure.Type                   FluxConst
pfset Patch.y-upper.BCPressure.Cycle                  "constant"
pfset Patch.y-upper.BCPressure.alltime.Value          0.0

pfset Patch.z-lower.BCPressure.Type                   DirEquilRefPatch
pfset Patch.z-lower.BCPressure.Cycle                  "constant"
pfset Patch.z-lower.BCPressure.RefGeom                domain
pfset Patch.z-lower.BCPressure.RefPatch               z-upper
pfset Patch.z-lower.BCPressure.alltime.Value          -3.0

pfset Patch.z-upper.BCPressure.Type                   OverlandFlow
pfset Patch.z-upper.BCPressure.Cycle                  "constant"
pfset Patch.z-upper.BCPressure.alltime.Value           0.0



#  TOPOGRAPHY & SLOPES IN
#  BOTH X- & Y- DIRECTIONS
#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------
pfset TopoSlopesX.Type                    "Constant"
pfset TopoSlopesX.GeomNames               "domain"
pfset TopoSlopesX.Geom.domain.Value       0.0

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------
pfset TopoSlopesY.Type                    "Constant"
pfset TopoSlopesY.GeomNames               "domain"
pfset TopoSlopesY.Geom.domain.Value       0.0

#---------------------------------------------------------
# Mannings coefficient
#---------------------------------------------------------
pfset Mannings.Type			              "Constant"
pfset Mannings.GeomNames		          "domain"
pfset Mannings.Geom.domain.Value	      5.52e-6

#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------
pfset ICPressure.Type                    __pfl_ICPpressureType__
pfset ICPressure.GeomNames               domain
pfset Geom.domain.ICPressure.Value       __pfl_ICPpressureValue__
pfset Geom.domain.ICPressure.FileName    "__pfl_ICPpressureFileName__"
pfset Geom.domain.ICPressure.RefGeom     domain
pfset Geom.domain.ICPressure.RefPatch    z-upper
#
#pfset ICPressure.Type                    HydroStaticPatch
#pfset ICPressure.GeomNames               domain
#pfset Geom.domain.ICPressure.Value       -2.0
#pfset Geom.domain.ICPressure.RefGeom     domain
#pfset Geom.domain.ICPressure.RefPatch    z-upper

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------
pfset PhaseSources.water.Type			     Constant
pfset PhaseSources.water.GeomNames		     domain
pfset PhaseSources.water.Geom.domain.Value	 0.0

#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------
pfset KnownSolution				             NoKnownSolution

# Set solver parameters
#-----------------------------------------------------------------------------
pfset Solver					                       Richards
pfset Solver.MaxIter				                   10000

pfset Solver.TerrainFollowingGrid                      True

pfset Solver.Nonlinear.MaxIter			               100
pfset Solver.Nonlinear.ResidualTol		               1e-5
pfset Solver.Nonlinear.EtaChoice		               Walker1
pfset Solver.Nonlinear.EtaChoice		               EtaConstant
pfset Solver.Nonlinear.EtaValue			               0.001
pfset Solver.Nonlinear.UseJacobian		               False
pfset Solver.Nonlinear.DerivativeEpsilon	           1e-16
pfset Solver.Nonlinear.StepTol			               1e-12
pfset Solver.Nonlinear.Globalization		           LineSearch
pfset Solver.Linear.KrylovDimension		               30
pfset Solver.Linear.MaxRestart			               2

#pfset Solver.Linear.Preconditioner                     PFMG
pfset Solver.Linear.Preconditioner			           MGSemi
pfset Solver.Linear.Preconditioner.MGSemi.MaxIter	   1
pfset Solver.Linear.Preconditioner.MGSemi.MaxLevels	   10
pfset Solver.PrintSubsurf				               False
pfset Solver.Drop					                   1E-20
pfset Solver.AbsTol					                   1E-12

pfset Solver.PrintSaturation                           True
pfset Solver.PrintPressure                             True
pfset Solver.PrintSubsurf                              False
pfset Solver.Nonlinear.PrintFlag                       LowVerbosity

pfset Solver.WriteSiloSubsurfData		               False
pfset Solver.WriteSiloPressure				           False
pfset Solver.WriteSiloSaturation		               False
pfset Solver.WriteSiloMask			                   False
pfset Solver.WriteCLMBinary			                   False

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------
pfwritedb idiurnal
