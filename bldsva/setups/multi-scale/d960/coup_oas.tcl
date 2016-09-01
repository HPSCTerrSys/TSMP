# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*
#
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
pfset ComputationalGrid.Lower.Z			 0.0

pfset ComputationalGrid.DX			 960.
pfset ComputationalGrid.DY		         960. 
pfset ComputationalGrid.DZ			 1.0 

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
pfset GeomInput.Names	                     "solidinput"
pfset GeomInput.solidinput.InputType         SolidFile 
pfset GeomInput.solidinput.GeomNames         domain
pfset GeomInput.solidinput.FileName          "__forcingdir__/geom.pfsol"

pfset Geom.domain.Patches	             "top bottom perimeter"

#-----------------------------------------------------------------------------
# VARIABLE dz ASSIGNMENTS
#-----------------------------------------------------------------------------
pfset Solver.Nonlinear.VariableDz            True 
pfset dzScale.GeomNames                      domain
pfset dzScale.Type                           nzList
pfset dzScale.nzListNumber                   [pfget ComputationalGrid.NZ] 
pfset Cell.0.dzScale.Value                   1.35
pfset Cell.1.dzScale.Value                   1.35
pfset Cell.2.dzScale.Value                  1.35
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

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names "constant rainrec"
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat             -1
# rainfall and recession time periods are defined here
# rain for 2 days, recession for 
pfset Cycle.rainrec.Names                 "rain rec"
pfset Cycle.rainrec.rain.Length           16
pfset Cycle.rainrec.rec.Length            222
pfset Cycle.rainrec.Repeat                -1

#  HYDROLOGICAL PARAMETERS
#Schaap and Leiz (1998), Soil Science
#  SETUP AND VALUES
#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------
pfset Geom.Perm.Names			 "domain"
pfset Geom.domain.Perm.Type		 Constant
pfset Geom.domain.Perm.Value		 0.0034
#pfset Geom.domain.Perm.Value		 0.0

#pfset Geom.domain.Perm.Type              "TurnBands"
#pfset Geom.domain.Perm.LambdaX           100
#pfset Geom.domain.Perm.LambdaY           10
#pfset Geom.domain.Perm.LambdaZ           1
#pfset Geom.domain.Perm.GeomMean          0.01
#pfset Geom.domain.Perm.Sigma             0.5

#pfset Geom.domain.Perm.NumLines          100
#pfset Geom.domain.Perm.RZeta             5.0
#pfset Geom.domain.Perm.KMax              10.0
#pfset Geom.domain.Perm.DelK              0.2
#pfset Geom.domain.Perm.Seed              333
#pfset Geom.domain.Perm.LogNormal         Log
#pfset Geom.domain.Perm.StratType         Bottom

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
pfset Geom.Porosity.GeomNames		 domain
pfset Geom.domain.Porosity.Type		 Constant
pfset Geom.domain.Porosity.Value	 0.4449

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------
pfset Phase.RelPerm.Type		 VanGenuchten
pfset Phase.RelPerm.GeomNames		 "domain"
pfset Geom.domain.RelPerm.Alpha            2.1 
pfset Geom.domain.RelPerm.N                2.0 
#---------------------------------------------------------
# Saturation
#---------------------------------------------------------
pfset Phase.Saturation.Type		 VanGenuchten
pfset Phase.Saturation.GeomNames	 domain
pfset Geom.domain.Saturation.Alpha	 2.1 
pfset Geom.domain.Saturation.N		 2.0
pfset Geom.domain.Saturation.SRes	 0.15
pfset Geom.domain.Saturation.SSat	 1.0

#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names				 ""

#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------

pfset BCPressure.PatchNames [pfget Geom.domain.Patches]

pfset Patch.top.BCPressure.Type                      OverlandFlow 
pfset Patch.top.BCPressure.Cycle                     "constant"
pfset Patch.top.BCPressure.alltime.Value             0.0

pfset Patch.bottom.BCPressure.Type                   FluxConst
pfset Patch.bottom.BCPressure.Cycle                  "constant"
pfset Patch.bottom.BCPressure.alltime.Value          0.0

pfset Patch.perimeter.BCPressure.Type                FluxConst
pfset Patch.perimeter.BCPressure.Cycle               "constant"
pfset Patch.perimeter.BCPressure.alltime.Value       0.0

#  TOPOGRAPHY & SLOPES IN
#  BOTH X- & Y- DIRECTIONS
#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

pfset TopoSlopesX.Type			 "PFBFile"
pfset TopoSlopesX.GeomNames		 "domain"
pfset TopoSlopesX.FileName		 "__forcingdir__/xslope960sc.pfb"

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------

pfset TopoSlopesY.Type			 "PFBFile"
pfset TopoSlopesY.GeomNames		 "domain"
pfset TopoSlopesY.FileName		 "__forcingdir__/yslope960sc.pfb"

#---------------------------------------------------------
# Mannings coefficient
#---------------------------------------------------------
pfset Mannings.Type			 "Constant"
pfset Mannings.GeomNames		 "domain"
pfset Mannings.Geom.domain.Value	 5.56e-5

#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------
pfset ICPressure.Type                    "PFBFile"
pfset ICPressure.GeomNames               "domain"
pfset Geom.domain.ICPressure.FileName    "__forcingdir__/ptb960sc_gwt05m.pfb"
#for restart purpose
pfset Geom.domain.ICPressure.FileName   "__pfl_ICPpressureFileName__"

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
pfset Solver                                     Richards
pfset Solver.MaxIter                             100000

pfset Solver.Nonlinear.PrintFlag                 LowVerbosity
pfset Solver.TerrainFollowingGrid                True

pfset Solver.Nonlinear.MaxIter                   500
pfset Solver.Nonlinear.ResidualTol               1e-4
pfset Solver.Nonlinear.EtaChoice                 Walker1
pfset Solver.Nonlinear.EtaChoice                 EtaConstant
pfset Solver.Nonlinear.EtaValue                  0.001
pfset Solver.Nonlinear.UseJacobian               False 
pfset Solver.Nonlinear.DerivativeEpsilon         1e-16
pfset Solver.Nonlinear.StepTol                   1e-20
pfset Solver.Nonlinear.Globalization             LineSearch
pfset Solver.Linear.KrylovDimension              30
pfset Solver.Linear.MaxRestart                   5 
#
pfset Solver.Linear.Preconditioner		         PFMGOctree	
pfset Solver.PrintSubsurf				 False
pfset Solver.Drop					 1E-20
pfset Solver.AbsTol					 1E-12

pfset Solver.PrintSaturation                            True 
pfset Solver.PrintSubsurf                               False
pfset Solver.PrintPressure                              True 
pfset Solver.PrintSubsurf                               False

pfset Solver.WriteSiloSubsurfData		        False	
pfset Solver.WriteSiloPressure				False
pfset Solver.WriteSiloSaturation		        False	
pfset Solver.WriteSiloMask			        False	
pfset Solver.WriteCLMBinary			        False	

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------
pfwritedb rurlaf 
