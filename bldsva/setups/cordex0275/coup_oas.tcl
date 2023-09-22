# Import the ParFlow TCL package
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*

#-------------------------------------------------------------------------------
pfset FileVersion 4
#-------------------------------------------------------------------------------
pfset Process.Topology.P __nprocx_pfl_bldsva__
pfset Process.Topology.Q __nprocy_pfl_bldsva__
pfset Process.Topology.R 1

# THE COMPUTATIONAL GRID IS A (BOX) THT CONTAINS THE MAIN PROBLEM. THIS CAN 
# EITHER BE EXACTLY THE SIZE OF THE PROBLEM OR LARGER. A BOX GEOMETRY IN PARFLOW 
# CAN BE ASIGNED BY EITHER SPECIFYING COORDINATES FOR TWO CORNERS OF THE BOX OR 
# GRID SIZE AND NUMBER OF CELLS IN X,Y, AND Z.
#-------------------------------------------------------------------------------
# Computational Grid: It Defines The Grid Resolutions within The Domain
#-------------------------------------------------------------------------------
pfset ComputationalGrid.Lower.X			 0.0
pfset ComputationalGrid.Lower.Y			 0.0
pfset ComputationalGrid.Lower.Z			 0.0

pfset ComputationalGrid.DX			 3000.
pfset ComputationalGrid.DY		         3000. 
pfset ComputationalGrid.DZ			 2.00

pfset ComputationalGrid.NX			 __ngpflx_bldsva__
pfset ComputationalGrid.NY			 __ngpfly_bldsva__
pfset ComputationalGrid.NZ			 15 

# DOMAIN GEOMETRY IS THE (EXACTLY) OUTER DOMAIN OR BOUNDARY OF THE MODEL 
# PROBLEM. IT HAS TO BE CONTAINED WITHIN THE COMPUTATIONAL GRID (i.e. OR HAS THE 
# SAME SIZE OF IT). THE DOMAIN GEOMETRY COULD BE A BOX OR IT COULD BE A 
# SOLID-FILE. BOUNDARY CONDITIONS ARE ASSIGNED TO THE DOMAIN SIDES WITH 
# SOMETHING CALLED (PATCHES) IN THE TCL-SCRIPT. A BOX HAS SIX (6) SIDES AND (6) 
# PATCHES WHILE A SOLID-FILE CAN HAVE ANY NUMBER OF PATCHES.
#-------------------------------------------------------------------------------
# Domain
#-------------------------------------------------------------------------------
pfset Domain.GeomName                            domain

#-------------------------------------------------------------------------------
# Domain Geometry Input 
#-------------------------------------------------------------------------------

 #pfset GeomInput.Names                 "solidinput"
 pfset GeomInput.Names                 "solidinput indi_input"
 pfset GeomInput.solidinput.InputType  SolidFile
 pfset GeomInput.solidinput.GeomNames  domain
 pfset GeomInput.solidinput.FileName   __pfl_solidinput_filename__
 pfset Geom.domain.Patches             "top bottom perimeter"

#-------------------------------------------------------------------------------
# VARIABLE dz ASSIGNMENTS
#-------------------------------------------------------------------------------
pfset Solver.Nonlinear.VariableDz    True
pfset dzScale.GeomNames              domain
pfset dzScale.Type                   nzList

# NOTE each cell depth is dz*dzScale!!!!!!
pfset dzScale.nzListNumber           15
pfset Cell.0.dzScale.Value           9.00
pfset Cell.1.dzScale.Value           7.50
pfset Cell.2.dzScale.Value           5.0
pfset Cell.3.dzScale.Value           5.0
pfset Cell.4.dzScale.Value           2.0
pfset Cell.5.dzScale.Value           0.5
pfset Cell.6.dzScale.Value           0.35
pfset Cell.7.dzScale.Value           0.25
pfset Cell.8.dzScale.Value           0.15
pfset Cell.9.dzScale.Value           0.10
pfset Cell.10.dzScale.Value          0.065
pfset Cell.11.dzScale.Value          0.035
pfset Cell.12.dzScale.Value          0.025
pfset Cell.13.dzScale.Value          0.015
pfset Cell.14.dzScale.Value          0.01


#-------------------------------------------------------------------------------
# Indicator Geometry Input
#-------------------------------------------------------------------------------

pfset GeomInput.indi_input.InputType IndicatorField
pfset GeomInput.indi_input.GeomNames "TC01 TC02 TC03 TC04 TC05 TC06 TC07 TC08 TC09 TC10 TC11 TC12 Allv BGR1 BGR2 BGR3 BGR4 BGR5 BGR6 Lake Sea"
pfset Geom.indi_input.FileName EUR-0275_TSMP_FZJ-IBG3_CLMPFLDomain_1592x1544_INDICATOR_regridded_rescaled_SoilGrids250-v2017_BGR3_alv.pfb

# USDA classes calculated based on SoilGrids
# 12 soil classes TC01-TC12 
# sand
pfset GeomInput.TC01.Value 1
# loamy sand
pfset GeomInput.TC02.Value 2
# sandy loam
pfset GeomInput.TC03.Value 3
# loam
pfset GeomInput.TC04.Value 4
# silt loam
pfset GeomInput.TC05.Value 5
# silt
pfset GeomInput.TC06.Value 6
# sandy clay loam
pfset GeomInput.TC07.Value 7
# clay loam
pfset GeomInput.TC08.Value 8
# silty clay loam
pfset GeomInput.TC09.Value 9
# sandy clay
pfset GeomInput.TC10.Value 10
# silty clay
pfset GeomInput.TC11.Value 11
# clay
pfset GeomInput.TC12.Value 12

# Alluvium burned in soil
pfset GeomInput.Allv.Value 13

# IHME plus GLHYMPS dataset
# Geology BGR1-BGR6 from 'highe' permeabilities to 'low' permabilities
pfset GeomInput.BGR1.Value 14
pfset GeomInput.BGR2.Value 15
pfset GeomInput.BGR3.Value 16
pfset GeomInput.BGR4.Value 17
pfset GeomInput.BGR5.Value 18
pfset GeomInput.BGR6.Value 19

# Lake
pfset GeomInput.Lake.Value 21
# Sea
pfset GeomInput.Sea.Value 22

#TIME SETUP
#-------------------------------------------------------------------------------
# Setup timing info
#-------------------------------------------------------------------------------
pfset TimingInfo.BaseUnit                __base_pfl__
pfset TimingInfo.StartCount              __start_cnt_pfl__
pfset TimingInfo.StartTime               0.0
pfset TimingInfo.StopTime                __stop_pfl_bldsva__
pfset TimeStep.Type                      Constant
pfset TimeStep.Value                     __dt_pfl_bldsva__
pfset TimingInfo.DumpInterval            __dump_pfl_interval__

# Time Cycles
#-------------------------------------------------------------------------------
pfset Cycle.Names "constant"
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat             -1

# HYDROLOGICAL PARAMETERS
# Schaap and Leiz (1998), Soil Science
# SETUP AND VALUES
#-------------------------------------------------------------------------------
# Perm
#-------------------------------------------------------------------------------
pfset Geom.Perm.Names              "domain TC01 TC02 TC03 TC04 TC05 TC06 TC07 TC08 TC09 TC10 TC11 TC12 Allv BGR1 BGR2 BGR3 BGR4 BGR5 BGR6 Lake Sea"

pfset Geom.domain.Perm.Type        Constant
pfset Geom.domain.Perm.Value       0.1

# USDA classes calculated based on SoilGrids
# Values from ROSATTA project:
# https://www.ars.usda.gov/pacific-west-area/riverside-ca/agricultural-water-efficiency-and-salinity-reSearch-unit/docs/model/rosetta-class-average-hydraulic-parameters/
pfset Geom.TC01.Perm.Type           Constant
pfset Geom.TC01.Perm.Value          0.267787

pfset Geom.TC02.Perm.Type           Constant
pfset Geom.TC02.Perm.Value          0.043832

pfset Geom.TC03.Perm.Type           Constant
pfset Geom.TC03.Perm.Value          0.015951

pfset Geom.TC04.Perm.Type           Constant
pfset Geom.TC04.Perm.Value          0.005021

pfset Geom.TC05.Perm.Type           Constant
pfset Geom.TC05.Perm.Value          0.0076

pfset Geom.TC06.Perm.Type           Constant
pfset Geom.TC06.Perm.Value          0.01823

pfset Geom.TC07.Perm.Type           Constant
pfset Geom.TC07.Perm.Value          0.005493

pfset Geom.TC08.Perm.Type           Constant
pfset Geom.TC08.Perm.Value          0.00341

pfset Geom.TC09.Perm.Type           Constant
pfset Geom.TC09.Perm.Value          0.004632

pfset Geom.TC10.Perm.Type           Constant
pfset Geom.TC10.Perm.Value          0.004729

pfset Geom.TC11.Perm.Type           Constant
pfset Geom.TC11.Perm.Value          0.004007

pfset Geom.TC12.Perm.Type           Constant
pfset Geom.TC12.Perm.Value          0.006149

# Alluvium burned in soil
pfset Geom.Allv.Perm.Type          Constant
pfset Geom.Allv.Perm.Value         0.1

# IHME plus GLHYMPS geology where permeability ranges are set according to
# ABE setting (from Wendy)
pfset Geom.BGR1.Perm.Type            Constant
pfset Geom.BGR1.Perm.Value           0.1

pfset Geom.BGR2.Perm.Type            Constant
pfset Geom.BGR2.Perm.Value           0.05

pfset Geom.BGR3.Perm.Type            Constant
pfset Geom.BGR3.Perm.Value           0.01

pfset Geom.BGR4.Perm.Type            Constant
pfset Geom.BGR4.Perm.Value           0.005

pfset Geom.BGR5.Perm.Type            Constant
pfset Geom.BGR5.Perm.Value           0.001

pfset Geom.BGR6.Perm.Type            Constant
pfset Geom.BGR6.Perm.Value           0.0005

# Lake
pfset Geom.Lake.Perm.Type          Constant
pfset Geom.Lake.Perm.Value         1.0e-5
# Sea
pfset Geom.Sea.Perm.Type           Constant
pfset Geom.Sea.Perm.Value          1.0e-5

pfset Perm.TensorType			 TensorByGeom
pfset Geom.Perm.TensorByGeom.Names	 "domain"
pfset Geom.domain.Perm.TensorValX	 1000.0 
pfset Geom.domain.Perm.TensorValY	 1000.0
pfset Geom.domain.Perm.TensorValZ	 1.0

#-------------------------------------------------------------------------------
# Specific Storage
#-------------------------------------------------------------------------------
pfset SpecificStorage.Type			 Constant
pfset SpecificStorage.GeomNames			 "domain"
pfset Geom.domain.SpecificStorage.Value		 1.0e-4

#-------------------------------------------------------------------------------
# Phases
#-------------------------------------------------------------------------------
pfset Phase.Names			 "water"
pfset Phase.water.Density.Type		 Constant
pfset Phase.water.Density.Value		 1.0
pfset Phase.water.Viscosity.Type	 Constant
pfset Phase.water.Viscosity.Value	 1.0

#-------------------------------------------------------------------------------
# Gravity
#-------------------------------------------------------------------------------
pfset Gravity				 1.0

#-------------------------------------------------------------------------------
# Contaminants
#-------------------------------------------------------------------------------
pfset Contaminants.Names		 ""

#-------------------------------------------------------------------------------
# Retardation
#-------------------------------------------------------------------------------
pfset Geom.Retardation.GeomNames	 ""

#-------------------------------------------------------------------------------
# Porosity
#-------------------------------------------------------------------------------
pfset Geom.Porosity.GeomNames          "domain TC01 TC02 TC03 TC04 TC05 TC06 TC07 TC08 TC09 TC10 TC11 TC12"

pfset Geom.domain.Porosity.Type        Constant
pfset Geom.domain.Porosity.Value       0.4

# USDA classes calculated based on SoilGrids
# Values from ROSATTA project:
# https://www.ars.usda.gov/pacific-west-area/riverside-ca/agricultural-water-efficiency-and-salinity-reSearch-unit/docs/model/rosetta-class-average-hydraulic-parameters/
pfset Geom.TC01.Porosity.Type            Constant
pfset Geom.TC01.Porosity.Value           0.375

pfset Geom.TC02.Porosity.Type            Constant
pfset Geom.TC02.Porosity.Value           0.39

pfset Geom.TC03.Porosity.Type            Constant
pfset Geom.TC03.Porosity.Value           0.387

pfset Geom.TC04.Porosity.Type            Constant
pfset Geom.TC04.Porosity.Value           0.399

pfset Geom.TC05.Porosity.Type            Constant
pfset Geom.TC05.Porosity.Value           0.439

pfset Geom.TC06.Porosity.Type            Constant
pfset Geom.TC06.Porosity.Value           0.489

pfset Geom.TC07.Porosity.Type            Constant
pfset Geom.TC07.Porosity.Value           0.384

pfset Geom.TC08.Porosity.Type            Constant
pfset Geom.TC08.Porosity.Value           0.442

pfset Geom.TC09.Porosity.Type            Constant
pfset Geom.TC09.Porosity.Value           0.482

pfset Geom.TC10.Porosity.Type            Constant
pfset Geom.TC10.Porosity.Value           0.385

pfset Geom.TC11.Porosity.Type            Constant
pfset Geom.TC11.Porosity.Value           0.481

pfset Geom.TC12.Porosity.Type            Constant
pfset Geom.TC12.Porosity.Value           0.459

#-------------------------------------------------------------------------------
# Relative Permeability
#-------------------------------------------------------------------------------
pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          "domain TC01 TC02 TC03 TC04 TC05 TC06 TC07 TC08 TC09 TC10 TC11 TC12"

pfset Geom.domain.RelPerm.Alpha        2.0
pfset Geom.domain.RelPerm.N            4.

# USDA classes calculated based on SoilGrids
# Values from ROSATTA project:
# https://www.ars.usda.gov/pacific-west-area/riverside-ca/agricultural-water-efficiency-and-salinity-reSearch-unit/docs/model/rosetta-class-average-hydraulic-parameters/
# rounded to second digit after comma
pfset Geom.TC01.RelPerm.Alpha          3.52
pfset Geom.TC01.RelPerm.N              4.18

pfset Geom.TC02.RelPerm.Alpha          3.48
pfset Geom.TC02.RelPerm.N              2.75

pfset Geom.TC03.RelPerm.Alpha          2.67
pfset Geom.TC03.RelPerm.N              2.45

pfset Geom.TC04.RelPerm.Alpha          1.11
pfset Geom.TC04.RelPerm.N              2.47

pfset Geom.TC05.RelPerm.Alpha          0.51
pfset Geom.TC05.RelPerm.N              2.66

pfset Geom.TC06.RelPerm.Alpha          0.66
pfset Geom.TC06.RelPerm.N              2.68
pfset Geom.TC07.RelPerm.Alpha          2.11
pfset Geom.TC07.RelPerm.N              2.33

pfset Geom.TC08.RelPerm.Alpha          1.58
pfset Geom.TC08.RelPerm.N              2.42

pfset Geom.TC09.RelPerm.Alpha          0.84
pfset Geom.TC09.RelPerm.N              2.52

pfset Geom.TC10.RelPerm.Alpha          3.34
pfset Geom.TC10.RelPerm.N              2.21

pfset Geom.TC11.RelPerm.Alpha          1.62
pfset Geom.TC11.RelPerm.N              2.32

pfset Geom.TC12.RelPerm.Alpha          1.50
pfset Geom.TC12.RelPerm.N              2.25

#-------------------------------------------------------------------------------
# Saturation
#-------------------------------------------------------------------------------
pfset Phase.Saturation.Type              VanGenuchten
pfset Phase.Saturation.GeomNames         "domain TC01 TC02 TC03 TC04 TC05 TC06 TC07 TC08 TC09 TC10 TC11 TC12"

pfset Geom.domain.Saturation.Alpha        2.0
pfset Geom.domain.Saturation.N            4.
pfset Geom.domain.Saturation.SRes         0.1
pfset Geom.domain.Saturation.SSat         1.0

# USDA classes calculated based on SoilGrids
# Values from ROSATTA project:
# https://www.ars.usda.gov/pacific-west-area/riverside-ca/agricultural-water-efficiency-and-salinity-reSearch-unit/docs/model/rosetta-class-average-hydraulic-parameters/
# rounded to second digit after comma
pfset Geom.TC01.Saturation.Alpha       3.52
pfset Geom.TC01.Saturation.N           4.18
pfset Geom.TC01.Saturation.SRes        0.14
pfset Geom.TC01.Saturation.SSat        1.00

pfset Geom.TC02.Saturation.Alpha       3.48
pfset Geom.TC02.Saturation.N           2.75
pfset Geom.TC02.Saturation.SRes        0.13
pfset Geom.TC02.Saturation.SSat        1.00

pfset Geom.TC03.Saturation.Alpha       2.67
pfset Geom.TC03.Saturation.N           2.45
pfset Geom.TC03.Saturation.SRes        0.10
pfset Geom.TC03.Saturation.SSat        1.00

pfset Geom.TC04.Saturation.Alpha       1.11
pfset Geom.TC04.Saturation.N           2.47
pfset Geom.TC04.Saturation.SRes        0.15
pfset Geom.TC04.Saturation.SSat        1.00

pfset Geom.TC05.Saturation.Alpha       0.51
pfset Geom.TC05.Saturation.N           2.66
pfset Geom.TC05.Saturation.SRes        0.15
pfset Geom.TC05.Saturation.SSat        1.00

pfset Geom.TC06.Saturation.Alpha       0.66
pfset Geom.TC06.Saturation.N           2.68
pfset Geom.TC06.Saturation.SRes        0.10
pfset Geom.TC06.Saturation.SSat        1.00

pfset Geom.TC07.Saturation.Alpha       2.11
pfset Geom.TC07.Saturation.N           2.33
pfset Geom.TC07.Saturation.SRes        0.16
pfset Geom.TC07.Saturation.SSat        1.00

pfset Geom.TC08.Saturation.Alpha       1.58
pfset Geom.TC08.Saturation.N           2.42
pfset Geom.TC08.Saturation.SRes        0.18
pfset Geom.TC08.Saturation.SSat        1.00

pfset Geom.TC09.Saturation.Alpha       0.84
pfset Geom.TC09.Saturation.N           2.52
pfset Geom.TC09.Saturation.SRes        0.19
pfset Geom.TC09.Saturation.SSat        1.00

pfset Geom.TC10.Saturation.Alpha       3.34
pfset Geom.TC10.Saturation.N           2.21
pfset Geom.TC10.Saturation.SRes        0.30
pfset Geom.TC10.Saturation.SSat        1.00

pfset Geom.TC11.Saturation.Alpha       1.62
pfset Geom.TC11.Saturation.N           2.32
pfset Geom.TC11.Saturation.SRes        0.23
pfset Geom.TC11.Saturation.SSat        1.00

pfset Geom.TC12.Saturation.Alpha       1.50
pfset Geom.TC12.Saturation.N           2.25
pfset Geom.TC12.Saturation.SRes        0.21
pfset Geom.TC12.Saturation.SSat        1.00

#-------------------------------------------------------------------------------
# Wells
#-------------------------------------------------------------------------------
pfset Wells.Names				 ""

#-------------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-------------------------------------------------------------------------------

pfset BCPressure.PatchNames                     [pfget Geom.domain.Patches]

pfset Patch.top.BCPressure.Type                 OverlandFlow
pfset Patch.top.BCPressure.Cycle                "constant"
pfset Patch.top.BCPressure.alltime.Value        0.0

pfset Patch.bottom.BCPressure.Type              FluxConst
pfset Patch.bottom.BCPressure.Cycle             "constant"
pfset Patch.bottom.BCPressure.alltime.Value     0.0

pfset Patch.perimeter.BCPressure.Type           FluxConst
pfset Patch.perimeter.BCPressure.Cycle          "constant"
pfset Patch.perimeter.BCPressure.alltime.Value  0.0


# Dirichlet BC
 pfset Patch.perimeter.BCPressure.Type          DirEquilRefPatch
 pfset Patch.perimeter.BCPressure.Cycle         "constant"
 pfset Patch.perimeter.BCPressure.RefGeom       domain
 pfset Patch.perimeter.BCPressure.RefPatch      top
 pfset Patch.perimeter.BCPressure.alltime.Value -0.05

#  TOPOGRAPHY & SLOPES IN
#  BOTH X- & Y- DIRECTIONS
#-------------------------------------------------------------------------------
# Topo slopes in x-direction
#-------------------------------------------------------------------------------
pfset TopoSlopesX.Type			 "PFBFile"
pfset TopoSlopesX.GeomNames		 "domain"
pfset TopoSlopesX.FileName		 "EUR-0275_TSMP_FZJ-IBG3_CLMPFLDomain_1592x1544_XSLOPE_TPS_HydroRIVER_sea_streams_corr.pfb"

#-------------------------------------------------------------------------------
# Topo slopes in y-direction
#-------------------------------------------------------------------------------
pfset TopoSlopesY.Type			 "PFBFile"
pfset TopoSlopesY.GeomNames		 "domain"
pfset TopoSlopesY.FileName		 "EUR-0275_TSMP_FZJ-IBG3_CLMPFLDomain_1592x1544_YSLOPE_TPS_HydroRIVER_sea_streams_corr.pfb"

#-------------------------------------------------------------------------------
# Mannings coefficient
#-------------------------------------------------------------------------------
pfset Mannings.Type			 "Constant"
pfset Mannings.GeomNames		 "domain"
pfset Mannings.Geom.domain.Value	 5.5e-5

#-------------------------------------------------------------------------------
# Phase sources:
#-------------------------------------------------------------------------------
pfset PhaseSources.water.Type			 Constant
pfset PhaseSources.water.GeomNames		 domain
pfset PhaseSources.water.Geom.domain.Value	 0.0

#-------------------------------------------------------------------------------
# Exact solution specification for error calculations
#-------------------------------------------------------------------------------
pfset KnownSolution				                      NoKnownSolution

#-------------------------------------------------------------------------------
# Set solver parameters
#-------------------------------------------------------------------------------
pfset Solver					                          Richards
pfset Solver.MaxIter				                    100000

pfset Solver.TerrainFollowingGrid               True

pfset Solver.Nonlinear.MaxIter			            400
pfset Solver.Nonlinear.ResidualTol		          1e-4
pfset Solver.Nonlinear.EtaChoice		            Walker1
pfset Solver.Nonlinear.UseJacobian		          True
pfset Solver.Nonlinear.DerivativeEpsilon	      1e-16
pfset Solver.Nonlinear.StepTol			            1e-16
pfset Solver.Nonlinear.Globalization		        LineSearch
pfset Solver.Linear.KrylovDimension		          30
pfset Solver.Linear.MaxRestart			            8
pfset Solver.MaxConvergenceFailures             8

pfset Solver.Linear.Preconditioner              PFMG

pfset Solver.PrintPressure                      False
pfset Solver.PrintSaturation                    False
pfset Solver.PrintSubsurfData                   True
pfset Solver.PrintLSMSink                       False
pfset Solver.PrintEvapTransSum                  False
pfset Solver.PrintOverlandSum                   False
#
pfset NetCDF.NumStepsPerFile                    96
pfset NetCDF.WritePressure                      True
pfset NetCDF.WriteSaturation                    True
pfset NetCDF.WriteMannings                      True
pfset NetCDF.WriteSlopes                        True
# Set NetCDF.WriteSubsurface and NetCDF.WriteMask to False, as those are coverd 
# by Solver.PrintSubsurfData and vanGenuchten parameter are coverd by 
# Solver.PrintSubsurfData as well, but not by NetCDf output.
pfset NetCDF.WriteSubsurface                    False
pfset NetCDF.WriteMask                          False
pfset NetCDF.WriteDZMultiplier                  True
pfset NetCDF.WriteEvapTrans                     True
#
pfset Solver.WriteSiloSubsurfData		            False
pfset Solver.WriteSiloPressure				          False
pfset Solver.WriteSiloSaturation		            False	
pfset Solver.WriteSiloMask			                False	
pfset Solver.WriteCLMBinary			                False	

#-------------------------------------------------------------------------------
# Initial conditions: water pressure
#-------------------------------------------------------------------------------
pfset ICPressure.Type                           HydroStaticPatch
pfset ICPressure.GeomNames                      domain
pfset Geom.domain.ICPressure.Value              -0.2
pfset Geom.domain.ICPressure.RefGeom            domain
pfset Geom.domain.ICPressure.RefPatch           top

#pfset ICPressure.Type                    NCFile
#pfset ICPressure.GeomNames               domain
#pfset Geom.domain.ICPressure.FileName    "__ICPressure__"
#pfdist    "__ICPressure__"
#pfset Geom.domain.ICPressure.RefGeom     domain
#pfset Geom.domain.ICPressure.RefPatch    top

#-------------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-------------------------------------------------------------------------------
pfwritedb cordex0.0275
