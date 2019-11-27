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
pfset ComputationalGrid.Lower.Z			 0.0

pfset ComputationalGrid.DX			 12500.
pfset ComputationalGrid.DY		         12500. 
pfset ComputationalGrid.DZ			 2.00

pfset ComputationalGrid.NX			 __ngpflx_bldsva__
pfset ComputationalGrid.NY			 __ngpfly_bldsva__
pfset ComputationalGrid.NZ			 15 

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

 #pfset GeomInput.Names                 "solidinput"
 pfset GeomInput.Names                 "solidinput indi_input"
 pfset GeomInput.solidinput.InputType  SolidFile
 pfset GeomInput.solidinput.GeomNames  domain
 pfset GeomInput.solidinput.FileName   __pfl_solidinput_filename__
 pfset Geom.domain.Patches             "top bottom perimeter"

#-----------------------------------------------------------------------------
# VARIABLE dz ASSIGNMENTS
#-----------------------------------------------------------------------------
pfset Solver.Nonlinear.VariableDz    True
pfset dzScale.GeomNames              domain
pfset dzScale.Type                   nzList


pfset dzScale.nzListNumber           15
pfset Cell.0.dzScale.Value           7.50
pfset Cell.1.dzScale.Value           7.50
pfset Cell.2.dzScale.Value           5.0
pfset Cell.3.dzScale.Value           5.0
pfset Cell.4.dzScale.Value           2.0
pfset Cell.5.dzScale.Value           0.50
#pfset Cell.5.dzScale.Value           1.00
pfset Cell.6.dzScale.Value           0.350
#pfset Cell.6.dzScale.Value           0.70
pfset Cell.7.dzScale.Value           0.25
#pfset Cell.7.dzScale.Value           0.50
pfset Cell.8.dzScale.Value           0.15
#pfset Cell.8.dzScale.Value           0.30
pfset Cell.9.dzScale.Value           0.10
#pfset Cell.9.dzScale.Value           0.20
pfset Cell.10.dzScale.Value          0.065
#pfset Cell.10.dzScale.Value          0.13
pfset Cell.11.dzScale.Value          0.035
#pfset Cell.11.dzScale.Value          0.07
pfset Cell.12.dzScale.Value          0.025
#pfset Cell.12.dzScale.Value          0.05
pfset Cell.13.dzScale.Value          0.015
#pfset Cell.13.dzScale.Value          0.03
pfset Cell.14.dzScale.Value          0.01
#pfset Cell.14.dzScale.Value          0.02


#----------------------------------------------------------------
# Indicator Geometry Input
#----------------------------------------------------------------

pfset GeomInput.indi_input.InputType IndicatorField
pfset GeomInput.indi_input.GeomNames "W1 W2 W3 W4 W6 W7 W8 W9 W11 W12 W13 W15 W16 W14"
#pfset GeomInput.indi_input.GeomNames "W1 W2 W3 W4 W5 W6 W7 W8 W9 W10 W11 W12 W13 water W15 W16"
pfset Geom.indi_input.FileName pfl_ind.pfb

# sand
pfset GeomInput.W1.Value 1
# loamy sand
pfset GeomInput.W2.Value 2
# sandy loam
pfset GeomInput.W3.Value 3
# silt loam
pfset GeomInput.W4.Value 4
# silt
#pfset GeomInput.W5.Value 5
# loam
pfset GeomInput.W6.Value 6
# sandy clay loam
pfset GeomInput.W7.Value 7
# silty clay loam
pfset GeomInput.W8.Value 8
# clay loam
pfset GeomInput.W9.Value 9
# sandy clay
#pfset GeomInput.W10.Value 10
# silty clay
pfset GeomInput.W11.Value 11
# clay
pfset GeomInput.W12.Value 12
# organic material
pfset GeomInput.W13.Value 13
# water
pfset GeomInput.W14.Value 14
#pfset GeomInput.water.Value 14
# bedrock
pfset GeomInput.W15.Value 15
# others
pfset GeomInput.W16.Value 16

#TIME SETUP
#-----------------------------------------------------------------------------
# Setup timing info
#-----------------------------------------------------------------------------
pfset TimingInfo.BaseUnit                __base_pfl__
pfset TimingInfo.StartCount              __start_cnt_pfl__
pfset TimingInfo.StartTime               0.0
pfset TimingInfo.StopTime                __stop_pfl_bldsva__
pfset TimeStep.Type                      Constant
pfset TimeStep.Value                     __dt_pfl_bldsva__
pfset TimingInfo.DumpInterval            __dump_pfl_interval__

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
#pfset Geom.Perm.Names             "domain"
pfset Geom.Perm.Names		   "domain W1 W2 W3 W4 W6 W7 W8 W9 W11 W12 W13 W15 W16 W14"
#pfset Geom.Perm.Names		   "domain water W1 W2 W3 W4 W5 W6 W7 W8 W9 W10 W11 W12 W13 W15 W16"

pfset Geom.domain.Perm.Type        Constant
pfset Geom.domain.Perm.Value       0.01

#pfset Geom.water.Perm.Type        Constant
#pfset Geom.water.Perm.Value       0.01
#pfset Geom.water.Perm.Value       0.000000001
pfset Geom.W1.Perm.Type           Constant
pfset Geom.W1.Perm.Value          0.269022595
pfset Geom.W2.Perm.Type           Constant
pfset Geom.W2.Perm.Value          0.043630356
pfset Geom.W3.Perm.Type           Constant
pfset Geom.W3.Perm.Value          0.015841225
pfset Geom.W4.Perm.Type           Constant
pfset Geom.W4.Perm.Value          0.007582087
#pfset Geom.W5.Perm.Type           Constant
#pfset Geom.W5.Perm.Value          0.01818816
pfset Geom.W6.Perm.Type           Constant
pfset Geom.W6.Perm.Value          0.026289889
pfset Geom.W7.Perm.Type           Constant
pfset Geom.W7.Perm.Value          0.005492736
pfset Geom.W8.Perm.Type           Constant
pfset Geom.W8.Perm.Value          0.004675077
pfset Geom.W9.Perm.Type           Constant
pfset Geom.W9.Perm.Value          0.003386794
#pfset Geom.W10.Perm.Type          Constant
#pfset Geom.W10.Perm.Value         0.004783973
pfset Geom.W11.Perm.Type          Constant
pfset Geom.W11.Perm.Value         0.003979136
pfset Geom.W12.Perm.Type          Constant
pfset Geom.W12.Perm.Value         0.006162952
pfset Geom.W13.Perm.Type          Constant
pfset Geom.W13.Perm.Value         0.01
#pfset Geom.W13.Perm.Value         0.001
pfset Geom.W15.Perm.Type          Constant
pfset Geom.W15.Perm.Value         0.5
pfset Geom.W16.Perm.Type          Constant
pfset Geom.W16.Perm.Value         0.1
pfset Geom.W14.Perm.Type          Constant
pfset Geom.W14.Perm.Value         0.1

pfset Perm.TensorType			 TensorByGeom
pfset Geom.Perm.TensorByGeom.Names	 "domain"
pfset Geom.domain.Perm.TensorValX	 1000.0
pfset Geom.domain.Perm.TensorValY	 1000.0
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
#pfset Geom.Porosity.GeomNames         "domain"
#pfset Geom.Porosity.GeomNames          "domain W1 W2 W3 W4 W5 W6 W7 W8 W9 W10 W11 W12 W13 W15 W16"
pfset Geom.Porosity.GeomNames          "domain W1 W2 W3 W4 W6 W7 W8 W9 W11 W12 W13 W15 W16 W14"

pfset Geom.domain.Porosity.Type        Constant
pfset Geom.domain.Porosity.Value       0.4

pfset Geom.W1.Porosity.Type            Constant
pfset Geom.W1.Porosity.Value           0.3693
pfset Geom.W2.Porosity.Type            Constant
pfset Geom.W2.Porosity.Value           0.3819
pfset Geom.W3.Porosity.Type            Constant
pfset Geom.W3.Porosity.Value           0.4071
pfset Geom.W4.Porosity.Type            Constant
pfset Geom.W4.Porosity.Value           0.4760
#pfset Geom.W5.Porosity.Type            Constant
#pfset Geom.W5.Porosity.Value           0.4760
pfset Geom.W6.Porosity.Type            Constant
pfset Geom.W6.Porosity.Value           0.4390
pfset Geom.W7.Porosity.Type            Constant
pfset Geom.W7.Porosity.Value           0.4040
pfset Geom.W8.Porosity.Type            Constant
pfset Geom.W8.Porosity.Value           0.4640
pfset Geom.W9.Porosity.Type            Constant
pfset Geom.W9.Porosity.Value           0.4386
#pfset Geom.W10.Porosity.Type           Constant
#pfset Geom.W10.Porosity.Value          0.4134
pfset Geom.W11.Porosity.Type           Constant
pfset Geom.W11.Porosity.Value          0.4789
pfset Geom.W12.Porosity.Type           Constant
pfset Geom.W12.Porosity.Value          0.4680
pfset Geom.W13.Porosity.Type           Constant
pfset Geom.W13.Porosity.Value          0.4
pfset Geom.W15.Porosity.Type           Constant
pfset Geom.W15.Porosity.Value          0.1
pfset Geom.W16.Porosity.Type           Constant
pfset Geom.W16.Porosity.Value          0.1
pfset Geom.W14.Porosity.Type           Constant
pfset Geom.W14.Porosity.Value          0.1

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

pfset Phase.RelPerm.Type               VanGenuchten
#pfset Phase.RelPerm.GeomNames         "domain"
#pfset Phase.RelPerm.GeomNames          "domain W1 W2 W3 W4 W5 W6 W7 W8 W9 W10 W11 W12 W13 W15 W16"
pfset Phase.RelPerm.GeomNames          "domain W1 W2 W3 W4 W6 W7 W8 W9 W11 W12 W13 W15 W16 W14"

pfset Geom.domain.RelPerm.Alpha        2.0
pfset Geom.domain.RelPerm.N            3.

pfset Geom.W1.RelPerm.Alpha            3.548134
pfset Geom.W1.RelPerm.N                3.162278 
pfset Geom.W1.RelPerm.NumSamplePoints  $VG_points
pfset Geom.W1.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W2.RelPerm.Alpha            3.467369
#pfset Geom.W2.RelPerm.N                1.737801
pfset Geom.W2.RelPerm.N                2.01
#pfset Geom.W2.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W2.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W3.RelPerm.Alpha            2.691535
#pfset Geom.W3.RelPerm.N                1.445440
pfset Geom.W3.RelPerm.N                2.01
#pfset Geom.W3.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W3.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W4.RelPerm.Alpha            0.501187
#pfset Geom.W4.RelPerm.N                1.659587
pfset Geom.W4.RelPerm.N                2.01
#pfset Geom.W4.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W4.RelPerm.MinPressureHead  $VG_pmin
#pfset Geom.W5.RelPerm.Alpha            0.660693
##pfset Geom.W5.RelPerm.N                1.659587
#pfset Geom.W5.RelPerm.N                2.01
#pfset Geom.W5.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W5.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W6.RelPerm.Alpha            1.122018
#pfset Geom.W6.RelPerm.N                1.318257
pfset Geom.W6.RelPerm.N                2.01
#pfset Geom.W6.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W6.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W7.RelPerm.Alpha            2.089296
#pfset Geom.W7.RelPerm.N                1.318257
pfset Geom.W7.RelPerm.N                2.01
#pfset Geom.W7.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W7.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W8.RelPerm.Alpha            0.831764
#pfset Geom.W8.RelPerm.N                1.513561
pfset Geom.W8.RelPerm.N                2.01
#pfset Geom.W8.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W8.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W9.RelPerm.Alpha            1.584893
#pfset Geom.W9.RelPerm.N                1.412538
pfset Geom.W9.RelPerm.N                2.01
#pfset Geom.W9.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W9.RelPerm.MinPressureHead  $VG_pmin
#pfset Geom.W10.RelPerm.Alpha           3.311311
##pfset Geom.W10.RelPerm.N               1.202264
#pfset Geom.W10.RelPerm.N               2.01
#pfset Geom.W10.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W10.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W11.RelPerm.Alpha           1.621810
#pfset Geom.W11.RelPerm.N               1.318257
pfset Geom.W11.RelPerm.N               2.01
#pfset Geom.W11.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W11.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W12.RelPerm.Alpha            1.513561
#pfset Geom.W12.RelPerm.N                1.258925
pfset Geom.W12.RelPerm.N                2.01
#pfset Geom.W12.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W12.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W13.RelPerm.Alpha            2.0
pfset Geom.W13.RelPerm.N                3. 
#pfset Geom.W13.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W13.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W15.RelPerm.Alpha            2.0
pfset Geom.W15.RelPerm.N                3.
#pfset Geom.W15.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W15.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W16.RelPerm.Alpha            2.0
pfset Geom.W16.RelPerm.N                3. 
#pfset Geom.W16.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W16.RelPerm.MinPressureHead  $VG_pmin
pfset Geom.W14.RelPerm.Alpha            2.0
pfset Geom.W14.RelPerm.N                3. 
#pfset Geom.W14.RelPerm.NumSamplePoints  $VG_points
#pfset Geom.W14.RelPerm.MinPressureHead  $VG_pmin

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type              VanGenuchten
#pfset Phase.Saturation.GeomNames         "domain"
#pfset Phase.Saturation.GeomNames          "domain W1 W2 W3 W4 W5 W6 W7 W8 W9 W10 W11 W12 W13 W15 W16"
pfset Phase.Saturation.GeomNames          "domain W1 W2 W3 W4 W6 W7 W8 W9 W11 W12 W13 W15 W16 W14"

pfset Geom.domain.Saturation.Alpha        2.0
pfset Geom.domain.Saturation.N            3.
pfset Geom.domain.Saturation.SRes         0.1
pfset Geom.domain.Saturation.SSat         1.0

pfset Geom.W1.Saturation.Alpha            3.548134
pfset Geom.W1.Saturation.N                3.162278
pfset Geom.W1.Saturation.SRes             0.141333
pfset Geom.W1.Saturation.SSat             1.0
pfset Geom.W2.Saturation.Alpha            3.467369
#pfset Geom.W2.Saturation.N                1.737801
pfset Geom.W2.Saturation.SRes             0.125641
pfset Geom.W2.Saturation.SSat             1.0
pfset Geom.W2.Saturation.N                2.01
pfset Geom.W3.Saturation.Alpha            2.691535
#pfset Geom.W3.Saturation.N                1.445440
pfset Geom.W3.Saturation.SRes             0.100775
pfset Geom.W3.Saturation.SSat             1.0
pfset Geom.W3.Saturation.N                2.01
pfset Geom.W4.Saturation.Alpha            0.501187
#pfset Geom.W4.Saturation.N                1.659587
pfset Geom.W4.Saturation.SRes             0.148064
pfset Geom.W4.Saturation.SSat             1.0
pfset Geom.W4.Saturation.N                2.01
#pfset Geom.W5.Saturation.Alpha            0.660693
##pfset Geom.W5.Saturation.N                1.659587
#pfset Geom.W5.Saturation.SRes             0.102249
#pfset Geom.W5.Saturation.SSat             1.0
#pfset Geom.W5.Saturation.N                2.01
pfset Geom.W6.Saturation.Alpha            1.122018
#pfset Geom.W6.Saturation.N                1.479108
pfset Geom.W6.Saturation.SRes             0.152882
pfset Geom.W6.Saturation.SSat             1.0
pfset Geom.W6.Saturation.N                2.01
pfset Geom.W7.Saturation.Alpha            2.089296
#pfset Geom.W7.Saturation.N                1.131826
pfset Geom.W7.Saturation.SRes             0.164063
pfset Geom.W7.Saturation.SSat             1.0
pfset Geom.W7.Saturation.N                2.01
pfset Geom.W8.Saturation.Alpha            0.831764
#pfset Geom.W8.Saturation.N                1.513561
pfset Geom.W8.Saturation.SRes             0.186722
pfset Geom.W8.Saturation.SSat             1.0
pfset Geom.W8.Saturation.N                2.01
pfset Geom.W9.Saturation.Alpha            1.584893
#pfset Geom.W9.Saturation.N                1.412538
pfset Geom.W9.Saturation.SRes             0.178733
pfset Geom.W9.Saturation.SSat             1.0
pfset Geom.W9.Saturation.N                2.01
#pfset Geom.W10.Saturation.Alpha           3.311311
##pfset Geom.W10.Saturation.N               1.202264
#pfset Geom.W10.Saturation.SRes            0.303896
#pfset Geom.W10.Saturation.SSat            1.0
#pfset Geom.W10.Saturation.N               2.01
pfset Geom.W11.Saturation.Alpha           1.621810
#pfset Geom.W11.Saturation.N               1.318257
pfset Geom.W11.Saturation.SRes            0.230769
pfset Geom.W11.Saturation.SSat            1.0
pfset Geom.W11.Saturation.N               2.01
pfset Geom.W12.Saturation.Alpha           1.513561
#pfset Geom.W12.Saturation.N               1.258925
pfset Geom.W12.Saturation.SRes            0.213508
pfset Geom.W12.Saturation.SSat            1.0
pfset Geom.W12.Saturation.N               2.01
pfset Geom.W13.Saturation.Alpha           2.0
pfset Geom.W13.Saturation.N               3.
pfset Geom.W13.Saturation.SRes            0.2
pfset Geom.W13.Saturation.SSat            1.0
pfset Geom.W15.Saturation.Alpha           2.0
pfset Geom.W15.Saturation.N               3.
pfset Geom.W15.Saturation.SRes            0.2
pfset Geom.W15.Saturation.SSat            1.0
pfset Geom.W16.Saturation.Alpha           2.0
pfset Geom.W16.Saturation.N               3.
pfset Geom.W16.Saturation.SRes            0.2
pfset Geom.W16.Saturation.SSat            1.0
pfset Geom.W14.Saturation.Alpha           2.0
pfset Geom.W14.Saturation.N               3.
pfset Geom.W14.Saturation.SRes            0.2
pfset Geom.W14.Saturation.SSat            1.0



#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names				 ""

#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------

pfset BCPressure.PatchNames                   [pfget Geom.domain.Patches]

pfset Patch.top.BCPressure.Type                     OverlandFlow
pfset Patch.top.BCPressure.Cycle                    "constant"
pfset Patch.top.BCPressure.alltime.Value            0.0

pfset Patch.bottom.BCPressure.Type                  FluxConst
pfset Patch.bottom.BCPressure.Cycle                 "constant"
pfset Patch.bottom.BCPressure.alltime.Value         0.0

pfset Patch.perimeter.BCPressure.Type               FluxConst
pfset Patch.perimeter.BCPressure.Cycle              "constant"
pfset Patch.perimeter.BCPressure.alltime.Value      0.0


# Dirichlet BC
 pfset Patch.perimeter.BCPressure.Type               DirEquilRefPatch
 pfset Patch.perimeter.BCPressure.Cycle              "constant"
 pfset Patch.perimeter.BCPressure.RefGeom            domain
 pfset Patch.perimeter.BCPressure.RefPatch           top
 pfset Patch.perimeter.BCPressure.alltime.Value      -0.05


# Dirichlet BC
#pfset Patch.perimeter.BCPressure.Type              PressureFile
#pfset Patch.perimeter.BCPressure.Cycle             "constant"
#set data [pfload -sa bc_hs.sa]
#pfsave $data -pfb bc_hs.pfb
#pfdist bc_hs.pfb
#pfset Patch.perimeter.BCPressure.alltime.FileName  bc_hs.pfb


#  TOPOGRAPHY & SLOPES IN
#  BOTH X- & Y- DIRECTIONS
#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------
pfset TopoSlopesX.Type			 "PFBFile"
pfset TopoSlopesX.GeomNames		 "domain"
pfset TopoSlopesX.FileName		 "slopex.pfb"

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------
pfset TopoSlopesY.Type			 "PFBFile"
pfset TopoSlopesY.GeomNames		 "domain"
pfset TopoSlopesY.FileName		 "slopey.pfb"

#---------------------------------------------------------
# Mannings coefficient
#---------------------------------------------------------
pfset Mannings.Type			 "Constant"
pfset Mannings.GeomNames		 "domain"
pfset Mannings.Geom.domain.Value	 5.5e-5
#pfset Mannings.Geom.domain.Value	 1.e-5

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
pfset Solver.MaxIter				 2500

pfset Solver.TerrainFollowingGrid                True

pfset Solver.Nonlinear.MaxIter			 400
pfset Solver.Nonlinear.ResidualTol		 1e-1
pfset Solver.Nonlinear.EtaChoice		 Walker1
#pfset Solver.Nonlinear.EtaChoice		 EtaConstant
#pfset Solver.Nonlinear.EtaValue			 1.0e-16
pfset Solver.Nonlinear.UseJacobian		 True
pfset Solver.Nonlinear.DerivativeEpsilon	 1e-16
pfset Solver.Nonlinear.StepTol			 1e-16
pfset Solver.Nonlinear.Globalization		 LineSearch
pfset Solver.Linear.KrylovDimension		 30
pfset Solver.Linear.MaxRestart			 8
pfset Solver.MaxConvergenceFailures              8

#pfset Solver.Linear.Preconditioner               PFMGOctree

pfset Solver.Linear.Preconditioner                      PFMG
#pfset Solver.Linear.Preconditioner			 MGSemi
#pfset Solver.Linear.Preconditioner.MGSemi.MaxIter	 1
#pfset Solver.Linear.Preconditioner.MGSemi.MaxLevels	 10
#pfset Solver.PrintSubsurf				 False
#pfset Solver.Drop					 1E-20
#pfset Solver.AbsTol					 1E-12

pfset Solver.PrintSaturation                            True 
pfset Solver.PrintPressure                              True 
#pfset Solver.PrintSubsurf                               False
#pfset Solver.Nonlinear.PrintFlag                        LowVerbosity

pfset Solver.WriteSiloSubsurfData		        False
pfset Solver.WriteSiloPressure				False
pfset Solver.WriteSiloSaturation		        False	
pfset Solver.WriteSiloMask			        False	
pfset Solver.WriteCLMBinary			        False	

#---------------------------------------------------------
# Initial conditions: water pressure
# HydroStaticPatch > PFBFile
#---------------------------------------------------------
#
pfset ICPressure.Type                    __pfl_ICPpressureType__
pfset ICPressure.GeomNames               domain
pfset Geom.domain.ICPressure.Value       __pfl_ICPpressureValue__
pfset Geom.domain.ICPressure.FileName    "__pfl_ICPpressureFileName__"
pfset Geom.domain.ICPressure.RefGeom     domain
pfset Geom.domain.ICPressure.RefPatch    top

#pfset ICPressure.Type                                   HydroStaticPatch
#pfset ICPressure.Type                                   Constant
#pfset ICPressure.GeomNames                              domain
#pfset Geom.domain.ICPressure.Value                      -0.2

#pfset Geom.domain.ICPressure.RefGeom                    domain
#pfset Geom.domain.ICPressure.RefPatch                   top


#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------
pfwritedb cordex0.11
