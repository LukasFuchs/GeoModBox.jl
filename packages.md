GeoModBox - Packages
--------------------
src
---
AdvectionEquation/1Dsolvers.jl:using Dierckx
AdvectionEquation/2Dsolvers.jl:using Interpolations
GeoModBox.jl:    using Statistics: mean
GeoModBox.jl:    using ExtendableSparse, LinearAlgebra
HeatEquation/1Dsolvers.jl:using ExtendableSparse
HeatEquation/AnalyticsDiffusion2D.jl:using ExactFieldSolutions, StaticArrays
HeatEquation/2Dsolvers.jl:using ExtendableSparse
InitialCondition/2Dini.jl:using Base.Threads
InitialCondition/2Dini.jl:using GeoModBox
MomentumEquation/1Dsolvers.jl:using ExtendableSparse
MomentumEquation/2Dsolvers.jl:using ExtendableSparse
Tracers/2Dsolvers.jl:using Base.Threads
Tracers/2Dsolvers.jl:using Statistics, Printf

examples
--------
AdvectionEquation/2D_Advection.jl:# Two-dimensional advection solver using different numerical methods.
AdvectionEquation/2D_Advection.jl:# using Statistics
AdvectionEquation/2D_Advection.jl:using Plots, Interpolations
AdvectionEquation/2D_Advection.jl:using GeoModBox.AdvectionEquation.TwoD, GeoModBox.Tracers.TwoD
AdvectionEquation/2D_Advection.jl:using GeoModBox.InitialCondition
AdvectionEquation/2D_Advection.jl:using Base.Threads
AdvectionEquation/2D_Advection.jl:using Printf, TimerOutputs, LaTeXStrings, Measures
AdvectionEquation/2D_Advection_ResolutionTest.jl:using Plots, Interpolations, Statistics
AdvectionEquation/2D_Advection_ResolutionTest.jl:using GeoModBox.AdvectionEquation.TwoD
AdvectionEquation/2D_Advection_ResolutionTest.jl:using GeoModBox.InitialCondition, GeoModBox.Tracers.TwoD
AdvectionEquation/2D_Advection_ResolutionTest.jl:using Base.Threads
AdvectionEquation/2D_Advection_ResolutionTest.jl:using Printf, TimerOutputs, LaTeXStrings, Measures
DiffusionEquation/2D/ForwardEuler.jl:using GeoModBox.HeatEquation.TwoD
DiffusionEquation/2D/ForwardEuler.jl:using Plots
DiffusionEquation/2D/ForwardEuler.jl:using TimerOutputs
DiffusionEquation/2D/ForwardEuler.jl:using ExactFieldSolutions
DiffusionEquation/2D/BackwardEuler.jl:using GeoModBox.HeatEquation.TwoD
DiffusionEquation/2D/BackwardEuler.jl:using ExactFieldSolutions, LinearAlgebra, Plots, Printf
DiffusionEquation/2D/BackwardEuler.jl:using TimerOutputs
DiffusionEquation/2D/GeneralSolverTest.jl:using Plots, GeoModBox.HeatEquation.TwoD, ExtendableSparse
DiffusionEquation/2D/GeneralSolverTest.jl:using Statistics, LinearAlgebra
DiffusionEquation/2D/GeneralSolverTest.jl:using TimerOutputs, LaTeXStrings, Measures
DiffusionEquation/2D/GeneralSolverTest.jl:using ExactFieldSolutions
DiffusionEquation/2D/GeneralSolverTest_variable_k.jl:using Plots, GeoModBox.HeatEquation.TwoD, ExtendableSparse
DiffusionEquation/2D/GeneralSolverTest_variable_k.jl:using Statistics, LinearAlgebra
DiffusionEquation/2D/GeneralSolverTest_variable_k.jl:using TimerOutputs, ExactFieldSolutions
DiffusionEquation/2D/Poisson_ResTest.jl:using GeoModBox.HeatEquation.TwoD, ExtendableSparse, Plots, Statistics
DiffusionEquation/2D/Poisson_ResTest.jl:using GLM, DataFrames
DiffusionEquation/2D/Poisson_ResTest.jl:using TimerOutputs
DiffusionEquation/2D/Poisson_variable_k.jl:using GeoModBox.HeatEquation.TwoD, ExtendableSparse, Plots
DiffusionEquation/2D/Poisson_variable_k.jl:using TimerOutputs
DiffusionEquation/2D/Gaussian_Diffusion.jl:using Plots, GeoModBox.HeatEquation.TwoD, ExtendableSparse
DiffusionEquation/2D/Gaussian_Diffusion.jl:using Statistics
DiffusionEquation/2D/Gaussian_Diffusion.jl:using TimerOutputs, LaTeXStrings, Measures
DiffusionEquation/2D/Gaussian_Diffusion.jl:using ExactFieldSolutions
DiffusionEquation/1D/Heat_1D_dc.jl:using Plots, Printf, LinearAlgebra, ExtendableSparse
DiffusionEquation/1D/Heat_1D_dc.jl:using GeoModBox.HeatEquation.OneD
DiffusionEquation/1D/Heat_1D_dc.jl:using TimerOutputs
DiffusionEquation/1D/OceanicGeotherm_1D.jl:using Plots, SpecialFunctions
DiffusionEquation/1D/OceanicGeotherm_1D.jl:using GeoModBox.HeatEquation.OneD
DiffusionEquation/1D/OceanicGeotherm_1D.jl:using TimerOutputs
DiffusionEquation/1D/ContinentalGeotherm_1D_dc.jl:using Plots, ExtendableSparse, LinearAlgebra
DiffusionEquation/1D/ContinentalGeotherm_1D_dc.jl:using GeoModBox.HeatEquation.OneD
DiffusionEquation/1D/ContinentalGeotherm_1D_dc.jl:using TimerOutputs
DiffusionEquation/1D/Heat_1D_discretization.jl:using Plots, Printf, LinearAlgebra, ExtendableSparse
DiffusionEquation/1D/Heat_1D_discretization.jl:using GeoModBox.HeatEquation.OneD
DiffusionEquation/1D/Heat_1D_discretization.jl:using TimerOutputs
DiffusionEquation/1D/ContinentalGeotherm_1D.jl:using Plots
DiffusionEquation/1D/ContinentalGeotherm_1D.jl:using GeoModBox.HeatEquation.OneD
DiffusionEquation/1D/ContinentalGeotherm_1D.jl:using TimerOutputs
DiffusionEquation/1D/OceanicGeotherm_1D_dc.jl:using Plots, SpecialFunctions
DiffusionEquation/1D/OceanicGeotherm_1D_dc.jl:using GeoModBox.HeatEquation.OneD
DiffusionEquation/1D/OceanicGeotherm_1D_dc.jl:using TimerOutputs, ExtendableSparse, Printf, LinearAlgebra
ShearHeating/2D/ShearHeatingShearBands.jl:using Plots, ExtendableSparse, Base.Threads
ShearHeating/2D/ShearHeatingShearBands.jl:using GeoModBox, GeoModBox.Tracers.TwoD
ShearHeating/2D/ShearHeatingShearBands.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
ShearHeating/2D/ShearHeatingShearBands.jl:using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
ShearHeating/2D/ShearHeatingShearBands.jl:using Statistics
ShearHeating/2D/ShearHeatingShearBands.jl:using Printf, LinearAlgebra, TimerOutputs, Interpolations, LsqFit
ShearHeating/2D/ShearHeatingShearBands.jl:using Measures
ShearHeating/2D/PlotData.jl:using DelimitedFiles
ShearHeating/2D/PlotData.jl:using Plots
ShearHeating/2D/PlotData.jl:using LaTeXStrings
ShearHeating/2D/PlotData.jl:using FileIO, ImageIO, Measures
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNM.jl:using Plots
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNM.jl:using ExtendableSparse
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNM.jl:using GeoModBox
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNM.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNM.jl:using GeoModBox.AdvectionEquation.TwoD
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNM.jl:using GeoModBox.Tracers.TwoD
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNM.jl:using Base.Threads
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNM.jl:using Printf, LinearAlgebra
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNM.jl:using TimerOutputs
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNM.jl:using Random
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNM.jl:    # properties using a mixing law, whereas the direct interpolation 
StokesEquation/2D/ViscousInclusion.jl:using GeoModBox, Printf, Plots
StokesEquation/2D/ViscousInclusion.jl:using GeoModBox.InitialCondition, GeoModBox.Tracers.TwoD
StokesEquation/2D/ViscousInclusion.jl:using GeoModBox.MomentumEquation.TwoD
StokesEquation/2D/ViscousInclusion.jl:using Base.Threads, LinearAlgebra, Statistics
StokesEquation/2D/RTI.jl:using Plots
StokesEquation/2D/RTI.jl:using ExtendableSparse
StokesEquation/2D/RTI.jl:using GeoModBox
StokesEquation/2D/RTI.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
StokesEquation/2D/RTI.jl:using GeoModBox.AdvectionEquation.TwoD
StokesEquation/2D/RTI.jl:using GeoModBox.Tracers.TwoD
StokesEquation/2D/RTI.jl:using Base.Threads
StokesEquation/2D/RTI.jl:using Printf, LinearAlgebra
StokesEquation/2D/RTI.jl:using TimerOutputs
StokesEquation/2D/RTI.jl:    # properties using a mixing law, whereas the direct interpolation 
StokesEquation/2D/VanKekenBenchmark_scaled.jl:using Plots
StokesEquation/2D/VanKekenBenchmark_scaled.jl:using ExtendableSparse
StokesEquation/2D/VanKekenBenchmark_scaled.jl:using GeoModBox
StokesEquation/2D/VanKekenBenchmark_scaled.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
StokesEquation/2D/VanKekenBenchmark_scaled.jl:using GeoModBox.AdvectionEquation.TwoD, GeoModBox.Scaling
StokesEquation/2D/VanKekenBenchmark_scaled.jl:using GeoModBox.Tracers.TwoD
StokesEquation/2D/VanKekenBenchmark_scaled.jl:using Base.Threads
StokesEquation/2D/VanKekenBenchmark_scaled.jl:using Printf, LinearAlgebra
StokesEquation/2D/VanKekenBenchmark_scaled.jl:using TimerOutputs
StokesEquation/2D/VanKekenBenchmark_scaled.jl:using Statistics
StokesEquation/2D/FallingBlockBenchmark_DC.jl:using Plots, ExtendableSparse, LinearAlgebra
StokesEquation/2D/FallingBlockBenchmark_DC.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
StokesEquation/2D/FallingBlockBenchmark_DC.jl:using GeoModBox.AdvectionEquation.TwoD
StokesEquation/2D/FallingBlockBenchmark_DC.jl:using GeoModBox.Tracers.TwoD
StokesEquation/2D/FallingBlockBenchmark_DC.jl:using Base.Threads
StokesEquation/2D/FallingBlockBenchmark_DC.jl:using Printf, TimerOutputs, LaTeXStrings, Measures
StokesEquation/2D/VanKekenBenchmark.jl:using Plots
StokesEquation/2D/VanKekenBenchmark.jl:using ExtendableSparse
StokesEquation/2D/VanKekenBenchmark.jl:using GeoModBox
StokesEquation/2D/VanKekenBenchmark.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
StokesEquation/2D/VanKekenBenchmark.jl:using GeoModBox.AdvectionEquation.TwoD
StokesEquation/2D/VanKekenBenchmark.jl:using GeoModBox.Tracers.TwoD
StokesEquation/2D/VanKekenBenchmark.jl:using Base.Threads
StokesEquation/2D/VanKekenBenchmark.jl:using Printf, LinearAlgebra
StokesEquation/2D/VanKekenBenchmark.jl:using TimerOutputs
StokesEquation/2D/VanKekenBenchmark.jl:using Statistics
StokesEquation/2D/RTI_GrowthRate.jl:using Plots
StokesEquation/2D/RTI_GrowthRate.jl:using ExtendableSparse
StokesEquation/2D/RTI_GrowthRate.jl:using GeoModBox
StokesEquation/2D/RTI_GrowthRate.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
StokesEquation/2D/RTI_GrowthRate.jl:using GeoModBox.AdvectionEquation.TwoD
StokesEquation/2D/RTI_GrowthRate.jl:using GeoModBox.Tracers.TwoD
StokesEquation/2D/RTI_GrowthRate.jl:using Base.Threads
StokesEquation/2D/RTI_GrowthRate.jl:using Printf, LinearAlgebra
StokesEquation/2D/RTI_GrowthRate.jl:using TimerOutputs
StokesEquation/2D/RTI_GrowthRate.jl:    # properties using a mixing law, whereas the direct interpolation 
StokesEquation/2D/FallingBlockConstEta_DC.jl:using Plots, ExtendableSparse
StokesEquation/2D/FallingBlockConstEta_DC.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
StokesEquation/2D/FallingBlockConstEta_DC.jl:using Printf, LinearAlgebra
StokesEquation/2D/FallingBlockConstEta_DC.jl:using TimerOutputs
StokesEquation/2D/FallingBlockConstEta_DC.jl:    # problem using the defect correction solution method.                #
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNC.jl:using Plots
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNC.jl:using ExtendableSparse
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNC.jl:using GeoModBox
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNC.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNC.jl:using GeoModBox.AdvectionEquation.TwoD
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNC.jl:using GeoModBox.Tracers.TwoD
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNC.jl:using Base.Threads
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNC.jl:using Printf, LinearAlgebra
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNC.jl:using TimerOutputs
StokesEquation/2D/RTI_Growth_Rate_Res_Test_CNC.jl:    # properties using a mixing law, whereas the direct interpolation 
StokesEquation/2D/FallingBlockBenchmark_direct.jl:using Plots, ExtendableSparse
StokesEquation/2D/FallingBlockBenchmark_direct.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
StokesEquation/2D/FallingBlockBenchmark_direct.jl:using GeoModBox.AdvectionEquation.TwoD
StokesEquation/2D/FallingBlockBenchmark_direct.jl:using GeoModBox.Tracers.TwoD
StokesEquation/2D/FallingBlockBenchmark_direct.jl:using Base.Threads
StokesEquation/2D/FallingBlockBenchmark_direct.jl:using Printf, TimerOutputs, LaTeXStrings, Measures
StokesEquation/2D/FallingBlockVarEta_DC.jl:using Plots
StokesEquation/2D/FallingBlockVarEta_DC.jl:using ExtendableSparse
StokesEquation/2D/FallingBlockVarEta_DC.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
StokesEquation/2D/FallingBlockVarEta_DC.jl:using GeoModBox.AdvectionEquation.TwoD
StokesEquation/2D/FallingBlockVarEta_DC.jl:using GeoModBox.Tracers.TwoD
StokesEquation/2D/FallingBlockVarEta_DC.jl:using Base.Threads
StokesEquation/2D/FallingBlockVarEta_DC.jl:using Printf, LinearAlgebra
StokesEquation/2D/FallingBlockVarEta_DC.jl:using TimerOutputs
StokesEquation/1D/ChannelFlow_1D.jl:using Plots, ExtendableSparse, Printf, LinearAlgebra
StokesEquation/1D/ChannelFlow_1D.jl:using GeoModBox.MomentumEquation.OneD
StokesEquation/1D/ChannelFlow_1D.jl:using TimerOutputs
ThermalConvection/BottomHeated.jl:using Plots, ExtendableSparse
ThermalConvection/BottomHeated.jl:using GeoModBox
ThermalConvection/BottomHeated.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
ThermalConvection/BottomHeated.jl:using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
ThermalConvection/BottomHeated.jl:using GeoModBox.Scaling
ThermalConvection/BottomHeated.jl:using Statistics, LinearAlgebra
ThermalConvection/BottomHeated.jl:using Printf, TimerOutputs, LaTeXStrings, Measures
ThermalConvection/README.md:- The momentum conservation equation is solved using the general solver.
ThermalConvection/README.md:- The energy conservation equation is solved using the general solver with a Crank–Nicolson discretization.
ThermalConvection/README.md:- Temperature advection is performed using a semi-Lagrangian scheme.
ThermalConvection/MixedHeated.jl:using Plots, ExtendableSparse
ThermalConvection/MixedHeated.jl:using GeoModBox
ThermalConvection/MixedHeated.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
ThermalConvection/MixedHeated.jl:using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
ThermalConvection/MixedHeated.jl:using GeoModBox.Scaling
ThermalConvection/MixedHeated.jl:using Statistics, LinearAlgebra
ThermalConvection/MixedHeated.jl:using Printf, TimerOutputs, LaTeXStrings, Measures
ThermalConvection/BottomHeated_VP.jl:using Plots, ExtendableSparse
ThermalConvection/BottomHeated_VP.jl:using GeoModBox
ThermalConvection/BottomHeated_VP.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
ThermalConvection/BottomHeated_VP.jl:using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
ThermalConvection/BottomHeated_VP.jl:using GeoModBox.Scaling
ThermalConvection/BottomHeated_VP.jl:# using GeoModBox.Tracers.TwoD
ThermalConvection/BottomHeated_VP.jl:# using Base.Threads
ThermalConvection/BottomHeated_VP.jl:using Statistics, LinearAlgebra
ThermalConvection/BottomHeated_VP.jl:using Printf, TimerOutputs
ThermalConvection/InternallyHeated.jl:using Plots, ExtendableSparse
ThermalConvection/InternallyHeated.jl:using GeoModBox
ThermalConvection/InternallyHeated.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
ThermalConvection/InternallyHeated.jl:using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
ThermalConvection/InternallyHeated.jl:using GeoModBox.Scaling
ThermalConvection/InternallyHeated.jl:using Statistics, LinearAlgebra
ThermalConvection/InternallyHeated.jl:using Printf, TimerOutputs, LaTeXStrings, Measures
ThermalConvection/BottomHeated_VarEta.jl:using Plots, ExtendableSparse
ThermalConvection/BottomHeated_VarEta.jl:using GeoModBox
ThermalConvection/BottomHeated_VarEta.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
ThermalConvection/BottomHeated_VarEta.jl:using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
ThermalConvection/BottomHeated_VarEta.jl:using GeoModBox.Scaling
ThermalConvection/BottomHeated_VarEta.jl:using Statistics, LinearAlgebra
ThermalConvection/BottomHeated_VarEta.jl:using Printf, TimerOutputs, LaTeXStrings, Measures
ThermalConvection/BottomHeated_VarEta.jl:# The momentum equation is solved using η/η₀,
ThermalConvection/Q_H_Relation.jl:using Printf
ThermalConvection/Blankenbach_var_eta.jl:using Plots, ExtendableSparse
ThermalConvection/Blankenbach_var_eta.jl:using GeoModBox
ThermalConvection/Blankenbach_var_eta.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
ThermalConvection/Blankenbach_var_eta.jl:using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
ThermalConvection/Blankenbach_var_eta.jl:using GeoModBox.Scaling
ThermalConvection/Blankenbach_var_eta.jl:using Statistics, LinearAlgebra
ThermalConvection/Blankenbach_var_eta.jl:using Printf, LaTeXStrings, Measures


exercises
---------
01_1D_Euler_Advection.ipynb:    "using Plots"
01_1D_Euler_Advection_en.ipynb:    "using Plots\n",
02_1D_Heat_explicit.ipynb:    "using Plots\n",
02_1D_Heat_explicit.ipynb:    "using GeoModBox.HeatEquation.OneD\n",
02_1D_Heat_explicit_en.ipynb:    "using Plots\n",
02_1D_Heat_explicit_en.ipynb:    "using GeoModBox.HeatEquation.OneD\n",
03_1D_Heat_implicit.ipynb:    "using Plots, ExtendableSparse, LinearAlgebra, Printf\n",
03_1D_Heat_implicit.ipynb:    "using GeoModBox.HeatEquation.OneD\n",
03_1D_Heat_implicit_en.ipynb:    "using Plots, ExtendableSparse, LinearAlgebra, Printf\n",
03_1D_Heat_implicit_en.ipynb:    "using GeoModBox.HeatEquation.OneD\n",
04_2D_Diffusion_Stationary.ipynb:    "using GeoModBox.HeatEquation.TwoD, ExtendableSparse, Plots\n",
04_2D_Diffusion_Stationary_en.ipynb:    "using GeoModBox.HeatEquation.TwoD, ExtendableSparse, Plots\n",
05_2D_Diffusion_TD_Plume.ipynb:    "using Plots, GeoModBox.HeatEquation.TwoD, ExtendableSparse\n",
05_2D_Diffusion_TD_Plume.ipynb:    "using Printf, LinearAlgebra\n",
05_2D_Diffusion_TD_Plume_en.ipynb:    "using Plots, GeoModBox.HeatEquation.TwoD, ExtendableSparse\n",
05_2D_Diffusion_TD_Plume_en.ipynb:    "using Printf, LinearAlgebra\n",
05_2D_Diffusion_TD_Sill.ipynb:    "using Plots, GeoModBox.HeatEquation.TwoD, ExtendableSparse\n",
05_2D_Diffusion_TD_Sill.ipynb:    "using Printf, LinearAlgebra\n",
05_2D_Diffusion_TD_Sill_en.ipynb:    "using Plots, GeoModBox.HeatEquation.TwoD, ExtendableSparse\n",
05_2D_Diffusion_TD_Sill_en.ipynb:    "using Printf, LinearAlgebra\n",
06_1D_Advection.ipynb:    "using Plots, Interpolations\n",
06_1D_Advection.ipynb:    "using GeoModBox.AdvectionEquation.OneD, GeoModBox.Tracers.OneD\n",
06_1D_Advection_en.ipynb:    "using Plots, Interpolations\n",
06_1D_Advection_en.ipynb:    "using GeoModBox.AdvectionEquation.OneD, GeoModBox.Tracers.OneD\n",
07_2D_Energy_Equation.ipynb:    "using Plots, Interpolations, ExtendableSparse, LinearAlgebra\n",
07_2D_Energy_Equation.ipynb:    "using GeoModBox # Enthält die Strukturen\n",
07_2D_Energy_Equation.ipynb:    "using GeoModBox.AdvectionEquation.TwoD, GeoModBox.InitialCondition\n",
07_2D_Energy_Equation.ipynb:    "using GeoModBox.HeatEquation.TwoD, GeoModBox.Tracers.TwoD \n",
07_2D_Energy_Equation.ipynb:    "using Base.Threads, Printf\n",
07_2D_Energy_Equation_en.ipynb:    "using Plots, Interpolations, ExtendableSparse, LinearAlgebra\n",
07_2D_Energy_Equation_en.ipynb:    "using GeoModBox\n",
07_2D_Energy_Equation_en.ipynb:    "using GeoModBox.AdvectionEquation.TwoD, GeoModBox.InitialCondition\n",
07_2D_Energy_Equation_en.ipynb:    "using GeoModBox.HeatEquation.TwoD, GeoModBox.Tracers.TwoD \n",
07_2D_Energy_Equation_en.ipynb:    "using Base.Threads, Printf\n",
08_1D_Stokes.ipynb:    "using Plots, ExtendableSparse\n",
08_1D_Stokes_en.ipynb:    "using Plots, ExtendableSparse\n",
09_2D_Falling_Block.ipynb:    "using Plots\n",
09_2D_Falling_Block.ipynb:    "using ExtendableSparse\n",
09_2D_Falling_Block.ipynb:    "using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD\n",
09_2D_Falling_Block_en.ipynb:    "using Plots\n",
09_2D_Falling_Block_en.ipynb:    "using ExtendableSparse\n",
09_2D_Falling_Block_en.ipynb:    "using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD"
10_2D_Falling_Block_td.ipynb:    "using Plots\n",
10_2D_Falling_Block_td.ipynb:    "using ExtendableSparse\n",
10_2D_Falling_Block_td.ipynb:    "using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD\n",
10_2D_Falling_Block_td.ipynb:    "using GeoModBox.AdvectionEquation.TwoD\n",
10_2D_Falling_Block_td.ipynb:    "using GeoModBox.Tracers.TwoD\n",
10_2D_Falling_Block_td.ipynb:    "using Base.Threads\n",
10_2D_Falling_Block_td.ipynb:    "using Printf\n",
10_2D_Falling_Block_td_en.ipynb:    "using Plots\n",
10_2D_Falling_Block_td_en.ipynb:    "using ExtendableSparse\n",
10_2D_Falling_Block_td_en.ipynb:    "using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD\n",
10_2D_Falling_Block_td_en.ipynb:    "using GeoModBox.AdvectionEquation.TwoD\n",
10_2D_Falling_Block_td_en.ipynb:    "using GeoModBox.Tracers.TwoD\n",
10_2D_Falling_Block_td_en.ipynb:    "using Base.Threads\n",
10_2D_Falling_Block_td_en.ipynb:    "using Printf\n",
11_2D_Thermal_Convection.ipynb:    "using Plots, ExtendableSparse, Base.Threads\n",
11_2D_Thermal_Convection.ipynb:    "using GeoModBox\n",
11_2D_Thermal_Convection.ipynb:    "using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD\n",
11_2D_Thermal_Convection.ipynb:    "using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD\n",
11_2D_Thermal_Convection.ipynb:    "using GeoModBox.Tracers.TwoD\n",
11_2D_Thermal_Convection.ipynb:    "using Statistics, LinearAlgebra\n",
11_2D_Thermal_Convection.ipynb:    "using Printf\n",
11_2D_Thermal_Convection_en.ipynb:    "using Plots, ExtendableSparse, Base.Threads\n",
11_2D_Thermal_Convection_en.ipynb:    "using GeoModBox\n",
11_2D_Thermal_Convection_en.ipynb:    "using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD\n",
11_2D_Thermal_Convection_en.ipynb:    "using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD\n",
11_2D_Thermal_Convection_en.ipynb:    "using GeoModBox.Tracers.TwoD\n",
11_2D_Thermal_Convection_en.ipynb:    "using Statistics, LinearAlgebra\n",
11_2D_Thermal_Convection_en.ipynb:    "using Printf\n",
12_2D_Thermal_Convection_scaled.ipynb:    "using Plots, ExtendableSparse\n",
12_2D_Thermal_Convection_scaled.ipynb:    "using GeoModBox\n",
12_2D_Thermal_Convection_scaled.ipynb:    "using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD\n",
12_2D_Thermal_Convection_scaled.ipynb:    "using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD\n",
12_2D_Thermal_Convection_scaled.ipynb:    "using GeoModBox.Scaling\n",
12_2D_Thermal_Convection_scaled.ipynb:    "using GeoModBox.Tracers.TwoD\n",
12_2D_Thermal_Convection_scaled.ipynb:    "using Base.Threads\n",
12_2D_Thermal_Convection_scaled.ipynb:    "using Statistics, LinearAlgebra\n",
12_2D_Thermal_Convection_scaled.ipynb:    "using Printf, LaTeXStrings, Measures\n",
12_2D_Thermal_Convection_scaled_en.ipynb:    "using Plots, ExtendableSparse\n",
12_2D_Thermal_Convection_scaled_en.ipynb:    "using GeoModBox\n",
12_2D_Thermal_Convection_scaled_en.ipynb:    "using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD\n",
12_2D_Thermal_Convection_scaled_en.ipynb:    "using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD\n",
12_2D_Thermal_Convection_scaled_en.ipynb:    "using GeoModBox.Scaling\n",
12_2D_Thermal_Convection_scaled_en.ipynb:    "using GeoModBox.Tracers.TwoD\n",
12_2D_Thermal_Convection_scaled_en.ipynb:    "using Base.Threads\n",
12_2D_Thermal_Convection_scaled_en.ipynb:    "using Statistics, LinearAlgebra\n",
12_2D_Thermal_Convection_scaled_en.ipynb:    "using Printf, LaTeXStrings, Measures\n",
13_Blankenbach_Benchmark.ipynb:    "using Plots, ExtendableSparse\n",
13_Blankenbach_Benchmark.ipynb:    "using GeoModBox\n",
13_Blankenbach_Benchmark.ipynb:    "using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD\n",
13_Blankenbach_Benchmark.ipynb:    "using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD\n",
13_Blankenbach_Benchmark.ipynb:    "using GeoModBox.Scaling\n",
13_Blankenbach_Benchmark.ipynb:    "using Statistics, LinearAlgebra\n",
13_Blankenbach_Benchmark.ipynb:    "using Printf, LaTeXStrings, Measures\n",
13_Blankenbach_Benchmark.jl:using Plots, ExtendableSparse
13_Blankenbach_Benchmark.jl:using GeoModBox
13_Blankenbach_Benchmark.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
13_Blankenbach_Benchmark.jl:using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
13_Blankenbach_Benchmark.jl:using GeoModBox.Scaling
13_Blankenbach_Benchmark.jl:using Statistics, LinearAlgebra
13_Blankenbach_Benchmark.jl:using Printf, LaTeXStrings, Measures
13_Blankenbach_Benchmark_en.ipynb:    "using Plots, ExtendableSparse\n",
13_Blankenbach_Benchmark_en.ipynb:    "using GeoModBox\n",
13_Blankenbach_Benchmark_en.ipynb:    "using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD\n",
13_Blankenbach_Benchmark_en.ipynb:    "using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD\n",
13_Blankenbach_Benchmark_en.ipynb:    "using GeoModBox.Scaling\n",
13_Blankenbach_Benchmark_en.ipynb:    "using Statistics, LinearAlgebra\n",
13_Blankenbach_Benchmark_en.ipynb:    "using Printf, LaTeXStrings, Measures\n",
Blankenbach_ResTest.jl:using Plots, ExtendableSparse
Blankenbach_ResTest.jl:using GeoModBox
Blankenbach_ResTest.jl:using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
Blankenbach_ResTest.jl:using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
Blankenbach_ResTest.jl:using GeoModBox.Scaling
Blankenbach_ResTest.jl:using Statistics, LinearAlgebra
Blankenbach_ResTest.jl:using Printf, LaTeXStrings, Measures











