using Plots
using ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD
using GeoModBox.Tracers.TwoD
using Base.Threads
using Printf

function RTI()
    # Define Initial Condition ========================================== #
    #   1) block
    Ini         =   (p=:RTI,) 
    λ           =   1.0         #   Perturbation wavelength[ km ]
    δA          =   500/15      #   Amplitude [ km ]
    # ------------------------------------------------------------------- #
    # Plot Settings ===================================================== #
    Pl  =   (
        qinc    =   4,
        mainc   =   2,
        qsc     =   100*(60*60*24*365.25)*5e1
    )
    # ------------------------------------------------------------------- #
    # Geometry ========================================================== #
    M       =   Geometry(
        ymin    =   -1.0e3,     # [ m ]
        ymax    =   0.0,
        xmin    =   0.0,
    )
    # -------------------------------------------------------------------- 
end

RTI()