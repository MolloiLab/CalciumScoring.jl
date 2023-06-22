module CalciumScoring

using Statistics
using Distributions
using DSP

include("./agatston.jl")
include("./integrated.jl")
include("./spatially_weighted.jl")
include("./volume_fraction.jl")

export score,
    # agatston.jl functions
    CalciumScore, Agatston,

    # integrated.jl functions
    Integrated,

    # material_decomposition functions,
    MaterialDecomposition,
    fit_calibration,

    # spatially_weighted.jl functions
    SpatiallyWeighted,

    # # utils.jl functions
    # mask_elements,

    # volume_fraction.jl
    VolumeFraction

end
