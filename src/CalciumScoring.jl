module CalciumScoring

using Images
using Statistics
using Distributions
using DSP

include("./agatston.jl")
include("./integrated.jl")
include("./spatially_weighted.jl")
include("./utils.jl")

export score,
    # agatston.jl functions
    Agatston,

    # integrated.jl functions
    Integrated,

    # spatially_weighted.jl functions
    SpatiallyWeighted

    # utils.jl functions
    mask_elements

end
