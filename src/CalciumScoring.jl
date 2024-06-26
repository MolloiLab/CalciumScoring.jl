module CalciumScoring
using DispatchDoctor: @stable
@stable default_mode="disable" begin

include("agatston.jl")
include("integrated.jl")
include("spatially_weighted.jl")
include("volume_fraction.jl")

export score,
    # agatston.jl functions
    CalciumScore, Agatston,

    # integrated.jl functions
    Integrated,

    # spatially_weighted.jl functions
    SpatiallyWeighted,

    # volume_fraction.jl
    VolumeFraction

end #stable

end  #module