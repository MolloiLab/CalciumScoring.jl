module CalciumScoring

using Images

include("./agatston.jl")
include("./integrated.jl")
include("./utils.jl")

export score,
    # agatston.jl functions
    Agatston,

    # integrated.jl functions
    Integrated,

    # utils.jl functions
    mask_elements

end
