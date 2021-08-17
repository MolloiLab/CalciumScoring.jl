module CalciumScoring

using Images

include("./Integrated.jl")
include("./utils.jl")

export
    score,

    # Integrated.jl functions
    Integrated,
    
    # utils.jl functions
    mask_elements

end
