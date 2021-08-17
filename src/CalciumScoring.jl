module CalciumScoring

using Images

include("./Integrated.jl")
include("./utils.jl")

export
    # Integrated.jl functions
    score,
    Integrated,
    
    # utils.jl functions
    mask_elements

end
