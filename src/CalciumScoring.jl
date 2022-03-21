module CalciumScoring

using Images
using Statistics
using DICOMUtils
using PhantomSegmentation

include("./agatston.jl")
include("./integrated.jl")
include("./utils.jl")

export score,
    # agatston.jl functions
    Agatston,

    # integrated.jl functions
    Integrated,

    # utils.jl functions
    mask_elements,
    mass_calibration

end
