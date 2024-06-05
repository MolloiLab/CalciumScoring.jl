using Preferences: set_preferences!
set_preferences!("CalciumScoring", "instability_check" => "error")

include("agatston.jl")
include("integrated.jl")
include("spatially_weighted.jl")
include("volume_fraction.jl")
