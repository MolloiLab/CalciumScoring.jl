using CalciumScoring
using DispatchDoctor: @stable

@stable default_mode="disable" begin

"""
## `Integrated`

```julia
function Integrated(
    vol, 
    I=sum(vol), 
    N=length(vol)
)
```

#### Inputs
- vol: region of interest
- I: integrated intensity of the entire `vol`
- N: total number of voxels contained in the entire `vol`

#### References
[Integrated intensity-based technique for coronary artery calcium mass measurement: A phantom study](https://doi.org/10.1002/mp.16326)

[Accurate quantification of vessel cross-sectional area using CT angiography: a simulation study](https://doi.org/10.1007/s10554-016-1007-9)
"""
struct Integrated{T1<:AbstractArray,T2<:AbstractFloat} <: CalciumScoring.CalciumScore
    vol::T1
    I::T2
    N::T2
end

function Integrated(
    vol::AbstractArray,
    I::AbstractFloat=Float64.(sum(vol)),
    N::AbstractFloat=Float64.(length(vol)),
)
    return Integrated(I, N)
end

"""
## `score` (Integrated)

```julia
score(S_Bkg, S_Obj, algorithm::Integrated)

score(S_Bkg, S_Obj, size, algorithm::Integrated)

score(S_Bkg, S_Obj, size, ρ, algorithm::Integrated)
```

Use the integrated intensity approach to find the total number of object
voxels/volume/mass contained in `vol`.

The function can be called with different parameters:

1. `S_Bkg`, `S_Obj`, `algorithm` - Returns the total number of object voxels.
2. `S_Bkg`, `S_Obj`, `size`, `algorithm` - Returns the total object volume.
3. `S_Bkg`, `S_Obj`, `size`, `ρ`, `algorithm` - Returns the total object mass.

#### Inputs
- `S_Bkg`: pure background signal intensity in the `vol`
- `S_Obj`: pure object signal intensity in the `vol`
- `size`: (optional) size of the voxels in the given `vol` (only for the 2nd and 3rd variations)
- `ρ`: (optional) density of the given object of interest contained within the `vol` (only for the 3rd variation)
- `algorithm`: the `Integrated` algorithm which accounts for noise and partial volume effect by integrating the intensity of the entire volume `vol`

#### References
[Integrated intensity-based technique for coronary artery calcium mass measurement: A phantom study](https://doi.org/10.1002/mp.16326)

[Accurate quantification of vessel cross-sectional area using CT angiography: a simulation study](https://doi.org/10.1007/s10554-016-1007-9)
"""
function score(S_Bkg, S_Obj, algorithm::Integrated)
    I = algorithm.I
    N = algorithm.N

    N_Obj = (I - (N * S_Bkg)) / (S_Obj - S_Bkg)
    return N_Obj
end

function score(S_Bkg, S_Obj, size, algorithm::Integrated)
    I = algorithm.I
    N = algorithm.N

    N_Obj = (I - (N * S_Bkg)) / (S_Obj - S_Bkg)
    Vol_Obj = N_Obj * size[1] * size[2] * size[3]
    return Vol_Obj
end

function score(S_Bkg, S_Obj, size, ρ, algorithm::Integrated)
    I = algorithm.I
    N = algorithm.N

    N_Obj = (I - (N * S_Bkg)) / (S_Obj - S_Bkg)
    Vol_Obj = N_Obj * size[1] * size[2] * size[3]
    Mass_Obj = Vol_Obj * ρ
    return Mass_Obj
end

end # stable