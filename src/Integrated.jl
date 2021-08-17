abstract type CalciumScore end

"""
    function Integrated(
        vol, 
        I=sum(vol), 
        N=length(vol)
        )

# Arguments
- vol: region of interest
- I: integrated intensity of the entire `vol`
- N: total number of voxels contained in the entire `vol`

# Citation
'Accurate quantification of vessel cross-sectional area using CT
angiography: a simulation study' [Molloi, Johnson, Ding, Lipinski] 
(DOI: 10.1007/s10554-016-1007-9)
"""
struct Integrated{T1<:AbstractArray,T2<:AbstractFloat} <: CalciumScore
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
    score(vol, S_Bkg, S_Obj, algorithm::Integrated)
    score(vol, S_Bkg, S_Obj, size, algorithm::Integrated)
    score(vol, S_Bkg, S_Obj, size, ρ, algorithm::Integrated)

Given a volume `vol` of interest, use the integrated intensity
approach to find the total number of object voxels/volume/mass 
contained in `vol`

# Arguments
- vol: region of interest
- S_Bkg: pure background signal intensity in the `vol`
- S_Obj: pure object signal intensity in the `vol`
- size: size of the voxels in the given `vol`
- ρ: density of the given object of interest contained within
    the `vol`
- algorithm: the `Integrated` algorithm which accounts for
    noise and partial volume effect by integrating the
    the intensity of the entire volume `vol`

# Citation
'Accurate quantification of vessel cross-sectional area using CT
angiography: a simulation study' [Molloi, Johnson, Ding, Lipinski] 
(DOI: 10.1007/s10554-016-1007-9)
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
