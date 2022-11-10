### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ ce82e0f9-e06b-48fe-8282-354a019ab89d
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate("..")
	using Revise
	using PlutoUI
	using CalciumScoring
end

# ╔═╡ 2d3c5074-1b69-472e-a3c6-5d4c5dbc8e65
TableOfContents()

# ╔═╡ aa43a020-2927-49df-9ade-778532dcd530
md"""
# Integrated Calcium Scoring
"""

# ╔═╡ 6d2364f9-744e-4a67-aed4-1df09a7d5879
md"""
## Integrated Score
"""

# ╔═╡ 1adb881e-b42e-4116-b622-c49ce5530a88
begin
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
end

# ╔═╡ 9f8b0e46-58f6-49f7-93e3-58841b28f93a
md"""
## Volume Score
"""

# ╔═╡ d526845c-0115-4ad1-a936-61e50a55078f
md"""
## Mass Score
"""

# ╔═╡ 422f4642-f1ca-4842-a364-a70f04b939fb
"""
```julia
    score(vol, S_Bkg, S_Obj, algorithm::Integrated)
```

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

# ╔═╡ f246bdbb-1273-4608-be7a-ad02ef0293cb
"""
```julia
    score(vol, S_Bkg, S_Obj, size, algorithm::Integrated)
```

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
function score(S_Bkg, S_Obj, size, algorithm::Integrated)
    I = algorithm.I
    N = algorithm.N

    N_Obj = (I - (N * S_Bkg)) / (S_Obj - S_Bkg)
    Vol_Obj = N_Obj * size[1] * size[2] * size[3]
    return Vol_Obj
end

# ╔═╡ 56312f94-8c5f-4b97-b79a-3049951f0fc7
"""
```julia
    score(vol, S_Bkg, S_Obj, size, ρ, algorithm::Integrated)
```

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
function score(S_Bkg, S_Obj, size, ρ, algorithm::Integrated)
    I = algorithm.I
    N = algorithm.N

    N_Obj = (I - (N * S_Bkg)) / (S_Obj - S_Bkg)
    Vol_Obj = N_Obj * size[1] * size[2] * size[3]
    Mass_Obj = Vol_Obj * ρ
    return Mass_Obj
end

# ╔═╡ Cell order:
# ╠═ce82e0f9-e06b-48fe-8282-354a019ab89d
# ╠═2d3c5074-1b69-472e-a3c6-5d4c5dbc8e65
# ╟─aa43a020-2927-49df-9ade-778532dcd530
# ╟─6d2364f9-744e-4a67-aed4-1df09a7d5879
# ╠═1adb881e-b42e-4116-b622-c49ce5530a88
# ╟─9f8b0e46-58f6-49f7-93e3-58841b28f93a
# ╠═422f4642-f1ca-4842-a364-a70f04b939fb
# ╟─d526845c-0115-4ad1-a936-61e50a55078f
# ╠═f246bdbb-1273-4608-be7a-ad02ef0293cb
# ╠═56312f94-8c5f-4b97-b79a-3049951f0fc7
