### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 3b4c5908-7035-11ed-21ef-f1b5cdb49f4f
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate("..")
    using Revise
    using PlutoUI
    using CalciumScoring
	using CairoMakie
	using Noise
	using ImageFiltering
end

# ╔═╡ a69f6bf2-c37b-481f-ae3e-6893cf7c7b5b
TableOfContents()

# ╔═╡ 8987307f-2021-4f2a-aa64-81059a46eca7
md"""
# Volume Fraction
"""

# ╔═╡ 2a78760b-7546-4acd-a7f4-9a0e3bb451a8
"""
```julia
struct VolumeFraction <: CalciumScore end
```
Lets Julia know (via multiple dispatch) that the algorithm of choice when calculating the calcium score should be the volume fraction algorithm
"""
struct VolumeFraction <: CalciumScore end

# ╔═╡ fc7679be-c412-47e3-8c3f-b6a36321f0ed
"""
```julia
_percentage_calcium(voxel_intensity, hu_calcium, hu_heart_tissue)
```

Find the percentage of calcium within one voxel given the voxel intensity `voxel_intensity`, the Hounsfield unit associated with a known calcium density `hu_calcium`, and the Hounsfield unit associated with background heart tissue `hu_heart_tissue`
"""
function _percentage_calcium(voxel_intensity, hu_calcium, hu_heart_tissue)
	return precentage_calcium = (voxel_intensity - hu_heart_tissue) / (hu_calcium - hu_heart_tissue)
end

# ╔═╡ 7afe5a72-abeb-41c0-9226-fea57ee7d23a
md"""
## 2D
"""

# ╔═╡ e7118615-2c84-4d2b-b0de-6aa3e1d756cb
"""
```julia
score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, calculate the number of voxels that are calcified
"""
function score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, alg::VolumeFraction)
	number_calcium_voxels = []
	for i in axes(vol, 1)
		for j in axes(vol, 2)
			m = vol[i, j]
			percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
			push!(number_calcium_voxels, percent)
		end
	end
	return sum(number_calcium_voxels)
end

# ╔═╡ bb65a2fd-577b-42f3-91bd-f20140c4c5d6
"""
```julia
score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, and the size of an individual pixel/voxel `voxel_size`, calculate the volume of the calcification
"""
function score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, alg::VolumeFraction)
	number_calcium_voxels = []
	for i in axes(vol, 1)
		for j in axes(vol, 2)
			m = vol[i, j]
			percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
			push!(number_calcium_voxels, percent)
		end
	end
	return sum(number_calcium_voxels) * voxel_size * density_calcium
end

# ╔═╡ 35516ea1-527c-4338-9aa0-bfff0ae214ed
"""
```julia
score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, the size of an individual pixel/voxel `voxel_size`, and the density of the previously mentioned calcium `density_calcium`, calculate the mass of the calcification
"""
function score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)
	number_calcium_voxels = []
	for i in axes(vol, 1)
		for j in axes(vol, 2)
			m = vol[i, j]
			percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
			push!(number_calcium_voxels, percent)
		end
	end
	return sum(number_calcium_voxels) * voxel_size * density_calcium
end

# ╔═╡ ad34b8f1-1f3b-4e0f-ad15-3a6ddd4d2c97
md"""
## 3D
"""

# ╔═╡ c4955bb6-eb32-48e6-9ac0-5c0067a0dd39
"""
```julia
score(vol::AbstractArray, hu_calcium, hu_heart_tissue, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, calculate the number of voxels that are calcified
"""
function score(vol::AbstractArray, hu_calcium, hu_heart_tissue, alg::VolumeFraction)
	number_calcium_voxels = []
	for i in axes(vol, 1)
		for j in axes(vol, 2)
			for k in axes(vol, 3)
				m = vol[i, j, k]
				percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
				push!(number_calcium_voxels, percent)
			end
		end
	end
	return sum(number_calcium_voxels)
end

# ╔═╡ 02bef673-e982-41f9-8781-0b11fb9477ac
"""
```julia
score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, and the size of an individual pixel/voxel `voxel_size`, calculate the volume of the calcification
"""
function score(vol::AbstractArray, hu_calcium, hu_heart_tissue, voxel_size, alg::VolumeFraction)
	number_calcium_voxels = []
	for i in axes(vol, 1)
		for j in axes(vol, 2)
			for k in axes(vol, 3)
				m = vol[i, j, k]
				percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
				push!(number_calcium_voxels, percent)
			end
		end
	end
	return sum(number_calcium_voxels) * voxel_size * density_calcium
end

# ╔═╡ 48f4660c-b6b9-402b-b6e1-a0207f1dc468
"""
```julia
score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, the size of an individual pixel/voxel `voxel_size`, and the density of the previously mentioned calcium `density_calcium`, calculate the mass of the calcification
"""
function score(vol::AbstractArray, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)
	number_calcium_voxels = []
	for i in axes(vol, 1)
		for j in axes(vol, 2)
			for k in axes(vol, 3)
				m = vol[i, j, k]
				percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
				push!(number_calcium_voxels, percent)
			end
		end
	end
	return sum(number_calcium_voxels) * voxel_size * density_calcium
end

# ╔═╡ Cell order:
# ╠═3b4c5908-7035-11ed-21ef-f1b5cdb49f4f
# ╠═a69f6bf2-c37b-481f-ae3e-6893cf7c7b5b
# ╟─8987307f-2021-4f2a-aa64-81059a46eca7
# ╟─2a78760b-7546-4acd-a7f4-9a0e3bb451a8
# ╟─fc7679be-c412-47e3-8c3f-b6a36321f0ed
# ╟─7afe5a72-abeb-41c0-9226-fea57ee7d23a
# ╠═e7118615-2c84-4d2b-b0de-6aa3e1d756cb
# ╠═bb65a2fd-577b-42f3-91bd-f20140c4c5d6
# ╠═35516ea1-527c-4338-9aa0-bfff0ae214ed
# ╟─ad34b8f1-1f3b-4e0f-ad15-3a6ddd4d2c97
# ╠═c4955bb6-eb32-48e6-9ac0-5c0067a0dd39
# ╠═02bef673-e982-41f9-8781-0b11fb9477ac
# ╠═48f4660c-b6b9-402b-b6e1-a0207f1dc468
