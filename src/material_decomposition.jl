### A Pluto.jl notebook ###
# v0.19.18

using Markdown
using InteractiveUtils

# ╔═╡ cf470008-87ba-11ed-3976-23a40998fc3a
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate("..")
    using Revise, PlutoUI, CalciumScoring, LsqFit
end

# ╔═╡ 78432aae-91d8-4d2c-8bed-315ac9e6effb
TableOfContents()

# ╔═╡ 63a2221a-bf1c-4920-be3b-9fc970d88f9f
md"""
# Material Decomposition
"""

# ╔═╡ 9e0cb45a-816b-4fd2-b678-82856c3d7329
struct MaterialDecomposition <: CalciumScore end

# ╔═╡ 044f521c-cc42-4ade-9d34-a2a60dabe546
md"""
## Calibration (Curve Fit)

Let's apply the calibration formula found in [An accurate method for direct dual-energy calibration and decomposition](https://www.researchgate.net/publication/20771282_An_accurate_method_for_direct_dual-energy_calibration_and_decomposition)

```math
\begin{aligned}	F = \frac{a_o + a_1x + a_2y + a_3x^2 + a_4xy + a_5y^2}{1 + b_1x + b_2y} 
\end{aligned}
\tag{1}
```

We can rename all the ``\alpha`` and ``b`` parameters to align more closely with the Julia code we will write. *Note: A minimum of 9 measurements (calibration points) is required for accurate calibration.*

```math
\begin{aligned}	F = \frac{p_{1} + p_{2} + p_{3} + p_{4}x^2 + p_{5}xy + p_{6}y^2}{1 + p_{7}x + p_{8}y} 
\end{aligned}
\tag{2}
```

"""

# ╔═╡ 9102343a-f606-4533-8724-f4f43c66507f
"""
	fit_calibration(calculated_calcium_intensities, calcium_densities, p0::AbstractVector=zeros(8))

Calibration equation, where `p0` is a vector of initial parameters that must be fit to underlying data (`calculated_calcium_intensities`) which contains the "high" and "low" energy measurements for various densities of calcium calibration rods (`known_calcium_densities`).

`LsqFit.curve_fit` then returns calibrated parameters `p`.

Note: `calculated_calcium_intensities` is expected to contain "low energy" intensity calculations in the first column of the `n x 2` array and "high energy" intensity calculations in the second column of the `n x 2` array.
"""
function fit_calibration(calculated_calcium_intensities::AbstractMatrix, known_calcium_densities::AbstractVector, p0::AbstractVector=zeros(8))
	F(x, p) = (p[1] .+ (p[2] .* x[:, 1]) .+ (p[3] .* x[:, 2]) .+ (p[4] .* x[:, 1].^2) .+ (p[5] .* x[:, 1] .* x[:, 2]) .+ (p[6] .* x[:, 2].^2)) ./ (1 .+ (p[7] .* x[:, 1]) + (p[8] .* x[:, 2]))
	
	p = LsqFit.curve_fit(F, calculated_calcium_intensities, known_calcium_densities, p0).param
	return p
end

# ╔═╡ 841a0e32-5ef4-4ce2-b9df-aaf1ab89f5de
md"""
## Score
"""

# ╔═╡ e06e98bb-afcf-411a-9161-2b1b970823ee
"""
	score(low_energy_intensity, high_energy_intensity, p, alg::MaterialDecomposition)

Given a dual-energy CT image. First, calculate the measured intensity of a region of interest (ROI) with suspected calcium for both low (`low_energy_intensity`) and high (`high_energy_intensity`) energy scans, then utilize previously calibrated parameters (`p`) (see `CalciumScoring.fit_calibration`) to calculate the density of the suspected calcium within the ROI.

"""
function score(low_energy_intensity, high_energy_intensity, p, alg::MaterialDecomposition)
	A = p[1] + (p[2] * low_energy_intensity) + (p[3] * high_energy_intensity) + (p[4] * low_energy_intensity^2) + (p[5] * low_energy_intensity * high_energy_intensity) + (p[6] * high_energy_intensity^2)
	B = 1 + (p[7] * low_energy_intensity) + (p[8] * high_energy_intensity)
	F = A / B
	return F
end

# ╔═╡ 7aa8559f-f099-41d4-bddd-8af6b012a4a5
"""
	score(low_energy_intensity, high_energy_intensity, p, vol_ROI, alg::MaterialDecomposition)

Given a dual-energy CT image. First, calculate the measured intensity of a region of interest (ROI) with suspected calcium for both low (`low_energy_intensity`) and high (`high_energy_intensity`) energy scans, then utilize previously calibrated parameters (`p`) (see `CalciumScoring.fit_calibration`) to calculate the density of the suspected calcium within the ROI. This can then be multiplied by the volume of the ROI (`vol_ROI`) to calculate the mass of the suspected calcium.

"""
function score(low_energy_intensity, high_energy_intensity, p, vol_ROI, alg::MaterialDecomposition)
	A = p[1] + (p[2] * low_energy_intensity) + (p[3] * high_energy_intensity) + (p[4] * low_energy_intensity^2) + (p[5] * low_energy_intensity * high_energy_intensity) + (p[6] * high_energy_intensity^2)
	B = 1 + (p[7] * low_energy_intensity) + (p[8] * high_energy_intensity)
	F = A / B
	mass = F * vol_ROI
	return mass
end

# ╔═╡ Cell order:
# ╠═cf470008-87ba-11ed-3976-23a40998fc3a
# ╠═78432aae-91d8-4d2c-8bed-315ac9e6effb
# ╟─63a2221a-bf1c-4920-be3b-9fc970d88f9f
# ╠═9e0cb45a-816b-4fd2-b678-82856c3d7329
# ╟─044f521c-cc42-4ade-9d34-a2a60dabe546
# ╠═9102343a-f606-4533-8724-f4f43c66507f
# ╟─841a0e32-5ef4-4ce2-b9df-aaf1ab89f5de
# ╠═e06e98bb-afcf-411a-9161-2b1b970823ee
# ╠═7aa8559f-f099-41d4-bddd-8af6b012a4a5
