using Revise, CalciumScoring, LsqFit

struct MaterialDecomposition <: CalciumScore end


""""
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

"""
	fit_calibration(calculated_calcium_intensities, calcium_densities, p0::AbstractVector=zeros(8))

Calibration equation, where `p0` is a vector of initial parameters that must be fit to underlying data (`calculated_calcium_intensities`) which contains the "high" and "low" energy measurements for various densities of calcium calibration rods (`known_calcium_densities`).

`LsqFit.curve_fit` then returns calibrated parameters `p`.

Note: `calculated_calcium_intensities` is expected to contain "low energy" intensity calculations in the first column of the `n x 2` array and "high energy" intensity calculations in the second column of the `n x 2` array.
"""
function fit_calibration(calculated_calcium_intensities::AbstractMatrix, known_calcium_densities::AbstractVector, p0::AbstractVector=zeros(8))
    F(x, p) = (p[1] .+ (p[2] .* x[:, 1]) .+ (p[3] .* x[:, 2]) .+ (p[4] .* x[:, 1] .^ 2) .+ (p[5] .* x[:, 1] .* x[:, 2]) .+ (p[6] .* x[:, 2] .^ 2)) ./ (1 .+ (p[7] .* x[:, 1]) + (p[8] .* x[:, 2]))

    p = LsqFit.curve_fit(F, calculated_calcium_intensities, known_calcium_densities, p0).param
    return p
end

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
