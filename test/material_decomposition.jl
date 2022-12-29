### A Pluto.jl notebook ###
# v0.19.18

using Markdown
using InteractiveUtils

# ╔═╡ 0411cbbe-87c4-11ed-1b3e-434ee6932f19
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate("..")
    using Revise, PlutoUI, CalciumScoring, Test
end

# ╔═╡ 1723575b-2a73-4c4e-a268-67e60a9651f7
TableOfContents()

# ╔═╡ 29d2496a-50d8-4d1b-ae5f-5566e90ab02e
md"""
# Material Decomposition
"""

# ╔═╡ d240a8fc-ddc2-4bb7-af80-63f3390dca7b
@testset "MaterialDecomposition()" begin
	calculated_intensities = [69.6768 55.2704; 113.7936 88.7968; 201.6592 152.9552; 288.9744 217.8288; 376.816 282.888; 463.416 345.3504; 548.5744 408.76; 633.9584 471.6464; 717.8224 534.9792; 803.4352 596.9072; 886.5104 657.928; 970.7024 719.7872; 1050.8064 780.0448; 1134.1088 841.7904; 1297.008 961.6592; 1379.3536 1022.3536]
	calcium_densities = [25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800]
	ps = fit_calibration(calculated_intensities, calcium_densities)


	predicted_densities = []
	for i in 1:length(calcium_densities)
		append!(
			predicted_densities, 
			score(
				calculated_intensities[i, 1], 
				calculated_intensities[i, 2], 
				ps,
				MaterialDecomposition()
			)
		)
	end
	@test calcium_densities == abs.(round.(predicted_densities, sigdigits=2))
end

# ╔═╡ 43b231cf-8480-4735-832d-527c93821665
@testset "MaterialDecomposition()" begin
	calculated_intensities = [69.6768 55.2704; 113.7936 88.7968; 201.6592 152.9552; 288.9744 217.8288; 376.816 282.888; 463.416 345.3504; 548.5744 408.76; 633.9584 471.6464; 717.8224 534.9792; 803.4352 596.9072; 886.5104 657.928; 970.7024 719.7872; 1050.8064 780.0448; 1134.1088 841.7904; 1297.008 961.6592; 1379.3536 1022.3536]

	calculated_intensities = calculated_intensities .* rand()
	
	calcium_densities = [25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800]
	ps = fit_calibration(calculated_intensities, calcium_densities)


	predicted_densities = []
	for i in 1:length(calcium_densities)
		append!(
			predicted_densities, 
			score(
				calculated_intensities[i, 1], 
				calculated_intensities[i, 2], 
				ps,
				MaterialDecomposition()
			)
		)
	end
	@test calcium_densities == abs.(round.(predicted_densities, sigdigits=2))
end

# ╔═╡ Cell order:
# ╠═0411cbbe-87c4-11ed-1b3e-434ee6932f19
# ╠═1723575b-2a73-4c4e-a268-67e60a9651f7
# ╟─29d2496a-50d8-4d1b-ae5f-5566e90ab02e
# ╠═d240a8fc-ddc2-4bb7-af80-63f3390dca7b
# ╠═43b231cf-8480-4735-832d-527c93821665
