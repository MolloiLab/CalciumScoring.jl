using Revise
using CalciumScoring
using Test

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
