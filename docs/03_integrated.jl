### A Pluto.jl notebook ###
# v0.19.37

#> [frontmatter]
#> title = "Integrated"
#> category = "Tutorials"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 1fe4bdc9-d58e-428b-bae7-2fbffec1e537
# ╠═╡ show_logs = false
using Pkg; Pkg.activate("."); Pkg.instantiate()

# ╔═╡ 44121e75-c7b6-4b6b-bd79-42c3869b4ab4
using PlutoUI: TableOfContents, bind, Slider

# ╔═╡ c478f1fd-fbdf-4f8e-af48-ffc21866af63
using CairoMakie: Figure, Axis, heatmap!, Colorbar, scatterlines!

# ╔═╡ 2459744e-f901-449d-8f82-d18f1b67dd24
using CalciumScoring: score, VolumeFraction

# ╔═╡ 696f3b0e-b120-4036-965b-091615d952a1
using Statistics: mean

# ╔═╡ 7e46973d-d11d-47e5-ae22-88a4b76318fa
using ImageMorphology: dilate, erode

# ╔═╡ 8473a26b-3b27-4c53-ba11-27d310409339
using DataFrames: DataFrame

# ╔═╡ 9869067f-a3fd-43a3-b408-9d7b08d375c9
include(joinpath(pwd(), "utils.jl")); # helper functions for creating phantoms

# ╔═╡ 2647a56b-6984-4819-a129-8cd498504342
md"""
# Overview
"""

# ╔═╡ 6859b0cf-e18f-4032-9297-b0198b6f4a17
html"""
<!DOCTYPE html>
<html lang="en">
<head>
    <script src="https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.min.js"></script>
    <script>
        mermaid.initialize({startOnLoad:true});
    </script>
</head>
<body>
    <div class="mermaid">
        graph TB
            classDef unique1 fill:#f9d71c,stroke:#333,stroke-width:2px;
            classDef unique2 fill:#8cd9f5,stroke:#333,stroke-width:2px;
            A["Overview: Four Different CAC Quantification Methods"]
            A --> AgatstonRegime
            A --> CalciumMassRegime
            subgraph AgatstonRegime["Agatston Scoring Regime"]
                B["Agatston Scoring"]
                E["Spatially Weighted Calcium Scoring"]
            end
            subgraph CalciumMassRegime["Calcium Mass Quantification Regime"]
                C["Volume Fraction Calcium Mass"]
                D["Integrated Calcium Mass"]
            end
            class D unique1;
            class B,C,E unique2;
    </div>
</body>
</html>
"""

# ╔═╡ 6b54bcdd-2db0-4b75-acd9-1f1262d2f779
md"""
!!! success "Overview"
	[Previously](02_volume_fraction.jl), we showcased the Volume Fraction Calcium Mass method. This notebook will examine how to use the Integrated Calium Mass method. To understand the theory behind this technique, please see [Coronary artery calcium mass measurement based on integrated intensity and volume fraction techniques](https://doi.org/10.1117/1.JMI.10.4.043502)
"""

# ╔═╡ 48b80de0-ec49-40a3-bd0c-f1a38d97629f
md"""
## Activate environment
"""

# ╔═╡ c92539e6-bbf1-4b47-b1ce-f0b0fa034db8
md"""
!!! info
	This documentation should work well as a downloadable tutorial, since it is just Pluto notebook. If you are using this notebook in an isolated environment, you can comment out the above cell 
	```julia
	# using Pkg; Pkg.activate("."); Pkg.instantiate()
	```
	and let Pluto.jl handle the notebook environment automatically.

	Alternatively, you must manually add each of these `import`ed (`using`) packages. For example: 
	
	```julia
	using Pkg
	Pkg.activate(temp = true)
	Pkg.add("CalciumScoring")
	...
	Pkg.add("DataFrames")
	
	using CalciumScoring
	...
	using DataFrames
	```

	You will need to manually locate and include the [utils.jl file](https://github.com/Dale-Black/CalciumScoring.jl/blob/master/docs/utils.jl) to automatically create the phantoms.
"""

# ╔═╡ 47a98c85-2bdc-4d5c-89a0-098421567108
md"""
## Import Packages
"""

# ╔═╡ c6ef9b15-2b10-492a-9b86-d5b693cd4e62
TableOfContents()

# ╔═╡ 1ba711ae-8779-4437-b211-203e5ae7a3b4
md"""
!!! info "Create Simulated Images"
	Let's quickly recreate the simulated images from the [previous tutorial]((02) Volume Fraction.jl).
"""

# ╔═╡ 01e45ba9-1ddc-4009-be74-61bcbcb957e5
md"""
# Integrated
"""

# ╔═╡ ca76a4c0-d71c-4a56-94e8-cc9f87e4c7b1
md"""
## Calibration Phantom
"""

# ╔═╡ 53f679ba-6b32-4566-aa0b-9bd5df83ea78
begin
    densities_cal = [0.025, 0.100, 0.200, 0.400] # mg/cm^3
	size_3d = (512, 512, 3)
	
    phantom_cal, masks_cal = create_calcium_calibration_phantom(
		size_3d, densities_cal; energy = 120, calcium_radius = 25
	)
end;

# ╔═╡ 7e61eafa-00c0-453b-a4d9-b6e87a745fa1
begin
	# Dilate the masks to avoid partial volume effect from surrounding tissue
    masks_cal_erode = Bool.(copy(masks_cal))
    for idx in axes(masks_cal_erode, 4)
        for _ in 1:5
            masks_cal_erode[:, :, :, idx] = erode(masks_cal_erode[:, :, :, idx])
        end
    end
end

# ╔═╡ 1c094dd4-7cd4-4534-aad4-9a099eba09f9
md"""
Select Slice: $(@bind z_cal Slider(axes(phantom_cal, 3), show_value = true))

Select Mask Number: $(@bind mask_num_cal Slider(axes(masks_cal_erode, 4), show_value = true))
"""

# ╔═╡ 89e52901-95f5-4e6c-9b3b-55df53e676ef
let
    mask = masks_cal_erode[:, :, :, mask_num_cal];
    
    f = Figure()
    ax = Axis(
        f[1, 1],
        title = "Measurement Phantom (120 kV)",
        titlesize = 20
    )   
    hm = heatmap!(phantom_cal[:, :, z_cal], colormap = :grays)
    Colorbar(f[1, 2], hm)

    ax = Axis(
        f[1, 3],
        title = "Dilated Masks",
        titlesize = 20
    )   
    hm = heatmap!(phantom_cal[:, :, z_cal], colormap = :grays)
    hm = heatmap!(mask[:, :, z_cal], colormap = (:jet, 0.5))
    f
end

# ╔═╡ 341299d8-6d6a-4374-bf5a-38640251b57e
ρ_rod = densities_cal[1] # mg/mm^3

# ╔═╡ b83d1f21-d107-4b95-82f5-5a63212b92a5
calibration_rod_intensity = mean(phantom_cal[masks_cal_erode[:, :, :, 1]])

# ╔═╡ d7ab71bf-75ca-422f-9cfc-c500dcce3dd9
md"""
## Measurement Phantom
"""

# ╔═╡ 162f9379-a269-4baa-809b-acdedca4018b
begin
	densities_meas = [0.025, 0.100, 0.250] # mg/cm^3
	phantom_meas, masks_meas = create_calcium_measurement_phantom(size_3d, densities_meas; energy = 120)

	# Dilate the masks to avoid partial volume effect from surrounding tissue
    masks_dil_meas = Bool.(copy(masks_meas))
    for idx in axes(masks_dil_meas, 4)
        for _ in 1:5
            masks_dil_meas[:, :, :, idx] = dilate(masks_dil_meas[:, :, :, idx])
        end
    end
end

# ╔═╡ 4eff5ef2-bac3-4fd7-b9af-63aa1dfaaaa8
md"""
For the Volume Fraction method, we not only need to segment the coronary artery (calcium insert), we also need to get a calculation of the surrounding background tissue. To do this, we will dilate the mask and subtract it, which will leave just a ring of background mask
"""

# ╔═╡ b00bd4f8-9910-4fb5-93c9-03843ed03a81
begin
	# Dilate the masks to avoid partial volume effect from surrounding tissue
    _bg_mask = Bool.(copy(masks_dil_meas))
    for idx in axes(_bg_mask, 4)
        for _ in 1:3
            _bg_mask[:, :, :, idx] = dilate(_bg_mask[:, :, :, idx])
        end
    end

	bg_mask = Bool.(_bg_mask .- masks_dil_meas)
end;

# ╔═╡ 43ea9eb5-d4ca-40d8-9543-7599c5e56358
md"""
Select Slice: $(@bind z_meas Slider(axes(phantom_meas, 3), show_value = true))

Select Mask Number: $(@bind mask_num_meas Slider(axes(masks_meas, 4), show_value = true))
"""

# ╔═╡ a7f00860-c6b5-4cc6-a531-71b78b33fd4b
let
    mask_dil = masks_dil_meas[:, :, :, mask_num_meas];
    bg_mask_num = bg_mask[:, :, :, mask_num_meas];
    
    f = Figure()
    ax = Axis(
        f[1, 1],
        title = "Measurement Phantom (120 kV)",
        titlesize = 20
    )   
    hm = heatmap!(phantom_meas[:, :, z_meas], colormap = :grays)
    Colorbar(f[1, 2], hm)

    ax = Axis(
        f[1, 3],
        title = "Masks",
        titlesize = 20
    )   
    hm = heatmap!(phantom_meas[:, :, z_meas], colormap = :grays)
    hm = heatmap!(mask_dil[:, :, z_meas], colormap = (:jet, 0.5))
    hm = heatmap!(bg_mask_num[:, :, z_meas], colormap = (:viridis, 0.2))
    f
end

# ╔═╡ 5cdc016e-c8de-4869-825a-8a1a4a2fdc2c
md"""
# Calibration Curve
Unlike the Volume Fraction method, Integrated Calcium Mass requires a calibration curve from various Calibration rods. This makes it more complex but allows for potentially better accuracy than the Volume Fraction Calcium Mass method. A calibration rod and soft tissue ROI can create the calibration curve if only one calibration rod is included.
"""

# ╔═╡ c3bc29b4-30f7-41e3-aadc-d1908f9d6856
xs = densities_cal

# ╔═╡ e97121fd-a762-4afd-bb43-cbedf538d43f
ys = [mean(phantom_cal[masks_cal_erode[:, :, :, y]]) for y in axes(masks_cal_erode, 4)]

# ╔═╡ dce61794-6faf-435d-b3ab-676d34df4f8a
let
	f = Figure()

	ax = Axis(
		f[1, 1],
		title = "Calibration Curve",
		titlesize = 20
	)
	scatterlines!(xs, ys, color = :red)
	f
end

# ╔═╡ 6c234d50-1830-4c78-bd79-19576fc25334
m = (ys[2] - ys[1]) / (xs[2] - xs[1])

# ╔═╡ 52dcddee-f62c-4a63-abde-52ba264028dd
b = ys[1] - m * xs[1]

# ╔═╡ f7d6c357-7c12-4db7-a369-bba0692d07b6
function calculate_intensity(x)
    return m * x + b
end

# ╔═╡ d1ee5640-6086-413a-afa8-5ac9654d4625
md"""
Now we can use any arbitrary calibration density as our input and calculate the intensity corresponding to the calibration density. With more calibration points, this becomes more accurate. Let's use `0.1` as the density and get the corresponding intensity.
"""

# ╔═╡ 76b91e86-97b8-4726-ae2e-bc66a7e40767
ρ_01 = 0.1

# ╔═╡ 3885184f-1045-4265-800c-aba7662125ee
calcium_intensity = calculate_intensity(ρ_01)

# ╔═╡ 5193c28f-b985-4542-9205-ec866d979904
gt_mass = 376.91 # mg

# ╔═╡ 9a657d01-4a8d-42d8-bcf2-66a987e8a223
md"""
# Calculations
"""

# ╔═╡ ac808783-5043-40a1-8270-dea1b310600d
voxel_size = spacing[1] * spacing[2] * spacing[3]

# ╔═╡ 2a050381-6b98-449a-b2f5-42d5be899324
rod_img = img_noisy[:, :, 4:6];

# ╔═╡ 04772e69-3b59-4f2a-81d8-ab1849f72c9e
alg = Integrated(artery_img[mask_large_3D])

# ╔═╡ de44e934-9aed-4be5-9446-2fe981b7d7a9
integrated_mass = score(bkg_intensity, calcium_intensity, spacing, ρ_01, alg)

# ╔═╡ 96dd09d6-fbb0-43f5-94d1-f1f0b0391074
md"""
!!! info 
	We can see that `score(..., Integrated())` returns a value of about `376.87` mg. Compare this to the true mass of `376.91` mg, and the calcium mass was estimated slightly more accurately than Volume Fraction and considerably more accurately than Agatston scoring. This has profound implications in diagnosing and preventing coronary artery disease, specifically when identifying calcifications that Agatston scoring missed.
"""

# ╔═╡ 7e965e51-494c-4796-b06b-f301b02538fb
md"""
# Next Steps
We just demonstrated how `score()` can be used with the `Integrated()` algorithm. Let's now look at the last algorithm [Spatially Weighted Calcium Scoring]((04) Spatially Weighted.jl).
"""

# ╔═╡ Cell order:
# ╟─2647a56b-6984-4819-a129-8cd498504342
# ╟─6859b0cf-e18f-4032-9297-b0198b6f4a17
# ╟─6b54bcdd-2db0-4b75-acd9-1f1262d2f779
# ╟─48b80de0-ec49-40a3-bd0c-f1a38d97629f
# ╠═1fe4bdc9-d58e-428b-bae7-2fbffec1e537
# ╟─c92539e6-bbf1-4b47-b1ce-f0b0fa034db8
# ╟─47a98c85-2bdc-4d5c-89a0-098421567108
# ╠═44121e75-c7b6-4b6b-bd79-42c3869b4ab4
# ╠═c478f1fd-fbdf-4f8e-af48-ffc21866af63
# ╠═2459744e-f901-449d-8f82-d18f1b67dd24
# ╠═696f3b0e-b120-4036-965b-091615d952a1
# ╠═7e46973d-d11d-47e5-ae22-88a4b76318fa
# ╠═8473a26b-3b27-4c53-ba11-27d310409339
# ╠═9869067f-a3fd-43a3-b408-9d7b08d375c9
# ╠═c6ef9b15-2b10-492a-9b86-d5b693cd4e62
# ╟─1ba711ae-8779-4437-b211-203e5ae7a3b4
# ╟─01e45ba9-1ddc-4009-be74-61bcbcb957e5
# ╟─ca76a4c0-d71c-4a56-94e8-cc9f87e4c7b1
# ╠═53f679ba-6b32-4566-aa0b-9bd5df83ea78
# ╠═7e61eafa-00c0-453b-a4d9-b6e87a745fa1
# ╟─1c094dd4-7cd4-4534-aad4-9a099eba09f9
# ╟─89e52901-95f5-4e6c-9b3b-55df53e676ef
# ╠═341299d8-6d6a-4374-bf5a-38640251b57e
# ╠═b83d1f21-d107-4b95-82f5-5a63212b92a5
# ╟─d7ab71bf-75ca-422f-9cfc-c500dcce3dd9
# ╠═162f9379-a269-4baa-809b-acdedca4018b
# ╟─4eff5ef2-bac3-4fd7-b9af-63aa1dfaaaa8
# ╠═b00bd4f8-9910-4fb5-93c9-03843ed03a81
# ╟─43ea9eb5-d4ca-40d8-9543-7599c5e56358
# ╟─a7f00860-c6b5-4cc6-a531-71b78b33fd4b
# ╟─5cdc016e-c8de-4869-825a-8a1a4a2fdc2c
# ╠═c3bc29b4-30f7-41e3-aadc-d1908f9d6856
# ╠═e97121fd-a762-4afd-bb43-cbedf538d43f
# ╠═dce61794-6faf-435d-b3ab-676d34df4f8a
# ╠═6c234d50-1830-4c78-bd79-19576fc25334
# ╠═52dcddee-f62c-4a63-abde-52ba264028dd
# ╠═f7d6c357-7c12-4db7-a369-bba0692d07b6
# ╟─d1ee5640-6086-413a-afa8-5ac9654d4625
# ╠═76b91e86-97b8-4726-ae2e-bc66a7e40767
# ╠═3885184f-1045-4265-800c-aba7662125ee
# ╠═5193c28f-b985-4542-9205-ec866d979904
# ╟─9a657d01-4a8d-42d8-bcf2-66a987e8a223
# ╠═ac808783-5043-40a1-8270-dea1b310600d
# ╠═2a050381-6b98-449a-b2f5-42d5be899324
# ╠═04772e69-3b59-4f2a-81d8-ab1849f72c9e
# ╠═de44e934-9aed-4be5-9446-2fe981b7d7a9
# ╟─96dd09d6-fbb0-43f5-94d1-f1f0b0391074
# ╟─7e965e51-494c-4796-b06b-f301b02538fb
