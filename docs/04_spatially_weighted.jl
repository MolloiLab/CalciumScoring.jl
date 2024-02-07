### A Pluto.jl notebook ###
# v0.19.37

#> [frontmatter]
#> title = "Spatially Weighted"
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

# ╔═╡ 51ac5d43-ecbf-4885-9570-ac447d31604c
# ╠═╡ show_logs = false
using Pkg; Pkg.activate("."); Pkg.instantiate()

# ╔═╡ a8eaa113-78f7-4348-a6ef-758f9868dfe6
using PlutoUI: TableOfContents, bind, Slider

# ╔═╡ c3d2bccd-f0f1-4353-adbd-b2579e4b75fa
using CairoMakie: Figure, Axis, heatmap!, Colorbar

# ╔═╡ a57ff4cd-ecde-4cc0-ae4a-e78a6fba641a
using CalciumScoring: score, SpatiallyWeighted

# ╔═╡ 33c09680-bd53-499d-8dc2-09545f1b1ea5
using Statistics: mean, std

# ╔═╡ 391a92bb-af9e-42be-ba7c-feb71012711c
using ImageMorphology: dilate, erode

# ╔═╡ ff986002-bf09-4966-bc54-cb803d2ff00c
using DataFrames: DataFrame

# ╔═╡ 5cc569e9-455e-4b95-8eed-7cbf3139dcb2
include(joinpath(pwd(), "utils.jl")); # helper functions for creating phantoms

# ╔═╡ 795190a6-5c79-47e7-936d-d041e8641a37
md"""
# Overview
"""

# ╔═╡ 57bc064f-fb1e-43f0-a380-f7c518bda319
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
            class B unique1;
            class C,D,E unique2;
    </div>
</body>
</html>
"""

# ╔═╡ 0a30f133-e587-407e-804f-bacd99c07865
md"""
!!! success "Overview"
	[Previously](03_integrated.jl), we showcased the Integrated Calcium Mass method. This notebook will examine how to use the Spatially Weighted Calcium Scoring method. To understand the theory behind this technique, please see [Spatially Weighted Coronary Artery Calcium Score and Coronary Heart Disease Events in the Multi-Ethnic Study of Atherosclerosis](https://doi.org/10.1161/CIRCIMAGING.120.011981)

"""

# ╔═╡ ffb5514e-ef5f-4fa5-ae04-2a860e6ed29d
md"""
## Activate environment
"""

# ╔═╡ 1a752fda-c11e-4a56-80c5-fe89b1be7c11
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

# ╔═╡ 227331bb-62ee-4a48-b604-1464ca7b9426
md"""
## Import Packages
"""

# ╔═╡ edd50476-7bb4-461c-8b5a-3dbbbbcde6a6
TableOfContents()

# ╔═╡ e1e94261-b503-4801-a71f-ef6951d85ce1
md"""
# Spatially Weighted

Unlike the Agatston scoring method, Spatially Weighted Calcium Scoring does not require thresholding but does require a calibration rod with a density of 0.1 ``\frac{mg}{cm^3}``. The rod's mean and standard deviation intensity are then utilized within the SWCS algorithm.

"""

# ╔═╡ a72d1057-a772-47a4-82b7-2124f6113bbb
md"""
## Calibration Phantom
"""

# ╔═╡ 70c28823-05c4-423c-9513-20af00fad03b
size_3d = (512, 512, 3)

# ╔═╡ a53b60d6-fc23-46ce-8344-1af54948d67d
begin
    densities_cal = [0.100] # mg/cm^3
	
    phantom_cal, masks_cal = create_calcium_calibration_phantom(
		size_3d, densities_cal; 
		energy = 120, calcium_radius = 25, noise_std = 25
	)
end;

# ╔═╡ 37780bc2-4dc4-472b-814d-b248e93393e4
begin
	# Dilate the masks to avoid partial volume effect from surrounding tissue
    masks_cal_erode = Bool.(copy(masks_cal))
    for idx in axes(masks_cal_erode, 4)
        for _ in 1:5
            masks_cal_erode[:, :, :, idx] = erode(masks_cal_erode[:, :, :, idx])
        end
    end
end

# ╔═╡ 76758ecd-ec36-4e6e-9d84-09695c15d0e4
md"""
Select Slice: $(@bind z_cal Slider(axes(phantom_cal, 3), show_value = true))

Select Mask Number: $(@bind mask_num_cal Slider(axes(masks_cal_erode, 4), show_value = true))
"""

# ╔═╡ a892147a-c4f6-4f96-9625-7e5ce3e1ed0f
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

# ╔═╡ d3ce5915-6ed4-4554-9aff-cf88324cb607
μ = mean(phantom_cal[masks_cal_erode[:, :, :, 1]])

# ╔═╡ 301581df-5600-4d80-b79d-8ca8e9715d1f
σ = std(phantom_cal[masks_cal_erode[:, :, :, 1]])

# ╔═╡ 6653fdeb-9ee7-432c-ae24-cf600c23c1eb
md"""
## Measurement Phantom
"""

# ╔═╡ 67687d36-bd37-4e8d-bf02-e7ce6d2da86e
begin
	densities_meas = [0.025, 0.100, 0.250] # mg/cm^3
	phantom_meas, masks_meas = create_calcium_measurement_phantom(size_3d, densities_meas; energy = 120)
end;

# ╔═╡ 05c53231-2a20-497d-b2cd-b5e4e0ed5114
begin
    # Dilate the masks to avoid partial volume effect from surrounding tissue
    masks_dil_meas = copy(masks_meas)
    for idx in axes(masks_dil_meas, 4)
        for _ in 1:5
            masks_dil_meas[:, :, :, idx] = dilate(masks_dil_meas[:, :, :, idx])
        end
    end
end

# ╔═╡ e4228825-f33b-41a4-bb9f-a1f7c8315e8b
md"""
Select Slice: $(@bind z_meas Slider(axes(phantom_meas, 3), show_value = true))

Select Mask Number: $(@bind mask_num_meas Slider(axes(masks_meas, 4), show_value = true))
"""

# ╔═╡ 6c13c0c0-177e-40f8-ad62-448ef64b7d72
let
    mask_dil = masks_dil_meas[:, :, :, mask_num_meas];
    
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
        title = "Dilated Masks",
        titlesize = 20
    )   
    hm = heatmap!(phantom_meas[:, :, z_meas], colormap = :grays)
    hm = heatmap!(mask_dil[:, :, z_meas], colormap = (:jet, 0.5))
    f
end

# ╔═╡ 9885f2d6-090f-4300-ac25-9659eba97ee7
md"""
# Calculations
"""

# ╔═╡ 2f9f8100-80f7-46c6-b727-7616a8b4cdad
md"""
## Low Density
"""

# ╔═╡ 326644c5-e94d-476f-ac44-c65216ac59bd
insert_low_density = phantom_meas .* masks_dil_meas[:, :, :, 1];

# ╔═╡ a8df8b9b-c30f-44ce-8697-9952ec39eeaf
swcs_low_density = score(insert_low_density, μ, σ, SpatiallyWeighted())

# ╔═╡ 40e82087-c778-41ab-8bac-fdbe13049e85
md"""
## Medium Density
"""

# ╔═╡ b2c2d8ab-d69d-49d4-a714-7830e193b121
insert_medium_density = phantom_meas .* masks_dil_meas[:, :, :, 2];

# ╔═╡ 0546770f-f8a0-4ba3-b3aa-e63e8a7cad04
swcs_medium_density = score(insert_medium_density, μ, σ, SpatiallyWeighted())

# ╔═╡ f1789af6-a7fb-4384-bbd4-37c77b5b67bc
md"""
## High Density
"""

# ╔═╡ d5c1846d-3168-43ab-a685-4a74b2dd518e
insert_high_density = phantom_meas .* masks_dil_meas[:, :, :, 3];

# ╔═╡ 066c72f7-bc33-48fd-9a07-908ad2e81e8d
swcs_high_density = score(insert_high_density, μ, σ, SpatiallyWeighted())

# ╔═╡ 65b5b84f-3c6e-460e-94df-a0a7b6fca06d
md"""
# Results
"""

# ╔═╡ 2261fbf7-8ac5-4aa0-b303-e9060625880f
md"""
!!! info "Recall"

	From the [Agatston Scoring tutorial](01_agatston_scoring.jl), we got pure scoring results (non mass) that were something like `[0.0, 41.25, 369.25]`. We can now see the benefit of spatially weighted calcium scoring below. Namely, the `SpatiallyWeighted` approach is correlated with `Agatston` scoring at higher densities, but at low densities, where `Agatston` scoring results in false-negatives, `SpatiallyWeighted` detecs some calcium.
"""

# ╔═╡ 82201e91-7110-4105-8543-6a418ff78735
agatston = [0.0, 41.25, 369.25]

# ╔═╡ 0f04f1e3-6810-4d7b-8f01-ae59c6ba7c08
swcs = [swcs_low_density, swcs_medium_density, swcs_high_density]

# ╔═╡ f07903f4-3fa0-4c60-a087-b5e12ef00862
df_mass = DataFrame(
	"Agatston Scoring" => agatston,
    "Spatially Weighted Calcium Score" => swcs
)

# ╔═╡ ebda2d13-d461-4736-9e45-817ffd6c4642
md"""
!!! warning "Limitations"

	One of the main limitations of the `SpatiallyWeighted` methods is the arbitrary score that it returns. Unlike `Agatston` scoring, this score is not straightforward to convert into a mass measurement. Although this approach improves upon Agatston scoring in various way it also has its own limitations.
"""

# ╔═╡ Cell order:
# ╟─795190a6-5c79-47e7-936d-d041e8641a37
# ╟─57bc064f-fb1e-43f0-a380-f7c518bda319
# ╟─0a30f133-e587-407e-804f-bacd99c07865
# ╟─ffb5514e-ef5f-4fa5-ae04-2a860e6ed29d
# ╠═51ac5d43-ecbf-4885-9570-ac447d31604c
# ╟─1a752fda-c11e-4a56-80c5-fe89b1be7c11
# ╟─227331bb-62ee-4a48-b604-1464ca7b9426
# ╠═a8eaa113-78f7-4348-a6ef-758f9868dfe6
# ╠═c3d2bccd-f0f1-4353-adbd-b2579e4b75fa
# ╠═a57ff4cd-ecde-4cc0-ae4a-e78a6fba641a
# ╠═33c09680-bd53-499d-8dc2-09545f1b1ea5
# ╠═391a92bb-af9e-42be-ba7c-feb71012711c
# ╠═ff986002-bf09-4966-bc54-cb803d2ff00c
# ╠═5cc569e9-455e-4b95-8eed-7cbf3139dcb2
# ╠═edd50476-7bb4-461c-8b5a-3dbbbbcde6a6
# ╟─e1e94261-b503-4801-a71f-ef6951d85ce1
# ╟─a72d1057-a772-47a4-82b7-2124f6113bbb
# ╠═70c28823-05c4-423c-9513-20af00fad03b
# ╠═a53b60d6-fc23-46ce-8344-1af54948d67d
# ╠═37780bc2-4dc4-472b-814d-b248e93393e4
# ╟─76758ecd-ec36-4e6e-9d84-09695c15d0e4
# ╟─a892147a-c4f6-4f96-9625-7e5ce3e1ed0f
# ╠═d3ce5915-6ed4-4554-9aff-cf88324cb607
# ╠═301581df-5600-4d80-b79d-8ca8e9715d1f
# ╟─6653fdeb-9ee7-432c-ae24-cf600c23c1eb
# ╠═67687d36-bd37-4e8d-bf02-e7ce6d2da86e
# ╠═05c53231-2a20-497d-b2cd-b5e4e0ed5114
# ╟─e4228825-f33b-41a4-bb9f-a1f7c8315e8b
# ╟─6c13c0c0-177e-40f8-ad62-448ef64b7d72
# ╟─9885f2d6-090f-4300-ac25-9659eba97ee7
# ╟─2f9f8100-80f7-46c6-b727-7616a8b4cdad
# ╠═326644c5-e94d-476f-ac44-c65216ac59bd
# ╠═a8df8b9b-c30f-44ce-8697-9952ec39eeaf
# ╟─40e82087-c778-41ab-8bac-fdbe13049e85
# ╠═b2c2d8ab-d69d-49d4-a714-7830e193b121
# ╠═0546770f-f8a0-4ba3-b3aa-e63e8a7cad04
# ╟─f1789af6-a7fb-4384-bbd4-37c77b5b67bc
# ╠═d5c1846d-3168-43ab-a685-4a74b2dd518e
# ╠═066c72f7-bc33-48fd-9a07-908ad2e81e8d
# ╟─65b5b84f-3c6e-460e-94df-a0a7b6fca06d
# ╟─2261fbf7-8ac5-4aa0-b303-e9060625880f
# ╠═82201e91-7110-4105-8543-6a418ff78735
# ╠═0f04f1e3-6810-4d7b-8f01-ae59c6ba7c08
# ╠═f07903f4-3fa0-4c60-a087-b5e12ef00862
# ╟─ebda2d13-d461-4736-9e45-817ffd6c4642
