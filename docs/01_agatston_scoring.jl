### A Pluto.jl notebook ###
# v0.19.37

#> [frontmatter]
#> title = "Agatston Scoring"
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

# ╔═╡ fd8dd444-ff11-491e-93bb-1660ef17593b
# ╠═╡ show_logs = false
using Pkg; Pkg.activate("."); Pkg.instantiate()

# ╔═╡ 56a705cd-3791-4404-a13d-9845bee40062
using PlutoUI: TableOfContents, bind, Slider

# ╔═╡ 84bc14d2-fb94-400d-9489-bd9f37594283
using CairoMakie: Figure, Axis, heatmap!, Colorbar

# ╔═╡ c9074763-f0ea-4941-aa3b-9a554658c337
# ╠═╡ show_logs = false
using CalciumScoring: score, Agatston

# ╔═╡ afcf6b72-8907-462a-a21c-1d523c647623
using Statistics: mean

# ╔═╡ ee4aa1d5-2b94-42a5-b85a-89534064b220
using ImageMorphology: dilate, erode

# ╔═╡ 61c966fb-1ba1-40b9-a783-2442f2a42871
using DataFrames: DataFrame

# ╔═╡ 9a9bb5c2-6a1b-4d47-9428-014bc2701ce8
include(joinpath(pwd(), "utils.jl")); # helper functions for creating phantoms

# ╔═╡ bee9b716-27d3-47ca-bc2f-c710422cfda2
md"""
# Overview
"""

# ╔═╡ 5bbe69ff-3326-4e21-b9ae-3cd61a847006
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

# ╔═╡ a402159d-1625-43d8-a908-d73da3e84282
md"""
!!! success "Overview"
	This notebook will examine how to use the Agatston scoring method.
"""

# ╔═╡ e1e7bc0a-7677-43c8-b71f-e76db0e5440f
md"""
## Activate environment
"""

# ╔═╡ 61ec4074-97a5-4920-af0d-5c7f8cb7da1e
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

# ╔═╡ 1926f907-c4f9-4b2c-bc30-d4b3a9be0743
md"""
## Import Packages
"""

# ╔═╡ 326eb4c9-462d-4fec-ad6b-0eeeb291e044
TableOfContents()

# ╔═╡ aa3585ed-2a4e-4ae7-9c02-7fd05cd91242
md"""
# Agatston Scoring

**Agatston Score**

Using the traditional Agatston score, we can quantify the amount of calcium contained within the region of interest. We will not need to assume any previous knowledge about calcium within the coronary artery.

**Agatston Mass**

However, the Agatston score returns an arbitrary number, and the direct correlation to physical quantities isn't clear. Luckly, the conversion into a mass score is straightforward forward and this is where the calibration rod comes in. If we measure the intensity of the calibration rod with a known calcium density, we can then correlate this to the unknown calcium in the coronary artery.

CalciumScoring.jl takes care of most of this for us, with the same `score()` function, simply including a `mass_cal_factor` in the arguments. The mass calibration factor is derived by finding the calibration rod's mean intensity and dividing the rod's density by this mean intensity.
"""

# ╔═╡ a3ba8309-06f1-4260-a3e7-9d5bb8952e04
md"""
## Calibration Phantom
"""

# ╔═╡ 4a9146c1-e4b9-41a9-93c3-887f2ce82c92
md"""
A calibration phantom is essentially a model or object with pre-known properties, predominantly used to calibrate imaging systems. Within this step, we'll be crafting a digital calibration phantom. This phantom will incorporate 1 rod, with a known density of 0.200 mg/cm^2. This calibration process becomes crucial for the subsequent phase when the task is to quantify calcium samples of unknown densities.
"""

# ╔═╡ 85723662-70ed-4fc6-b438-8c418ac63f1b
size_3d = (512, 512, 3)

# ╔═╡ 0ca6a5a9-83a6-4d06-8604-1cb6f7b484f9
begin
    densities_cal = [0.200] # mg/cm^3
	
    phantom_cal, masks_cal = create_calcium_calibration_phantom(
		size_3d, densities_cal; energy = 120, calcium_radius = 25
	)
end;

# ╔═╡ a4eeb30d-86fd-4ff0-9ec3-689f7d36335b
begin
	# Dilate the masks to avoid partial volume effect from surrounding tissue
    masks_cal_erode = Bool.(copy(masks_cal))
    for idx in axes(masks_cal_erode, 4)
        for _ in 1:5
            masks_cal_erode[:, :, :, idx] = erode(masks_cal_erode[:, :, :, idx])
        end
    end
end

# ╔═╡ aad1f9a8-3ff8-492f-b037-3cc9c46c8736
md"""
Select Slice: $(@bind z_cal Slider(axes(phantom_cal, 3), show_value = true))

Select Mask Number: $(@bind mask_num_cal Slider(axes(masks_cal_erode, 4), show_value = true))
"""

# ╔═╡ 74f656ec-54cb-4ac9-b6d1-c1b03c21d0a5
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

# ╔═╡ e6b713c6-c7d9-43f8-819f-2e38a0a1c086
ρ_rod = densities_cal[1] # mg/mm^3

# ╔═╡ 63636479-d1f9-4705-9e84-4b6ded3edad2
calibration_rod_intensity = mean(phantom_cal[masks_cal_erode[:, :, :, 1]])

# ╔═╡ fea38227-450a-49ed-893a-1d74f55c9198
md"""
## Measurement Phantom

The measurement phantom contains calcium inserts of "unknown" mass that we are going to try to measure using various calcium algorithms.

**Dilation of Masks for Measurement Phantom**

In the context of image processing, dilation is a mathematical morphological operation that enlarges the boundaries of the foreground object. 

In this specific step, we dilate the masks associated with the measurement phantom. The primary motivation behind this dilation is to circumvent the partial volume effect that may arise from the surrounding tissue. The partial volume effect can cause blurring and if we don't include the surrounding tissue, some of the calcium intensity can be blurred into the surrounding background tissue and be lost in our measurements
"""

# ╔═╡ 82297cd0-66b1-45e8-be88-f29c677be5e2
begin
	densities_meas = [0.025, 0.100, 0.250] # mg/cm^3
	phantom_meas, masks_meas = create_calcium_measurement_phantom(size_3d, densities_meas; energy = 120)
end;

# ╔═╡ f8547b8c-53ee-46b1-94a9-3b2fc5a28bc0
begin
    # Dilate the masks to avoid partial volume effect from surrounding tissue
    masks_dil_meas = copy(masks_meas)
    for idx in axes(masks_dil_meas, 4)
        for _ in 1:5
            masks_dil_meas[:, :, :, idx] = dilate(masks_dil_meas[:, :, :, idx])
        end
    end
end

# ╔═╡ 540c1c19-d5dc-431f-8936-55a8f6cac269
md"""
Select Slice: $(@bind z_meas Slider(axes(phantom_meas, 3), show_value = true))

Select Mask Number: $(@bind mask_num_meas Slider(axes(masks_meas, 4), show_value = true))
"""

# ╔═╡ bc513431-a289-4dc4-a1bb-6889fe4085df
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

# ╔═╡ 70a5dc1e-febd-4ed1-8d18-93653a94f317
md"""
# Calculations
"""

# ╔═╡ df7ded9f-5550-429d-a93f-17355e52f376
md"""
!!! warning "Agatston Scoring Threshold"

	One important note is that the threshold of 130 HU is energy dependent. The traditional technique was defined at 120 kV, but recent studies have shown how to adapt the threshold for use in other tube voltages. This is important for reducing patient dose. 

	`score()` accounts for this providing an optional argument `kv = 120`. This assumes that the input image was acquired at a tube voltage of 130, but the user can adjust this and the algorithm automatically updates the threshold and weighting factors accordingly.
"""

# ╔═╡ 8fb9dc5a-6a8b-4c68-856b-31950e61a2f2
md"""
## Ground Truth Mass
"""

# ╔═╡ ef0e3fdd-7a57-4d6d-a8c3-c2ceb889de6a
spacing = [0.5, 0.5, 0.5] # mm

# ╔═╡ f3459701-0440-45e7-b204-a2017855d565
begin
	voxel_vol = 0.5^3
	vol = sum(masks_meas[:, :, :, 1]) * voxel_vol
	gt_mass = densities_meas .* vol # mg
end

# ╔═╡ f60562c3-4837-4a64-94a3-b6993ae268ef
md"""
## Low Density
"""

# ╔═╡ 009280e2-04a8-4a55-9ae7-d6d551d8ac01
mass_cal_factor = ρ_rod / calibration_rod_intensity

# ╔═╡ d3750403-27bd-4055-9378-1c54472fa936
insert_low_density = phantom_meas .* masks_dil_meas[:, :, :, 1];

# ╔═╡ 982ffa5f-9b63-4b66-8404-059e7900bd0d
agatson_score_low_density, volume_score_low_density, mass_score_low_density = score(insert_low_density, spacing, mass_cal_factor, Agatston())

# ╔═╡ 2d169fba-fa50-46ef-8d3a-a3ca875bdaa7
md"""
## Medium Density
"""

# ╔═╡ 51a3b79c-38cc-4e54-a221-decb81b127a8
insert_medium_density = phantom_meas .* masks_dil_meas[:, :, :, 2];

# ╔═╡ 3ab6f3fe-1cf7-4846-8ce9-90333085598f
agatson_score_medium_density, volume_score_medium_density, mass_score_medium_density = score(insert_medium_density, spacing, mass_cal_factor, Agatston())

# ╔═╡ 624cd45b-c9b8-4995-9e70-06c0f898182f
md"""
## High Density
"""

# ╔═╡ 64076971-fd4e-4698-a54d-27f1855eeea7
insert_high_density = phantom_meas .* masks_dil_meas[:, :, :, 3];

# ╔═╡ 4e4c6ddf-3b98-4eff-a454-1eb7f471c250
agatson_score_high_density, volume_score_high_density, mass_score_high_density = score(insert_high_density, spacing, mass_cal_factor, Agatston())

# ╔═╡ bcdbee5e-f742-4244-b81a-046e99b15399
mass_pred = [mass_score_low_density, mass_score_medium_density, mass_score_high_density]

# ╔═╡ df284b6b-aefb-4651-93b8-9102eaca8795
md"""
# Results
"""

# ╔═╡ 81badee4-e40e-4a07-bfca-cafebe18a8c1
df_mass = DataFrame(
    "Ground Truth Mass (mg)" => gt_mass,
    "Predicted Mass (mg)" => mass_pred
)

# ╔═╡ e934660e-a1d1-4f74-b1f0-0070a373556c
md"""
!!! info "Agatston Scoring Limitations"
	We can see that the mass score returns a value that is highly accurate as long as the density is high enough. Compare this to the lowest density mass score with the ground truth mass that labeled this region as no calcium (`CAC = 0`). This is not a surprise and this is one of the limitations of Agatston scoring. The Agatston scoring approach applies an arbitrary threshold to the input array, and any voxel with an intensity below ``130`` (Hounsfield Units, HU) is removed from the calculation.

	This is one of the main limitations of Agatston scoring and some of the motivations behind the recent "calcium mass" approaches.
"""

# ╔═╡ 4327c28c-473d-40a5-8530-b3b7dd86b6f8
md"""
# Next Steps
We just demonstrated how `score()` can be used with the `Agatston()` algorithm. This is the traditional calcium scoring algorithm. Check out the [Volume Fraction Calcium Mass](02_volume_fraction.jl) tutorial to see how to implement a more recent approach with various benefits.
"""

# ╔═╡ Cell order:
# ╟─bee9b716-27d3-47ca-bc2f-c710422cfda2
# ╟─5bbe69ff-3326-4e21-b9ae-3cd61a847006
# ╟─a402159d-1625-43d8-a908-d73da3e84282
# ╟─e1e7bc0a-7677-43c8-b71f-e76db0e5440f
# ╠═fd8dd444-ff11-491e-93bb-1660ef17593b
# ╟─61ec4074-97a5-4920-af0d-5c7f8cb7da1e
# ╟─1926f907-c4f9-4b2c-bc30-d4b3a9be0743
# ╠═56a705cd-3791-4404-a13d-9845bee40062
# ╠═84bc14d2-fb94-400d-9489-bd9f37594283
# ╠═c9074763-f0ea-4941-aa3b-9a554658c337
# ╠═afcf6b72-8907-462a-a21c-1d523c647623
# ╠═ee4aa1d5-2b94-42a5-b85a-89534064b220
# ╠═61c966fb-1ba1-40b9-a783-2442f2a42871
# ╠═9a9bb5c2-6a1b-4d47-9428-014bc2701ce8
# ╠═326eb4c9-462d-4fec-ad6b-0eeeb291e044
# ╟─aa3585ed-2a4e-4ae7-9c02-7fd05cd91242
# ╟─a3ba8309-06f1-4260-a3e7-9d5bb8952e04
# ╟─4a9146c1-e4b9-41a9-93c3-887f2ce82c92
# ╠═85723662-70ed-4fc6-b438-8c418ac63f1b
# ╠═0ca6a5a9-83a6-4d06-8604-1cb6f7b484f9
# ╠═a4eeb30d-86fd-4ff0-9ec3-689f7d36335b
# ╟─aad1f9a8-3ff8-492f-b037-3cc9c46c8736
# ╟─74f656ec-54cb-4ac9-b6d1-c1b03c21d0a5
# ╠═e6b713c6-c7d9-43f8-819f-2e38a0a1c086
# ╠═63636479-d1f9-4705-9e84-4b6ded3edad2
# ╟─fea38227-450a-49ed-893a-1d74f55c9198
# ╠═82297cd0-66b1-45e8-be88-f29c677be5e2
# ╠═f8547b8c-53ee-46b1-94a9-3b2fc5a28bc0
# ╟─540c1c19-d5dc-431f-8936-55a8f6cac269
# ╟─bc513431-a289-4dc4-a1bb-6889fe4085df
# ╟─70a5dc1e-febd-4ed1-8d18-93653a94f317
# ╟─df7ded9f-5550-429d-a93f-17355e52f376
# ╟─8fb9dc5a-6a8b-4c68-856b-31950e61a2f2
# ╠═ef0e3fdd-7a57-4d6d-a8c3-c2ceb889de6a
# ╠═f3459701-0440-45e7-b204-a2017855d565
# ╟─f60562c3-4837-4a64-94a3-b6993ae268ef
# ╠═009280e2-04a8-4a55-9ae7-d6d551d8ac01
# ╠═d3750403-27bd-4055-9378-1c54472fa936
# ╠═982ffa5f-9b63-4b66-8404-059e7900bd0d
# ╟─2d169fba-fa50-46ef-8d3a-a3ca875bdaa7
# ╠═51a3b79c-38cc-4e54-a221-decb81b127a8
# ╠═3ab6f3fe-1cf7-4846-8ce9-90333085598f
# ╟─624cd45b-c9b8-4995-9e70-06c0f898182f
# ╠═64076971-fd4e-4698-a54d-27f1855eeea7
# ╠═4e4c6ddf-3b98-4eff-a454-1eb7f471c250
# ╠═bcdbee5e-f742-4244-b81a-046e99b15399
# ╟─df284b6b-aefb-4651-93b8-9102eaca8795
# ╠═81badee4-e40e-4a07-bfca-cafebe18a8c1
# ╟─e934660e-a1d1-4f74-b1f0-0070a373556c
# ╟─4327c28c-473d-40a5-8530-b3b7dd86b6f8
