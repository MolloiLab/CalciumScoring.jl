### A Pluto.jl notebook ###
# v0.19.32

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
using CalciumScoring: score, Agatston

# ╔═╡ afcf6b72-8907-462a-a21c-1d523c647623
using Statistics: mean

# ╔═╡ ee4aa1d5-2b94-42a5-b85a-89534064b220
using ImageMorphology: dilate

# ╔═╡ 61c966fb-1ba1-40b9-a783-2442f2a42871
using DataFrames: DataFrame

# ╔═╡ 9a9bb5c2-6a1b-4d47-9428-014bc2701ce8
include(joinpath(pwd(), "utils.jl")); # helper functions for creating phantoms

# ╔═╡ bee9b716-27d3-47ca-bc2f-c710422cfda2
md"""
# Set up
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
	[Previously]((00) Getting Started.jl), we introduced the CalciumScoring.jl package. This notebook will examine how to use the Agatston scoring method.
"""

# ╔═╡ a662a8cb-002f-4f00-904e-684576029fa9
md"""
## Create Simulated CT Images
This phantom contains calcium rods of known densities, but the goal is to test our calcium scoring accuracy by quantifying these inserts without directly referencing their true densities. After creating the phantom, its visualizations are presented, and the intensity of each calcium insert is measured. Using our `score()`ing function we can then predict the calcium in each insert and compare that to the ground truth.
"""

# ╔═╡ 85723662-70ed-4fc6-b438-8c418ac63f1b
size_3d = (512, 512, 3)

# ╔═╡ 82297cd0-66b1-45e8-be88-f29c677be5e2
begin
	densities_meas = [0.125, 0.275, 0.425] # mg/cm^3
	phantom_meas, masks_meas = create_calcium_measurement_phantom(size_3d, densities_meas; energy = 120)
end;

# ╔═╡ 171a80e4-e587-4378-8c2e-3aa79b582bd8
md"""
**Dilation of Masks for Measurement Phantom**

In the context of image processing, dilation is a mathematical morphological operation that enlarges the boundaries of the foreground object. 

In this specific step, we dilate the masks associated with the measurement phantom. The primary motivation behind this dilation is to circumvent the partial volume effect that may arise from the surrounding tissue. The partial volume effect can cause blurring and if we don't include the surrounding tissue, some of the calcium intensity can be blurred into the surrounding background tissue and be lost in our measurements
"""

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

# ╔═╡ 70a5dc1e-febd-4ed1-8d18-93653a94f317
md"""
# Agatston (Mass) Score

**Agatston Score**

Using the traditional Agatston score, we can now quantify the amount of calcium contained within the region of interest. We will not need to assume any previous knowledge about calcium within the coronary artery.

**Agatston Mass**

However, the Agatston score returns an arbitrary number, and the direct correlation to physical quantities isn't clear. Luckly, the conversion into a mass score is straightforward forward and this is where the calibration rod comes in. If we measure the intensity of the calibration rod with a known calcium density, we can then correlate this to the unknown calcium in the coronary artery.

CalciumScoring.jl takes care of most of this for us, with the same `score()` function, simply including a `mass_cal_factor` in the arguments. The mass calibration factor is derived by finding the calibration rod's mean intensity and dividing the rod's density by this mean intensity.
"""

# ╔═╡ df7ded9f-5550-429d-a93f-17355e52f376
md"""
!!! warning "Agatston Scoring Threshold"

	One important note is that the threshold of 130 HU is energy dependent. The traditional technique was defined at 120 kV, but recent studies have shown how to adapt the threshold for use in other tube voltages. This is important for reducing patient dose. 

	`score()` accounts for this providing an optional argument `kv = 120`. This assumes that the input image was acquired at a tube voltage of 130, but the user can adjust this and the algorithm automatically updates the threshold and weighting factors accordingly.
"""

# ╔═╡ f60562c3-4837-4a64-94a3-b6993ae268ef
md"""
## Low Density
"""

# ╔═╡ e1ac238b-8446-48ef-8f6a-4f3390a56ca3
ρ_rod = 0.2 # mg/mm^3

# ╔═╡ 009280e2-04a8-4a55-9ae7-d6d551d8ac01
mass_cal_factor = ρ_rod / mean(200)

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
	We can see that the mass score returns a value that is highly accurate as long as the density is high enough. Compare this to the lowest density mass score with the ground truth mass and we see nearly 50% underestimation. This is not a surprise and this is one of the limitations of Agatston scoring. The Agatston scoring approach applies an arbitrary threshold to the input array, and any voxel with an intensity below ``130`` (Hounsfield Units, HU) is removed from the calculation.

	This is one of the main limitations of Agatston scoring and some of the motivations behind the recent "calcium mass" approaches.
"""

# ╔═╡ 4327c28c-473d-40a5-8530-b3b7dd86b6f8
md"""
# Next Steps
We just demonstrated how `score()` can be used with the `Agatston()` algorithm. This is the traditional calcium scoring algorithm. Check out the [Volume Fraction Calcium Mass]((02) Volume Fraction.jl) tutorial to see how to implement a more recent approach with various benefits.
"""

# ╔═╡ Cell order:
# ╟─bee9b716-27d3-47ca-bc2f-c710422cfda2
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
# ╟─5bbe69ff-3326-4e21-b9ae-3cd61a847006
# ╟─a402159d-1625-43d8-a908-d73da3e84282
# ╟─a662a8cb-002f-4f00-904e-684576029fa9
# ╠═85723662-70ed-4fc6-b438-8c418ac63f1b
# ╠═82297cd0-66b1-45e8-be88-f29c677be5e2
# ╟─171a80e4-e587-4378-8c2e-3aa79b582bd8
# ╠═f8547b8c-53ee-46b1-94a9-3b2fc5a28bc0
# ╟─540c1c19-d5dc-431f-8936-55a8f6cac269
# ╟─bc513431-a289-4dc4-a1bb-6889fe4085df
# ╟─8fb9dc5a-6a8b-4c68-856b-31950e61a2f2
# ╠═ef0e3fdd-7a57-4d6d-a8c3-c2ceb889de6a
# ╠═f3459701-0440-45e7-b204-a2017855d565
# ╟─70a5dc1e-febd-4ed1-8d18-93653a94f317
# ╟─df7ded9f-5550-429d-a93f-17355e52f376
# ╟─f60562c3-4837-4a64-94a3-b6993ae268ef
# ╠═e1ac238b-8446-48ef-8f6a-4f3390a56ca3
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
