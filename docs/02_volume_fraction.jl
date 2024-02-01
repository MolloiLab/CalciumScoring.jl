### A Pluto.jl notebook ###
# v0.19.37

#> [frontmatter]
#> title = "Volume Fraction"
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

# ╔═╡ 7ae3a705-bc9b-4514-a179-4a82f08df9e6
# ╠═╡ show_logs = false
using Pkg; Pkg.activate("."); Pkg.instantiate()

# ╔═╡ efe79ac3-fe84-47a7-9ada-7357f4ed97d7
using PlutoUI: TableOfContents, bind, Slider

# ╔═╡ 4bfb1bee-9a91-41c0-ae1b-096a3a0207fe
using CairoMakie: Figure, Axis, heatmap!, Colorbar

# ╔═╡ a0733e33-4b7e-4494-9bb3-13206e3b83af
using CalciumScoring: score, VolumeFraction

# ╔═╡ e85a8141-30db-48ff-b84e-4e3f68841d90
using Statistics: mean

# ╔═╡ 3bebec59-265d-4a28-aa25-4894af7075cc
using ImageMorphology: dilate, erode

# ╔═╡ 46604fc3-6cb0-4c7c-8b61-fe8baae5d6a6
using DataFrames: DataFrame

# ╔═╡ 34a80e1a-2a95-4b8d-89bc-e73d6ad200c3
include(joinpath(pwd(), "utils.jl")); # helper functions for creating phantoms

# ╔═╡ ff9845e8-ce4b-44a5-90b4-038335d36fd8
md"""
# Overview
"""

# ╔═╡ df71ccb2-4062-4d56-a887-d1de363beffa
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
            class C unique1;
            class B,D,E unique2;
    </div>
</body>
</html>
"""

# ╔═╡ cce7e28d-98ed-473b-8a15-add461821b60
md"""
!!! success "Overview"
	[Previously](01_agatston_scoring.jl), we showcased the Agatston scoring method. This notebook will examine how to use the Volume Fraction Calium Mass method. To understand the theory behind this technique, please see [Coronary artery calcium mass measurement based on integrated intensity and volume fraction techniques](https://doi.org/10.1117/1.JMI.10.4.043502)
"""

# ╔═╡ 7bcef3ec-fc9f-4430-be69-01d17c9ec403
md"""
## Activate environment
"""

# ╔═╡ 059ceb3f-6ad0-4414-9063-0632623ca392
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

# ╔═╡ f6379941-0639-438f-8eac-1d926775bdde
md"""
## Import Packages
"""

# ╔═╡ c4e1c8b0-633c-4c3e-b6da-da8183e7783f
TableOfContents()

# ╔═╡ fe57f9cc-f298-45a6-8e09-327099b955d5
md"""
# Volume Fraction
"""

# ╔═╡ 2c8e78b9-b8fb-42e7-a577-bbc636ff571f
md"""
!!! info
	Let's quickly recreate the simulated images from the [previous tutorial]((01) Agatston Scoring.jl) along with the ground truth mass calculations
"""

# ╔═╡ 627dca2d-ce48-4a0f-8cbd-7482f53cd21f
md"""
## Calibration Phantom
"""

# ╔═╡ 1388c359-5d41-46ee-acfd-4d159936953c
begin
    densities_cal = [0.200] # mg/cm^3
	size_3d = (512, 512, 3)
	
    phantom_cal, masks_cal = create_calcium_calibration_phantom(
		size_3d, densities_cal; energy = 120, calcium_radius = 25
	)
end;

# ╔═╡ 2921276a-a15a-4eae-9147-b6d5dec2672d
begin
	# Dilate the masks to avoid partial volume effect from surrounding tissue
    masks_cal_erode = Bool.(copy(masks_cal))
    for idx in axes(masks_cal_erode, 4)
        for _ in 1:5
            masks_cal_erode[:, :, :, idx] = erode(masks_cal_erode[:, :, :, idx])
        end
    end
end

# ╔═╡ fb13404c-fe9b-44e5-ab7b-8f3349a7a082
md"""
Select Slice: $(@bind z_cal Slider(axes(phantom_cal, 3), show_value = true))

Select Mask Number: $(@bind mask_num_cal Slider(axes(masks_cal_erode, 4), show_value = true))
"""

# ╔═╡ 2b5c4ab4-92d9-4353-90cd-a5ad9261432a
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

# ╔═╡ 5fe2138f-2170-4933-9ef0-fd1a921411bd
ρ_rod = densities_cal[1] # mg/mm^3

# ╔═╡ 98111958-9a44-4097-a998-9ef1cc364a6d
calibration_rod_intensity = mean(phantom_cal[masks_cal_erode[:, :, :, 1]])

# ╔═╡ 8c51ece9-de81-446e-9232-d2dc71326dfd
md"""
## Measurement Phantom
"""

# ╔═╡ 081220be-1d71-48e2-8f4f-d89cfe6d2ba2
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

# ╔═╡ f253a3ec-988d-4627-bbe4-da6a023901cc
md"""
For the Volume Fraction method, we not only need to segment the coronary artery (calcium insert), we also need to get a calculation of the surrounding background tissue. To do this, we will dilate the mask and subtract it, which will leave just a ring of background mask
"""

# ╔═╡ 6ddd6b81-e048-4ed6-9bb5-c6d2adad3ac5
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

# ╔═╡ a79d2779-2563-4b83-8e35-d97dbe5faad7
md"""
Select Slice: $(@bind z_meas Slider(axes(phantom_meas, 3), show_value = true))

Select Mask Number: $(@bind mask_num_meas Slider(axes(masks_meas, 4), show_value = true))
"""

# ╔═╡ 3b1744a3-2434-42f2-96ad-eb0cd6845866
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

# ╔═╡ a4f4db07-adf3-49f5-b9ef-464f6249dc41
md"""
# Calculations
"""

# ╔═╡ f408e8ed-21f5-4c27-8e0a-459fdd9c95db
md"""
## Ground Truth Mass
"""

# ╔═╡ 1ac1a1bb-d469-457f-be90-fe1b92275b90
begin
	spacing = [0.5, 0.5, 0.5] # mm

	voxel_vol = 0.5^3
	vol = sum(masks_meas[:, :, :, 1]) * voxel_vol
	gt_mass = densities_meas .* vol # mg
end

# ╔═╡ a6cdb21c-42c7-44eb-84a1-f34bb67fc831
voxel_size = spacing[1] * spacing[2] * spacing[3]

# ╔═╡ 0c00cc8c-acf2-434c-94f6-0cea88106c81
md"""
## Low Density
"""

# ╔═╡ ab4d7abc-83b8-43d4-ad2a-03ef11119fcb
insert_low_density = phantom_meas[masks_dil_meas[:, :, :, 1]]

# ╔═╡ ead0c005-1327-4e36-a0e2-05c697e781dd
bkg_intensity_low_density = mean(phantom_meas[bg_mask[:, :, :, 1]])

# ╔═╡ a94859b0-7158-4542-b8ce-bac63635e003
vf_mass_low_density = score(insert_low_density, calibration_rod_intensity, bkg_intensity_low_density, voxel_size, ρ_rod, VolumeFraction())

# ╔═╡ 5a14252c-b42f-4a18-a6ec-9d4b1034c286
md"""
## Medium Density
"""

# ╔═╡ dac9434a-80ef-46b8-877e-d6ee00433ebc
insert_medium_density = phantom_meas[masks_dil_meas[:, :, :, 2]]

# ╔═╡ 3717ed21-75ed-4154-be33-dc8e3f5a430a
bkg_intensity_medium_density = mean(phantom_meas[bg_mask[:, :, :, 2]])

# ╔═╡ 531628e9-54bb-4e55-9867-c156b53a2d47
vf_mass_medium_density = score(insert_medium_density, calibration_rod_intensity, bkg_intensity_medium_density, voxel_size, ρ_rod, VolumeFraction())

# ╔═╡ 5cca98a4-d988-4f8f-8cb6-4abc667e04bb
md"""
## High Density
"""

# ╔═╡ a39d293e-0ba9-4c5f-829d-9904ff75dd25
insert_high_density = phantom_meas[masks_dil_meas[:, :, :, 3]]

# ╔═╡ 1d897d4b-7cf3-4673-8fdb-8198197c7c7f
bkg_intensity_high_density = mean(phantom_meas[bg_mask[:, :, :, 3]])

# ╔═╡ 2be68f1e-6d4f-45ff-9f8c-dfde8ddb6992
vf_mass_high_density = score(insert_high_density, calibration_rod_intensity, bkg_intensity_high_density, voxel_size, ρ_rod, VolumeFraction())

# ╔═╡ ece2c7e4-b52b-4d06-b429-629373771462
md"""
# Results
"""

# ╔═╡ 19e69e03-ac1c-46ef-a813-4c6bb25bd3ca
mass_pred = [vf_mass_low_density, vf_mass_medium_density, vf_mass_high_density]

# ╔═╡ e4a90130-4499-4cf3-be7b-870bd66d10f3
df_mass = DataFrame(
    "Ground Truth Mass (mg)" => gt_mass,
    "Predicted Mass (mg)" => mass_pred
)

# ╔═╡ c267e5de-afb8-4c42-8db1-77edbae9e669
md"""
!!! info 
	We can see that `score(..., VolumeFraction())` returns a value that is much more accurate than the corresponding `Agatston` scoring approach, especially for the low density regions.
"""

# ╔═╡ ba100057-22e8-45be-a2d7-71d43a921928
md"""
# Next Steps
We demonstrated how `score()` can be used with the `VolumeFraction()` algorithm. This is the most well-tested calcium scoring algorithm in the library, but check out the [Integrated Calcium Mass](03_integrated.jl) tutorial to see how to implement another approach.
"""

# ╔═╡ Cell order:
# ╟─ff9845e8-ce4b-44a5-90b4-038335d36fd8
# ╟─df71ccb2-4062-4d56-a887-d1de363beffa
# ╟─cce7e28d-98ed-473b-8a15-add461821b60
# ╟─7bcef3ec-fc9f-4430-be69-01d17c9ec403
# ╠═7ae3a705-bc9b-4514-a179-4a82f08df9e6
# ╟─059ceb3f-6ad0-4414-9063-0632623ca392
# ╟─f6379941-0639-438f-8eac-1d926775bdde
# ╠═efe79ac3-fe84-47a7-9ada-7357f4ed97d7
# ╠═4bfb1bee-9a91-41c0-ae1b-096a3a0207fe
# ╠═a0733e33-4b7e-4494-9bb3-13206e3b83af
# ╠═e85a8141-30db-48ff-b84e-4e3f68841d90
# ╠═3bebec59-265d-4a28-aa25-4894af7075cc
# ╠═46604fc3-6cb0-4c7c-8b61-fe8baae5d6a6
# ╠═34a80e1a-2a95-4b8d-89bc-e73d6ad200c3
# ╠═c4e1c8b0-633c-4c3e-b6da-da8183e7783f
# ╟─fe57f9cc-f298-45a6-8e09-327099b955d5
# ╟─2c8e78b9-b8fb-42e7-a577-bbc636ff571f
# ╟─627dca2d-ce48-4a0f-8cbd-7482f53cd21f
# ╠═1388c359-5d41-46ee-acfd-4d159936953c
# ╠═2921276a-a15a-4eae-9147-b6d5dec2672d
# ╟─fb13404c-fe9b-44e5-ab7b-8f3349a7a082
# ╟─2b5c4ab4-92d9-4353-90cd-a5ad9261432a
# ╠═5fe2138f-2170-4933-9ef0-fd1a921411bd
# ╠═98111958-9a44-4097-a998-9ef1cc364a6d
# ╟─8c51ece9-de81-446e-9232-d2dc71326dfd
# ╠═081220be-1d71-48e2-8f4f-d89cfe6d2ba2
# ╟─f253a3ec-988d-4627-bbe4-da6a023901cc
# ╠═6ddd6b81-e048-4ed6-9bb5-c6d2adad3ac5
# ╟─a79d2779-2563-4b83-8e35-d97dbe5faad7
# ╟─3b1744a3-2434-42f2-96ad-eb0cd6845866
# ╟─a4f4db07-adf3-49f5-b9ef-464f6249dc41
# ╟─f408e8ed-21f5-4c27-8e0a-459fdd9c95db
# ╠═a6cdb21c-42c7-44eb-84a1-f34bb67fc831
# ╠═1ac1a1bb-d469-457f-be90-fe1b92275b90
# ╟─0c00cc8c-acf2-434c-94f6-0cea88106c81
# ╠═ab4d7abc-83b8-43d4-ad2a-03ef11119fcb
# ╠═ead0c005-1327-4e36-a0e2-05c697e781dd
# ╠═a94859b0-7158-4542-b8ce-bac63635e003
# ╟─5a14252c-b42f-4a18-a6ec-9d4b1034c286
# ╠═dac9434a-80ef-46b8-877e-d6ee00433ebc
# ╠═3717ed21-75ed-4154-be33-dc8e3f5a430a
# ╠═531628e9-54bb-4e55-9867-c156b53a2d47
# ╟─5cca98a4-d988-4f8f-8cb6-4abc667e04bb
# ╠═a39d293e-0ba9-4c5f-829d-9904ff75dd25
# ╠═1d897d4b-7cf3-4673-8fdb-8198197c7c7f
# ╠═2be68f1e-6d4f-45ff-9f8c-dfde8ddb6992
# ╟─ece2c7e4-b52b-4d06-b429-629373771462
# ╠═19e69e03-ac1c-46ef-a813-4c6bb25bd3ca
# ╠═e4a90130-4499-4cf3-be7b-870bd66d10f3
# ╟─c267e5de-afb8-4c42-8db1-77edbae9e669
# ╟─ba100057-22e8-45be-a2d7-71d43a921928
