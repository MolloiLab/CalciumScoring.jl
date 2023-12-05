### A Pluto.jl notebook ###
# v0.19.32

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
using ImageMorphology: dilate

# ╔═╡ 46604fc3-6cb0-4c7c-8b61-fe8baae5d6a6
using DataFrames: DataFrame

# ╔═╡ 34a80e1a-2a95-4b8d-89bc-e73d6ad200c3
include(joinpath(pwd(), "utils.jl")); # helper functions for creating phantoms

# ╔═╡ ff9845e8-ce4b-44a5-90b4-038335d36fd8
md"""
# Set up
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
	[Previously]((01) Agatston Scoring.jl), we showcased the Agatston scoring method. This notebook will examine how to use the Volume Fraction Calium Mass method. To understand the theory behind this technique, please see [Coronary artery calcium mass measurement based on integrated intensity and volume fraction techniques](https://doi.org/10.1117/1.JMI.10.4.043502)
"""

# ╔═╡ 2c8e78b9-b8fb-42e7-a577-bbc636ff571f
md"""
!!! info
	Let's quickly recreate the simulated images from the [previous tutorial]((01) Agatston Scoring.jl) along with the ground truth mass calculations
"""

# ╔═╡ 081220be-1d71-48e2-8f4f-d89cfe6d2ba2
begin
	size_3d = (512, 512, 3)

	densities_meas = [0.125, 0.275, 0.425] # mg/cm^3
	phantom_meas, masks_meas = create_calcium_measurement_phantom(size_3d, densities_meas; energy = 120)

	# Dilate the masks to avoid partial volume effect from surrounding tissue
    masks_dil_meas = copy(masks_meas)
    for idx in axes(masks_dil_meas, 4)
        for _ in 1:5
            masks_dil_meas[:, :, :, idx] = dilate(masks_dil_meas[:, :, :, idx])
        end
    end
end

# ╔═╡ a79d2779-2563-4b83-8e35-d97dbe5faad7
md"""
Select Slice: $(@bind z_meas Slider(axes(phantom_meas, 3), show_value = true))

Select Mask Number: $(@bind mask_num_meas Slider(axes(masks_meas, 4), show_value = true))
"""

# ╔═╡ 3b1744a3-2434-42f2-96ad-eb0cd6845866
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

# ╔═╡ 1ac1a1bb-d469-457f-be90-fe1b92275b90
begin
	spacing = [0.5, 0.5, 0.5] # mm

	voxel_vol = 0.5^3
	vol = sum(masks_meas[:, :, :, 1]) * voxel_vol
	gt_mass = densities_meas .* vol # mg
end

# ╔═╡ a4f4db07-adf3-49f5-b9ef-464f6249dc41
md"""
## Calculation
"""

# ╔═╡ a6cdb21c-42c7-44eb-84a1-f34bb67fc831
voxel_size = spacing[1] * spacing[2] * spacing[3]

# ╔═╡ 53e784fa-2604-4165-a09d-34bc1cc32ac2
volume_fraction_mass = score(artery_img[mask_large_3D], calibration_rod_intensity, bkg_intensity, voxel_size, ρ_rod, VolumeFraction())

# ╔═╡ c267e5de-afb8-4c42-8db1-77edbae9e669
md"""
!!! info 
	We can see that `score(..., VolumeFraction())` returns a value of about `375.60` mg. Compare this to the true mass of `376.91` mg and we see that the calcium mass was estimated more accurately than Agatston scoring. This has profound implications in the diagnosis and prevention of coronary artery disease, specifically when identifying calcifications that Agatston scoring missed.
"""

# ╔═╡ ba100057-22e8-45be-a2d7-71d43a921928
md"""
# Next Steps
We demonstrated how `score()` can be used with the `VolumeFraction()` algorithm. This is the most well-tested calcium scoring algorithm in the library, but check out the [Integrated Calcium Mass]((03) Integrated.jl) tutorial to see how to implement another approach.
"""

# ╔═╡ Cell order:
# ╟─ff9845e8-ce4b-44a5-90b4-038335d36fd8
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
# ╟─df71ccb2-4062-4d56-a887-d1de363beffa
# ╟─cce7e28d-98ed-473b-8a15-add461821b60
# ╟─2c8e78b9-b8fb-42e7-a577-bbc636ff571f
# ╠═081220be-1d71-48e2-8f4f-d89cfe6d2ba2
# ╟─a79d2779-2563-4b83-8e35-d97dbe5faad7
# ╟─3b1744a3-2434-42f2-96ad-eb0cd6845866
# ╠═1ac1a1bb-d469-457f-be90-fe1b92275b90
# ╟─a4f4db07-adf3-49f5-b9ef-464f6249dc41
# ╠═a6cdb21c-42c7-44eb-84a1-f34bb67fc831
# ╠═53e784fa-2604-4165-a09d-34bc1cc32ac2
# ╟─c267e5de-afb8-4c42-8db1-77edbae9e669
# ╟─ba100057-22e8-45be-a2d7-71d43a921928
