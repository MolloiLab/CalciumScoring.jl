### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 2407ead4-e7d3-4885-b3ac-d7ef1da1454a
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(temp = true)
	Pkg.add(url = "https://github.com/Dale-Black/CalciumScoring.jl")
	Pkg.add(["PlutoUI", "CairoMakie"])
	Pkg.add(["ImageFiltering", "ImageMorphology", "Noise"])
	Pkg.add("Statistics")
	
	using CalciumScoring
	using PlutoUI, CairoMakie
	using ImageFiltering, ImageMorphology, Noise
	using Statistics
end

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
	[Previously](/docs/(01) Agatston Scoring.jl), we showcased the Agatston scoring method. This notebook will examine how to use the Volume Fraction Calium Mass method. To understand the theory behind this technique, please see [Coronary artery calcium mass measurement based on integrated intensity and volume fraction techniques](https://doi.org/10.1117/1.JMI.10.4.043502)
"""

# ╔═╡ d77638af-8af0-435c-b81b-734ad859cd67
md"""
## Import Packages
First, let's import CalciumScoring.jl, along with CairoMakie.jl for graphs, and PlutoUI.jl for some interactivity. Be aware this can take a long time, especially if this is the first time being downloaded. Future work on this package will focus on improving this.
"""

# ╔═╡ c4e1c8b0-633c-4c3e-b6da-da8183e7783f
TableOfContents()

# ╔═╡ 2c8e78b9-b8fb-42e7-a577-bbc636ff571f
md"""
!!! info "Create Simulated Images"
	Let's quickly rectreate the simulated images from the [previous tutorial](/docs/(01) Agatston Scoring.jl).
"""

# ╔═╡ 081220be-1d71-48e2-8f4f-d89cfe6d2ba2
begin
    img = zeros(300, 300, 6)
    h, k = size(img)[1:2] .// 2
    r_small = 20
    r_large = 50
    mask_small = zeros(Bool, size(img)[1:2])
    mask_large = zeros(Bool, size(img)[1:2])
    for i in axes(img, 1)
        for j in axes(img, 2)
            for z in axes(img, 3)
                if z <= 3
                    if (i-h)^2 + (j-k)^2 <= r_small^2
                        img[i, j, z] = 130
                        mask_small[i, j] = true
                    end
                else
                    if (i-h)^2 + (j-k)^2 <= r_large^2
                        img[i, j, z] = 260
                        mask_large[i, j] = true
                    end
                end
            end
        end
	end
    
    img_noisy = copy(img)
    for z in axes(img, 3)
        img_noisy[:,:,z] = imfilter(img[:,:,z], Kernel.gaussian((10, 10)))
        img_noisy[:,:,z] = mult_gauss(img_noisy[:,:,z])
    end

	for _ in 1:15
		dilate!(mask_small)
	end
	for _ in 1:10
		erode!(mask_large)
	end
end;

# ╔═╡ c606a885-82ee-43d0-9fa7-e3d1bc36a2d4
mask_small_3D = cat(mask_small, mask_small, mask_small; dims = 3);

# ╔═╡ d15fecd3-79b1-4fbf-9df8-8fa45aa6e9da
mask_large_3D = cat(mask_large, mask_large, mask_large; dims = 3);

# ╔═╡ 62e468a3-718a-4adf-93c6-9e15e79b9230
md"""
!!! info

	Remember, the calibration rod has a known calcium density of ``0.1 \frac{mg}{mm^3}``, the voxel spacing is ``1 \ mm \times 1 \ mm \times 1 \ mm``, and the ground truth calcium mass in the coronary artery is ``376.91 \ mg``
"""

# ╔═╡ b7161fa3-4000-4983-900e-89c2121f6857
ρ_rod = 0.2 # mg/mm^3

# ╔═╡ 632bdc3b-8422-4555-aa13-bc9f002e4851
spacing = [1, 1, 1] # mm

# ╔═╡ 0ea05245-a7de-4c98-ab9e-67f9ef196af9
gt_mass = 376.91 # mg

# ╔═╡ ae93019d-66a3-486a-b272-ad9d25a6ed1b
md"""
## Background Mask
This calcium quantification technique does not require any intensity-based thresholding. The mean intensity of the background material is needed to calculate the calcium contained within a region of interest (ROI).

To find the background intensity, we will first create a background ring mask by dilating the dilated mask again and then subtracting it.
"""

# ╔═╡ d51c034a-1f35-443a-b862-aa95a7be18ce
function dilate_recursively(mask, n)
    dilated_mask = copy(mask)
    for _ in 1:n
        dilated_mask = dilate(dilated_mask)
    end
    return dilated_mask
end

# ╔═╡ 52a86bdc-5f06-4ba9-b193-356a24a6b847
background_mask = dilate_recursively(mask_small, 5) - mask_small;

# ╔═╡ 44baf015-d340-4092-8ed1-74ca1690de52
artery_img = img_noisy[:, :, 1:3];

# ╔═╡ 250b430a-fe29-4db3-b6f8-85c1323c10b2
let
	f = Figure()
	z = 1
	mask = getindex.(findall(isone, mask_small), [1 2])
	bkg_mask = getindex.(findall(isone, background_mask), [1 2])

	ax = CairoMakie.Axis(
		f[1, 1],
		title = "Masks Overlayed (z = $(z))"
	)
	heatmap!(artery_img[:, :, z], colormap = :grays)
	scatter!(mask[:, 1], mask[:, 2], color = (:red, 0.1), markersize = 5, label = "object mask")
	scatter!(bkg_mask[:, 1], bkg_mask[:, 2], color = (:blue, 0.1), markersize = 5, label = "background mask")
	axislegend(ax)
	f
end

# ╔═╡ 8acd5e7b-bd9a-4a40-b786-65940efb9962
background_mask_3D = Bool.(cat(background_mask, background_mask, background_mask; dims = 3));

# ╔═╡ df691f58-a424-43b7-a211-750387436f67
bkg_intensity = mean(artery_img[background_mask_3D])

# ╔═╡ df75acb8-46c9-4445-a9f4-fabcdf6a29fc
md"""
!!! info

	The voxel size is simply `1` since `1 x 1 x 1 = 1`
"""

# ╔═╡ a4f4db07-adf3-49f5-b9ef-464f6249dc41
md"""
## Calculation
"""

# ╔═╡ a6cdb21c-42c7-44eb-84a1-f34bb67fc831
voxel_size = spacing[1] * spacing[2] * spacing[3]

# ╔═╡ d56954da-dbed-4392-8796-5943c8f85fd1
rod_img = img_noisy[:, :, 4:6];

# ╔═╡ a9158ed6-b383-438f-86c3-1f4f8aef65a9
calibration_rod_intensity = mean(rod_img[mask_large_3D])

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
We demonstrated how `score()` can be used with the `VolumeFraction()` algorithm. This is the most well-tested calcium scoring algorithm in the library, but check out the [Integrated Calcium Mass](/docs/(03) Integrated.jl) tutorial to see how to implement another approach.
"""

# ╔═╡ Cell order:
# ╟─fe57f9cc-f298-45a6-8e09-327099b955d5
# ╟─df71ccb2-4062-4d56-a887-d1de363beffa
# ╟─cce7e28d-98ed-473b-8a15-add461821b60
# ╟─d77638af-8af0-435c-b81b-734ad859cd67
# ╠═2407ead4-e7d3-4885-b3ac-d7ef1da1454a
# ╠═c4e1c8b0-633c-4c3e-b6da-da8183e7783f
# ╟─2c8e78b9-b8fb-42e7-a577-bbc636ff571f
# ╠═081220be-1d71-48e2-8f4f-d89cfe6d2ba2
# ╠═c606a885-82ee-43d0-9fa7-e3d1bc36a2d4
# ╠═d15fecd3-79b1-4fbf-9df8-8fa45aa6e9da
# ╟─62e468a3-718a-4adf-93c6-9e15e79b9230
# ╠═b7161fa3-4000-4983-900e-89c2121f6857
# ╠═632bdc3b-8422-4555-aa13-bc9f002e4851
# ╠═0ea05245-a7de-4c98-ab9e-67f9ef196af9
# ╟─ae93019d-66a3-486a-b272-ad9d25a6ed1b
# ╠═d51c034a-1f35-443a-b862-aa95a7be18ce
# ╠═52a86bdc-5f06-4ba9-b193-356a24a6b847
# ╠═44baf015-d340-4092-8ed1-74ca1690de52
# ╟─250b430a-fe29-4db3-b6f8-85c1323c10b2
# ╠═8acd5e7b-bd9a-4a40-b786-65940efb9962
# ╠═df691f58-a424-43b7-a211-750387436f67
# ╟─df75acb8-46c9-4445-a9f4-fabcdf6a29fc
# ╟─a4f4db07-adf3-49f5-b9ef-464f6249dc41
# ╠═a6cdb21c-42c7-44eb-84a1-f34bb67fc831
# ╠═d56954da-dbed-4392-8796-5943c8f85fd1
# ╠═a9158ed6-b383-438f-86c3-1f4f8aef65a9
# ╠═53e784fa-2604-4165-a09d-34bc1cc32ac2
# ╟─c267e5de-afb8-4c42-8db1-77edbae9e669
# ╟─ba100057-22e8-45be-a2d7-71d43a921928
