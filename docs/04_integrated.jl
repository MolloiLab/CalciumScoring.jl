### A Pluto.jl notebook ###
# v0.19.26

#> [frontmatter]
#> title = "Integrated"
#> category = "Tutorials"

using Markdown
using InteractiveUtils

# ╔═╡ 129354eb-7262-4f2d-b189-1610f02d2396
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	
	using CalciumScoring
	using PlutoUI, CairoMakie
	using ImageFiltering, ImageMorphology, Noise
	using Statistics
end

# ╔═╡ 33590140-3daa-4b1e-9a5a-5dde7b27d384
md"""
# Integrated
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
	[Previously]((02) Volume Fraction.jl), we showcased the Volume Fraction Calcium Mass method. This notebook will examine how to use the Integrated Calium Mass method. To understand the theory behind this technique, please see [Coronary artery calcium mass measurement based on integrated intensity and volume fraction techniques](https://doi.org/10.1117/1.JMI.10.4.043502)
"""

# ╔═╡ a1ae748b-fdf1-4a3d-b4c8-9acc23601e4b
md"""
## Import Packages
First, let's import CalciumScoring.jl, along with CairoMakie.jl for graphs, and PlutoUI.jl for some interactivity. Be aware this can take a long time, especially if this is the first time being downloaded. Future work on this package will focus on improving this.
"""

# ╔═╡ 5e7df225-3622-4ee0-90d2-f54f16f6e0f1
TableOfContents()

# ╔═╡ 1ba711ae-8779-4437-b211-203e5ae7a3b4
md"""
!!! info "Create Simulated Images"
	Let's quickly recreate the simulated images from the [previous tutorial]((02) Volume Fraction.jl).
"""

# ╔═╡ 5a3e5de5-d1a4-4c90-8de6-e2e1cdf61f71
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

# ╔═╡ 6dfad6c7-85e0-439e-9801-5c2709d63313
mask_small_3D = cat(mask_small, mask_small, mask_small; dims = 3);

# ╔═╡ a08b8a77-7936-489b-8274-f37e21783484
mask_large_3D = cat(mask_large, mask_large, mask_large; dims = 3);

# ╔═╡ 5ef95123-7217-4d73-815a-ecfe58f6ffb4
md"""
!!! info

	Remember, the calibration rod has a known calcium density of ``0.1 \frac{mg}{mm^3}``, the voxel spacing is ``1 \ mm \times 1 \ mm \times 1 \ mm``, and the ground truth calcium mass in the coronary artery is ``376.91 \ mg``
"""

# ╔═╡ 53a8c163-98c4-4149-92be-c2d33785166f
ρ_rod = 0.2 # mg/mm^3

# ╔═╡ 5f2e9e76-2109-403d-b9da-86eb80d96503
spacing = [1, 1, 1] # mm

# ╔═╡ 5193c28f-b985-4542-9205-ec866d979904
gt_mass = 376.91 # mg

# ╔═╡ 5cdc016e-c8de-4869-825a-8a1a4a2fdc2c
md"""
## Calibration Curve
Unlike the Volume Fraction method, Integrated Calcium Mass requires a calibration curve from various Calibration rods. This makes it more complex but allows for potentially better accuracy than the Volume Fraction Calcium Mass method. A calibration rod and soft tissue ROI can create the calibration curve if only one calibration rod is included.
"""

# ╔═╡ c9d14cfd-0837-4ae3-9ebf-4871dbf58dbc
calibration_img = img_noisy[:, :, 4:end];

# ╔═╡ 714c4020-d45a-4b3d-85bb-90c9fee62074
function offset_mask(mask, offset)
    mask_large2 = zeros(Bool, size(mask))
    for idx in findall(isone, mask)
        new_i = idx[1] - offset[1]
        new_j = idx[2] - offset[2]
        
        # Check if the new indices are within the bounds
        if 1 <= new_i <= size(mask, 1) && 1 <= new_j <= size(mask, 2)
            mask_large2[new_i, new_j] = true
        end
    end
    return mask_large2
end



# ╔═╡ 9df11db6-ad71-401a-900f-3a9b480bdc16
begin
	offset = [80, 80]
	mask_large2 = offset_mask(mask_large, offset)
end;

# ╔═╡ 281cd28a-91fd-4fb8-aecd-da388625002e
mask_large2_3D = cat(mask_large2, mask_large2, mask_large2; dims = 3);

# ╔═╡ c036a512-41fd-4dc2-922f-5ebd4501b748
let
	z = 1
	title = "Cross Section of Calibration Rod"
	mask = getindex.(findall(isone, mask_large), [1 2])
	mask2 = getindex.(findall(isone, mask_large2), [1 2])
	
	f = Figure()

	ax = CairoMakie.Axis(
		f[1, 1],
		title = "Mask Overlayed (z = $(z))",
		titlesize = 20
	)
	heatmap!(calibration_img[:, :, z], colormap = :grays)
	scatter!(mask[:, 1], mask[:, 2], color = :red, markersize = 1.3, label = "calibration mask")
	scatter!(mask2[:, 1], mask2[:, 2], color = :orange, markersize = 1.3, label = "soft tissue mask")

	axislegend(ax)
	f
end

# ╔═╡ c3bc29b4-30f7-41e3-aadc-d1908f9d6856
xs = [0, ρ_rod]

# ╔═╡ e97121fd-a762-4afd-bb43-cbedf538d43f
ys = [mean(calibration_img[mask_large2_3D]), mean(calibration_img[mask_large_3D])]

# ╔═╡ dce61794-6faf-435d-b3ab-676d34df4f8a
let
	f = Figure()

	ax = CairoMakie.Axis(
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

# ╔═╡ 60f08fb4-3268-43c0-9a49-bf30e6cde8b8
md"""
## Background Mask
This calcium quantification technique does not require any intensity-based thresholding. The mean intensity of the background material is needed to calculate the calcium contained within a region of interest (ROI). Similar to the Volume Fraction method

To find the background intensity, we will first create a background ring mask by dilating the dilated mask again and then subtracting it.
"""

# ╔═╡ 2aef781e-04ed-4177-a4dc-c9609129d4cf
function dilate_recursively(mask, n)
    dilated_mask = copy(mask)
    for _ in 1:n
        dilated_mask = dilate(dilated_mask)
    end
    return dilated_mask
end

# ╔═╡ 1352ee5a-0c77-4155-8fd6-c4a372507688
background_mask = dilate_recursively(mask_small, 5) - mask_small;

# ╔═╡ 1c13630e-0ced-492d-a82c-2be6a1080183
artery_img = img_noisy[:, :, 1:3];

# ╔═╡ 4c27fcab-a299-4e30-ab2d-ce0394bff9a9
let
	z = 1
	f = Figure()
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

# ╔═╡ 1fe40c15-d24d-4fb6-8c4b-c6d05d0d3de4
background_mask_3D = Bool.(cat(background_mask, background_mask, background_mask; dims = 3));

# ╔═╡ 7c267a32-24a3-486d-b95d-039470b6e1a7
bkg_intensity = mean(artery_img[background_mask_3D])

# ╔═╡ d3071e3e-2df0-4877-9418-414685c7a023
md"""
!!! info

	The voxel size is simply `1` since `1 x 1 x 1 = 1`
"""

# ╔═╡ 9a657d01-4a8d-42d8-bcf2-66a987e8a223
md"""
## Calculation
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
# ╟─33590140-3daa-4b1e-9a5a-5dde7b27d384
# ╟─6859b0cf-e18f-4032-9297-b0198b6f4a17
# ╟─6b54bcdd-2db0-4b75-acd9-1f1262d2f779
# ╟─a1ae748b-fdf1-4a3d-b4c8-9acc23601e4b
# ╠═129354eb-7262-4f2d-b189-1610f02d2396
# ╠═5e7df225-3622-4ee0-90d2-f54f16f6e0f1
# ╟─1ba711ae-8779-4437-b211-203e5ae7a3b4
# ╠═5a3e5de5-d1a4-4c90-8de6-e2e1cdf61f71
# ╠═6dfad6c7-85e0-439e-9801-5c2709d63313
# ╠═a08b8a77-7936-489b-8274-f37e21783484
# ╟─5ef95123-7217-4d73-815a-ecfe58f6ffb4
# ╠═53a8c163-98c4-4149-92be-c2d33785166f
# ╠═5f2e9e76-2109-403d-b9da-86eb80d96503
# ╠═5193c28f-b985-4542-9205-ec866d979904
# ╟─5cdc016e-c8de-4869-825a-8a1a4a2fdc2c
# ╠═c9d14cfd-0837-4ae3-9ebf-4871dbf58dbc
# ╠═714c4020-d45a-4b3d-85bb-90c9fee62074
# ╠═9df11db6-ad71-401a-900f-3a9b480bdc16
# ╠═281cd28a-91fd-4fb8-aecd-da388625002e
# ╟─c036a512-41fd-4dc2-922f-5ebd4501b748
# ╠═c3bc29b4-30f7-41e3-aadc-d1908f9d6856
# ╠═e97121fd-a762-4afd-bb43-cbedf538d43f
# ╟─dce61794-6faf-435d-b3ab-676d34df4f8a
# ╠═6c234d50-1830-4c78-bd79-19576fc25334
# ╠═52dcddee-f62c-4a63-abde-52ba264028dd
# ╠═f7d6c357-7c12-4db7-a369-bba0692d07b6
# ╟─d1ee5640-6086-413a-afa8-5ac9654d4625
# ╠═76b91e86-97b8-4726-ae2e-bc66a7e40767
# ╠═3885184f-1045-4265-800c-aba7662125ee
# ╟─60f08fb4-3268-43c0-9a49-bf30e6cde8b8
# ╠═2aef781e-04ed-4177-a4dc-c9609129d4cf
# ╠═1352ee5a-0c77-4155-8fd6-c4a372507688
# ╠═1c13630e-0ced-492d-a82c-2be6a1080183
# ╟─4c27fcab-a299-4e30-ab2d-ce0394bff9a9
# ╠═1fe40c15-d24d-4fb6-8c4b-c6d05d0d3de4
# ╠═7c267a32-24a3-486d-b95d-039470b6e1a7
# ╟─d3071e3e-2df0-4877-9418-414685c7a023
# ╟─9a657d01-4a8d-42d8-bcf2-66a987e8a223
# ╠═ac808783-5043-40a1-8270-dea1b310600d
# ╠═2a050381-6b98-449a-b2f5-42d5be899324
# ╠═04772e69-3b59-4f2a-81d8-ab1849f72c9e
# ╠═de44e934-9aed-4be5-9446-2fe981b7d7a9
# ╟─96dd09d6-fbb0-43f5-94d1-f1f0b0391074
# ╟─7e965e51-494c-4796-b06b-f301b02538fb
