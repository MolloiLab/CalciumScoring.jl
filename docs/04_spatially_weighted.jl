### A Pluto.jl notebook ###
# v0.19.26

#> [frontmatter]
#> title = "Spatially Weighted"
#> category = "Tutorials"

using Markdown
using InteractiveUtils

# ╔═╡ 5356d60d-6bea-4995-b62b-43580455f160
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

# ╔═╡ d06d5523-aa89-417b-8109-ed32cc5c0002
md"""
# Spatially Weighted
"""

# ╔═╡ 6d1e1483-9c73-4815-8b14-b3572ca4a8da
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
            class E unique1;
            class B,C,D unique2;
    </div>
</body>
</html>
"""

# ╔═╡ ac4cfdee-4564-4232-b5b9-5129e077467f
md"""
!!! success "Overview"
	[Previously]((03) Integrated.jl), we showcased the Integrated Calcium Mass method. This notebook will examine how to use the Spatially Weighted Calcium Scoring method. This method attempts to calculate scores equivalent to the Agatston method without thresholding.
"""

# ╔═╡ 56f4299d-630c-48bd-b124-31fc81718a45
md"""
## Import Packages
First, let's import CalciumScoring.jl, along with CairoMakie.jl for graphs, and PlutoUI.jl for some interactivity. Be aware this can take a long time, especially if this is the first time being downloaded. Future work on this package will focus on improving this.
"""

# ╔═╡ 4f56aa03-d299-4b1b-80a3-207e7d4fce99
TableOfContents()

# ╔═╡ f49dd3ea-c6f8-409a-9897-71c7a4bb33cc
md"""
!!! info "Create Simulated Images"
	Let's quickly recreate the simulated images from the [previous tutorial]((03) Integrated.jl).
"""

# ╔═╡ 38083227-3246-45c0-93f8-235f59241c2b
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

# ╔═╡ 9c4c3bb5-0f77-4054-a326-5bac955a3432
mask_small_3D = cat(mask_small, mask_small, mask_small; dims = 3);

# ╔═╡ 24fb41ca-c7cc-4264-90f9-240c4c5587b5
mask_large_3D = cat(mask_large, mask_large, mask_large; dims = 3);

# ╔═╡ ccbda8e1-1d84-4a8c-acb8-b1e2a4242cb0
md"""
!!! info

	Remember, the calibration rod has a known calcium density of ``0.1 \frac{mg}{mm^3}``, the voxel spacing is ``1 \ mm \times 1 \ mm \times 1 \ mm``, and the ground truth calcium mass in the coronary artery is ``376.91 \ mg``
"""

# ╔═╡ 2a16e67d-41db-4032-b723-400146a0ab7f
ρ_rod = 0.2 # mg/mm^3

# ╔═╡ 6f98b16b-b0ce-4284-a036-a2bc76741689
gt_mass = 376.91 # mg

# ╔═╡ 81f3b156-8723-4ee4-aef8-3563652c38fb
md"""
## Calibration Distribution
Unlike the Agatston scoring method, Spatially Weighted Calcium Scoring does not require thresholding but does require a calibration rod with a density of 0.1 ``\frac{mg}{mm^3}``. The rod's mean and standard deviation intensity are then utilized within the SWCS algorithm. Since our calibration rod is at 0.2 ``\frac{mg}{mm^3}``, let's divide the mean in half
"""

# ╔═╡ 3344e286-0600-46b7-bfd4-d3fab6015a51
calibration_img = img_noisy[:, :, 4:end];

# ╔═╡ afb35dc0-377e-474a-8656-67badecffea8
let
	z = 1
	title = "Cross Section of Calibration Rod"
	mask = getindex.(findall(isone, mask_large), [1 2])
	
	f = Figure()

	ax = CairoMakie.Axis(
		f[1, 1],
		title = "Mask Overlayed (z = $(z))",
		titlesize = 20
	)
	heatmap!(calibration_img[:, :, z], colormap = :grays)
	scatter!(mask[:, 1], mask[:, 2], color = :red, markersize = 1.3, label = "calibration mask")

	axislegend(ax)
	f
end

# ╔═╡ 8e89655a-ccab-4b05-b072-b04c5169403d
μ, σ = mean(calibration_img[mask_large_3D]) / 2, std(calibration_img[mask_large_3D])

# ╔═╡ d46f7abf-f605-4c73-8dab-d8a9c0cab518
md"""
## Calculation
"""

# ╔═╡ d278b0b6-0876-49b2-814b-38586089d943
artery_img = img_noisy[:, :, 1:3] .* mask_small_3D;

# ╔═╡ ccab840e-3d04-4c18-80fd-51e55fc6664a
swcs = score(artery_img, μ, σ, SpatiallyWeighted())

# ╔═╡ 70a448fa-909e-4131-a87a-3c0c47bbf046
md"""
!!! info 
	If we remember from the [Agatston Scoring]((01) Agatston Scoring.jl) we see that the Agatston score was `77`, which correlated to a mass score of `204` mg. This greatly underestimated the ground truth mass of `376.91` mg. Spatially weighted calcium scoring cannot be used to calculate mass directly, but the spatially weighted calcium score of 393 is much higher than the corresponding Agatston score. This shows the benefit of the spatially weighted approach when dealing with low-density calcifications compared to the Agatston scoring approach. However, it seems like both volume fraction and integrated methods are more robust.
"""

# ╔═╡ Cell order:
# ╟─d06d5523-aa89-417b-8109-ed32cc5c0002
# ╟─6d1e1483-9c73-4815-8b14-b3572ca4a8da
# ╟─ac4cfdee-4564-4232-b5b9-5129e077467f
# ╟─56f4299d-630c-48bd-b124-31fc81718a45
# ╠═5356d60d-6bea-4995-b62b-43580455f160
# ╠═4f56aa03-d299-4b1b-80a3-207e7d4fce99
# ╟─f49dd3ea-c6f8-409a-9897-71c7a4bb33cc
# ╠═38083227-3246-45c0-93f8-235f59241c2b
# ╠═9c4c3bb5-0f77-4054-a326-5bac955a3432
# ╠═24fb41ca-c7cc-4264-90f9-240c4c5587b5
# ╟─ccbda8e1-1d84-4a8c-acb8-b1e2a4242cb0
# ╠═2a16e67d-41db-4032-b723-400146a0ab7f
# ╠═6f98b16b-b0ce-4284-a036-a2bc76741689
# ╟─81f3b156-8723-4ee4-aef8-3563652c38fb
# ╠═3344e286-0600-46b7-bfd4-d3fab6015a51
# ╟─afb35dc0-377e-474a-8656-67badecffea8
# ╠═8e89655a-ccab-4b05-b072-b04c5169403d
# ╟─d46f7abf-f605-4c73-8dab-d8a9c0cab518
# ╠═d278b0b6-0876-49b2-814b-38586089d943
# ╠═ccab840e-3d04-4c18-80fd-51e55fc6664a
# ╟─70a448fa-909e-4131-a87a-3c0c47bbf046
