### A Pluto.jl notebook ###
# v0.19.26

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

# ╔═╡ 2ed53026-7897-465a-a388-3d800ac2de96
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

# ╔═╡ 24de8bc6-e495-49d4-a1b2-673f714bdf9f
md"""
# Agatston Scoring
"""

# ╔═╡ 1e656925-3b85-498d-b141-2066de514aed
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

# ╔═╡ aa529c31-0503-46c9-a027-b911797a7306
md"""
!!! success "Overview"
	[Previously]((01) Getting Started.jl), we introduced the CalciumScoring.jl package. This notebook will examine how to use the Agatston scoring method.
"""

# ╔═╡ 3c866da2-8263-4d4b-93f8-f178d243130d
md"""
## Import Packages
First, let's import CalciumScoring.jl, along with CairoMakie.jl for graphs, and PlutoUI.jl for some interactivity. Be aware this can take a long time, especially if this is the first time being downloaded. Future work on this package will focus on improving this.
"""

# ╔═╡ 9ee13ae7-71b0-447b-8ce3-60e9649439a9
TableOfContents()

# ╔═╡ 5ddaf3d2-7dd6-49ac-934d-cb7bdfd37a43
md"""
## Create Simulated CT Images
Below we will create a simulated 3D CT image. The first 3 slices correspond to a simulated coronary artery containing some unknown amount of calcium. The last 3 slices correspond to a simulated calibration rod of calcium with a known amount of calcium.

We will also create two masks. One mask that contains the entire coronary artery including the region along the boundary of the artery that is affected by partial volume effect (blurring). For the other mask, we will contain only the inner-most part of the calibration rod, avoiding any voxels potentially affected by the partial volume effect.
"""

# ╔═╡ aee04655-a30d-430d-a732-b3897938d740
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

# ╔═╡ f5c7fb00-ee43-44bf-84de-e270f77a0930
@bind z1 PlutoUI.Slider(axes(img_noisy, 3), show_value = true)

# ╔═╡ 58760da5-0b94-4c01-9024-a46fe2e7832a
let
	if z1 in 1:3
		title = "Cross Section of Simulated Coronary Artery"
		mask = getindex.(findall(isone, mask_small), [1 2])
	else
		title = "Cross Section of Calibration Rod"
		mask = getindex.(findall(isone, mask_large), [1 2])
	end
	
	
	f = Figure(resolution = (1200, 800))
	ax = CairoMakie.Axis(
		f[1, 1],
		title = title,
		titlesize = 20
	)
	heatmap!(img_noisy[:, :, z1], colormap = :grays)

	ax = CairoMakie.Axis(
		f[1, 2],
		title = "Mask Overlayed",
		titlesize = 20
	)
	heatmap!(img_noisy[:, :, z1], colormap = :grays)
	scatter!(mask[:, 1], mask[:, 2], color = :red, markersize = 1.3)
	f
end

# ╔═╡ 6616ca9c-5de2-4615-8707-8b7138b9f2e7
md"""
!!! success "Calculate Ground Truth Calcium Mass"
	If assume that the calcium in the coronary artery (slices 1-3) has a density of ``0.1 \frac{mg}{mm^3}``, we can then calculate the mass of this calcium within the entire 3D image. We first assume the voxels are of size ``1 \ mm \times 1 \ mm \times 1 \ mm`` and calculate the volume of the calcium.
	
	```math
	\begin{align*}
	\text{Volume}_\text{calcium} &= \pi \times r^2 \times \text{height} \\
	&= \pi (20 mm)^2(3 mm) \\
	&= 3769.91 mm^3 \\ \\
	
	\text{Mass}_\text{calcium} &= \text{Volume}_\text{calcium} \times \text{Density}_\text{calcium} \\
	&= (3769.91 mm^3)(0.1 \frac{mg}{mm^3}) \\
	&= 376.91 mg
	\end{align*}
	```
"""

# ╔═╡ 46b33f08-05b9-4c26-a9e4-49b0d1dbd0b3
gt_mass = 376.91 # mg

# ╔═╡ db66f0c7-4256-4ffa-9cd3-90766e08381d
md"""
## Agatston Score
We can now quantify the amount of calcium contained within the region of interest, using the traditional Agatston score. We will not need to assume any previous knowledge about calcium within the coronary artery.
"""

# ╔═╡ 98b7ed8c-97d8-4966-a65a-27aabf51ec55
spacing = [1, 1, 1] # mm

# ╔═╡ dbe72760-be53-4f72-a8ef-6b46d7a0b824
mask_small_3D = cat(mask_small, mask_small, mask_small; dims = 3);

# ╔═╡ 611a3d4e-18ec-45db-998c-be40e7f59cff
artery_img = img_noisy[:, :, 1:3] .* mask_small_3D;

# ╔═╡ 8d883827-995a-433d-9143-03b2326c0610
@bind z2 PlutoUI.Slider(axes(artery_img, 3), show_value = true)

# ╔═╡ 161b9ba5-013e-4d7f-86b8-c3aaa53d3206
heatmap(artery_img[:, :, z2], colormap = :grays)

# ╔═╡ 5b8e3935-3dc5-4f7b-b83e-6192f1e285cb
agatson_score, volume_score = score(artery_img, spacing, Agatston())

# ╔═╡ aeb4d61b-c091-47f5-8ed5-bddee40960d2
md"""
## Mass Score
We can see that some amount of calcium was detected. But the Agatston score returns an arbitray number and the direct correlation to physical quantities isn't clear. Luckly, the conversion into a mass score is straight forward and this is where the calibration rod comes in. If we measure the intensity of the calibration rod with a known calcium density, we can then correlate this to the unknown calcium in the coronary artery.

CalciumScoring.jl takes care of most of this for us, with the same `score()` function, simply including a `mass_cal_factor` in the arguments. The mass calibration factor is derived by finding the mean intensity of the calibration rod and dividing the density of the rod by this mean intensity.
"""

# ╔═╡ 1525560b-1ef9-4c9f-b557-0d8f3f37d0cd
mask_large_3D = cat(mask_large, mask_large, mask_large; dims = 3);

# ╔═╡ 0f63ae63-d8cd-4b38-99d0-e53d0d68ea11
rod_img = img_noisy[:, :, 4:6] .* mask_large_3D;

# ╔═╡ dcee2b99-1b32-4945-bd4b-5ad5c07cd8b5
ρ_rod = 0.2 # mg/mm^3

# ╔═╡ 07eccab2-6f71-4eaf-b98a-00bc5eaea693
mass_cal_factor = ρ_rod / mean(rod_img)

# ╔═╡ 0231b3f4-8e2b-4f68-be5f-553714069082
_, _, mass_score = score(artery_img, spacing, mass_cal_factor, Agatston())

# ╔═╡ 3f05439b-afdb-43dc-9246-7c00cdedffb3
md"""
!!! info "Agatston Scoring Limitations"
	We can see that the mass score returns a value of about `204.01` mg. Compare this to the true mass of `376.91` mg and we see that the mass was underestimated by almost half. This is not a surprise and this is one of the limitations of Agatston scoring. The Agatston scoring approach applies an arbitrary threshold to the input array, and any voxel with an intensity below ``130`` (Hounsfield Units, HU) is removed from the calculation.

	The simulated calcium image we created set the voxels to an intensity of exactly 130 HU. After noise was added, a number of those voxels were then below the 130 HU threshold and, therefore not included in the Agatston scoring calculation.

	This is one of the main limitations of Agatston scoring and some of the motivations behind the recent calcium scoring approaches.
"""

# ╔═╡ 48d50aa7-cdf6-48d4-97f5-ed9713d3b325
md"""
!!! warning "Agatston Scoring Threshold"

	One important note is that the threshold of 130 HU is energy dependent. The traditional technique was defined at 120 kV, but recent studies have shown how to adapt the threshold for use in other tube voltages. This is important for reducing patient dose. 

	`score()` accounts for this providing an optional argument `kv = 120`. This assumes that the input image was acquired at a tube voltage of 130, but the user adjust this and the algorithm automatically updates the threshold and weighting factors accordingly.
"""

# ╔═╡ ba707bc7-6b36-4131-97ad-ea757d17d646
md"""
# Next Steps
We just demonstrated how `score()` can be used with the `Agatston()` algorithm. This is the traditional calcium scoring algorithm. Check out the [Volume Fraction Calcium Mass]((02) Volume Fraction.jl) tutorial to see how to implement a more recent approach with various benefits.
"""

# ╔═╡ Cell order:
# ╟─24de8bc6-e495-49d4-a1b2-673f714bdf9f
# ╟─1e656925-3b85-498d-b141-2066de514aed
# ╟─aa529c31-0503-46c9-a027-b911797a7306
# ╟─3c866da2-8263-4d4b-93f8-f178d243130d
# ╠═2ed53026-7897-465a-a388-3d800ac2de96
# ╠═9ee13ae7-71b0-447b-8ce3-60e9649439a9
# ╟─5ddaf3d2-7dd6-49ac-934d-cb7bdfd37a43
# ╠═aee04655-a30d-430d-a732-b3897938d740
# ╟─f5c7fb00-ee43-44bf-84de-e270f77a0930
# ╟─58760da5-0b94-4c01-9024-a46fe2e7832a
# ╟─6616ca9c-5de2-4615-8707-8b7138b9f2e7
# ╠═46b33f08-05b9-4c26-a9e4-49b0d1dbd0b3
# ╟─db66f0c7-4256-4ffa-9cd3-90766e08381d
# ╠═98b7ed8c-97d8-4966-a65a-27aabf51ec55
# ╠═dbe72760-be53-4f72-a8ef-6b46d7a0b824
# ╠═611a3d4e-18ec-45db-998c-be40e7f59cff
# ╟─8d883827-995a-433d-9143-03b2326c0610
# ╠═161b9ba5-013e-4d7f-86b8-c3aaa53d3206
# ╠═5b8e3935-3dc5-4f7b-b83e-6192f1e285cb
# ╟─aeb4d61b-c091-47f5-8ed5-bddee40960d2
# ╠═1525560b-1ef9-4c9f-b557-0d8f3f37d0cd
# ╠═0f63ae63-d8cd-4b38-99d0-e53d0d68ea11
# ╠═dcee2b99-1b32-4945-bd4b-5ad5c07cd8b5
# ╠═07eccab2-6f71-4eaf-b98a-00bc5eaea693
# ╠═0231b3f4-8e2b-4f68-be5f-553714069082
# ╟─3f05439b-afdb-43dc-9246-7c00cdedffb3
# ╟─48d50aa7-cdf6-48d4-97f5-ed9713d3b325
# ╟─ba707bc7-6b36-4131-97ad-ea757d17d646
