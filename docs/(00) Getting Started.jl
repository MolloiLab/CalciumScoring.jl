### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 19b78bb6-ba8c-4e43-85c2-ba9148535f4b
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(temp = true)
	Pkg.add("PlutoUI")
	
	using PlutoUI
end

# ╔═╡ e0ee3cd9-f889-48de-99d8-dd1ac0ad4dc2
md"""
## Overview
"""

# ╔═╡ 0dd245c9-eb79-4c95-b649-108fa372fe37
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
            class B,C,D,E unique2;
    </div>
</body>
</html>
"""

# ╔═╡ 41898031-4baa-47a8-95ae-bb42b5483193
md"""
!!! success "Overview"
	This package contains four state-of-the-art coronary artery calcium (CAC) quantification methods.
	- Agatston scoring [1, 2]
	- Volume fraction calcium mass [3]
	- Integrated calcium mass [3, 4]
	- Spatially weighted calcium scoring [5, 6]
	
	This collection of notebooks will briefly showcase how to use CalciumScoring.jl to quantify calcium using these techniques.
	
	**References**
	1. [Quantification of coronary artery calcium using ultrafast computed tomography](https://doi.org/10.1016/0735-1097(90)90282-t)
	2. [Ultra-low-dose coronary artery calcium scoring using novel scoring thresholds for low tube voltage protocols—a pilot study ](https://doi.org/10.1093/ehjci/jey019)
	3. [Coronary artery calcium mass measurement based on integrated intensity and volume fraction techniques](https://doi.org/10.1117/1.JMI.10.4.043502)
	4. [Integrated intensity-based technique for coronary artery calcium mass measurement: A phantom study](https://doi.org/10.1002/mp.16326)
	5. [An alternative method for quantifying coronary artery calcification: the multi-ethnic study of atherosclerosis (MESA)](https://doi.org/10.1186/1471-2342-12-14)
	6. [Spatially Weighted Coronary Artery Calcium Score and Coronary Heart Disease Events in the Multi-Ethnic Study of Atherosclerosis](https://doi.org/10.1161/CIRCIMAGING.120.011981)

"""

# ╔═╡ cf35a8b5-a427-4c7d-9a54-32c4243995f4
md"""
## TLDR
"""

# ╔═╡ 6301aaa6-de72-4c12-9aa8-9ededce85ace
md"""
The general approach is always the same. First, you must prepare the region of interest within the CT scan containing potential calcium. Then, you must gather important scan information like voxel size, spacing information, calibration intensity, etc. Lastly, you calculate the calcium using `score()` with the appropriate algorithm and arguments.
"""

# ╔═╡ a2676713-b1c6-4b6a-9803-8601bf889e29
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
            A["1. Prepare Region of Interest"] --> B["2. Gather CT Scan Info"]
            B --> C["3. Calculate Calcium"]
    </div>
</body>
</html>
"""

# ╔═╡ b89820d7-1f4d-4a67-9df8-51a7ce1657d7
md"""
This is an example of the process for Agatston scoring:
"""

# ╔═╡ 7cef5506-04a6-45cb-a685-d21888b28b76
md"""
```julia
# 1
region_of_interest = ct_scan .* mask

# 2
spacing = [1, 1, 1]
mass_cal_factor = density_calibration / mean(intensity_calibration)

# 3
agatson_score, volume_score, mass_score = score(region_of_interest, spacing, mass_cal_factor, Agatston())
```
"""

# ╔═╡ 17a5a766-2795-4b4b-b553-84c0d3d1ca6e
md"""
## Next Steps

The core function, `score()`, implements all four types of CAC quantification methods. Agatston scoring is the traditional method, but recent advances in the field of medical physics have introduced better approaches to coronary artery calcium mass quantification. We will explore these as we go. 

Check out the [Agatston scoring](docs/(01) Agatston Scoring.jl) tutorial to see how to implement the traditional approach.
"""

# ╔═╡ 92d11a7d-0cc9-42bb-80d9-226b61a3949b
TableOfContents()

# ╔═╡ Cell order:
# ╟─e0ee3cd9-f889-48de-99d8-dd1ac0ad4dc2
# ╟─0dd245c9-eb79-4c95-b649-108fa372fe37
# ╟─41898031-4baa-47a8-95ae-bb42b5483193
# ╟─cf35a8b5-a427-4c7d-9a54-32c4243995f4
# ╟─6301aaa6-de72-4c12-9aa8-9ededce85ace
# ╟─a2676713-b1c6-4b6a-9803-8601bf889e29
# ╟─b89820d7-1f4d-4a67-9df8-51a7ce1657d7
# ╟─7cef5506-04a6-45cb-a685-d21888b28b76
# ╟─17a5a766-2795-4b4b-b553-84c0d3d1ca6e
# ╠═19b78bb6-ba8c-4e43-85c2-ba9148535f4b
# ╠═92d11a7d-0cc9-42bb-80d9-226b61a3949b
