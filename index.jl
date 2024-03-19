### A Pluto.jl notebook ###
# v0.19.40

#> [frontmatter]
#> title = "Home"

using Markdown
using InteractiveUtils

# ╔═╡ 1f656df7-f84c-42f7-8b55-1b415474ce5b
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(joinpath(pwd(), "docs"))
	Pkg.instantiate()

	using HTMLStrings: to_html, head, link, script, divv, h1, img, p, span, a, figure, hr
	using PlutoUI
end

# ╔═╡ 4c391e35-5623-4903-89b4-17e48ae4d7af
md"""
# Quickstart
"""

# ╔═╡ cef64d20-56c9-4478-bcdc-3a965884e4ad
md"""
## TLDR
"""

# ╔═╡ f0b3cbeb-9132-441b-8ff2-090d0502804c
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

# ╔═╡ b4f6e7e3-bfbf-46ca-8ea4-4cc80fad07c3
md"""
!!! success ""
	This package contains four state-of-the-art coronary artery calcium (CAC) quantification methods.
	- Agatston scoring [1, 2]
	- Volume fraction calcium mass [3]
	- Integrated calcium mass [3, 4]
	- Spatially weighted calcium scoring [5, 6]
	
	This collection of notebooks will briefly showcase how to use CalciumScoring.jl to quantify calcium using these techniques.
"""

# ╔═╡ 24d80989-0109-45b7-8d59-188b29186037
md"""
The general approach is always the same. First, you must prepare the region of interest within the CT scan containing potential calcium. Then, you must gather important scan information like voxel size, spacing information, calibration intensity, etc. Lastly, you calculate the calcium using `score()` with the appropriate algorithm and arguments.
"""

# ╔═╡ 93926ede-e23d-4f2e-b485-daae86ed5cbf
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

# ╔═╡ d0befdf6-462a-4058-85c6-78f92995ce90
md"""
This is an example of the process for Agatston scoring:
"""

# ╔═╡ 1324b7d2-2dad-4313-bdf5-2c349bdd4a91
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

# ╔═╡ a92cbf7e-0b09-457c-a461-2c649e8395e6
md"""
## References
"""

# ╔═╡ 7ab48343-2827-4f32-ab5b-a086981e62a1
md"""
1. [Quantification of coronary artery calcium using ultrafast computed tomography](https://doi.org/10.1016/0735-1097(90)90282-t)
2. [Ultra-low-dose coronary artery calcium scoring using novel scoring thresholds for low tube voltage protocols—a pilot study ](https://doi.org/10.1093/ehjci/jey019)
3. [Coronary artery calcium mass measurement based on integrated intensity and volume fraction techniques](https://doi.org/10.1117/1.JMI.10.4.043502)
4. [Integrated intensity-based technique for coronary artery calcium mass measurement: A phantom study](https://doi.org/10.1002/mp.16326)
5. [An alternative method for quantifying coronary artery calcification: the multi-ethnic study of atherosclerosis (MESA)](https://doi.org/10.1186/1471-2342-12-14)
6. [Spatially Weighted Coronary Artery Calcium Score and Coronary Heart Disease Events in the Multi-Ethnic Study of Atherosclerosis](https://doi.org/10.1161/CIRCIMAGING.120.011981)

"""

# ╔═╡ ab941a39-c28c-46dd-ab69-9d244216f3d1
md"""
# Next Steps

The core function, `score()`, implements all four types of CAC quantification methods. Agatston scoring is the traditional method, but recent advances in the field of medical physics have introduced better approaches to coronary artery calcium mass quantification. We will explore these as we go. 

Check out the tutorials sections to see how to implement the traditional Agatston approach or more recent improved approaches, like the Volume Fraction method.
"""

# ╔═╡ 4cf29548-8ecd-40ff-a573-3d3475ccf382
md"""
## Tutorials
"""

# ╔═╡ e01cc6df-da72-44c4-904c-2e2740df5a00
md"""
## API
"""

# ╔═╡ a2ada594-28d5-4981-87bf-0c369b260f57
to_html(hr())

# ╔═╡ d02cda6f-7ff0-4fc5-a5bb-5af4b2043c33
TableOfContents()

# ╔═╡ a1af303f-6214-4763-bd17-d989e5c0ab2d
data_theme = "synthwave";

# ╔═╡ 5c69ba9d-1b4d-46b5-872d-ada6035624ab
function index_title_card(title::String, subtitle::String, image_url::String; data_theme::String = "pastel", border_color::String = "primary")
	return to_html(
	    divv(
	        head(
				link(:href => "https://cdn.jsdelivr.net/npm/daisyui@3.7.4/dist/full.css", :rel => "stylesheet", :type => "text/css"),
	            script(:src => "https://cdn.tailwindcss.com")
	        ),
			divv(:data_theme => "$data_theme", :class => "card card-bordered flex justify-center items-center border-$border_color text-center w-full dark:text-[#e6e6e6]",
				divv(:class => "card-body flex flex-col justify-center items-center",
					img(:src => "$image_url", :class => "h-24 w-24 md:h-40 md:w-40 rounded-md", :alt => "$title Logo"),
					divv(:class => "text-3xl md:text-5xl font-bold bg-gradient-to-r from-accent to-primary inline-block text-transparent bg-clip-text py-10", "$title"),
					p(:class => "card-text text-md font-serif", "$subtitle"
					)
				)
			)
	    )
	)
end;

# ╔═╡ 11d6cc2b-f30b-4d4f-8346-738650489b21
index_title_card(
	"CalciumScoring.jl",
	"Comprehensive Suite of Coronary Artery Calcium Quantification Algorithms",
	"https://img.freepik.com/free-vector/ct-scan-concept-illustration_114360-7073.jpg";
	data_theme = data_theme
)

# ╔═╡ 0c7647e5-033a-4b11-9305-bcd274b96bc7
begin
	struct Article
		title::String
		path::String
		image_url::String
	end

	function article_card(article::Article, color::String; data_theme = "pastel")
		a(:href => article.path, :class => "w-full md:w-1/2 p-2",
			divv(:data_theme => "$data_theme", :class => "card card-bordered border-$color text-center dark:text-[#e6e6e6]",
				divv(:class => "card-body justify-center items-center h-32 md:h-40",
					p(:class => "text-lg md:text-2xl", article.title),
					p(:class => "text-sm md:text-base", "Click to open the notebook")
				),
				figure(
					img(:class => "w-full h-40 md:h-48 object-cover", :src => article.image_url, :alt => article.title)
				)
			)
		)
	end

	article_list_tutorials = Article[
		Article("Agatston Scoring", "docs/01_agatston_scoring.jl", "https://img.freepik.com/free-vector/human-heart-disease-symbol_1308-107392.jpg"),
		Article("Volume Fraction", "docs/02_volume_fraction.jl", "https://img.freepik.com/free-vector/heart-with-heartbeat_78370-2027.jpg"),
		Article("Integrated", "docs/03_integrated.jl", "https://img.freepik.com/free-vector/high-blood-pressure-abstract-concept-illustration_335657-4603.jpg"),
		Article("Spatially Weighted", "docs/04_spatially_weighted.jl", "https://img.freepik.com/free-vector/heartbeat-with-heart-rate-graph_1308-108566.jpg"),
	]

	article_list_api = Article[
		Article("API", "docs/99_api.jl", "https://img.freepik.com/free-photo/modern-technology-workshop-creativity-innovation-communication-development-generated-by-ai_188544-24548.jpg"),
	]
end;

# ╔═╡ 1ba43e44-7839-46e1-bb8e-cb36b5cbe3a1
to_html(
    divv(:class => "flex flex-wrap justify-center items-start",
        [article_card(article, "secondary"; data_theme = data_theme) for article in article_list_tutorials]...
    )
)

# ╔═╡ e27e3573-f422-4dc0-bc38-21720dcf8f89
to_html(
    divv(:class => "flex flex-wrap justify-center items-start",
        [article_card(article, "secondary"; data_theme = data_theme) for article in article_list_api]...
    )
)

# ╔═╡ Cell order:
# ╟─11d6cc2b-f30b-4d4f-8346-738650489b21
# ╟─4c391e35-5623-4903-89b4-17e48ae4d7af
# ╟─cef64d20-56c9-4478-bcdc-3a965884e4ad
# ╟─f0b3cbeb-9132-441b-8ff2-090d0502804c
# ╟─b4f6e7e3-bfbf-46ca-8ea4-4cc80fad07c3
# ╟─24d80989-0109-45b7-8d59-188b29186037
# ╟─93926ede-e23d-4f2e-b485-daae86ed5cbf
# ╟─d0befdf6-462a-4058-85c6-78f92995ce90
# ╟─1324b7d2-2dad-4313-bdf5-2c349bdd4a91
# ╟─a92cbf7e-0b09-457c-a461-2c649e8395e6
# ╟─7ab48343-2827-4f32-ab5b-a086981e62a1
# ╟─ab941a39-c28c-46dd-ab69-9d244216f3d1
# ╟─4cf29548-8ecd-40ff-a573-3d3475ccf382
# ╟─1ba43e44-7839-46e1-bb8e-cb36b5cbe3a1
# ╟─e01cc6df-da72-44c4-904c-2e2740df5a00
# ╟─e27e3573-f422-4dc0-bc38-21720dcf8f89
# ╟─a2ada594-28d5-4981-87bf-0c369b260f57
# ╟─d02cda6f-7ff0-4fc5-a5bb-5af4b2043c33
# ╟─1f656df7-f84c-42f7-8b55-1b415474ce5b
# ╟─a1af303f-6214-4763-bd17-d989e5c0ab2d
# ╟─5c69ba9d-1b4d-46b5-872d-ada6035624ab
# ╟─0c7647e5-033a-4b11-9305-bcd274b96bc7
