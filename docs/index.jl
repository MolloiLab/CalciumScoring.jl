### A Pluto.jl notebook ###
# v0.19.26

#> [frontmatter]
#> title = "docs/CalciumScoring.jl"
#> sidebar = "false"

using Markdown
using InteractiveUtils

# ╔═╡ 1f656df7-f84c-42f7-8b55-1b415474ce5b
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")
	Pkg.instantiate()

	using HTMLStrings: to_html, head, link, script, divv, h1, img, p, span, a, figure, hr
	using PlutoUI
end

# ╔═╡ e099ee7c-b797-4aa0-8a2a-209415a72c09
md"""
## QuickStart
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
					img(:src => "$image_url", :class => "h-24 w-24 md:h-52 md:w-52 rounded-md", :alt => "$title Logo"),
					divv(:class => "text-5xl font-bold bg-gradient-to-r from-accent to-primary inline-block text-transparent bg-clip-text py-10", "$title"),
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

# ╔═╡ 1c9a2ea8-8037-4592-8506-83659b1c92be
struct Article
	title::String
	path::String
	image_url::String
end

# ╔═╡ ab3e7613-3926-4216-a2ea-8bef708a9e20
article_list_quickstart = Article[
	Article("Getting Started", "01_getting_started.jl", "https://img.freepik.com/free-photo/futuristic-spaceship-takes-off-into-purple-galaxy-fueled-by-innovation-generated-by-ai_24640-100023.jpg"),
];

# ╔═╡ 201962ba-6c3d-43bb-8593-56a46fe1e9a4
article_list_tutorials = Article[
	Article("Agatston Scoring", "02_agatston_scoring.jl", "https://img.freepik.com/free-vector/human-heart-disease-symbol_1308-107392.jpg"),
	Article("Volume Fraction", "03_volume_fraction.jl", "https://img.freepik.com/free-vector/heart-with-heartbeat_78370-2027.jpg"),
	Article("Integrated", "04_integrated.jl", "https://img.freepik.com/free-vector/high-blood-pressure-abstract-concept-illustration_335657-4603.jpg"),
	Article("Spatially Weighted", "05_spatially_weighted.jl", "https://img.freepik.com/free-vector/heartbeat-with-heart-rate-graph_1308-108566.jpg"),
];

# ╔═╡ 1f7202ca-5612-487e-9fd8-718fe3245fab
article_list_api = Article[
	Article("API", "06_api.jl", "https://img.freepik.com/free-photo/modern-technology-workshop-creativity-innovation-communication-development-generated-by-ai_188544-24548.jpg"),
];

# ╔═╡ 9fe7717e-0d5a-4238-8834-ea7608b27790
function article_card(article::Article, color::String; data_theme = "pastel")
    a(:href => article.path, :class => "w-1/2 p-2",
		divv(:data_theme => "$data_theme", :class => "card card-bordered border-$color text-center dark:text-[#e6e6e6]",
			divv(:class => "card-body justify-center items-center",
				p(:class => "card-title", article.title),
				p("Click to open the notebook")
			),
			figure(
				img(:class =>"w-full", :src => article.image_url, :alt => article.title)
			)
        )
    )
end;

# ╔═╡ 36db6fb8-e010-44c7-9d68-5da260b0da51
to_html(
    divv(:class => "flex flex-wrap justify-center items-start",
        [article_card(article, "accent"; data_theme = data_theme) for article in article_list_quickstart]...
    )
)

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
# ╟─e099ee7c-b797-4aa0-8a2a-209415a72c09
# ╟─36db6fb8-e010-44c7-9d68-5da260b0da51
# ╟─4cf29548-8ecd-40ff-a573-3d3475ccf382
# ╟─1ba43e44-7839-46e1-bb8e-cb36b5cbe3a1
# ╟─e01cc6df-da72-44c4-904c-2e2740df5a00
# ╟─e27e3573-f422-4dc0-bc38-21720dcf8f89
# ╟─a2ada594-28d5-4981-87bf-0c369b260f57
# ╟─d02cda6f-7ff0-4fc5-a5bb-5af4b2043c33
# ╟─1f656df7-f84c-42f7-8b55-1b415474ce5b
# ╟─a1af303f-6214-4763-bd17-d989e5c0ab2d
# ╟─5c69ba9d-1b4d-46b5-872d-ada6035624ab
# ╟─1c9a2ea8-8037-4592-8506-83659b1c92be
# ╟─ab3e7613-3926-4216-a2ea-8bef708a9e20
# ╟─201962ba-6c3d-43bb-8593-56a46fe1e9a4
# ╟─1f7202ca-5612-487e-9fd8-718fe3245fab
# ╟─9fe7717e-0d5a-4238-8834-ea7608b27790
