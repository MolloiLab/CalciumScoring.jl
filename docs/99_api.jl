### A Pluto.jl notebook ###
# v0.19.26

#> [frontmatter]
#> title = "API"
#> category = "API"

using Markdown
using InteractiveUtils

# ╔═╡ 4430a6d2-6e99-416a-8afa-27c3c049dc08
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")
	Pkg.instantiate()

	using PlutoUI, CalciumScoring
end;

# ╔═╡ 4ba7105a-4d49-447c-b053-3c265ec38139
md"""
# Overview
"""

# ╔═╡ 70da4548-6a59-4b6b-8885-0dbab4dbffb4
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

# ╔═╡ fbc3ffb1-414f-472b-a04d-f98d43d49a1c
md"""
# References
"""

# ╔═╡ ce457dad-300c-430b-9fb3-b2f0048d99f3
TableOfContents()

# ╔═╡ 3e0d2872-7d4d-449b-816e-85859b7159d2
all_names = [name for name in names(CalciumScoring)];

# ╔═╡ 5f465635-4f44-4509-80ca-0742b65bd810
exported_functions = filter(x -> x != :CalciumScoring, all_names);

# ╔═╡ d69659cc-0cd9-40f8-8902-55ed87c0d322
function generate_docs(exported_functions)
    PlutoUI.combine() do Child
        md"""
        $([md" $(Docs.doc(eval(name)))" for name in exported_functions])
        """
    end
end;

# ╔═╡ 26d2e031-7a40-424b-8196-afb675cdf157
generate_docs(exported_functions)

# ╔═╡ Cell order:
# ╟─4ba7105a-4d49-447c-b053-3c265ec38139
# ╟─70da4548-6a59-4b6b-8885-0dbab4dbffb4
# ╟─fbc3ffb1-414f-472b-a04d-f98d43d49a1c
# ╟─26d2e031-7a40-424b-8196-afb675cdf157
# ╟─4430a6d2-6e99-416a-8afa-27c3c049dc08
# ╟─ce457dad-300c-430b-9fb3-b2f0048d99f3
# ╟─3e0d2872-7d4d-449b-816e-85859b7159d2
# ╟─5f465635-4f44-4509-80ca-0742b65bd810
# ╟─d69659cc-0cd9-40f8-8902-55ed87c0d322
