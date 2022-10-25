### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ f914d5dc-9d93-49df-b05d-2a752ef8bc60
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate("..")
	using Revise
	using PlutoUI
	using Test
	using CalciumScoring
end

# ╔═╡ 6788b122-50b1-11ed-3603-c30ea99a0af2
@testset "Agatston" begin
	v1 = ones((4, 2, 2))
	v2 = zeros((4, 2, 2))
	vol = hcat(v1, v2) * 400
	spacing = [0.5, 0.5]
	alg = Agatston()
	answer = 16
	@test score(vol, spacing, alg) == answer
end

# ╔═╡ Cell order:
# ╠═f914d5dc-9d93-49df-b05d-2a752ef8bc60
# ╠═6788b122-50b1-11ed-3603-c30ea99a0af2
