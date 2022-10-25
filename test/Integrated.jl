### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ 1b309749-f59a-484f-a83d-01cc1c57a1cd
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate("..")
	using Revise
	using PlutoUI
	using Test
	using CalciumScoring
end

# ╔═╡ 4f93f1ec-13c3-4bfc-8d96-29021bd060ef
@testset "number" begin
	v1 = [
		-1024 -1024 -1024 -1024
		-1024 800 800 -1024
		-1024 800 800 -1024
		-1024 -1024 -1024 -1024
	]
	vol = cat(v1, v1; dims=3)
	alg = Integrated(vol)
	S_Bkg = -1024
	S_Obj = 800
	answer = 8
	@test score(S_Bkg, S_Obj, alg) == answer
end

# ╔═╡ 5b911c45-a69f-4fdf-b4b5-00fdfc560946
@testset "volume" begin
	v1 = [
		-1024 -1024 -1024 -1024
		-1024 800 800 -1024
		-1024 800 800 -1024
		-1024 -1024 -1024 -1024
	]
	vol = cat(v1, v1; dims=3)
	alg = Integrated(vol)
	S_Bkg = -1024
	S_Obj = 800
	size = [0.5, 0.5, 0.5]
	answer = 8 * 0.5 * 0.5 * 0.5
	@test score(S_Bkg, S_Obj, size, alg) == answer
end

# ╔═╡ c6ff2b47-c371-452c-b764-52ef566ed583
@testset "mass" begin
	v1 = [
		-1024 -1024 -1024 -1024
		-1024 800 800 -1024
		-1024 800 800 -1024
		-1024 -1024 -1024 -1024
	]
	vol = cat(v1, v1; dims=3)
	alg = Integrated(vol)
	S_Bkg = -1024
	S_Obj = 800
	size = [0.5, 0.5, 0.5]
	ρ = 10.0
	answer = 8 * 0.5 * 0.5 * 0.5 * 10
	@test score(S_Bkg, S_Obj, size, ρ, alg) == answer
end

# ╔═╡ Cell order:
# ╠═1b309749-f59a-484f-a83d-01cc1c57a1cd
# ╠═4f93f1ec-13c3-4bfc-8d96-29021bd060ef
# ╠═5b911c45-a69f-4fdf-b4b5-00fdfc560946
# ╠═c6ff2b47-c371-452c-b764-52ef566ed583
