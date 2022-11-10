### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 8409fd11-e215-4982-a183-0d1621a8264d
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate("..")
	using Revise
	using PlutoUI
	using Test
	using CalciumScoring
end

# ╔═╡ 72565940-19f0-4e23-aeff-4bf47b15f03f
@testset "mask_elements" begin
	img = [
		1 10 1 1
		-2 -2 -2 -2
		1 1 1 1
		1 1 -2 1
		1 1 1 1
		-2 -2 -2 -2
	]
	answer = Bool.([
		1 1 1 1
		0 0 0 0
		1 1 1 1
		1 1 1 1
		1 1 1 1
		0 0 0 0
	])
	@test mask_elements(img, 0, 2) == answer
end

# ╔═╡ 5b2a38ac-6124-11ed-1ccc-2194d81ab664
@testset "mask_elements" begin
	img = [
		1 -2 -2 -2
		-2 -2 -2 -2
		1 -2 -2 -2
		-2 -2 -2 -2
		1 -2 -2 -2
		-2 -2 -2 -2
	]
	answer = Bool.([
		1 0 0 0
		0 0 0 0
		1 0 0 0
		0 0 0 0
		1 0 0 0
		0 0 0 0
	])
	@test mask_elements(img, 0, 2) == answer
end

# ╔═╡ Cell order:
# ╠═8409fd11-e215-4982-a183-0d1621a8264d
# ╠═72565940-19f0-4e23-aeff-4bf47b15f03f
# ╠═5b2a38ac-6124-11ed-1ccc-2194d81ab664
