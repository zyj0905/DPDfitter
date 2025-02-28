### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 19487fd8-f5aa-11ef-20ca-f99cfee6a036
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/JuliaHEP/LorentzVectorBase.jl"),
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/FourVectors.jl"),
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/DecayAngles.jl"),
		Pkg.PackageSpec("Parameters"),
		Pkg.PackageSpec("DataFrames")])
	# 
	using DecayAngles
	using FourVectors
	using DataFrames
	using Parameters
end

# ╔═╡ fb89cd09-4921-471c-981a-f0702224f22e
md"""
# Kimenatics of Dˣ D π

In this notebook, we take a single kinematic point,

```
Dst: 2.0085299,0.0570074,-0.026685,0.0479813
D: 1.8967276,-0.108349,-0.056907,-0.295222
Pi: 0.3153471, 0.1034199, 0.0873037,0.2482869
```

and compute all helicity angles

"""

# ╔═╡ b40d123e-dedc-445b-95a5-ed043541e0d8
function pure_B(p::FourVector, p_ref::FourVector)
	@unpack cosθ, ϕ = spherical_coordinates(p_ref)
	θ = acos(cosθ)
	γ = boost_gamma(p_ref)
	p |> Rz(-ϕ) |> Ry(-θ) |> Bz(-γ) |> Ry(θ) |> Rz(ϕ)
end

# ╔═╡ ce9d24dc-eaac-4403-8d35-4db79c0f503e
function pure_B(system::NamedTuple)
	ptot = collect(system) |> sum
	map(system) do p
		pure_B(p, ptot)
	end |> NamedTuple{keys(system)}
end

# ╔═╡ 361f84d9-0093-4580-bf6d-cd9d8f8cfb90
function helicity_angles(four_vectors_rf, topology)
	momenta_dict = Dict(pairs(four_vectors_rf))
	tree_empty = DecayNode(topology);
	tree_with_particle_order = add_indices_order(tree_empty);
	tree_with_four_vectors = add_transform_through(
		HelicityTransformation, tree_with_particle_order, momenta_dict);
	# 
	decay_angles(tree_with_four_vectors)
end

# ╔═╡ 3e46401a-c950-49d6-8717-1255ced02ef7
four_vectors_nt = (
	Dst = FourVector(0.0570074,-0.026685,0.0479813; E=2.0085299),
	D = FourVector(-0.108349,-0.056907,-0.295222; E=1.8967276),
	Pi = FourVector(0.1034199, 0.0873037,0.2482869; E=0.3153471)
)

# ╔═╡ e01f7bd5-2089-4d0c-a0b4-93c6917ed7ae
## test if masses are reasonable
four_vectors_nt |> collect .|> mass

# ╔═╡ 91154f6e-4469-4de2-bae9-a7bcd867be38
# sqrt(s)
collect(four_vectors_nt) |> sum |> mass

# ╔═╡ ac528202-f867-41f3-8aca-3ca4993f2987
# pure boost to the rest frame
four_vectors_rf = pure_B(four_vectors_nt)

# ╔═╡ a84fbe38-6d26-4738-a9bc-aa479d32f22a
begin
	vcat(
		helicity_angles(four_vectors_rf, ((:Pi, :D), :Dst)),
		helicity_angles(four_vectors_rf, ((:D, :Dst), :Pi)),
		helicity_angles(four_vectors_rf, ((:Dst, :Pi), :D))
	) |> DataFrame
end

# ╔═╡ Cell order:
# ╟─fb89cd09-4921-471c-981a-f0702224f22e
# ╠═19487fd8-f5aa-11ef-20ca-f99cfee6a036
# ╠═b40d123e-dedc-445b-95a5-ed043541e0d8
# ╠═ce9d24dc-eaac-4403-8d35-4db79c0f503e
# ╠═361f84d9-0093-4580-bf6d-cd9d8f8cfb90
# ╠═3e46401a-c950-49d6-8717-1255ced02ef7
# ╠═e01f7bd5-2089-4d0c-a0b4-93c6917ed7ae
# ╠═91154f6e-4469-4de2-bae9-a7bcd867be38
# ╠═ac528202-f867-41f3-8aca-3ca4993f2987
# ╠═a84fbe38-6d26-4738-a9bc-aa479d32f22a
