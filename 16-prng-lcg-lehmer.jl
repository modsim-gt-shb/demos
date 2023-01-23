### A Pluto.jl notebook ###
# v0.18.2

using Markdown
using InteractiveUtils

# ╔═╡ 17c99eca-f8ac-4ba4-88ce-bdb5bc0e0754
md"""
# Linear congruential generators #

This notebook accompanies the slides on pseudorandom number generation.
"""

# ╔═╡ 6569c6c1-aaba-4378-9f2d-0b3f3bb4641d
"""
    g(x; m, a, c)

Returns the value of the discrete map corresponding to a linear congruential generator, ``x_{i+1} = (a x_i + c) \\mod m``.
"""
function g(x; m, a, c)
	(a*x + c) % m
end

# ╔═╡ 028b8696-20d1-44cb-be37-348c21925a72
"""
    generate_lcg(x0; kwargs...)

Generates a sequence of values ``x_0, x_1, \\ldots, x_p`` such that

- ``x_0 = `` `x0`
- ``x_{i+1} = g(x_i; m, a, c)``, that is, the sequence of values produced by the linear congruential generator ``g(x) = (ax + c) \\mod m``
- ``p =`` the period of the sequence, that is, ``x_0 = x_p``.

This function returns the tuple ``(p, x)`` where ``x`` is the sequence (including both ``x_0`` and ``x_p``).

> The `kwargs...` are the same as those of the LCG discrete-map function, `g`.
"""
function generate_lcg(x0; kwargs...)
	X = [x0]
	while !(last(X) in X[1:length(X)-1])
		push!(X, g(last(X); kwargs...))
	end
	(length(X)-1, X)
end

# ╔═╡ e70aff7a-f416-4425-8f73-14bab4b3cd69
"""
    generate_lehmer(x0; kwargs...)

Generates a Lehmer sequence, i.e., the linear congruential sequence produced by `generate_lcg` with `c=0`.

> The `kwargs...` are the same as those of the LCG discrete-map function, `g`.
"""
function generate_lehmer(x0; kwargs...)
	generate_lcg(x0; c=0, kwargs...)
end

# ╔═╡ bb2b2eba-1353-4fbe-8a64-b425c36486cb
# Demo:
generate_lehmer(2; m=2^6, a=13) # try a=2

# ╔═╡ 156f0a1c-ce21-4749-bfae-f48767fcf1aa
md"""
# All sequences of a Lehmer generator #
"""

# ╔═╡ 6d9b3d10-c891-49d1-a181-70539af41ac5
md"""
Let's pick some parameters for a Lehmer generator. If you are playing with this notebook, try these settings:
- `m=2^6`, `a=13`
- `m=2^6`, `a=2`
- `m=2^5-1`, `a=13` [_Note:_ `2^5-1=31` is a [Mersenne prime](https://en.wikipedia.org/wiki/Mersenne_prime).]
- `m=2^5-1`, `a=2`
"""

# ╔═╡ 1d246900-b43f-11ec-3f93-d3eb44bf8338
begin
	m = 2^6 # also try: 2^5-1, which is a Mersenne Prime
	a = 13
	m, a
end

# ╔═╡ fd62ed94-7b3f-473a-b343-f192898502fe
"""
    generate_all_subsequences(m, a)

Given a Lehmer generator with modulus `m` and multiplier `a`, this function calculates and returns all of its subsequences.

More specifically, it returns a pair (`n`, `χ`) where
- `n` is the number of subsequences
- Let `1 <= i <= n` denote the `i`-th subsequence. Then a value `k` is a part of that subsequence if `χ[k] == i`.
"""
function generate_all_subsequences(m, a)
	χ = repeat([0], m-1)
	num_subsequences = 0
	for k in 1:m-1
		if χ[k] == 0
			num_subsequences += 1
			_, X_k = generate_lehmer(k; m=m, a=a)
			χ[X_k] .= num_subsequences
		end
	end
	(num_subsequences, χ)
end

# ╔═╡ bd536a3b-7848-49f2-8d9d-09e1ad387709
begin
	num_subsequences, χ = generate_all_subsequences(m, a)
	for i in 1:num_subsequences
		p_i, χ_i = generate_lehmer(findall(χ .== i)[1]; m=m, a=a)
		println(i, " (p=", p_i, "): ", χ_i)
	end
end

# ╔═╡ Cell order:
# ╟─17c99eca-f8ac-4ba4-88ce-bdb5bc0e0754
# ╠═6569c6c1-aaba-4378-9f2d-0b3f3bb4641d
# ╟─028b8696-20d1-44cb-be37-348c21925a72
# ╟─e70aff7a-f416-4425-8f73-14bab4b3cd69
# ╠═bb2b2eba-1353-4fbe-8a64-b425c36486cb
# ╠═156f0a1c-ce21-4749-bfae-f48767fcf1aa
# ╟─6d9b3d10-c891-49d1-a181-70539af41ac5
# ╠═1d246900-b43f-11ec-3f93-d3eb44bf8338
# ╠═fd62ed94-7b3f-473a-b343-f192898502fe
# ╠═bd536a3b-7848-49f2-8d9d-09e1ad387709
