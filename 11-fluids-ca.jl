### A Pluto.jl notebook ###
# v0.17.5

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

# ╔═╡ f80f5adc-97e5-11ec-341c-337dc6192b37
begin
	using LaTeXStrings
	using Plots
	using OffsetArrays
	using StatsBase
	using Random
	using PlutoUI
end

# ╔═╡ a5dbdc32-d7a6-4e2b-8c47-0de9b90d7235
md"""
# A "fluidic" cellular automaton #

This notebook implements the "FHP" cellular automaton model of a fluid discussed in the reading [G20, Chapter 10].
"""

# ╔═╡ 915ee6cf-ff20-48b0-8803-1ff77f8bde93
md"""
## Geometry ##

The domain is tiled into equilateral triangles whose vertices are "sites." Each of the ``2M \times 2M`` sites is assigned a _logical_ coordinate pair ``(i, j)`` where ``-M \leq i, j < M``.

Each site also has _physical_ coordinates. Suppose the edge-length of each triangle is ``s``. Then each triangle's "height" is ``r = s\sqrt{\frac{3}{4}}.`` _(The height is the distance from any vertex to the center of the edge opposite from it.)_ Then we will assign each site ``(i, j)`` to the physical coordinate position ``(x_{i,j}, y_{i,j})``:

``x_{i,j} = \left\{ \begin{array}{ll} s\cdot i & \mbox{if } j \mbox{ is even} \\ s \cdot \left(i + \frac{1}{2}\right) & \mbox{otherwise} \end{array} \right.``

and

``y_{i,j} = r \cdot j``.

Words, words, words. 😉. Here is a picture, for ``M=3`` and ``s=1``:
"""

# ╔═╡ 52d94124-6cc3-4af1-b010-63c668fd0378
STEPSIZE = 1

# ╔═╡ fb6e0d76-7c3f-4f1d-a4f0-a78f39b8b7ca
md"""
As also shown above, each site is centered within the hexagon formed from the convex hull of incident triangles. Each site's six nearest neighbors constitute its neighborhood.
"""

# ╔═╡ 473ca08e-eff6-48a5-8981-afad0a9186d5
md"""
## Simulator "architecture" ##

The function `simulate`, below, implements the simulation framework for our CA model.
"""

# ╔═╡ 9b1088a7-adbf-47da-9c78-032b55e394e8
md"""
## States and update rules ##

The _state_ ``x`` of any site is a six-bit vector, where each bit corresponds to the edge that connects that site to one of its neighbors. We denote each bit by an integer ``\alpha`` numbered from 0 to 5.
"""

# ╔═╡ aa90a3c7-6c7a-4f2a-a08f-36cc7ede1b30
md"""
If a bit is set, then there is a particle at the site on that edge. At the start of each simulation time step, our convention will be that the particle is traveling toward the site.
"""

# ╔═╡ 6392b141-0548-48ac-91c8-b5589279aff3
md"""
The state is updated according to one of the following three cases:
1. Two-body collisions
2. Three-body collisions
3. Passthrough (all other cases)
"""

# ╔═╡ 34c88f6b-97e3-4e9b-9301-d44269b0ded7
md"""
### Case 1: Two-body collisions ###

This case occurs when there are _exactly_ two particles at the site and they collide "head on." That means they live on opposite parallel edges of the site.
"""

# ╔═╡ 3c4da1eb-15fe-4565-a3c1-cf5db155efe6
"""
    update_2body!(Step, Sites, i, j)

If a two-body collision exists at `Sites[i, j, :]`, sets the new state at `Step[i, j, :]`, which represents a random rotation of the bodies by one unit, and returns `true`. Otherwise, returns `false` and leaves `Step[i, j, :]` unmodified.
"""
function update_2body!(Step, Site, i, j)
	active = findall(Site[i, j, :])
	if diff(active) == [3]
		dir = (rand() < 0.5) ? 1 : -1
		Step[i, j, active] .= 0
		Step[i, j, (active .+ 6 .+ dir) .% 6] .= 1
		return true
	else
		return false
	end
end

# ╔═╡ e9c572e9-63ee-4e71-88a6-a16f2d079fdd
md"""
### Case 2: Three-body collisions ###

This case occurs when there are _exactly_ three particles at the site and they collide "symmetrically," meaning they occur on every other edge.
"""

# ╔═╡ 4e699d98-2786-4aa5-90e9-b8a094ea7263
"""
    update_3body!(Step, Sites, i, j)

If a three-body collision exists at `Sites[i, j, :]`, sets the new state at `Step[i, j, :]`, which is just a copy of the existing state. Otherwise, returns `false` and leaves `Step[i, j, :]` unmodified.
"""
function update_3body!(Step, Site, i, j)
	###
	### YOUR CODE HERE ###
	###
	### BEGIN SOLUTION
	active = findall(Site[i, j, :])
	if diff(active) == [2, 2]
		Step[i, j, :] = Site[i, j, :]
		return true
	else
		return false
	end
	### END SOLUTION
end

# ╔═╡ bf939585-832b-463c-8752-4ab79bbcffab
md"""
### Case 3: Passthrough ###

For all remaining cases, the particles should "pass-through." That means they move to the opposite edge.
"""

# ╔═╡ 0e9d40e8-a5e3-4c3a-9f8e-b485be8554e7
"""
    update_passthrough!(Step, Sites, i, j)

Moves each particle at `Sites[i, j, a]` to its opposite edge in `Step[i, j, b]`. Always returns `true`.
"""
function update_passthrough!(Step, Site, i, j)
	###
	### YOUR CODE HERE ###
	###
	### BEGIN SOLUTION
	active = findall(Site[i, j, :])
	Step[i, j, :] .= 0
	Step[i, j, (active .+ 3) .% 6] .= 1
	true
	### END SOLUTION
end

# ╔═╡ ce6284b2-4fef-4ad5-998e-449db4a7a943
md"""
### Transport ###

The last operation is to move all outgoing particles to their neighboring sites.

To help you out, note the following helper function, which determines the site and edge that "mirrors" a given site and edge.
"""

# ╔═╡ 0e90eb77-d635-4190-abfd-b2ad2198a520
md"""# A simulation #"""

# ╔═╡ 4d668f04-0b5e-4e0f-9d9a-dabef42d855a
begin
	M = 10
	T = 0:100
	
	L"M=%$(M),\ s=%$(STEPSIZE) \quad \Longrightarrow \quad r\ \simeq\ s\sqrt{\frac{3}{4}}\ \simeq\ %$(round(STEPSIZE*sqrt(3/4), sigdigits=4))"
end

# ╔═╡ c7bff650-7801-4d88-9a27-9153b7eca61f
md"**Pick a time from the range $(T).**"

# ╔═╡ 4075dd9c-88b8-46fb-99eb-d3e03e151a05
@bind t_example Slider(T)

# ╔═╡ dc9c9d76-5967-43a2-958a-9c61ee9101bc
md"# Helper functions #"

# ╔═╡ df529624-20c5-4252-897c-aaaed4db673e
function site_apply(Sites, fun)
	mapslices(fun, Sites, dims=[3])[:, :, 0]
end

# ╔═╡ 8a7b06e0-5cdb-49b8-9658-c94db5399c6f
function αs(i, j)
	if (j % 2) == 0 # even
		OffsetArray([(1, 0), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1)], 0:5)
	else
		OffsetArray([(1, 0), (1, 1), (0, 1), (-1, 0), (0, -1), (1, -1)], 0:5)
	end
end

# ╔═╡ 4b7b3eeb-8fc7-4cb3-9b33-2ec2687ab219
function wrap(i, I)
	if i == maximum(I)+1
		i = minimum(I)
	elseif i == minimum(I)-1
		i = maximum(I)
	end
	i
end

# ╔═╡ 1eaae315-4502-4b4a-9646-84ce085232cb
function αs((i, j))
	neighbors(i, j)
end

# ╔═╡ d70fcd31-6dfd-4dfd-a7f2-9807fc6b7201
αs(0, 0)

# ╔═╡ 73c8a825-dfc8-435f-8911-29c47370e80f
function α(k, i, j)
	αs(i, j)[k]
end

# ╔═╡ e8033a15-4993-4991-8740-92084467ecb3
function _αs(i, j)
	αs(i, j)[3, 4, 5, 0, 1, 2]
end

# ╔═╡ 8488470a-9583-4c71-8779-98558bf9f4c2
function _αs((i, j))
	_αs(i, j)
end

# ╔═╡ e3bd66bd-67bb-4a60-b90f-e7c3aca5fbe9
function _α(k, i, j)
	_αs(i, j)[k]
end

# ╔═╡ 0b475664-8915-4cc4-89c5-90788fc10091
function α(k, (i, j))
	α(k, i, j)
end

# ╔═╡ 35c4b183-9e11-4076-827c-902139733e02
function mirror(k, i, j; I, J)
	i, j = wrap(i, I), wrap(j, J)
	
	di, dj = α(k, i, j)
	i2, j2 = wrap(i+di, I), wrap(j+dj, J)
	
	_k = [3, 4, 5, 0, 1, 2][k+1]
	(i2, j2, _k)
end

# ╔═╡ 3ec5f04b-5b04-4a92-80fa-9ab4f4fec596
"""
    max_abs(x)

Returns the maximum absolute-value of the collection `x`. This function works on ranges (e.g., `UnitRange`, `IdOffsetRange`) as well as regular array types.
"""
function max_abs(x)
	maximum([abs(minimum(x)) abs(maximum(x))])
end

# ╔═╡ 4920893e-3a35-4c2c-abac-4ae723950d51
"""
    indices(OA, dim=nothing)

Returns the indices of `OA`, an `OffsetArray`, along the given dimension `dim`.

The returned object is a `UnitRange` if `dim` is given or a `Tuple` of `UnitRange` objects if `dim` is `nothing`.
"""
function indices(Sites, dim=nothing)
	if dim == nothing
		(minimum(a):maximum(a) for a in axes(Sites))
	else
		minimum(axes(Sites, dim)):maximum(axes(Sites, dim))
	end
end

# ╔═╡ 2ef1a611-9251-4b8d-9d77-2060ef9393da
function calc_outgoing!(Sites_out, Sites)
	I, J = indices(Sites)
	for i in I
		for j in J
			update_2body!(Sites_out, Sites, i, j) ||
			    update_3body!(Sites_out, Sites, i, j) ||
				update_passthrough!(Sites_out, Sites, i, j)
		end
	end
end

# ╔═╡ 9561a486-cab9-46d1-a15a-e0742bc98001
"""
    transport!(Sites, Step)

Let `Step` be a lattice whose particles are all outgoing. This function moves them to their destination sites, storing the result in `Sites`.
"""
function transport!(Sites, Step)
	I, J = indices(Step)
	Sites[:, :, :] .= 0

	###
	### YOUR CODE HERE ###
	###
	### BEGIN SOLUTION
	for i in I
		for j in J
			active = findall(Step[i, j, :])
			for a in active
				i2, j2, a2 = mirror(a, i, j; I=I, J=J)
				Sites[i2, j2, a2] = 1
			end
		end
	end
	### END SOLUTION
end

# ╔═╡ 10e7c3a8-ebab-4897-b457-ff7eb78916f9
"""
    initialize_simple_case!(Sites; case=0)

Some cases:
1. Boundary pass-through
2. Two-body collision
3. Three-body collision
4. Three-body pass-through
5. Wraparound "swap"
"""
function initialize_simple_case!(Sites; case=0)
	I, J = indices(Sites)
	M = maximum([max_abs(I) max_abs(J)])
	if case == 1 # boundary pass-through
		Sites[-M, 0, 0] = 1
	elseif case == 2 # 2-body
		Sites[0, 0, 0] = 1
		Sites[0, 0, 3] = 1
	elseif case == 3 # 3-body
		Sites[0, 0, 0] = 1
		Sites[0, 0, 2] = 1
		Sites[0, 0, 4] = 1
	elseif case == 4 # rando
		Sites[-M+1, 0, 0] = 1
		Sites[-M+1, 0, 2] = 1
		Sites[-M+1, 0, 3] = 1
	elseif case == 5 # rando
		Sites[-M, -M+1, [1 3]] .= 1
		Sites[-M+1, -M, 5] = 1
		Sites[ M-1, -M, 4] = 1
	else # default: illustrate where the edges are
		Sites[0, 0, :] .= 1
	end
end

# ╔═╡ 3407de4d-8916-4d80-ae95-6921fa1054ba
function simulate_step!(Sites)
	I, J = indices(Sites, 1), indices(Sites, 2)
	Step = copy(Sites)
	for i in I
		for j in J
			step_site!(Step, Sites, i, j)
		end
	end
	Sites[:, :, :] = Step
	transport!(Sites, Step)
end

# ╔═╡ 49e3115e-56d0-483b-b93f-9e11f08629db
function plot_lattice3(Sites; s=STEPSIZE, title="")
	I, J = indices(Sites, 1), indices(Sites, 2)
	grid = site_apply(Sites, sum)
	plt = heatmap(I, J, collect(transpose(grid)),
	              aspect_ratio=:equal, clim=(0, 6),
	              c=:Blues_7, title=title)
	hline!([0], color=:grey, linestyle=:dot, legend=nothing)
	vline!([0], color=:grey, linestyle=:dot, legend=nothing)
	plt
end

# ╔═╡ 3c7d92ad-6297-4bf5-8d6c-9dca1f3b7d9c
function generate_frame(Sites, t)
	I, J = indices(Sites, 2), indices(Sites, 3)
	num_particles = sum(Sites[t, I, J, :])
	plot_lattice3(Sites[t, :, :, :];
				  title=L"t=%$(t);\ n=%$(num_particles)")
end

# ╔═╡ efdbf204-9abe-496f-9c74-b76ab1b33111
function generate_animation(Sites; T=T)
	anim = @animate for t in T
		generate_frame(Sites, t)
	end
	gif(anim, fps=5)
end

# ╔═╡ a4ed6980-5868-487d-9eb2-772c989e49e9
function ind2coord(i, j; which=:both, s=STEPSIZE)
	if s == nothing
		x, y = i, j
	else
		r = s * sqrt(3/4)
		y = j * r
		if j%2 == 0 # even
			x = i * s
		else
			x = (i + 0.5) * s
		end
	end
	(x, y)
end

# ╔═╡ 936a23a5-6167-4147-90a7-87196ab1552d
function ind2coord((i, j); kwargs...)
	ind2coord(i, j; kwargs...)
end

# ╔═╡ 6a003973-f857-40f8-9fc1-c977ffa55d25
function plot_site(i, j, Sites=nothing; s=STEPSIZE,
	               add=true, label_coords=true, label_dirs=true,
	               kwargs...)
	x, y = ind2coord(i, j; s=s)
	plot!([x], [y], seriestype=:scatter, color=:red; kwargs...)
	if Sites == nothing
		alphas = 0:5
	else
		alphas = findall(Sites[i, j, :])
	end
	if label_coords
		annotate!(x, y, text(L"(%$(i),%$(j))", :left, :top, :red, 10))
	end
	for a in alphas
		di, dj = α(a, i, j)
		i2, j2 = i+di, j+dj
		x2, y2 = ind2coord(i2, j2; s=s)
		if label_coords
			annotate!(x2, y2, text(L"(%$(i2),%$(j2))", :left, :top, :black, 10))
		end
		if label_dirs
			annotate!(x2, y2, text(L"\alpha=%$(a)", :left, :bottom, :black, 10))
		end
	end
end

# ╔═╡ 28c59939-4ec1-41c9-b5dc-90e545e63f36
function create_sites(I, J; T=nothing)
	K = 0:length(αs(0, 0))-1
	if T == nothing
		Values = zeros(Bool, length(I), length(J), length(K))
		Sites = OffsetArray(Values, I, J, K)
	else
		Values = zeros(Bool, length(T), length(I), length(J), length(K))
		Sites = OffsetArray(Values, T, I, J, K)
	end
	Sites
end

# ╔═╡ 078e7283-b583-4328-bb68-1a7dd9e676c2
"""
    simulate(Sites0, T)

Inputs:
- `Sites0`, an initial condition for all lattice sites; `Sites0[i, j, a]` refers to bit-`a` of the site at logical position `(i, j)`.
- `T` is a unit-range variable indicating which time steps to evaluate.

Then, for each time step `t in T`, this function
- analyzes incoming particles and determine there outgoing directions;
- and then transports outgoing particles to the next site.

It returns `Sites[t, :, :, :]`, the sequence of simulation grids.
"""
function simulate(Sites0, T)
	# Create output grid: Sites[t, i, j, alpha]
	I, J = indices(Sites0, 1), indices(Sites0, 2)
	Sites = create_sites(I, J; T=T)

	# Initial condition
	t0, tmax = minimum(T), maximum(T)
	Sites[t0, :, :, :] = Sites0

	# Main simulation loop
	for t in t0:tmax-1
		# Incoming -> outgoing
		Sites_outgoing = Sites[t, :, :, :]
		calc_outgoing!(Sites_outgoing, Sites[t, :, :, :])

		# Transport outgoing particles to the next site
		Sites_next = copy(Sites_outgoing)
		transport!(Sites_next, Sites_outgoing)
		Sites[t+1, :, :, :] = Sites_next
	end
	Sites
end

# ╔═╡ ada051d1-49d4-40e1-a837-5c418014336b
function pair(i, j)
	(i, j)
end

# ╔═╡ 88ecd24a-dd16-4d07-b7f9-8a5a30cc00e6
function group_coords(I, J)
	I_coords = repeat(I, 1, length(J))
	J_coords = repeat(transpose(J), length(I), 1)
	OffsetArray(pair.(I_coords, J_coords), I, J)
end

# ╔═╡ 65a1f149-bfc2-49e6-9549-c0e3ce423976
function calc_coords(I, J; kwargs...)
	I_coords = repeat(I, 1, length(J))
	J_coords = repeat(transpose(J), length(I), 1)
	OffsetArray(ind2coord.(I_coords, J_coords; kwargs...), I, J)	
end

# ╔═╡ 8dff633c-2512-49e5-bc59-9aa1780c7e82
function plot_lattice(Sites;
	                  s=STEPSIZE, incoming=true, label_coords=false,
	                  title=nothing, arrows=nothing, markersize=nothing)
	I, J = indices(Sites)
	M = maximum([max_abs(I) max_abs(J)])

	if arrows == nothing
		arrows = maximum(size(Sites)) <= 21
	end
	if markersize == nothing
		if maximum(size(Sites)) <= 21
			markersize = 2
		else
			markersize = 0.5
		end
	end
	if title == nothing
		title = L"n=%$(sum(Sites))"
	end

	r = s * sqrt(3/4)
	XY = calc_coords(I, J; s=s)
	Xs = [xy[1] for xy in reshape(XY, length(XY))]
	Ys = [xy[2] for xy in reshape(XY, length(XY))]
	plt = plot(Xs, Ys, seriestype=:scatter, color=:grey,
		       size=(800, 600), aspect=1,
		       legend=nothing, markersize=markersize, title=title,
	           xlims=(minimum(Xs)-s*1.25, maximum(Xs)+s*1.25),
		       ylims=(minimum(Ys)-r*1.25, maximum(Ys)+r*1.25))
	hline!([0], color=:grey, linestyle=:dot)
	vline!([0], color=:grey, linestyle=:dot)

	xb_min, xb_max = minimum(Xs)-s/2, maximum(Xs)+s/2
	yb_min, yb_max = minimum(Ys)-r/2, maximum(Ys)+r/2
	hline!([yb_min, yb_max], color=:grey, linestyle=:dash)
	vline!([xb_min, xb_max], color=:grey, linestyle=:dash)

	if label_coords
		for i in I
			for j in J
				x, y = ind2coord(i, j; s=s)
				annotate!(x, y, text(L"(%$(i),%$(j))", :left, :top, 9))
			end
		end
	end

	edges = findall(Sites)
	edges_x = zeros(2, length(edges))
	edges_y = zeros(2, length(edges))
	dir_x = zeros(length(edges))
	dir_y = zeros(length(edges))
	for (e, k) in enumerate(edges)
		i, j, k = Tuple(k)
		αi, αj = α(k, i, j)
		edges_x[1, e], edges_y[1, e] = ind2coord(i, j; s=s)
		edges_x[2, e], edges_y[2, e] = ind2coord(i+αi, j+αj; s=s)
		(source, sink) = incoming ? (2, 1) : (1, 2)
		dir_x[e] = edges_x[sink, e] - edges_x[source, e]
		dir_y[e] = edges_y[sink, e] - edges_y[source, e]
	end

	if arrows
		source = incoming ? 2 : 1
		quiver!(edges_x[source, :], edges_y[source, :],
			    quiver=(dir_x, dir_y),
				color=:grey,
				arrow=arrow(:closed, 0.1))
	else
		plot!(edges_x, edges_y, color=:black)
	end

	plt
end

# ╔═╡ dc991c4a-11c9-4b25-9e50-8f27a8b8c1f2
begin
	S2 = create_sites(-3:2, -3:2)
	initialize_simple_case!(S2; case=2) # Try 1-5
	plt2 = plot_lattice(S2; arrows=true, incoming=true)
	plot_site(0, 0, S2)
	plt2
end

# ╔═╡ c729b897-bb84-4e80-b4e2-3f7480c20651
# Example: Mirror of bit α=2 at site (i=-2, j=0) => (i', j', α')
mirror(2, -2, 0; I=indices(S2, 1), J=indices(S2, 2))

# ╔═╡ 2b14df96-0721-45a9-9501-27afad4af458
begin
	S2_after = copy(S2)
	update_2body!(S2_after, S2, 0, 0)
	plt2_after = plot_lattice(S2_after; arrows=true, incoming=false)
	plot_site(0, 0, S2_after)
	plt2_after
end

# ╔═╡ 095312d0-81eb-4812-b61b-964d0aaaf5d9
# Three-body demo
begin
	S3 = create_sites(-3:2, -3:2)
	initialize_simple_case!(S3; case=3) # Try 1-5
	plt3 = plot_lattice(S3; arrows=true, incoming=true)
	plot_site(0, 0, S3)
	plt3
end

# ╔═╡ 996f93d9-edf3-4331-8735-88f76a5cbebb
begin
	S3_after = copy(S3)
	update_3body!(S3_after, S3, 0, 0)
	plt3_after = plot_lattice(S3_after; arrows=true, incoming=false)
	plot_site(0, 0, S3_after)
	plt3_after	
end

# ╔═╡ 983b744e-8e94-44a0-9392-786e58a3a66e
# Passthrough demo
begin
	S4 = create_sites(-3:2, -3:2)
	initialize_simple_case!(S4; case=4) # Try 1-5
	plt4 = plot_lattice(S4; arrows=true, incoming=true)
	plot_site(-2, 0, S4)
	plt4
end

# ╔═╡ cdb2122d-7d7a-4177-beeb-44776caa58a1
begin
	S4_after = copy(S4)
	update_passthrough!(S4_after, S4, 0, 0)
	plt4_after = plot_lattice(S4_after; arrows=true, incoming=false)
	plot_site(-2, 0, S4_after)
	plt4_after	
end

# ╔═╡ b8f34552-65ca-419d-9940-8013e2978262
# Transport demo: Incoming
begin
	S5 = create_sites(-3:2, -3:2)
	initialize_simple_case!(S5; case=5) # Try 1-5
	plt5 = plot_lattice(S5; arrows=true, incoming=true)
	plot_site(-2, -3, S5)
	plot_site(-3, -2, S5)
	plot_site(2, -3, S5)
	plt5
end

# ╔═╡ 60b3d814-c279-4ba4-95d2-2a93ccd199c9
# Transport demo: Outgoing
begin
	S5_outgoing = copy(S5)
	calc_outgoing!(S5_outgoing, S5)
	plt5_outgoing = plot_lattice(S5_outgoing; arrows=true, incoming=false)
	plot_site(-2, -2, S5_outgoing)
	plot_site(-2, -3, S5_outgoing)
	plot_site(-3, -2, S5_outgoing)
	plot_site(2, -3, S5_outgoing)
	plt5_outgoing
end

# ╔═╡ f97212e1-27ca-4ca3-99f1-ae4225537a20
# Transport demo: End
begin
	S5_end = copy(S5_outgoing)
	transport!(S5_end, S5_outgoing)
	plt5_end = plot_lattice(S5_end; arrows=true, incoming=true)
	plot_site(-2, -2, S5_end)
	plot_site(-2, -3, S5_end)
	plot_site(-3, -2, S5_end)
	plot_site(2, -3, S5_end)
	plt5_end
end

# ╔═╡ 1ebc624a-1253-48c0-b3c1-f56eec97c7c6
function illustrate_sites(M=3; case=0, s=STEPSIZE, label_coords=false)
	S = create_sites(-M:(M-1), -M:(M-1))
	initialize_simple_case!(S; case=case)
	plt = plot_lattice(S; s=s, arrows=false, label_coords=label_coords)
	S, plt
end

# ╔═╡ 2b848d43-c63a-49ca-b3d0-6e660ceb5273
begin
	_, plt0 = illustrate_sites(3; label_coords=true)
	plt0
end

# ╔═╡ f95cde71-1e7b-4021-a8bf-a1c3ef146f8e
begin
	_, plt1 = illustrate_sites(2)
	plot_site(0, 0)
	plt1
end

# ╔═╡ e7532f00-ce2f-40ad-801d-cb3c77841b87
function sample_edges(p, I, J; A=αs(0, 0))
	N = length(I) * length(J) * length(A)
	K = 0:length(A)-1
	Items = reshape([item for item in Iterators.product(I, J, K)], N)
	sample(Items, Int(round(p*N)))
end

# ╔═╡ 9a9d975b-3146-4622-bf29-b856ae560df0
function randomize_sites!(Sites, p; seed=nothing, kwargs...)
	I, J = indices(Sites, 1), indices(Sites, 2)
	if seed != nothing
		Random.seed!(3)
	end
	
	# Pick random edges
	edges = sample_edges(p, I, J; kwargs...)
	for (i, j, k) in edges
		Sites[i, j, k] = 1
	end
end

# ╔═╡ 624221b1-c7eb-4588-b040-717a99c07da7
begin
	Sites = create_sites(-M:M-1, -M:M-1)
	randomize_sites!(Sites, 0.5; seed=3)
end

# ╔═╡ 40565068-46f1-4dad-9006-c3d1b2f687ac
begin
	Sites_sim = simulate(Sites, T)
	size(Sites_sim)
end

# ╔═╡ e074d7c8-71b8-4dc2-b277-28bac17c8d62
generate_animation(Sites_sim)

# ╔═╡ 2ff103f1-35a4-4b9c-a817-54cab4f182dc
plot_lattice(Sites_sim[t_example, :, :, :];
	         arrows=false,
	         title=L"t=%$(t_example)" *
			       L",\ n=%$(sum(Sites_sim[t_example, :, :, :]))")

# ╔═╡ 6f495ae4-4a62-4c5f-aa4a-4613e77eaf0d
function check_sites(Sites; kwargs...)
	active = findall(Sites)
	for a in active
		i, j, k = Tuple(a)
		i_prime, j_prime, k_prime = mirror(k, i, j; kwargs...)
		@assert Sites[i_prime, j_prime, k_prime] == false "Eep! Check: (i,j,k)=($(i),$(j),$(k)) == $(Sites[i, j, k]) <== (i', j', k')=($(i_prime), $(j_prime), $(k_prime)) == $(Sites[i_prime, j_prime, k_prime])"
	end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
LaTeXStrings = "~1.3.0"
OffsetArrays = "~1.10.8"
Plots = "~1.25.11"
PlutoUI = "~0.7.35"
StatsBase = "~0.33.16"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a6552bfeab40de157a297d84e03ade4b8177677f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.12"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "13468f237353112a01b2d6b32f3d0f80219944aa"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "5c907bdee5966a9adb8a106807b7c387e51e4d6c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.11"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "85bf3e4bd279e405f91489ce518dedb1e32119cb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.35"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "de893592a221142f3db370f48290e3a2ef39998f"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "995a812c6f7edea7527bb570f0ac39d0fb15663c"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "6354dfaf95d398a1a70e0b28238321d5d17b2530"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═f80f5adc-97e5-11ec-341c-337dc6192b37
# ╟─a5dbdc32-d7a6-4e2b-8c47-0de9b90d7235
# ╟─915ee6cf-ff20-48b0-8803-1ff77f8bde93
# ╟─52d94124-6cc3-4af1-b010-63c668fd0378
# ╠═2b848d43-c63a-49ca-b3d0-6e660ceb5273
# ╟─fb6e0d76-7c3f-4f1d-a4f0-a78f39b8b7ca
# ╟─473ca08e-eff6-48a5-8981-afad0a9186d5
# ╠═078e7283-b583-4328-bb68-1a7dd9e676c2
# ╠═2ef1a611-9251-4b8d-9d77-2060ef9393da
# ╟─9b1088a7-adbf-47da-9c78-032b55e394e8
# ╠═f95cde71-1e7b-4021-a8bf-a1c3ef146f8e
# ╟─aa90a3c7-6c7a-4f2a-a08f-36cc7ede1b30
# ╠═dc991c4a-11c9-4b25-9e50-8f27a8b8c1f2
# ╟─6392b141-0548-48ac-91c8-b5589279aff3
# ╟─34c88f6b-97e3-4e9b-9301-d44269b0ded7
# ╟─3c4da1eb-15fe-4565-a3c1-cf5db155efe6
# ╠═2b14df96-0721-45a9-9501-27afad4af458
# ╟─e9c572e9-63ee-4e71-88a6-a16f2d079fdd
# ╠═4e699d98-2786-4aa5-90e9-b8a094ea7263
# ╟─095312d0-81eb-4812-b61b-964d0aaaf5d9
# ╟─996f93d9-edf3-4331-8735-88f76a5cbebb
# ╟─bf939585-832b-463c-8752-4ab79bbcffab
# ╠═0e9d40e8-a5e3-4c3a-9f8e-b485be8554e7
# ╟─983b744e-8e94-44a0-9392-786e58a3a66e
# ╠═cdb2122d-7d7a-4177-beeb-44776caa58a1
# ╟─ce6284b2-4fef-4ad5-998e-449db4a7a943
# ╠═c729b897-bb84-4e80-b4e2-3f7480c20651
# ╠═9561a486-cab9-46d1-a15a-e0742bc98001
# ╠═b8f34552-65ca-419d-9940-8013e2978262
# ╠═60b3d814-c279-4ba4-95d2-2a93ccd199c9
# ╠═f97212e1-27ca-4ca3-99f1-ae4225537a20
# ╟─0e90eb77-d635-4190-abfd-b2ad2198a520
# ╠═4d668f04-0b5e-4e0f-9d9a-dabef42d855a
# ╠═624221b1-c7eb-4588-b040-717a99c07da7
# ╠═40565068-46f1-4dad-9006-c3d1b2f687ac
# ╠═e074d7c8-71b8-4dc2-b277-28bac17c8d62
# ╟─c7bff650-7801-4d88-9a27-9153b7eca61f
# ╠═4075dd9c-88b8-46fb-99eb-d3e03e151a05
# ╠═2ff103f1-35a4-4b9c-a817-54cab4f182dc
# ╟─dc9c9d76-5967-43a2-958a-9c61ee9101bc
# ╠═1ebc624a-1253-48c0-b3c1-f56eec97c7c6
# ╟─10e7c3a8-ebab-4897-b457-ff7eb78916f9
# ╠═6a003973-f857-40f8-9fc1-c977ffa55d25
# ╠═d70fcd31-6dfd-4dfd-a7f2-9807fc6b7201
# ╠═efdbf204-9abe-496f-9c74-b76ab1b33111
# ╠═3c7d92ad-6297-4bf5-8d6c-9dca1f3b7d9c
# ╠═3407de4d-8916-4d80-ae95-6921fa1054ba
# ╠═df529624-20c5-4252-897c-aaaed4db673e
# ╠═8a7b06e0-5cdb-49b8-9658-c94db5399c6f
# ╠═73c8a825-dfc8-435f-8911-29c47370e80f
# ╠═e8033a15-4993-4991-8740-92084467ecb3
# ╠═e3bd66bd-67bb-4a60-b90f-e7c3aca5fbe9
# ╠═35c4b183-9e11-4076-827c-902139733e02
# ╠═4b7b3eeb-8fc7-4cb3-9b33-2ec2687ab219
# ╟─1eaae315-4502-4b4a-9646-84ce085232cb
# ╟─8488470a-9583-4c71-8779-98558bf9f4c2
# ╟─0b475664-8915-4cc4-89c5-90788fc10091
# ╟─3ec5f04b-5b04-4a92-80fa-9ab4f4fec596
# ╟─4920893e-3a35-4c2c-abac-4ae723950d51
# ╠═8dff633c-2512-49e5-bc59-9aa1780c7e82
# ╠═49e3115e-56d0-483b-b93f-9e11f08629db
# ╠═a4ed6980-5868-487d-9eb2-772c989e49e9
# ╠═936a23a5-6167-4147-90a7-87196ab1552d
# ╠═28c59939-4ec1-41c9-b5dc-90e545e63f36
# ╠═ada051d1-49d4-40e1-a837-5c418014336b
# ╠═88ecd24a-dd16-4d07-b7f9-8a5a30cc00e6
# ╠═65a1f149-bfc2-49e6-9549-c0e3ce423976
# ╠═e7532f00-ce2f-40ad-801d-cb3c77841b87
# ╠═9a9d975b-3146-4622-bf29-b856ae560df0
# ╠═6f495ae4-4a62-4c5f-aa4a-4613e77eaf0d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
