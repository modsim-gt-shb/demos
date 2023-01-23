### A Pluto.jl notebook ###
# v0.17.7

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

# ╔═╡ 49b1ec9f-00ce-4c4b-8416-80e09844107b
begin
	using LaTeXStrings
	using PlutoUI
	using CairoMakie
	using LinearAlgebra
	using DifferentialEquations
	
	md"""
	> This notebook is based on Steven Strogatz's, _Love affairs and differential equations (1988)_.
	>
	> - [[PDF link](http://ai.stanford.edu/~rajatr/articles/SS_love_dEq.pdf)]
	> - [[Official citation link](https://doi.org/10.1080/0025570X.1988.11977342)]
	>
	> Expand this cell to see the Julia modules that this notebook uses.
	"""
end

# ╔═╡ 50a510fe-8082-11ec-0138-934922c3e62b
md"""
# Star-crossed Lovers #

Valentine's Day is nearly upon us. It's time to build mathematical models of love, to maximize your chances of success!
"""

# ╔═╡ 70ff7b68-c5aa-4028-9aad-06844fae534a
begin
	# ========== Pick your celebrity couple here ========== #

	❤1, ❤2 = "Beyoncé", "Jay-Z"
	#❤1, ❤2 = "Shahid Kapoor", "Mira Rajput" # https://www.cosmopolitan.in/celebrity/news/a11161/9-most-controversial-bollywood-marriages-all-time
	#❤1, ❤2 = "Marie Antoinette", "Louis XVI" # https://inspirelle.com/11-famous-french-couples-love-better-worse/
	#❤1, ❤2 = "신민아", "김우빈" # "Shin Min-ah", "Kim Woo-bin" — https://www.scmp.com/magazines/style/celebrity/article/3157724/7-korean-celebrity-couples-who-overcame-tough-times-k-pop
	#❤1, ❤2 = "范冰冰", "李晨" # "Fan Bingbing", "Li Chen" — https://www.scmp.com/magazines/style/celebrity/article/3156705/5-chinese-celebrity-couples-who-dated-screen-then-real
	#❤1, ❤2 = "The Royal Formerly Known as Prince Henry", "Meghan Markle" # https://www.wonderwall.com/celebrity/couples/controversial-celebrity-marriages-386831.gallery
	
	# ========== Pick your celebrity couple here ========== #

	md"""
	## Pick your poison ##
	_(Or, pick your lovers... or your poisoned lovers?)_
	
	Let's look to a celebrity couple for inspiration. By default, we'll use [Beyoncé and Jay-Z](https://www.usmagazine.com/entertainment/pictures/beyonce-jay-zs-sexy-romance-2011113/), but feel free to change this default (by inspecting the code for this cell and editing it appropriately).

	Your celebrity couple: ❤1 ``\!\equiv\!`` **$(❤1)**, ❤2 ``\!\equiv\!`` **$(❤2)**
	"""
end

# ╔═╡ 5eaea167-90eb-4e3e-bb0b-718a5cb66187
md"""
## A formal model of love ##
_(Linear is best, obvi)_

Let ``x_1=x_1(t)`` be a real-valued function of time that measures how much $(❤1) loves $(❤2) at time ``t``. Define ``x_2=x_2(t)`` similarly for $(❤2). Suppose that positive values indicate love and negative values mean loathing. Zero values presumably mean indifference.

> You might argue instead that love is complex, that is, ``x_1`` and ``x_2`` should be complex-valued with both real and imaginary parts. That certainly rings true. But more on that later...
"""

# ╔═╡ f7b1d204-c49c-4e52-a651-a011ccb9a564
md"""Let's further suppose love's evolution through time can be modeled by a linear dynamical system:"""

# ╔═╡ 46a86e66-dff9-49d9-83e7-3d1edf9d479f
md"""
Pause for a moment and think through what different coefficient values might mean.

Here is one model of a couple as suggested by Strogatz. The more $(❤2) loves $(❤1), the more $(❤1) loves $(❤2) back. But $(❤2) is a fickle lover: the more $(❤1) loves $(❤2), the more $(❤2) dislikes $(❤1). These might correspond to the system where ``a_{12} = 1`` whereas ``a_{21} = -1`` (and ``a_{11} = a_{22} = 0``).

> Strogatz has several other entertaining suggestions for lover types. For instance, ``a_{21} > 0`` and ``a_{22} > 0`` might correspond to $(❤2) being an "eager beaver." Or ``a_{11} < 0`` and ``a_{12} > 0`` might correspond to $(❤1) being a "cautious lover."
"""

# ╔═╡ 6aecb035-6b9c-4983-87fd-689f38b34ef7
md"""
## Numerical simulation ##
"""

# ╔═╡ a943ebae-1563-4bfd-a396-9f6ed97a74e4
begin
	coefs = -1:0.25:1
#	slider_a11 = @bind a11 PlutoUI.Slider(-1:0.25:1; default=0)
#	slider_a12 = @bind a12 PlutoUI.Slider(-1:0.25:1; default=1)
#	slider_a21 = @bind a21 PlutoUI.Slider(-1:0.25:1; default=-1)
#	slider_a22 = @bind a22 PlutoUI.Slider(-1:0.25:1; default=0)
	slider_a11 = PlutoUI.Slider(coefs; default=0)
	slider_a12 = PlutoUI.Slider(coefs; default=1)
	slider_a21 = PlutoUI.Slider(coefs; default=-1)
	slider_a22 = PlutoUI.Slider(coefs; default=0)
	@bind A_entries confirm(PlutoUI.combine() do Child
	md"""
	|       | $(❤1)         | $(❤2)         |
	|:-----:|:-------------:|:-------------:|
	| $(❤1) | $(Child("a11", slider_a11)) | $(Child("a12", slider_a12)) |
	| $(❤2) | $(Child("a21", slider_a21)) | $(Child("a22", slider_a22)) |
	"""
	end)
end

# ╔═╡ ad85b5af-86dd-435a-bc3e-283057013ad4
md"""
Use the sliders below to set the coefficients that best describe your celebrity couple.

Sliders correspond to the limits and steps, $(coefs).
"""

# ╔═╡ dfdd3e0a-9ed3-4e09-a006-584b588e9127
begin
	slider_x1_0 = PlutoUI.Slider(-1:0.2:1; default=1)
	slider_x2_0 = PlutoUI.Slider(-1:0.2:1; default=1)
	@bind x_0 confirm(PlutoUI.combine() do Child
		md"""
		Use the sliders to pick an initial condition, and then let's solve this system to see the "trajectory of love" for our star-crossed celebrity lovers.
		
		|       |                               |
		|:-----:|:-----------------------------:|
		| $(❤1) | $(Child("x1_0", slider_x1_0)) |
		| $(❤2) | $(Child("x2_0", slider_x2_0)) |
		"""
	end)
end

# ╔═╡ 801699c7-027e-4156-b81c-3d5aece5fada
md"""
Use the next slider to pick a maximum simulation time.

|                            |                                                   |
|:--------------------------:|:-------------------------------------------------:|
| ``t_{\small\mathrm{max}}`` | $(@bind t_max confirm(PlutoUI.Slider(1:1:20; default=10))) |
"""

# ╔═╡ eac57cbc-f690-46a3-8ee3-ab663f6b22d2
L"""t_{\small\textrm{max}} = %$(t_max)"""

# ╔═╡ 784c4d9e-009b-4853-ab69-fff007fa482d
md"""Here is the trajectory of our celebrity lovers over time."""

# ╔═╡ 4c93a57c-651d-4a5d-b0c8-4b66d5f32a31
md"""
## Eigenpairs Analysis ##
"""

# ╔═╡ 8c4ec357-df54-4003-8fca-aa87ba77ccdd
md"""The eigenvalues of ``A`` are"""

# ╔═╡ cdfbcbde-286e-415f-8cc6-6107afae7488
md"""and the eigenvectors are the columns of"""

# ╔═╡ d7a22bdc-3a79-4f57-bb2e-de7d9342b8eb
md"""
## Love Flows ##
"""

# ╔═╡ 9c2a069f-09b5-43de-8b3a-52f2c53b90a3
md"""**Postscript: Auxiliary functions.**"""

# ╔═╡ d7e9af87-d749-405e-80e7-d6f3668d7955
function real_to_string(x; r=3)
	if x == trunc(x)
		string(trunc(Int, x))
	elseif r != nothing
		string(round(x, digits=r))
	else
		string(x)
	end
end

# ╔═╡ bdab2ce9-2f77-4dc4-9677-1f2e500dfa44
function complex_to_string(x)
	re, im = real(x), imag(x)
	parts = []
	if re != 0
		append!(parts, real_to_string(re))
		if im == 1
			append!(parts, " + i")
		elseif im == -1
			append!(parts, " - i")
		elseif im > 0
			append!(parts, " + " * real_to_string(im) * "i")
		elseif im < 0
			append!(parts, " - " * real_to_string(abs(im)) * "i")
		end
	elseif im != 0
		if im == 1
			append!(parts, "i")
		elseif im == -1
			append!(parts, "-i")
		else
			append!(parts, real_to_string(im) * "i")
		end
	end
	s = join(parts)
	if s == ""
		"0"
	else
		s
	end
end

# ╔═╡ d8613c51-3646-4941-93ae-e1cfd0d23fb1
function complex_to_string(s::String)
	s
end

# ╔═╡ 452529c3-9194-4c8b-b2ca-9fdae3781ef7
function matrix_to_string(A)
	rows = []
	for i in 1:size(A, 1)
		row = []
		for j in 1:size(A, 2)
			append!(row, [complex_to_string(A[i, j])])
		end
		row_str = join(row, " & ")
		append!(rows, [row_str])
	end
	rows_str = join(rows, " \\\\ ")
	open_str = "\\left[\\begin{matrix}"
    close_str = "\\end{matrix}\\right]"
	join([open_str, rows_str, close_str], " ")
end

# ╔═╡ 452b1cf9-06d6-42c1-bc42-c7471bc5c376
begin
	x1_0, x2_0 = x_0.x1_0, x_0.x2_0
	latexstring("x(0) = " * matrix_to_string([x1_0; x2_0]))
end

# ╔═╡ 79d3e426-708e-4df7-bc43-ad3b04baa268
function symbolic_mat(m, n; base="a", indcomma=false)
	indcomma_str = indcomma ? "," : ""
	A = nothing
	for i in 1:m
		a_i = Array{String}(undef, 1, n)
		if n > 1
			for j in 1:n
				if m > 1
					a_i[j] = base * "_{$(i)$(indcomma_str)$(j)}"
				else
					a_i[j] = base * "_{$(j)}"
				end
			end
		else
			a_i = base * "_{$(i)}"
		end
		if A == nothing
			A = a_i
		else
			A = [A; a_i]
		end
	end
	A
end

# ╔═╡ e8e4b37f-4fde-45b9-86d2-8a04e2913cc3
begin
	x_str = matrix_to_string(["x_1";
		                      "x_2"])
	F_str = matrix_to_string(["a_{11} x_1 + a_{12} x_2";
	                          "a_{21} x_1 + a_{22} x_2"])
	A_str = matrix_to_string(symbolic_mat(2, 2; base="a"))
	
	latexstring("D" * x_str * " = " * F_str * " = " * A_str * x_str * "\\Longrightarrow Dx = Ax")
end

# ╔═╡ 9f276336-3da4-4ce2-b78e-338a6e995c3f
begin
	function build_term(c; v="x")
		if c != 0
			if c == 1
				x = v
			elseif c == -1
				x = "-" * v
			else
				x = complex_to_string(c) * " " * v
			end
		else
			x = ""
		end
		x
	end
	
	function build_terms(coefs; base="x")
		terms = Array{String,1}()
		for (k, c) in enumerate(coefs)
			t = build_term(c; v=base * "_" * string(k))
			if t != ""
				append!(terms, [t])
			end
		end
		terms, join(terms, " + ")
	end	

	function system_latex(a11, a12, a21, a22)
		_, a1_terms_str = build_terms([a11, a12])
		_, a2_terms_str = build_terms([a21, a22])
		a11_str = complex_to_string(a11)
		a12_str = complex_to_string(a12)
		a21_str = complex_to_string(a21)
		a22_str = complex_to_string(a22)

		x_str = matrix_to_string(symbolic_mat(2, 1; base="x"))
		F_str = matrix_to_string([a1_terms_str; a2_terms_str])
		final_str = "D" * x_str * " = " * F_str * "\\quad \\Longrightarrow \\quad Dx = Ax, \\quad A = " * matrix_to_string([a11_str a12_str; a21_str a22_str])
		latexstring(final_str)
	end
end

# ╔═╡ 2a99271f-f9a1-46d9-9a4e-843c1b25657c
system_latex(0, 1, -1, 0)

# ╔═╡ 77cfe6d5-622b-43bc-8b4a-1c7ed447c3fe
begin
	a11 = A_entries.a11
	a12 = A_entries.a12
	a21 = A_entries.a21
	a22 = A_entries.a22
	system_latex(a11, a12, a21, a22)
end

# ╔═╡ e7f03841-4b20-42e3-b4e9-b21aa1a7fcf6
A = [a11 a12; a21 a22]

# ╔═╡ 72bb7e9f-b485-4724-9f2e-4f53621d4d2b
(λ1, λ2), Q = eigen(A);

# ╔═╡ c898725c-7cc5-43d6-8ba0-0d3abec3b6c5
begin
	λ1_str = complex_to_string(λ1)
	λ2_str = complex_to_string(λ2)
	L"""\lambda_1 = %$(λ1_str) \qquad \lambda_2 = %$(λ2_str)"""
end

# ╔═╡ e704ecd0-53c3-4e32-9de9-20684f272e49
begin
	prefix_str = "\\begin{align}"
	Q_str = "Q &= " * matrix_to_string(Q)
	q1_str = "q_1 &= " * matrix_to_string(Q[:, 1])
	q2_str = "q_2 = " * matrix_to_string(Q[:, 2])
	suffix_str = "\\end{align}"
	final_str = prefix_str * Q_str * "\\\\ \\quad \\implies \\ " * q1_str * ", \\ " * q2_str * suffix_str
	latexstring(final_str)
end

# ╔═╡ bb733454-f2f9-47e6-9ace-7fa2bd69aed7
begin
	xv = [0; 0]
	yv = [0; 0]
	q1, q2 = real(Q[:, 1]), real(Q[:, 2])
	if q1 == q2
		q2[:] = [cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)] * q1
	end
	uv = [q1[1]; q2[1]]
	vv = [q1[2]; q2[2]]

	fig_phase = Figure(fontsize=20, resolution=(640, 480))
	ax_phase = Axis(fig_phase[1, 1],
					title=L"Phase field, $Dx = Ax$",
					xlabel=❤1, ylabel=❤2)
	streamplot!((x) -> Point2(A*x), -2..2, -2..2)
	hlines!(ax_phase, 0, color=:gray, linestyle=:dot)
	vlines!(ax_phase, 0, color=:gray, linestyle=:dot)
	arrows!(ax_phase, xv, yv, uv, vv, color=:darkgray, linewidth=3, arrowsize=15)
	text!(ax_phase, L"q_1", position=(q1[1], q1[2]))
	text!(ax_phase, L"q_2", position=(q2[1], q2[2]))
	Colorbar(fig_phase[1, 2])
	ax_phase.aspect = AxisAspect(1)
	fig_phase
end

# ╔═╡ 8d7d8e3b-f42d-48d5-a034-bcbe61474ad0
md"""**ODE solver results.**"""

# ╔═╡ 59be2873-d6fa-4986-81f4-6c02ca2f03d2
begin
	function F!(dx, x, p, t)
		dx[:] = p[1] * x
	end
	x0 = [x1_0; x2_0]
	tspan = (0, t_max)
	prob = ODEProblem(F!, x0, tspan, [A], saveat=0.1)
	sol = solve(prob);
end

# ╔═╡ c33e5bbc-29c8-4241-99c8-7bacb28a5456
begin
	xlim_traj = ceil(maximum(hcat([e[1] for e in sol.u], [e[2] for e in sol.u]))*10)/10
	trajectory = Observable(Point2f[(sol.u[1][1], sol.u[1][2])])
	fig_traj, ax_traj = scatter(Point2f[(sol.u[1][1], sol.u[1][2])], marker=:star5, markersize=40)
	scatter!(ax_traj, trajectory)
	hlines!(ax_traj, 0, color=:gray, linestyle=:dashdot)
	vlines!(ax_traj, 0, color=:gray, linestyle=:dashdot)
	ax_traj.xlabel = ❤1
	ax_traj.ylabel = ❤2
	ax_traj.aspect = AxisAspect(1)
	limits!(ax_traj, -xlim_traj, xlim_traj, -xlim_traj, xlim_traj)
	frames_traj = 2:size(sol.t, 1)
#	record(fig_traj, "x.mp4", frames_traj; framerate=30) do k
	CairoMakie.Makie.Record(fig_traj, frames_traj) do k
		new_point = Point2f(sol.u[k][1], sol.u[k][2])
		trajectory[] = push!(trajectory[], new_point)
		ax_traj.title = "t = " * real_to_string(sol.t[k])
	end
end

# ╔═╡ 6b153d9a-04f7-43be-9bce-ec8558102646
md"""
**Postscript.** This Julia/Pluto notebook was prepared for use in the Georgia Tech course on Modeling and Computer Simulation, CSE 6730 / CX 4230, Spring 2022.
"""

# ╔═╡ Cell order:
# ╟─50a510fe-8082-11ec-0138-934922c3e62b
# ╟─49b1ec9f-00ce-4c4b-8416-80e09844107b
# ╟─70ff7b68-c5aa-4028-9aad-06844fae534a
# ╟─5eaea167-90eb-4e3e-bb0b-718a5cb66187
# ╟─f7b1d204-c49c-4e52-a651-a011ccb9a564
# ╟─e8e4b37f-4fde-45b9-86d2-8a04e2913cc3
# ╟─46a86e66-dff9-49d9-83e7-3d1edf9d479f
# ╟─2a99271f-f9a1-46d9-9a4e-843c1b25657c
# ╟─6aecb035-6b9c-4983-87fd-689f38b34ef7
# ╟─ad85b5af-86dd-435a-bc3e-283057013ad4
# ╟─a943ebae-1563-4bfd-a396-9f6ed97a74e4
# ╟─77cfe6d5-622b-43bc-8b4a-1c7ed447c3fe
# ╟─dfdd3e0a-9ed3-4e09-a006-584b588e9127
# ╟─452b1cf9-06d6-42c1-bc42-c7471bc5c376
# ╟─801699c7-027e-4156-b81c-3d5aece5fada
# ╟─eac57cbc-f690-46a3-8ee3-ab663f6b22d2
# ╟─784c4d9e-009b-4853-ab69-fff007fa482d
# ╟─c33e5bbc-29c8-4241-99c8-7bacb28a5456
# ╟─4c93a57c-651d-4a5d-b0c8-4b66d5f32a31
# ╟─e7f03841-4b20-42e3-b4e9-b21aa1a7fcf6
# ╠═72bb7e9f-b485-4724-9f2e-4f53621d4d2b
# ╟─8c4ec357-df54-4003-8fca-aa87ba77ccdd
# ╟─c898725c-7cc5-43d6-8ba0-0d3abec3b6c5
# ╟─cdfbcbde-286e-415f-8cc6-6107afae7488
# ╟─e704ecd0-53c3-4e32-9de9-20684f272e49
# ╟─d7a22bdc-3a79-4f57-bb2e-de7d9342b8eb
# ╟─bb733454-f2f9-47e6-9ace-7fa2bd69aed7
# ╟─9c2a069f-09b5-43de-8b3a-52f2c53b90a3
# ╟─d7e9af87-d749-405e-80e7-d6f3668d7955
# ╟─bdab2ce9-2f77-4dc4-9677-1f2e500dfa44
# ╟─d8613c51-3646-4941-93ae-e1cfd0d23fb1
# ╟─452529c3-9194-4c8b-b2ca-9fdae3781ef7
# ╟─79d3e426-708e-4df7-bc43-ad3b04baa268
# ╟─9f276336-3da4-4ce2-b78e-338a6e995c3f
# ╟─8d7d8e3b-f42d-48d5-a034-bcbe61474ad0
# ╟─59be2873-d6fa-4986-81f4-6c02ca2f03d2
# ╟─6b153d9a-04f7-43be-9bce-ec8558102646
