### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ b5204be0-e2ad-4645-8440-ea0e82f9e205
begin
	using Plots
	using MAT
	using LinearAlgebra
	using LaTeXStrings
	using TSVD
end

# ╔═╡ 0f352c21-f8ec-4d60-9f89-acfb061d6350
md"""
# Dynamical systems theory meets data analysis: dynamic mode decomposition (DMD) #

This notebook is a Julia adaptation of an example taken from the book,

* _Dynamic Mode Decomposition: Data-driven Modeling of Complex Systems,_ by N. Kutz, S. L. Brunton, B. W. Brunton, J. L. Proctor (2016). http://www.dmdbook.com/
"""

# ╔═╡ 5be85952-37d6-4dba-a859-9b4f4e3e1e14
md"""
# Example: Fluid flowing around a cylinder #

The plot and movie below are the output of a computer simulation of a fluid flowing around a cylinder from the left side of the plot toward the right. Ripples form, oscillate, and dissipate over time.
"""

# ╔═╡ 9996d46a-881d-4350-bd31-5c49df04228c
md"""
A simulation, where the (nonlinear) governing equations are given, generated these outputs. But what if the scenario were reversed? Suppose we have these outputs—perhaps as observations of a system—but don't know the governing equations. Can we still learn something about the dynamics or physics of the system? The dynamic mode decomposition, or DMD, was invented for this setting.
"""

# ╔═╡ f7c19246-5f95-47e0-9fa5-fea212046ee0
md"""
# Learning a linear model from data #

**Problem statement.** Consider a nonlinear dynamical system,

$$\frac{dx}{dt} = f(x),$$

where $x = x(t)$ is a vector of length $m$ that represents the state of the system at time $t$. Suppose we observe $x$ at a sequence of time points, $x_k = x(t_k)$, where $t_k = k \Delta t$ for $1 \leq k \leq n$.

If we knew $f(x)$, we could try to generate one observation $x_{k+1}$ from the previous one $x_k$ by integrating the system:

$$x_{k+1} = x_k + \int_{t_k}^{t_{k+1}} f(x(\tau)) d\tau.$$

But suppose we don't know $f(x)$. How can we use the observations to learn something about it? Let's apply a familiar idea from earlier in the class: approximate an observation $x_{k+1}$ from a previous one, $x_k$, using a _locally linear_ approximation, $x_{k+1} \approx A x_k$, for some unknown $A$. The problem is to determine $A$.
"""

# ╔═╡ 3a09cfee-5a6f-411b-9236-32fa02351948
md"""
**An idea.** If such an $A$ exists and holds for all observations, we can try to use those observations to estimate it.

Let $X \equiv [\begin{matrix} x_1 & x_2 & \cdots & x_{n-1} \end{matrix}]$ be the data matrix containing the first $n-1$ samples and let $Y \equiv [\begin{matrix} x_2 & x_3 & \cdots & x_n \end{matrix}]$ be the last $n-1$ samples, i.e., the shift of $X$ by 1 to the left (plus the last observation). The action of $A$ is to evolve the system step-by-step (observation-to-observation), i.e.,

$$Y \approx A X.$$

One way to pick $A$ is to solve a _linear least-squares regression problem_, where we choose a matrix $A_*$ to minimize the residual sum-of-squares:

$$A_* = \underset{A}{\arg\!\min} \|Y - AX\|_F^2.$$

This problem has a known analytic solution, which is $A_* = Y X^{\dagger}$, where $X^{\dagger}$ is the [_Moore-Penrose pseudoinverse_](https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse). There are different ways to compute it, including via $QR$ decomposition and the [_singular value decomposition_ (SVD)](https://en.wikipedia.org/wiki/Singular_value_decomposition).
"""

# ╔═╡ 295d1dd8-983d-443f-ab6a-28cce3125f37
md"""
**Refining the idea using the SVD.** Let's use the SVD to calculate the pseudoinverse and solve the linear least-squares regression problem.

Let $X$ be an $m \times n$ matrix, and let $s = \min(m, n)$ be the smaller of these two dimensions. The SVD of $X$ is the factorization $X = U \Sigma V^T$, where
- the matrix $U$ is an $m \times s$ matrix such that $U^T U = I$;
- the matrix $V$ is an $n \times s$ matrix such that $V^T V = I$; and
- the matrix $\Sigma = \mathrm{diag}(\sigma_1, \sigma_2, \ldots, \Sigma_s)$ is a $s \times s$ _diagonal_ matrix, where $\sigma_k \geq 0$.

The columns of $U$ are called the _left singular vectors_; the columns of $V$ are the _right singular vectors_; and the entries of $\Sigma$ are the _singular values_. By convention, the singular values are ordered from largest to smallest, $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_s$.

> _Note 1:_ This definition is the [_thin_ or _economy_ SVD](https://en.wikipedia.org/wiki/Singular_value_decomposition#Thin_SVD). There are other variations.
>
> _Note 2:_ If you don't remember or appreciate the SVD, take a moment now to do so. It is truly magical. The special structure of the SVD admits many manipulations, making it the "Swiss Army Knife" of linear algebra.

Armed with the SVD, observe the following:

$$\begin{eqnarray}
Y & \approx & A X = A (U \Sigma V^T) \\
\implies \quad A & \approx & Y V \Sigma^{-1} U^T.
\end{eqnarray}$$

Thus, given the first $n-1$ observations $X$ and last $n-1$ observations $Y$, our algorithm is simply as follows.

1. Calculate the SVD, $X = U \Sigma V^T$.
2. Form $A = Y V \Sigma^{-1} U^T$.
3. Predict any observation $x_{k+1} = A x_k$.

Voilá! 
"""

# ╔═╡ 804ab926-f054-4faf-8b00-63e1e8bf6250
md"""
**The catch(es).** From a computational perspective, our approach so far suffers from two potential issues.

1. Calculating the SVD could be expensive.
2. The matrix $A_*$ might be huge!

To see why these are issues, recall that there are $n$ observations of length $m$ each, or $\mathcal{O}(m n)$ elements in total. Typically, $m \gg n$. In the motivating movie, each frame is about $500 \times 200$ pixels, meaning $m=500 \cdot 200 = 100,000$.

Calculating the full SVD costs $\mathcal{O}(\min\{m n^2, m^2 n\})$ operations. Moreover, the size of $A$ is $m \times m$, thereby requiring $\mathcal{O}(m^2)$ storage. It doesn't seem like our approach will scale.

Luckily, there are two insights that will rescue us.
"""

# ╔═╡ 14a4e646-b775-4b4f-8989-ba0aa3603c30
md"""
**Insight 1: We might not need the full SVD.**

The special form of the SVD is not its only magic. A second is its "optimal approximation" property.

Suppose $X = U \Sigma V^T$ be the SVD of $X$. Let $\Sigma_r$ be the diagonal matrix with the top $r$ singular values of $\Sigma$ ($r \leq \min\{m, n\}$), and let $U_r$ and $V_r$ denote the associated singular vectors, i.e., $U_r$ holds the first $r$ columns of $U$ and $V_r$ the first $r$ columns of $V$. The _$r$-truncated SVD_ is the product $U_r \Sigma_r V_r^T$.

Consider any rank-$r$ decomposition of $X$. That is, suppose $X = FG^T$ where $F$ is $m \times r$ and $G$ is $n \times r$ with both $F$ and $G$ being of full rank (rank $r$). Then,

$$\min_{F, G} \|X - FG^T\|_F^2 = \|X - U_r \Sigma_r V_r^T\|_F^2 = \sigma_{r+1}^2 + \sigma_{r+2}^2 + \cdots + \sigma_s^2.$$

This statement succinctly summarizes two remarkable facts:

1. Keeping the top $r$ singular vectors of the SVD—the ones associated with the $r$ largest singular values—leads to the [best rank-$r$ approximation of $X$](https://en.wikipedia.org/wiki/Low-rank_approximation).
2. The sum of squares of the last $s-r$ singular values measures the approximation error.

So, instead of calculating the _exact_ SVD, $X = U \Sigma V^T$, we can approximate it by $X \approx U_r \Sigma_r V_r^T$. Its cost will be $\mathcal{O}(\max\{m, n\} r^2)$ operations. This approximation will work well if $r$ can be kept small. A small $r$ will be possible when the singular values decay rapidly (exponentially)—recall that the approximation error is the sum of the squares of the truncated singular values.
"""

# ╔═╡ c77885c5-c6d2-4ce7-9b84-d49d7a06028f
md"""
**Insight 2: We don't need $A$.** Recall how we analyzed these systems earlier in our class. Given $f(x)$, we could determine a linear approximation $A x$ near a point $x$. But our analysis didn't use $A$ directly; rather, we turned to $A$'s eigenvalues and eigenvectors. For a linear system—or a linear approximation to a nonlinear system—the eigenvectors and eigenvalues revealed the fundamental behavior of solution trajectories.

So what we need is not $A$ but insight into the _action_ of $A$ via its eigendecomposition.
"""

# ╔═╡ a6386cb2-e244-44de-b7a6-698c2d458dee
md"""
# The Dynamic Mode Decomposition (DMD) algorithm #

The dynamic mode decomposition (DMD) algorithm estimates the action of $A$. The key idea is to approximate $A$ by a much smaller matrix $A_r$, which we will refer to as $A$'s "reduced-order model." Our choice of $A_r$ will be judicious: we will ensure that the eigenvalues of $A_r$ are also eigenvalues of $A$. We will also show how to construct the corresponding eigenvectors of $A$ from those of $A_r$. Thus, we will be able to get the information about $A$ that we need _without_ having to construct it explicitly.
"""

# ╔═╡ af98ac0a-7975-4f1b-886e-ae255a85c36e
md"""
**Starting point: the truncated SVD.** Per "Insight 1," assume we have the truncated SVD of $X$, so that $X \approx U_r \Sigma_r V_r^T$. Thus, $A \approx Y V_r \Sigma_r^{-1} U_r^T$. The value $r$ is a tuning parameter we need to choose.
"""

# ╔═╡ b35fbda6-473f-4f0b-a90c-74a4c52982b8
md"""
**An eigenvalue-preserving projection.** However, we don't want to form $A$—nor should we need to, per "Insight 2." Instead, consider a reduced-order model, $A_r$, of $A$, formed by

$$A_r = U_r^T A U_r \approx U_r^T Y V_r \Sigma_r^{-1}.$$

The matrix $A_r$ is $r \times r$, and thus much smaller than the $m \times m$ matrix $A$. More importantly, while we don't want to form $A$ itself, we can still calculate $A_r$ without it: since $A_r$ by projecting $A$ onto $U_r$, we no longer need to form an intermediate $m \times m$ object.

This project is only useful if $A_r$ preserves the important information we need from $A$. To check that, suppose we calculate an eigendecomposition of $U_r^T Y V_r \Sigma_r^{-1}$:

$$\begin{eqnarray}
               U_r^T Y V_r \Sigma_r^{-1} & = & W_r \Lambda_r W_r^{-1} \\
\implies \quad U_r^T Y V_r \Sigma_r^{-1} W_r & = & W_r \Lambda_r.
\end{eqnarray}$$
"""

# ╔═╡ b9758b50-a41b-4eb1-a5d3-f7b9f0cdf096
md"""
We would like for the eigenvalues $\Lambda_r$ to be eigenvalues of $A$, too. If that were true, then there would be eigenvectors $\Phi_r$ corresponding to $\Lambda_r$ such that

$$\Phi_r \Lambda_r = A \Phi_r.$$

Recall that $A \approx Y V_r \Sigma_r^{-1} U_r^T.$ If you stare long enough at this relation and the approximate eigendecomposition of $A_r$, you might try the following. Let $\Phi_r \equiv Y V_r \Sigma_r^{-1} W_r$. Then,

$$\begin{eqnarray}
Y V_r \Sigma_r^{-1} U_r^T \Phi_r &    =    & Y V_r \Sigma_r^{-1} \underbrace{U_r^T (Y V_r \Sigma_r^{-1}}_{\approx A_r} W_r) \\
                                 & \approx & Y V_r \Sigma_r^{-1} \underline{A_r W_r} \\
                                 & \approx & \underline{Y V_r \Sigma_r^{-1} W_r} \Lambda_r \\
                                 &    =    & \Phi_r \Lambda_r.
\end{eqnarray}$$

Thus, $\Phi_r \equiv Y V_r \Sigma_r^{-1} W_r$ are approximate eigenvectors of $A$ corresponding to the (approximate) eigenvalues $\Lambda_r$.
"""

# ╔═╡ e7e11aa5-95c1-4569-9d39-2218ab9d6a0a
md"""
**From modes to dynamics.** The eigenvectors $\Phi_r$ will be (approximate) fundamental modes of solutions to the dynamical system. Thus, any solution will be a linear combination of these modes. That is, each observation $x_k$ may be written as

$$x_k \approx b_1 \lambda_1^{k-1} \phi_1 + b_2 \lambda_2^{k-1} \phi_2 + \cdots + b_r \lambda_r^{k-1} \phi_r \equiv \Phi_r b,$$

where the vector $$b$$ are coefficients to be determined. We can determine them from the first observation as follows.

$$\begin{eqnarray}
        x_1 & \approx & \Phi_r b \\
        x_1 & \approx & Y V_r \Sigma_r^{-1} W_r b \\
  U_r^T x_1 & \approx & \underline{U_r Y V_r \Sigma_r^{-1}} W_r b \\
  U_r^T x_1 & \approx & A_r W_r b \\
  U_r^T x_1 & \approx & W_r \Lambda_r b \\
\implies  b & \approx & (W_r \Lambda_r)^{-1} (U_r^T x_1).
\end{eqnarray}$$
"""

# ╔═╡ bec98e48-3286-4048-b383-9e661557bc07
md"""
The $k$-th observation will be approximated as $x_k \approx \Phi_r \Lambda_r^{k-1} b$. Since we assumed uniformly spaced time points, $t_k = k \Delta t$, then a simple way to interpolate to any arbitrary $t$ is

$$x(t) \approx \Phi_r \Lambda_r^{\frac{t}{\Delta t} - 1} b.$$
"""

# ╔═╡ 28f179c8-3796-4b27-a929-8ccec3d9396d
md"""
**Algorithm (summary).** In summary, here are all the steps of the DMD algorithm.

- **Step 1:** Calculate the $r$-truncated singular value decomposition (SVD) of $X$: $U_r \Sigma_r V_r^T \approx X$
- **Step 2:** Form the approximate projected system operator, $A_r \approx U_r^T Y V_r \Sigma_r^{-1}$. _Claim:_ The eigenvalues of $A_r$ will be approximate eigenvalues of $A$, too.
- **Step 3:** Calculate $A_r \approx W_r \Lambda_r W_r^{-1}$, the approximate eigendecomposition of $A_r$.
- **Step 4:** Calculate $\Phi_r = Y V_r \Sigma_r^{-1} W_r$. These will be the approximate eigenvectors of $A$ corresponding to the approximate eigenvalues $\Lambda_r$.
- **Step 5:** Calculate $b = (W_r \Lambda_r)^{-1} (U_r^T x_1)$. These will be the coefficients of the modes characterizing the system dynamics.
"""

# ╔═╡ 2888f19d-1706-4311-a607-e7c38ed8a194
function dmd(observations; r=50)
	m, n = size(observations)
	r = min(r, n)
	X = observations[:, 1:end-1]
	Y = observations[:, 2:end]
	Ur, σr, Vr = tsvd(X, r) # `σr` is a vector (not a diag matrix)
	Ar = transpose(Ur) * Y * Vr ./ transpose(σr)
	λr, Wr = eigen(Ar) # `λr` is a vector (not a diag matrix)
	Φr = Y * (Vr ./ transpose(σr)) * Wr
	return (Φr, Wr, λr, Ur, σr, Vr)
end

# ╔═╡ 00aff0d3-aefa-459b-9b37-2ea648c3dedf
function dmd_fit(x, Ur, Wr, λr)
	b = (Wr .* λr) \ (transpose(Ur) * x)
end

# ╔═╡ 9274a087-943e-4c12-96ed-db26e76b3aae
md"""
# Experiments #

Let's apply DMD to the vortex flow observations.
"""

# ╔═╡ ad41f62c-bf60-4e9c-b300-f0ce837d97d2
md"""
**Singular values.** For the approximation to be a good one, we need the observations matrix $X$ to have low rank. We can check that by inspecting its singular values, to see that they decay "quickly" (ideally, exponentially). Indeed, that appears to be the case!
"""

# ╔═╡ a95ac675-3829-4c77-82a4-808cd939e3fc
md"""
**Modes.** Let's also inspect the eigenvalues and eigenvectors of $A_r$, the reduced-order model. Since our vortex had oscillatory and decaying behaviors, we should see complex eigenvalues. And it did not appear that the flows were "blowing up." So we might also see that the eigenvalues are bounded in magnitude by 1.
"""

# ╔═╡ 7b4b73e6-6b8e-42a7-9062-d7d226980e46
md"""
Here is a mode. Since it has the same dimensions as an input observation, we can draw the mode in the same spatial coordinates.
"""

# ╔═╡ 47f26dba-c61f-43ad-a5dd-44b4ce0676fb
md"""
# Reconstructing the signal #

As the last step, let's reconstruct the flow from just the first observation, using the eigendecomposition.
"""

# ╔═╡ dd616c1d-969a-4739-a192-4eb2ce7b17fc
function reconstruct(Φr, λr, b, t_max; x1=nothing, Ur=nothing, Wr=nothing)
	if b == nothing
		b = dmd_fit(x1, Ur, Wr, λr)
	end
	Z = zeros(ComplexF64, (size(Φr, 1), t_max))
	for t in 1:t_max
		if t == 1
			Z[:, t] = Φr * b
		else
			Z[:, t] = Φr .* transpose(λr.^(t-1))*b
		end
	end
	return Z
end

# ╔═╡ be181401-7a68-4b96-8b91-46ecd7ef37c4
md"""
# Synthetic input #

For fun, suppose we generate a random initial condition and use the model to calculate the flow. Here is a result. (I have no idea if this picture makes any sense.)
"""

# ╔═╡ 60e450a1-c6d5-4411-8c3e-6e21fa50b5cb
md"""
# Helper functions #
"""

# ╔═╡ d981e923-21fe-4139-9bfa-9f442da486a0
"""
    draw_frame(t, frames; lims=nothing, colormap=:vik, kwargs...)

Renders as a sequence of movie frames.

- `frames[1:m, 1:n, 1:T]` — A sequence of `T` movie frames, each of size `m x n`.
- `t` — The frame to draw, i.e., `frames[:, :, t]`.
- (OPTIONAL) `lims` — Range of values to decide the pixel colors.

Returns a handle to a heatmap object.
"""
function draw_frame(t, frames; lims=nothing, colormap=:vik, kwargs...)
	if lims == nothing
		lims = (minimum(frames), maximum(frames))
	end
	heatmap(frames[:, :, t], c=colormap, clim=lims; kwargs...)
end

# ╔═╡ af700794-ce71-48f9-854b-48c20950cdae
begin
	ALL_mat = matopen("/Users/richie/teaching/cse6730/sp22/src.git/DATA/FLUIDS/CYLINDER_ALL.mat", "r")
	VORTDIMS = (199, 449, 151)
	VORTFRAMES = reshape(read(ALL_mat, "VORTALL"), VORTDIMS)
	VORTLIMS = (minimum(VORTFRAMES), maximum(VORTFRAMES))
	draw_frame(1, VORTFRAMES)
end

# ╔═╡ f762064a-fbc1-4a7e-82f6-4ed4cb0cddd4
begin
	Data = reshape(VORTFRAMES, (VORTDIMS[1]*VORTDIMS[2], VORTDIMS[3]))
	r = 50
	Φr, Wr, λr, Ur, σr, Vr = dmd(Data; r=r)
	σr
end

# ╔═╡ d0cd2f28-f5ec-4f98-a18e-86ff8dad44fd
begin
	scatter(σr, label=nothing, markersize=1.5, yaxis=:log,
	        xlabel=L"r",
		    title=md"``\sigma_r`` (singular values)", titlelocation=:left)
	plot!(σr[1] ./ (1:r), label=L"r^{-1}", linestyle=:dash)
	plot!(σr[1] * exp.(-(0:r-1)./3), label=L"e^{-r/3}", linestyle=:dash)
	plot!(σr[1] * exp.(-(0:r-1)./2), label=L"e^{-r/2}", linestyle=:dash)
	vline!([43], label=nothing, color=:grey, linestyle=:dot)
end

# ╔═╡ 554f13eb-9f26-429c-86ec-83465ba424c2
scatter(real.(λr), imag.(λr), legend=nothing, aspect_ratio=:equal,
	    title="Eigenvalues of \$A_r\$")

# ╔═╡ 7abcdba4-62f5-4675-a915-004f3b1fad41
begin
	Modes = reshape(Φr, (VORTDIMS[1], VORTDIMS[2], size(Φr, 2)))
	k = 1
	mode_k_lims = (minimum(abs.(Modes[:, :, k])), maximum(abs.(Modes[:, :, k])))
	mode_k_real = draw_frame(k, real.(Modes); lims=mode_k_lims,
	                         title="Mode $(k) (real part)")
#	mode_k_imag = draw_frame(k, imag.(Modes); lims=mode_k_lims)
end

# ╔═╡ 9fe39a00-4375-45c6-8575-367eb58f7870
"""
    animate(frames; skip=0, fps=10, filename="anim.gif", kwargs...)

Renders an animation. Uses the `draw_frame` subroutine to render each frame, supplying any keyword arguments via `kwargs...`.

Other options:
- `skip`: Frames to skip in between rendering (0=render every frame, 1=render every other frame, ...)
- `fps`: Speed of the final animation in frames per second
- `filename`: Name of the output file to hold the rendered movie
"""
function animate(frames; skip=0, fps=10, filename="anim.gif", kwargs...)
	step = 1 + minimum([0, skip])
	anim = @animate for t in 1:step:size(frames, 3)
		draw_frame(t, frames; kwargs...)
	end
	gif(anim, filename, fps=fps)
end

# ╔═╡ b63f438e-4741-470d-8aff-ecbe4207abb9
animate(VORTFRAMES; filename="anim_vortices2.gif", lims=VORTLIMS)

# ╔═╡ a45724b7-b413-49f2-b66e-e75e037daf45
begin
	Z = -reconstruct(Φr, λr, nothing, VORTDIMS[3]; x1=Data[:, 1], Ur=Ur, Wr=Wr)
	F = reshape(Z, VORTDIMS)
	animate(real.(F); filename="anim_vortices_approx.gif", lims=VORTLIMS)
end

# ╔═╡ 4959d354-d4e4-4932-821f-5220ff8ad841
begin
	ΔF = real.(F) - VORTFRAMES
	animate(ΔF; filename="anim_vortices_approx_error.gif", # lims=VORTLIMS,
		        title="Approximation error in the reconstruction")
end

# ╔═╡ 55f4fc09-d1f2-4e76-9e77-866e33d883c2
begin
	x1 = (rand(VORTDIMS[1]*VORTDIMS[2]) .- 0.5) .* 100
	Z2 = -reconstruct(Φr, λr, nothing, 25; x1=x1, Ur=Ur, Wr=Wr)
	F2 = reshape(Z2, (VORTDIMS[1], VORTDIMS[2], size(Z2, 2)))
	animate(real.(F2); filename="anim_vortices_approx2.gif") #, lims=VORTLIMS)
end

# ╔═╡ 48659cd7-0589-4fd1-995a-dbf84275a849


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MAT = "23992714-dd62-5051-b70f-ba57cb901cac"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
TSVD = "9449cd9e-2762-5aa3-a617-5413e99d722e"

[compat]
LaTeXStrings = "~1.3.0"
MAT = "~0.10.3"
Plots = "~1.27.6"
TSVD = "~0.4.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.3"
manifest_format = "2.0"
project_hash = "2eabd07635d8c60bf41d1757e99c96da41b9b6ab"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BufferedStreams]]
deps = ["Compat", "Test"]
git-tree-sha1 = "5d55b9486590fdda5905c275bb21ce1f0754020f"
uuid = "e1450e63-4bb3-523b-b2a4-4ffa8c0fd77d"
version = "1.0.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

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
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

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
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

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
git-tree-sha1 = "af237c08bda486b74318c8070adb96efa6952530"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.2"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "cd6efcf9dc746b06709df14e462f0a3fe0786b1e"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.2+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

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

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "Mmap", "Random", "Requires"]
git-tree-sha1 = "cdd249512de03cbf8370365a0a08b9a24955dca9"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.16.6"

[[deps.HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "bab67c0d1c4662d2c4be8c6007751b0b6111de5c"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.1+0"

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

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

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

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

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
git-tree-sha1 = "6f14549f7760d84b2db7a9b10b88cd3cc3025730"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.14"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

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
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

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
git-tree-sha1 = "a970d55c2ad8084ca317a4658ba6ce99b7523571"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.12"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MAT]]
deps = ["BufferedStreams", "CodecZlib", "HDF5", "SparseArrays"]
git-tree-sha1 = "971be550166fe3f604d28715302b58a3f7293160"
uuid = "23992714-dd62-5051-b70f-ba57cb901cac"
version = "0.10.3"

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
version = "2.28.0+0"

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
version = "2022.2.1"

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

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
git-tree-sha1 = "3b429f37de37f1fc603cc1de4a799dc7fbe4c0b6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "6f2dd1cf7a4bbf4f305a0d8750e351cb46dfbe80"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.6"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

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
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

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
version = "0.7.0"

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
git-tree-sha1 = "4f6ec5d99a28e1a749559ef7dd518663c5eca3d5"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "8d7530a38dbd2c397be7ddd01a424e4f411dcc41"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.2"

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
version = "1.0.0"

[[deps.TSVD]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "61cd1ce64b4ffb69e2d156ff7166a8eb796d699a"
uuid = "9449cd9e-2762-5aa3-a617-5413e99d722e"
version = "0.4.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

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
version = "1.2.12+3"

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
version = "5.1.1+0"

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
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

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
# ╟─0f352c21-f8ec-4d60-9f89-acfb061d6350
# ╠═b5204be0-e2ad-4645-8440-ea0e82f9e205
# ╟─5be85952-37d6-4dba-a859-9b4f4e3e1e14
# ╠═af700794-ce71-48f9-854b-48c20950cdae
# ╠═b63f438e-4741-470d-8aff-ecbe4207abb9
# ╟─9996d46a-881d-4350-bd31-5c49df04228c
# ╟─f7c19246-5f95-47e0-9fa5-fea212046ee0
# ╟─3a09cfee-5a6f-411b-9236-32fa02351948
# ╟─295d1dd8-983d-443f-ab6a-28cce3125f37
# ╟─804ab926-f054-4faf-8b00-63e1e8bf6250
# ╟─14a4e646-b775-4b4f-8989-ba0aa3603c30
# ╟─c77885c5-c6d2-4ce7-9b84-d49d7a06028f
# ╟─a6386cb2-e244-44de-b7a6-698c2d458dee
# ╟─af98ac0a-7975-4f1b-886e-ae255a85c36e
# ╟─b35fbda6-473f-4f0b-a90c-74a4c52982b8
# ╟─b9758b50-a41b-4eb1-a5d3-f7b9f0cdf096
# ╟─e7e11aa5-95c1-4569-9d39-2218ab9d6a0a
# ╟─bec98e48-3286-4048-b383-9e661557bc07
# ╟─28f179c8-3796-4b27-a929-8ccec3d9396d
# ╠═2888f19d-1706-4311-a607-e7c38ed8a194
# ╠═00aff0d3-aefa-459b-9b37-2ea648c3dedf
# ╟─9274a087-943e-4c12-96ed-db26e76b3aae
# ╠═f762064a-fbc1-4a7e-82f6-4ed4cb0cddd4
# ╟─ad41f62c-bf60-4e9c-b300-f0ce837d97d2
# ╠═d0cd2f28-f5ec-4f98-a18e-86ff8dad44fd
# ╟─a95ac675-3829-4c77-82a4-808cd939e3fc
# ╟─554f13eb-9f26-429c-86ec-83465ba424c2
# ╟─7b4b73e6-6b8e-42a7-9062-d7d226980e46
# ╠═7abcdba4-62f5-4675-a915-004f3b1fad41
# ╟─47f26dba-c61f-43ad-a5dd-44b4ce0676fb
# ╠═dd616c1d-969a-4739-a192-4eb2ce7b17fc
# ╠═a45724b7-b413-49f2-b66e-e75e037daf45
# ╠═4959d354-d4e4-4932-821f-5220ff8ad841
# ╟─be181401-7a68-4b96-8b91-46ecd7ef37c4
# ╠═55f4fc09-d1f2-4e76-9e77-866e33d883c2
# ╟─60e450a1-c6d5-4411-8c3e-6e21fa50b5cb
# ╠═d981e923-21fe-4139-9bfa-9f442da486a0
# ╟─9fe39a00-4375-45c6-8575-367eb58f7870
# ╠═48659cd7-0589-4fd1-995a-dbf84275a849
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
