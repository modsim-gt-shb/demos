### A Pluto.jl notebook ###
# v0.18.2

using Markdown
using InteractiveUtils

# ╔═╡ 67eaf394-8ea7-4bc0-b006-1d6124c51e61
using PlutoUI

# ╔═╡ a57ff680-56f6-496e-9249-ce2dc3776d63
html"<button onclick='present()'>present</button>"

# ╔═╡ 5d9e752d-ae7a-4f84-ae1b-c4535e6afc63
md"""
# Tasking in Julia #

To implement processes for our process-oriented discrete-event simulator, we will use [Julia's asynchronous tasking](https://docs.julialang.org/en/v1/manual/asynchronous-programming/) module. Here is a quick tour of how it works.
"""

# ╔═╡ 90ad4c40-cdcb-4e3b-8f76-e5fd097d0c1f
md"""
# Goal: Suspending and resuming processes #

To support a process-oriented world view, we need a way to run a process, suspend it, and later resume it. Consider this cartoon of such a process, which has a unique identifier (`my_id`) and a predicate-like object on which it wants to suspend:
"""

# ╔═╡ 9410a588-45b7-405e-bea8-f625815341c4
md"""
## Demo: `new_process` without suspending ##

Here is a sample run, where the "predicate" `c` assumes its default value, `nothing`, which we'll take to mean, "immediately resume."
"""

# ╔═╡ dfe7916f-87dd-4cb4-a767-e10207170a26
md"""
# Tasking basics #

The first concept we need is that of a **task**, which is an _asynchronous unit of execution_.

- Unit of execution: A task is typically a function (method) call, like an invocation of `new_process`. However, it can also be an arbitrary block of code (`begin ... end`).
- Asynchronous: A task may run concurrently with other tasks.

Tasks will enable us to achieve a picture like the one below from [our class notes on the process-oriented world view](https://gatech.instructure.com/courses/233272/files/31041077/download).
"""

# ╔═╡ 3da381b3-9dd2-46c6-bc03-93e9bfefec61
LocalResource("./process-oriented-cartoon.png")

# ╔═╡ 78f3fae9-6a99-45b6-855c-29661142e78d
md"""
## The lifecycle of a `Task` ##

The "life" of a task has three "acts."

1. Creating the task using the `@task` macro, which returns a handle as a `Task` object
2. Declaring the task ready for execution, also called _scheduling_ the task, using `schedule`
3. Waiting for the task to complete, using `wait`

Here is an example.
"""

# ╔═╡ 1875ead1-d815-4849-91c5-74ceb7eea378
function nap(name, duration)
	@info "Hi, I'm $(name), and I'm going to nap for $(duration) seconds..."
	sleep(duration)
	@info "It's $(name) again, and I'm awake!"
end

# ╔═╡ 6b7f97c5-96ae-45ae-aa7a-6ad4064c00c1
sleepy_task = @task nap("Sally", 5)

# ╔═╡ 9a86a1e3-a39c-42df-8415-045ee372e197
schedule(sleepy_task)

# ╔═╡ e2ee6ddc-031d-438d-b002-3dbe82efd2e2
wait(sleepy_task)

# ╔═╡ cc66fa69-9081-4580-8bc3-3263cb9e3ce1
md"""
## Multiple tasks ##

Of course, the most important property of tasking is concurrency!
"""

# ╔═╡ 78dd4816-2a13-43bd-b14c-f48d478b721b
begin
	nappers = Dict{String, Task}()
	people = ["Alice", "Anne", "Ann", "Alfredo"]
end

# ╔═╡ 47903565-8013-4adf-81d8-8259890b8575
md"""
In this case, let's use the `@async` macro to create and schedule the task in one shot, for each of three tasks.
"""

# ╔═╡ 214ba8cc-a62b-4387-b24b-e90092431340
for person in people
	nappers[person] = @async nap(person, length(person))
end

# ╔═╡ 9510d7e0-296c-4ac2-a47e-a06fdded9a46
sort(people, by=length, rev=true)

# ╔═╡ d51a3bb4-0139-484e-9671-e215e13bb015
for person in sort(people, by=length, rev=true)
	@debug person
	wait(nappers[person])
end

# ╔═╡ 1cbc333e-6ccf-494f-8c94-cc33177c57f1
md"""
# Channels: A abstraction for producer-consumer communication #

One way that tasks can talk to each other is through **channels** (`Channel` objects). 

A channel expresses a producer-consumer pattern: one task places an object into the channel while another task waits for an object to appear on the channel, subsequently extracting it.
"""

# ╔═╡ 09f9e766-eec1-4ab0-a055-da45e2d9062c
md"""
## `wait_until` using channels ##

For example, consider this implementation of `wait_until`, which simply takes an object from the channel.
"""

# ╔═╡ ff11ee92-7139-4f11-9192-c1b02799e36b
function wait_until(c::Union{Nothing,Channel})
	if c != nothing
		while !take!(c)
			; # "spin loop"
		end
	end
end

# ╔═╡ 03afdf8a-a6c1-11ec-3693-cf8c233309db
function new_process(my_id; c=nothing)
	@info md"# [$(my_id)] #"
	@info md"[$(my_id)] 1. work work work..."
	
	wait_until(c)
	
	@info md"[$(my_id)] 2. work work work..."
	
	wait_until(c)
	
	@info md"[$(my_id)] 3. work work work..."
	
	wait_until(c)
	
	@info md"# [$(my_id)] done! #"
end

# ╔═╡ 628e2281-2e45-4783-953e-bdd3645f46dc
new_process("A")

# ╔═╡ 4d9732f2-fd28-4af0-8ef0-9ba0a7c00773
md"""
## Demo: Channels in action ##

Let's now create a `new_process` task, associating it with a channel.
"""

# ╔═╡ 10ecdb47-d3b2-4622-b7a5-1efb481da84a
begin
	c0 = Channel(Inf)  # `Inf` indicates the capacity of this channel (infinite)
	t0 = @async new_process(0; c=c0) # Start a task
end

# ╔═╡ 0df18b65-1bce-4c15-a122-4be77b0f1fce
md"""
Now repeatedly run the cell below. When you want the task to proceed from the `wait_until`, change the "token" placed in the channel to `true`.
"""

# ╔═╡ e9fe56ec-0e0c-45a8-bf7e-d624e8a67ef8
istaskdone(t0), put!(c0, false)

# ╔═╡ 26b7633f-9a0b-4911-a0af-b2bc9876fc59
md"""
## Variation on a theme ##

You can also bind channels to tasks, such that the channel closes automatically when the task has finished.
"""

# ╔═╡ 21e8f274-ad47-4fc5-9b4f-f161bde458df
begin
	t1 = Ref{Task}()
	f1(c) = new_process(1; c)
	c1 = Channel(f1, Inf; taskref=t1, spawn=true)
end

# ╔═╡ 9730487f-680b-46f8-a02d-6a795cc904b2
md"""
> Try changing logical-or, `||`, to logical-and, `&&`, to see what happens when you try to put an object onto `c1` _after_ task `t1[]` is complete.
"""

# ╔═╡ 5bae1b7b-6add-41c1-856c-6eb2fe4ec0e3
istaskdone(t1[]) || put!(c1, false)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.37"
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

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "bf0a1121af131d9974241ba53f601211e9303a9e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.37"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─a57ff680-56f6-496e-9249-ce2dc3776d63
# ╠═67eaf394-8ea7-4bc0-b006-1d6124c51e61
# ╟─5d9e752d-ae7a-4f84-ae1b-c4535e6afc63
# ╟─90ad4c40-cdcb-4e3b-8f76-e5fd097d0c1f
# ╠═03afdf8a-a6c1-11ec-3693-cf8c233309db
# ╟─9410a588-45b7-405e-bea8-f625815341c4
# ╠═628e2281-2e45-4783-953e-bdd3645f46dc
# ╟─dfe7916f-87dd-4cb4-a767-e10207170a26
# ╟─3da381b3-9dd2-46c6-bc03-93e9bfefec61
# ╟─78f3fae9-6a99-45b6-855c-29661142e78d
# ╠═1875ead1-d815-4849-91c5-74ceb7eea378
# ╠═6b7f97c5-96ae-45ae-aa7a-6ad4064c00c1
# ╠═9a86a1e3-a39c-42df-8415-045ee372e197
# ╠═e2ee6ddc-031d-438d-b002-3dbe82efd2e2
# ╟─cc66fa69-9081-4580-8bc3-3263cb9e3ce1
# ╠═78dd4816-2a13-43bd-b14c-f48d478b721b
# ╟─47903565-8013-4adf-81d8-8259890b8575
# ╠═214ba8cc-a62b-4387-b24b-e90092431340
# ╠═9510d7e0-296c-4ac2-a47e-a06fdded9a46
# ╠═d51a3bb4-0139-484e-9671-e215e13bb015
# ╟─1cbc333e-6ccf-494f-8c94-cc33177c57f1
# ╟─09f9e766-eec1-4ab0-a055-da45e2d9062c
# ╠═ff11ee92-7139-4f11-9192-c1b02799e36b
# ╟─4d9732f2-fd28-4af0-8ef0-9ba0a7c00773
# ╠═10ecdb47-d3b2-4622-b7a5-1efb481da84a
# ╟─0df18b65-1bce-4c15-a122-4be77b0f1fce
# ╠═e9fe56ec-0e0c-45a8-bf7e-d624e8a67ef8
# ╟─26b7633f-9a0b-4911-a0af-b2bc9876fc59
# ╠═21e8f274-ad47-4fc5-9b4f-f161bde458df
# ╟─9730487f-680b-46f8-a02d-6a795cc904b2
# ╠═5bae1b7b-6add-41c1-856c-6eb2fe4ec0e3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
