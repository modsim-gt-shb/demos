### A Pluto.jl notebook ###
# v0.18.2

using Markdown
using InteractiveUtils

# ╔═╡ 893b562c-a08f-416a-89b9-29e3647fd4f7
begin
	using LaTeXStrings
	using DataStructures
	using PlutoUI
end

# ╔═╡ 7f27e8ef-e052-4769-88dd-6fea09405d0a
html"<button onclick='present()'>present</button>"

# ╔═╡ 89d8b9ca-3ded-4017-b074-5eb715561452
md"""
# A process-oriented discrete-event simulator #

In a **process-oriented** world view of a discrete-event simulation, we focus on the behavior of each object of the simulated system. The [class notes](https://gatech.instructure.com/courses/233272/files/31041077/download) show one example, specifically for a simple single-runway airport system.

Here is what that example looks like:
"""

# ╔═╡ 7d10cc24-dcb2-41fd-bdbb-edc73dd6e679
LocalResource("./airport-activities.png")

# ╔═╡ 51b8e3ac-c2e6-4227-98fd-2372e2df61aa
md"""
## Airport state variables ##
"""

# ╔═╡ c4f2b0ce-2f4e-4b12-b682-3b15266a47a5
"""
An airport is characterized by behavioral parameters and observable states. The parameters are:

- `R::Float64`: The time that an aircraft spends on the runway while landing
- `G::Float64`: The time that an aircraft spends at the gate after landing but prior to departure

The state variables representing the airport are:

- `runway_free::Ref{Bool}`: A reference to a boolean variable indicating whether the runway is free or not.
- `in_air::Int`: The number of aircraft that are either in the air or occupying the runway during a landing.
- `on_ground::Int`: The number of aircraft occupying a gate prior to departure from the airport.

> We'll explain why `runway_free` is a `Ref{Bool}` rather than simple a `Bool` when we discuss the simulator implementation.

Lastly, we'll augment the primary state variables with some additional "derived observables," used to track statistics about the state variables.

- `in_air_max::Int`: Maximum observed value of `in_air`
- `on_ground_max::Int`: Maximum observed value of `on_ground`
"""
mutable struct AirportState
	# Primary state variables
	runway_free::Ref{Bool}
	in_air::Int
	on_ground::Int
	
	# Parameters
	R::Float64  # Runway time
	G::Float64  # Gate time

	# Derived observables
	in_air_max::Int
	on_ground_max::Int

	# Constructor
    function AirportState(; kwargs...)
        S = new(true, 0, 0, 1.0, 1.0, 0, 0)
        for (key, value) in kwargs
            setfield!(S, key, value)
        end
        return S
    end
end

# ╔═╡ 2b4482cd-b366-47b2-b100-fe8b06377e25
md"""
## Aircraft process ##

When an aircraft arrives at the airport, the simulator will need to create and begin executing a new process for it. The process will interact with the simulator through the primitives `advance_time` and `wait_until`.

Recall the conceptual model of a process, below, and compare it to the process code, implemented in the function `aircraft`.
"""

# ╔═╡ 904595f8-3b00-4623-84ed-673ccf9db0a4
LocalResource("./airport-activities.png")

# ╔═╡ 4541530f-37da-4938-9717-6c75dc2de7cc
md"""
# Airport simulation demo #

Here is a sample execution given the following scenario:
- A runway landing time of `R=3` time units
- An on-ground (at the gate) time of `G=6` time units
- Two initial arrival events, one at time `t=1` and another at time `t=3`.
"""

# ╔═╡ 45d0e179-e250-4e98-a18b-e0b9617de140
LocalResource("./airport-sample-trace.png")

# ╔═╡ a6fd7d7b-8c1c-439a-99a2-ec8340e70425
md"""
In our basic simulation framework, each aircraft is associated with its own simulation process. A simulation is instantiated by

- creating a `simulation` object of type `Simulator`, which is a wrapper around our state variables;
- creating initial events to kickstart each aircraft process;
- and creating a new process for the simulator itself.

In our implementation, each process is implemented as an [asynchronous Julia task](https://docs.julialang.org/en/v1/manual/asynchronous-programming/), including a "main" asynchronous task to run the simulator itself.

> The main task is `t_main`. We are launching it as its own asynchronous task mainly for debugging purposes—since it returns right away, if there are any task synchronization errors that cause the tasks to deadlock, we can still terminate `t_main` interactively.
"""

# ╔═╡ 5654b04d-5aae-4e8f-a299-b0efa3c6fd40
md"""
The following code cell will log the final simulation state.
"""

# ╔═╡ 891289dd-f864-407c-a3d1-e05e4a532d8c
# Uncomment if needed for debugging
begin
	#[terminate(p) for p in simulation.processes if !istaskdone(p)]
	#terminate(t_main)
	#simulation
end

# ╔═╡ b1d87f6c-5dee-4d89-8140-8b3085003307
md"# Simulator infrastructure #"

# ╔═╡ dff2370b-a079-4e4e-8643-a1c0906ef9f5
md"""
The infrastructure for a process-oriented simulator extends our earlier [event-driven framework](https://vuduc.org/cse6730/sp22/discrete-events-from-scratch2.jl.html).
"""

# ╔═╡ 38016d5d-ba91-4757-a39e-8b74fb6f5cfe
md"""
## Unmodified bits ##

First, let's recall some elements of the event-driven framework that can be carried over to the process-oriented framework _without_ modification, namely, `Event` objects and the `FutureEventsList` type.
"""

# ╔═╡ d536d59d-6328-483a-a0c3-62f4afcc6615
"""
An event is a pair consisting of

* `timestamp`, which indicates when the event occurs in logical simulation time; and
* `handler!`, which is a callback function that the simulator should invoke to "execute" the event.

The event handler should have a signature of the form:
```
    function handler!(sim::Simulator)
        ...
    end
```
"""
struct Event
	timestamp
	handler!::Function
end

# ╔═╡ cb8bb339-0bbd-416f-bb1d-8788e0f23fa0
function Base.isless(a::Event, b::Event)
	a.timestamp < b.timestamp
end

# ╔═╡ 35801cdb-715a-47a5-bfbc-54e808f84eb9
# Demo: isless(a::Event, b::Event)
begin
	a = Event(7, (s, fev) -> ())
	b = Event(3, (s, fev) -> ())
	a < b
end

# ╔═╡ 6d14ddce-66ab-4180-b147-5902220aee16
@doc """
A min-heap container for `Event` objects.
"""
FutureEventsList = BinaryMinHeap{Event}

# ╔═╡ 8089a1cf-a18b-4f65-a279-f123b6108f4a
md"""
## Extensions for process management ##

The first extension we need is the addition of new simulation-state variables that can be used to coordinate processes, i.e., cooperatively schedule them.
"""

# ╔═╡ 10c04213-651a-4c0f-aa2b-8a2096459a8a
md"""
## Waiting lists ##

Recall that a process might need to wait until a predicate becomes true. To implement this functionality, we'll use the abstraction of a **wait list.**

A waiting list is a dictionary that maps a predicate reference to a list of all processes waiting on it.

We are using references to boolean variables (`Ref{Bool}`) instead of boolean values (`Bool`) so that multiple asynchronous processes can observe or update the same predicate.
"""

# ╔═╡ 1520cc22-ceb2-47be-9474-0a8c5b5c1a4d
"""
A dictionary that maps predicate variables (as references, `Ref{Bool}`) to a list of waiting process IDs.
"""
WaitList = Dict{Ref{Bool}, Vector{Int}}

# ╔═╡ d2b931d6-d94a-40d9-9859-4dbba62b1f8e
"""
    Simulator(app; now=0, future_events=FutureEventsList(), ...)

The event-driven framework included the following components:

- `app`: The application's model-specific state variables
- `now`: The current value of the global simulation clock
- `future_events`: The future events list

In the process-oriented framework, we need to track all running processes, their run conditions (running or suspended), and any predicates they might be waiting on. For this purpose, let's extend the event-driven simulation state with the following components:

- `processes`: A list of processes (`Task` objects). The process whose ID is `k` is associated with the task `processes[k]`.
- `suspended`: A list of channels used to suspend and resume processes; `suspended[k]` is the channel associated with the process whose ID is `k`
- `resume_scheduler` is a channel used to let the scheduler know it should resume its own execution
- `predicates`: A dictionary mapping each predicate variable to a list of the process IDs currently waiting on it
"""
mutable struct Simulator
	# Event-driven
	now
	future_events::FutureEventsList
	app

	# Additional fields for managing processes
	processes::Vector{Task}
	suspended::Vector{Channel}
	resume_scheduler::Channel
	predicates::WaitList
	
	function Simulator(app=nothing; kwargs...)
		S = new(0, FutureEventsList(), app,
			    Task[], Channel[], Channel(Inf), WaitList())
        for (key, value) in kwargs
            setfield!(S, key, value)
        end
        return S
    end
end

# ╔═╡ ba9bb382-dc2e-475d-addf-05bac66510a8
md"""
## Protocol for suspending and resuming processes ##

Julia's asynchronous tasking mechanism creates concurrently executing tasks. In principle, these are all executing at the same time and do not "stop." Therefore, to emulate the suspending of a process or thread, we need a coordination protocol of some kind.

The protocol implemented here is as follows.
- Each process has an integer ID, `k`, and is associated with the `Task` object, `processes[k]`.
- Each process has a channel, `suspended[k]`. A process waits to run by "taking" a token from its channel, `suspended[k]`.
- The main scheduling task is the first one to run and decides which process gets to execute. If it picks process `k`, then it will
    1. put a token onto channel `suspended[k]`, and
    2. suspend itself by taking (waiting on) a token from `resume_scheduler`.
- If process `k` is running and wishes to suspend, then it puts a token onto `resume_scheduler` and takes (waits on) a token from `suspended[k]`.
"""

# ╔═╡ 02f22473-6828-4637-ba9d-52512164ede3
md"""
## Events ##

We can use the `Event` abstraction of an event-oriented framework to implement two kinds of functionality in a process-oriented framework: creating processes and advancing time.
"""

# ╔═╡ 1d19bcfa-a708-4f1e-a526-c26fabc1e239
md"""
## Creating a process ##
"""

# ╔═╡ eabbab9b-9553-457f-af13-01cf09933730
md"""
The simulation developer should call this first function, `add_process!`, to create an `Event` for kickstarting a new process.

For example, in the airport demo, the developer might create two new aircraft processes to start at times 1.0 and 3.0 via these two calls:
```julia
    add_process!(simulation_state, 1.0, aircraft)
	add_process!(simulation_state, 3.0, aircraft)
```
"""

# ╔═╡ 4dd5ef63-4f1e-41a8-b975-a21033b1941a
md"""
This second function is the handler for new-process events and is _internal_ to the simulator. It creates an asynchronous task and an associated channel so that it can coordinate execution with the main scheduling task.
"""

# ╔═╡ 37e3d5a0-e22d-413a-ab69-1c9efa1f3b4e
"""
    start_process(sim::Simulator; new_process::Function)

This function starts a new process by invoking the user-defined function `new_process` in the context of the simulation state `sim`.

The function `new_process` should have the signature,

```julia
function new_process(app)
    ...
```

where `app` is the user-defined application state.

> The developer's interface hides the _simulator_ state. However, that state is available to each process via Julia's [task-local storage abstraction](https://docs.julialang.org/en/v1/base/parallel/#Base.task_local_storage-Tuple{Any}). This functionality is needed to support the simple interfaces for new processes and both the `advance_time` and `wait_until` primitives.
"""
function start_process(sim::Simulator; new_process::Function)
	id = length(sim.processes) + 1
	task = @task begin
		    task_local_storage("sim", sim)
			task_local_storage("id", id)
			new_process(sim.app)
			put!(sim.resume_scheduler, true)
	end
	push!(sim.processes, task)
	push!(sim.suspended, Channel(Inf))
	bind(sim.suspended[id], sim.processes[id])
	schedule(sim.processes[id])
	take!(sim.resume_scheduler)
end

# ╔═╡ f84f10e6-5116-4205-89b3-7109bf43fb12
"""
    add_process!(sim::Simulator, timestamp, new_process::Function)

Adds a simulation process, `new_process(app)`, to begin at simulation time `timestamp` in the context of the simulation `sim`.
Adds `event` to the future event list of `simulation`.
"""
function add_process!(sim::Simulator, timestamp, new_process::Function)
	handler! = (s) -> start_process(s; new_process=new_process)
	event = Event(timestamp, handler!)
	push!(sim.future_events, event)
end

# ╔═╡ 69921892-3c48-4975-b1e3-c9acc499f884
md"""
## Advancing time (`advance_time`) ##

The developer should use this next function, `advance_time(t)`, to advance the simulation clock by `t` time units. Internally, the function will suspend the process schedule a new event to resume the process at that time.
"""

# ╔═╡ 4a1201a0-5c9d-4f34-886c-8e06d3bed1c5
md"""
The event handler for resuming an event is the following _internal_ function.
"""

# ╔═╡ 3a8c9aba-834d-4c68-ba45-1a52cdb37b0c
"""
    resume_process(id::Int, sim::Simulator)

Resumes a given process and blocks until it suspends again.
"""
function resume_process(id::Int, sim::Simulator)
	put!(sim.suspended[id], true)
	take!(sim.resume_scheduler)
end

# ╔═╡ af30372b-7ce5-4d32-b616-861a53dc35bf
"""
    advance_time(t::Float64)

Suspends the current process until `t` simulation time units have advanced.
"""
function advance_time(t::Float64)
	id = task_local_storage("id")
	sim = task_local_storage("sim")
	resume_handler = (s) -> resume_process(id, s)
	resume_event = Event(sim.now + t, resume_handler)
	push!(sim.future_events, resume_event)
	put!(sim.resume_scheduler, true)
	take!(sim.suspended[id])
end

# ╔═╡ d01a82c6-a315-4bf6-a5b6-2be9b3c1b082
md"""
## Predicates (and `wait_until`) ##

The developer can use the `wait_until(p)` primitive to wait until a particular predicate variable becomes `true`.

The implementation of this primitive tracks the predicate and suspends the process, returning control to the scheduler per the communication protocol explained previously.
"""

# ╔═╡ 3ded01ae-0257-4dae-be87-4e91296a6b7c
md"""
The next two helper functions are internal to the simulator. They allow adding a new predicate to the waiting list and resuming any wait-listed tasks.
"""

# ╔═╡ 4bb42014-9ba9-4c19-af26-2f0df03e00ec
"""
    add_predicate(p::Ref{Bool}, w::WaitList)

Adds process `id` to the waiting list `w` under predicate `p`.
"""
function add_predicate(id::Int, p::Ref{Bool}, w::WaitList)
	if !haskey(w, p)
		w[p] = Int[]
	end
	push!(w[p], id)
end

# ╔═╡ 3572b01d-a56b-4a13-816d-5317bd009ad5
"""
    wait_until(p::Ref{Bool})

Suspends the current process until the predicate `p` is `true`.
"""
function wait_until(p::Ref{Bool})
	id = task_local_storage("id")
	sim = task_local_storage("sim")
	add_predicate(id, p, sim.predicates)
	put!(sim.resume_scheduler, true)
	take!(sim.suspended[id])
end

# ╔═╡ eb27da59-35ab-49c5-8d95-8fca3897ccdd
"""
    aircraft(app::AirportState)

Implements the behavior of an aircraft, which is the following sequence of activities:

1. Arriving at the airport and waiting for the runway to be free.
2. Occupying the runway for `R` time units.
3. Sitting on the ground at a gate for `G` time units.
4. Departing the airport.
"""
function aircraft(app::AirportState)
	# Arrive
	app.in_air += 1
	app.in_air_max = max(app.in_air_max, app.in_air)
	wait_until(app.runway_free)

	# Landing
	app.runway_free[] = false
	advance_time(app.R)
	app.runway_free[] = true

	# Waiting to depart
	app.in_air -= 1
	app.on_ground += 1
	app.on_ground_max = max(app.on_ground_max, app.on_ground)
	advance_time(app.G)

	# Departing
	app.on_ground -= 1
end

# ╔═╡ 09804f46-2e65-4f03-b042-aaf3a6f2f025
"""
    resume_waitlisted(sim)

Resumes any processes waiting on predicates that are now true.
"""
function resume_waitlisted(sim)
	for (p, w) in sim.predicates
		if p[] && length(w) > 0
			# resume process `id`:
			id = popfirst!(w)
			put!(sim.suspended[id], true)
			take!(sim.resume_scheduler)
		end
	end
end

# ╔═╡ ce7783b5-f204-409c-85f5-a93ffa77b1a8
md"""
## Simulation loop ##

With the preceding infrastructure, we are now ready to extend the main simulation loop. This loop is the same as its event-oriented counterpart with the addition of code to resume wait-listed processes.
"""

# ╔═╡ ab83e453-2f27-4067-8b50-34d257251ee4
md"""
## Miscellaneous helper functions ##
"""

# ╔═╡ d907b146-08f5-4ba2-9464-c9a29f3dbb55
"""
    any_events(fev)

Returns `true` if there are any pending events.
"""
any_events(fev::FutureEventsList) = length(fev) > 0

# ╔═╡ 9956c490-a6c6-11ec-28f2-3b4516e927a8
"""
    run!(sim, fev::FutureEventsList; t_max=Inf)

Runs a discrete-event simulation starting from an initial state (`state`) and future event list (`fev`), updating both as each event is executed. Ends when the future event list is empty or the simulation time is greater than or equal to `t_max`. Returns the final simulation time and state.
"""
function run!(sim::Simulator; t_max=Inf)
	while any_events(sim.future_events) && (sim.now < t_max)
		e = pop!(sim.future_events) # next event
		sim.now = e.timestamp # advance clock
		e.handler!(sim) # execute event
		@debug sim # log event	
		resume_waitlisted(sim)
	end
	sim.now
end

# ╔═╡ 4fdefcb2-8d8e-4dba-bfc5-5bdcb4d87053
begin
	simulation = Simulator(AirportState(R=3.0, G=6.0))
	add_process!(simulation, 1.0, aircraft)
	add_process!(simulation, 3.0, aircraft)
	t_main = @async run!(simulation; t_max=15)
end

# ╔═╡ b34c62e5-db00-4a8f-9df9-6fe9b1e177c6
begin
	timedwait(() -> istaskdone(t_main), 10)
	@info simulation
	simulation.now, istaskdone(t_main)
end

# ╔═╡ fda04866-5946-4e3a-a241-9569bbbab1d7
"""
    terminate(task::Task)

Terminates a task (by sending an interrupt exception to it).
"""
function terminate(task::Task)
	if !istaskdone(task)
		@async Base.throwto(task, InterruptException())
	end
end

# ╔═╡ b20d7824-13ef-4a53-a685-76c99a1c3a5f
md"""
# Postscript: Coroutines in other programming languages #

If you are looking for comparable functionality in other languages, you will discover a variety of options.

The process-oriented implementation described in this notebook does not exploit full concurrency. It effectively implements [coroutines](https://en.wikipedia.org/wiki/Coroutine) (or cooperative multitasking, or non-prëemptive multitasking): in reality, only one thread is active at a time, and each thread schedules the next one. Julia does not support this style explicitly, but other languages do:
- In Python, you can use the `yield` and generator constructs, e.g., see [this tutorial on generators and `yield`](https://realpython.com/introduction-to-python-generators/) or [PEP 342](https://peps.python.org/pep-0342/).
- In C, there are the `setjmp` and `longjmp` functions, which can be used as demonstrated [here](https://fanf.livejournal.com/105413.html?).
- C++20 supports tasks, generators, and yield constructs: [C++ coroutines](https://en.cppreference.com/w/cpp/language/coroutines)

Check the [wiki page on coroutines](https://en.wikipedia.org/wiki/Coroutine) for examples in other languages.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataStructures = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DataStructures = "~0.18.11"
LaTeXStrings = "~1.3.0"
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

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

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

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

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

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

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

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

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
# ╟─7f27e8ef-e052-4769-88dd-6fea09405d0a
# ╠═893b562c-a08f-416a-89b9-29e3647fd4f7
# ╟─89d8b9ca-3ded-4017-b074-5eb715561452
# ╟─7d10cc24-dcb2-41fd-bdbb-edc73dd6e679
# ╟─51b8e3ac-c2e6-4227-98fd-2372e2df61aa
# ╠═c4f2b0ce-2f4e-4b12-b682-3b15266a47a5
# ╟─2b4482cd-b366-47b2-b100-fe8b06377e25
# ╟─904595f8-3b00-4623-84ed-673ccf9db0a4
# ╠═eb27da59-35ab-49c5-8d95-8fca3897ccdd
# ╟─4541530f-37da-4938-9717-6c75dc2de7cc
# ╟─45d0e179-e250-4e98-a18b-e0b9617de140
# ╟─a6fd7d7b-8c1c-439a-99a2-ec8340e70425
# ╠═4fdefcb2-8d8e-4dba-bfc5-5bdcb4d87053
# ╟─5654b04d-5aae-4e8f-a299-b0efa3c6fd40
# ╠═b34c62e5-db00-4a8f-9df9-6fe9b1e177c6
# ╠═891289dd-f864-407c-a3d1-e05e4a532d8c
# ╟─b1d87f6c-5dee-4d89-8140-8b3085003307
# ╟─dff2370b-a079-4e4e-8643-a1c0906ef9f5
# ╟─38016d5d-ba91-4757-a39e-8b74fb6f5cfe
# ╠═d536d59d-6328-483a-a0c3-62f4afcc6615
# ╠═cb8bb339-0bbd-416f-bb1d-8788e0f23fa0
# ╠═35801cdb-715a-47a5-bfbc-54e808f84eb9
# ╠═6d14ddce-66ab-4180-b147-5902220aee16
# ╟─8089a1cf-a18b-4f65-a279-f123b6108f4a
# ╠═d2b931d6-d94a-40d9-9859-4dbba62b1f8e
# ╟─10c04213-651a-4c0f-aa2b-8a2096459a8a
# ╠═1520cc22-ceb2-47be-9474-0a8c5b5c1a4d
# ╟─ba9bb382-dc2e-475d-addf-05bac66510a8
# ╟─02f22473-6828-4637-ba9d-52512164ede3
# ╟─1d19bcfa-a708-4f1e-a526-c26fabc1e239
# ╟─eabbab9b-9553-457f-af13-01cf09933730
# ╠═f84f10e6-5116-4205-89b3-7109bf43fb12
# ╟─4dd5ef63-4f1e-41a8-b975-a21033b1941a
# ╠═37e3d5a0-e22d-413a-ab69-1c9efa1f3b4e
# ╟─69921892-3c48-4975-b1e3-c9acc499f884
# ╠═af30372b-7ce5-4d32-b616-861a53dc35bf
# ╟─4a1201a0-5c9d-4f34-886c-8e06d3bed1c5
# ╠═3a8c9aba-834d-4c68-ba45-1a52cdb37b0c
# ╟─d01a82c6-a315-4bf6-a5b6-2be9b3c1b082
# ╠═3572b01d-a56b-4a13-816d-5317bd009ad5
# ╟─3ded01ae-0257-4dae-be87-4e91296a6b7c
# ╠═4bb42014-9ba9-4c19-af26-2f0df03e00ec
# ╠═09804f46-2e65-4f03-b042-aaf3a6f2f025
# ╟─ce7783b5-f204-409c-85f5-a93ffa77b1a8
# ╠═9956c490-a6c6-11ec-28f2-3b4516e927a8
# ╟─ab83e453-2f27-4067-8b50-34d257251ee4
# ╠═d907b146-08f5-4ba2-9464-c9a29f3dbb55
# ╠═fda04866-5946-4e3a-a241-9569bbbab1d7
# ╟─b20d7824-13ef-4a53-a685-76c99a1c3a5f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
