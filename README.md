# Assignment10

A simple "tank" model for gas flow from a well controlled by pressure is

$$
\frac{V_p T_s}{P_s T}\frac{d}{dt}\left(\frac{P(t)}{z\left(P(t), T)\right)}\right) + J \left(P(t)^2
- P_{wf}^2\right)^n
$$

where $P(t)$ is the tank pressure, $T$ is the tank temperature (constant),
$V_p$ is the pore volume of the tank, $P_s$ and $T_s$ are standard temperature
and pressure, $J$ is the well productivity index, $n$ is a fitting parameter
that captures inertial effects near the well, $P_{wf}$ is the well
pressure, and $z$ is the gas expansion factor.

Use [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) to solve this
differential equation.  Implement your solution inside the function
`gas_solver()` in the [assignment10.jl](src/assignment10.jl) file.  The
argument list for `gas_solver` should be obvious from the equation above aside
from `Pₒ` which is the initial tank pressure and `tmax` which is the maximum time
you want the solver to run.

The function should return the full solution composite type (the thing that the
`solve` function from DifferentialEquations.jl returns).  A set of arguments to
`gas_solver` that should return a decent solution is available in the
[`runtests.jl`](test/runtests.jl) file.

Just use the default solver, do not specify any additional parameters to
`solve`.  The $z$-factor function has already been coded up for your use.

## Testing

To see if you answers are correct, run the following command at the Terminal
command line from the repository's root directory

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

the tests will run and report if passing or failing.
