using DrWatson
@quickactivate :GFDM

# load point cloud from data directory
# the square domain has its boundaries labeled as "left", "right", "top", and "bottom"
pc = load_csv(datadir("pointclouds", "Square_h=0.025_type=meshfree.csv"))
labelpoints!(pc, Val(:square))

# define heat equation of the form c * dT/dt = div(lambda * grad(T)) + q
# heat conductivity discontinuous
λ(p, t) = prod(cospi, p.x/2) > 3/4 ? 10.0 : 1.0
q(p, t) = pi^2/2 * prod(cospi, p.x/2)
model = HeatEquation(
    heat_capacity       = (p, t) -> 1.0,
    heat_conductivity   = λ,
    heat_source         = q,
    initial_temperature = (p) -> 0.0
)

# define boundary conditions for each boundary segment
bcon = (
    left   = Dirichlet(0.0),
    right  = Dirichlet(0.0),
    top    = Dirichlet(0.0),
    bottom = Dirichlet(0.0)
)

# define method
laplace_wlsq = Laplace_WLSQ(order=2, dd=OneDimensionalCorrection_Default())
diffusion_method = DiffusionSystem(
    # do not smooth diffusivity
    smoothing = Smoothing(),
    # no domain decomposition
    operator = DiffusionOperator_Single(
        # scaled WLSQ-based Lapalce operator with arithmetic averaging of diffusivity
        DivEtaGrad_ScaledLaplacian(:arithmetic_mean, laplace_wlsq)
    )
)

# solve heat equation
solve!(pc, model, bcon;
    method   = diffusion_method,
    timespan = (0.0, 0.8),
    Δt_max   = 1/60,
    filename = datadir("heat_equation_discontinuous","results")  #save result files for paraview
)
