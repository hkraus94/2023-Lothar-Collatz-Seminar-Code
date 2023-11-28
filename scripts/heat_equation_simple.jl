using DrWatson
@quickactivate :GFDM

# load point cloud from data directory
# the square domain has its boundaries labeled as "left", "right", "top", and "bottom"
pc = load_csv(datadir("pointclouds", "Square_h=0.025_type=meshfree.csv"))

# define heat equation of the form c * dT/dt = div(lambda * grad(T)) + q
# each function depends on the point p (see Particle) and time t
model = HeatEquation(
    heat_capacity       = (p, t) -> 1.0,
    heat_conductivity   = (p, t) -> 0.5,
    heat_source         = (p, t) -> 0.0,
    # initial temperature does not depend on time, T(x, y) = cos(π/2*x) * cos(π/2*y)
    initial_temperature = (p) -> prod(cospi, p.x/2)
)

# define boundary conditions for each boundary segment
bcon = (
    left   = Dirichlet(0.0),
    right  = Dirichlet(0.0),
    top    = Dirichlet(0.0),
    bottom = Dirichlet(0.0)
)

# define method
diffusion_method = DiffusionSystem(
    # do not smooth diffusivity
    smoothing = Smoothing(),
    # no domain decomposition
    operator = DiffusionOperator_Single(
        # weighted least squares optimization for each point with a
        # one-dimensional diagonal dominance correction
        # to turn off correction replace OneDimensionalCorrection_Default with DD_Off()
        DivEtaGrad_WLSQ(2, false, OneDimensionalCorrection_Default())
    )
)

# solve heat equation
solve!(pc, model, bcon;
    method   = diffusion_method,
    timespan = (0.0, 2.0),
    Δt_max   = 1/60,
    filename = datadir("heat_equation_simple","results")  #save result files for paraview
)
