include("LocalOperators/abstract.jl")
include("LocalOperators/approximation.jl")
include("LocalOperators/directionalderivative.jl")
include("LocalOperators/laplace.jl")
include("LocalOperators/smoothing.jl")
include("LocalOperators/div_eta_grad.jl")
include("LocalOperators/neumann.jl")

include("LocalOperators/boundaryconditions.jl")

include("GlobalOperators/diffusion_system.jl")

function default_operators!(pc::PointCloud{dim, Tv, Data}) where {dim, Tv, Data}
    names = coordnames(Val(dim))

    for i = eachindex(pc)
        pc.p.c_gradient[i] = gradient_row(pc.p.x[i], pc.p.h[i], pc.p.neighbors[i], 2, pc, DD_Off())
        pc.p.c_laplace[i]  = isboundary(pc.p[i]) ?
            directional_derivative_row(pc.p.n[i], pc.p.x[i], pc.p.h[i], pc.p.neighbors[i], 1, pc, OneDimensionalCorrection_Default()) :
            laplace_row(i, pc, Laplace_WLSQ())
    end

end