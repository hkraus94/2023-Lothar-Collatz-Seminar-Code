"""
Sets up poisson equation \n
    -∇⋅(λ∇T) = q.
Dirichlet and neumann boundary conditions can be used.
"""
struct PoissonEquation{Λ, Q} <: AbstractPartialDifferentialEquation
    λ::Λ
    q::Q
end


function prepare_timestep!(
    pc::PointCloud;

    compute_triangulation = true,
    checks = false, αtol = deg2rad(15),
    compute_interfaces = true,
    compute_neighborhoods = true,
    compute_sparsitypattern = true,
    compute_default_operators = true,
    kwargs...
)
    if compute_neighborhoods
        neighbors!(pc)
    end

    if compute_triangulation
        delaunay!(pc; checks, αtol)
        volumes!(pc; kwargs...)
    end

    if compute_interfaces && isseparated(pc)
        interfaces!(pc)
        neighbors!(pc)
    end

    if compute_sparsitypattern
        pc.sparsity = COOPattern(pc.p.neighbors)
    end

    if compute_default_operators
        default_operators!(pc)
    end
end

function setvars!(pc::PointCloud, model::PoissonEquation)
    pc.p.λ .= model.λ.(particles(pc))
    pc.p.q .= model.q.(particles(pc))
end

function solve!(
    pc::PointCloud{2},
    model::PoissonEquation,
    boundaryconditions;
    method = DiffusionSystem(
        Smoothing(w -> exp(-3w), 1),
        DiffusionOperator_Single(
            DivEtaGrad_WLSQ(2, true, OneDimensionalCorrection_Default())
        )
    ),
    prepared_timestep = false
)
    check_consistency_boundary(pc, boundaryconditions)

    # bfunT = extract_boundary(boundaryconditions, :T)

    prepared_timestep || prepare_timestep!(pc)

    setvars!(pc, model)

    M = zero(pc.sparsity)
    B = zero(pc.sparsity)

    rm = zeros(length(pc))
    rb = zeros(length(pc))

    assemble_diffusion_system!(M, rm, B, rb, model.λ, model.q, boundaryconditions, method, pc)

    A = B - M
    r = rm + rb

    Pl = Diagonal(ones(length(pc)))

    for i = eachindex(pc)
        # λ = 1 / A[i,i]
        # λ = 1 / norm(view(A, i, pc.p.neighbors[i]))
        λ = 1 / maximum(view(A, i, pc.p.neighbors[i]))
        Pl[i,i] = λ

        if A[i,i] < 0
            A[i, :] *= -1
            r[i] *= -1
        end
    end

    A = Pl * A
    r = Pl * r

    LU = IncompleteLU.ilu(A)

    (pc.p.T[:], log) = IterativeSolvers.bicgstabl(A, r, Pl = LU, verbose = false, log = true, abstol = 1e-14, reltol = 1e-14)
    println(log)
    println("Residual = $(norm(A * pc.p.T - r))")
    return log

    # pc.p.T[:] = A \ r
    # println("Residual = $(norm(A * pc.p.T - r))")
    # return (isconverged = true, iters = 1)
end