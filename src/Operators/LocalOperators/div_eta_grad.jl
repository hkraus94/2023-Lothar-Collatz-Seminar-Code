abstract type AbstractDiffusionOperator end

struct DivEtaGrad_WLSQ{DD <: AbstractDiagonalDominance} <: AbstractDiffusionOperator
    order::Int
    scaling::Bool
    dd::DD
end

struct DivEtaGrad_ScaledLaplacian{LAP <: AbstractLaplaceOperator} <: AbstractDiffusionOperator
    etaMean::Symbol
    laplacian::LAP
end

"""
    Computes strong form div eta grad operator using a point cloud.
    scaling: \n
    0 -> no scaling
    1 -> scaling μ = log(η) with exp(μ) rescaling
"""
function div_eta_grad_row(i, pc, method::DivEtaGrad_WLSQ)
    neighbors = pc.p.neighbors[i]
    cgrad = pc.p.c_gradient[i]
    eta  = pc.p.eta[neighbors]
    mue  = pc.p.mue[neighbors]

    rhs  = zeros(nfuns(method.order, 2))
    K, W = localmatrices(pc.p.x[i], pc.p.h[i], neighbors, method.order, pc)

    if method.scaling == false
        rhs[2:3] = cgrad * eta
        rhs[4] = 2eta[1]
        rhs[6] = 2eta[1]
    elseif method.scaling == true
        rhs[2:3] = cgrad * mue
        rhs[4] = 2.0
        rhs[6] = 2.0
    end

    ρ = if method.scaling == false
        1.0
    elseif method.scaling == true
        exp(mue[1])
    end
    return ρ * leastsquares(K, W, rhs, method.dd)

end

# function div_eta_grad_row(i, pc, method::DivEtaGradRow{:flux_conservation})
#     neighbors = pc.p.neighbors[i]

#     cx = pc.operators[:x][i, neighbors]
#     cy = pc.operators[:y][i, neighbors]

#     eta = pc.p.eta[neighbors]

#     dir = normalize(SVector(dot(cx, eta), dot(cy, eta)))

#     ptsL = Int[]
#     ptsR = Int[]

#     ptsLloc = Int[]
#     ptsRloc = Int[]

#     for (k,j) = enumerate(neighbors)
#         dv = pc.p.x[j] - pc.p.x[i]
#         dp = dot(dv, dir)
#         dp <= 0 && (push!(ptsL, j), push!(ptsLloc, k))
#         dp >= 0 && (push!(ptsR, j), push!(ptsRloc, k))
#     end

#     etaLs = eta[ptsLloc]
#     etaRs = eta[ptsRloc]

#     nL = length(ptsL)
#     nR = length(ptsR)

#     etaL = exp(sum(log, etaLs) ./ nL)
#     etaR = exp(sum(log, etaRs) ./ nR)

#     # println("i = $i: etaL = $etaL, etaR = $etaR")

#     cxL = directional_derivative_row(dir, pc.p.x[i], pc.p.h[i], ptsL, method.order, pc, method.dd)
#     cxR = directional_derivative_row(dir, pc.p.x[i], pc.p.h[i], ptsR, method.order, pc, method.dd)

#     v = zeros(precision(pc), length(neighbors))

#     v[ptsLloc] += etaL * cxL
#     v[ptsRloc] -= etaR * cxR

#     return v[1] < zero(v[1]) ? v : -v
# end

# function div_eta_grad_row(i, pc, ::DivEtaGradRow{:weak_voronoi})
#     edgelength, area = voronoi(i, pc)
#     neighbors = pc.p.neighbors[i]
#     eta = pc.p.eta[neighbors]

#     v = zeros(precision(pc), length(neighbors))

#     for (k, j) = enumerate(neighbors)
#         edgelength[k] == 0.0 && continue
#         # etaij = (eta[1] + eta[k]) / 2
#         etaij = 2 * (eta[1] * eta[k]) / (eta[1] + eta[k])
#         # etaij = sqrt(eta[1] * eta[k])
#         W = etaij * ( edgelength[k] / norm(pc.p.x[i] - pc.p.x[j]) )
#         v[1] -= W
#         v[k] += W
#     end

#     return v / area
# end

eta_mean(i, pc, ::Val{:minimum}) = map(x -> min(pc.p.eta[i], x), pc.p.eta[pc.p.neighbors[i]])
eta_mean(i, pc, ::Val{:maximum}) = map(x -> max(pc.p.eta[i], x), pc.p.eta[pc.p.neighbors[i]])
eta_mean(i, pc, ::Val{:arithmetic_mean}) = map(x -> (x + pc.p.eta[i]) / 2, pc.p.eta[pc.p.neighbors[i]])
eta_mean(i, pc, ::Val{:harmonic_mean})   = map(x -> 2 / (1 / pc.p.eta[i] + 1 / x), pc.p.eta[pc.p.neighbors[i]])
eta_mean(i, pc, ::Val{:geometric_mean})  = map(x -> sqrt(x * pc.p.eta[i]), pc.p.eta[pc.p.neighbors[i]])
eta_mean(i, pc, ::Val{:center})  = pc.p.eta[i]

function gradient_weno(i, U, pc; ϵ = 1e-5, r = 2)
    ∇U = SVector{dimension(pc), GFDM.precision(pc)}[]
    neighbors = pc.p.neighbors[i]
    xi = pc.p.x[i]
    Ui = U[i]
    A = zero(MMatrix{2, 2, GFDM.precision(pc)})
    b = zero(MVector{2,    GFDM.precision(pc)})
    for j = 2:length(neighbors)
        for k = j+1:length(neighbors)
            xj = pc.p.x[neighbors[j]]
            xk = pc.p.x[neighbors[k]]

            A[1,:] = xj - xi
            A[2,:] = xk - xi

            det(A) == 0 && continue

            b[1] = U[neighbors[j]] - Ui
            b[2] = U[neighbors[k]] - Ui

            ∇Ujk = A \ b

            push!(∇U, ∇Ujk)
        end
    end

    weights = (norm.(∇U) .+ ϵ).^(-r)
    sumweights = sum(weights)

    w = weights / sumweights

    return sum(w .* ∇U)
end

function eta_mean(i, pc, ::Val{:gradient_reconstruction})
    neighbors = pc.p.neighbors[i]

    # gradients = map(i -> gradient_weno(i, pc.p.eta, pc), neighbors)
    gradients = map(j -> pc.p.c_gradient[j] * pc.p.eta[pc.p.neighbors[j]], neighbors)
    DX = pc.p.x[neighbors] .- (pc.p.x[i], )

    return map((k, j) -> (pc.p.eta[i] + pc.p.eta[j]) / 2 + 1/8 * dot(gradients[1] - gradients[k], DX[k]), 1:length(neighbors), neighbors)
end

function div_eta_grad_row(i, pc, method::DivEtaGrad_ScaledLaplacian)
    etamean = eta_mean(i, pc, Val(method.etaMean))
    clap    = laplace_row(i, pc, method.laplacian)

    cdeg    = etamean .* clap
    cdeg[1] = -sum(cdeg[2:end])

    return cdeg
end

# function div_eta_grad_row(i, pc, ::DivEtaGradRow{:weak_centroid})
#     neighbors   = pc.p.neighbors[i]
#     n_neighbors = length(neighbors)

#     eta = pc.p.eta[neighbors]

#     area = pc.p.dV[i]
#     delaunay = pc.cells[pc.p.simplices[i]]

#     edgelength = zeros(precision(pc), n_neighbors)

#     for tri = delaunay
#         centroid = sum(pc.p.x[[tri]]) / 3
#         edge_out = setdiff(tri, i)

#         for j = edge_out
#             pmid = (pc.p.x[i]+ pc.p.x[j]) / 2

#             k = findfirst(==(j), neighbors)
#             edgelength[k] += norm(pmid - centroid)
#         end
#     end

#     v = zeros(precision(pc), n_neighbors)

#     for k = eachindex(neighbors)
#         edgelength[k] == 0.0 && continue
#         j = neighbors[k]
#         etaij = (eta[1] + eta[k]) / 2
#         W = (etaij * edgelength[k]) / norm(pc.p.x[i] - pc.p.x[j])
#         v[1] -= W
#         v[k] += W
#     end

#     return v / area
# end

# function div_eta_grad_row(i, pc, method::DivEtaGradRow{:gfvm})
#     edgelength, area = voronoi(i, pc)
#     neighbors = pc.p.neighbors[i]
#     eta = pc.p.eta[neighbors]

#     c = zeros(length(neighbors))

#     for (k, j) = enumerate(neighbors)
#         (edgelength[k] ≈ 0.0 || j == k) && continue
#         etaij = 2 * (eta[1] * eta[k]) / (eta[1] + eta[k])

#         nij = normalize(pc.p.x[j] - pc.p.x[i])
#         xij = (pc.p.x[i] + pc.p.x[j]) / 2
#         γik = directional_derivative_row(nij, xij, pc.p.h[i], pc.p.neighbors[i], method.order, pc, method.dd)

#         if γik[1] >= 0 || any(<(0), γik[2:end])
#             error("Found non diagonally dominant γik for point $i and neighbor $j")
#         end

#         c += edgelength[k] * etaij * γik
#     end

#     if area ≈ 0
#         error("Area for point $i almost zero, area = $area")
#     end

#     return c / area
# end