function localmatrices(p, h, neighbors, order, pc; scale = false)

    dim   = dimension(pc)
    ncols = length(neighbors)
    nrows = nfuns(order, dim)

    K = zeros(nrows, ncols)
    w = zeros(ncols)

    nrows > ncols && @warn "More conditions than neighbors!"

    cnt = 1
    if dim == 2
        for j = neighbors
            dp = pc.p.x[j] - p

            w[cnt] = exp(-2 * norm(dp) / (h + pc.p.h[j]))

            if scale
                dp *= 1 / h
            end

            dx = dp[1]
            dy = dp[2]

            K[1, cnt] = 1

            if order > 0
                K[2, cnt] = dx
                K[3, cnt] = dy
            end

            if order > 1
                K[4, cnt] = dx^2
                K[5, cnt] = dx * dy
                K[6, cnt] = dy^2
            end

            if order > 2
                K[7, cnt] = dx^3
                K[8, cnt] = dx^2 * dy
                K[9, cnt] = dx * dy^2
                K[10, cnt] = dy^3
            end

            if order > 3
                K[11, cnt] = dx^4
                K[12, cnt] = dx^3 * dy
                K[13, cnt] = dx^2 * dy^2
                K[14, cnt] = dx * dy^3
                K[15, cnt] = dy^4
            end

            cnt += 1
        end
    else
        error("Dimension $dim not implemented!")
    end

    return K, Diagonal(w)
end

function localmatrices(i, order, pc; scale = false)
    return localmatrices(pc.p.x[i], pc.p.h[i], pc.p.neighbors[i], order, pc; scale)
end

nfuns(order, dim) = binomial(order + dim, dim)