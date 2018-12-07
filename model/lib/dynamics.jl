function bdiffx(fn, nx, ny, dx)
    # backward difference
    fnp1 = zeros((ny, nx))
    
    for j in 2:(nx-1)
        fnp1[:, j] = (fn[:, j] - fn[:, j-1]) / dx
    end

    return fnp1
end

function fdiffx(fn, nx, ny, dx)
    # backward difference
    fnp1 = zeros((ny, nx))
    
    for j in 1:nx
        fnp1[:, j] = (fn[:, j+1] - fn[:, j]) / dx
    end

    return fnp1
end

function bdiffy(fn, nx, ny, dy)
    # backward difference
    fnp1 = zeros((ny, nx))
    
    for i in 2:(ny-1)
        fnp1[i, :] = (fn[i, :] - fn[i-1, :]) / dy
    end

    return fnp1
end

function fdiffy(fn, nx, ny, dy)
    # backward difference
    fnp1 = zeros((ny, nx))
    
    for i in 1:ny
        fnp1[i, :] = (fn[i+1, :] - fn[i, :]) / dy
    end

    return fnp1
end

function ucoriolis(lats, yy, vn, nx, ny)
    omega = 2 * pi / 86400.
    R = 6.371e6
    beta = zeros((ny+1, nx))
    unp1 = zeros((ny, nx+1))

    for i in 1:(ny+1)
        for j in 1:nx
            beta[i, j] = 2 * omega / R
        end
    end

    f = beta .* yy

    for i in 1:ny
        for j in 2:nx
            unp1[i, j] = (f[i, j-1] * vn[i, j-1] + f[i, j] * vn[i, j] + 
                f[i+1, j-1] * vn[i+1, j-1] + f[i+1, j] * vn[i+1, j]) / 4. 
        end
    end

    return unp1
end

function vcoriolis(lats, yy, un, nx, ny)
    omega = 2 * pi / 86400.
    R = 6.371e6
    beta = zeros((ny, nx+1))
    vnp1 = zeros((ny+1, nx))

    for i in 1:ny
        for j in 1:(nx+1)
            beta[i, j] = 2 * omega / R
        end
    end

    f = beta .* yy

    for i in 2:ny
        for j in 1:nx
            vnp1[i, j] = (f[i-1, j] * un[i-1, j] + f[i-1, j+1] * un[i-1, j+1] + 
                f[i, j] * un[i, j] + f[i, j+1] * un[i, j+1]) / 4. 
        end
    end

    return vnp1
end