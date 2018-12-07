include("./dynamics.jl")

function u_rhs(parms, un, phin, vn; bc=1)
    vg, ulats, vlats, uyy, vyy, nx, ny, dx, dy, dt = parms
    rhs = - bdiffx(phin, nx+1, ny, dx) + ucoriolis(vlats, vyy, vn, nx, ny)

    # x - radiational boundary conditions
    rhs[:, nx+1] = - bc * vg * (un[:, nx+1] - un[:, nx]) / dx
    rhs[:, 1] = + bc * vg * (un[:, 2] - un[:, 1]) / dx    

    return rhs
end

function phi_rhs(parms, phin, un, vn, fphi)
    vg, ulats, vlats, uyy, vyy, nx, ny, dx, dy, dt = parms
    rhs = - (vg * vg) .* (fdiffx(un, nx, ny, dx) + fdiffy(vn, nx, ny, dy)) .+ fphi

    return rhs
end

function v_rhs(parms, vn, phin, un; bc=1)
    vg, ulats, vlats, uyy, vyy, nx, ny, dx, dy, dt = parms
    rhs = - bdiffy(phin, nx, ny+1, dy) - vcoriolis(ulats, uyy, un, nx, ny)

    # y - radiational boundary conditions
    rhs[ny+1, :] = - bc * vg * (vn[ny+1, :] - vn[ny, :]) / dy
    rhs[1, :] = + bc * vg * (vn[2, :] - vn[1, :]) / dy

    return rhs
end

function rk4(parms, dt, un, phin, vn, fphi)
    vg, ulats, vlats, uyy, vyy, nx, ny, dx, dy, dt = parms
    scaler = [0.5, 0.5, 1]

    ku = zeros((4, ny, nx+1))
    kv = zeros((4, ny+1, nx))
    kphi = zeros((4, ny, nx))

    uns = un * 1.
    vns = vn * 1.
    phins = phin * 1.

    for i in 1:3
        ku[i, :, :] = u_rhs(parms, uns, phins, vns)
        kv[i, :, :] = v_rhs(parms, vns, phins, uns)
        kphi[i, :, :] = phi_rhs(parms, phins, uns, vns, fphi)

        uns = un .+ scaler[i] * dt * ku[i, :, :]
        vns = vn .+ scaler[i] * dt * kv[i, :, :]
        phins = phin .+ scaler[i] * dt * kphi[i, :, :]
    end

    ku[4, :, :] = u_rhs(parms, uns, phins, vns)
    kv[4, :, :] = v_rhs(parms, vns, phins, uns)
    kphi[4, :, :] = phi_rhs(parms, phins, uns, vns, fphi)

    unp1 = un + (dt / 6.) * (ku[1, :, :] .+ 2 * ku[2, :, :] .+ 2 * ku[3, :, :] .+ ku[4, :, :])
    vnp1 = vn + (dt / 6.) * (kv[1, :, :] .+ 2 * kv[2, :, :] .+ 2 * kv[3, :, :] .+ kv[4, :, :])
    phinp1 = phin + (dt / 6.) * (kphi[1, :, :] .+ 2 * kphi[2, :, :] .+ 2 * kphi[3, :, :] .+ kphi[4, :, :])

    return unp1, phinp1, vnp1
end
