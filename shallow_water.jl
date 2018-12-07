include("./lib/methods.jl")
include("./lib/utils.jl")


function shallow_water(dx, dy, H, dt, nt, output, x0, y0, width)

    # define constant variables
    g = 9.8
    vg = sqrt(H * g)

    # number of timesteps and interval to
    # save new data
    nsteps = trunc(Int, (nt * 3600.) / dt)
    nskip = trunc(Int, (output * 3600.) / dt)

    lons, lats, yy, xx = gridgen(-40, 40, -80, 80, dx=dx, dy=dy)
    ulons, ulats, uyy, uxx = ugridgen(-40, 40, -80, 80, dx=dx, dy=dy)
    vlons, vlats, vyy, vxx = vgridgen(-40, 40, -80, 80, dx=dx, dy=dy)

    ny, nx = size(lons)
    uny, unx = size(ulons)
    vny, vnx = size(vlons)

    # time grid
    tmin = 0
    tmax = (nsteps - 1) * dt

    t = tmin:(dt * nskip):(tmax + dt * nskip)

    re = width * dx
    fphi = 0.01 * gaussian(xx, yy, re, re, nx, ny, x0=x0, y0=y0)

    # initial values
    phin = zeros((ny, nx))
    un = zeros((uny, unx))
    vn = zeros((vny, vnx))

    uall = reshape(un, (1, uny, unx))
    vall = reshape(vn, (1, vny, vnx))
    phiall = reshape(phin, (1, ny, nx))

    parms = (vg, ulats, vlats, uyy, vyy, nx, ny, dx, dy, dt)

    for n in 1:nsteps

        if n % nskip == 0 
            uall = vcat(uall, reshape(un, (1, uny, unx)))
            phiall = vcat(phiall, reshape(phin, (1, ny, nx)))
            vall = vcat(vall, reshape(vn, (1, vny, vnx)))
        end

        unp1, phinp1, vnp1 = rk4(parms, dt, un, phin, vn, fphi)

        un = unp1 * 1.
        phin = phinp1 * 1.
        vn = vnp1 * 1.
    end

    # unstagger
    uall = (uall[:, :, 1:end-1] + uall[:, :, 2:end]) / 2.
    vall = (vall[:, 1:end-1, :] + vall[:, 2:end, :]) / 2.

    _div = divergence(uall, vall, dx, dy)
    _curl = curl(uall, vall, dx, dy)


    return phiall, uall, vall, fphi, _div, _curl, lons, lats, t
end

#shallow_water(100000, 100000, 250, 400, 72, 1, -1.5e6, -1.5e6, 10)