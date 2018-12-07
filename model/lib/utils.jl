function gaussian(xx, yy, rex, rey, nx, ny; x0=0, y0=0)
    
    f0 = zeros((ny, nx))

    for i in 1:ny
        for j in 1:nx
            f0[i, j] = exp(-((xx[i, j] - x0)^2 / rex^2 + (yy[i, j] - y0)^2 / rey^2));
        end
    end

    return f0
end

function gridgen(lati, latf, loni, lonf; dx=100000, dy=100000)
    # creates latitude and longitude grid based on initial
    # anf final lat/lon and grid spacing

    R = 6.371e6
    rlati = lati * pi / 180.
    rlatf = latf * pi / 180.
    rloni = loni * pi / 180.
    rlonf = lonf * pi / 180.

    dlat = dy / R
    dlon = dx / R
    lat = rlati:dlat:(rlatf + dlat)
    lon = rloni:dlon:(rlonf + dlon)

    DY = (rlatf - rlati) * R
    DX = (rlonf - rloni) * R
    y = (-DY / 2.):dy:(DY / 2. + dy)
    x = (-DX / 2.):dx:(DX / 2. + dx)

    nx = size(lon)[1]
    ny = size(lat)[1]

    lats = zeros((ny, nx))
    lons = zeros((ny, nx))
    yy = zeros((ny, nx))
    xx = zeros((ny, nx))

    for i in 1:ny
        for j in 1:nx
            lats[i, j] = lat[i]
            lons[i, j] = lon[j]
            yy[i, j] = y[i]
            xx[i, j] = x[j]
        end
    end

    lats = lats * 180 / pi
    lons = lons * 180 / pi

    return lons, lats, yy, xx
end

function ugridgen(lati, latf, loni, lonf; dx=100000, dy=100000)
    # creates latitude and longitude grid based on initial
    # anf final lat/lon and grid spacing

    R = 6.371e6
    rlati = lati * pi / 180.
    rlatf = latf * pi / 180.
    rloni = loni * pi / 180.
    rlonf = lonf * pi / 180.

    dlat = dy / R
    dlon = dx / R
    lat = rlati:dlat:(rlatf + dlat)
    lon = (rloni - dlon / 2.):dlon:(rlonf + dlon + dlon / 2.)

    DY = (rlatf - rlati) * R
    DX = (rlonf - rloni) * R
    y = (-DY / 2.):dy:(DY / 2. + dy)
    x = (-DX / 2. - dx/2.):dx:(DX / 2. + dx + dx/2.)

    nx = size(lon)[1]
    ny = size(lat)[1]

    lats = zeros((ny, nx))
    lons = zeros((ny, nx))
    yy = zeros((ny, nx))
    xx = zeros((ny, nx))

    for i in 1:ny
        for j in 1:nx
            lats[i, j] = lat[i]
            lons[i, j] = lon[j]
            yy[i, j] = y[i]
            xx[i, j] = x[j]
        end
    end

    lats = lats * 180 / pi
    lons = lons * 180 / pi

    return lons, lats, yy, xx
end

function vgridgen(lati, latf, loni, lonf; dx=100000, dy=100000)
    # creates latitude and longitude grid based on initial
    # anf final lat/lon and grid spacing

    R = 6.371e6
    rlati = lati * pi / 180.
    rlatf = latf * pi / 180.
    rloni = loni * pi / 180.
    rlonf = lonf * pi / 180.

    dlat = dy / R
    dlon = dx / R
    lat = (rlati - dlat / 2.):dlat:(rlatf + dlat + dlat / 2.)
    lon = rloni:dlon:(rlonf + dlon)

    DY = (rlatf - rlati) * R
    DX = (rlonf - rloni) * R
    y = (-DY / 2. - dy/2.):dy:(DY / 2. + dy + dy/2.)
    x = (-DX / 2.):dx:(DX / 2. + dx)

    nx = size(lon)[1]
    ny = size(lat)[1]

    lats = zeros((ny, nx))
    lons = zeros((ny, nx))
    yy = zeros((ny, nx))
    xx = zeros((ny, nx))

    for i in 1:ny
        for j in 1:nx
            lats[i, j] = lat[i]
            lons[i, j] = lon[j]
            yy[i, j] = y[i]
            xx[i, j] = x[j]
        end
    end

    lats = lats * 180 / pi
    lons = lons * 180 / pi

    return lons, lats, yy, xx
end

function divergence(u, v, dx, dy)

    nt, ny, nx = size(u)
    dudx = zeros((nt, ny, nx))
    dvdy = zeros((nt, ny, nx))

    for n in 1:nt
        for i in 1:ny
            for j in 2:(nx-1)
                dudx[n, i, j] = (u[n, i, j+1] - u[n, i, j-1]) / (2 * dx)
            end
            dudx[n, i, nx] = (u[n, i, nx] - u[n, i, nx-1]) / dx
            dudx[n, i, 1] = (u[n, i, 2] - u[n, i, 1]) / dx
        end
    end

    for n in 1:nt
        for j in 1:nx
            for i in 2:(ny-1)
                dvdy[n, i, j] = (v[n, i+1, j] - v[n, i-1, j]) / (2 * dy)
            end
            dvdy[n, ny, j] = (v[n, ny, j] - v[n, ny-1, j]) / dy
            dvdy[n, 1, j] = (v[n, 2, j] - v[n, 1, j]) / dy
        end
    end

    _div = dudx + dvdy

    return _div
end

function curl(u, v, dx, dy)

    nt, ny, nx = size(u)
    dudy = zeros((nt, ny, nx))
    dvdx = zeros((nt, ny, nx))

    for n in 1:nt
        for i in 1:ny
            for j in 2:(nx-1)
                dvdx[n, i, j] = (v[n, i, j+1] - v[n, i, j-1]) / (2 * dx)
            end
            dvdx[n, i, nx] = (v[n, i, nx] - v[n, i, nx-1]) / dx
            dvdx[n, i, 1] = (v[n, i, 2] - v[n, i, 1]) / dx
        end
    end

    for n in 1:nt
        for j in 1:nx
            for i in 2:(ny-1)
                dudy[n, i, j] = (u[n, i+1, j] - u[n, i-1, j]) / (2 * dy)
            end
            dudy[n, ny, j] = (u[n, ny, j] - u[n, ny-1, j]) / dy
            dudy[n, 1, j] = (u[n, 2, j] - u[n, 1, j]) / dy
        end
    end

    _curl = dvdx - dudy

    return _curl
end