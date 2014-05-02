from math import pi, sqrt, sin, cos, fabs, atan


def cartesian2geodetic(x, y, z):
    ''' ECEF2WGS84 converts earth-centered, earth-fixed coordinates to geodetic
    latitude, longitude, and HAE (WGS84) using the improved Bowring method.

    [LAT, LON, ELE] = ECEF2WGS84(X, Y, Z) converts ECEF cartesian coordinates
    X, Y, Z in meters into an array of coordinates. Both LAT and LON are in
    degrees. ELE is in meters.

    References:
    [1] Toms, R. M., "An Efficient Algorithm for Geocentric to Geodetic
    Coordinate Conversion", Lawrence Livermore National Laboratory, 1995,
    http://www.osti.gov/energycitations/servlets/purl/110235-MLuPMu/
    webviewable/110235.pdf
    '''
    a = 6378137.0       # WGS84 equatorial radius (m)
    c = 6356752.3142    # WGS84 polar radius (m)
    aD_c = 1.0026        # only valid below 2000 km HAE

    f = (a - c) / a
    e2 = (a * a - c * c) / (a * a)
    E2 = (a * a - c * c) / (c * c)

    w = sqrt(x * x + y * y)
    t_0 = z * aD_c
    s_0 = sqrt((z * aD_c) ** 2 + w ** 2)

    sin_b_0 = t_0 / s_0
    cos_b_0 = w / s_0

    t_1 = z + c * E2 * sin_b_0 ** 3
    s_1 = sqrt(t_1 ** 2 + (w - a * e2 * cos_b_0 ** 3) ** 2)

    sin_p_1 = t_1 / s_1
    cos_p_1 = (w - a * e2 * cos_b_0 ** 3) / s_1

    r_n = a / sqrt(1 - e2 * sin_p_1 ** 2)

    if fabs(cos_p_1) >= cos(67.5 * (pi / 180)):
        ele = w / fabs(cos_p_1) - r_n
    else:
        ele = z / sin_b_0 + r_n * (e2 - 1)

    lat = atan(sin_p_1 / cos_p_1) * (180 / pi)
    lon = atan(y / x) * (180 / pi)

    return lat, lon, ele


def geodetic2cartesian(lat, lon, ele):
    '''
    References:
    [1] https://geodesy.curtin.edu.au/local/docs/1Featherstone295.pdf
    '''
    a = 6378137.0       # WGS84 equatorial radius in meters (semi-major axis)
    b = 6356752.3142    # WGS84 polar radius in meters (semi-minor axis)

    f = (a - b) / a     # Ellipsoid flattening factor
    e = sqrt(2 * f - f * f) # First numerical eccentricity of ellipsoid

    lat = lat * (pi / 180)
    lon = lon * (pi / 180)

    v = a / sqrt(1 - e * e * sin(lat) * sin(lat)) # Radius of curvature in the prime vertical of the surface of the geodetic ellipsoid

    x = (v + ele) * cos(lat) * cos(lon)
    y = (v + ele) * cos(lat) * sin(lon)
    z = (v * (1 - e * e) + ele) * sin(lat)

    return x, y, z

