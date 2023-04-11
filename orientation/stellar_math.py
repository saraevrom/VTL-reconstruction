import numpy as np
import numba as nb

'''
Astronomy math module
Glossary:
ECI -- Earth-Centered Inertial coordinate system. 
    Origin -- center of earth
    X - points to vernal equinox. 
    Z - aligned with Earth rotation axis
    Y - chosen to make XYZ left-handed coordinate system
ECEF --  Earth-Centered Earth-Fixed coordinate system
    Origin -- center of earth
    X - lies in a plane of equator, extending to the prime meridian.
    Z - aligned with Earth rotation axis
    Y - chosen to make XYZ left-handed coordinate system
OCEF -- Observer-Centered Earth-Fixed coordinate system
    Origin -- observer
    X - points to zenith
    Y - points to east
    Z - points to north
'''

@nb.njit()
def rotate_yz(x, y, z, ang):
    '''
    Rotates given vector in plane YZ (around X)
    :param x:
    :param y:
    :param z:
    :param ang: angle in radians
    :return:
    '''
    cos_era = np.cos(ang)
    sin_ang = np.sin(ang)
    x1 = x
    y1 = y * cos_era - z * sin_ang
    z1 = y * sin_ang + z * cos_era
    return x1, y1, z1


@nb.njit()
def rotate_zy(x, y, z, ang):
    return rotate_yz(x, y, z, -ang)


@nb.njit()
def rotate_xz(x, y, z, ang):
    '''
    Rotates given vector in plane XZ (around Y)
    :param x:
    :param y:
    :param z:
    :param ang:
    :return:
    '''
    cos_era = np.cos(ang)
    sin_ang = np.sin(ang)
    y1 = y
    x1 = x * cos_era - z * sin_ang
    z1 = x * sin_ang + z * cos_era
    return x1, y1, z1


@nb.njit()
def rotate_zx(x, y, z, ang):
    return rotate_xz(x, y, z, -ang)

@nb.njit()
def rotate_xy(x,y,z, ang):
    '''

    Rotates given vector in plane XY (around Z)
    :param x:
    :param y:
    :param z:
    :param ang:
    :return:
    '''
    cos_era = np.cos(ang)
    sin_ang = np.sin(ang)
    z1 = z
    x1 = x * cos_era - y * sin_ang
    y1 = x * sin_ang + y * cos_era
    return x1, y1, z1


@nb.njit()
def rotate_yx(x, y, z, ang):
    return rotate_xy(x, y, z, -ang)


@nb.njit()
def radec_to_eci(ra, dec):
    '''
    Transform direction from stellar coordinates (RA, dec) to ECI cartesian coordinates
    :param ra: Right ascension in radians
    :param dec: Declination in radians
    :return: 3D vector coordinates with length 1.0
    '''
    cos_dec = np.cos(dec)
    z = np.sin(dec)
    x = cos_dec * np.cos(ra)
    y = cos_dec * np.sin(ra)
    return x, y, z


@nb.njit()
def eci_to_ecef(x_eci, y_eci, z_eci, era):
    '''
    Transform direction from ECI to ECEF cartesian coordinates
    :param x_eci: x coordinate in ECI
    :param y_eci: y coordinate in ECI
    :param z_eci: z coordinate in ECI
    :param era: earth rotation angle in radians
    :return: 3D vector coordinates with length 1.0
    '''
    return rotate_xy(x_eci, y_eci, z_eci, -era)


@nb.njit()
def eci_to_ocef(x_eci, y_eci, z_eci, era, lat, lon):
    '''
    Transform direction from ECI to OCEF cartesian coordinates
    :param x_eci: x coordinate in ECI
    :param y_eci: y coordinate in ECI
    :param z_eci: z coordinate in ECI
    :param era: Earth rotation angle
    :param lat: Latitude, radians
    :param lon: Longitude, radians
    :return:
    '''
    x_0, y_0, z_0 = rotate_xy(x_eci, y_eci, z_eci, - era - lon)  # longitude offsets ERA
    # Meridian is aligned
    return rotate_xz(x_0, y_0, z_0, -lat)


@nb.njit()
def ocef_to_detector_plane(x_local, y_local, z_local, focal_distance):
    '''
    Transform direction from OCEF to detector coordinates
    :param x_local: x coordinate OCEF
    :param y_local: y coordinate OCEF
    :param z_local: z coordinate OCEF
    :param focal_distance: Focal distance of lens in detector
    :return: (x,y,v) x,y - image position on detector. v - if image is visible
    '''
    v = x_local > 0
    x_p = -y_local*focal_distance/x_local
    y_p = z_local*focal_distance/x_local
    return x_p, y_p, v


def detector_plane_to_ocef(x_pdf, y_pdf, focal_distance):
    x_ocef = np.ones(x_pdf.shape)
    y_ocef = -x_pdf/focal_distance
    z_ocef = y_pdf/focal_distance
    norm = np.sqrt(x_ocef**2+y_ocef**2+z_ocef**2)
    return x_ocef/norm, y_ocef/norm, z_ocef/norm

@nb.njit()
def ocef_to_altaz(x_local, y_local, z_local):
    '''
    Transform direction from OCEF to altitude/azimuth
    :param x_local: x coordinate OCEF
    :param y_local: y coordinate OCEF
    :param z_local: z coordinate OCEF
    :return: altitude, azimuth in radians. Azimuth is counted from north towards east
    '''
    horizontal = np.sqrt(y_local**2 + z_local**2)
    alt = np.arctan(x_local/horizontal)
    az = (np.arctan2(y_local, z_local) + 2*np.pi) % (2*np.pi)
    return alt, az

@nb.njit()
def unixtime_to_era(unixtime):
    ut1 = unixtime/86400.0 - 10957.5
    return np.pi*2*(0.7790572732640 + 1.00273781191135448*ut1) % (2*np.pi)