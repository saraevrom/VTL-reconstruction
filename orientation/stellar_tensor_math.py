import numpy as np
import pymc as pm
import pytensor.tensor as pt

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

def rotate_yz_pt(x, y, z, ang):
    '''
    Rotates given vector in plane YZ (around X)
    :param x:
    :param y:
    :param z:
    :param ang: angle in radians
    :return:
    '''
    cos_era = pt.cos(ang)
    sin_ang = pt.sin(ang)
    x1 = x
    y1 = y * cos_era - z * sin_ang
    z1 = y * sin_ang + z * cos_era
    return x1, y1, z1


def rotate_zy_pt(x, y, z, ang):
    return rotate_yz_pt(x, y, z, -ang)


def rotate_xz_pt(x, y, z, ang):
    '''
    Rotates given vector in plane XZ (around Y)
    :param x:
    :param y:
    :param z:
    :param ang:
    :return:
    '''
    cos_era = pt.cos(ang)
    sin_ang = pt.sin(ang)
    y1 = y
    x1 = x * cos_era - z * sin_ang
    z1 = x * sin_ang + z * cos_era
    return x1, y1, z1


def rotate_zx_pt(x, y, z, ang):
    return rotate_xz_pt(x, y, z, -ang)


def rotate_xy_pt(x,y,z, ang):
    '''

    Rotates given vector in plane XY (around Z)
    :param x:
    :param y:
    :param z:
    :param ang:
    :return:
    '''
    cos_era = pt.cos(ang)
    sin_ang = pt.sin(ang)
    z1 = z
    x1 = x * cos_era - y * sin_ang
    y1 = x * sin_ang + y * cos_era
    return x1, y1, z1


def rotate_yx_pt(x, y, z, ang):
    return rotate_xy_pt(x, y, z, -ang)


def eci_to_ecef_pt(x_eci, y_eci, z_eci, era):
    '''
    Transform direction from ECI to ECEF cartesian coordinates
    :param x_eci: x coordinate in ECI
    :param y_eci: y coordinate in ECI
    :param z_eci: z coordinate in ECI
    :param era: earth rotation angle in radians
    :return: 3D vector coordinates with length 1.0
    '''
    return rotate_xy_pt(x_eci, y_eci, z_eci, -era)



def eci_to_ocef_pt(x_eci, y_eci, z_eci, era, lat, lon):
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
    x_0, y_0, z_0 = rotate_xy_pt(x_eci, y_eci, z_eci, - era - lon)  # longitude offsets ERA
    # Meridian is aligned
    return rotate_xz_pt(x_0, y_0, z_0, -lat)


def ecef_to_ocef_pt(x_ecef, y_ecef, z_ecef, lat, lon):
    '''
    Transform direction from ECEF to OCEF cartesian coordinates
    :param x_ecef:
    :param y_ecef:
    :param z_ecef:
    :param lat:
    :param lon:
    :return:
    '''
    return eci_to_ocef_pt(x_ecef, y_ecef, z_ecef, 0, lat, lon) # When ERA=0 ECI and ECEF are same

def ocef_to_detector_plane_pt(x_local, y_local, z_local, focal_distance):
    '''
    Transform direction from OCEF to detector coordinates
    :param x_local: x coordinate OCEF
    :param y_local: y coordinate OCEF
    :param z_local: z coordinate OCEF
    :param focal_distance: Focal distance of lens in detector
    :return: (x,y,v) x,y - image position on detector. v - if image is visible
    '''
    v = pt.switch(x_local > 0, 1.0, 0.0)
    x_p = -y_local*focal_distance/x_local
    y_p = z_local*focal_distance/x_local
    return x_p, y_p, v


def ocef_to_altaz_pt(x_local, y_local, z_local):
    '''
    Transform direction from OCEF to altitude/azimuth
    :param x_local: x coordinate OCEF
    :param y_local: y coordinate OCEF
    :param z_local: z coordinate OCEF
    :return: altitude, azimuth in radians. Azimuth is counted from north towards east
    '''
    horizontal = pt.sqrt(y_local**2 + z_local**2)
    alt = pt.arctan(x_local/horizontal)
    az = (pt.arctan2(y_local, z_local) + 2*np.pi) % (2*np.pi)
    return alt, az
