from datetime import datetime
import numpy as np
from .slowmath import Vector3, Quaternion, Vector2



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
    return Vector3(x, y, z)

def eci_to_radec(v:Vector3):
    dec = np.arcsin(v.z)
    ra = np.arctan2(v.y, v.x)

    return ra, dec


def eci_to_ecef(era):
    '''
    Transform direction from ECI to ECEF cartesian coordinates
    :param x_eci: x coordinate in ECI
    :param y_eci: y coordinate in ECI
    :param z_eci: z coordinate in ECI
    :param era: earth rotation angle in radians
    :return: 3D vector coordinates with length 1.0
    '''
    return Quaternion.rotate_xy(-era)

def eci_to_ocef(era, lat, lon, self_rot=0):
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
    q0 = Quaternion.rotate_xy(- era - lon)  # longitude offsets ERA
    # Meridian is aligned
    q1 = Quaternion.rotate_xz(-lat)*q0
    return Quaternion.rotate_yz(self_rot)*q1


def ocef_to_eci(era, lat, lon):
    q0 = Quaternion.rotate_xz(lat)
    return Quaternion.rotate_xy(era + lon)*q0


def detector_to_local(ra, ha, eigen,lat,lon):
    norot = Quaternion.rotate_yz(-eigen)
    eci = ocef_to_eci(0, ra, ha)
    return eci_to_ocef(0,lat,lon)*eci*norot

def local_to_detector(ra, ha, eigen, lat, lon):
    return detector_to_local(ra, ha, eigen, lat, lon).inverse()


def ocef_to_detector_plane(v_local:Vector3, focal_distance):
    '''
    Transform direction from OCEF to detector coordinates as it was pointed in zenith
    :param v_local: coordinates in OCEF
    :param focal_distance: Focal distance of lens in detector
    :return: (x,y,v) x,y - image position on detector. v - if image is visible
    '''
    x_local = v_local.x
    y_local = v_local.y
    z_local = v_local.z
    v = x_local > 0
    x_p = -y_local*focal_distance/x_local
    y_p = z_local*focal_distance/x_local
    return Vector2(x_p, y_p), v


def ocef_to_altaz(v_local:Vector3, backend=None):
    '''
    Transform direction from OCEF to altitude/azimuth
    :param x_local: x coordinate OCEF
    :param y_local: y coordinate OCEF
    :param z_local: z coordinate OCEF
    :return: altitude, azimuth in radians. Azimuth is counted from north towards east
    '''
    x_local = v_local.x
    y_local = v_local.y
    z_local = v_local.z
    if backend is None:
        backend = np
    horizontal = backend.sqrt(y_local**2 + z_local**2)
    alt = backend.arctan2(x_local, horizontal)
    az = (backend.arctan2(y_local, z_local) + 2*np.pi) % (2*np.pi)
    return alt, az

def altaz_to_ocef(alt, az):
    x = np.sin(alt)
    z = np.cos(alt)*np.cos(az)
    y = np.cos(alt)*np.sin(az)
    return Vector3(x, y, z)


def unixtime_to_era(unixtime):
    ut1 = unixtime/86400.0 - 10957.5
    return np.pi*2*(0.7790572732640 + 1.00273781191135448*ut1) % (2*np.pi)


def datetime_to_era(dt:datetime):
    ref_dt = datetime(1970, 1, 1)
    unixtime = (dt - ref_dt).total_seconds()
    return unixtime_to_era(unixtime)


def radec_to_ocef(ra, dec, lat, lon, self_rotation, era) -> Vector3:
    v_eci = radec_to_eci(ra, dec)
    q_ocef = eci_to_ocef(era, lat, lon)
    return Quaternion.rotate_yz(self_rotation)*q_ocef*v_eci

def ocef_to_radec(ocef:Vector3,lat,lon,self_rotation,era):
    q = Quaternion.rotate_yz(self_rotation)*eci_to_ocef(era, lat, lon)
    q_rev = q.inverse()
    v_eci:Vector3 = q_rev*ocef
    return eci_to_radec(v_eci)

def ecef_to_ocef(lat, lon, self_rot=0):
    '''
    Transform direction from ECEF to OCEF cartesian coordinates
    :param lat:
    :param lon:
    :param self_rot:
    :return:
    '''
    return eci_to_ocef(0, lat, lon, self_rot=self_rot)


def detector_plane_to_ocef(v_pdf:Vector2, focal_distance):
    x_ocef = np.ones(v_pdf.x.shape)
    y_ocef = -v_pdf.x/focal_distance
    z_ocef = v_pdf.y/focal_distance
    return Vector3(x_ocef, y_ocef, z_ocef).normalized()


def detector_plane_to_ocef_f(v_pdf:Vector2, focal_distance):
    x_ocef = 1
    y_ocef = -v_pdf.x/focal_distance
    z_ocef = v_pdf.y/focal_distance
    return Vector3(x_ocef, y_ocef, z_ocef).normalized()
