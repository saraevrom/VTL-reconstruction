from datetime import datetime

import numpy
import numpy as np
from .slowmath import Vector3, Quaternion, Vector2



'''
Alternative astronomy math module
Glossary:
ECI -- Earth-Centered Inertial coordinate system. 
    Origin -- center of earth
    Z - points to vernal equinox. 
    Y - aligned with Earth rotation axis
    X - chosen to make XYZ left-handed coordinate system
ECEF --  Earth-Centered Earth-Fixed coordinate system
    Origin -- center of earth
    Z - lies in a plane of equator, extending to the prime meridian.
    Y - aligned with Earth rotation axis
    X - chosen to make XYZ left-handed coordinate system
OCEF -- Observer-Centered Earth-Fixed coordinate system
    Origin -- observer
    Z - points to zenith
    X - points to east
    Y - points to north
'''


def radec_to_eci(ra, dec,backend=numpy):
    '''
    Transform direction from stellar coordinates (RA, dec) to ECI cartesian coordinates
    :param ra: Right ascension in radians
    :param dec: Declination in radians
    :return: 3D vector coordinates with length 1.0
    '''
    cos_dec = backend.cos(dec)
    y = backend.sin(dec)
    z = cos_dec * backend.cos(ra)
    x = cos_dec * backend.sin(ra)
    return Vector3(x, y, z)


def eci_to_radec(v:Vector3,backend=numpy):
    dec = backend.arcsin(v.y)
    ra = backend.arctan2(v.x, v.z)
    return ra, dec


def latlon_quaternion(lat,lon, backend=numpy):
    return Quaternion.rotate_zx(lon, backend)*Quaternion.rotate_zy(lat, backend)


def eci_to_ocef(era, dec, gha, self_rot=0, backend=numpy):
    '''
    Transform direction from ECI to OCEF cartesian coordinates
    :param backend: math backend
    :param era: Earth rotation angle
    :param dec: Declination, radians
    :param gha: Greenwich hour angle, radians
    :param self_rot: Own eastward rotation
    :return:
    '''
    ra = gha+era
    return Quaternion.rotate_xy(self_rot, backend=backend)*latlon_quaternion(dec, ra, backend=backend).conj()


def ocef_to_eci(era, dec, gha, self_rot=0, backend=numpy):
    return eci_to_ocef(era, dec, gha, self_rot, backend=backend).conj()


def detector_to_local(dec, gha, self_rot, lat,lon, backend=numpy):
    eci = ocef_to_eci(0, dec, gha, self_rot, backend=backend)
    return eci_to_ocef(0,lat,lon, backend=backend)*eci


def local_to_detector(dec, gha, eigen, lat, lon, backend=numpy):
    return detector_to_local(dec, gha, eigen, lat, lon, backend=backend).conj()


def ocef_to_detector_plane(v_local:Vector3, focal_distance):
    '''
    Transform direction from OCEF to detector coordinates as it was pointed in zenith
    :param v_local: coordinates in OCEF
    :param focal_distance: Focal distance of lens in detector
    :return: (pos,v) pos - image position on detector. v - if image is visible
    '''
    x_local = v_local.x
    y_local = v_local.y
    z_local = v_local.z
    v = z_local > 0
    x_p = -x_local*focal_distance/z_local
    y_p = y_local*focal_distance/z_local
    return Vector2(x_p, y_p), v


def ocef_to_altaz(v_local:Vector3, backend=None, allow_neg=False):
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
    horizontal = backend.sqrt(x_local**2 + y_local**2)
    alt = backend.arctan2(z_local, horizontal)
    az = backend.arctan2(x_local, y_local)
    if not allow_neg:
        az = (az+2*np.pi) % (2*np.pi)
    return alt, az


def altaz_to_ocef(alt, az, backend=numpy):
    z = backend.sin(alt)
    y = backend.cos(alt)*backend.cos(az)
    x = backend.cos(alt)*backend.sin(az)
    return Vector3(x, y, z)


def unixtime_to_era(unixtime):
    ut1 = unixtime/86400.0 - 10957.5
    return np.pi*2*(0.7790572732640 + 1.00273781191135448*ut1) % (2*np.pi)


def datetime_to_era(dt:datetime):
    ref_dt = datetime(1970, 1, 1)
    unixtime = (dt - ref_dt).total_seconds()
    return unixtime_to_era(unixtime)


def radec_to_ocef(ra, dec, lat, lon, self_rotation, era, backend=numpy) -> Vector3:
    v_eci = radec_to_eci(ra, dec)
    q_ocef = eci_to_ocef(era, lat, lon, self_rotation, backend=backend)
    return q_ocef*v_eci


def ocef_to_radec(ocef:Vector3,lat,lon,self_rotation,era, backend=numpy):
    q = eci_to_ocef(era, lat, lon, self_rotation, backend=backend)
    q_rev = q.conj()
    v_eci:Vector3 = q_rev*ocef
    return eci_to_radec(v_eci)

def ecef_to_ocef(lat, lon, self_rot=0, backend=numpy):
    '''
    Transform direction from ECEF to OCEF cartesian coordinates
    :param lat:
    :param lon:
    :param self_rot:
    :return:
    '''
    return eci_to_ocef(0, lat, lon, self_rot=self_rot, backend=backend)


def detector_plane_to_ocef(v_pdf:Vector2, focal_distance):
    z_ocef = np.ones(v_pdf.x.shape)
    x_ocef = -v_pdf.x/focal_distance
    y_ocef = v_pdf.y/focal_distance
    return Vector3(x_ocef, y_ocef, z_ocef).normalized()


def detector_plane_to_ocef_f(v_pdf:Vector2, focal_distance):
    z_ocef = 1
    x_ocef = -v_pdf.x/focal_distance
    y_ocef = v_pdf.y/focal_distance
    return Vector3(x_ocef, y_ocef, z_ocef).normalized()
