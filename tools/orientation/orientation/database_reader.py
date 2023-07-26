import pandas as pd
import os.path as ospath
import gc
from vtl_common.parameters import MAX_STAR_MAGNITUDE
from .stellar_math import radec_to_eci
import numpy as np
from .stellar_math import eci_to_ocef, rotate_yz, ocef_to_detector_plane
from vtl_common.parameters import APERTURE # in mm^2


LIGHTSPEED = 299792458. # m/s
NU1 = LIGHTSPEED/300e-9 # Hz
NU2 = LIGHTSPEED/400e-9 # Hz
DELTA_NU = NU1-NU2
EFFICIENCY = 0.15

#PIXEL_AREA = (PIXEL_SIZE*1e-3)**2 # in m^2
JY_TO_W = 1e-23 * 1e-7 * (APERTURE/100) * DELTA_NU * EFFICIENCY
JY_TO_PICO_W = JY_TO_W*1e12


class DatabaseError(Exception):
    def __init__(self, key):
        super().__init__(f"Could not find entry for {key}")


def name_a_star(row, name_limit=3):
    data = []
    if not pd.isnull(row["proper"]):
        data.append(row["proper"])
    if not pd.isnull(row["bf"]):
        data.append(row["bf"])
    if not pd.isnull(row["hip"]) and row["hip"]!=0:
        data.append(f"HIP {int(row['hip'])}")
    if not pd.isnull(row["hr"]):
        data.append(f"HR {int(row['hr'])}")
    if not pd.isnull(row["gl"]):
        data.append(f"Gliese {row['gl']}")
    if data:
        if len(data) > name_limit:
            data = data[:name_limit]
        if len(data) > 1:
            return f"{data[0]} ({', '.join(data[1:])})"
        else:
            return data[0]
    else:
        return f"UNKNOWN (StarID {row['id']})"

DATABASE = None

def get_database():
    global DATABASE
    if DATABASE is None:
        CWD = ospath.dirname(__file__)
        DATABASE = pd.read_csv(CWD+"/hygdata_v3_photometric_mod.csv", sep=",")
        DATABASE = DATABASE[(DATABASE["mag"] < MAX_STAR_MAGNITUDE) & (DATABASE["id"] != 0)]  # Remove Sun
        x,y,z = radec_to_eci(DATABASE["ra"].to_numpy()*np.pi/12, DATABASE["dec"].to_numpy()*np.pi/180)

        DATABASE["eci_x"] = x
        DATABASE["eci_y"] = y
        DATABASE["eci_z"] = z
        DATABASE["star_name"] = DATABASE.apply(name_a_star, axis=1)
        gc.collect()
    return DATABASE

class StarEntry(object):
    def __init__(self, record):
        self.record = record

    def magnitude(self):
        return self.record["mag"]

    def magnitude_b(self):
        '''
        Stellar magnitude in U band
        :return:
        '''
        ci = None
        if "BV" in self.record.keys():
            ci = self.record["BV"]
        if "ci" in self.record.keys():
            ci = self.record["ci"]
        else:
            raise DatabaseError("BV/ci")

        return self.magnitude()+ci   # B = V + BV

    def magnitude_u(self):
        b = self.magnitude_b()
        ub = None
        if "UB" in self.record.keys():
            ub = self.record["UB"]
        else:
            raise DatabaseError("UB")

        return b + ub  # U = B + UB = V + BV + UB

    def analyzable(self):
        return ("BV" in self.record.keys() or "ci" in self.record.keys()) and "UB" in self.record.keys()

    def energy(self):
        return 3640 * JY_TO_PICO_W * 10**(-self.magnitude()/2.5)

    def energy_u(self):
        return 1810 * JY_TO_PICO_W * 10**(-self.magnitude_u()/2.5)


    def __eq__(self, other):
        return self.record["id"] == other.record["id"]

    def name(self):
        return name_a_star(self.record)

    def primary_name(self):
        return name_a_star(self.record, 1)

    def get_eci(self):
        return [self.record[i] for i in ["eci_x", "eci_y", "eci_z"]]

    def position_on_plane(self, params, era):
        x_eci, y_eci, z_eci = self.get_eci()
        x_ocef, y_ocef, z_ocef = eci_to_ocef(x_eci, y_eci, z_eci, era,
                                             lat=params["VIEW_LATITUDE"] * np.pi / 180,
                                             lon=params["VIEW_LONGITUDE"] * np.pi / 180)
        x_ocef, y_ocef, z_ocef = rotate_yz(x_ocef, y_ocef, z_ocef, params["SELF_ROTATION"] * np.pi / 180)
        x_pdf, y_pdf, visible = ocef_to_detector_plane(x_ocef, y_ocef, z_ocef, params["FOCAL_DISTANCE"])
        return x_pdf, y_pdf, visible

