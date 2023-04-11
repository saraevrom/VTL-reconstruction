import pandas as pd
import os.path as ospath
import gc
from vtl_common.parameters import MAX_STAR_MAGNITUDE
from .stellar_math import radec_to_eci
import numpy as np
from .stellar_math import eci_to_ocef, rotate_yz, ocef_to_detector_plane


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
        DATABASE = pd.read_csv(CWD+"/hygdata_v3.csv", sep=",")
        DATABASE = DATABASE[(DATABASE["mag"] < MAX_STAR_MAGNITUDE) & (DATABASE["id"] != 0)]  # Remove Sun
        x,y,z = radec_to_eci(DATABASE["ra"].to_numpy()*np.pi/12, DATABASE["dec"].to_numpy()*np.pi/180)

        DATABASE["eci_x"] = x
        DATABASE["eci_y"] = y
        DATABASE["eci_z"] = z
        DATABASE["star_name"] = DATABASE.apply(name_a_star, axis=1)
        gc.collect()
    return DATABASE

VEGA_MAGNITUDE = 0.03
VEGA_LUM = 2.54e-6

class StarEntry(object):
    def __init__(self, record):
        self.record = record

    def magnitude(self):
        return self.record["mag"]

    def energy(self):
        return VEGA_LUM * 10**((VEGA_MAGNITUDE - self.magnitude())/2.5)

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
        x_ocef, y_ocef, z_ocef = rotate_yz(x_ocef, y_ocef, z_ocef, params["SELF_ROTATION"])
        x_pdf, y_pdf, visible = ocef_to_detector_plane(x_ocef, y_ocef, z_ocef, params["FOCAL_DISTANCE"])
        return x_pdf, y_pdf, visible