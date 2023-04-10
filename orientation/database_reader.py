import pandas as pd
import os.path as ospath
import gc
from vtl_common.parameters import MAX_STAR_MAGNITUDE
from .stellar_math import radec_to_eci
import numpy as np


def name_a_star(row):
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
        if len(data) > 3:
            data = data[:3]
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
