import pandas as pd
import os.path as ospath
import gc
from vtl_common.parameters import MAX_STAR_MAGNITUDE

CWD = ospath.dirname(__file__)
DATABASE = pd.read_csv(CWD+"/hygdata_v3.csv", sep=",")
DATABASE = DATABASE[DATABASE["mag"]<MAX_STAR_MAGNITUDE | DATABASE["id"] != 0]  # Remove Sun
gc.collect()
