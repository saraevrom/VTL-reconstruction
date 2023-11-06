from .base import GUIList
from ..orientation.database_reader import get_database
import re
from tkinter.simpledialog import askstring
from vtl_common.localization import get_locale
from ..orientation.database_reader import StarEntry


BAYER_LIST = "Alp Bet Gam Del Eps Zet Eta The Iot Kap Lam Mu Nu Xi Omi Pi Rho Sig Tau Ups Phi Chi Psi Ome".split(" ")
BAYER_LETTER_PART = "(?:"+"|".join(BAYER_LIST)+")"
BAYER_REGEX = rf"{BAYER_LETTER_PART}(?:-\d)?"

def match_by_group(database, match, field):
    return database[database[field] == match.group()]

def match_by_first_group(database, match, field, cast_type=int):
    target = cast_type(match.groups()[0])
    #print("MATCH", field, target)
    res = database[database[field] == target]
    #print(res)
    return res

def match_by_dual(database, match, field1, field2, cast1=str, cast2=str, grp1=0, grp2=1):
    #print(f"{field1}-{field2} match")
    grps = match.groups()
    target1 = cast1(grps[grp1])
    target2 = cast2(grps[grp2])
    #print(f"targets {target1}; {target2}")
    m1 = database[field1] == target1
    m2 = database[field2] == target2
    res = database[m1 & m2]
    return res


filters = [
    [r"\s+", None],
    [r"Gliese\s*(\d+)", lambda x, y: match_by_first_group(x, y, "gl")],
    [r"HIP\s*(\d+)", lambda x, y: match_by_first_group(x, y, "hip")],
    [r"HR\s*(\d+)", lambda x, y: match_by_first_group(x, y, "hr")],
    [r"Gl\s*(\d+)", lambda x, y: match_by_first_group(x, y, "gl")],
    [rf'({BAYER_REGEX})\s+(\w+)', lambda x, y: match_by_dual(x,y,"bayer","con")], # Bayer match
    [r'(\d+)\s+(\w+)', lambda x, y: match_by_dual(x,y,"flam","con", cast1=int)], # Flamsteed match
    [r"\w+", lambda x, y: match_by_group(x, y, "proper")],
]

for i in range(len(filters)):
    filters[i][0] = re.compile(filters[i][0])

def get_ids(entries):
    return [item["id"] for item in entries]





class StarList(GUIList):
    def __init__(self, master):
        super().__init__(master, "orientation.selection.stars")

    def obtain_items(self):
        asked_string = askstring(title=get_locale("orientation.selection.stars.dialog.title"),
                                 prompt=get_locale("orientation.selection.stars.dialog.prompt"))
        if not asked_string:
            return None
        ptr = 0
        items = []
        database = get_database()
        while ptr<len(asked_string):
            mat_pair = None
            for filt in filters:
                regex: re.Pattern
                regex, processor = filt
                mat = regex.match(asked_string[ptr:])
                if mat:
                    mat_pair = mat,processor
                    break
            if mat_pair is None:
                print("Unknown pattern on", ptr)
                return None
            else:
                mat, proc = mat_pair
                ptr += mat.end()
                if proc is not None:
                    filtered = proc(database,mat)
                    if filtered.shape[0] == 1:
                        star = StarEntry(filtered.to_dict('records')[0])
                        if star.analyzable():
                            items.append(star)
                            print("Found",items[-1].name())
                        else:
                            print("NOT ANALYZABLE", items[-1].name())
        print("PARSE COMPLETED")
        return items

    @staticmethod
    def represent_item(item):
        return item.name()


