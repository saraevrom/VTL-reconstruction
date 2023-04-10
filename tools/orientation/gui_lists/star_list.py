from .base import GUIList
from orientation.database_reader import get_database, name_a_star
import re
from tkinter.simpledialog import askstring
from vtl_common.localization import get_locale


def match_by_group(database, match, field):
    return database[database[field] == match.group()]

def match_by_first_group(database, match, field, cast_type=int):
    target = cast_type(match.groups()[0])
    #print("MATCH", field, target)
    res = database[database[field] == target]
    #print(res)
    return res

filters = [
    [r"\s+", None],
    [r"Gliese\s*(\d+)", lambda x, y: match_by_first_group(x, y, "gl")],
    [r"HIP\s*(\d+)", lambda x, y: match_by_first_group(x, y, "hip")],
    [r"HR\s*(\d+)", lambda x, y: match_by_first_group(x, y, "hr")],
    [r"Gl\s*(\d+)", lambda x, y: match_by_first_group(x, y, "gl")],
    [r"\d+\w+\s+\w+", lambda x, y: match_by_group(x, y, "bf")],
    [r"\w+", lambda x, y: match_by_group(x, y, "proper")],
]

for i in range(len(filters)):
    filters[i][0] = re.compile(filters[i][0])

def get_ids(entries):
    return [item["id"] for item in entries]


class StarEntry(object):
    def __init__(self, record):
        self.record = record

    def __eq__(self, other):
        return self.record["id"] == other.record["id"]

    def name(self):
        return name_a_star(self.record)

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
                        items.append(StarEntry(filtered.to_dict('records')[0]))
                        print("Found",items[-1].name())
        print("PARSE COMPLETED")
        return items

    @staticmethod
    def represent_item(item):
        return item.name()

    @staticmethod
    def modify_contents(container, new_item):
        if new_item not in container:
            container.append(new_item)