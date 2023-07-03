from copy import deepcopy, copy

class Unchanged(object):

    def __init__(self,v):
        self.value = v

    def __repr__(self):
        return f"<UNCHANGED {self.value}>"

    def __str__(self):
        return "<UNCHANGED {self.value}>"


class Deleted(object):
    def __repr__(self):
        return "<DELETED>"

    def __str__(self):
        return "<DELETED>"


DELETED = Deleted()

def diff(struct_a, struct_b):
    '''
    Takes two JSON serializable objects and makes a new one that represents difference
    :param struct_a: start object
    :param struct_b: end object
    :return: diff object
    '''
    if isinstance(struct_a, dict) and isinstance(struct_b, dict):
        set_a = set(struct_a.keys())
        set_b = set(struct_b.keys())

        deleted = set_a - set_b
        new = set_b - set_a
        inner = set_a.intersection(set_b)
        res_dict = dict()
        for i in deleted:
            res_dict[i] = DELETED

        for i in inner:
            res_dict[i] = diff(struct_a[i], struct_b[i])

        for i in new:
            res_dict[i] = struct_b[i]

        return res_dict

    elif isinstance(struct_a, list) and isinstance(struct_b, list):
        # List-list delta case
        res_list = []
        for i, item_b in enumerate(struct_b):
            if i>=len(struct_a):
                res_list.append(item_b)
            else:
                res_list.append(diff(struct_a[i], item_b))
        return res_list
    else:
        if type(struct_a) == type(struct_b) and struct_a == struct_b:
            return Unchanged(struct_b)
        else:
            return struct_b

def patch(struct_a, diff_obj):
    '''
    Applies changes obtained by diff() function to struct
    :param struct_a: start object
    :param diff_obj: diff object
    :return: object with applied patch
    '''

    if isinstance(struct_a, dict) and isinstance(diff_obj, dict):
        res = deepcopy(struct_a)
        for k in diff_obj.keys():
            if diff_obj[k] == DELETED and k in res:
                del res[k]
            else:
                if k in res.keys():
                    res[k] = patch(res[k], diff_obj[k])
                elif isinstance(diff_obj[k], Unchanged):
                    res[k] = diff_obj[k].value
                else:
                    res[k] = diff_obj[k]
        return res
    elif isinstance(struct_a, list) and isinstance(diff_obj, list):
        res = []
        for i, diff_i in enumerate(diff_obj):
            if i >= len(struct_a):
                if isinstance(diff_i, Unchanged):
                    res.append(diff_i.value)
                else:
                    res.append(diff_i)
            else:
                res.append(patch(struct_a[i], diff_obj))
        return res
    elif isinstance(diff_obj, Unchanged):
        return struct_a
    else:
        return copy(diff_obj)
