from datetime import datetime


def datetime_to_copystr(dt:datetime):
    s = dt.strftime("%Y-%m-%d %H:%M:%S:%f")
    s = s[:-3]
    return s

def ut0_to_copystr(ut0):
    dt = datetime.utcfromtimestamp(ut0)
    return datetime_to_copystr(dt)