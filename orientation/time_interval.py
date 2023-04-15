import datetime
from vtl_common.datetime_parser import datetime_to_unixtime
from vtl_common.parameters import DATETIME_FORMAT

class TimeRange(object):
    def __init__(self, start:datetime.datetime, end:datetime.datetime, stride=1):
        self.start = start
        self.end = end
        self.stride = stride

    @staticmethod
    def from_unixtime(ut_start, ut_end,stride=1):
        start = datetime.datetime.utcfromtimestamp(ut_start)
        end = datetime.datetime.utcfromtimestamp(ut_end)
        return TimeRange(start,end,stride)

    def unixtime_intervals(self):
        return datetime_to_unixtime(self.start), datetime_to_unixtime(self.end)


    def intersects(self, other):
        '''
        Checks if given interval intersects with another one
        :param other:
        :return:
        '''
        if other.start<=self.start<=other.end:
            return True
        if other.start <= self.end <= other.end:
            return True
        if self.start<=other.start<=self.end:
            return True
        if self.start <= other.end <= self.end:
            return True
        return False

    def name(self):
        if self.start == self.end:
            s = self.start.strftime(DATETIME_FORMAT)
        else:
            s = self.start.strftime(DATETIME_FORMAT)+" : "+self.end.strftime(DATETIME_FORMAT)

        if self.stride!=1:
            s += " : "+str(self.stride)
        return s

    def unify(self, other):
        return TimeRange(min(self.start, other.start), max(self.end, other.end), max(self.stride, other.stride))

    def intersection(self, other):
        '''
        Returns intersection interval in case of intersection with another one
        :param other:
        :return:
        '''
        if self.intersects(other):
            return TimeRange(max(self.start, other.start), min(self.end, other.end), max(self.stride, other.stride))

        return None

    def __lt__(self, other):
        return self.start < other.start

    def __le__(self, other):
        return self.start <= other.start

    def __gt__(self, other):
        return self.start > other.start

    def __ge__(self, other):
        return self.start >= other.start

    def __eq__(self, other):
        return self.start == other.start
