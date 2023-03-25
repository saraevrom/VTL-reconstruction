

class Cutter(object):
    def cut(self, plot_data):
        raise NotImplementedError()


class WholeCutter(Cutter):
    def cut(self, plot_data):
        return plot_data

class RangeCutter(Cutter):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def cut(self, plot_data):
        return plot_data[self.start: self.end]
