import numpy as np
import numba as nb

class Cutter(object):
    def cut(self, plot_data):
        raise NotImplementedError()


class WholeCutter(Cutter):
    def cut(self, plot_data):
        return None

class RangeCutter(Cutter):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def cut(self, plot_data):
        return self.start, self.end


@nb.njit()
def find_first_trigger(curve, thresh):
    n = len(curve)
    for i in range(n):
        if curve[i] >= thresh:
            return i
    return n-1


@nb.njit()
def find_last_trigger(curve, thresh):
    n = len(curve)
    for i in range(n-1, -1, -1):
        if curve[i] >= thresh:
            return i
    return 0


class ThresholdCutter(Cutter):
    def __init__(self, threshold):
        self.threshold = threshold

    def cut(self, plot_data):
        lightcurve = np.max(plot_data, axis=(1,2))
        start = find_first_trigger(lightcurve, self.threshold)
        end = find_last_trigger(lightcurve, self.threshold)
        if end>start:
            print(f"CUTTING FROM {start} to {end}")
            return start, end