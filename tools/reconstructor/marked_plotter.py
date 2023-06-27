import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Rectangle
from vtl_common.localized_GUI import GridPlotter
from vtl_common.parameters import HALF_PIXELS, PIXEL_SIZE, HALF_GAP_SIZE

from orientation.stellar_math import radec_to_eci, eci_to_ocef, ocef_to_detector_plane, unixtime_to_era, rotate_yz

SPAN = HALF_PIXELS*PIXEL_SIZE+HALF_GAP_SIZE


def _set_mask(patch, value):
    if value:
        patch.set_alpha(0.0)
    else:
        patch.set_alpha(0.5)


class HighlightingPlotter(GridPlotter):
    def __init__(self, master):
        super().__init__(master)
        self.bottom_left = Rectangle((0, 0), -SPAN, -SPAN, color="gray", alpha=0.0)
        self.bottom_right = Rectangle((0, 0), SPAN, -SPAN, color="gray", alpha=0.0)
        self.top_left = Rectangle((0, 0), -SPAN, SPAN, color="gray", alpha=0.0)
        self.top_right = Rectangle((0, 0), SPAN, SPAN, color="gray", alpha=0.0)

        self.axes.add_patch(self.bottom_left)
        self.axes.add_patch(self.bottom_right)
        self.axes.add_patch(self.top_left)
        self.axes.add_patch(self.top_right)
        self.added_patches = []
        self._origin = (0, 0)  # Where to start draw arrow for Direction pick
        self._point_target=(1,1)
        self._pointer_time_index = 0  # Frame from start for ut0
        self._direction_arrow = self.axes.arrow(0, 0, 1, 1, width=0.2,length_includes_head=True,
                                                edgecolor="green", facecolor="green")
        self._direction_arrow.set_visible(False)
        self._not_visible = self.axes.text(0,0, "N/V", color="red",ha='center', va='center')
        self._not_visible.set_visible(False)
        self._pointer_data = None
        self._times = None


    def set_point_direction(self, formdata, ut0_arr):
        self._pointer_data = formdata
        self._times = ut0_arr
        self.update_point()


    def update_point(self):
        if self._pointer_data is not None:
            formdata = self._pointer_data
            if formdata["show"]:
                ra = formdata["ra"]
                dec = formdata["dec"]
                print("RA=", ra)
                orientation = formdata["orientation"]
                ut0 = self._times[self._pointer_time_index]
                era = unixtime_to_era(ut0)
                lat = orientation["VIEW_LATITUDE"] * np.pi / 180
                lon = orientation["VIEW_LONGITUDE"] * np.pi / 180
                self_rotation = orientation["SELF_ROTATION"] * np.pi / 180
                f = orientation["FOCAL_DISTANCE"]

                x_eci, y_eci, z_eci = radec_to_eci(ra, dec)
                x_ocef, y_ocef, z_ocef = eci_to_ocef(x_eci, y_eci, z_eci, era, lat, lon)
                x_ocef, y_ocef, z_ocef = rotate_yz(x_ocef, y_ocef, z_ocef, self_rotation)
                x, y, v = ocef_to_detector_plane(x_ocef, y_ocef, z_ocef, f)
                if v:
                    self._not_visible.set_visible(False)
                    self.point_to(x, y)
                else:
                    self.set_nv()
            else:
                self.hide_pointer()

    def set_nv(self):
        self._direction_arrow.set_visible(False)
        self._not_visible.set_visible(True)
        self.draw()

    def hide_pointer(self):
        self._direction_arrow.set_visible(False)
        self._not_visible.set_visible(False)
        self.draw()

    def set_origin(self, x, y, k0=None, view_enable=False):
        self._origin = (x, y)
        if k0 is not None:
            self._pointer_time_index = k0
        self.update_point()

    def point_to(self,x,y, view_enable=True):
        self._point_target = (x,y)
        self.update_point_arrow(view_enable)

    def update_point_arrow(self, view_enable=False):
        x0 = self._origin[0]
        y0 = self._origin[1]
        x = self._point_target[0]
        y = self._point_target[1]
        self._direction_arrow.set_data(x=x0, y=y0, dx=x-x0, dy=y-y0)
        if view_enable:
            self._direction_arrow.set_visible(True)
        self.draw()

    def plot_lines(self, *args, **kwargs):
        lines = self.axes.plot(*args,**kwargs)
        for l in lines:
            self.added_patches.append(l)

    def plot_arrow(self, *args, **kwargs):
        arrow = self.axes.arrow(*args, **kwargs)
        self.added_patches.append(arrow)

    def plot_circle(self, *args, **kwargs):
        patch = plt.Circle(*args, **kwargs)
        self.axes.add_patch(patch)
        self.added_patches.append(patch)

    def clear_added_patches(self):
        for x in self.added_patches:
            x.remove()
        self.added_patches.clear()

    def set_mask(self, source):
        _set_mask(self.bottom_left, source["bottom_left"])
        _set_mask(self.bottom_right, source["bottom_right"])
        _set_mask(self.top_left, source["top_left"])
        _set_mask(self.top_right, source["top_right"])

    def set_mask_vars(self, bl, br, tl, tr):
        _set_mask(self.bottom_left, bl.get())
        _set_mask(self.bottom_right, br.get())
        _set_mask(self.top_left, tl.get())
        _set_mask(self.top_right, tr.get())
