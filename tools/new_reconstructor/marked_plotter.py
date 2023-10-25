import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Rectangle

from fixed_rotator import Vector2
from vtl_common.localized_GUI import GridPlotter
from vtl_common.parameters import HALF_PIXELS, PIXEL_SIZE, HALF_GAP_SIZE
from vtl_common.parameters import MAIN_LATITUDE, MAIN_LONGITUDE

# from tools.orientation.orientation.stellar_math import radec_to_eci, eci_to_ocef
# from tools.orientation.orientation.stellar_math import ocef_to_detector_plane, unixtime_to_era, rotate_yz
# from tools.orientation.orientation.stellar_math import ocef_to_altaz

from fixed_rotator.astro_math_z_aligned import radec_to_eci, eci_to_ocef, ocef_to_detector_plane, unixtime_to_era, Quaternion
from fixed_rotator.astro_math_z_aligned import ocef_to_altaz, radec_to_ocef, Vector3
from fixed_rotator.astro_math_z_aligned import ocef_to_eci, ocef_to_radec, detector_plane_to_ocef_f
from specific_ui.data_output import DataOutput


SPAN = HALF_PIXELS*PIXEL_SIZE+HALF_GAP_SIZE
NA_TEXT = "N/A"+"\n"*4

def _set_mask(patch, value):
    if value:
        patch.set_alpha(0.0)
    else:
        patch.set_alpha(0.5)


class PlotProxy(object):
    def __init__(self, axes):
        self.axes = axes
        self.__proxy_origin = None

    def append_patch(self, patch):
        pass

    def plot_lines(self, *args, **kwargs):
        lines = self.axes.plot(*args,**kwargs)
        for l in lines:
            self.append_patch(l)

    def plot_arrow(self, *args, **kwargs):
        arrow = self.axes.arrow(*args, **kwargs)
        self.append_patch(arrow)

    def plot_circle(self, *args, **kwargs):
        patch = plt.Circle(*args, **kwargs)
        self.axes.add_patch(patch)
        self.append_patch(patch)

    def set_origin(self, x, y, k0=None, view_enable=False, writer=None):
        pass


ARROW_ARGS = dict(width=0.2, length_includes_head=True, edgecolor="green", facecolor="green")


def detector_point_to_radec(orientation, direction_ocef, era):
    v_lat = orientation["VIEW_LATITUDE"] * np.pi / 180
    v_lon = orientation["VIEW_LONGITUDE"] * np.pi / 180
    self_rotation = orientation["SELF_ROTATION"] * np.pi / 180

    ra,dec = ocef_to_radec(direction_ocef,v_lat, v_lon,self_rotation,era)
    return ra,dec


def write_to(writer:DataOutput,label,text):
    if writer is not None:
        writer.add_entry(label,text)


def add_sep(writer:DataOutput,label,size=14):
    if writer is not None:
        writer.add_separator(label,size)


def write_radec(writer, ra,dec, era, orientation):
    add_sep(writer, "SKY", 10)
    write_to(writer, "DEC [°]", f"{dec * 180 / np.pi:.2f}")
    write_to(writer, "RA [°]", f"{ra * 180 / np.pi:.2f}")
    write_to(writer, "RA [h]", f"{ra * 12 / np.pi:.2f}")

    v_lat = orientation["VIEW_LATITUDE"] * np.pi / 180
    v_lon = orientation["VIEW_LONGITUDE"] * np.pi / 180
    self_rotation = orientation["SELF_ROTATION"] * np.pi / 180
    f = orientation["FOCAL_DISTANCE"]

    v_ocef = radec_to_ocef(ra, dec, v_lat, v_lon, self_rotation, era)
    xy, v = ocef_to_detector_plane(v_ocef, f)
    x, y = xy.x, xy.y

    angsep = np.arctan((x * x + y * y) ** 0.5 / f) * 180 / np.pi
    if not v:
        angsep = 180 - angsep

    add_sep(writer, "DEVICE", 10)
    write_to(writer, "γ [°]", f"{angsep:.2f}")

    psi = np.arctan2(y, x) * 180 / np.pi
    write_to(writer,"ψ [°]",f"{psi:.2f}")
    v_obs = radec_to_ocef(ra, dec,
                          MAIN_LATITUDE * np.pi / 180,
                          MAIN_LONGITUDE * np.pi / 180,
                          0.0, era)

    alt, az = ocef_to_altaz(v_obs)
    alt *= 180 / np.pi
    az *= 180 / np.pi

    zang = 90 - alt
    add_sep(writer, "HORIZON", 10)
    write_to(writer,"θ [°]", f"{zang:.2f}")
    write_to(writer,"φ [°]", f"{az:.2f}")

class HighlightingPlotter(GridPlotter, PlotProxy):
    def __init__(self, master, controller):
        GridPlotter.__init__(self,master)
        PlotProxy.__init__(self,self.axes)
        self._controller = controller
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
        self._point_target = (1, 1)
        self._pointer_time_index = 0  # Frame from start for ut0
        self._direction_arrow = self.axes.arrow(0, 0, 1, 1, **ARROW_ARGS)
        self._direction_arrow.set_visible(False)
        self._not_visible = self.axes.text(0, 0, "N/V", color="red",ha='center', va='center')
        self._not_visible.set_visible(False)
        self._pointer_data = None
        self._times = None


    def set_point_direction(self, formdata, ut0_arr, writer=None):
        self._pointer_data = formdata
        self._times = ut0_arr
        return self.update_point(writer)

    def _write_origin(self, writer,origin):
        formdata = self._pointer_data
        if formdata["show"]:
            ra = formdata["ra"]
            dec = formdata["dec"]
            print("RA=", ra)
            orientation = formdata["orientation"]
            print(self._pointer_time_index)
            ut0 = self._times[self._pointer_time_index]
            era = unixtime_to_era(ut0)
            v_lat = orientation["VIEW_LATITUDE"] * np.pi / 180
            v_lon = orientation["VIEW_LONGITUDE"] * np.pi / 180
            self_rotation = orientation["SELF_ROTATION"] * np.pi / 180
            f = orientation["FOCAL_DISTANCE"]

            origin_dir = Vector3(origin[0], origin[1], f).normalized()
            origin_ocef = detector_plane_to_ocef_f(Vector2(origin[0], origin[1]), f)
            ra_origin, dec_origin = detector_point_to_radec(orientation, origin_ocef, era)
            if ra_origin < 0:
                ra_origin += np.pi * 2


            write_radec(writer, ra_origin, dec_origin, era, orientation)

            v_ocef = radec_to_ocef(ra, dec, v_lat, v_lon, self_rotation, era)
            xy, v = ocef_to_detector_plane(v_ocef, f)
            x, y = xy.x, xy.y

            tgtdir = Vector3(x, y, f).normalized()
            if not v:
                tgtdir = -tgtdir

            x1 = x - self._origin[0]
            y1 = y - self._origin[1]

            angsep1 = np.arccos(tgtdir.dot(origin_dir)) * 180 / np.pi
            psi1 = np.arctan2(y1, x1) * 180 / np.pi

            write_to(writer, "γ (relative) [°]", f"{round(angsep1, 2)}")
            write_to(writer, "ψ (relative) [°]", f"{round(psi1, 2)}")
            return psi1

    def update_point(self, writer=None):
        if self._pointer_data is not None:
            formdata = self._pointer_data
            if formdata["show"]:
                ra = formdata["ra"]
                dec = formdata["dec"]
                print("RA=", ra)
                orientation = formdata["orientation"]
                print(self._pointer_time_index)
                ut0 = self._times[self._pointer_time_index]
                era = unixtime_to_era(ut0)
                v_lat = orientation["VIEW_LATITUDE"] * np.pi / 180
                v_lon = orientation["VIEW_LONGITUDE"] * np.pi / 180
                self_rotation = orientation["SELF_ROTATION"] * np.pi / 180
                f = orientation["FOCAL_DISTANCE"]

                v_ocef = radec_to_ocef(ra, dec, v_lat, v_lon, self_rotation, era)
                xy, v = ocef_to_detector_plane(v_ocef, f)
                x,y = xy.x, xy.y
                print(self._origin)

                if v:
                    self._not_visible.set_visible(False)
                    self.point_to(x, y)
                else:
                    self.set_nv()

                add_sep(writer,"STELLAR")
                #s = "STELLAR\n"
                #s += display_radec(ra,dec,era, orientation)
                write_radec(writer,ra,dec,era, orientation)

                add_sep(writer,"START POINT")
                psi1 = self._write_origin(writer,self._origin)
                recos = self._controller.get_current_reconstruction()
                for reco in recos:
                    mode, obj = reco
                    add_sep(writer,f"RECO {mode}")
                    phi0 = obj.ask_parameter("phi0")
                    x0 = obj.ask_parameter("x0")
                    y0 = obj.ask_parameter("y0")
                    if x0 is None or y0 is None:
                        psi1_use = psi1
                    else:
                        psi1_use = self._write_origin(writer,(x0,y0))

                    if phi0 is not None:
                        delta_psi = (180+psi1_use)-phi0
                        if delta_psi>360:
                            delta_psi = delta_psi%360
                        if delta_psi>180:
                            delta_psi -= 360
                        if phi0 is not None:
                            write_to(writer, f"Δψ [°]", f"{delta_psi:.2f}")
                return
            else:
                self.hide_pointer()
        add_sep(writer,NA_TEXT)

    def set_nv(self):
        self._direction_arrow.set_visible(False)
        self._not_visible.set_visible(True)
        self.draw()

    def hide_pointer(self):
        self._direction_arrow.set_visible(False)
        self._not_visible.set_visible(False)
        self.draw()

    def set_origin(self, x, y, k0=None, view_enable=False, writer=None):
        x = float(x)
        y = float(y)
        self._origin = (x, y)
        if k0 is not None:
            self._pointer_time_index = int(k0)
        return self.update_point(writer)

    def point_to(self,x,y, view_enable=True):
        self._point_target = (x,y)
        self.update_point_arrow(view_enable)

    def get_pointer_data(self):
        return self._origin, self._point_target

    def update_point_arrow(self, view_enable=False):
        x0 = self._origin[0]
        y0 = self._origin[1]
        x = self._point_target[0]
        y = self._point_target[1]
        self._direction_arrow.set_data(x=x0, y=y0, dx=x-x0, dy=y-y0)
        if view_enable:
            self._direction_arrow.set_visible(True)
        self.draw()

    def append_patch(self, patch):
        self.added_patches.append(patch)

    def clear_added_patches(self):
        for x in self.added_patches:
            x.remove()
        self.added_patches.clear()
        self.set_origin(0, 0, 0)

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

    def mirror_arrow_direction(self, proxy:PlotProxy):
        x0 = self._origin[0]
        y0 = self._origin[1]
        x = self._point_target[0]
        y = self._point_target[1]
        if self._direction_arrow.get_visible():
            proxy.axes.arrow(x0, y0, x-x0, y-y0, **ARROW_ARGS)