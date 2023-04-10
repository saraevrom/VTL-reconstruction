from vtl_common.localized_GUI import Plotter
import numpy as np
from orientation.database_reader import get_database
from vtl_common.parameters import MAIN_LATITUDE, MAIN_LONGITUDE, MAX_STAR_MAGNITUDE
from vtl_common.parameters import HALF_PIXELS, PIXEL_SIZE, HALF_GAP_SIZE
from orientation.stellar_math import unixtime_to_era, eci_to_ocef, ocef_to_altaz, detector_plane_to_ocef
from orientation.stellar_math import rotate_yz, rotate_xz, rotate_xy

DETECTOR_SPAN = HALF_PIXELS*PIXEL_SIZE+HALF_GAP_SIZE

MIN_STAR_MAGNITUDE =-1.5


def params_rotate(xs,ys,  params):
    focus = params["FOCAL_DISTANCE"]
    x_ocef, y_ocef, z_ocef = detector_plane_to_ocef(xs, ys, focus)
    self_rotation = params["SELF_ROTATION"] * np.pi / 180
    lat = params["VIEW_LATITUDE"]
    lon = params["VIEW_LONGITUDE"]

    x_ocef, y_ocef, z_ocef = rotate_yz(x_ocef, y_ocef, z_ocef, -self_rotation)
    x_ocef, y_ocef, z_ocef = rotate_xz(x_ocef, y_ocef, z_ocef, lat * np.pi / 180)
    x_ocef, y_ocef, z_ocef = rotate_xy(x_ocef, y_ocef, z_ocef, (lon - MAIN_LONGITUDE) * np.pi / 180)
    x_ocef, y_ocef, z_ocef = rotate_xz(x_ocef, y_ocef, z_ocef, -MAIN_LATITUDE * np.pi / 180)
    alt, az = ocef_to_altaz(x_ocef, y_ocef, z_ocef)
    rs = 90 - alt * 180 / np.pi
    return az, rs

class SkyPlotter(Plotter):
    def __init__(self, master):
        super().__init__(master, polar=True)
        self.axes.set_theta_zero_location('N')
        self.axes.set_xticks([0, np.pi/2, np.pi, np.pi*1.5])
        self.axes.set_xticklabels(["N","E","S","W"])
        self.axes.set_ylim(0, 90.0)
        self._stars = None
        self._fov = None
        self._xarrow = None
        self._yarrow = None
        self._starnames = None
        self._annot = self.axes.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        self._annot.set_visible(False)
        self.figure.canvas.mpl_connect("motion_notify_event", self._hover)
        self.figure.canvas.mpl_connect("axes_leave_event", self._leave)

    def _update_annot(self, ind):
        #print("IND:", ind)
        pos = self._stars.get_offsets()[ind["ind"][0]]
        self._annot.xy = pos
        text = self._starnames.iloc[ind["ind"][0]]
        self._annot.set_text(text)
        #self._annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        self._annot.get_bbox_patch().set_alpha(0.4)

    def _hover(self, event):
        if self._stars is not None:
            vis = self._annot.get_visible()
            if event.inaxes == self.axes:
                cont, ind = self._stars.contains(event)
                if cont:
                    self._update_annot(ind)
                    self._annot.set_visible(True)
                    self.figure.canvas.draw_idle()
                else:
                    if vis:
                        self._annot.set_visible(False)
                        self.figure.canvas.draw_idle()

    def _leave(self, event):
        self._annot.set_visible(False)
        self.figure.canvas.draw_idle()

    def plot_stars(self, unixtime):
        era = unixtime_to_era(unixtime)
        database = get_database()
        eci_x = database["eci_x"].to_numpy()
        eci_y = database["eci_y"].to_numpy()
        eci_z = database["eci_z"].to_numpy()
        x,y,z = eci_to_ocef(eci_x, eci_y, eci_z, era, MAIN_LATITUDE*np.pi/180, MAIN_LONGITUDE*np.pi/180)
        alt, az = ocef_to_altaz(x,y,z)
        if self._stars is None:
            mag = database["mag"]*1.0
            y1 = -2.5
            y2 = 0.0
            x1 = MIN_STAR_MAGNITUDE
            x2 = 5.0
            k = (y2-y1)/(x2-x1)
            b = y1 - k*x1
            self._starnames = database["star_name"]
            self._stars = self.axes.scatter(az, 90 - 180 * alt / np.pi, s=10 ** -(k * mag + b))
        else:
            offsets = np.vstack([az,90-180*alt/np.pi]).T
            self._stars.set_offsets(offsets)
        self.draw()

    def draw_fov(self, params):
        spanned_range = np.linspace(-DETECTOR_SPAN, DETECTOR_SPAN, 100)
        constant_range = np.full(shape=(100,), fill_value=DETECTOR_SPAN)

        dir_i_x = np.array([0.0, DETECTOR_SPAN])
        dir_i_y = np.array([0.0, 0.0])

        dir_j_x = np.array([0.0, 0.0])
        dir_j_y = np.array([0.0, DETECTOR_SPAN])

        xs = np.concatenate([constant_range, -spanned_range, -constant_range, spanned_range])
        ys = np.concatenate([spanned_range, constant_range, -spanned_range, -constant_range])

        az, rs = params_rotate(xs,ys,params)
        if self._fov is None:
            self._fov, = self.axes.plot(az, rs, "-", color="black")
        else:
            self._fov.set_xdata(az)
            self._fov.set_ydata(rs)

        azx, rsx = params_rotate(dir_i_x, dir_i_y, params)
        if self._xarrow is None:
            self._xarrow = self.axes.arrow(x=azx[0], y=rsx[0], dx=azx[1]-azx[0], dy=rsx[1]-rsx[0], color="red", width=0.1)
        else:
            self._xarrow.set_data(x=azx[0], y=rsx[0], dx=azx[1]-azx[0], dy=rsx[1]-rsx[0])

        azy, rsy = params_rotate(dir_j_x, dir_j_y, params)
        if self._yarrow is None:
            self._yarrow = self.axes.arrow(x=azy[0], y=rsy[0], dx=azy[1]-azy[0], dy=rsy[1]-rsy[0], color="blue", width=0.1)
        else:
            self._yarrow.set_data(x=azy[0], y=rsy[0], dx=azy[1]-azy[0], dy=rsy[1]-rsy[0])
        self.draw()


