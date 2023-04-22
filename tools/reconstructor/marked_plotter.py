
from vtl_common.localized_GUI import GridPlotter
from vtl_common.parameters import HALF_PIXELS, PIXEL_SIZE, HALF_GAP_SIZE
from matplotlib.patches import Rectangle

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
