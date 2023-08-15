#!/usr/bin/env python3
import os, sys, json, argparse
import h5py
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from matplotlib.patches import Rectangle
import matplotlib
import pandas as pd

matplotlib.rc('font', size=30)


HALF_PIXELS = 8
PIXEL_SIZE = 2.85
HALF_GAP_SIZE = 2.0

LOWER_EDGES = np.arange(HALF_PIXELS)*PIXEL_SIZE+HALF_GAP_SIZE
LOWER_EDGES = np.concatenate([-np.flip(LOWER_EDGES)-PIXEL_SIZE, LOWER_EDGES])

SPAN = HALF_PIXELS * PIXEL_SIZE + HALF_GAP_SIZE

RECO_POOL = [
    dict(linestyle="--", color="red"),
    dict(linestyle="dashdot", color="red"),
]

parser = argparse.ArgumentParser()
parser.add_argument("source_dir", type=str)
parser.add_argument("dest_dir", type=str)
parser.add_argument('--stack_plot', action='store_true')

def load_conf(srcdir):
    srcfile = os.path.join(srcdir,"config.json")
    if os.path.isfile(srcfile):
        with open(srcfile,"r") as fp:
            obj = json.load(fp)
    else:
        obj = dict()
    return obj

def postprocess_conf(conf_dict, key, fig, ax):
    if key in conf_dict.keys():
        entry:dict
        entry = conf_dict[key]
        if "xlim" in entry.keys():
            ax.set_xlim(*entry["xlim"])
        if "ylim" in entry.keys():
            ax.set_ylim(*entry["ylim"])

def optional_value(conf_dict, main_key, aux_key, default):
    if main_key in conf_dict.keys():
        entry: dict
        entry = conf_dict[key]
        return entry.get(aux_key, default)
    return default

def create_grid(axes):
    axes.vlines(LOWER_EDGES, -SPAN, -HALF_GAP_SIZE, colors="black")
    axes.vlines(LOWER_EDGES, SPAN, HALF_GAP_SIZE, colors="black")
    axes.vlines([-HALF_GAP_SIZE, SPAN], -SPAN, -HALF_GAP_SIZE, colors="black")
    axes.vlines([-HALF_GAP_SIZE, SPAN], SPAN, HALF_GAP_SIZE, colors="black")
    axes.hlines(LOWER_EDGES, -SPAN, -HALF_GAP_SIZE, colors="black")
    axes.hlines(LOWER_EDGES, SPAN, HALF_GAP_SIZE, colors="black")
    axes.hlines([-HALF_GAP_SIZE, SPAN], -SPAN, -HALF_GAP_SIZE, colors="black")
    axes.hlines([-HALF_GAP_SIZE, SPAN], SPAN, HALF_GAP_SIZE, colors="black")
    axes.set_xlim(-SPAN,SPAN)
    axes.set_ylim(-SPAN,SPAN)

def x0_from_pmt(alive_matrix):
    low = HALF_GAP_SIZE
    high = -HALF_GAP_SIZE
    fullsize=False
    if alive_matrix[:8,:].any():
        low = -SPAN
    if alive_matrix[8:,:].any():
        high = SPAN
        fullsize = low<0
    if high<low:
        fullsize = True
        low = -SPAN
        high = SPAN
    return low, high, fullsize

def y0_from_pmt(alive_matrix):
    low = HALF_GAP_SIZE
    high = -HALF_GAP_SIZE
    fullsize=False
    if alive_matrix[:, :8].any():
        low = -SPAN
    if alive_matrix[:, 8:].any():
        fullsize = low<0
        high = SPAN
    if high < low:
        fullsize = True
        low = -SPAN
        high = SPAN
    return low, high, fullsize

# GRID BASED COLORATION
def hsv_to_rgb(h,s,v):
    '''
    0<=h<=360
    0<=s<=1
    0<=v<=1
    '''
    h = h % 360
    c = v*s
    x = c*(1-abs((h/60)%2-1))
    m = v-c
    if h<60:
        _r = c; _g = x; _b = 0
    elif h<120:
        _r = x
        _g = c
        _b = 0
    elif h<180:
        _r = 0
        _g = c
        _b = x
    elif h<240:
        _r = 0
        _g = x
        _b = c
    elif h<300:
        _r = x
        _g = 0
        _b = c
    else:
        _r = c
        _g = 0
        _b = x
    return _r+m, _g+m, _b+m


def h_color(i, hue_shift=0.0,s_shift = 0.0, v_shift = 0.0):
    h = (i)/8*360+hue_shift
    s = 1-s_shift
    v = 1-v_shift
    return hsv_to_rgb(h,s,v)

WIDTH = 2*HALF_PIXELS
HEIGHT = 2*HALF_PIXELS


def floormod(x,y):
    pivot = int(np.floor(x/y))*y
    return x-pivot

def get_color(i,j):
    if i%2==0:
        j1 = j
    else:
        j1 = j + 1
    shift_id = floormod(floormod(i-j1*WIDTH//4,WIDTH),WIDTH)
    gray_shift = 0.0
    if j%2==0 and (i-j//2)%2==0:
        gray_shift = 1.0
    return h_color(shift_id,j/HEIGHT*180,
                   v_shift=gray_shift*0.5,
                   s_shift=gray_shift*0.3)

def get_parameter_estimation(summary, key):
    return (summary[summary["parameter"]==key]["median"]).iloc[0]

def intersect_test(iter_1, iter_2):
    return len(set(iter_1).intersection(set(iter_2)))>1


def argsort_metric(x:np.ndarray):
    slide = sliding_window_view(x, axis=1,window_shape=3)
    smoothed = slide.mean(axis=-1)
    return np.argmax(smoothed, axis=1)

if __name__=="__main__":
    args = parser.parse_args()
    src_dir = args.source_dir
    dest_dir = args.dest_dir
    if not os.path.exists(src_dir):
        print("Source directory does not exist", file=sys.stderr)
        exit(1)
    conf = load_conf(src_dir)
    os.makedirs(dest_dir, exist_ok=True)
    pdfs = os.path.join(dest_dir, "pdf")
    os.makedirs(pdfs, exist_ok=True)
    pngs = os.path.join(dest_dir, "png")
    os.makedirs(pngs, exist_ok=True)
    source_filenames = filter(lambda x:x.endswith(".h5"),os.listdir(src_dir))
    for filename in source_filenames:
        filepath = os.path.join(src_dir,filename)
        with h5py.File(filepath,"r") as h5file:
            print("Processing", filepath)
            xs = np.array(h5file["time_s"])
            origin_xs = np.array(h5file["x_data"])
            ys = np.array(h5file["y_data"])
            selection = np.array(h5file["selection"])

            fig, ax = plt.subplots(figsize=(16, 12))

            ax.set_xlabel("Time, s")
            ax.set_ylabel("Signal, a.u.")

            key = filename[:-3]
            inv_legend = optional_value(conf, key, "inv_legend", False)
            if inv_legend:
                axin1 = ax.inset_axes([0.70, 0.65, 0.35, 0.35])
            else:
                axin1 = ax.inset_axes([0.00, 0.65, 0.35, 0.35])
            axin1.get_xaxis().set_visible(False)
            axin1.get_yaxis().set_visible(False)
            axin1.set_aspect("equal")
            create_grid(axin1)
            xmin, xmax, fs1 = x0_from_pmt(selection)
            ymin, ymax, fs2 = y0_from_pmt(selection)
            if fs1 or fs2:
                xmin, xmax = -SPAN, SPAN
                ymin, ymax = -SPAN, SPAN

            axin1.set_xlim(xmin, xmax)
            axin1.set_ylim(ymin, ymax)

            stackplot_data = []
            stackplot_colors = []
            # PLOT PIXELMAP
            for i, x in enumerate(LOWER_EDGES):
                for j, y in enumerate(LOWER_EDGES):
                    if selection[i,j]:
                        if not args.stack_plot:
                            line, = ax.plot(xs, ys[:, i, j], color=get_color(i,j), linewidth=1.0)
                        stackplot_data.append(ys[:, i, j])
                        stackplot_colors.append(get_color(i,j))
                        rect = Rectangle((x, y), PIXEL_SIZE, PIXEL_SIZE, color=get_color(i,j))
                        axin1.add_patch(rect)

            if args.stack_plot:
                stackplot_colors = np.array(stackplot_colors)
                stackplot_data = np.vstack(stackplot_data)

                # sort stack histogram data
                argmaxed = argsort_metric(stackplot_data)
                order = np.argsort(argmaxed)
                stackplot_data = stackplot_data[order]
                stackplot_colors = stackplot_colors[order]

                # PLOT stack histogram
                ax.stackplot(xs, stackplot_data, colors= stackplot_colors)

            # PLOT LC
            accum = np.sum(ys[:,selection], axis=1)

            ax.plot(xs,accum, color="black", label="LC", linewidth=2)

            reco_i = 0
            x_com = 0
            y_com = 0
            com_cnt = 0
            # PLOT RECO
            reco_keys = ["RECO_"+mode for mode in "MABCD"]
            needs_diff = intersect_test(reco_keys, h5file.keys())
            for mode in "MABCD":
                reco_key = "RECO_"+mode
                if reco_key in h5file.keys():
                    reco = h5file[reco_key]
                    reco_xs = np.array(reco["lc_xs"]) # in frames
                    reco_start = reco_xs[0]
                    reco_end = reco_xs[-1]
                    print("RECO", reco_start, reco_end)
                    k0 = (reco_start+reco_end)/2
                    reco_xs = np.interp(reco_xs, origin_xs, xs) # Transform to seconds

                    reco_ys = np.array(reco["lc_ys"])
                    if needs_diff:
                        legend_label = f"Reco LC "+mode
                    else:
                        legend_label = f"Reco LC"
                    ax.plot(reco_xs, reco_ys, label=legend_label, linewidth=2, **RECO_POOL[reco_i])
                    reco_i +=1
                    summary = pd.read_hdf(filepath, mode="r",key=reco_key+"/summary")
                    print(filename, mode)
                    print(summary)
                    phi0 = get_parameter_estimation(summary, "Phi0")*np.pi/180
                    x0 = get_parameter_estimation(summary, "X0")
                    y0 = get_parameter_estimation(summary, "Y0")
                    u0 = get_parameter_estimation(summary, "U0")
                    sigma_psf = get_parameter_estimation(summary, "SigmaPSF")
                    ts = np.array([reco_start, reco_end])-k0
                    reco_overlay_x = PIXEL_SIZE*u0*ts*np.cos(phi0) + x0
                    reco_overlay_y = PIXEL_SIZE*u0*ts*np.sin(phi0) + y0
                    dx = reco_overlay_x[1]-reco_overlay_x[0]
                    dy = reco_overlay_y[1]-reco_overlay_y[0]
                    axin1.arrow(x=reco_overlay_x[0], dx=dx,
                                y=reco_overlay_y[0], dy=dy,
                                color="red", width=sigma_psf*PIXEL_SIZE,  length_includes_head=True)
                    astro_tgt_x = h5file.attrs["astro_target_x"]
                    astro_tgt_y = h5file.attrs["astro_target_y"]
                    x_com += x0
                    y_com += y0
                    com_cnt += 1
            if com_cnt>0:
                x_com /= com_cnt
                y_com /= com_cnt
                dx = astro_tgt_x - x_com
                dy = astro_tgt_y - y_com
                axin1.arrow(x=x_com, dx=dx,
                            y=y_com, dy=dy,
                            color="green", width=sigma_psf * PIXEL_SIZE/4, length_includes_head=True)


            if inv_legend:
                ax.legend(loc="upper center")
            else:
                ax.legend(loc="upper right")
            postprocess_conf(conf,key,fig,ax)
            fig.tight_layout()
            fig.savefig(os.path.join(pngs, filename[:-3])+".png", dpi=250)
            fig.savefig(os.path.join(pdfs, filename[:-3])+".pdf")
