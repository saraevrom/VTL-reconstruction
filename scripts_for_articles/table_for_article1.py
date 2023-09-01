#!/usr/bin/env python3
import argparse, os, json, sys
import re
import h5py
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("source_dir", type=str)

FILENAME_REGEX = re.compile(r"(F\d+)-(\w+)\.h5")
FOCUS = 154.3293791839621

def intersect_test(iter_1, iter_2):
    return len(set(iter_1).intersection(set(iter_2)))>1


def get_parameter_estimation(summary, key):
    return (summary[summary["parameter"]==key]["median"]).iloc[0]

def load_conf(srcdir):
    srcfile = os.path.join(srcdir,"config.json")
    if os.path.isfile(srcfile):
        with open(srcfile,"r") as fp:
            obj = json.load(fp)
    else:
        obj = dict()
    return obj


TABLE_MATCH = {
    "LIN":"lin",
    "CONST":"const",
    "EXP":"expo",
    "EXPOLINEAR":"expo-lin",
    "MIXED": "FIXME!!",
    "GAUSS":"gauss"
}

def sort_key_for_filename(fn):
    mat = FILENAME_REGEX.match(fn)
    if mat:
        return int(mat.groups()[0][1:])
    return -1


if __name__=="__main__":
    args = parser.parse_args()
    src_dir = args.source_dir
    conf = load_conf(src_dir)
    if not os.path.exists(src_dir):
        print("Source directory does not exist", file=sys.stderr)
        exit(1)


    # Prepare table header
    table = "\\begin{tabular}{|l|*{8}{c|}}\n"
    table += "Event & LC & UTC & ($\\gamma$, $\\gamma_\mathrm{sh}$), $\\,^\\circ$  & $\\alpha$, $\\,^\\circ$ & $\\Delta\\Phi$, $\\,^\\circ$ & $U_0$, px/fr & $\\sigma_\\mathrm{psf}$, px & $\\tau$, fr \\\\\n"
    table += "\\hline\n"
    rows = []

    source_filenames = list(filter(lambda x: x.endswith(".h5"), os.listdir(src_dir)))
    source_filenames.sort(key=sort_key_for_filename)
    for filename in source_filenames:
        fnmatch = FILENAME_REGEX.match(filename)
        if fnmatch:
            event, lc = fnmatch.groups()
            num = int(event[1:])
            filepath = os.path.join(src_dir,filename)
            with h5py.File(filepath, "r") as h5file:
                # Shower radiant position
                astro_tgt_x = h5file.attrs["astro_target_x"]
                astro_tgt_y = h5file.attrs["astro_target_y"]

                reco_keys = ["RECO_" + mode for mode in "MABCD"]
                needs_diff = intersect_test(reco_keys, h5file.keys())
                for mode in "MABCD":
                    reco_key = "RECO_" + mode

                    # Get event label
                    if needs_diff:
                        event_disp = f"{event}-{mode}"
                    else:
                        event_disp = event

                    if reco_key in h5file.keys():
                        # Get gamma/phi angles
                        summary = pd.read_hdf(filepath, mode="r", key=reco_key + "/summary")
                        x0 = get_parameter_estimation(summary, "X0")
                        y0 = get_parameter_estimation(summary, "Y0")
                        phi0 = get_parameter_estimation(summary, "Phi0")
                        u0 = get_parameter_estimation(summary, "U0")
                        sigma_psf = get_parameter_estimation(summary, "SigmaPSF")
                        parameters = np.array(summary["parameter"])
                        print(parameters)
                        if "τ_LC" in parameters:
                            tau = get_parameter_estimation(summary, "τ_LC")
                            tau = f"{tau:.3f}"
                        elif "τ_L" in parameters:
                            tau1 = get_parameter_estimation(summary, "τ_L")
                            tau2 = get_parameter_estimation(summary, "τ_R")
                            tau = f"{tau1:.3f}/{tau2:.3f}"
                        else:
                            tau = ""

                        phi_sh = np.arctan2(y0-astro_tgt_y,x0-astro_tgt_x)*180/np.pi

                        dot = (x0*astro_tgt_x+y0*astro_tgt_y+FOCUS**2)
                        norm0 = (x0**2+y0**2+FOCUS**2)**0.5
                        norm_astro = (astro_tgt_x**2+astro_tgt_y**2+FOCUS**2)**0.5
                        cos_alpha = dot/(norm0*norm_astro)
                        alpha = np.arccos(cos_alpha) * 180 / np.pi
                        alpha = round(alpha, 3)

                        cos_gamma = FOCUS/norm0
                        gamma = np.arccos(cos_gamma)*180/np.pi

                        cos_gamma_sh = FOCUS/norm_astro
                        gamma_sh = np.arccos(cos_gamma_sh)*180/np.pi


                        delta_phi = phi_sh-phi0
                        # delta_phi = round(delta_phi,3)
                        # u0 = round(u0,3)
                        # sigma_psf = round(sigma_psf,3)
                        utc = conf[filename[:-3]]["utc"]
                        # gamma = round(gamma,3)
                        # gamma_sh = round(gamma_sh,3)

                        # Add new row
                        table += f"{event_disp} & {TABLE_MATCH[lc]} & {utc} & "
                        table += f"({gamma:.3f},{gamma_sh:.3f}) & {alpha:.3f} & {delta_phi:.3f} & {u0:.3f} & "
                        table += f"{sigma_psf:.3f} & "
                        if tau:
                            table += tau+" "
                        table += f"\\\\\n"
            table += "\\hline\n"

    # Finish table
    for row in rows:
        table+= row[1]
    table += "\\end{tabular}"
    print(table)