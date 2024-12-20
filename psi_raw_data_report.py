#!/usr/bin/env python3
"""

called from psi.R so as to present figures and create a raw data report

"""

__author__ = "Andy Schauer"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-12-19"
__version__ = "1.0"
__copyright__ = "Copyright 2024, Andy Schauer"
__license__ = "Apache 2.0"


# -------------------- Imports --------------------
import argparse
from bokeh.plotting import figure, output_file, save, show
from bokeh.embed import file_html
from bokeh.io import export_svgs
from bokeh import palettes
from bokeh.resources import CDN  # , INLINE
from bokeh.models import Legend
#import cairosvg
import csv
import datetime as dt
import json
import isolab_lib as lab 
import matplotlib.pyplot as pplt
from natsort import natsorted
import numpy as np
import os
import shutil
import sys
import time
import warnings
import webbrowser


# -------------------- Functions ------------------




# -------------------- Get Session ------------------
parser = argparse.ArgumentParser()
parser.add_argument('current_session', type=str, help="This allows R to pass the session name to python.")
args = parser.parse_args()
# print(f"\n\n    {args.current_session}")



# -------------------- paths --------------------
project_path = lab.get_path("psi", "project")



# ---------- Identify session to be loaded ----------
session = args.current_session
session_path = f'{project_path}/{session}/'
report_path = f"{session_path}raw_data_report/"
archive_path = f"{session_path}archive/"

session_samplelog = f"psi_{session}_samplelog.csv"
# session_samplelog = f"psi_{session}_HEG_samplelog.csv"

headers, data = lab.read_file(os.path.join(session_path, session_samplelog), ',')

original_data = data.copy()
trust0_indices = [i for i, e in enumerate(data['trust']) if int(e) == 0]
trust1_indices = [i for i, e in enumerate(data['trust']) if int(e) == 1]
for header in headers:
    data[header] = [data[header][index] for index in trust1_indices]

numlist = ['UID', 'trust', 'd13C', 'd13Csd', 'd13Cn', 'd13Cslope', 'd13Cr2', 'd13Cpval', 'd18O', 'd18Osd', 'd18On', 'd18Oslope', 'd18Or2', 'd18Opval',
           'd45', 'd45sd', 'd45n', 'd45slope', 'd45r2', 'd45pval', 'd46', 'd46sd', 'd46n', 'd46slope', 'd46r2', 'd46pval',
           'd47', 'd47sd', 'd47n', 'd47slope', 'd47r2', 'd47pval', 'D47raw', 'D47sd', 'D47n', 'D47slope', 'D47r2', 'D47pval',
           'd48', 'd48sd', 'd48n', 'd48slope', 'd48r2', 'd48pval', 'D48raw', 'D48sd', 'D48n', 'D48slope', 'D48r2', 'D48pval',
           'd49', 'd49sd', 'd49n', 'd49slope', 'd49r2', 'd49pval', 'D49raw', 'D49sd', 'D49n', 'D49slope', 'D49r2', 'D49pval',
           'R45sam_mean', 'R45sam_sd', 'R45wg_mean', 'R45wg_sd', 'R45sam_wg', 'R46sam_mean', 'R46sam_sd', 'R46wg_mean', 'R46wg_sd', 'R46sam_wg',
           'R47sam_mean', 'R47sam_sd', 'R47wg_mean', 'R47wg_sd', 'R47sam_wg', 'R48sam_mean', 'R48sam_sd', 'R48wg_mean', 'R48wg_sd', 'R48sam_wg',
           'R49sam_mean', 'R49sam_sd', 'R49wg_mean', 'R49wg_sd', 'R49sam_wg', 
           'bkgnd_44', 'bkgnd_45', 'bkgnd_46', 'bkgnd_47', 'bkgnd_48', 'bkgnd_49', 'unix_time', 'mass', 'sample_pressure', 'yield_mbar_ug', 'residual_gas',
           'nchops', 'min_ref_beam', 'max_ref_beam', 'total_blocks', 'total_cycles', 'initial_sample_beam', 'pre_balance_sample_beam',
           'nu_balance', 'sam_max44nA', 'sam_mean44nA', 'sam_min44nA', 'sam_max47nA', 'sam_mean47nA', 'sam_min47nA', 'min_balance', 'mean_balance', 'max_balance',
           'sam_decay', 'wg_decay']

strlist = set(headers) - set(numlist)

for i in numlist:
    data[i] = [x if x != "NA" else np.nan for x in data[i]]
    globals()[i] = np.asarray(data[i], dtype=float)

for i in strlist:
    globals()[i] = np.asarray(data[i])



ETH1 = [i for i, e in enumerate(Sample) if e=='ETH-1']
ETH2 = [i for i, e in enumerate(Sample) if e=='ETH-2']
ETH3 = [i for i, e in enumerate(Sample) if e=='ETH-3']
ETH4 = [i for i, e in enumerate(Sample) if e=='ETH-4']
IAEAC1 = [i for i, e in enumerate(Sample) if e=='IAEA-C1']
IAEAC2 = [i for i, e in enumerate(Sample) if e=='IAEA-C2']
Merck = [i for i, e in enumerate(Sample) if e=='MERCK']
GU1 = [i for i, e in enumerate(Sample) if e=='GU1']


carbs = [i for i, e in enumerate(material_type) if e=="CarbStd"]
pressure_vs_mass_line = "transducer_pressure = 0.03791 * mass - 3.99629"
#mass_for_line = np.array([np.nanmin(mass), np.nanmax(mass)])
mass_for_line = np.array([0, 10000])
predicted_pressure = 0.03791 * mass_for_line - 3.99629




# sample summary
samples_summary = {}
samples_summary['name'] = list(set(Sample))
[samples_summary['name'].remove(i) for i in ['ETH-1', 'ETH-2', 'ETH-3', 'ETH-4', 'IAEA-C1', 'IAEA-C2', 'MERCK', 'GU1']]
samples_summary['n'] = np.zeros(len(samples_summary['name']))
samples_summary['mass'] = np.zeros(len(samples_summary['name']))
samples_summary['mass_sd'] = np.zeros(len(samples_summary['name']))
samples_summary['mass_max'] = np.zeros(len(samples_summary['name']))
samples_summary['proportional_carbonate'] = np.zeros(len(samples_summary['name']))
samples_summary['proportional_carbonate_sd'] = np.zeros(len(samples_summary['name']))
samples_summary['target_mass'] = np.zeros(len(samples_summary['name']))
samples_summary['too_small'] = np.zeros(len(samples_summary['name']))

for i,sample in enumerate(samples_summary['name']):
    samples_summary['n'][i] = len(mass[Sample == sample])
    samples_summary['mass'][i] = np.round(np.mean(mass[Sample == sample]))
    samples_summary['mass_sd'][i] = np.round(2*np.std(mass[Sample == sample]))
    samples_summary['mass_max'][i] = np.max(mass[Sample == sample])
    samples_summary['proportional_carbonate'][i] = np.round(np.mean((sample_pressure[Sample == sample] + 3.99629) / 0.03791 / mass[Sample == sample]),3)
    samples_summary['proportional_carbonate_sd'][i] = np.round(np.std((sample_pressure[Sample == sample] + 3.99629) / 0.03791 / mass[Sample == sample]),3)
    samples_summary['target_mass'][i] = np.round(4500 / samples_summary['proportional_carbonate'][i])
    if samples_summary['mass_max'][i] < samples_summary['target_mass'][i]:
        samples_summary['too_small'][i] = 1

from tabulate import tabulate

# Assuming samples_summary has been defined as in your code
# Prepare data for tabulation
table_data = []
for i, name in enumerate(samples_summary['name']):
    table_data.append([
        name,
        samples_summary['n'][i],
        samples_summary['mass'][i],
        samples_summary['mass_sd'][i],
        samples_summary['mass_max'][i],
        samples_summary['proportional_carbonate'][i],
        samples_summary['proportional_carbonate_sd'][i],
        samples_summary['target_mass'][i],
        samples_summary['too_small'][i],
    ])

# Define headers
headers = ["Sample Name", "n", "Mass", "Mass SD", "Mass Max", "Prop. Carb.", "Pro Carb SD", "Target Mass", "Too Small"]

# Print the table
print(tabulate(table_data, headers=headers, tablefmt="grid"))














# -------------------- figures --------------------
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    print('\n\n        Making figures...')

    figures = {}
    font_size = '24pt'
    ms = 15
    fig_n = 1


    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. CO2 pressure vs material mass - Here are all recognized carbonate standards from the current session. The line was generated from a range of IAEA-C2 masses 
                                over several sessions. Any points that fall significantly below the line indicate a problem, perhaps explained by a clog or lack of acid.
                                Points above 2000 ug are expected to be below the line due to a non-linear response of pressure gauge. 
                                Use this equation to calculate the amount of carbonate in your samples."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="mass (ug)", y_axis_label="Sample pressure (mbar)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(mass[carbs], sample_pressure[carbs], size=10, line_color='black', fill_color='black')
    figures[fig_n]['fig'].line(mass_for_line, predicted_pressure, line_color='black')
    figures[fig_n]['fig'].text(x=0, y=60, text=[pressure_vs_mass_line], text_font_size='16pt')
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size

    fig_n += 1


    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. Sample Pressure - If we are weighing an consistent and appropriate amount of material, all of our sample pressures should be about the same. They
                                should also be chopped only once (chop = 1). """
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="Sample Pressure",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(UID[np.where(nchops==3)], sample_pressure[np.where(nchops==3)], legend_label="chops = 3",  size=10, line_color='black', fill_color='white')
    figures[fig_n]['fig'].circle(UID[np.where(nchops==2)], sample_pressure[np.where(nchops==2)], legend_label="chops = 2",  size=10, line_color='black', fill_color='cyan')
    figures[fig_n]['fig'].circle(UID[np.where(nchops==1)], sample_pressure[np.where(nchops==1)], legend_label="chops = 1",  size=10, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].circle(UID[np.where(nchops==0)], sample_pressure[np.where(nchops==0)], legend_label="chops = 0", size=10, line_color='black', fill_color='black')
    figures[fig_n]['fig'].circle(UID[np.where(nchops==-1)], sample_pressure[np.where(nchops==-1)], legend_label="chops = -1", size=10, line_color='black', fill_color='green')
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].legend.label_text_font_size = font_size
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')


    fig_n += 1


    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. Residual Gas - This is an estimate of the gases leftover after the reaction has completed and all of the CO2 and H2O have been frozen. Incondensibles
                                or leaks may cause this value to be high."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="Residual Gas",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(UID, residual_gas, size=10, line_color='black', fill_color='black')
    
    
    
    fig_n += 1



    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. CO2 Yield - This figure show the amount of CO2 obtained per microgram of carbonate. Only the standards are shown. """
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="Yield (mbar / ug)",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(UID[ETH1], yield_mbar_ug[ETH1], legend_label="ETH-1", size=ms, line_color='black', fill_color=(255, 51, 230))
    figures[fig_n]['fig'].circle(UID[ETH2], yield_mbar_ug[ETH2], legend_label="ETH-2", size=ms, line_color='black', fill_color=(255, 102, 153))
    figures[fig_n]['fig'].circle(UID[ETH3], yield_mbar_ug[ETH3], legend_label="ETH-3", size=ms, line_color='black', fill_color=(255, 153, 77))
    figures[fig_n]['fig'].circle(UID[ETH4], yield_mbar_ug[ETH4], legend_label="ETH-4", size=ms, line_color='black', fill_color=(255, 204, 0))
    figures[fig_n]['fig'].circle(UID[IAEAC1], yield_mbar_ug[IAEAC1], legend_label="IAEA-C1", size=ms, line_color='black', fill_color='white')
    figures[fig_n]['fig'].circle(UID[IAEAC2], yield_mbar_ug[IAEAC2], legend_label="IAEA-C2", size=ms, line_color='black', fill_color='tan')
    figures[fig_n]['fig'].circle(UID[Merck], yield_mbar_ug[Merck], legend_label="Merck", size=ms, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(UID[GU1], yield_mbar_ug[GU1], legend_label="GU1", size=ms, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')


    fig_n += 1



    bkgnd_44_index = bkgnd_44 / np.mean(bkgnd_44)
    bkgnd_45_index = bkgnd_45 / np.mean(bkgnd_45)
    bkgnd_46_index = bkgnd_46 / np.mean(bkgnd_46)
    bkgnd_47_index = bkgnd_47 / np.mean(bkgnd_47)
    bkgnd_48_index = bkgnd_48 / np.mean(bkgnd_48)

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. Normalized Background - This is the background on each faraday cup divided by the mean of all backgounds for that cup. """
    figures[fig_n]['fig'] = figure(y_range=(0,2), width=1400, height=700, x_axis_label="UID", y_axis_label="Backgrounds",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(UID, bkgnd_44_index, legend_label="m/z 44", size=10, line_color='black', fill_color='black')
    figures[fig_n]['fig'].circle(UID, bkgnd_45_index, legend_label="m/z 45", size=10, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].circle(UID, bkgnd_46_index, legend_label="m/z 46", size=10, line_color='black', fill_color='pink')
    figures[fig_n]['fig'].circle(UID, bkgnd_47_index, legend_label="m/z 47", size=10, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(UID, bkgnd_48_index, legend_label="m/z 48", size=10, line_color='black', fill_color='yellow')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')


    fig_n += 1


    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. Beam Size - This is the background on each faraday cup divided by the mean of all backgounds for that cup. """
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="m/z 44 signal (Amp)",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(UID, initial_sample_beam, legend_label="Initial Sample Beam", size=10, line_color='black', fill_color='black')
    figures[fig_n]['fig'].circle(UID, sam_max44nA, legend_label="Max Sample Beam during measurement", size=10, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(UID, sam_mean44nA, legend_label="Mean Sample Beam during measurement", size=10, line_color='black', fill_color='yellow')
    figures[fig_n]['fig'].circle(UID, sam_min44nA, legend_label="Min Sample Beam during measurement", size=10, line_color='black', fill_color='pink')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')



    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="D47 standard deviation",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(UID, D47sd, size=10, line_color='black', fill_color='black')

    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="Mean beam size (Amp)", y_axis_label="D47 standard deviation",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(sam_mean44nA, D47sd, size=10, line_color='black', fill_color='yellow')



    fig_n += 1


    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="D48 standard deviation",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(UID, D48sd, size=10, line_color='black', fill_color='black')


    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. D47 vs d47"""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="d47", y_axis_label="D47",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(d47[ETH1], D47raw[ETH1], legend_label="ETH-1", size=ms, line_color='black', fill_color=(255, 51, 230))
    figures[fig_n]['fig'].circle(d47[ETH2], D47raw[ETH2], legend_label="ETH-2", size=ms, line_color='black', fill_color=(255, 102, 153))
    figures[fig_n]['fig'].circle(d47[ETH3], D47raw[ETH3], legend_label="ETH-3", size=ms, line_color='black', fill_color=(255, 153, 77))
    figures[fig_n]['fig'].circle(d47[ETH4], D47raw[ETH4], legend_label="ETH-4", size=ms, line_color='black', fill_color=(255, 204, 0))
    figures[fig_n]['fig'].circle(d47[IAEAC1], D47raw[IAEAC1], legend_label="IAEA-C1", size=ms, line_color='black', fill_color='white')
    figures[fig_n]['fig'].circle(d47[IAEAC2], D47raw[IAEAC2], legend_label="IAEA-C2", size=ms, line_color='black', fill_color='tan')
    figures[fig_n]['fig'].circle(d47[Merck], D47raw[Merck], legend_label="MERCK", size=ms, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(d47[GU1], D47raw[GU1], legend_label="GU1", size=ms, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')

    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. D48 vs d48"""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="d48", y_axis_label="D48",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(d48[ETH1], D48raw[ETH1], legend_label="ETH-1", size=ms, line_color='black', fill_color=(255, 51, 230))
    figures[fig_n]['fig'].circle(d48[ETH2], D48raw[ETH2], legend_label="ETH-2", size=ms, line_color='black', fill_color=(255, 102, 153))
    figures[fig_n]['fig'].circle(d48[ETH3], D48raw[ETH3], legend_label="ETH-3", size=ms, line_color='black', fill_color=(255, 153, 77))
    figures[fig_n]['fig'].circle(d48[ETH4], D48raw[ETH4], legend_label="ETH-4", size=ms, line_color='black', fill_color=(255, 204, 0))
    figures[fig_n]['fig'].circle(d48[IAEAC1], D48raw[IAEAC1], legend_label="IAEA-C1", size=ms, line_color='black', fill_color='white')
    figures[fig_n]['fig'].circle(d48[IAEAC2], D48raw[IAEAC2], legend_label="IAEA-C2", size=ms, line_color='black', fill_color='tan')
    figures[fig_n]['fig'].circle(d48[Merck], D48raw[Merck], legend_label="MERCK", size=ms, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(d48[GU1], D48raw[GU1], legend_label="GU1", size=ms, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')

    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="d13C", y_axis_label="d18O",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(d13C[ETH1], d18O[ETH1], legend_label="ETH-1", size=ms, line_color='black', fill_color=(255, 51, 230))
    figures[fig_n]['fig'].circle(d13C[ETH2], d18O[ETH2], legend_label="ETH-2", size=ms, line_color='black', fill_color=(255, 102, 153))
    figures[fig_n]['fig'].circle(d13C[ETH3], d18O[ETH3], legend_label="ETH-3", size=ms, line_color='black', fill_color=(255, 153, 77))
    figures[fig_n]['fig'].circle(d13C[ETH4], d18O[ETH4], legend_label="ETH-4", size=ms, line_color='black', fill_color=(255, 204, 0))
    figures[fig_n]['fig'].circle(d13C[IAEAC1], d18O[IAEAC1], legend_label="IAEA-C1", size=ms, line_color='black', fill_color='white')
    figures[fig_n]['fig'].circle(d13C[IAEAC2], d18O[IAEAC2], legend_label="IAEA-C2", size=ms, line_color='black', fill_color='tan')
    figures[fig_n]['fig'].circle(d13C[Merck], d18O[Merck], legend_label="Merck", size=ms, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(d13C[GU1], d18O[GU1], legend_label="GU1", size=ms, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')

    fig_n += 1


    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="Residual d13C",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(UID[ETH1], d13C[ETH1]-np.nanmean(d13C[ETH1]), legend_label="ETH-1", size=ms, line_color='black', fill_color=(255, 51, 230))
    figures[fig_n]['fig'].circle(UID[ETH2], d13C[ETH2]-np.nanmean(d13C[ETH2]), legend_label="ETH-2", size=ms, line_color='black', fill_color=(255, 102, 153))
    figures[fig_n]['fig'].circle(UID[ETH3], d13C[ETH3]-np.nanmean(d13C[ETH3]), legend_label="ETH-3", size=ms, line_color='black', fill_color=(255, 153, 77))
    figures[fig_n]['fig'].circle(UID[ETH4], d13C[ETH4]-np.nanmean(d13C[ETH4]), legend_label="ETH-4", size=ms, line_color='black', fill_color=(255, 204, 0))
    figures[fig_n]['fig'].circle(UID[IAEAC1], d13C[IAEAC1]-np.nanmean(d13C[IAEAC1]), legend_label="IAEA-C1", size=ms, line_color='black', fill_color='white')
    figures[fig_n]['fig'].circle(UID[IAEAC2], d13C[IAEAC2]-np.nanmean(d13C[IAEAC2]), legend_label="IAEA-C2", size=ms, line_color='black', fill_color='tan')
    figures[fig_n]['fig'].circle(UID[Merck], d13C[Merck]-np.nanmean(d13C[Merck]), legend_label="Merck", size=ms, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(UID[GU1], d13C[GU1]-np.nanmean(d13C[GU1]), legend_label="GU1", size=ms, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')

    fig_n += 1


    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="Residual d18O",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(UID[ETH1], d18O[ETH1]-np.nanmean(d18O[ETH1]), legend_label="ETH-1", size=ms, line_color='black', fill_color=(255, 51, 230))
    figures[fig_n]['fig'].circle(UID[ETH2], d18O[ETH2]-np.nanmean(d18O[ETH2]), legend_label="ETH-2", size=ms, line_color='black', fill_color=(255, 102, 153))
    figures[fig_n]['fig'].circle(UID[ETH3], d18O[ETH3]-np.nanmean(d18O[ETH3]), legend_label="ETH-3", size=ms, line_color='black', fill_color=(255, 153, 77))
    figures[fig_n]['fig'].circle(UID[ETH4], d18O[ETH4]-np.nanmean(d18O[ETH4]), legend_label="ETH-4", size=ms, line_color='black', fill_color=(255, 204, 0))
    figures[fig_n]['fig'].circle(UID[IAEAC1], d18O[IAEAC1]-np.nanmean(d18O[IAEAC1]), legend_label="IAEA-C1", size=ms, line_color='black', fill_color='white')
    figures[fig_n]['fig'].circle(UID[IAEAC2], d18O[IAEAC2]-np.nanmean(d18O[IAEAC2]), legend_label="IAEA-C2", size=ms, line_color='black', fill_color='tan')
    figures[fig_n]['fig'].circle(UID[Merck], d18O[Merck]-np.nanmean(d18O[Merck]), legend_label="Merck", size=ms, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(UID[GU1], d18O[GU1]-np.nanmean(d18O[GU1]), legend_label="GU1", size=ms, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')

    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="Residual D47raw",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(UID[ETH1], D47raw[ETH1]-np.nanmean(D47raw[ETH1]), legend_label="ETH-1", size=ms, line_color='black', fill_color=(255, 51, 230))
    figures[fig_n]['fig'].circle(UID[ETH2], D47raw[ETH2]-np.nanmean(D47raw[ETH2]), legend_label="ETH-2", size=ms, line_color='black', fill_color=(255, 102, 153))
    figures[fig_n]['fig'].circle(UID[ETH3], D47raw[ETH3]-np.nanmean(D47raw[ETH3]), legend_label="ETH-3", size=ms, line_color='black', fill_color=(255, 153, 77))
    figures[fig_n]['fig'].circle(UID[ETH4], D47raw[ETH4]-np.nanmean(D47raw[ETH4]), legend_label="ETH-4", size=ms, line_color='black', fill_color=(255, 204, 0))
    figures[fig_n]['fig'].circle(UID[IAEAC1], D47raw[IAEAC1]-np.nanmean(D47raw[IAEAC1]), legend_label="IAEA-C1", size=ms, line_color='black', fill_color='white')
    figures[fig_n]['fig'].circle(UID[IAEAC2], D47raw[IAEAC2]-np.nanmean(D47raw[IAEAC2]), legend_label="IAEA-C2", size=ms, line_color='black', fill_color='tan')
    figures[fig_n]['fig'].circle(UID[Merck], D47raw[Merck]-np.nanmean(D47raw[Merck]), legend_label="Merck", size=ms, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(UID[GU1], D47raw[GU1]-np.nanmean(D47raw[GU1]), legend_label="GU1", size=ms, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')

    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="Residual D48raw",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(UID[ETH1], D48raw[ETH1]-np.nanmean(D48raw[ETH1]), legend_label="ETH-1", size=ms, line_color='black', fill_color=(255, 51, 230))
    figures[fig_n]['fig'].circle(UID[ETH2], D48raw[ETH2]-np.nanmean(D48raw[ETH2]), legend_label="ETH-2", size=ms, line_color='black', fill_color=(255, 102, 153))
    figures[fig_n]['fig'].circle(UID[ETH3], D48raw[ETH3]-np.nanmean(D48raw[ETH3]), legend_label="ETH-3", size=ms, line_color='black', fill_color=(255, 153, 77))
    figures[fig_n]['fig'].circle(UID[ETH4], D48raw[ETH4]-np.nanmean(D48raw[ETH4]), legend_label="ETH-4", size=ms, line_color='black', fill_color=(255, 204, 0))
    figures[fig_n]['fig'].circle(UID[IAEAC1], D48raw[IAEAC1]-np.nanmean(D48raw[IAEAC1]), legend_label="IAEA-C1", size=ms, line_color='black', fill_color='white')
    figures[fig_n]['fig'].circle(UID[IAEAC2], D48raw[IAEAC2]-np.nanmean(D48raw[IAEAC2]), legend_label="IAEA-C2", size=ms, line_color='black', fill_color='tan')
    figures[fig_n]['fig'].circle(UID[Merck], D48raw[Merck]-np.nanmean(D48raw[Merck]), legend_label="Merck", size=ms, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(UID[GU1], D48raw[GU1]-np.nanmean(D48raw[GU1]), legend_label="GU1", size=ms, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')

    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="Sample Pressure", y_axis_label="Residual D47raw",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(sample_pressure[ETH1], D47raw[ETH1]-np.nanmean(D47raw[ETH1]), legend_label="ETH-1", size=ms, line_color='black', fill_color=(255, 51, 230))
    figures[fig_n]['fig'].circle(sample_pressure[ETH2], D47raw[ETH2]-np.nanmean(D47raw[ETH2]), legend_label="ETH-2", size=ms, line_color='black', fill_color=(255, 102, 153))
    figures[fig_n]['fig'].circle(sample_pressure[ETH3], D47raw[ETH3]-np.nanmean(D47raw[ETH3]), legend_label="ETH-3", size=ms, line_color='black', fill_color=(255, 153, 77))
    figures[fig_n]['fig'].circle(sample_pressure[ETH4], D47raw[ETH4]-np.nanmean(D47raw[ETH4]), legend_label="ETH-4", size=ms, line_color='black', fill_color=(255, 204, 0))
    figures[fig_n]['fig'].circle(sample_pressure[IAEAC1], D47raw[IAEAC1]-np.nanmean(D47raw[IAEAC1]), legend_label="IAEA-C1", size=ms, line_color='black', fill_color='white')
    figures[fig_n]['fig'].circle(sample_pressure[IAEAC2], D47raw[IAEAC2]-np.nanmean(D47raw[IAEAC2]), legend_label="IAEA-C2", size=ms, line_color='black', fill_color='tan')
    figures[fig_n]['fig'].circle(sample_pressure[Merck], D47raw[Merck]-np.nanmean(D47raw[Merck]), legend_label="Merck", size=ms, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(sample_pressure[GU1], D47raw[GU1]-np.nanmean(D47raw[GU1]), legend_label="GU1", size=ms, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')

    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="Sample Pressure", y_axis_label="Residual D48raw",
                               tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(sample_pressure[ETH1], D48raw[ETH1]-np.nanmean(D48raw[ETH1]), legend_label="ETH-1", size=ms, line_color='black', fill_color=(255, 51, 230))
    figures[fig_n]['fig'].circle(sample_pressure[ETH2], D48raw[ETH2]-np.nanmean(D48raw[ETH2]), legend_label="ETH-2", size=ms, line_color='black', fill_color=(255, 102, 153))
    figures[fig_n]['fig'].circle(sample_pressure[ETH3], D48raw[ETH3]-np.nanmean(D48raw[ETH3]), legend_label="ETH-3", size=ms, line_color='black', fill_color=(255, 153, 77))
    figures[fig_n]['fig'].circle(sample_pressure[ETH4], D48raw[ETH4]-np.nanmean(D48raw[ETH4]), legend_label="ETH-4", size=ms, line_color='black', fill_color=(255, 204, 0))
    figures[fig_n]['fig'].circle(sample_pressure[IAEAC1], D48raw[IAEAC1]-np.nanmean(D48raw[IAEAC1]), legend_label="IAEA-C1", size=ms, line_color='black', fill_color='white')
    figures[fig_n]['fig'].circle(sample_pressure[IAEAC2], D48raw[IAEAC2]-np.nanmean(D48raw[IAEAC2]), legend_label="IAEA-C2", size=ms, line_color='black', fill_color='tan')
    figures[fig_n]['fig'].circle(sample_pressure[Merck], D48raw[Merck]-np.nanmean(D48raw[Merck]), legend_label="Merck", size=ms, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(sample_pressure[GU1], D48raw[GU1]-np.nanmean(D48raw[GU1]), legend_label="GU1", size=ms, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')

# -------------------- create html content for report page --------------------
header = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <!-- py by Andrew Schauer -->
        <meta http-equiv="Content-Type" content="text/html charset=UTF-8" />
        <meta name="viewport" content="width=device-width,initial-scale=1">
        <link rel="stylesheet" type="text/css" href="py_report_style.css">
        <title>{session} Raw Data Report</title>
    </head>"""

body = f"""
    <body>
    <div class="entire_page">
    <h2>Psi Raw Data Report - {session}</h2>
    <div class="created-date">Created - {str(dt.datetime.now())}</div>
    <p>Use these figures to help you identify those samples or standards that should be distrusted before running psi_calibrate.py.</p><br><br><hr>"""


figure_block = [f"""<div class="clear-both">{file_html(figures[i]['fig'], CDN)}{figures[i]['cap']}<br><br><hr></div>""" for i in figures.keys()]


footer = f"""
    </body></div></html>"""

# -------------------- raw data report --------------------

if os.path.exists(archive_path)==False:
    os.mkdir(archive_path)

if os.path.exists(report_path)==False:
    os.mkdir(report_path)
else:
    shutil.move(report_path, os.path.join(archive_path, f"report_{int(dt.datetime.utcnow().timestamp())}"))
    os.mkdir(report_path)


log_summary_page = os.path.join(report_path, 'raw_data_report.html')
with open(log_summary_page, 'w') as report:
    report.write(header)
    report.write(body)

    [report.write(i) for i in figure_block]

    report.write(footer)

webbrowser.open(log_summary_page)
