#!/usr/bin/env python3
"""

Used to assess long term trends across all sessions in Psi. 

"""

__author__ = "Andy Schauer"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-12-19"
__version__ = "1.0"
__copyright__ = "Copyright 2025, Andy Schauer"
__license__ = "Apache 2.0"


from D47crunch import *
from bokeh.plotting import figure, output_file, save, show
from bokeh.embed import file_html
from bokeh.io import export_svgs
from bokeh import palettes
from bokeh.resources import CDN  # , INLINE
from bokeh.models import Span
import csv
import datetime as dt
import json
import isolab_lib
import matplotlib.pyplot as pplt
from natsort import natsorted
import numpy as np
import os
import shutil
import sys
import time
import warnings
import webbrowser


project_path = isolab_lib.get_path("psi", "project")


sessions = {
            "session_230629": {},
            "session_230922": {},
            "session_231016": {},
            "session_231104": {},
            "session_231227": {},
            "session_240122": {},
            "session_240219": {},
            "session_240409": {},
            "session_240618": {},
            "session_240813_dolomite": {},
            "session_241004": {},
            "session_241109": {},
            "session_241203": {},
            "session_241231": {}
            }

starting_UID = [1064, 1809, 2089, 2306, 2889, 3198, 3634, 4156, 4817, 5508, 5934, 6356, 6565, 6739]
carbonate_standards = ['ETH-1', 'ETH-2', 'ETH-3', 'ETH-4', 'IAEA-C1', 'IAEA-C2', 'MERCK', 'GU1']
symbol_colors = [(255, 51, 230), (255, 102, 153), (255, 153, 77), (255, 204, 0), 'white', 'tan', 'green', 'blue']


session_all = {}
for carbonate in carbonate_standards:
    session_all[carbonate] = {}
    session_all[carbonate]['UID'] = []
    session_all[carbonate]['d13C'] = []
    session_all[carbonate]['d18O'] = []
    session_all[carbonate]['D47raw'] = []
    session_all[carbonate]['D48raw'] = []


figures = {}
font_size = '24pt'
ms = 15
fig_n=0

fig_n +=1
figures[fig_n] = {}
figures[fig_n]['name'] = 'NuCarb-ResidualGas__vs__UID'
figures[fig_n]['cap'] = f"""Figure {fig_n}. NuCarb-ResidualGas vs UID for all analyses."""
figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="Residual Gas",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")

fig_n +=1
figures[fig_n] = {}
figures[fig_n]['name'] = 'Background__vs__UID'
figures[fig_n]['cap'] = f"""Figure {fig_n}. Background vs UID for all sessions that contain unknown samples."""
figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="Background (Amps)",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")

fig_n +=1
figures[fig_n] = {}
figures[fig_n]['name'] = 'Balace__vs__UID'
figures[fig_n]['cap'] = f"""Figure {fig_n}. Signal Balance vs UID for all sessions that contain unknown samples."""
figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="Signal Balance",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")

fig_n +=1
figures[fig_n] = {}
figures[fig_n]['name'] = 'D47sd__vs__UID'
figures[fig_n]['cap'] = f"""Figure {fig_n}. D47sd vs UID for all analyses."""
figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="D47sd",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")

fig_n +=1
figures[fig_n] = {}
figures[fig_n]['name'] = 'D48sd__vs__UID'
figures[fig_n]['cap'] = f"""Figure {fig_n}. D48sd vs UID for all analyses."""
figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="D48sd",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")

fig_n +=1
figures[fig_n] = {}
figures[fig_n]['name'] = 'd13C_residual__vs__UID'
figures[fig_n]['cap'] = f"""Figure {fig_n}. d13C residual vs UID for all sessions that contain unknown samples."""
figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="d13C residual",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")

fig_n +=1
figures[fig_n] = {}
figures[fig_n]['name'] = 'd18O_residual__vs__UID'
figures[fig_n]['cap'] = f"""Figure {fig_n}. d18O residual vs UID for all sessions that contain unknown samples."""
figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="d18O residual",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")

fig_n +=1
figures[fig_n] = {}
figures[fig_n]['name'] = 'D47raw_residual__vs__UID'
figures[fig_n]['cap'] = f"""Figure {fig_n}. D47raw residual vs UID for all sessions that contain unknown samples."""
figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="D47raw residual",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")

fig_n +=1
figures[fig_n] = {}
figures[fig_n]['name'] = 'D48raw_residual__vs__UID'
figures[fig_n]['cap'] = f"""Figure {fig_n}. D48raw residual vs UID for all sessions that contain unknown samples."""
figures[fig_n]['fig'] = figure(width=1400, height=700, x_axis_label="UID", y_axis_label="D48raw residual",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")


# loop through all sessions and load calibrated data
for session in sessions:
    fig_n = 0
    print(session)
    headers, data = isolab_lib.read_file(os.path.join(project_path, session, f'psi_{session}_samplelog.csv'), ',')

    trust0_indices = [i for i, e in enumerate(data['trust']) if int(e) == 0]
    trust1_indices = [i for i, e in enumerate(data['trust']) if int(e) == 1]
    for header in headers:
        data[header] = [data[header][index] for index in trust1_indices]

    sessions[session] = {'headers':headers, 'data':data}

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
        sessions[session]['data'][i] = [x if x != "NA" else np.nan for x in sessions[session]['data'][i]]
        sessions[session]['data'][i] = np.asarray(sessions[session]['data'][i], dtype=float)
    else:
        sessions[session]['data'][i] = np.asarray(sessions[session]['data'][i])

    UID = sessions[session]['data']['UID']
    d13C = sessions[session]['data']['d13C']
    d18O = sessions[session]['data']['d18O']
    D47raw = sessions[session]['data']['D47raw']
    D48raw = sessions[session]['data']['D48raw']

    for carbonate in carbonate_standards:
        sessions[session][carbonate] = {}
        sessions[session][carbonate]['index'] = [i for i, e in enumerate(sessions[session]['data']['Sample']) if e==carbonate]
        sessions[session][carbonate]['d13C'] = d13C[sessions[session][carbonate]['index']]
        sessions[session][carbonate]['d18O'] = d18O[sessions[session][carbonate]['index']]
        sessions[session][carbonate]['D47raw'] = D47raw[sessions[session][carbonate]['index']]
        sessions[session][carbonate]['D48raw'] = D48raw[sessions[session][carbonate]['index']]
        session_all[carbonate]['UID'].append(UID[sessions[session][carbonate]['index']])
        session_all[carbonate]['d13C'].append(d13C[sessions[session][carbonate]['index']])
        session_all[carbonate]['d18O'].append(d18O[sessions[session][carbonate]['index']])
        session_all[carbonate]['D47raw'].append(D47raw[sessions[session][carbonate]['index']])
        session_all[carbonate]['D48raw'].append(D48raw[sessions[session][carbonate]['index']])

    
    fig_n += 1
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['residual_gas'], size=10, line_color='black', fill_color='black')    
    
    fig_n += 1
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['bkgnd_44'], legend_label="m/z 44", size=10, line_color='black', fill_color='black')
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['bkgnd_45'], legend_label="m/z 45", size=10, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['bkgnd_46'], legend_label="m/z 46", size=10, line_color='black', fill_color='pink')
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['bkgnd_47'], legend_label="m/z 47", size=10, line_color='black', fill_color='green')
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['bkgnd_48'], legend_label="m/z 48", size=10, line_color='black', fill_color='yellow')

    fig_n += 1
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['nu_balance'], legend_label="Nu calculated", size=10, line_color='black', fill_color='black')
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['min_balance'], legend_label="IsoLab min", size=10, line_color='black', fill_color='blue')
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['mean_balance'], legend_label="IsoLab mean", size=10, line_color='black', fill_color='pink')
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['max_balance'], legend_label="IsoLab max", size=10, line_color='black', fill_color='green')

    fig_n += 1
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['D47sd'], size=10, line_color='black', fill_color='black')    
    
    fig_n += 1
    figures[fig_n]['fig'].circle(sessions[session]['data']['UID'], sessions[session]['data']['D48sd'], size=10, line_color='black', fill_color='black')    



for i,carbonate in enumerate(carbonate_standards):
    session_all[carbonate]['UID_combined'] = np.concatenate(session_all[carbonate]['UID'])
    session_all[carbonate]['d13C_combined'] = np.concatenate(session_all[carbonate]['d13C'])
    session_all[carbonate]['d13C_mean'] = np.mean(session_all[carbonate]['d13C_combined'])
    session_all[carbonate]['d18O_combined'] = np.concatenate(session_all[carbonate]['d18O'])
    session_all[carbonate]['d18O_mean'] = np.mean(session_all[carbonate]['d18O_combined'])
    session_all[carbonate]['D47raw_combined'] = np.concatenate(session_all[carbonate]['D47raw'])
    session_all[carbonate]['D47raw_mean'] = np.mean(session_all[carbonate]['D47raw_combined'])
    session_all[carbonate]['D48raw_combined'] = np.concatenate(session_all[carbonate]['D48raw'])
    session_all[carbonate]['D48raw_mean'] = np.mean(session_all[carbonate]['D48raw_combined'])

    fig_n = 6
    figures[fig_n]['fig'].circle(session_all[carbonate]['UID_combined'], session_all[carbonate]['d13C_combined'] - session_all[carbonate]['d13C_mean'], legend_label=carbonate, size=10, line_color='black', fill_color=symbol_colors[i])

    fig_n += 1
    figures[fig_n]['fig'].circle(session_all[carbonate]['UID_combined'], session_all[carbonate]['d18O_combined'] - session_all[carbonate]['d18O_mean'], legend_label=carbonate, size=10, line_color='black', fill_color=symbol_colors[i])

    fig_n += 1
    figures[fig_n]['fig'].circle(session_all[carbonate]['UID_combined'], session_all[carbonate]['D47raw_combined'] - session_all[carbonate]['D47raw_mean'], legend_label=carbonate, size=10, line_color='black', fill_color=symbol_colors[i])

    fig_n += 1
    figures[fig_n]['fig'].circle(session_all[carbonate]['UID_combined'], session_all[carbonate]['D48raw_combined'] - session_all[carbonate]['D48raw_mean'], legend_label=carbonate, size=10, line_color='black', fill_color=symbol_colors[i])



fig_n = 1
vline = Span(location=5793, dimension='height', line_color='red', line_width=1)
figures[fig_n]['fig'].add_layout(vline)
for i in starting_UID:
    vline = Span(location=i, dimension='height', line_color='green', line_width=1)
    figures[fig_n]['fig'].add_layout(vline)

fig_n += 1
vline = Span(location=5793, dimension='height', line_color='red', line_width=1)
figures[fig_n]['fig'].add_layout(vline)
figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')
for i in starting_UID:
    vline = Span(location=i, dimension='height', line_color='green', line_width=1)
    figures[fig_n]['fig'].add_layout(vline)

fig_n += 1
vline = Span(location=5793, dimension='height', line_color='red', line_width=1)
figures[fig_n]['fig'].add_layout(vline)
figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')
for i in starting_UID:
    vline = Span(location=i, dimension='height', line_color='green', line_width=1)
    figures[fig_n]['fig'].add_layout(vline)

fig_n += 1
vline = Span(location=5793, dimension='height', line_color='red', line_width=1)
figures[fig_n]['fig'].add_layout(vline)
for i in starting_UID:
    vline = Span(location=i, dimension='height', line_color='green', line_width=1)
    figures[fig_n]['fig'].add_layout(vline)

fig_n += 1
vline = Span(location=5793, dimension='height', line_color='red', line_width=1)
figures[fig_n]['fig'].add_layout(vline)
for i in starting_UID:
    vline = Span(location=i, dimension='height', line_color='green', line_width=1)
    figures[fig_n]['fig'].add_layout(vline)

fig_n += 1
vline = Span(location=5793, dimension='height', line_color='red', line_width=1)
figures[fig_n]['fig'].add_layout(vline)
figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')
for i in starting_UID:
    vline = Span(location=i, dimension='height', line_color='green', line_width=1)
    figures[fig_n]['fig'].add_layout(vline)

fig_n += 1
vline = Span(location=5793, dimension='height', line_color='red', line_width=1)
figures[fig_n]['fig'].add_layout(vline)
figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')
for i in starting_UID:
    vline = Span(location=i, dimension='height', line_color='green', line_width=1)
    figures[fig_n]['fig'].add_layout(vline)

fig_n += 1
vline = Span(location=5793, dimension='height', line_color='red', line_width=1)
figures[fig_n]['fig'].add_layout(vline)
figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')
for i in starting_UID:
    vline = Span(location=i, dimension='height', line_color='green', line_width=1)
    figures[fig_n]['fig'].add_layout(vline)

fig_n += 1
vline = Span(location=5793, dimension='height', line_color='red', line_width=1)
figures[fig_n]['fig'].add_layout(vline)
figures[fig_n]['fig'].add_layout(figures[fig_n]['fig'].legend[0], 'right')
for i in starting_UID:
    vline = Span(location=i, dimension='height', line_color='green', line_width=1)
    figures[fig_n]['fig'].add_layout(vline)


# -------------------- create html content for report page --------------------
header = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <!-- py by Andrew Schauer -->
        <meta http-equiv="Content-Type" content="text/html charset=UTF-8" />
        <meta name="viewport" content="width=device-width,initial-scale=1">
        <link rel="stylesheet" type="text/css" href="py_report_style.css">
        <title>Psi - all sessions - figures</title>
    </head>"""

body = f"""
    <body>
    <div class="entire_page">
    <h2></h2>
    <div class="created-date">Created - {str(dt.datetime.now())}</div>
    <p></p>"""


figure_block = [f"""<div class="clear-both">{file_html(figures[i]['fig'], CDN)}{figures[i]['cap']}<br><br><hr></div>""" for i in figures.keys()]


footer = f"""
    </body></div></html>"""


html_figures_page = os.path.join(project_path, 'sessions_figures.html')
with open(html_figures_page, 'w') as html_page:
    html_page.write(header)
    html_page.write(body)

    [html_page.write(i) for i in figure_block]

    html_page.write(footer)

webbrowser.open(html_figures_page)
