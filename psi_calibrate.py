"""
This script started out as dressing for M. Daeron's D47crunch carbonate clumped isotope python library. As of version
1.0, I am trying to have it write a proper report in addition to making full use of D47crunch.

Version 0.1 - 2023.05.26 => created
Version 0.5 - 2023.06.02 => now we can set a samples flag to 0 (ignore or exclude) in the raw data file and pass
    only flag 1 data to D47crunch.
Version 0.6 - 2024.02.06 => this script is now called from psi_calibrate.Rmd, an Rmarkdown that oversees the calibration
    of psi data followed by creation of a report.
Version 0.7 - 2024.04.30 => created way to choose a sample log file to process after dividing sessions into chop2 and chop3
Version 1.0 - 2024.05.31 => started creating proper report
Version 1.1 - 2024.07.10 => added GU1 to D48 calibration; added repeated removal of standard as an estimate of accuracy
Version 1.2 - 2024.08.09 => changed flag to trust
Version 1.3 - 2024.08.14 => updating to get session read in as variable from data file and used accordingly
Version 1.4 - 2024.08.16 => can filter based on nchops at the command line; included Emma's addition of equilibrated gases
Version 1.5 - 2024.09.16 => added Kim et al 2015 18O acid fractionation equation
Version 2.0 - 2024.12.06 => trying to finish this off, or at least make it a final looking report, rather than the snark it is currently populated with
"""

__author__ = "Andy Schauer"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-12-06"
__version__ = "2.0"
__copyright__ = "Copyright 2024, Andy Schauer"
__license__ = "Apache 2.0"
__acknowledgements__ = "K. Huntington, E. Heitmann, C. Saenger, S. Mat, V. Ravi, M. Leite. Special thanks to M. Daeron for D47crunch!"



# -------------------- IMPORTS --------------------
import argparse
import csv
from D47crunch import *
import isolab_lib
import json
import matplotlib.pyplot as pplt
from natsort import natsorted
import os
import shutil
import sys
import time as t
import webbrowser



# ---------- GET SESSION ----------
project_path = isolab_lib.get_path("psi", "project")
session_list = natsorted(os.listdir(project_path))
exclude_from_session_list = ['archive', 'compilations']
session_list = [session for session in session_list if session not in exclude_from_session_list]
print('\nChoose a session from the list below:')
[print(f'    {i}') for i in session_list]
identified_session = 0
while identified_session == 0:
    session_search = input('\nEnter the session you wish to process: ')
    isdir = [session_search[0: len(session_search)] in x for x in session_list]
    if len(np.where(isdir)[0]) == 1:
        identified_session = 1
        session = session_list[np.where(isdir)[0][0]]
        print(f'    Processing session {session}...')
    else:
        print('\n** More than one session found. **\n')

session_path = f'{project_path}{session}/'
report_path = f"{session_path}report/"


session_file_list = natsorted([file for file in os.listdir(session_path) if "samplelog.csv" in file])
if len(session_file_list) == 1:
    identified_session_file = 1
    session_samplelog = session_file_list[0]
    print(f'    Processing session file {session_samplelog}...')
else:
    print('\nChoose a session log file from the list below:')
    [print(f'    {i}') for i in session_file_list]
    identified_session_file = 0
    while identified_session_file == 0:
        session_file_search = input('\n\nEnter the session log file you wish to process: ')
        isdir = [session_file_search[0: len(session_file_search)] in x for x in session_file_list]
        if len(np.where(isdir)[0]) == 1:
            identified_session_file = 1
            session_samplelog = session_file_list[np.where(isdir)[0][0]]
            print(f'    Processing session file {session_samplelog}...')
        else:
            print('\n** More than one session file found. **\n')

D47crunch_input_file = f"psi_{session}_D47crunch_input.csv"
D47crunch_output_analyses_file = f"{session}_D47crunch_output_all_analyses.csv"
D47crunch_output_samples_file = f"{session}_D47crunch_output_summarized.csv"
D47crunch_remove_one_error =  f"{session}_D47crunch_remove_one_error_summarized.csv"
psi_calibrated_session_file = f"psi_{session}_calibrated.csv"



# ------------------ CREATE REPORT DIRECTORY -------------------------------
import datetime as dt
archive_path = f"{session_path}archive/"
if os.path.exists(archive_path)==False:
    os.mkdir(archive_path)

if os.path.exists(report_path):
    shutil.move(report_path, os.path.join(archive_path, f"report_{int(dt.datetime.utcnow().timestamp())}"))

os.mkdir(report_path)
os.mkdir(f"{report_path}data/")
os.mkdir(f"{report_path}figures/")
os.mkdir(f"{report_path}scripts/")

shutil.copy2(os.path.join(session_path, session_samplelog), os.path.join(report_path, 'data/', session_samplelog))




# -------------------- SCRIPTS --------------------
script_path = isolab_lib.get_path("psi", "python")
scripts = {'psi.R': '', 'isolab_lib.py': '', 'psi_raw_data_report.py': '', 'psi_calibrate.py': ''}
scripts = {key: (t.strftime('%Y-%m-%d %H:%M:%S', t.localtime(os.path.getmtime(f'{script_path}{key}')))) for key, value in scripts.items()}

[shutil.copy2(os.path.join(script_path, script), os.path.join(report_path, f"scripts/{script}_REPORT_COPY")) for script in scripts]
shutil.copy2(os.path.join(script_path, 'py_report_style.css'), os.path.join(report_path, 'py_report_style.css'))



# -------------------- PRE-LOAD RAW DATA FILE AND REMOVE TRUST=0 --------------------
with open(f"{session_path}{session_samplelog}", 'r') as datafile:
    reader = csv.reader(datafile)
    headers = next(reader)

    n_everything = 0
    rows = []
    for row in reader:
        n_everything += 1
        if int(row[4]) == 1:
            rows.append(row)

headers.append('Session')
for row in rows:
    row.append(session)

with open(f"{report_path}data/{D47crunch_input_file}", 'w', newline='') as datafile:
    writer = csv.writer(datafile)
    writer.writerow(headers)
    writer.writerows(rows)



# -------------------- REFERENCE MATERIALS --------------------

# CARBONATES
#    - read in json file, make dictionaries for each carbonate, and create an index list within each dictionary
with open(isolab_lib.get_path("psi", "standards"), 'r') as refmat_file:
    refmat = json.load(refmat_file)

refmat_keys = refmat['carbonates'].keys()
for i in refmat_keys:
    globals()[i] = refmat['carbonates'][i]
    globals()[i]['index'] = np.empty(0, dtype="int16")


# GASES
# Theoretical Values
t4 = 4 + 273.15
t25 = 23 + 273.15
t60 = 57 + 273.15
t1000 = 1000 + 273.15
        
# “True D47” calculated using Eq. 28 from Wang et al. (2004) and IUPAC parameters ***Preferred equation to replace Dennis et al. (2011) Eq. A2***
t4_D47rf = 0.00050479*(1000/t4)**7 - 0.00885734*(1000/t4)**6 + 0.06385048*(1000/t4)**5 - 0.23891768*(1000/t4)**4 + 0.46854990*(1000/t4)**3 - 0.34158204*(1000/t4)**2 + 0.12940422*(1000/t4) - 0.01752753;
t25_D47rf = 0.00050479*(1000/t25)**7 - 0.00885734*(1000/t25)**6 + 0.06385048*(1000/t25)**5 - 0.23891768*(1000/t25)**4 + 0.46854990*(1000/t25)**3 - 0.34158204*(1000/t25)**2 + 0.12940422*(1000/t25) - 0.01752753;
t60_D47rf = 0.00050479*(1000/t60)**7 - 0.00885734*(1000/t60)**6 + 0.06385048*(1000/t60)**5 - 0.23891768*(1000/t60)**4 + 0.46854990*(1000/t60)**3 - 0.34158204*(1000/t60)**2 + 0.12940422*(1000/t60) - 0.01752753;
t1000_D47rf = 0.00050479*(1000/t1000)**7 - 0.00885734*(1000/t1000)**6 + 0.06385048*(1000/t1000)**5 - 0.23891768*(1000/t1000)**4 + 0.46854990*(1000/t1000)**3 - 0.34158204*(1000/t1000)**2 + 0.12940422*(1000/t1000) - 0.01752753;



# USE THESE CARBONATES AND GASES TO CALIBRATE D47 and D48
d13C_calibration_standards = ['Merck', 'ETH4', 'ETH1']
d18O_calibration_standards = ['Merck', 'ETH4', 'ETH1']
D47_calibration_standards = ['ETH1', 'ETH2', 'ETH4', 'Merck', 'IAEAC1', 'IAEAC2']
D48_calibration_standards = ['ETH1', 'ETH2', 'ETH4', 'GU1']



# -------------------- D47crunch --------------------
psi47 = D47data(session=session)

rxnT = 70
psi47.ALPHA_18O_ACID_REACTION = np.exp((3.48 * (1000/(rxnT+273.15)) - 1.47)/1000) # Kim et al 2015

psi47.Nominal_d13C_VPDB = {eval(i)['names'][0]: eval(i)['d13C'] for i in d13C_calibration_standards}
psi47.Nominal_d18O_VPDB = {eval(i)['names'][0]: eval(i)['d18O_vpdb'] for i in d18O_calibration_standards}
psi47.Nominal_D47 = {eval(i)['names'][0]: eval(i)['D47'] for i in D47_calibration_standards}

psi47.read(f"{report_path}data/{D47crunch_input_file}")
psi47.wg()
psi47.crunch()
psi47.standardize()
psi47.summary(verbose=True, save_to_file=True)


psi48 = D48data(psi47, session=session)
psi48.Nominal_D48 = {eval(i)['names'][0]: eval(i)['D48'] for i in D48_calibration_standards}
psi48.Nominal_d13C_VPDB = {eval(i)['names'][0]: eval(i)['d13C'] for i in d13C_calibration_standards}
psi48.Nominal_d18O_VPDB = {eval(i)['names'][0]: eval(i)['d18O_vpdb'] for i in d18O_calibration_standards}
psi48.crunch()
psi48.standardize()
psi48.summary(verbose=True)


table_of_samples(data47 = psi47,
                  data48 = psi48,
                  dir = f"{report_path}data/",
                  filename = D47crunch_output_samples_file,
                  save_to_file = True,
                  print_out = True,
                  output = None)


table_of_analyses(data47 = psi47,
                  data48 = psi48,
                  dir = f"{report_path}data/",
                  filename = D47crunch_output_analyses_file,
                  save_to_file = True,
                  print_out = False,
                  output = None)




# ------------------ remove-one error analysis -------------------------------
# psi47_remove_one_error_list = ['ETH-1', 'ETH-2', 'ETH-3', 'ETH-4', 'IAEA-C1', 'MERCK', 'IAEA-C2']
# psi47_removed_D47 = {}

# for removed in psi47_remove_one_error_list:

    # psi47_remove_one_error = D47data(session=session)
    # [j for j in psi47_remove_one_error_list if j != removed]

    # psi47_remove_one_error.Nominal_D47 = {key: value for key, value in psi47.Nominal_D47.items() if key != removed}

    # psi47_remove_one_error.read(f"{session_path}{D47crunch_input_file}")
    # psi47_remove_one_error.wg()
    # psi47_remove_one_error.crunch()
    # psi47_remove_one_error.standardize()

    # table_of_samples(data47 = psi47_remove_one_error,
                      # dir = f"{session_path}",
                      # filename = D47crunch_remove_one_error,
                      # save_to_file = True,
                      # print_out = True,
                      # output = None)

    # herror, derror = isolab_lib.read_file(os.path.join(session_path, D47crunch_remove_one_error), ',')

    # removed_index = [j for j, e in enumerate(derror['Sample']) if str(e).lower() == removed.lower()]
    # psi47_removed_D47[removed] = {'mean': derror['D47'][removed_index[0]], 'se': derror['SE'][removed_index[0]], 'accuracy': float(derror['D47'][removed_index[0]])-psi47.Nominal_D47[removed]}

    # y = {outer_key: inner_dict['accuracy'] for outer_key, inner_dict in psi47_removed_D47.items()}
    # y = [key for value, key in y.items()]
    # x = [value for key, value in psi47.Nominal_D47.items()]



# -------------------- combine Psi meta data with D47crunch output data --------------------
h1, d1 = isolab_lib.read_file(os.path.join(report_path, 'data/', D47crunch_input_file), ',')
h2, d2 = isolab_lib.read_file(os.path.join(report_path, 'data/', D47crunch_output_analyses_file), ',')

d3 = {**d1, **d2}
h3=list(d3.keys())
h3[-2] = 'D47_iCDES'
# h3[-1] = 'D48_iCDES'

with open(f"{session_path}{session_samplelog}", 'r') as datafile:
    reader = csv.reader(datafile)
    headers0 = next(reader)
    rows0 = [row for row in reader if int(row[4]) == 0]


with open(os.path.join(report_path, 'data/', psi_calibrated_session_file), 'w', newline='') as cal_sesh_file:
        datawriter = csv.writer(cal_sesh_file)
        datawriter.writerow(d3.keys())
        datawriter.writerows(zip(*d3.values()))
        datawriter.writerow("")
        datawriter.writerows(rows0)

psi47.plot_sessions(dir=f"{report_path}figures/")
psi48.plot_sessions(dir=f"{report_path}figures/")



# ------------------ read in final calibrated file for interactive plotting etc -------------------------------

headers, data = isolab_lib.read_file(os.path.join(report_path, 'data/', psi_calibrated_session_file), ',')

original_data = data.copy()
trust0_indices = [i for i, e in enumerate(data['trust']) if int(e) == 0]
trust1_indices = [i for i, e in enumerate(data['trust']) if int(e) == 1]
for header in headers[:-1]:
    data[header] = [data[header][index] for index in trust1_indices]

numlist = ['UID', 'trust', 'd13C', 'd13Csd', 'd13Cn', 'd13Cslope', 'd13Cr2', 'd13Cpval', 'd18O', 'd18Osd', 'd18On', 'd18Oslope', 'd18Or2', 'd18Opval',
           'd45', 'd45sd', 'd45n', 'd45slope', 'd45r2', 'd45pval', 'd46', 'd46sd', 'd46n', 'd46slope', 'd46r2', 'd46pval',
           'd47', 'd47sd', 'd47n', 'd47slope', 'd47r2', 'd47pval', 'D47raw', 'D47sd', 'D47n', 'D47slope', 'D47r2', 'D47pval',
           'd48', 'd48sd', 'd48n', 'd48slope', 'd48r2', 'd48pval', 'D48raw', 'D48sd', 'D48n', 'D48slope', 'D48r2', 'D48pval',
           'd49', 'd49sd', 'd49n', 'd49slope', 'd49r2', 'd49pval', 'D49raw', 'D49sd', 'D49n', 'D49slope', 'D49r2', 'D49pval',
           'R45sam_mean', 'R45sam_sd', 'R45wg_mean', 'R45wg_sd', 'R45sam_wg', 'R46sam_mean', 'R46sam_sd', 'R46wg_mean', 'R46wg_sd', 'R46sam_wg',
           'R47sam_mean', 'R47sam_sd', 'R47wg_mean', 'R47wg_sd', 'R47sam_wg', 'R48sam_mean', 'R48sam_sd', 'R48wg_mean', 'R48wg_sd', 'R48sam_wg',
           'R49sam_mean', 'R49sam_sd', 'R49wg_mean', 'R49wg_sd', 'R49sam_wg', 'unix_time', 'mass', 'sample_pressure', 'yield_mbar_ug', 'residual_gas',
           'nchops', 'min_ref_beam', 'max_ref_beam', 'total_blocks', 'total_cycles', 'initial_sample_beam', 'pre_balance_sample_beam',
           'nu_balance', 'sam_max44nA', 'sam_mean44nA', 'sam_min44nA', 'sam_max47nA', 'sam_mean47nA', 'sam_min47nA', 'min_balance', 'mean_balance', 'max_balance',
           'sam_decay', 'wg_decay', 'd13C_VPDB', 'd18O_VSMOW', 'D47', 'D48']

strlist = set(headers) - set(numlist)

for i in numlist:
    data[i] = [x if x != "NA" else np.nan for x in data[i]]
    globals()[i] = np.asarray(data[i], dtype=float)

for i in strlist:
    globals()[i] = np.asarray(data[i])

ETH1['index'] = [i for i, e in enumerate(Sample) if e=='ETH-1']
ETH2['index'] = [i for i, e in enumerate(Sample) if e=='ETH-2']
ETH3['index'] = [i for i, e in enumerate(Sample) if e=='ETH-3']
ETH4['index'] = [i for i, e in enumerate(Sample) if e=='ETH-4']
IAEAC1['index'] = [i for i, e in enumerate(Sample) if e=='IAEA-C1']
IAEAC2['index'] = [i for i, e in enumerate(Sample) if e=='IAEA-C2']
Merck['index'] = [i for i, e in enumerate(Sample) if e=='MERCK']
GU1['index'] = [i for i, e in enumerate(Sample) if e=='GU1']

n_trusted = len(set(UID))
n_standards = len(ETH1['index']) + len(ETH2['index']) + len(ETH3['index']) + len(ETH4['index']) + len(IAEAC1['index']) + len(IAEAC2['index']) + len(Merck['index']) + len(GU1['index'])
n_samples = n_trusted - n_standards

report_page = os.path.join(report_path, f'{session}_calibration_summary.html')


# ---------- REPORT BITS ---------- 
print('Making html page...')
header = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <!-- py by Andrew Schauer -->
        <meta http-equiv="Content-Type" content="text/html charset=UTF-8" />
        <meta name="viewport" content="width=device-width,initial-scale=1">
        <link rel="stylesheet" type="text/css" href="py_report_style.css">
        <title>Psi Calibration Report</title>
    </head>"""


calculation_notes_block = """<li>An R script called psi.R was used to take high frequency faraday-collector current measurements and convert them to isotope ratios and raw measured delta values.</li>
                             <li>Analyses were screened for quality based on a variety of measures. Distrusted analyses are given a "0" value in the "trust" column.</li>
                             <li>Psi_calibrate.py is used as an entry into D47crunch.py, Mathieu Daeron's clumped isotope python library.</li>"""


calibration_standards = list(set(d13C_calibration_standards + d18O_calibration_standards + D47_calibration_standards + D48_calibration_standards))

for i in calibration_standards:
    eval(i)['purpose'] = ''
    if i in d13C_calibration_standards:
        eval(i)['purpose'] += 'd13C-'

for i in calibration_standards:
    if i in d18O_calibration_standards:
        eval(i)['purpose'] += 'd18O-'

for i in calibration_standards:
    if i in D47_calibration_standards:
        eval(i)['purpose'] += 'D47-'

for i in calibration_standards:
    if i in D48_calibration_standards:
        eval(i)['purpose'] += 'D48-'


refmat_block = str([f"""<tr>
                             <td>{eval(i)['names'][0]}</td>
                             <td>{eval(i)['material']}</td>
                             <td>{eval(i)['d13C']}</td>
                             <td>{eval(i)['d18O_vpdb']}</td>
                             <td>{eval(i)['D47']}</td>
                             <td>{eval(i)['D48']}</td>
                             <td>{eval(i)['purpose']}</td>
                         </tr>""" for i in calibration_standards]).replace("[", "").replace("'", "").replace("]", "").replace(", ", "").replace("\\n", "")


D47_summary_table = "<table>\n<tr><th></th><th></th></tr>\n"
with open("output/D47_summary.csv", "r") as f:
    reader = csv.reader(f)
    for i, row in enumerate(reader):
        D47_summary_table += "<tr>" + "".join(f"<td>{col}</td>" for col in row) + "</tr>\n"
D47_summary_table += "</table>"


D48_summary_table = "<table>\n<tr><th></th><th></th></tr>\n"
with open("output/D48_summary.csv", "r") as f:
    reader = csv.reader(f)
    for i, row in enumerate(reader):
        D48_summary_table += "<tr>" + "".join(f"<td>{col}</td>" for col in row) + "</tr>\n"
D48_summary_table += "</table>"


replicate_summary_table = "<table>\n"
with open(f"{report_path}data/{D47crunch_output_samples_file}", "r") as f:
    reader = csv.reader(f)
    for i, row in enumerate(reader):
        if i == 0:
            replicate_summary_table += "<tr>" + "".join(f"<th>{col}</th>" for col in row) + "</tr>\n"
        else:
            replicate_summary_table += "<tr>" + "".join(f"<td>{col}</td>" for col in row) + "</tr>\n"
replicate_summary_table += "</table>"

body = f"""
    <body>
    <div class="entire_page">
    <h2>Psi Calibration Report - {session}</h2>
    <div class="created-date">Created - {str(dt.datetime.now())}</div>
    <h2>Introduction</h2>
    <div class="text-indent">
        <p>This report is meant to be a stand-alone collection of methods,
        data, scripts, and notes related to calibration of carbonate clumped
        isotope data generated with Psi, a Nu Perspective / NuCarb isotope
        ratio mass spectrometer / automated carbonate digestion system.
        You can read more about our implementation of this method on our website
        <a href="https://isolab.ess.washington.edu/SOPs/psi.php">
        https://isolab.ess.washington.edu/SOPs/psi.php</a>.</p>

        <p>The data and python scripts used to generate this page are linked
        and described in the <a href="#refs">References</a> section below. If you wish
        to save this report, <a href="report.zip">save the zip file</a> or copy and paste
        or download the entire 'report' directory to a place of your choosing and all
        html, images, data files, and python scripts will be saved. <strong>
        &iexcl; <a href="report.zip">Save a copy if you are finished analyzing your samples</a> !</strong></p>
    </div>

    <h2>My data</h2>
    <div class="text-indent">
        <p>This technical stuff is fine and all but where are <strong><a href="data/{psi_calibrated_session_file}">my data</a></strong>?
        This calibrated data file contains sample IDs, dates of analyses, unique analysis numbers, total mass weighed for analysis,
        blah blah and blah. Two sections of data are separated by an empty row. The first section of data are trusted; the second section
        of data are not trusted. Under the "trust" heading, "1" indicates good, trusted data while "0" indicates poor
        quality data that should probably be distrusted. Up to you if you want to use it. Untrusted data are given the reason for
        distrust. If you want to negate the comments and trust and distrusted analysis, find the raw data samplelog in the main session folder
        and change the trust value to "1"; then rerun psi_calibrate.py. <strong>&iexcl; <a href="report.zip">If you are done analyzing samples,
        save a copy of the entire report directory elsewhere, not just a copy of your data file</a> !</strong></p>
    </div>

    <h2>Data operations</h2>
    <div class="text-indent">
        <p>A suite of mathmatical operations were completed on these data prior to claiming they are final. Here are notes associated with
        these calculations:
        <ul>
            {calculation_notes_block}
        </ul>
        </p>
    </div>

    <h2>Run inventory</h2>
    <div class="text-indent">
        <table>
            <tr><th></th><th></th></tr>
            <tr><td>Total number of trusted analyses in session</td><td>{n_trusted}</td></tr>
            <tr><td>Total number of trusted standards in session</td><td>{n_standards}</td></tr>
            <tr><td>Total number of trusted samples in session</td><td>{n_samples}</td></tr>
            <tr><td><br></td></tr>
            <tr><td>Number of <a href="#distrusted">distrusted analyses in session</a></td><td>{n_everything - n_trusted}</td></tr>
        </table>
    </div>

    <h2>Reference materials</h2>
    <div class="text-indent">
        <p>Table showing all of our reference carbonates, their accepted values, and their calibration purpose. Accepted values for &Delta;<sub>47</sub> come from Bernasconi et al. 2021. &Delta;<sub>48</sub> accepted values for GU1, ETH-1, ETH-2, and ETH-4 are from Fiebig et al. 2024 (recommended by Jens Fiebig over email exchange with Andy on 2024-12-06). Given the nitrate contamination in ETH-3, I suggest not using it as a calibration standard until we reach consensus on how to proceed.</p>
        <table>
            <tr><th>Reference<br>name</th><th>Reference<br>material</th><th>&delta;<sup>13</sup>C<br>accepted<br>(permil)</th><th>&delta;<sup>18</sup>O<br>accepted<br>(permil)</th><th>&Delta;<sub>47</sub></th><th>&Delta;<sub>48</sub></th><th>Purpose</th></tr>
            {refmat_block}
            <tr><td> </td><td></td><td> </td><td> </td><td> </td><td> </td><td> </td></tr>
            <tr><td>ETH-3</td><td></td><td>{ETH3['d13C']}</td><td>{ETH3['d18O_vpdb']}</td><td>{ETH3['D47']}</td><td>{ETH3['D48']}</td><td>** Excluded until we deal with nitrate issue</td></tr>
        </table>
    </div>

    <h2>Data quality</h2>
    <div class="text-indent">
        <div class="figThumb"><strong>&Delta;<sub>47</sub> Repeatability</strong>
        {D47_summary_table}</div>

        <div class="figThumb"><strong>&Delta;<sub>48</sub> Repeatability</strong>
        {D48_summary_table}</div>

        <div class="clear-both"></div>

        <h3>Replicate Summary Table</h3>
        <p>{replicate_summary_table}</p>
    </div>



    <h2>Figures</h2>"""

figure_block = f"""<div class="clear-both"><img src="figures/D47_plot_{session}.png"><hr></div>
                   <div class="clear-both"><img src="figures/D48_plot_{session}.png"><hr></div>
                   """

# python_scripts_block = str([f'<li><a href="python/{key}_REPORT_COPY">{key}</a> - {value}</li>' for key, value in python_scripts.items()]).replace("[", "").replace("'", "").replace("]", "").replace(", ", "")
python_scripts_block = 'python scripts here'

footer = f"""
    <h2 id="refs">References</h2>
    <div class="references">
        <ul>
            <li>Bernasconi et al. 2021. InterCarb: A Community Effort to Improve Interlaboratory Standardization of the Carbonate Clumped Isotope Thermometer Using Carbonate Standards. <em>Geochemistry, Geophysics, Geosystems</em>. DOI: <a href="https://doi.org/10.1029/2020GC009588">10.1029/2020GC009588</a></li>
            <li>Fiebig et al. 2024. Carbonate clumped isotope values compromised by nitrate-derived NO2 interferent. <em>Chemical Geology</em>. DOI: <a href="https://doi.org/10.1016/j.chemgeo.2024.122382">10.1016/j.chemgeo.2024.122382</a></li>
            <li>Python scripts - modification date:
                <ul>
                    {python_scripts_block}
                </ul>
            <li></li>
            <li><a href="https://github.com/andyschauer/shrekCN">github repository</a></li>
            <li>Data files - <a href="data/{psi_calibrated_session_file}">{psi_calibrated_session_file}</a></li>
            <li><a href="https://isolab.ess.washington.edu/SOPs/psi.php">Running carbonates on Psi.</a></li>
            <li><a href="report.zip">Zip file of entire report directory.</a></strong>.</li>
        </ul>
    </div>
    </body></div></html>"""



# -------------------- WRITE REPORT --------------------
with open(report_page, 'w') as report:
    report.write(header)
    report.write(body)
    report.write(figure_block)
    # [report.write(i) for i in figure_block]
    report.write(footer)
    report.close()
webbrowser.open(report_page)



# -------------------- REPORT ZIP --------------------
print('Archiving report as zip file...')
try:
    shutil.make_archive('report', 'zip', report_path)
    shutil.move('report.zip', os.path.join(report_path, 'report.zip'))
except PermissionError as e:
    print(f"Failed to make zip file due to error - {e}")




