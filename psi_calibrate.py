"""
This script reads in analysis log files from Psi, a Nu Perspective with an attached NuCarb carbonate device, for 
    calibration and reporting. The brunt of the calibration is via D47crunch, a python library written by M.
    Daeron; everything else in the present script is dressing, presentation, and exporting of salient data. The
    raw data from Psi are initially processed via psi.R with the product being the analysis log files read in here.


Change Log:
    2023.05.26 => created
    2023.06.02 => now we can set a samples flag to 0 (ignore or exclude) in the raw data file and pass
        only flag 1 data to D47crunch.
    2024.02.06 => this script is now called from psi_calibrate.Rmd, an Rmarkdown that oversees the calibration
        of psi data followed by creation of a report.
    2024.04.30 => created way to choose a sample log file to process after dividing sessions into chop2 and chop3
    2024.05.31 => started creating proper report
    2024.07.10 => added GU1 to D48 calibration; added repeated removal of standard as an estimate of accuracy
    2024.08.09 => changed flag to trust
    2024.08.14 => updating to get session read in as variable from data file and used accordingly
    2024.08.16 => can filter based on nchops at the command line; included Emma's addition of equilibrated gases
    2024.09.16 => added Kim et al 2015 18O acid fractionation equation
    2024.12.06 => trying to finish this off, or at least make it a final looking report, rather than the snark it is currently populated with
    2025.01.13 => playing with True values of gases
    2025.01.31 => most of the work needed to get multiple sessions calibrated under a pooled standardization approach works as of this date
    2025.02.06 => finishing up pooled standardization approach; 
    2025.02.07 => created while in Psi time, Emma and I editing
    2025.02.25 => added text above an output table to clarify D47crunch d18O VSMOW values and that they are VSMOW_CO2; I also included a quick blurb about how to convert then to d18O_VSMOW_Carbonate.
    2025.02.27 => added d18O_VSMOWmineral, d18O_VPDBmineral columns to the calibrated data file
"""

__author__ = "Andy Schauer"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2025.02.27"
__copyright__ = "Copyright 2025, Andy Schauer"
__license__ = "Apache 2.0"



# -------------------- IMPORTS --------------------
import argparse
import csv
import datetime
from D47crunch import *
import isolab_lib
import json
import matplotlib.pyplot as pplt
from natsort import natsorted
import os
from pathlib import Path
import shutil
import sys
import time
import webbrowser



# -------------------- REFERENCE MATERIALS --------------------

# CARBONATES
#    - read in json file, make dictionaries for each carbonate, and create an index list within each dictionary
with open(isolab_lib.get_path("psi", "standards"), 'r') as refmat_file:
    refmat = json.load(refmat_file)

refmat_keys = refmat['carbonates'].keys()
for i in refmat_keys:
    globals()[i] = refmat['carbonates'][i]
    globals()[i]['index'] = np.empty(0, dtype="int16")

refmat_keys = refmat['gases']['CDES'].keys()
for i in refmat_keys:
    globals()[i] = refmat['gases']['CDES'][i]
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

# Andy and Kate typing this in from Frontiers review lookup table
t5_D48rf = 0.418308748 # NOTE - this is 5 *C, not 4
t1000_D48rf = -0.002023846

FCEN4['D47'] = t4_D47rf
FFSP4['D47'] = t4_D47rf
FCEN4['D48'] = t5_D48rf
FFSP4['D48'] = t5_D48rf
FC1000['D47'] = t1000_D47rf
FF1000['D47'] = t1000_D47rf
FC1000['D48'] = t1000_D48rf
FF1000['D48'] = t1000_D48rf


# USE THESE CARBONATES AND GASES TO CALIBRATE D47 and D48
d13C_calibration_standards = ['Merck', 'ETH4', 'ETH1']
d18O_calibration_standards = ['Merck', 'ETH4', 'ETH1']
D47_calibration_standards = ['ETH1', 'ETH2', 'ETH3', 'ETH4', 'Merck', 'IAEAC1', 'IAEAC2', 'FC1000', 'FF1000', 'FCEN4', 'FFSP4']
D48_calibration_standards = ['ETH1', 'ETH2', 'ETH4','GU1', 'FC1000', 'FF1000', 'FCEN4','FFSP4']



# -------------------- SETUP and SCRIPTS --------------------
project_path = isolab_lib.get_path("psi", "project")
results_path = f"{project_path}Results/"

script_path = isolab_lib.get_path("psi", "python")
scripts = {'psi.R': '', 'isolab_lib.py': '', 'psi_raw_data_report.py': '', 'psi_calibrate.py': ''}
scripts = {key: (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(f'{script_path}{key}')))) for key, value in scripts.items()}



# ------------------ GET / CREATE REPORT DIRECTORY -------------------------------

report_list = natsorted(os.listdir(f"{project_path}Reports/"))
print('\nChoose your Report directory from the list below:')
[print(f'    {i}') for i in report_list]

report_found = False

while not report_found:
    search_string = input('\nEnter the Report you wish to work on today (or leave blank and hit enter to create a new one): ')
    
    if not search_string:
        search_string = input('\nEnter the name of the Report you wish to create: ')
        Path(f"{project_path}Reports/{search_string}").mkdir(parents=True, exist_ok=True)
        report = search_string
        report_found = True
        print(f"\nReport '{search_string}' created successfully.")
    else:
        # Allow for partial matches in the report list
        matches = [x for x in report_list if search_string in x]

        if len(matches) == 0:
            print(f"\n    '{search_string}' was not found. Please re-enter your Report again.")
        elif len(matches) > 1:
            print(f"\n    More than one Report directory contains '{search_string}'. Please be more specific.")
        else:
            report = matches[0]  # Select the uniquely matched report
            report_found = True
            print(f"\nUsing report directory: {report}")

report_path = f"{project_path}Reports/{report}/report_{int(datetime.datetime.now(datetime.UTC).timestamp())}"
os.mkdir(report_path)
os.mkdir(f"{report_path}/data/")
os.mkdir(f"{report_path}/figures/")
os.mkdir(f"{report_path}/scripts/")



# -------------------- GET SESSION(S) --------------------
session_list = natsorted(os.listdir(results_path))
exclude_from_session_list = ['archive', 'compilations']
session_list = [session for session in session_list if session not in exclude_from_session_list]
print('\nChoose your session(s) from the list below:')
[print(f'    {i}') for i in session_list]
all_sessions_found = False
while all_sessions_found is False:
    search_string = input('\nEnter the session(s) you wish to process (separated by commas): ')
    parsed_string = search_string.split(',')
    sessions_to_process = []

    for session_search in parsed_string:
        is_session = [session_search[0: len(session_search)] in x for x in session_list]
        if len(np.where(is_session)[0]) == 0:
            print(f"\n    {session_search} was not found. Please re-enter your session(s) again.")
            break
        else:
            sessions_to_process.append(session_search)

    if len(parsed_string) != len(sessions_to_process):
        continue

    all_sessions_found = True
    print(f"\n    {len(sessions_to_process)} sessions found to be calibrated.")
    print("\n\n    Now we will find the samplelog files.")



# -------------------- GET SAMPLELOG FILE(S) --------------------
samplelogs_to_process = []
for session in sessions_to_process:
    session_path = f'{results_path}{session}/'

    time.sleep(0.5)
    print(f"    ...looking for all files in {session} that have 'samplelog.csv' in the filename...")
    
    session_file_list = natsorted([file for file in os.listdir(session_path) if "samplelog.csv" in file])
    if len(session_file_list) == 1:
        identified_session_file = 1
        session_samplelog = session_file_list[0]
        print("\n    Got it.")
        time.sleep(0.5)
        print(f'\n    Processing session file {session_samplelog}...')
    else:
        print('\nMore than one session log file found. Choose one from the list below:')
        [print(f'    {i}') for i in session_file_list]
        identified_session_file = 0
        while identified_session_file == 0:
            session_file_search = input('\nEnter the session log file you wish to process: ')
            isdir = [session_file_search[0: len(session_file_search)] in x for x in session_file_list]
            if len(np.where(isdir)[0]) == 1:
                identified_session_file = 1
                session_samplelog = session_file_list[np.where(isdir)[0][0]]
                print(f'    Processing session file {session_samplelog}...')
            else:
                print('\n** More than one session file found. **\n')
                print(f'\n**     Make sure {isdir[0]} does not exist in any other file names. **\n')
                print(f'\n**     The best way to name your custom log files is psi_{session}_MY-SESSION-NAME_samplelog.csv. **\n')
    samplelogs_to_process.append(f"{session_path}{session_samplelog}")

    if len(session_file_list)>1:
        print('\n\n    Next session...')



# -------------------- PRE-LOAD RAW DATA FILE(S) AND REMOVE TRUST=0 --------------------
n_everything = 0  # used in sample accounting below
for i,samplelog in enumerate(samplelogs_to_process):
    with open(samplelog, 'r') as datafile:
        reader = csv.reader(datafile)
        headers = next(reader)

        rows = []
        for row in reader:
            n_everything += 1
            if int(row[4]) == 1:
                rows.append(row)

    headers.append('Session')
    for row in rows:
        row.append(sessions_to_process[i])

    with open(f"{report_path}/data/{sessions_to_process[i]}_D47crunch_input_file.csv", 'w', newline='') as datafile:
        writer = csv.writer(datafile)
        writer.writerow(headers)
        writer.writerows(rows)



# -------------------- LOAD INPUT FILES AND CONVERT THEM TO DICTIONARIES --------------------
data_as_dictionary = {}
for i, session in enumerate(sessions_to_process):
    data_as_dictionary[sessions_to_process[i]] = D47data()
    data_as_dictionary[sessions_to_process[i]].read(f"{report_path}/data/{sessions_to_process[i]}_D47crunch_input_file.csv")

all_sessions = []
for lst in data_as_dictionary.values():
    all_sessions += lst



# -------------------- D47crunch --------------------
psi47 = D47data(all_sessions)

rxnT = 70
psi47.ALPHA_18O_ACID_REACTION = np.exp((3.48 * (1000/(rxnT+273.15)) - 1.47)/1000) # Kim et al 2015
psi47.Nominal_d13C_VPDB = {eval(i)['names'][0]: eval(i)['d13C'] for i in d13C_calibration_standards}
psi47.Nominal_d18O_VPDB = {eval(i)['names'][0]: eval(i)['d18O_vpdb'] for i in d18O_calibration_standards}
psi47.Nominal_D47 = {eval(i)['names'][0]: eval(i)['D47'] for i in D47_calibration_standards}
psi47.LEVENE_REF_SAMPLE = 'IAEA-C2'

psi47.wg()
psi47.crunch()
psi47.standardize()
psi47.summary(verbose=True, save_to_file=True)


psi48 = D48data(all_sessions)
psi48.ALPHA_18O_ACID_REACTION = np.exp((3.48 * (1000/(rxnT+273.15)) - 1.47)/1000) # Kim et al 2015
psi48.Nominal_d13C_VPDB = {eval(i)['names'][0]: eval(i)['d13C'] for i in d13C_calibration_standards}
psi48.Nominal_d18O_VPDB = {eval(i)['names'][0]: eval(i)['d18O_vpdb'] for i in d18O_calibration_standards}
psi48.Nominal_D48 = {eval(i)['names'][0]: eval(i)['D48'] for i in D48_calibration_standards}
psi48.LEVENE_REF_SAMPLE = 'IAEA-C2'

psi48.wg()
psi48.crunch()
psi48.standardize()
psi48.summary(verbose=True)

table_of_sessions(data47 = psi47,
                  data48 = psi48,
                  dir = f"{report_path}/data/",
                  filename = "D47crunch_table_of_sessions.csv",
                  save_to_file = True,
                  print_out = True,
                  output = None)

table_of_samples(data47 = psi47,
                 data48 = psi48,
                 dir = f"{report_path}/data/",
                 filename = "D47crunch_table_of_samples.csv",
                 save_to_file = True,
                 print_out = True,
                 output = None)

table_of_analyses(data47 = psi47,
                  data48 = psi48,
                  dir = f"{report_path}/data/",
                  filename = "D47crunch_table_of_analyses.csv",
                  save_to_file = True,
                  print_out = False,
                  output = None)


psi47.plot_sessions(dir=f"{report_path}/figures/", figsize=(8, 8), filetype='png', dpi=100)
psi48.plot_sessions(dir=f"{report_path}/figures/", figsize=(8, 8), filetype='png', dpi=100)

psi47.plot_distribution_of_analyses(dir=f"{report_path}/figures/", filename='time_distribution', figsize=(8, 8), subplots_adjust=(0.02, 0.13, 0.85, 0.8), dpi=100)



# ------------------------- COMBINE SESSIONS INPUT FILES INTO SINGLE FILE --------------------------
combined_sessions_file = f"{report_path}/data/D47crunch_session_input_files_combined.csv"
with open(combined_sessions_file, mode="w", newline="") as output:
    writer = csv.writer(output)
    
    headers = None
    for i, e in enumerate(sessions_to_process):
        inputfile = f"{report_path}/data/{sessions_to_process[i]}_D47crunch_input_file.csv"
        with open(inputfile, mode="r") as file:
            reader = csv.reader(file)
            if headers is None:
                headers = next(reader)
                writer.writerow(headers)
            else:
                next(reader)  # Skip the header row for subsequent files
            for row in reader:
                writer.writerow(row)



# -------------------- combine Psi input data with D47crunch output data --------------------

h1, d1 = isolab_lib.read_file(os.path.join(report_path, "data/", combined_sessions_file), ',')
h2, d2 = isolab_lib.read_file(os.path.join(report_path, "data/", "D47crunch_table_of_analyses.csv"), ',')
psi_calibrated_file = 'psi_calibrated_file.csv'

# D47crunch leaves d18O_VSMOW as the d18O of the CO2 relative to VSMOW. Some prefer to think about their carbonate
#    and not the resulting CO2. As such, the d18O of the carbonate mineral relative to both VSMOW and VPDB are added
#    to the final calibrated file.
d2['d18O_VSMOWmineral'] = [str((np.double(i)-8.709)/1.008709) for i in d2['d18O_VSMOW']] # assumes 70*C acid and equation 6 from Kim et al 2015: acid.alpha70 = exp((3.48 * (10^3/(70+273.15)) - 1.47)/1000) which is equal to 1.00870904256074
d2['d18O_VPDBmineral'] = [str(np.double(i)*0.97001 - 29.99) for i in d2['d18O_VSMOWmineral']]


d3 = {**d1, **d2}
h3=list(d3.keys())
h3[-2] = 'D47_iCDES'
# h3[-1] = 'D48_iCDES'

with open(f"{session_path}{session_samplelog}", 'r') as datafile:
    reader = csv.reader(datafile)
    headers0 = next(reader)
    rows0 = [row for row in reader if int(row[4]) == 0]


with open(os.path.join(report_path, 'data/', psi_calibrated_file), 'w', newline='') as cal_sesh_file:
        datawriter = csv.writer(cal_sesh_file)
        datawriter.writerow(d3.keys())
        datawriter.writerows(zip(*d3.values()))
        datawriter.writerow("")
        datawriter.writerows(rows0)



# ------------------ read in final calibrated file for interactive plotting etc -------------------------------

headers, data = isolab_lib.read_file(os.path.join(report_path, 'data/', psi_calibrated_file), ',')

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




# ---------- REPORT BITS ---------- 

print('Making html page...')
report_page = os.path.join(report_path, '_calibration_summary.html')

style = """
    header {
        border: 0px solid black;
        font-size:1.3rem;
        height:4rem;
        width:100%;
    }

    body{
        background-color: white;
        font-family: Arial, Helvetica, sans-serif;
        height: 100vh;
    }

    .entire_page{
        margin: 10px;
    }

    img{
        max-width: 100%;
        padding:2rem;
    }

    h2,h3{
        margin-bottom: 0rem;
        padding-top: 1.0rem;
        padding-bottom: 0rem;
    }

    .created-date{
        float: right;
        font-size: 0.8rem;
    }

    .text-indent{
        margin-left:2rem;
    }


    /* tables */
    table{
        border-bottom: 1px solid black;
    }
    table th{
        border-bottom: 1px solid black;
        text-align: center;
    }
    table td{
        padding-left: 1rem;
        padding-right: 1rem;
        text-align: center;
    }
    table tr:hover {
        background-color: #ffff99;
    }


    .figThumb{
        float:left;
        margin: 1rem;
        padding:1rem;
    }

    .figThumb img{
        height:auto;
        margin: 0;
        padding: 0;
        width:200px;
    }

    .clear-both{
        clear: both;
        margin-bottom: 2rem;
    }

    .largeFigs{
        border:  0px solid black;
        display:block;

    }

    .largeFigs img{
        border:  0px solid black;
        margin: 0;
        max-width: 100%;
        padding: 0;
    }

    .caption{
        font-size: 1.2rem;
        margin-left: 2rem;
        width: 1200px;
    }

    .thumbcaption{
        margin-left: 0;
        padding: 0;
        text-align: center;
        width: 200px;
    }

    .references li{
        padding-bottom:0.4rem;
    }
"""

header = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <!-- py by Andrew Schauer -->
        <meta http-equiv="Content-Type" content="text/html charset=UTF-8" />
        <meta name="viewport" content="width=device-width,initial-scale=1">
        <style>{style}</style>
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

summary_table_of_sessions = "<table>\n"
with open(f"{report_path}/data/D47crunch_table_of_sessions.csv", "r") as f:
    reader = csv.reader(f)
    for i, row in enumerate(reader):
        if i == 0:
            summary_table_of_sessions += "<tr>" + "".join(f"<th>{col}</th>" for col in row) + "</tr>\n"
        else:
            summary_table_of_sessions += "<tr>" + "".join(f"<td>{col}</td>" for col in row) + "</tr>\n"
summary_table_of_sessions += "</table>"

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
with open(f"{report_path}/data/D47crunch_table_of_samples.csv", "r") as f:
    reader = csv.reader(f)
    for i, row in enumerate(reader):
        if i == 0:
            replicate_summary_table += "<tr>" + "".join(f"<th>{col}</th>" for col in row) + "</tr>\n"
        else:
            replicate_summary_table += "<tr>" + "".join(f"<td>{col}</td>" for col in row) + "</tr>\n"
replicate_summary_table += "</table>"


figure_block_1 = [f"""<div class="clear-both"><img src="figures/D47_plot_{session}.png"><hr></div>""" for session in sessions_to_process]
figure_block_2 = [f"""<div class="clear-both"><img src="figures/D48_plot_{session}.png"><hr></div>""" for session in sessions_to_process]
figure_block_3 = f"""<div class="clear-both"><img src="figures/time_distribution.png"><hr></div>"""

body = f"""
    <body>
    <div class="entire_page">
    <h2>Psi Calibration Report - {', '.join(sessions_to_process)}</h2>
    <div class="created-date">Created - {str(datetime.datetime.now())}</div>
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
        <p>This technical stuff is fine and all but where are <strong><a href="data/psi_calibrated_file">my data</a></strong>?
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

        <div class="clear-both"></div>
        <h3>Summary Table of Sessions</h3>
        {summary_table_of_sessions}
        <div class="figThumb"><strong>&Delta;<sub>47</sub> Repeatability</strong>
        {D47_summary_table}</div>

        <div class="figThumb"><strong>&Delta;<sub>48</sub> Repeatability</strong>
        {D48_summary_table}</div>

        <div class="clear-both"></div>

        <h3>Replicate Summary Table</h3>
        <p>It is important to note that the d18O_vsmow values in the below table, and in all of D47crunch, are that of the CO2, not the carbonate. To convert these d18O_vsmow CO2 values to d18O_vsmow carbonate values
        you need to convert it based on the acid temperature fractionation factor. Since Psi is at 70*C, we are going to use equation 6 from Kim et al 2015:
        acid.alpha70 = exp((3.48 * (10^3/(70+273.15)) - 1.47)/1000) which is equal to 1.00870904256074. d18O_VSMOW_Carbonate = (d18O_VSMOW_CO2-8.709)/1.008709.</p>
        <p>{replicate_summary_table}</p>
    </div>



    <h2>Figures</h2>"""



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
            <li>Data files - <a href="data/psi_calibrated_file">psi_calibrated_file</a></li>
            <li><a href="https://isolab.ess.washington.edu/SOPs/psi.php">Running carbonates on Psi.</a></li>
            <li><a href="report.zip">Zip file of entire report directory.</a></strong>.</li>
        </ul>
    </div>
    </body></div></html>"""



# -------------------- WRITE REPORT --------------------
with open(report_page, 'w') as report:
    report.write(header)
    report.write(body)
    [report.write(i) for i in figure_block_1]
    [report.write(i) for i in figure_block_2]
    report.write(figure_block_3)
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




