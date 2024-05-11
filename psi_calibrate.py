"""
This script is really just dressing for M. Daeron's D47crunch carbonate clumped isotope python library.

Version 0.1 - 2023.05.26 => created
Version 0.5 - 2023.06.02 => now we can set a samples flag to 0 (ignore or exclude) in the raw data file and pass
    only flag 1 data to D47crunch.
Version 0.6 - 2024.02.06 => this script is now called from psi_calibrate.Rmd, an Rmarkdown that oversees the calibration
    of psi data followed by creation of a report.
Version 0.7 - 2024.04.30 => created way to choose a sample log file to process after dividing sessions into chop2 and chop3
"""

__author__ = "Andy Schauer"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-04-30"
__version__ = "0.7"
__copyright__ = "Copyright 2024, Andy Schauer"
__license__ = "Apache 2.0"
__acknowledgements__ = "K. Huntington, E. Heitmann, S. Mat, Vigash. Special thanks to M. Daeron for D47crunch!"



# -------------------- Imports --------------------
import csv
from D47crunch import *
import lab
import matplotlib.pyplot as pplt
from natsort import natsorted
import os
import shutil


# -------------------- Functions --------------------
def get_path(desired_path):
    """Make your life easier with this section. These are the paths that seem to change depending on the computer we are working on."""
    psi_path_file = os.path.join(os.getcwd(), 'psi_path.txt')
    if os.path.isfile(psi_path_file):
        # print(' :-) Using existing psi path file for a warm and fuzzy experience. (-:')
        with open(psi_path_file, 'r') as ppf:
            home, python_path, project_path, standards_path = ppf.readline().split(',')
            python_path = home + python_path
            project_path = home + project_path
            standards_path = home + standards_path

    else:
        python_path_check = False
        project_path_check = False
        standards_path_check = False
        print(' )-: Picarro path file does not exist yet. :-(')
        print(" Let's make one... :-| ")
        while python_path_check is False:
            python_path = input(f'Enter the current path to the psi python scripts. Perhaps it is {os.getcwd()}. ')
            if os.path.isdir(python_path):
                python_path_check = True
                if python_path[-1] != '/':
                    python_path += '/'
            else:
                print(f'oops, try typing that in again (you typed {python_path}): ')

        while project_path_check is False:
            project_path = input('Enter the current path to your projects: ')
            if os.path.isdir(project_path):
                project_path_check = True
                if project_path[-1] != '/':
                    project_path += '/'
            else:
                print(f'oops, try typing that in again (you typed {project_path}): ')

        while standards_path_check is False:
            standards_path = input('Enter the current path and filename to your reference materials file: ')
            if os.path.isfile(standards_path):
                standards_path_check = True
            else:
                print(f'oops, try typing that in again (you typed {standards_path}): ')

        with open(psi_path_file, 'w') as ppf:
            ppf.write(f'{python_path},{project_path},{standards_path}')

    if desired_path == "project":
        return project_path
    elif desired_path == "python":
        return python_path
    elif desired_path == "standards":
        return standards_path
    else:
        unknown_path = input('Enter the path to your project: ')
        return unknown_path



# -------------------- paths --------------------
project_path = get_path("project")

# ---------- Identify session to be loaded ----------
session_list = natsorted(os.listdir(project_path))
print('\nChoose a session from the list below:')
[print(f'    {i}') for i in session_list]
identified_session = 0
while identified_session == 0:
    session_search = input('Enter the session you wish to process: ')
    isdir = [session_search[0: len(session_search)] in x for x in session_list]
    if len(np.where(isdir)[0]) == 1:
        identified_session = 1
        session = session_list[np.where(isdir)[0][0]]
        print(f'    Processing session {session}...')
    else:
        print('\n** More than one session found. **\n')

session_path = f'{project_path}/{session}/'
report_path = f"{session_path}report/"


session_file_list = natsorted([file for file in os.listdir(session_path) if "samplelog.csv" in file])
print('\nChoose a session log file from the list below:')
[print(f'    {i}') for i in session_file_list]
identified_session_file = 0
while identified_session_file == 0:
    session_file_search = input('Enter the session log file you wish to process: ')
    isdir = [session_file_search[0: len(session_file_search)] in x for x in session_file_list]
    if len(np.where(isdir)[0]) == 1:
        identified_session_file = 1
        session_samplelog = session_file_list[np.where(isdir)[0][0]]
        print(f'    Processing session file {session_samplelog}...')
    else:
        print('\n** More than one session file found. **\n')

D47crunch_input_file = f"psi_{session}_D47crunch_input.csv"
D47crunch_output_file = f"{session}_D47crunch_output.csv"
psi_calibrated_session_file = f"psi_{session}_calibrated.csv"



# -------------------- preload rawdata file and remove flag0 samples --------------------
with open(f"{session_path}{session_samplelog}", 'r') as datafile:
    reader = csv.reader(datafile)
    headers = next(reader)
    rows = [row for row in reader if int(row[4]) == 1]

with open(f"{session_path}{D47crunch_input_file}", 'w', newline='') as datafile:
    writer = csv.writer(datafile)
    writer.writerow(headers)
    writer.writerows(rows)


# -------------------- define standards to use --------------------
D47data.Nominal_D4x = {
    a: D47data.Nominal_D4x[a]
    for a in ['ETH-1', 'ETH-2', 'ETH-3', 'ETH-4', 'IAEA-C1', 'MERCK']
    }

# to get used standards type  psi47.Nominal_D47

D48data.Nominal_D4x = {
    a: D48data.Nominal_D4x[a]
    for a in ['ETH-1', 'ETH-2', 'ETH-3', 'ETH-4']
    }

# -------------------- main interface with D47crunch --------------------
psi47 = D47data(session=session)
psi47.read(f"{session_path}{D47crunch_input_file}")
psi47.wg()
psi47.crunch()
psi47.standardize()

psi48 = D48data(psi47, session=session)
psi48.crunch()
psi48.standardize()


psi47.summary(verbose=True)
psi48.summary(verbose=True)

table_of_samples(psi47, psi48)

table_of_analyses(data47 = psi47,
                  data48 = psi48,
                  dir = f"{session_path}",
                  filename = D47crunch_output_file,
                  save_to_file = True,
                  print_out = False,
                  output = None)


# -------------------- combine Psi meta data with D47crunch output data --------------------
h1, d1 = lab.read_file(os.path.join(session_path, D47crunch_input_file), ',')
h2, d2 = lab.read_file(os.path.join(session_path, D47crunch_output_file), ',')

d3 = {**d1, **d2}
h3=list(d3.keys())
h3[-2] = 'D47_iCDES'
h3[-1] = 'D48_iCDES'


with open(f"{session_path}{session_samplelog}", 'r') as datafile:
    reader = csv.reader(datafile)
    headers0 = next(reader)
    rows0 = [row for row in reader if int(row[4]) == 0]


with open(os.path.join(session_path, psi_calibrated_session_file), 'w', newline='') as cal_sesh_file:
        datawriter = csv.writer(cal_sesh_file)
        datawriter.writerow(d3.keys())
        datawriter.writerows(zip(*d3.values()))
        datawriter.writerow("")
        datawriter.writerows(rows0)




# ------------------ create report directory -------------------------------
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

psi47.plot_sessions(dir=f"{report_path}figures/")
psi48.plot_sessions(dir=f"{report_path}figures/")





# ------------------ read in final calibrated file for interactive plotting etc -------------------------------

headers, data = lab.read_file(os.path.join(session_path, psi_calibrated_session_file), ',')

original_data = data.copy()
flag0_indices = [i for i, e in enumerate(data['flag']) if int(e) == 0]
flag1_indices = [i for i, e in enumerate(data['flag']) if int(e) == 1]
for header in headers[:-1]:
    data[header] = [data[header][index] for index in flag1_indices]

numlist = ['UID', 'flag', 'd13C', 'd13Csd', 'd13Cn', 'd13Cslope', 'd13Cr2', 'd13Cpval', 'd18O', 'd18Osd', 'd18On', 'd18Oslope', 'd18Or2', 'd18Opval',
           'd45', 'd45sd', 'd45n', 'd45slope', 'd45r2', 'd45pval', 'd46', 'd46sd', 'd46n', 'd46slope', 'd46r2', 'd46pval',
           'd47', 'd47sd', 'd47n', 'd47slope', 'd47r2', 'd47pval', 'D47raw', 'D47sd', 'D47n', 'D47slope', 'D47r2', 'D47pval',
           'd48', 'd48sd', 'd48n', 'd48slope', 'd48r2', 'd48pval', 'D48raw', 'D48sd', 'D48n', 'D48slope', 'D48r2', 'D48pval',
           'd49', 'd49sd', 'd49n', 'd49slope', 'd49r2', 'd49pval', 'D49raw', 'D49sd', 'D49n', 'D49slope', 'D49r2', 'D49pval',
           'R45sam_mean', 'R45sam_sd', 'R45wg_mean', 'R45wg_sd', 'R45sam_wg', 'R46sam_mean', 'R46sam_sd', 'R46wg_mean', 'R46wg_sd', 'R46sam_wg',
           'R47sam_mean', 'R47sam_sd', 'R47wg_mean', 'R47wg_sd', 'R47sam_wg', 'R48sam_mean', 'R48sam_sd', 'R48wg_mean', 'R48wg_sd', 'R48sam_wg',
           'R49sam_mean', 'R49sam_sd', 'R49wg_mean', 'R49wg_sd', 'R49sam_wg', 'unix_time', 'mass', 'sample_pressure', 'yield_mbar_ug', 'residual_gas',
           'nchops', 'min_ref_beam', 'max_ref_beam', 'total_blocks', 'total_cycles', 'initial_sample_beam', 'pre_balance_sample_beam',
           'nu_balance', 'sam_max44nA', 'sam_mean44nA', 'sam_min44nA', 'sam_max47nA', 'sam_mean47nA', 'sam_min47nA', 'min_balance', 'mean_balance', 'max_balance',
           'sam_decay', 'wg_decay']

strlist = set(headers) - set(numlist)

for i in numlist:
    globals()[i] = np.asarray(data[i], dtype=float)

for i in strlist:
    globals()[i] = np.asarray(data[i])

