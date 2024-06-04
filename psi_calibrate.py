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
"""

__author__ = "Andy Schauer"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-05-31"
__version__ = "1.0"
__copyright__ = "Copyright 2024, Andy Schauer"
__license__ = "Apache 2.0"
__acknowledgements__ = "K. Huntington, E. Heitmann, S. Mat, V. Ravi, M. Leite. Special thanks to M. Daeron for D47crunch!"



# -------------------- Imports --------------------
import csv
from D47crunch import *
import isolab_lib
import matplotlib.pyplot as pplt
from natsort import natsorted
import os
import shutil
import webbrowser



# -------------------- Functions --------------------




# ---------- Identify session to be loaded ----------
project_path = isolab_lib.get_path("psi", "project")
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
D47crunch_output_analyses_file = f"{session}_D47crunch_output_all_analyses.csv"
D47crunch_output_samples_file = f"{session}_D47crunch_output_summarized.csv"
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

table_of_samples(data47 = psi47,
                  data48 = psi48,
                  dir = f"{session_path}",
                  filename = D47crunch_output_samples_file,
                  save_to_file = True,
                  print_out = True,
                  output = None)


table_of_analyses(data47 = psi47,
                  data48 = psi48,
                  dir = f"{session_path}",
                  filename = D47crunch_output_analyses_file,
                  save_to_file = True,
                  print_out = False,
                  output = None)


# -------------------- combine Psi meta data with D47crunch output data --------------------
h1, d1 = isolab_lib.read_file(os.path.join(session_path, D47crunch_input_file), ',')
h2, d2 = isolab_lib.read_file(os.path.join(session_path, D47crunch_output_analyses_file), ',')

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

headers, data = isolab_lib.read_file(os.path.join(session_path, psi_calibrated_session_file), ',')

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



report_page = os.path.join(session_path, f'{session}_calibration_summary.html')


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

calculation_notes_block = "<li>magic</li>"
refmat_block = "<tr><td>wand</td><td>potions</td><td></td><td></td><td></td><td></td><td></td></tr>"
data_quality_block = "<tr><td>OWLS</td></tr>"
# calculation_notes_block = str([f"<li>{i}</li>" for i in calculation_notes]).replace("[", "").replace("'", "").replace("]", "").replace(", ", "")
# refmat_block = str([f"<tr><td>{eval(i)['names'][0]}</td><td>{eval(i)['material']}</td><td>{eval(i)['d15N']}</td><td>{eval(i)['fractionN']}</td><td>{eval(i)['d13C']}</td><td>{eval(i)['fractionC']}</td><td>{eval(i)['purpose']}</td></tr>" for i in included_isotope_standards]).replace("[", "").replace("'", "").replace("]", "").replace(", ", "")
# data_quality_block = f"""<tr><td>&delta;<sup>15</sup>N</td><td>{np.round(np.std(d15N_AirN2[eval(d15N_std_3)['index']])*2,3)} &permil;</td><td>{np.round(np.mean(d15N_AirN2[eval(d15N_std_3)['index']]) - eval(d15N_std_3)['d15N'], 3)} &permil;</td><td>{d15N_std_3}</td></tr>
                         # <tr><td>&delta;<sup>13</sup>C</td><td>{np.round(np.std(d13C_VPDB[eval(d13C_std_3)['index']])*2,3)} &permil;</td><td>{np.round(np.mean(d13C_VPDB[eval(d13C_std_3)['index']]) - eval(d13C_std_3)['d13C'],3)} &permil;</td><td>{d13C_std_3}</td></tr>
                         # <tr><td><br></td><td> </td><td> </td><td> </td></tr>
                         # <tr><td>&delta;<sup>15</sup>N</td><td>{np.round(d15N_AirN2_residual_std*2, 3)} &permil;</td><td> </td><td>all isotope reference materials</td></tr>
                         # <tr><td>&delta;<sup>13</sup>C</td><td>{np.round(d13C_VPDB_residual_std*2, 3)} &permil;</td><td> </td><td>all isotope reference materials</td></tr>
                         # <tr><td>N quantity</td>
                         #     <td>{np.round(np.nanstd(np.hstack([Nqty[GA1['index']] - GA1['fractionN']*Amount[GA1['index']], Nqty[GA2['index']] - GA2['fractionN']*Amount[GA2['index']], Nqty[SA['index']] - SA['fractionN']*Amount[SA['index']]])*1000)*2, 1)} &micro;g</td>
                         #     <td>{np.round(np.nanmean(np.hstack([Nqty[GA1['index']] - GA1['fractionN']*Amount[GA1['index']], Nqty[GA2['index']] - GA2['fractionN']*Amount[GA2['index']], Nqty[SA['index']] - SA['fractionN']*Amount[SA['index']]]))*1000, 1)} &micro;g</td>
                         #     <td>all isotope reference materials</td></tr>
                         # <tr><td>C quantity</td>
                         #     <td>{np.round(np.nanstd(np.hstack([Cqty[GA1['index']] - GA1['fractionC']*Amount[GA1['index']], Cqty[GA2['index']] - GA2['fractionC']*Amount[GA2['index']], Cqty[SA['index']] - SA['fractionC']*Amount[SA['index']]])*1000)*2, 1)} &micro;g</td>
                         #     <td>{np.round(np.nanmean(np.hstack([Cqty[GA1['index']] - GA1['fractionC']*Amount[GA1['index']], Cqty[GA2['index']] - GA2['fractionC']*Amount[GA2['index']], Cqty[SA['index']] - SA['fractionC']*Amount[SA['index']]]))*1000, 1)} &micro;g</td>
                         #     <td>all isotope reference materials</td></tr>"""

body = f"""
    <body>
    <div class="entire_page">
    <h2>Psi Calibration Report</h2>
    <div class="created-date">Created - {str(dt.datetime.now())}</div>
    <h2>Introduction</h2>
    <div class="text-indent">
        <p>This report is meant to be a stand-alone collection of methods,
        data, scripts, and notes related to calibration of carbonate clumped
        isotope data generated with Psi a Nu Perspective / NuCarb isotope
        ratio mass spectrometer / automated carbonate digestion system.
        You can read more about our implementation of this method on our website
        <a href="https://isolab.ess.washington.edu/laboratory/carbonate-D47.php">
        https://isolab.ess.washington.edu/laboratory/carbonate-D47.php</a>.</p>

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
        blah blah and blah. Each section of data is separated
        by an empty row. The first section of data are the trusted reference materials; the second section of data are trusted samples;
        the third section of data are untrusted. Under the "trust" heading, "1" indicates good, trusted data while "0" indicates poor
        quality data that should probably be distrusted. Up to you if you want to use it. Untrusted data are given the reason for
        distrust. <strong>&iexcl; <a href="report.zip">If you are done analyzing samples, save a copy of the entire report directory elsewhere,
        not just a copy of your data file</a> !</strong></p>
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
            <tr><td>Total number of analyses in session</td><td>{len(set(UID))}</td></tr>
            <tr><td>Total number of standards in session</td><td>Nothing but standards</td></tr>
            <tr><td>Total number of samples in run</td><td>0</td></tr>
            <tr><td><br></td></tr>
            <tr><td>Number of <a href="#distrusted">distrusted analyses in session</a></td><td>0</td></tr>
        </table>
    </div>

    <h2>Reference materials</h2>
    <div class="text-indent">
        <p>describe Bernasconi et al. 2021, Fiebig for GU1</p>
        <table>
            <tr><th>Reference<br>name</th><th>Reference<br>material</th><th>d13C<br>accepted<br>(permil)</th><th>d18O<br>accepted<br>(permil)</th><th>D47</th><th>D48</th><th>Purpose</th></tr>
            {refmat_block}
        </table>
    </div>

    <h2>Data quality</h2>
    <div class="text-indent"><p>Top quality of course.</p>
        <table>
            <tr><th>Parameter</th><th>Precision</th><th>Accuracy</th><th>Reference<br>Material(s)<br>Used</th></tr>
            {data_quality_block}
        </table><br>
    </div>

    <h2>Figures</h2>"""

figure_block = '<div class="clear-both"><img src="/home/aschauer/projects/psi/Results/session_240409/report/figures/D47_plot_mySession.png"><hr></div>'

# python_scripts_block = str([f'<li><a href="python/{key}_REPORT_COPY">{key}</a> - {value}</li>' for key, value in python_scripts.items()]).replace("[", "").replace("'", "").replace("]", "").replace(", ", "")
python_scripts_block = 'python scripts here'

footer = f"""
    <h2 id="refs">References</h2>
    <div class="references">
    <ul>
        <li>Python scripts - modification date:
            <ul>
                {python_scripts_block}
            </ul>
        <li></li>
        <li><a href="https://github.com/andyschauer/shrekCN">github repository</a></li>
        <li>Data files - <a href="data/{psi_calibrated_session_file}">{psi_calibrated_session_file}</a></li>
        <li><a href="https://isolab.ess.washington.edu/laboratory/solid-CN.php">IsoLab's carbon and nitrogen analysis overiew.</a></li>
        <li><a href="https://isolab.ess.washington.edu/SOPs/shrek-cn.php">IsoLab's water carbon and nitrogen analysis method.</a></li>
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
shutil.make_archive('report', 'zip', session_path)
shutil.move('report.zip', os.path.join(session_path, 'report.zip'))




