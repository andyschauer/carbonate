#!/usr/bin/Rscript
#
# psi.R
# 
# aschauer@uw.edu - 2023-05-21
#
# psi.R reads in the raw data files from psi and writes them to a csv file to be processed with psi_calibrate.R
#
# Change Log:
#    2023-05-21 - saved archive of previous version in an effort to simplify psi.R and make it streamlined with D47crunch.py
#    2023-05-23 - adding back in items that were witheld from the D47crunch file while I thought about data work flow
#    2023-05-24 - added psi_path.txt, added more quality control flagging
#    2023-06-01 - fixed issue regarding keeping track of processed files across platforms
#    2023-06-09 - removing separate file for flag0 because now the psi_daeron.py script takes care of flag 0 analyses;
#        also changed misbalance flagging thresholds from +0.56 and -0.24 to +/- 1.0.
#    2023-06-13 - added analysis time to display summary
#    2023-07-21 - added D47 and D48 calculation, added session selection prompt
#    2023-07-23 - delta and permil symbols, ensured data output to csv file, uploaded to bob data server
#    2023-08-30 - created reprocess and reprocess_all, combined outlier and linear model functions into 'summarize', 
#    2023-09-12 - continued development, added min and max ref beam to track reference refills
#    2023-09-21 - making a workflow that is R only, coming back to Daeron code later...
#    2023-09-26 - combined session_230629 specific correction code with mainstream psi.R
#    2023-09-27 - updated chop threshold and ensured bellows mode is still output
#    2023-10-18 - updated residual_gas threshold from 0.003 to 0.005 based on historical data; added messaging when a data file has .fail tag
#    2023-12-05 - updated while trying to get ref gases
#    2024-01-16 - no more flag=0 when residual gas is above certain threshold, converting all masses to micrograms, 
#        added delta value comments when standard deviation is high, using reference_materials.json file,
#    2024-01-18 - found bug in standardizing nomenclature of percent carbonate standards, removed for now but will have to make sure names are entered correctly manually
#    2024-02-12 - sort of combining psi_calibrate.R with psi.R, not to calibrate the data but to read in the current sample log and make plots so folks can
#        assess outliers and such before running psi_calibrate.py, the python script that uses M. Daeron's D47crunch python library.
#    2024-03-07 - abandoning Rmd for now in favor of making R run a python code to create a report. I know more python and
#        need to get this working for now. R and Rmd will happen later. Also removed psi_lib.R since this is now the only R script and pasted the contents of psi_lib
#        in here.
#    2024-03-18 - removed integration_summary and simply take the mean, and now the standard deviation; also, renames high_res to all_cycles
#    2024-03-19 - uploaded to bob, all_cycles saved as hdf5 file
#    2024-04-02 - making sure it works on bob and zax
#    2024-04-26 - found bug in all_cycles block where it errors out if not cycles exist, fixed with if cycles exist
#    2024-04-30 - updated residual gas threshold from 0.005 to 0.010
#    2024-05-18 - added bulk vs clumped conditional parts so psi.R can be used for both clumped and bulk runs
#    2024-08-06 - sometimes a Data file did not "fail" but does not contain any data, which used to cause the script to crash, this update should take care of that edge case 
#    2024-08-07 - added Session column to samplelog file, which should make D47crunch happier; changed wg_d18O_vpdb to 1.7515 from -7.4; changed 'sn' to UID throughout;
#    2024-08-09
#    2024-11-18 - updating after Ashley visit; bellows mode now works, so i have removed the flagging related to bellows mode,
#                 I also modified first two figures to make both on scale when the wg and same are severely misbalanced.
#    2025-01-13 - had to update the date string format from the Data files due to Nu change
#    2025-02-27 - updated values of reference gas as part of our recent VPDB-VSMOW discussions; ensuring values are consistent with the most thorough measurement from 2022-12-27; and in consideration with D47crunch estimates
#
#
# ToDo:
#    - outlier function needs work
#        - if one thing has an outlier, are all things pitched?
#    - save figures and such to a report for each sample
#    - while we can reprocess a single sample, we end up with duplicates in the log file, need to deal with this.



# -------------------- Introduction ------------------------

psi_R_version <- "psi.R - v2025.02.27"
message("\n")
message("Running ", psi_R_version)
message("\n")



# -------------------- Libraries ------------------------
library(rhdf5)
library(jsonlite)
library(lubridate)
library(plotly)
library(pracma)
library(rmarkdown)
library(sys)



# -------------------- Functions -------------------- 

calculate_ratios <- function(R45, R46){
    # Calculate component ratios (R13, R17, R18)
    # If vpdb.who is set to Daeron_ChemGeol_2016 in polly_standards.m, the calculate exactly
    #    as Mathieu's D47lib.py, mostly under the 'bulk_and_clumping_deltas' function.
    # copy and paste from Daeron_ChemGeol_2016 python script
    A <- -3 * K^2 * vpdb_R18^(2*a)
    B <- 2 * K * R45 * vpdb_R18^a
    C <- 2 * vpdb_R18
    D <- - R46
    aa <- A*a*(2*a-1) + B*a*(a-1)/2
    bb <- 2*A*a + B*a + C
    cc <- A+B+C+D
    d18O_PDBCO2 <- 1000 * ( -bb + ( bb^2 - 4 * aa * cc )^.5 ) / ( 2 * aa )
    R18 <- ( 1 + d18O_PDBCO2 / 1000 ) * vpdb_R18
    R17 <- K * R18 ^ a
    R13 <- R45 - 2 * R17
    return(c(R13, R17, R18))
}


get_session <- function(){
    session_list <- list.dirs(paste(project_dir, results_dir, sep='/'), full.names = FALSE, recursive = FALSE)
    found_session <- FALSE
    while (found_session == FALSE) {
      message("Here is a list of all sessions in the Results directory:")
      for (i in session_list[session_list != "archive"]){
        message(paste('    ', i, sep=' '))
      }
      message("\n")
      session_name <- readline("Enter the session you are working with: ")

      if (length(which(session_list==session_name))){
        found_session <- TRUE
      } else {
        message("\n")
        message(" ** ooops, the directory name you entered is not in the list, try again... ")
        message("\n")
        Sys.sleep(0.75)
      }
    }
    samplelog <- paste(project_dir, results_dir, session_name, paste("psi_", session_name, "_samplelog.csv", sep=""), sep="/")

    return(list(name=session_name, samplelog=samplelog))
}


summarize <- function(data, outlier_sigma){
    # Outlier test and removal
    # z = (x - m) / s     x = individual datum, m = mean of data set, s = standard deviation of dataset
    check_for_outliers <- TRUE
    adi <- c(1:length(data)) # all data index
    pdi <- adi # population data index
    while (check_for_outliers){
        m <- mean(data[pdi])
        s <- sd(data[pdi])
        z <- (data[pdi] - m) / s
        odi <- which(abs(z) > outlier_sigma) # outlier data index
        if (length(odi)>0){
          pdi <- pdi[-odi]
          # go again
        } else {
          check_for_outliers = FALSE
          odi <- setdiff(adi,pdi)
        }
    }
    # Fit a linear model of data to the cycle index, which is in chronological order
    model <- lm(data[pdi]~pdi)

    # returned items are:
        # m=mean(data[pdi])                =>          mean of data with outliers removed
        # s=sd(data[pdi])                  =>          standard deviation of data with outliers removed
        # n=length(data[pdi])              =>          number of cycles with outliers removed
        # data=data                        =>          all data
        # adi=adi                          =>          all data index (1:n)
        # pdi=pdi                          =>          index of data with outliers removed
        # odi=odi                          =>          index of outliers
        # slope=summary(model)[[4]][2,1]   =>          slope of data with index (which is cycles in chronological order)
        # R2=summary(model)[[8]]           =>          R-squared of the data with outliers removed
        # P=summary(model)[[4]][2,4]       =>          P value of slope of data with outliers removed
    return(list(m=mean(data[pdi]), s=sd(data[pdi]), n=length(data[pdi]), data=data, adi=adi, pdi=pdi, odi=odi, slope=summary(model)[[4]][2,1], R2=summary(model)[[8]], P=summary(model)[[4]][2,4]))
}




#  -------------------- Directory and File Setup -------------------- 

paths <- strsplit(readLines('psi_path_R.txt'), ",")
R_dir <- paths[[1]][1]
project_dir <- paths[[1]][2]
refmat_file <- paths[[1]][3]
py_dir <- paths[[1]][4]
results_dir <- "Results"


# -------------------- Constants ------------------------

# from Daeron_ChemGeol_2016
a <- 0.528
vpdb_R13 <- 0.01118
vsmow_R18 <- 0.00200520
vsmow_R17 <- 0.00038475
vpdb_R18 <- vsmow_R18 * 1.03092 * 1.01025  # These two assume an acid temperature of 25 *C
vpdb_R17 <- 0.00039310
# vpdb_R18 <- vsmow_R18 * 1.03092 * 1.008709
K <- 0.010224511345758351
# vpdb_R17 <- K/(vpdb_R18^-a)
Kcalc_vpdb <- vpdb_R17*vpdb_R18^-a
Kcalc_vsmow <- vsmow_R17*vsmow_R18^-a

# Psi working gas
#    Working gas accepted values 
wg_d13C_vpdb <- -10.16  # VPDB; Measured using polly against polly working gas on 2022.12.27; polly wg had been calibrated directly to NBS19 and VSMOW as part of the Saenger et al 17O effort
wg_d18O_vpdb <- 2.15  # VPDB_CO2; Measured using polly against polly working gas on 2022.12.27; polly wg had been calibrated directly to NBS19 and VSMOW as part of the Saenger et al 17O effort
wg_R13 <- (wg_d13C_vpdb/1000+1)*vpdb_R13
wg_R18 <- (wg_d18O_vpdb/1000+1)*vpdb_R18
wg_R17 = vpdb_R17*(wg_R18/vpdb_R18)^a

# Stochastic working gas
#    R45 and R46 are not stochastic, but in the interest of common nomenclature, the "_sto" tag is maintained.
wg_R45sto <- wg_R13 + 2*wg_R17
wg_R46sto <- 2*wg_R18 + 2*wg_R13*wg_R17 + (wg_R17^2)
wg_R47sto <- 2*wg_R13*wg_R18 + 2*wg_R17*wg_R18 + wg_R13*(wg_R17^2)
wg_R48sto <- wg_R18^2 + 2*wg_R13*wg_R17*wg_R18
wg_R49sto <- wg_R13*wg_R18^2


vpdb_vsmowslap <- 1.03092  # Kim_GCA_2015
acid <- list()
acid$d18O_alpha25 <- exp((3.48 * (10^3/(25+273.15)) - 1.47)/1000)  # 1.01025419478517
acid$d18O_alpha70 <- exp((3.48 * (10^3/(70+273.15)) - 1.47)/1000)  # 1.00870904256074
acid$d18O_alpha90 <- exp((3.48 * (10^3/(90+273.15)) - 1.47)/1000)  # 1.00814581509316



# -------------------- Read in reference materials -------------------- 
refmats <- fromJSON(refmat_file)

CarbStd_names <- list(c(refmats$carbonates$ETH1$names, refmats$carbonates$ETH2$names,
                        refmats$carbonates$ETH3$names, refmats$carbonates$ETH4$names,
                        refmats$carbonates$IAEAC1$names, refmats$carbonates$IAEAC2$names,
                        refmats$carbonates$Merck$names, refmats$carbonates$Coral$names,
                        refmats$carbonates$GU1$names))

RefGas_names <- list(c(refmats$gases$CDES$FC1000$names, refmats$gases$CDES$FF1000$names))



# -------------------- Setup -------------------- 

session <- get_session()

comment <- ""


samplelog_headers_clumped <- c("UID", "analysis_time", "Sample", "material_type", "trust","comment", "Session",
                               "d13C", "d13Csd", "d13Cn", "d13Cslope", "d13Cr2", "d13Cpval",
                               "d18O", "d18Osd", "d18On", "d18Oslope", "d18Or2", "d18Opval",
                               "d45", "d45sd", "d45n", "d45slope", "d45r2", "d45pval",
                               "d46", "d46sd", "d46n", "d46slope", "d46r2", "d46pval",
                               "d47", "d47sd", "d47n", "d47slope", "d47r2", "d47pval",
                               "D47raw", "D47sd", "D47n", "D47slope", "D47r2", "D47pval",
                               "d48", "d48sd", "d48n", "d48slope", "d48r2", "d48pval",
                               "D48raw", "D48sd", "D48n", "D48slope", "D48r2", "D48pval",
                               "d49", "d49sd", "d49n", "d49slope", "d49r2", "d49pval",
                               "D49raw", "D49sd", "D49n", "D49slope", "D49r2", "D49pval",
                               "R45sam_mean", "R45sam_sd", "R45wg_mean", "R45wg_sd", "R45sam_wg",
                               "R46sam_mean", "R46sam_sd", "R46wg_mean", "R46wg_sd", "R46sam_wg",
                               "R47sam_mean", "R47sam_sd", "R47wg_mean", "R47wg_sd", "R47sam_wg",
                               "R48sam_mean", "R48sam_sd", "R48wg_mean", "R48wg_sd", "R48sam_wg",
                               "R49sam_mean", "R49sam_sd", "R49wg_mean", "R49wg_sd", "R49sam_wg",
                               "bkgnd_44", "bkgnd_45", "bkgnd_46", "bkgnd_47", "bkgnd_48", "bkgnd_49",
                               "unix_time", "batch_name", "Nu_software_Version","psi_R_versions",
                               "mass", "sample_pressure", "yield_mbar_ug", "residual_gas", "nchops", "min_ref_beam", "max_ref_beam",
                               "total_blocks", "total_cycles",
                               "initial_sample_beam", "pre_balance_sample_beam", "nu_balance",
                               "sam_max44nA", "sam_mean44nA", "sam_min44nA", "sam_max47nA", "sam_mean47nA", "sam_min47nA",
                               "min_balance", "mean_balance", "max_balance", "sam_decay", "wg_decay")

samplelog_headers_bulk <- c("UID", "analysis_time", "Sample", "material_type", "trust","comment", "Session",
                           "d13C", "d13Csd", "d13Cn", "d13Cslope", "d13Cr2", "d13Cpval",
                           "d18O", "d18Osd", "d18On", "d18Oslope", "d18Or2", "d18Opval",
                           "d45", "d45sd", "d45n", "d45slope", "d45r2", "d45pval",
                           "d46", "d46sd", "d46n", "d46slope", "d46r2", "d46pval",
                           "R45sam_mean", "R45sam_sd", "R45wg_mean", "R45wg_sd", "R45sam_wg",
                           "R46sam_mean", "R46sam_sd", "R46wg_mean", "R46wg_sd", "R46sam_wg",
                           "unix_time", "batch_name", "Nu_software_Version","psi_R_versions",
                           "mass", "sample_pressure", "yield_mbar_ug", "residual_gas", "nchops", "min_ref_beam", "max_ref_beam",
                           "total_blocks", "total_cycles",
                           "initial_sample_beam", "pre_balance_sample_beam", "nu_balance",
                           "sam_max44nA", "sam_mean44nA", "sam_min44nA",
                           "min_balance", "mean_balance", "max_balance", "sam_decay", "wg_decay")



# -------------------- Organize all Raw Data Files -------------------- 

# Make a list of all files in the current session directory
if (dir.exists(paste(project_dir, results_dir, session$name, sep="/"))){
    # get list of all directories that may contain data
    all_dirs <- list.dirs(paste(project_dir, results_dir, session$name, sep="/"))
    if (length(all_dirs) == 0){
        all_dirs <- paste(project_dir, results_dir, session$name, sep="/")
    }
    all_dirs <- all_dirs[!grepl("archive", all_dirs)]
    # get list of all "Data_*.txt" files
    all_session_data_paths <- c()
    for (i in all_dirs){
        # list files with "Data_", at least one digit, any number of any character, any number of word characters, and .txt
        all_session_data_paths <- append(all_session_data_paths, list.files(i, pattern="Data_\\d+.*\\w*.txt", full.names=TRUE))
    }
    message("    The current session - ", session$name, " - has ", length(all_session_data_paths), " total analyses so far.\n")
} else {
    stop("    Path to data files does not exist - ", paste(project_dir, results_dir, session$name, sep="/"))
}
all_session_basenames <- basename(all_session_data_paths)


# Keep track of all raw data files we have already read in and written to the samplelog file
processed_files_log <- paste(project_dir, results_dir, session$name, paste("psi_", session$name, "_processed_files.csv", sep=""), sep="/")
if (exists('reprocess')){
    print(all_session_basenames)
    sample_to_reprocess_flag <- FALSE
    while (sample_to_reprocess_flag==FALSE){
        sample_to_reprocess <- readline("Enter the complete filename from the list above for the sample you wish to reprocess: ")
        if (sample_to_reprocess %in% all_session_basenames){
            sample_to_reprocess_flag <- TRUE
            unprocessed_basenames <- sample_to_reprocess
            new_session_paths <- all_session_data_paths[match(unprocessed_basenames, all_session_basenames)]
            comment <- paste(comment, 'reprocessed; ')
            message("    REPROCESSING file ", sample_to_reprocess)
        } else {
            print('ooops. Filename not recognized. Try again. Make sure you are copying the complete name without the quotes.')
        }
    }
    rm('reprocess')
} else if (exists('reprocess_all')) {
    if (file.exists(processed_files_log)){
        file.remove(processed_files_log)
    }
    new_session_paths <- all_session_data_paths
    file.rename(session$samplelog, paste(project_dir, results_dir, session$name, "archive", paste("psi_", session$name, "_samplelog_archive_", round(as.numeric(as.POSIXct(now()))), ".csv", sep=""), sep="/"))
} else {
    if (file.exists(processed_files_log)){
        processed_paths <- readLines(processed_files_log)
        processed_basenames <- basename(processed_paths)
        unprocessed_basenames <- subset(all_session_basenames, !(all_session_basenames %in% processed_basenames))
        new_session_paths <- all_session_data_paths[match(unprocessed_basenames, all_session_basenames)]
    } else {
        new_session_paths <- all_session_data_paths
    }
}



# -------------------- Process Raw Data Files -------------------- 
if (length(new_session_paths)>0){

    for (data_file in new_session_paths){

        if (grepl('.fail', data_file)){
            message("        File ", data_file, " failed ***** ")
        } else {

            tryCatch({
                # read in data as data frame
                raw_data <- read.table(data_file, header=FALSE, fill=TRUE, sep="", col.names=paste('col', 1:21, sep=""))

                if (dim(raw_data)[1] > 0){
                    # get run meta data
                    if (substr(raw_data[1,1],1,7)=='Clumped'){
                        analysis_type <- 'clumped'    
                    } else {
                        analysis_type <- 'bulk'
                    }
                    batch_name <- substr(data_file, regexpr(session$name,data_file)[1]+nchar(session$name), regexpr('Data_',data_file)[1]-2)
                    batch_name <- gsub("/", "", batch_name)
                    CF_status <- gsub(",|\"|\"", "", raw_data[which(raw_data == "Cold finger status", arr.ind = TRUE)[1],2])
                    pre_sn <- substr(data_file, regexpr('Data_',data_file)[1]+5, nchar(data_file))
                    UID <- substr(pre_sn, 1, regexpr(' ', pre_sn)[1]-1) # Unique IDentifier (UID)
                    sample_name <- gsub(",\"|\"","",raw_data[which(raw_data == "Sample Name is ", arr.ind = TRUE)[1],2])
                    sample_mass <- gsub(",","",raw_data[which(raw_data == "Sample Weight: ", arr.ind = TRUE)[1],2])
                    sample_pressure <- gsub(",","",raw_data[which(raw_data == "Inlet transducer pressure", arr.ind = TRUE)[1],2])
                    nu_software_version <- gsub("\"","",raw_data[which(raw_data == "Software Version", arr.ind = TRUE)[1],3])
                    nchops <- gsub(",","",raw_data[which(raw_data == "Number of sample chops", arr.ind = TRUE)[1],2])
                    residual_gas <- gsub(",|\"|\"","",raw_data[which(raw_data == "Max sample pump-over pressure (mBar)", arr.ind = TRUE)[1],2])
                    # analysis_string_time <- strptime(raw_data[2,1], format="Started analysis at %I:%M:%S %p on the %A, %B %d, %Y") # this is good for all data before 2025-01-10 - it used to look like this "Started analysis at 9:52:44 PM on the Wednesday, January 8, 2025"
                    analysis_string_time <- strptime(raw_data[2,1], format="Started analysis at %H:%M:%S on the %d %B %Y") # after 2025-01-10 the date format looks like this # "Started analysis at 10:13:26 on the 10 January 2025"
                    analysis_unix_time <- as.numeric(analysis_string_time)
                    initial_sample_beam <- gsub(",|\"","",raw_data[which(raw_data == "Initial Sam Beam", arr.ind = TRUE)[1],2])
                    pre_balance_sample_beam <- gsub(",|\"","",raw_data[which(raw_data == "Pre-balance Sam Beam", arr.ind = TRUE)[1],2])
                    nu_balance <- gsub(",|\"","",raw_data[which(raw_data == "Balance %", arr.ind = TRUE)[1],2])
                    min_ref_beam <- gsub(",|\"","",raw_data[which(raw_data == "Min ref beam", arr.ind = TRUE)[1],2])
                    max_ref_beam <- gsub(",|\"","",raw_data[which(raw_data == "Max ref beam", arr.ind = TRUE)[1],2])


                    # Identify the analysis as a standard or sample and normalize nomenclature
                    found_sample <- FALSE
                    for (carb in refmats$carbonates){
                        if (tolower(sample_name) %in% tolower(carb$names)){
                            if (grepl("sand", carb$material)){
                                # leave sample_name alone, standardize manually for now
                                sample_type <- "PercentCarbStd"
                            } else {
                                sample_name <- carb$names[1]
                                sample_type <- "CarbStd"
                            }
                            found_sample <- TRUE
                        }
                    }
                    for (gas in refmats$gases$CDES){
                        if (tolower(sample_name) %in% tolower(gas$names)) {
                            sample_name <- gas$names[1]
                            sample_type <- "RefGas"
                            found_sample <- TRUE
                        }
                    }
                    if (found_sample == FALSE) {
                        sample_type <- "sample"
                        found_sample <- TRUE
                    }


                    # Did the sample produce any data at all?
                    top_of_blocks <- which(raw_data == "Individual", arr.ind = TRUE)
                    message("\n")
                    if (length(top_of_blocks) == 0) {
                        message("        Sample ", sample_name, " (UID ", UID, ") did not produce any data")

                    } else {
                        # get number of cycles and number of blocks
                        integration <- as.numeric(gsub(",","",raw_data[which(raw_data == "Cycle_Length", arr.ind = TRUE)[1],2]))  
                        ncycles <- as.numeric(gsub(",","",raw_data[which(raw_data == "No_C_O_Cycles", arr.ind = TRUE)[1],2]))
                        nblocks <- as.numeric(gsub(",","",raw_data[which(raw_data == "Num Blocks", arr.ind = TRUE)[1],2]))
                        total_cycles <- ncycles * nblocks
                        l <- integration*total_cycles
                        if (CF_status == 'Bel') {
                            nchops <- -1
                        }

                        message("        Sample ", sample_name, " (UID ", UID, ") has ", total_cycles, " cycles in ", nblocks, " blocks")

                        if (analysis_type=='clumped'){
                            wg_cycles_49 <- vector(mode="numeric", length=total_cycles)
                            wg_cycles_48 <- vector(mode="numeric", length=total_cycles)
                            wg_cycles_47 <- vector(mode="numeric", length=total_cycles)
                            wg_cycles_49sd <- vector(mode="numeric", length=total_cycles)
                            wg_cycles_48sd <- vector(mode="numeric", length=total_cycles)
                            wg_cycles_47sd <- vector(mode="numeric", length=total_cycles)
                            sam_cycles_49 <- vector(mode="numeric", length=total_cycles)
                            sam_cycles_48 <- vector(mode="numeric", length=total_cycles)
                            sam_cycles_47 <- vector(mode="numeric", length=total_cycles)
                            sam_cycles_49sd <- vector(mode="numeric", length=total_cycles)
                            sam_cycles_48sd <- vector(mode="numeric", length=total_cycles)
                            sam_cycles_47sd <- vector(mode="numeric", length=total_cycles)
                        }

                        wg_cycles_46 <- vector(mode="numeric", length=total_cycles)
                        wg_cycles_45 <- vector(mode="numeric", length=total_cycles)
                        wg_cycles_44 <- vector(mode="numeric", length=total_cycles)
                        wg_cycles_46sd <- vector(mode="numeric", length=total_cycles)
                        wg_cycles_45sd <- vector(mode="numeric", length=total_cycles)
                        wg_cycles_44sd <- vector(mode="numeric", length=total_cycles)
                        sam_cycles_46 <- vector(mode="numeric", length=total_cycles)
                        sam_cycles_45 <- vector(mode="numeric", length=total_cycles)
                        sam_cycles_44 <- vector(mode="numeric", length=total_cycles)
                        sam_cycles_46sd <- vector(mode="numeric", length=total_cycles)
                        sam_cycles_45sd <- vector(mode="numeric", length=total_cycles)
                        sam_cycles_44sd <- vector(mode="numeric", length=total_cycles)

                        offset <- 6  # this is the offset in rows from the top of the current block
                        cycle_num_offset <- 0  # this is the offset in cycles from the first cycle -- should be renamed to cycle_counter
                        number_of_columns <- ncol(raw_data)
                        raw_i <- 1


                        for (block_num in 1:dim(top_of_blocks)[1]){
                            if (block_num==1 && raw_data[top_of_blocks[block_num]+5,2]=="Bla"){
                                if (analysis_type == 'clumped'){
                                    bkgnd_49 <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+0, (2:number_of_columns)]))
                                    bkgnd_49sd <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+0, (2:number_of_columns)]))
                                    bkgnd_48 <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)]))
                                    bkgnd_48sd <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)]))
                                    bkgnd_47 <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)]))
                                    bkgnd_47sd <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)]))
                                    bkgnd_46 <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+3, (2:number_of_columns)]))
                                    bkgnd_46sd <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+3, (2:number_of_columns)]))
                                    bkgnd_45 <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+4, (2:number_of_columns)]))
                                    bkgnd_45sd <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+4, (2:number_of_columns)]))
                                    bkgnd_44 <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+5, (2:number_of_columns)]))
                                    bkgnd_44sd <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+5, (2:number_of_columns)]))    
                                    offset <- 12
                                } else {
                                    bkgnd_46 <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+0, (2:number_of_columns)]))
                                    bkgnd_46sd <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+0, (2:number_of_columns)]))
                                    bkgnd_45 <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)]))
                                    bkgnd_45sd <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)]))
                                    bkgnd_44 <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)]))
                                    bkgnd_44sd <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)]))    
                                    offset <- 9
                                }
                            }

                            # summarize cycles
                            for (i in 1:ncycles){
                                cycle_num <- i + cycle_num_offset

                                # average adjacent, in time, reference values - this is an effort to predict what the ref signal would be if it were measured exactly when the sample signal is measured
                                if (analysis_type == 'clumped'){
                                    wg_cycles_49[cycle_num] <- mean(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1+14, (2:number_of_columns)]))))
                                    wg_cycles_49sd[cycle_num] <- sd(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1+14, (2:number_of_columns)]))))
                                    wg_cycles_48[cycle_num] <- mean(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2+14, (2:number_of_columns)]))))
                                    wg_cycles_48sd[cycle_num] <- sd(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2+14, (2:number_of_columns)]))))
                                    wg_cycles_47[cycle_num] <- mean(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+3, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+3+14, (2:number_of_columns)]))))
                                    wg_cycles_47sd[cycle_num] <- sd(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+3, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+3+14, (2:number_of_columns)]))))
                                    wg_cycles_46[cycle_num] <- mean(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+4, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+4+14, (2:number_of_columns)]))))
                                    wg_cycles_46sd[cycle_num] <- sd(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+4, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+4+14, (2:number_of_columns)]))))
                                    wg_cycles_45[cycle_num] <- mean(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+5, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+5+14, (2:number_of_columns)]))))
                                    wg_cycles_45sd[cycle_num] <- sd(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+5, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+5+14, (2:number_of_columns)]))))
                                    wg_cycles_44[cycle_num] <- mean(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+6, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+6+14, (2:number_of_columns)]))))
                                    wg_cycles_44sd[cycle_num] <- sd(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+6, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+6+14, (2:number_of_columns)]))))
                                    offset <- offset + 7
                                } else {
                                    wg_cycles_46[cycle_num] <- mean(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1+8, (2:number_of_columns)]))))
                                    wg_cycles_46sd[cycle_num] <- sd(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1+8, (2:number_of_columns)]))))
                                    wg_cycles_45[cycle_num] <- mean(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2+8, (2:number_of_columns)]))))
                                    wg_cycles_45sd[cycle_num] <- sd(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2+8, (2:number_of_columns)]))))
                                    wg_cycles_44[cycle_num] <- mean(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+3, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+3+8, (2:number_of_columns)]))))
                                    wg_cycles_44sd[cycle_num] <- sd(c(mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+3, (2:number_of_columns)])), mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+3+8, (2:number_of_columns)]))))
                                    offset <- offset + 4
                                }

                                if (analysis_type == 'clumped'){
                                    sam_cycles_49[cycle_num] <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)]))
                                    sam_cycles_49sd[cycle_num] <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)]))
                                    sam_cycles_48[cycle_num] <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)]))
                                    sam_cycles_48sd[cycle_num] <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)]))
                                    sam_cycles_47[cycle_num] <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+3, (2:number_of_columns)]))
                                    sam_cycles_47sd[cycle_num] <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+3, (2:number_of_columns)]))
                                    sam_cycles_46[cycle_num] <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+4, (2:number_of_columns)]))
                                    sam_cycles_46sd[cycle_num] <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+4, (2:number_of_columns)]))
                                    sam_cycles_45[cycle_num] <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+5, (2:number_of_columns)]))
                                    sam_cycles_45sd[cycle_num] <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+5, (2:number_of_columns)]))
                                    sam_cycles_44[cycle_num] <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+6, (2:number_of_columns)]))
                                    sam_cycles_44sd[cycle_num] <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+6, (2:number_of_columns)]))

                                    offset <- offset + 7
                                } else {
                                    sam_cycles_46[cycle_num] <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)]))
                                    sam_cycles_46sd[cycle_num] <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+1, (2:number_of_columns)]))
                                    sam_cycles_45[cycle_num] <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)]))
                                    sam_cycles_45sd[cycle_num] <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+2, (2:number_of_columns)]))
                                    sam_cycles_44[cycle_num] <- mean(as.numeric(raw_data[top_of_blocks[block_num]+offset+3, (2:number_of_columns)]))
                                    sam_cycles_44sd[cycle_num] <- sd(as.numeric(raw_data[top_of_blocks[block_num]+offset+3, (2:number_of_columns)]))

                                    offset <- offset + 4                            
                                }

                                raw_i <- raw_i + ncycles
                            }

                            if (analysis_type == 'clumped'){
                                offset <- 5
                            } else {
                                offset <- 2
                            }
                            cycle_num_offset <- cycle_num_offset + ncycles
                        
                        }

                        # signal related statistics
                        wg_cycles_44_diff <- diff(wg_cycles_44)
                        sam_cycles_44_diff <- diff(sam_cycles_44)
                        wg_initial_decay <- mean(wg_cycles_44_diff[1:5])
                        sam_initial_decay <- mean(sam_cycles_44_diff[1:5])
                        balance <- (sam_cycles_44 - wg_cycles_44) / sam_cycles_44 * 100


                        # Begin ratio and delta calculations

                        R45wg_raw <- wg_cycles_45/wg_cycles_44
                        R45sam_raw <- sam_cycles_45/sam_cycles_44

                        R46wg_raw <- wg_cycles_46/wg_cycles_44
                        R46sam_raw <- sam_cycles_46/sam_cycles_44

                        if (analysis_type == 'clumped'){
                            R47wg_raw <- wg_cycles_47/wg_cycles_44
                            R47sam_raw <- sam_cycles_47/sam_cycles_44

                            R48wg_raw<- wg_cycles_48/wg_cycles_44
                            R48sam_raw <- sam_cycles_48/sam_cycles_44

                            R49wg_raw <- wg_cycles_49/wg_cycles_44
                            R49sam_raw <- sam_cycles_49/sam_cycles_44
                        }

                        # session specific corrections, so far only session_230629
                        if (session$name=='session_230629'){

                            source(paste(project_dir, results_dir, session$name, 'psi_session_230629_CustomCorrection.R', sep="/"))

                        } else {
                            # Begin ratio and delta calculations
                            R45wg <- summarize(R45wg_raw, NULL)
                            R45sam <- summarize(R45sam_raw, NULL)

                            R46wg <- summarize(R46wg_raw, NULL)
                            R46sam <- summarize(R46sam_raw, NULL)

                            if (analysis_type == 'clumped'){
                                R47wg <- summarize(R47wg_raw, NULL)
                                R47sam <- summarize(R47sam_raw, NULL)

                                R48wg <- summarize(R48wg_raw, NULL)
                                R48sam <- summarize(R48sam_raw, NULL)

                                R49wg <- summarize(wg_cycles_49/wg_cycles_44, NULL)
                                R49sam <- summarize(sam_cycles_49/sam_cycles_44, NULL)
                            }
                        }

                        # Sample ratios normalized to stochastic working gas            
                        R45sam_wg <- R45sam$data / R45wg$data * wg_R45sto
                        R46sam_wg <- R46sam$data / R46wg$data * wg_R46sto

                        if (analysis_type == 'clumped'){
                            R47sam_wg <- R47sam$data / R47wg$data * wg_R47sto
                            R48sam_wg <- R48sam$data / R48wg$data * wg_R48sto
                            R49sam_wg <- R49sam$data / R49wg$data * wg_R49sto
                        }

                        # Calculate atomic ratios
                        R13sam <- vector(mode="numeric", length=total_cycles)
                        R17sam <- vector(mode="numeric", length=total_cycles)
                        R18sam <- vector(mode="numeric", length=total_cycles)
                        for (i in 1:total_cycles){
                          ratios <- calculate_ratios(R45sam_wg[i], R46sam_wg[i])
                          R13sam[i] <- ratios[1]
                          R17sam[i] <- ratios[2]
                          R18sam[i] <- ratios[3] 
                        }
                        
                        # stochastic sample gas
                        R45sam_sto <- R13sam + 2*R17sam
                        R46sam_sto <- 2*R18sam + 2*R13sam*R17sam + R17sam^2

                        if (analysis_type == 'clumped'){
                            R47sam_sto <- 2*R13sam*R18sam + 2*R17sam*R18sam + R13sam*R17sam^2
                            R48sam_sto <- R18sam^2 + 2*R13sam*R17sam*R18sam
                            R49sam_sto <- R13sam*R18sam^2
                        }

                        
                        # conduct an outlier test along the way 
                        nsigma <- 3.3

                        # molecular deltas
                        d45 <- summarize((R45sam$data / R45wg$data -1)*1000, nsigma)
                        d46 <- summarize((R46sam$data / R46wg$data -1)*1000, nsigma)

                        if (analysis_type == 'clumped'){
                            d47 <- summarize((R47sam$data / R47wg$data -1)*1000, nsigma)
                            d48 <- summarize((R48sam$data / R48wg$data -1)*1000, nsigma)
                            d49 <- summarize((R49sam$data / R49wg$data -1)*1000, nsigma)
                        
                            # raw capital Deltas
                            D45 <- summarize(((R45sam_wg/R45sam_sto)-1)*1000, nsigma)
                            D46 <- summarize(((R46sam_wg/R46sam_sto)-1)*1000, nsigma)
                            D47 <- summarize(((R47sam_wg/R47sam_sto)-1)*1000, nsigma)
                            D48 <- summarize(((R48sam_wg/R48sam_sto)-1)*1000, nsigma)
                            D49 <- summarize(((R49sam_wg/R49sam_sto)-1)*1000, nsigma)
                        }

                        # elemental deltas
                        d13C <- summarize((R13sam/vpdb_R13-1)*1000, nsigma)
                        d18O <- summarize((R18sam/vpdb_R18-1)*1000, nsigma)


                        # Thresholds, Flagging, and Comments on each analysis
                        trust <- 1
                        balance_flag <- FALSE
                        residual_gas_flag <- FALSE
                        d13Csd_flag <- FALSE
                        d18Osd_flag <- FALSE
                        outlier_flag <- FALSE

                        if (analysis_type == 'clumped'){
                            D47sd_flag <- FALSE
                            D48sd_flag <- FALSE
                        }

                        if (exists('reprocess')){
                            comment <- "reprocessed; "
                        } else if (exists('reprocess_all')) {
                            comment <- 'reprocess_all; '
                        } else {
                            comment <- ""
                        }

                        if (as.numeric(sample_mass)<10) {
                            if (sample_type != 'RefGas') {
                                sample_mass <- as.numeric(sample_mass)*1000
                                comment <- paste(comment, "sample mass originally entered in milligrams (mg) now in micrograms (ug)", sep="; ")
                            } else {
                                sample_mass <- NA
                                comment <- paste(comment, "mass set to NA for RefGas", sep="; ")
                            }
                        }
                        if (as.numeric(residual_gas) > 0.010){
                            residual_gas_flag <- TRUE
                            comment <- paste(comment, "residual gas high", sep="; ")
                        }
                        if (mean(balance[0:5]) > 1.0 || mean(balance[0:5]) < -1.0) {
                            balance_flag <- TRUE
                            comment <- paste(comment, "minimum balance is outside of threshold +/- 1.0", sep="; ")
                        }
                        if (d13C$s > 0.1) {
                            d13Csd_flag <- TRUE
                            comment <- paste(comment, "high d13C standard deviation", sep="; ")
                        }
                        if (d18O$s > 0.2) {
                            d18Osd_flag <- TRUE
                            comment <- paste(comment, "high d18O standard deviation", sep="; ")
                        }
                        if (analysis_type == 'clumped'){
                            if (D47$s > 0.2) {
                                D47sd_flag <- TRUE
                                comment <- paste(comment, "high D47 standard deviation", sep="; ")
                                trust <- 0
                            }
                            if (D48$s > 1.0) {
                                D48sd_flag <- TRUE
                                comment <- paste(comment, "high D48 standard deviation", sep="; ")
                            }
                        }
                        if (any(total_cycles>c(length(d45$pdi), length(d46$pdi), length(d47$pdi), length(d48$pdi), length(D47$pdi), length(D48$pdi)) )){
                            outlier_flag <-TRUE
                            comment <- paste(comment, "outliers were identified", sep="; ")
                        }

                        if (exists('reprocess_all')){
                            # don't plot or print summary info, just do the math and get on with it
                        } else {
                            # Figures and Summary
                            lfs <- 1.4 # label font size

                            layout(matrix(c(1, 2, 3,
                                            4, 5, 6,
                                            7, 8, 9), nrow=3, ncol=3, byrow=TRUE))
                            all_44 <- c(wg_cycles_44, sam_cycles_44)
                            y_limits <- range(all_44)
                            plot(wg_cycles_44, col="red",
                                 xlab="Cycle Number",
                                 ylab="m/z 44 signal (Amps)",
                                 ylim=y_limits,
                                 cex.lab=lfs)
                            points(sam_cycles_44, col="blue")

                            all_44_diff <- c(wg_cycles_44_diff, sam_cycles_44_diff)
                            y_limits <- range(all_44_diff)
                            plot(wg_cycles_44_diff,
                                 col="red",
                                 main = paste(UID, sample_name, sep=" - "),
                                 xlab="Cycle Number",
                                 ylab="signal decay (Amps / cycle)",
                                 ylim=y_limits,
                                 cex.lab=lfs)
                            points(sam_cycles_44_diff, col="blue")
                            
                            plot(wg_cycles_44-sam_cycles_44,
                                 xlab="Cycle Number",
                                 ylab="signal balance (Amps)",
                                 cex.lab=lfs)
                            
                            plot(d45$data,
                                 xlab="Cycle Number",
                                 ylab="d45 (permil)",
                                 cex.lab=lfs)
                            points(d45$odi, d45$data[d45$odi], pch=8, col='red')

                            if (analysis_type == 'clumped'){
                                plot(d47$data,
                                     xlab="Cycle Number",
                                     ylab="d47 (permil)",
                                     cex.lab=lfs)
                                points(d47$odi, d47$data[d47$odi], pch=8, col='red')
                                
                                plot(d48$data,
                                     xlab="Cycle Number",
                                     ylab="d48 (permil)",
                                     cex.lab=lfs)
                                points(d48$odi, d48$data[d48$odi], pch=8, col='red')
                            }

                            
                            plot(d46$data,
                                 xlab="Cycle Number",
                                 ylab="d46 (permil)",
                                 cex.lab=lfs)
                            points(d46$odi, d46$data[d46$odi], pch=8, col='red')
                            
                            if (analysis_type == 'clumped'){
                                plot(D47$data,
                                     xlab="Cycle Number",
                                     ylab="D47 (permil)",
                                     cex.lab=lfs)
                                points(D47$odi, D47$data[D47$odi], pch=8, col='red')
                                
                                plot(D48$data,
                                     xlab="Cycle Number",
                                     ylab="D48 (permil)",
                                     cex.lab=lfs)
                                points(D48$odi, D48$data[D48$odi], pch=8, col='red')
                            }
                            
                            message("            Analysis time: ", analysis_string_time)
                            message("            Mass: ", sample_mass)
                            message("            Residual Gas: ", round(as.numeric(residual_gas), 5), " mbar", if (residual_gas_flag){"    ** UPDATED !! - FLAGGED WITH HIGH RESIDUAL GAS **"})
                            message("            Chops: ", nchops)
                            message("            Signal Balance: ", if (balance_flag){"    ** FLAGGED AS MISBALANCED **"}, if (outlier_flag){"    ** OUTLIERS NOT REMOVED FROM BALANCE DATA **"})
                            message("                Starting Balance: ", round(balance[1], 3), " %")
                            message("                Mean Balance: ", round(mean(balance), 3), " %")
                            message("                Ending Balance: ", round(balance[length(balance)], 3), " %")
                            message("            Variance estimates: ")
                            message("                ",paste("d","45 sd",sep=''),": ", round(d45$s, 4), " \u2030", if (length(d45$pdi)<total_cycles){paste("    ** ", length(d45$odi)," outlier(s) removed **")})
                            message("                ",paste("d","46 sd",sep=''),": ", round(d46$s, 4), " \u2030", if (length(d46$pdi)<total_cycles){paste("    ** ", length(d46$odi)," outlier(s) removed **")})
                            if (analysis_type == 'clumped'){
                                message("                ",paste("d","47 sd",sep=''),": ", round(d47$s, 4), " \u2030", if (length(d47$pdi)<total_cycles){paste("    ** ", length(d47$odi)," outlier(s) removed **")})
                                message("                ",paste("d","48 sd",sep=''),": ", round(d48$s, 4), " \u2030", if (length(d48$pdi)<total_cycles){paste("    ** ", length(d48$odi)," outlier(s) removed **")})
                                message("                ",paste("D","47 sd",sep=''),": ", round(D47$s, 4), " \u2030", if (length(D47$pdi)<total_cycles){paste("    ** ", length(D47$odi)," outlier(s) removed **")})
                                message("                ",paste("D","48 sd",sep=''),": ", round(D48$s, 4), " \u2030", if (length(D48$pdi)<total_cycles){paste("    ** ", length(D48$odi)," outlier(s) removed **")})
                            }
                            message("            Finally, data I care about :) (versus working gas, which approximates VPDB_CO2)")
                            message("                ",paste("d","13C",sep=''),": ", round(d13C$m, 3), " \u2030", if(d13Csd_flag){"    ** HIGH STANDARD DEVIATION **"})
                            message("                ",paste("d","18O",sep=''),": ", round(d18O$m, 3), " \u2030", if(d18Osd_flag){"    ** HIGH STANDARD DEVIATION **"})
                            if (analysis_type == 'clumped'){
                                message("                ",paste("D","47",sep=''),": ", round(D47$m, 3), " \u2030", if(D47sd_flag){"    ** HIGH STANDARD DEVIATION **"})
                                message("                ",paste("D","48",sep=''),": ", round(D48$m, 3), " \u2030", if(D48sd_flag){"    ** HIGH STANDARD DEVIATION **"})
                            }
                            comment <- paste(comment, readline(prompt="            Enter comments and/or press [enter]: "), sep="; ")
                            comment <- gsub(",", ";", comment)
                        }

                        
            #           output headers and data
                        if (analysis_type == 'clumped'){

                            samplelog_headers <- samplelog_headers_clumped

                            samplelog_output <- c(UID, as.character(analysis_string_time), sample_name, sample_type, trust, comment, session$name,
                                                   d13C$m, d13C$s, d13C$n, d13C$slope, d13C$R2, d13C$P,
                                                   d18O$m, d18O$s, d18O$n, d18O$slope, d18O$R2, d18O$P,
                                                   d45$m, d45$s, d45$n, d45$slope, d45$R2, d45$P,
                                                   d46$m, d46$s, d46$n, d46$slope, d46$R2, d46$P,
                                                   d47$m, d47$s, d47$n, d47$slope, d47$R2, d47$P,
                                                   D47$m, D47$s, D47$n, D47$slope, D47$R2, D47$P,
                                                   d48$m, d48$s, d48$n, d48$slope, d48$R2, d48$P,
                                                   D48$m, D48$s, D48$n, D48$slope, D48$R2, D48$P,
                                                   d49$m, d49$s, d49$n, d49$slope, d49$R2, d49$P,
                                                   D49$m, D49$s, D49$n, D49$slope, D49$R2, D49$P,
                                                   R45sam$m, R45sam$s, R45wg$m, R45wg$s, mean(R45sam_wg), 
                                                   R46sam$m, R46sam$s, R46wg$m, R46wg$s, mean(R46sam_wg),
                                                   R47sam$m, R47sam$s, R47wg$m, R47wg$s, mean(R47sam_wg),
                                                   R48sam$m, R48sam$s, R48wg$m, R48wg$s, mean(R48sam_wg),
                                                   R49sam$m, R49sam$s, R49wg$m, R49wg$s, mean(R49sam_wg),
                                                   bkgnd_44, bkgnd_45, bkgnd_46, bkgnd_47, bkgnd_48, bkgnd_49,
                                                   analysis_unix_time, batch_name, nu_software_version, psi_R_version, 
                                                   sample_mass, sample_pressure, as.numeric(sample_pressure)/as.numeric(sample_mass), residual_gas, nchops, min_ref_beam, max_ref_beam,
                                                   nblocks, total_cycles,
                                                   initial_sample_beam, pre_balance_sample_beam, nu_balance,
                                                   max(sam_cycles_44), mean(sam_cycles_44), min(sam_cycles_44), max(sam_cycles_47), mean(sam_cycles_47), min(sam_cycles_47),
                                                   min(balance), mean(balance), max(balance), sam_initial_decay, wg_initial_decay)
                        } else {
                            samplelog_headers <- samplelog_headers_bulk
                            samplelog_output <- c(UID, as.character(analysis_string_time), sample_name, sample_type, trust, comment, session$name,
                                                   d13C$m, d13C$s, d13C$n, d13C$slope, d13C$R2, d13C$P,
                                                   d18O$m, d18O$s, d18O$n, d18O$slope, d18O$R2, d18O$P,
                                                   d45$m, d45$s, d45$n, d45$slope, d45$R2, d45$P,
                                                   d46$m, d46$s, d46$n, d46$slope, d46$R2, d46$P,
                                                   R45sam$m, R45sam$s, R45wg$m, R45wg$s, mean(R45sam_wg), 
                                                   R46sam$m, R46sam$s, R46wg$m, R46wg$s, mean(R46sam_wg),
                                                   analysis_unix_time, batch_name, nu_software_version, psi_R_version, 
                                                   sample_mass, sample_pressure, as.numeric(sample_pressure)/as.numeric(sample_mass), residual_gas, nchops, min_ref_beam, max_ref_beam,
                                                   nblocks, total_cycles,
                                                   initial_sample_beam, pre_balance_sample_beam, nu_balance,
                                                   max(sam_cycles_44), mean(sam_cycles_44), min(sam_cycles_44),
                                                   min(balance), mean(balance), max(balance), sam_initial_decay, wg_initial_decay)
                        }

                        if (file.exists(session$samplelog) == FALSE) {
                            write(samplelog_headers, session$samplelog, ncolumns=length(samplelog_headers), sep=",", append=TRUE)
                        }

                        write(samplelog_output, session$samplelog, ncolumns=length(samplelog_headers), sep=",", append=TRUE)

                    }

                    write(data_file, processed_files_log, sep="", append=TRUE)


                    if (exists('ncycles')){

                        # save all high res data into single data frame
                        if (analysis_type == 'clumped'){
                            all_cycles <- data.frame(c(1:ncycles), wg_cycles_44, sam_cycles_44, wg_cycles_45, sam_cycles_45, wg_cycles_46, sam_cycles_46,
                                                     wg_cycles_47, sam_cycles_47, wg_cycles_48, sam_cycles_48, wg_cycles_44sd, sam_cycles_44sd,
                                                     wg_cycles_45sd, sam_cycles_45sd, wg_cycles_46sd, sam_cycles_46sd, wg_cycles_47sd, sam_cycles_47sd,
                                                     wg_cycles_48sd, sam_cycles_48sd, R45wg$data, R45sam$data, R46wg$data, R46sam$data, R47wg$data,
                                                     R47sam$data, R48wg$data, R48sam$data, d45$data, d46$data, d47$data, d48$data, d13C$data, d18O$data,
                                                     D47$data, D48$data, d45$data - d45$m, d46$data - d46$m, d47$data - d47$m, d48$data - d48$m,
                                                     d13C$data - d13C$m, d18O$data - d18O$m, D47$data - D47$m, D48$data - D48$m)

                            all_cycles_column_names <- c('cycle', 'wg_cycles_44', 'sam_cycles_44', 'wg_cycles_45', 'sam_cycles_45', 'wg_cycles_46', 'sam_cycles_46',
                                                     'wg_cycles_47', 'sam_cycles_47', 'wg_cycles_48', 'sam_cycles_48', 'wg_cycles_44sd', 'sam_cycles_44sd',
                                                     'wg_cycles_45sd', 'sam_cycles_45sd', 'wg_cycles_46sd', 'sam_cycles_46sd', 'wg_cycles_47sd', 'sam_cycles_47sd',
                                                     'wg_cycles_48sd', 'sam_cycles_48sd', 'R45wg', 'R45sam', 'R46wg', 'R46sam', 'R47wg',
                                                     'R47sam', 'R48wg', 'R48sam', 'd45', 'd46', 'd47', 'd48', 'd13C', 'd18O',
                                                     'D47', 'D48', 'd45residual', 'd46residual', 'd47residual', 'd48residual',
                                                     'd13Cresidual', 'd18Oresidual', 'D47residual', 'D48residual')
                        } else {
                            all_cycles <- data.frame(c(1:ncycles), wg_cycles_44, sam_cycles_44, wg_cycles_45, sam_cycles_45, wg_cycles_46, sam_cycles_46,
                                                     wg_cycles_44sd, sam_cycles_44sd, wg_cycles_45sd, sam_cycles_45sd, wg_cycles_46sd, sam_cycles_46sd,
                                                     R45wg$data, R45sam$data, R46wg$data, R46sam$data, d45$data, d46$data, d13C$data, d18O$data,
                                                     d45$data - d45$m, d46$data - d46$m, d13C$data - d13C$m, d18O$data - d18O$m)

                            all_cycles_column_names <- c('cycle', 'wg_cycles_44', 'sam_cycles_44', 'wg_cycles_45', 'sam_cycles_45', 'wg_cycles_46', 'sam_cycles_46',
                                                     'wg_cycles_44sd', 'sam_cycles_44sd', 'wg_cycles_45sd', 'sam_cycles_45sd', 'wg_cycles_46sd', 'sam_cycles_46sd',
                                                     'R45wg', 'R45sam', 'R46wg', 'R46sam', 'd45', 'd46', 'd13C', 'd18O', 'd45residual', 'd46residual', 'd13Cresidual',
                                                     'd18Oresidual')
                        }
                        colnames(all_cycles) <- all_cycles_column_names

                        psi_all_cycles_file <- paste(project_dir, results_dir, session$name, paste("psi_", session$name, "_all_cycles.hdf5", sep=""), sep="/")
                        h5dataset <- paste('/UID',UID,sep="")
                        tryCatch({
                            h5write(all_cycles, psi_all_cycles_file, h5dataset)
                        }, error = function(e) {
                            message(paste('        Sample UID ',UID,' is already in the hdf5 file.', sep=""))
                        })
                    }

                } else {
                    message("        File ", data_file, " had no data ***** ")
                } 

            },  warning = function(w) {
                # Handle warnings
                print("** WARNING **")
                print(w)
                return(NA)  # Return NA if a warning occurs
            }, error = function(e) {
                # Handle errors
                print("** ERROR **")
                print(e)
                return(NA)  # Return NA if an error occurs
            }, finally = {
                # Code to run regardless of success or error
                # print("meh")
            }) 

        }
    }


} else {
    message("    All analyses in the ", session$name, " directory have already been processed. Go to the")
    message("        psi computer and double click on 'rclone_psi_results.bat' in the 'c:\\psi\\' directory.")
    message("        If new analyses exist, this will transfer them to the data server. Then come back and")
    message("        rerun psi.R.")
}

if (exists('reprocess')){
    rm(reprocess)
}
if (exists('reprocess_all')){
    rm(reprocess_all)
}


setwd(py_dir)
system(paste("py ", py_dir, "psi_raw_data_report.py ", session$name, sep=""))

setwd(R_dir)
message(" ")
message(" Done")

message(" ")
message(" ")
message("HOW DO I CALIBRATE?")
message("    Run the python script 'psi_calibrate.py'.")
message("        - You should see a python logo shortcut on the desktop and on the bottom menu bar.")
message("        - type 'py psi_calibrate.py' without the single quotes and hit enter")
message("        - You can also do the above from any Windows command prompt that has been set to S:/Data/python/")
