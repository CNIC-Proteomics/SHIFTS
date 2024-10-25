# USAGE

Every SHIFTS module uses configuration files (INI), but some of the parameters it contains can also be specified through the command line. Run each module with the help option, -h, (e.g., python DMcalibrator.py -h) for a list of parameters that are available through the command line and their description.

## 1. SHIFTSadapter

This module adapts a Comet-PTM or Recom file so that it can be analyzed with SHIFTS. It removes extra lines before the column headers and adds a column containing the input file name. If the input file is of a different type that already does not have any extra lines before the column headers and contains an input file name column, this step can be skipped. Every other module following this one assumes that the input file’s first line contains the column names.

  + Input:
    - A tab-separated file from Comet-PTM or Recom output (with header).

  + Output:
    - A tab-separated file (without header).
    - A log file containing the header information.


## 2. DuplicateRemover

This module will remove scan duplicates, keeping only the scan candidate that meets the following criteria: highest rank, highest score (if there are several with the same rank), highest sp_score (if there are several with the same score).

  + Input:
    - A file from SHIFTSadapter output.
    - Command-line parameters:
    ```sh
      -s (--scan) = ScanID column name (case sensitive, column must contain an unique ID for each scan)
      -n (--num) = Rank column name (case sensitive)
      -x (--score) = Score column name (case sensitive) 
      -p (--spscore) = SpScore column name (case sensitive)
    ```

  + Output:
    - The same file, with duplicates removed.
    - A log file.


## 3. DMcalibrator

This module calculates the calibrated values for experimental masses and deltamasses. The dataset is filtered by a user-defined minimum score and maximum ppm, to obtain a subset of high-quality identifications that can be used for calibration. The median relative error (median of the absolute error column divided by the MZ column) is calculated to obtain the ‘alpha’ value:
alpha = (df_filtered[abs_error]/df_filtered[mzcolumn]).median()
Then, this value is used to calibrate the experimental MZ values of the entire dataset as follows:

```
cal_exp_mz = exp_mz * (1  - alpha)
```
  + Input:
    - A file from SHIFTSadapter or DuplicateRemover output.
    - A configuration file (INI). There is a default INI in the “config” folder, called “SHIFTS.ini”. The following parameters are used:
    ```sh
			score_min: Minimum score to filter by
			ppm_max: Maximum PPM error to filter by
			scorecolumn: Name of column containing score (case-sensitive)
			zcolumn: Name of column containing charge (case-sensitive)
			mzcolumn: Name of column containing experimental M/Z (case-sensitive)
			seqcolumn: Name of column containing sequence (case-sensitive)
			seqdmcolumn: Name of column containing sequence with deltamass within square brackets (case-sensitive)
			proteincolumn: Name of column containing protein IDs (case-sensitive)
			decoyprefix: Prefix used for decoy protein IDs (case-sensitive)
			calseqcolumn: Name of output column containing sequence with calibrated deltamass (case-sensitive)
			decimal_places: Number of decimal places to use in sequence+deltamass output columns
    ```
    - A mass file (INI). There is a default INI in the “config” folder, called “MassMod.ini”. The values it contains will be used to calculate theoretical mass values for calibration. There are three sections:
    ```sh
			Aminoacid masses
			Fixed modifications for each aminoacid (if there are any)
			Other masses (Hydrogen, Oxygen, proton)
    ```

  + Output:
    - The same file with additional columns for theoretical MH and MZ, calibrated MH and MZ, calibrated deltamass MH and MZ, absolute error, ppm error, sequence with calibrated deltamass.
    - A log file containing systematic error before and after calibration, alpha, StdDevMAD_ppm (the mean absolute deviation for the ppm error), number of PSMs before and after filtering (size of the filtered, high quality dataset that was used for calibration)


## 4. PeakModeller

This module concatenates a group of files, generates a histogram grouped by deltamass bins of a user-specified width, and calculates the frequency and slopes (first and second derivatives) for each bin, using linear regression.

  + Input:
    - A text file containing a list of calibrated files to be modelled together
    - A configuration file (INI). There is a default INI in the “config” folder. The following parameters are used:
    ```sh
			bins: Width of the bins
			slope_points: Number of points (bins) to use for slope calculation
			smooth_points: Number of points (bins) to use for smoothing
    ```

  + Output:
    - A tab-separated file containing all the input files together (DMTable), including a new ‘Filename’ column that tracks the origin of each PSM.
    - A tab-separated file containing the histogram (DMHistogram). This table contains columns with the bin, bin midpoint, frequency, slope1 (first derivative) and slope2 (second derivative).
    - A log file.


## 5. PeakSelector

This module filters a DMHistogram table (as generated by PeakModeller) according to user-specified thresholds for the columns containing slope and frequency. This filtered dataset is used to calculate a list of apexes. For each bin (row) where the value of the first derivative (slope) changes from positive to negative, the exact DM value where the first derivative is 0 is interpolated using a user-defined number of points surrounding the apex bin.
Taking the x (bin midpoint) and y (first derivative) values for those points it will identify the apex as the x value where y equals 0.

  + Input:
    - A DMHistogram table (as generated by PeakModeller) to be filtered.
    - A configuration file (INI). There is a default INI in the “config” folder. The following parameters are used:
    ```sh
			slope: Threshold for slope of DM peak
			frequency: Threshold for number of PSMs
			apex_points: Number of points (bins) to use for apex calculation
    ```

  + Output:
    - A text file containing the apex list.
    - A log file containing the number of apexes that were calculated.


## 6. PeakAssignator

This module will assign every PSM to the closest deltamass peak found in the provided apex list and identify it as either belonging to that peak, or as an orphan. To do so, it will calculate the absolute distance between the deltamass value and the assigned peak and then the error in ppm:
```
distance = abs(assigned_peak - delta_MH)
distance_ppm = (distance / (theoretical_mass + assigned_peak)) * 1e6
```
If that value is larger than the user-defined threshold, it will be considered an orphan.

  + Input:
    - A DMtable.
    - An apex list.
    - A configuration file (INI). There is a default INI in the “config” folder. The following parameters are used:
    ```sh
			ppm_max: Maximum ppm difference for peak assignation
			peak_label: Label for peaks
			orphan_label: Label for orphans
			caldeltamh_column: Name of column containing calibrated Delta MH
			theomh_column: Name of column containing theoretical MH for ppm error calculation
			closestpeak_column: Output column that will contain the closest peak
			peak_column: Output column that will contain the peak/orphan labels
			deltamass_column: Output column that will contain the assigned deltamass
			ppm_column: Output column that will contain the ppm error
			mod_peptide_column: Name of column containing sequence with deltamass in XXX[DM]XXX or XXXXXX_DM format (case-sensitive)
			assign_seq_column: Name of output column containing sequence with assigned deltamass (case-sensitive)
			decimal_places: Number of decimal places to use in sequence+deltamass output columns
    ```

  + Output:
    - A DMtable with additional columns for the closest peak, assignation as peak or orphan, assigned deltamass, ppm error, sequence with assigned deltamass.
    - A log file.


## 7. PeakFDRer

This module will calculate global, local, and peak FDR values for a DMtable subdivided by experiment. Global FDR is calculated per experiment, local and peak FDR is calculated for the whole dataset. Experiment groups are created using the information in the “Filename” column (created by the DMcalibrator module) of the DMTable, and a user-provided list that assigns each filename to a group. For the global FDR, rather than taking an entire experiment, this module will separate it in two deltamass regions (defined by the parameter dm_region_limit, default value -56) and calculate a global FDR for each region.

  + Input:
    - A DMtable.
    - A tab-separated file containing, in order:
    ```sh
			A column with the batch name
			A column with the experiment name
			A column with the file path (must match the file paths in the Filename column of your DMtable).
    ```
    - A configuration file (INI). There is a default INI in the “config” folder. The following parameters are used:
    ```sh
			score_column: Name of column containing score (case-sensitive)
			dm_region_limit: Deltamass region limit for Global FDR. Two regions will be created, DM equal to or above and DM below this value
			dm_column: Name of column containing deltamass for region limits (case-sensitive)
			peak_outlier_value: Peak FDR value to be assigned to orphans (default=1)
			peak_label: Label for peaks
			peak_column: Output column that will contain the peak/orphan labels
			caldeltamh_column: Name of column containing calibrated Delta MH
			closestpeak_column: Output column that will contain the closest peak
    ```
    - An output path.

  + Output:
    - One DMtable for each batch with additional columns for the global, local and peak FDR ranks and values.
    - A log file.

