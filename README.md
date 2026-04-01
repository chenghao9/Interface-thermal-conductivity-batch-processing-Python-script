# Interface-thermal-conductivity-batch-processing-Python-script
A Python-based batch processing toolkit for calculating Interface Thermal Conductance (ITC) from LAMMPS simulation outputs. 
# ITC Batch Processing Toolkit

A Python-based batch processing toolkit for calculating **Interface Thermal Conductance (ITC)** from LAMMPS simulation outputs. This project reads temperature profile data and heat flux exchange data from multiple case directories, performs linear fitting, computes the interfacial temperature drop and average heat power, and finally reports ITC values in both **W/m²·K** and **MW/m²·K**. It also generates figures, logs, per-case summaries, and a batch summary table. 

## Features

- Batch processing of multiple simulation cases
- Automatic or manual discovery of case directories
- Extraction of averaged temperature profiles from chunked LAMMPS output
- Linear fitting of hot/cold energy exchange curves
- Linear fitting of temperature profiles on both sides of the interface
- Automatic calculation of:
  - average heat power
  - interfacial temperature difference
  - interface thermal conductance (ITC)
- Export of figures, logs, case summaries, and batch summary CSV files
- Flexible configuration through a plain-text `rc_itc.in` file with per-case overrides 

## Repository Structure

A typical working directory should look like this:

```text
project_root/
├─ deal_ITC_batch_rc.py
├─ rc_itc.in
└─ cases/
   ├─ case01/
   │  ├─ temp.CBNT.txt
   │  └─ heat_flux_exchange.txt
   ├─ case02/
   │  ├─ temp.CBNT.txt
   │  └─ heat_flux_exchange.txt
   └─ case03/
      ├─ temp.CBNT.txt
      └─ heat_flux_exchange.txt
```

In practice, the script uses the `root_dir` defined in `rc_itc.in` and processes the case folders inside that directory. 

## Requirements

- Python 3.9+
- NumPy
- Pandas
- Matplotlib
- SciPy

Install dependencies with:

```bash
pip install numpy pandas matplotlib scipy
```

## Input Files

Each case directory should contain the following input files:

- `temp.CBNT.txt` — temperature profile output
- `heat_flux_exchange.txt` — heat exchange data used for linear fitting of energy flow 

## Configuration

The program is controlled by a configuration file, for example `rc_itc.in`.

### Example configuration

```ini
# Root directory containing all case folders
root_dir = C:/Users/chenghao/OneDrive/Desktop/21

# list: only process folders listed in case_dirs
# auto: automatically process all subfolders under root_dir
scan_mode = auto

# Used only when scan_mode = list
case_dirs = case01, case02, case03

# Batch summary output name
batch_summary_name = batch_itc_results_summary.csv

[defaults]
# Input files
temp_input_name = temp.CBNT.txt
heat_flux_input_name = heat_flux_exchange.txt

# Output files
extracted_temp_name = extracted_temperature_data.txt
mean_temp_name = mean_temperature_data.txt
mean_temp_figure_name = mean_temperature_vs_position.png
heat_flux_figure_name = energy_flux_fit.png
temp_fit_figure_name = temperature_fit.png
summary_name = result_summary.csv
log_name = ITC_results_log.txt

# Calculation parameters
num_blocks = 3
num_chunks = 50
position_factor = 197.39605679
time_ps_factor = 0.0005
fit_range1 = 35, 92
fit_range2 = 112, 166
mid_x = 99
area_A2 = 261.6881

[case:case03]
fit_range1 = 36, 90
fit_range2 = 114, 168
mid_x = 100
area_A2 = 270.0
```

### Key parameters

- `root_dir`: root directory containing all case folders
- `scan_mode`: folder scanning mode
  - `auto`: process all subdirectories under `root_dir`
  - `list`: process only folders listed in `case_dirs`
- `num_blocks`: number of temperature blocks used for averaging
- `num_chunks`: number of chunks in the temperature profile
- `position_factor`: scaling factor for chunk position
- `time_ps_factor`: conversion factor from timestep to ps
- `fit_range1`, `fit_range2`: fitting ranges on the left and right sides of the interface
- `mid_x`: interface position used to evaluate the temperature jump
- `area_A2` or `area_m2`: interface area used in ITC calculation

The configuration system also supports per-case overrides using sections such as `[case:case03]`.

## Usage

Run with the default configuration file:

```bash
python deal_ITC_batch_rc.py
```

Run with a custom configuration file:

```bash
python deal_ITC_batch_rc.py rc_itc.in
```

The script will scan the target case folders, process them one by one, and write a batch summary CSV file into the configured `root_dir`. 

## Workflow

The script follows the steps below for each case:

1. Read the raw temperature file and extract chunk-wise temperature data.
2. Compute the mean temperature profile and save it.
3. Read heat exchange data and perform linear fitting for hot and cold regions.
4. Convert the fitted slopes into heat power.
5. Fit the temperature profile on both sides of the interface.
6. Evaluate the temperature discontinuity at the interface position.
7. Calculate ITC using:

```text
ITC = P_avg / (ΔT × A)
```

where:

- `P_avg` is the average heat power
- `ΔT` is the interfacial temperature drop
- `A` is the interface area 

## Output Files

For each case, the program generates:

- `extracted_temperature_data.txt` — extracted raw temperature data
- `mean_temperature_data.txt` — averaged temperature profile
- `mean_temperature_vs_position.png` — mean temperature profile plot
- `energy_flux_fit.png` — heat flux linear fitting plot
- `temperature_fit.png` — temperature fitting plot near the interface
- `result_summary.csv` — per-case summary
- `ITC_results_log.txt` — detailed calculation log

For batch execution, the root directory also contains:

- `batch_itc_results_summary.csv` — combined summary for all processed cases 

## Notes

- Make sure the input file names in each case folder match those defined in the configuration file.
- At least one of `area_A2` or `area_m2` must be provided.
- The fitting ranges should be selected carefully to avoid non-linear regions near thermostats or noisy interface regions.
- If a case fails during processing, the batch run continues and records the failure in the final summary table. 

## Citation and Attribution

- **Author**: Cheng Hao
- **Email**: chenghao8425@163.com
