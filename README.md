# BoltzTraP2 Post-Processing Tools

Welcome to the **BoltzTraP2 Post-Processing** repository! This project provides Python scripts to streamline the post-processing of results generated by [BoltzTraP2](https://boltztrap2y.readthedocs.io/en/latest/BoltzTraP2.html).

## Table of Contents

- [Script Overview](#ScriptOverview)
- [Citations](#citations)
- [License](#license)
- [Commit Sign-off Requirement](#CommitSign-offRequirement)
- [Running the Script](#RunningtheScript)<br>
    [Command-Line Arguments](#Command-LineArguments)<br>
    [Interactive Prompts](#InteractivePrompts)<br>
- [Functionalities](#Functionalities)
    [1. Unit Conversion](#UnitConvertion)<br>
    [2. Fermi Energy Adjustment](#FermiEnergyAdjustment)<br>
    [3. Plotting Options](#PlottingOptions)<br>
      [Option 1: Seebeck Coefficient vs. (μ - E_F)](#Seebeckvsmu)<br>
      [Option 2: Seebeck Coefficient vs. Temperature](#SeebeckvsT)<br>
      [Option 3: Electrical Conductivity vs. Temperature](#ElecCondvsT)<br>
      [Option 4: Electronic Thermal Conductivity vs. Temperature](#ElecThermCondvsT)<br>
- [Output](#Output)
- [Examples](#Examples)
- [Troubleshooting](#Troubleshooting)

## Script Overview

The script provided (BTP2-extract.py) is designed to process BoltzTraP2 output files (*.trace), perform unit conversions, adjust the chemical potential by subtracting the Fermi energy, and generate various plots for data analysis.

Key Functionalities:
	•	Load and Process Data: Reads the input data file and processes it into a usable format.
	•	Unit Conversion: Automatically detects and converts the unit of Ef (Rydberg, Hartree, or eV) to eV.
	•	Fermi Energy Adjustment: Optionally subtracts the provided Fermi energy from the chemical potential.
	•	Data Conversion: Converts the Seebeck coefficient from V/K to µV/K.
	•	Plot Generation: Offers multiple plotting options to visualize different physical properties.

## Citations

1. **BoltzTraP2**
   
Georg K.H. Madsen, Jesús Carrete, Matthieu J. Verstraete, BoltzTraP2, a program for interpolating band structures and calculating semi-classical transport coefficients, Computer Physics Communications, Volume 231, 2018, Pages 140-145, ISSN 0010-4655, https://doi.org/10.1016/j.cpc.2018.05.010.

2. **This Repository**
   
@software{ChT_BTP2-PP,
  author = {Tantardini, Christian},
  doi = {10.5281/zenodo.14058295},
  month = {11},
  title = {BoltzTraP2-post-processing},
  url = {https://github.com/Christian48596/BoltzTraP2-post-processing},
  year = {2024}
}

## License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

The full text of the GPL General Public License can be found in file LICENSE.

## Commit Sign-off Requirement

This repository requires all contributors to sign off on their commits to affirm compliance with the Developer Certificate of Origin (DCO). Signing off confirms that:

The contribution is your original work or you have permission to contribute it.
You agree to the repository’s contribution terms.
To sign off on a commit, include the following line in your commit message:
```
Signed-off-by: Your Name <your.email@example.com>
```
If using Git locally, you can add the sign-off automatically with:
```
git commit --signoff
```
## Prerequisites

Before using the script, ensure that you have the following installed on your system:
- Python 3.6 or higher
- Required Python Packages:<br>
• pandas<br>
• numpy<br>
• matplotlib<br>
• seaborn<br>
• argparse<br>

You can install the required packages using pip:
```
pip install pandas numpy matplotlib seaborn argparse
```

## Installation

1. **Clone the Repository**

```
git clone https://github.com/Christian48596/BoltzTraP2-post-processing.git
cd BoltzTraP2-post-processing
```

2.	**(Optional) Create a Virtual Environment**

It’s recommended to use a virtual environment to manage dependencies.

```bash
python3 -m venv env
source env/bin/activate  # On Windows, use `env\Scripts\activate`
```

## Running the Script

**Command-Line Arguments**

The script accepts the following command-line argument:<br>
• -i, --input: (Required) Path to the input .trace file containing the data.

Example:

```
python BTP2-extract.py -i results.trace
```

**Interactive Prompts**

After executing the script with the required input file, you will be guided through a series of interactive prompts:

**Contact Information Display**:
  The script begins by displaying the contact information of Dr. Christian Tantardini.

**Fermi Energy Subtraction**:
• Prompt: Do you want to subtract the Fermi energy from μ? (yes/no)<br>
• Action: If you choose yes, you will be prompted to enter the Fermi energy value in eV.<br>

**Plotting Options**:
  You will be presented with a menu to select the type of plot you wish to generate:<br>

Select a plotting option:
1. Plot S vs (μ - E_F) for chosen T(s) and save .csv file.<br>
2. Plot S vs T for chosen (μ - E_F) value(s) and save .csv file.<br>
3. Plot σ/τ₀ vs T for chosen (μ - E_F) value(s) and save .csv file.<br>
4. Plot κₑ/τ₀ vs T for chosen (μ - E_F) value(s) and save .csv file.<br>
5. Exit.<br>
Input: Enter the number corresponding to your choice (1-5).<br>

**Additional Inputs Based on Plot Selection**:
• Option 1:<br>
• Prompt: Enter temperature(s) in K separated by commas (e.g., 300,400,500)<br>
• Options 2-4:<br>
• Prompt: Enter (μ - E_F) value(s) in eV separated by commas (e.g., 1.000,2.000,3.000)<br>
After entering the required values, the script will generate the corresponding plots and save the data used for plotting as .csv files in the plots directory.<br>

**Continuation Prompt**:
After completing an action, you will be asked whether you want to perform another action:<br>
• Prompt: Do you want to perform another action? (yes/no)<br>
• Action: If yes, the plotting menu will be displayed again. If no, the script will exit.<br>

*Exiting the Script*<br>

At any point, you can choose to exit the script by selecting the appropriate option from the menu or by declining to perform another action when prompted.<br>

## Functionalities

1. **Unit Conversion**
   
• Automatic Detection: The script reads the header of the input file to detect the unit of Ef. Supported units are Rydberg (Ry), Hartree (Ha), and electron volts (eV).<br>
• Conversion to eV: If Ef is not in eV, it is converted using the following factors:<br>
• 1 Rydberg (Ry) = 13.605698 eV<br>
• 1 Hartree (Ha) = 27.211386 eV<br>

2. **Fermi Energy Adjustment**

• Optional Subtraction: Users can choose to subtract a specified Fermi energy value from the chemical potential (μ). This is useful for aligning the chemical potential with the Fermi level in your analysis.<br>

3. **Plotting Options**

The script offers four plotting options:<br>

**Option 1: Seebeck Coefficient vs. (μ - E_F)**
	• Description: Plots the Seebeck coefficient (S) against the adjusted chemical potential (μ - E_F) for selected temperatures.<br>
	• Inputs Required:<br>
	• Temperature(s) in Kelvin (K).<br>
	• Output:<br>
	• Plot image files (.png and .pdf) saved in the plots directory.<br>
	• Corresponding data used for plotting saved as .csv files.<br>

**Option 2: Seebeck Coefficient vs. Temperature**
	• Description: Plots the Seebeck coefficient (S) against temperature (T) for selected (μ - E_F) values.<br>
	• Inputs Required:<br>
	• (μ - E_F) value(s) in electron volts (eV).<br>
	• Output:<br>
	• Plot image files (.png and .pdf) saved in the plots directory.<br>
	• Corresponding data used for plotting saved as .csv files.<br>

**Option 3: Electrical Conductivity vs. Temperature**
	• Description: Plots the electrical conductivity (σ/τ₀) against temperature (T) for selected (μ - E_F) values.<br>
	• Inputs Required:<br>
	• (μ - E_F) value(s) in electron volts (eV).<br>
	• Output:<br>
	• Plot image files (.png and .pdf) saved in the plots directory.<br>
	• Corresponding data used for plotting saved as .csv files.<br>

**Option 4: Electronic Thermal Conductivity vs. Temperature**
	• Description: Plots the electronic thermal conductivity (κₑ/τ₀) against temperature (T) for selected (μ - E_F) values.<br>
	• Inputs Required:<br>
	• (μ - E_F) value(s) in electron volts (eV).<br>
	• Output:<br>
	• Plot image files (.png and .pdf) saved in the plots directory.<br>
	• Corresponding data used for plotting saved as .csv files.<br>

## Output
  • Plots Directory: All generated plots and corresponding .csv data files are saved in the plots directory within the script’s root folder.<br>
  • Plot Files: Each plot is saved in both .png and .pdf formats for versatility in usage.<br>
  • CSV Files: The data used for each plot is saved as a .csv file, allowing for further analysis or record-keeping.<br>

## Examples

**Example 1: Basic Usage**

Process a .trace file without subtracting Fermi energy and generate a Seebeck coefficient vs. temperature plot.<br>
```
python BTP2-extract.py -i results.trace
```

*Interactive Steps*:<br>
	1.	Fermi Energy Subtraction:<br>
	•	Input: no<br>
	2.	Select Plotting Option:<br>
	•	Input: 2<br>
	3.	Enter (μ - E_F) Values:<br>
	•	Input: 1.000,2.000,3.000<br>
	4.	Save and Exit:<br>
	•	Input: no<br>

*Outcome*:<br>
	•	Generates S_vs_T.png and S_vs_T.pdf in the plots directory.
	•	Saves corresponding .csv data files.

**Example 2: Advanced Usage with Fermi Energy Subtraction**

Process a .trace file, subtract a Fermi energy of 5.0 eV, and generate multiple plots.

python BTP2-extract.py -i results.trace

*Interactive Steps*:
	1.	Fermi Energy Subtraction:
	•	Input: yes
	•	Fermi Energy Value: 5.0
	2.	Select Plotting Option:
	•	Input: 1
	3.	Enter Temperature Values:
	•	Input: 300,400,500
	4.	Choose to Perform Another Action:
	•	Input: yes
	5.	Select Another Plotting Option:
	•	Input: 3
	6.	Enter (μ - E_F) Values:
	•	Input: 1.500,2.500
	7.	Choose to Perform Another Action:
	•	Input: no

*Outcome*:
	•	Generates S_vs_mu_Ef.png, S_vs_mu_Ef.pdf, sigma_vs_T.png, and sigma_vs_T.pdf in the plots directory.
	•	Saves corresponding .csv data files for each plot.

## Troubleshooting

**Missing Dependencies**:
*Error*: ModuleNotFoundError: No module named 'pandas'
*Solution*: Install the missing package using pip install pandas.

**Invalid Input File**:
*Error*: Error loading data: [Error Details]
*Solution*: Ensure the input file path is correct and the file is formatted as expected.

**No Data Found for Selected Parameters**:
*Warning*: Warning: No data found for T = 600 K.
*Solution*: Verify that the entered temperature or (μ - E_F) values exist in the dataset.

**Unit Detection Failure**:
*Message*: Ef unit not specified in the header. Assuming 'Ef[eV]'.
*Solution*: Ensure that the input file header includes the unit of Ef. If not, the script defaults to eV.

**Plot Not Generated**:
*Message*: No plots were generated for S vs μ - E_F due to missing data.
*Solution*: Check the entered parameters and ensure they match the data within the specified tolerance.

 
