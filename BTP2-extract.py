import pandas as pd
import argparse

print("Dr. Christian Tantardini")
print("Email: christiantantardini@ymail.com")
print("https://scholar.google.com/citations?user=QCZdlUQAAAAJ")

# Conversion factors
rydberg_to_ev = 13.605698  # Rydberg to eV
hartree_to_ev = 27.211386  # Hartree to eV

# Conversion factor for Seebeck coefficient to µV/K
s_scale_factor = 1e6  # To display S in units of µV/K

# Set up argument parser
parser = argparse.ArgumentParser(
    description="Process a data file to extract and convert specific columns for analysis.\n\n"
                "The script performs the following:\n"
                "1. Checks the unit of Ef (either Rydberg, Hartree, or eV) and converts to eV if needed.\n"
                "2. Optionally subtracts the provided Fermi energy from Ef[eV].\n"
                "3. Provides an option to use standard extraction or choose specific columns to extract.\n\n"
                "Input File Format: The input file should be a space-separated text file containing columns:\n"
                "['Ef[unit]', 'T[K]', 'N[e/uc]', 'DOS(ef)[1/(Ha*uc)]', 'S[V/K]',\n"
                "'sigma/tau0[1/(ohm*m*s)]', 'RH[m**3/C]', 'kappae/tau0[W/(m*K*s)]',\n"
                "'cv[J/(mol*K)]', 'chi[m**3/mol]']\n\n"
                "Note: Lines beginning with '#' will be ignored as comments."
)

parser.add_argument(
    '-i', '--input', required=True, help='Path to the input file containing the data (e.g., data.txt)'
)
parser.add_argument(
    '-f', '--fermi', type=float, required=True, help='Fermi energy in eV to optionally subtract from Ef[eV]'
)

# Parse the arguments
args = parser.parse_args()

# Load the data from the input file into a pandas DataFrame
data = pd.read_csv(args.input, sep=r'\s+', comment='#', header=None,
                   names=['Ef[unit]', 'T[K]', 'N[e/uc]', 'DOS(ef)[1/(Ha*uc)]', 'S[V/K]',
                          'sigma/tau0[1/(ohm*m*s)]', 'RH[m**3/C]', 'kappae/tau0[W/(m*K*s)]',
                          'cv[J/(mol*K)]', 'chi[m**3/mol]'])

# Inspect the first few lines of the file to detect the unit of Ef
with open(args.input, 'r') as file:
    for line in file:
        if 'Ef[Ry]' in line:
            unit_conversion = 'Ef[Ry]'
            break
        elif 'Ef[Ha]' in line:
            unit_conversion = 'Ef[Ha]'
            break
        elif 'Ef[eV]' in line:
            unit_conversion = 'Ef[eV]'
            break
    else:
        raise ValueError("Ef unit not recognized. Make sure it is either 'Ef[Ry]', 'Ef[Ha]', or 'Ef[eV]'.")

# Convert Ef to eV based on the detected unit
if unit_conversion == 'Ef[Ry]':
    data['Ef[eV]'] = data['Ef[unit]'] * rydberg_to_ev
elif unit_conversion == 'Ef[Ha]':
    data['Ef[eV]'] = data['Ef[unit]'] * hartree_to_ev
elif unit_conversion == 'Ef[eV]':
    data['Ef[eV]'] = data['Ef[unit]']  # Already in eV, no conversion needed

# Ask the user if they want to subtract the Fermi energy
subtract_fermi = input(f"Do you want to subtract the Fermi energy ({args.fermi} eV) from Ef[eV]? (yes/no): ").strip().lower()
if subtract_fermi == 'yes':
    data['Ef[eV]'] -= args.fermi
    print(f"Fermi energy of {args.fermi} eV subtracted from Ef[eV].")

# Convert S from V/K to µV/K
data['S[µV/K]'] = data['S[V/K]'] * s_scale_factor

# Ask the user whether to use standard extraction or choose columns
user_choice = input("Do you want to use standard extraction (SvsT.dat, sigmavsT.dat, kappaevsT.dat) or choose columns? (Enter 'standard' or 'choose'): ").strip().lower()

if user_choice == 'standard':
    # Perform standard extraction
    svst_data = data[['T[K]', 'S[µV/K]']]
    svst_data.to_csv('SvsT.dat', index=False, sep='\t', float_format="%.6f")

    sigmavst_data = data[['T[K]', 'sigma/tau0[1/(ohm*m*s)]']]
    sigmavst_data.to_csv('sigmavsT.dat', index=False, sep='\t', float_format="%.6f")

    kappaevst_data = data[['T[K]', 'kappae/tau0[W/(m*K*s)]']]
    kappaevst_data.to_csv('kappaevsT.dat', index=False, sep='\t', float_format="%.6f")

    print("Standard extraction complete. Data saved to 'SvsT.dat', 'sigmavsT.dat', and 'kappaevsT.dat'.")

elif user_choice == 'choose':
    # List of available columns for extraction
    available_columns = [
        'Ef[eV]', 'T[K]', 'N[e/uc]', 'DOS(ef)[1/(Ha*uc)]', 'S[µV/K]',
        'sigma/tau0[1/(ohm*m*s)]', 'RH[m**3/C]', 'kappae/tau0[W/(m*K*s)]',
        'cv[J/(mol*K)]', 'chi[m**3/mol]'
    ]

    # Display available columns
    print("Available columns for extraction:")
    for idx, col in enumerate(available_columns, start=1):
        print(f"{idx}. {col}")

    # Get user input for the columns to extract
    selected_indices = input("Enter the numbers of the columns you want to extract, separated by commas (e.g., 1, 3, 5): ")
    selected_indices = [int(idx.strip()) - 1 for idx in selected_indices.split(",")]

    # Get user input for the filename to save the extracted data
    filename = input("Enter the filename to save the extracted data (e.g., output.dat): ")

    # Extract the selected columns
    selected_columns = [available_columns[idx] for idx in selected_indices]
    extracted_data = data[selected_columns]

    # Manually format the data for aligned output
    with open(filename, 'w') as file:
        # Write the headers
        file.write("  ".join(f"{col:>15}" for col in selected_columns) + "\n")
        # Write each row of data
        for row in extracted_data.itertuples(index=False):
            file.write("  ".join(f"{value:>15.6f}" for value in row) + "\n")

    print(f"Custom extraction complete. Data saved to {filename}.")

else:
    print("Invalid choice. Please run the script again and enter either 'standard' or 'choose'.")