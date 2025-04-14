# Code to parse HepMC3 data and extract particles from events, with the ability to save modified files

def parse_and_modify_hepmc3_event(file, output_file):
    """
    Parse HepMC3 data from a file object, modify statuses as specified, and save to a new file.

    Args:
        file (TextIO): File object for the HepMC3 formatted file.
        output_file (TextIO): File object for writing the modified HepMC3 data.
    """
    current_event = None

    for line in file:
        line = line.strip()
        if line.startswith("E "):  # Start of an event
            current_event = int(line.split()[1])  # Extract event number
            output_file.write(line + "\n")
        elif line.startswith("P "):  # Particle line
            # Extract particle details
            columns = line.split()
            particle_id = columns[1]  # Particle ID
            pdg_code = columns[2]  # PDG Code
            status = int(columns[8])  # Status

            # Modify status if it is 4 and the particle is not proton (2212) or electron (11)
            if status == 4 and pdg_code not in ["11", "2212"]:
                columns[8] = "999"  # Example: Change status to 999

            # Write the modified particle line to the output file
            output_file.write(" ".join(columns) + "\n")
        else:
            # Write all other lines unchanged
            output_file.write(line + "\n")

# Example usage
if __name__ == "__main__":
    import sys
    import os

    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    output_path = os.path.splitext(file_path)[0] + "_status_corrected" + os.path.splitext(file_path)[1]

    # Open the input file for reading and output file for writing
    with open(file_path, 'r') as file, open(output_path, 'w') as output_file:
        parse_and_modify_hepmc3_event(file, output_file)

    print(f"Modified file saved as: {output_path}")

