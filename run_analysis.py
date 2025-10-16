import os
import subprocess
import argparse
import glob

def find_dcd_files(search_dir):
    """Find all DCD files in the specified directory and its subdirectories."""
    return glob.glob(os.path.join(search_dir, '**', '*.dcd'), recursive=True)

def run_ifcps_analysis(topology, trajectory, output, protein_selection, ligand_selection):
    """Run the IFCPS.py script with the given parameters."""
    command = [
        "python", "IFCPS.py",
        "--topology", topology,
        "--trajectory", trajectory,
        "--output", output,
        "--protein_selection", protein_selection,
        "--ligand_selection", ligand_selection
    ]
    print(f"Running analysis for {trajectory}...")
    subprocess.run(command, check=True)
    print(f"Analysis complete. Output saved to {output}")

def main():
    parser = argparse.ArgumentParser(description="Run IFCPS analysis on multiple DCD files.")
    parser.add_argument("--topology", required=True, help="Path to the topology file (e.g., prmtop).")
    parser.add_argument("--search_dir", default=".", help="Directory to search for DCD files.")
    parser.add_argument("--ligand_selection", default="resname LIG", help="Ligand selection string for IFCPS.py.")
    parser.add_argument("--protein_selection", default="protein", help="Protein selection string for IFCPS.py.")
    args = parser.parse_args()

    dcd_files = find_dcd_files(args.search_dir)

    if not dcd_files:
        print("No DCD files found.")
        return

    print(f"Found {len(dcd_files)} DCD file(s):")
    for dcd_file in dcd_files:
        output_file = os.path.splitext(dcd_file)[0] + "_ifcps.csv"
        run_ifcps_analysis(args.topology, dcd_file, output_file, args.protein_selection, args.ligand_selection)

if __name__ == "__main__":
    main()