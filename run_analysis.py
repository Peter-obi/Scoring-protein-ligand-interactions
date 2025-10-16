import os
import argparse
import glob
import MDAnalysis as mda
from IFCPS import calculate_similarity_for_trajectory

def find_dcd_files(search_dir):
    """Find all DCD files in the specified directory and its subdirectories."""
    return glob.glob(os.path.join(search_dir, '**', '*.dcd'), recursive=True)

def get_protein_pocket_selection(u, ligand_selection, distance):
    """
    Get the protein selection string for the pocket around the ligand in the first frame.
    """
    u.trajectory[0]
    ligand_atoms = u.select_atoms(ligand_selection)
    protein_atoms_near_ligand = u.select_atoms(f"(protein) and (around {distance} group ligand_atoms)", updating=False, ligand_atoms=ligand_atoms)

    if len(protein_atoms_near_ligand) == 0:
        raise ValueError("No protein atoms found within the specified distance of the ligand.")

    resids = ' '.join(map(str, sorted(set(protein_atoms_near_ligand.resids))))
    return f"protein and (resid {resids})"

def main():
    parser = argparse.ArgumentParser(description="Run IFCPS analysis on multiple DCD files.")
    parser.add_argument("--topology", required=True, help="Path to the topology file (e.g., prmtop).")
    parser.add_argument("--search_dir", default=".", help="Directory to search for DCD files.")
    parser.add_argument("--ligand_selection", default="resname LIG", help="Ligand selection string.")
    parser.add_argument("--distance", type=float, default=6.0, help="Distance for protein pocket selection (in Angstroms).")
    args = parser.parse_args()

    dcd_files = find_dcd_files(args.search_dir)

    if not dcd_files:
        print("No DCD files found.")
        return

    print(f"Found {len(dcd_files)} DCD file(s):")
    for dcd_file in dcd_files:
        print(f"Processing {dcd_file}...")
        u = mda.Universe(args.topology, dcd_file)

        try:
            protein_pocket_selection = get_protein_pocket_selection(u, args.ligand_selection, args.distance)

            protein = u.select_atoms(protein_pocket_selection)
            ligand = u.select_atoms(args.ligand_selection)

            output_file = os.path.splitext(dcd_file)[0] + "_ifcps.csv"

            calculate_similarity_for_trajectory(u, protein, ligand, output_file)
            print(f"Analysis complete. Output saved to {output_file}")
        except ValueError as e:
            print(f"Error processing {dcd_file}: {e}")


if __name__ == "__main__":
    main()