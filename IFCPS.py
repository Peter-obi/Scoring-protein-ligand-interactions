import MDAnalysis as mda
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from oddt.fingerprints import InteractionFingerprint
from oddt import toolkit
from oddt.fingerprints import sparse_to_dense
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
import argparse
import tempfile
import os

def pdb_to_sdf(pdb_file, sdf_file):
    """Convert PDB file to SDF format with proper bond connectivity."""
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not read molecule from {pdb_file}")
    
    with Chem.SDWriter(sdf_file) as writer:
        writer.write(mol)

def get_interaction_fingerprint(protein_pdb, ligand_pdb, ligand_sdf, protein_residues):
    """Generate the interaction fingerprint for a protein-ligand pair."""
    pdb_to_sdf(ligand_pdb, ligand_sdf)
    
    protein = next(toolkit.readfile('pdb', protein_pdb))
    protein.protein = True
    ligand = next(toolkit.readfile('sdf', ligand_sdf))
    
    ifp = InteractionFingerprint(ligand, protein, strict=True)
    return sparse_to_dense(ifp, len(protein_residues) * 8)

def process_trajectory(u, protein, ligand, output):
    """Process the trajectory to calculate and write similarity scores."""
    with tempfile.TemporaryDirectory() as temp_dir:
        protein_ref_pdb = os.path.join(temp_dir, "protein_ref.pdb")
        ligand_ref_pdb = os.path.join(temp_dir, "ligand_ref.pdb")
        ligand_ref_sdf = os.path.join(temp_dir, "ligand_ref.sdf")

        # Generate reference interaction fingerprint from the first frame
        u.trajectory[0]
        protein.write(protein_ref_pdb)
        ligand.write(ligand_ref_pdb)
        reference_ifp_dense = get_interaction_fingerprint(protein_ref_pdb, ligand_ref_pdb, ligand_ref_sdf, protein.residues)

        with open(output, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Frame', 'Similarity'])
            writer.writerow([u.trajectory[0].frame, 1.0])

            protein_pdb = os.path.join(temp_dir, "protein.pdb")
            ligand_pdb = os.path.join(temp_dir, "ligand.pdb")
            ligand_sdf = os.path.join(temp_dir, "ligand.sdf")

            for ts in u.trajectory[1:]:
                try:
                    protein.write(protein_pdb)
                    ligand.write(ligand_pdb)
                    current_ifp_dense = get_interaction_fingerprint(protein_pdb, ligand_pdb, ligand_sdf, protein.residues)

                    similarity = cosine_similarity([reference_ifp_dense], [current_ifp_dense])[0][0]
                    writer.writerow([ts.frame, similarity])
                except Exception as e:
                    print(f"Error processing frame {ts.frame}: {e}")
                    writer.writerow([ts.frame, 0.0])

def main():
    parser = argparse.ArgumentParser(description="Calculate interaction fingerprint similarity for a MD trajectory.")
    parser.add_argument("--topology", default="step5_input.pdb", help="Topology file (e.g., PDB)")
    parser.add_argument("--trajectory", default="traj.xtc", help="Trajectory file (e.g., XTC)")
    parser.add_argument("--output", default="output.csv", help="Output CSV file")
    parser.add_argument("--protein_selection", default="protein", help="Protein selection string")
    parser.add_argument("--ligand_selection", default="resname LIG", help="Ligand selection string")
    args = parser.parse_args()

    u = mda.Universe(args.topology, args.trajectory)
    protein = u.select_atoms(args.protein_selection)
    ligand = u.select_atoms(args.ligand_selection)

    process_trajectory(u, protein, ligand, args.output)

if __name__ == "__main__":
    main()