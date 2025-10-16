import MDAnalysis as mda
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from oddt.fingerprints import InteractionFingerprint
from oddt import toolkit
from oddt.fingerprints import sparse_to_dense
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
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

def calculate_similarity_for_trajectory(u, protein, ligand, output_file):
    """
    Process the trajectory to calculate and write similarity scores.

    Args:
        u (MDAnalysis.Universe): The universe containing the trajectory.
        protein (MDAnalysis.AtomGroup): The protein atom group.
        ligand (MDAnalysis.AtomGroup): The ligand atom group.
        output_file (str): Path to the output CSV file.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        protein_ref_pdb = os.path.join(temp_dir, "protein_ref.pdb")
        ligand_ref_pdb = os.path.join(temp_dir, "ligand_ref.pdb")
        ligand_ref_sdf = os.path.join(temp_dir, "ligand_ref.sdf")

        # Generate reference interaction fingerprint from the first frame
        u.trajectory[0]
        protein.write(protein_ref_pdb)
        ligand.write(ligand_ref_pdb)
        reference_ifp_dense = get_interaction_fingerprint(protein_ref_pdb, ligand_ref_pdb, ligand_ref_sdf, protein.residues)

        with open(output_file, 'w', newline='') as f:
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