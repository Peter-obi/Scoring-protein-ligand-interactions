import MDAnalysis as mda
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from oddt.fingerprints import InteractionFingerprint
from oddt import toolkit
from oddt.fingerprints import sparse_to_dense
import csv
from rdkit import Chem
from rdkit.Chem import AllChem

def pdb_to_sdf(pdb_file, sdf_file):
    """Convert PDB file to SDF format with proper bond connectivity."""
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not read molecule from {pdb_file}")
    
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()

def calculate_interaction_similarity():
    """
    Calculate interaction fingerprint similarity across MD trajectory frames.
    
    Compares protein-ligand interactions in each frame to a reference frame
    and outputs cosine similarity scores to CSV file.
    """
    u = mda.Universe("step5_input.pdb", "traj.xtc")
    protein = u.select_atoms("protein")
    ligand = u.select_atoms("resname LIG")
    u.trajectory[0]
    reference_protein = protein.positions.copy()
    reference_ligand = ligand.positions.copy()
    
    # Generate reference interaction fingerprint
    protein.write("protein_ref.pdb")
    ligand.write("ligand_ref.pdb")
    pdb_to_sdf("ligand_ref.pdb", "ligand_ref.sdf")
    
    protein_ref = next(toolkit.readfile('pdb', 'protein_ref.pdb'))
    protein_ref.protein = True
    ligand_ref = next(toolkit.readfile('sdf', 'ligand_ref.sdf'))
    
    reference_ifp = InteractionFingerprint(ligand_ref, protein_ref, strict=True)
    reference_ifp_dense = sparse_to_dense(reference_ifp, len(protein.residues) * 8)
    
    with open('output.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Frame', 'Similarity'])
        
        for ts in u.trajectory[1:]:
            try:
                protein.write("protein.pdb")
                ligand.write("ligand.pdb")
                pdb_to_sdf("ligand.pdb", "ligand.sdf")
                
                protein_curr = next(toolkit.readfile('pdb', 'protein.pdb'))
                protein_curr.protein = True
                ligand_curr = next(toolkit.readfile('sdf', 'ligand.sdf'))
                
                current_ifp = InteractionFingerprint(ligand_curr, protein_curr, strict=True)
                current_ifp_dense = sparse_to_dense(current_ifp, len(protein.residues) * 8)
                
                similarity = cosine_similarity([reference_ifp_dense], [current_ifp_dense])[0][0]
                writer.writerow([ts.frame, similarity])
                
            except Exception as e:
                print(f"Error processing frame {ts.frame}: {e}")
                writer.writerow([ts.frame, 0.0])

if __name__ == "__main__":
    calculate_interaction_similarity()
