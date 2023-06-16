import MDAnalysis as mda
import numpy as np
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity
from oddt.fingerprints import InteractionFingerprint
from oddt import toolkit
from oddt.fingerprints import sparse_to_dense
import csv
from iodata import load_one, dump_one

def convert_pdb_to_sdf(input_file, output_file):
    """Converts the format of a molecular data file.
    
    Parameters
    ----------
    input_file : str
        Path to the input file.
    output_file : str
        Path where the output file will be saved.
    """
    molecule = load_one(input_file)
    dump_one(molecule, output_file)

# load the trajectory
u = mda.Universe("step5_input.pdb", "traj.xtc")
protein = u.select_atoms("protein")
ligand = u.select_atoms("resname LIG")

# reference frame
u.trajectory[0]
reference_protein = protein.positions.copy()
reference_ligand = ligand.positions.copy()

# calculate reference interaction fingerprint
protein.write("protein_ref.pdb")
ligand.write("ligand_ref.pdb")
convert_pdb_to_sdf("ligand_ref.pdb", "lig.sdf")
protein_ref = next(toolkit.readfile('pdb', 'protein_ref.pdb'))
protein_ref.protein = True
ligand_ref = next(toolkit.readfile('sdf', 'lig.sdf'))
reference_ifp = InteractionFingerprint(ligand_ref, protein_ref, strict=True)
reference_ifp_dense = sparse_to_dense(reference_ifp, len(protein.residues) * 8)

# prepare CSV output
with open('output.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Frame', 'Similarity'])

    # loop through trajectory
    for ts in u.trajectory[1:]:
        # save current frame to pdb
        protein.write("protein.pdb")
        ligand.write("ligand.pdb")
        convert_pdb_to_sdf("ligand.pdb", "ligand.sdf")
    
        # calculate current interaction fingerprint
        protein_curr = next(toolkit.readfile('pdb', 'protein.pdb'))
        protein_curr.protein = True
        ligand_curr = next(toolkit.readfile('sdf', 'ligand.sdf'))
        current_ifp = InteractionFingerprint(ligand_curr, protein_curr, strict=True)
        current_ifp_dense = sparse_to_dense(current_ifp, len(protein.residues) * 8)

        # calculate cosine similarity
        similarity = cosine_similarity([reference_ifp_dense], [current_ifp_dense])[0][0]

        # write result to CSV file
        writer.writerow([ts.frame, similarity])

