import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# ----------- USER SETTINGS ------------
gro_file = "md.gro"
xtc_file = "md.xtc"
stride = 10          # Analyze every 10th frame to reduce memory use
chunk_size = 5       # Number of frames to process at once
output_fig = "interaction_bargraph.png"
ligand_resname = "UNL"
total_ns = 100.0     # Total simulation time in nanoseconds
# --------------------------------------

def is_hydrophobic(res):
    return res.name in ['ALA','VAL','LEU','ILE','MET','PHE','TRP','PRO']

def is_aromatic(res):
    return res.name in ['PHE','TYR','TRP','HIS']

def get_residue_id(res):
    return f"{res.name}{res.index+1}"

def is_cation(atom):
    return (atom.name.startswith("NZ") or atom.name.startswith("NH") or (atom.residue.name == ligand_resname and atom.element.symbol == "N"))

def is_anion(atom):
    return (atom.name.startswith("OD") or atom.name.startswith("OE") or (atom.residue.name == ligand_resname and atom.element.symbol == "O"))

for f in [gro_file, xtc_file]:
    if not os.path.isfile(f):
        print(f"ERROR: File not found: {f}")
        sys.exit(1)

print("Parsing topology...")
top = md.load_frame(xtc_file, 0, top=gro_file).topology

# Get time per frame and total frames using iterload (no .open())
print("Calculating time per frame and total frames...")
frame_times = []
total_frame_count = 0
for chunk in md.iterload(xtc_file, top=gro_file, chunk=1, stride=1):
    total_frame_count += chunk.n_frames
    frame_times.extend(list(chunk.time))
    if len(frame_times) >= 2:
        break

if len(frame_times) < 2:
    print("Could not determine time per frame. XTC may be too short.")
    sys.exit(1)
time_per_frame_ps = frame_times[1] - frame_times[0]
time_per_frame_ns = time_per_frame_ps / 1000.0

print(f"Total frames in trajectory: {total_frame_count}")
print(f"Time per frame: {time_per_frame_ns:.3f} ns")

# Select ligand and protein
ligand_atom_indices = top.select(f"resname {ligand_resname}")
protein_atom_indices = top.select("protein")
if len(ligand_atom_indices) == 0 or len(protein_atom_indices) == 0:
    print(f"ERROR: Could not find {ligand_resname} or protein atoms in topology.")
    sys.exit(1)

ligand_atoms = [atom for atom in top.atoms if atom.residue.name == ligand_resname]
protein_residues = [res for res in top.residues if res.is_protein]
protein_res_map = {get_residue_id(res): res for res in protein_residues}

interaction_types = ['hbond', 'salt_bridge', 'pi_pi', 'cation_pi', 'hydrophobic']
interaction_colors = {
    'hbond': '#1f77b4',
    'salt_bridge': '#d62728',
    'pi_pi': '#9467bd',
    'cation_pi': '#ff7f0e',
    'hydrophobic': '#2ca02c',
}
counts = {get_residue_id(res): {itype: 0 for itype in interaction_types} for res in protein_residues}

# ----------- CHUNKED ANALYSIS FOR LOW CONFIG SYSTEMS ------------
n_analyzed = 0
print("Analyzing trajectory in chunks...")
for chunk in md.iterload(xtc_file, top=gro_file, chunk=chunk_size, stride=stride):
    n_analyzed += chunk.n_frames

    hbond_tuples = md.baker_hubbard(chunk, periodic=False)
    for hbond in hbond_tuples:
        donor_idx, hydrogen_idx, acceptor_idx = hbond
        donor = top.atom(donor_idx)
        acceptor = top.atom(acceptor_idx)
        residues = [donor.residue, acceptor.residue]
        if any(res.name == ligand_resname for res in residues) and any(res.is_protein for res in residues):
            for res in residues:
                if res.is_protein:
                    res_id = get_residue_id(res)
                    counts[res_id]['hbond'] += 1

    for res in protein_residues:
        res_id = get_residue_id(res)
        cation_atoms = [a for a in res.atoms if is_cation(a)]
        anion_atoms = [a for a in res.atoms if is_anion(a)]
        if res.name in ['ASP','GLU']:
            for atom in anion_atoms:
                for lig_a in ligand_atoms:
                    if is_cation(lig_a):
                        dist = md.compute_distances(chunk, [[atom.index, lig_a.index]])
                        if np.any(dist < 0.4):
                            counts[res_id]['salt_bridge'] += 1
        if res.name in ['ARG','LYS','HIS']:
            for atom in cation_atoms:
                for lig_a in ligand_atoms:
                    if is_anion(lig_a):
                        dist = md.compute_distances(chunk, [[atom.index, lig_a.index]])
                        if np.any(dist < 0.4):
                            counts[res_id]['salt_bridge'] += 1

    for res in protein_residues:
        res_id = get_residue_id(res)
        if is_aromatic(res):
            prot_ring_atoms = [a.index for a in res.atoms if a.element.symbol == 'C']
            for lres in set([a.residue for a in ligand_atoms]):
                if is_aromatic(lres):
                    lig_ring_atoms = [a.index for a in lres.atoms if a.element.symbol == 'C']
                    pairs = np.array([[i, j] for i in prot_ring_atoms for j in lig_ring_atoms])
                    if len(pairs) > 0 and np.any(md.compute_distances(chunk, pairs) < 0.55):
                        counts[res_id]['pi_pi'] += 1
            cationic_lig_atoms = [a for a in ligand_atoms if is_cation(a)]
            for cat_a in cationic_lig_atoms:
                pairs = np.array([[cat_a.index, i] for i in prot_ring_atoms])
                if len(pairs) > 0 and np.any(md.compute_distances(chunk, pairs) < 0.6):
                    counts[res_id]['cation_pi'] += 1

    for res in protein_residues:
        res_id = get_residue_id(res)
        if is_hydrophobic(res):
            hydrophobic_atoms = [a.index for a in res.atoms if a.element.symbol == 'C']
            for lig_a in ligand_atoms:
                if lig_a.element.symbol == 'C':
                    pairs = np.array([[lig_a.index, i] for i in hydrophobic_atoms])
                    if len(pairs) > 0 and np.any(md.compute_distances(chunk, pairs) < 0.45):
                        counts[res_id]['hydrophobic'] += 1

print(f"Done. Analyzed {n_analyzed} frames.")

# ----------- REPORTING AND PLOTTING ------------
print("Generating plot...")

df = pd.DataFrame([
    {'residue': res_id, **counts[res_id]}
    for res_id in counts
])
df = df.set_index('residue')
df = df[(df.T != 0).any()]

if len(df) == 0:
    print("No interactions detected between ligand and protein residues.")
    sys.exit(0)

# --- Calculate percent of total simulation time covered by each detection ---
# Each detection is for one analyzed frame, each frame covers time_per_frame_ns * stride ns
frame_coverage_ns = n_analyzed * time_per_frame_ns * stride
scaling_factor = total_ns / frame_coverage_ns

for col in df.columns:
    # Cap at 100% and ensure no negative values
    df[col] = np.clip(100.0 * df[col] * scaling_factor / n_analyzed, 0, 100)

ax = df.plot(kind='bar', stacked=True,
             color=[interaction_colors[c] for c in interaction_types],
             figsize=(14, 6))
plt.title(f"Non-covalent interactions between {ligand_resname} and Protein (as % of 100 ns)")
plt.ylabel("Contact Frequency (% of 100 ns)")
plt.xlabel("Protein Residues")
plt.legend([c.replace('_', ' ').title() for c in interaction_types])
plt.tight_layout()
plt.savefig(output_fig, dpi=300)
plt.show()

print(f"Bar graph saved as: {output_fig}")
