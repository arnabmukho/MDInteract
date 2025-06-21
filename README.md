**MDInteract** is a Python tool for analyzing non-covalent interactions between a protein and a ligand in molecular dynamics (MD) simulations. It parses GROMACS trajectory files and generates a frequency bar graph of contacts (hydrogen bonds, salt bridges, π-π stacking, cation-π, hydrophobic) between your ligand and each protein residue.

## Features

- **Efficient analysis** of large trajectories using chunked and strided processing.
- **Detection of five interaction types:** hydrogen bonds, salt bridges, π-π stacking, cation-π, and hydrophobic contacts.
- **Customizable settings** for file names, stride, chunk size, ligand name, and simulation time.
- **Bar graph output** showing interaction frequency as a percentage of total simulation time, ready for publication or reports.
- **Automatic error-checks** for missing files or mismatches in selection.

## Requirements

- Python 3.x
- [MDTraj](http://mdtraj.org/)
- NumPy
- Pandas
- Matplotlib

Install these dependencies with:
```bash
pip install mdtraj numpy pandas matplotlib
```

## Usage

1. **Place Files in Directory:**  
   Ensure your GROMACS `.gro` (topology) and `.xtc` (trajectory) files are accessible.

2. **Edit User Settings in the Script:**  
   At the top of the script, adjust parameters as needed:
   ```python
   gro_file = "md.gro"          # Topology file (GROMACS .gro)
   xtc_file = "md.xtc"    # Trajectory file (.xtc)
   stride = 10                      # Analyze every 10th frame
   chunk_size = 5                   # Number of frames per chunk
   output_fig = "interaction_bargraph.png"  # Output plot filename
   ligand_resname = "UNL"           # Residue name of ligand
   total_ns = 100.0                 # Total simulation time for normalization (ns)
   ```

3. **Run the Tool:**
   ```bash
   python MDInteract.py
   ```

4. **View Output:**  
   The tool generates a bar graph (`interaction_bargraph.png`) summarizing the interaction frequency (as % of total simulation time) for each protein residue and interaction type. The plot is also displayed interactively.

## How it Works

- **Selection:** Finds ligand atoms by residue name, and all protein atoms.
- **Interaction Detection:**  
  - *Hydrogen bonds* via MDTraj's `baker_hubbard` algorithm.
  - *Salt bridges* by distance checks between charged protein and ligand groups.
  - *π-π stacking* and *cation-π* by aromatic and cationic group proximity.
  - *Hydrophobic contacts* by close carbon-carbon contacts.
- **Chunked Processing:** Processes the trajectory in configurable chunks and strides for memory efficiency.
- **Normalization:** Frequencies are scaled to represent their occurrence over the total simulation time.

## Customization

- Change `ligand_resname` to match your ligand’s residue name in your topology file.
- Adjust `stride` and `chunk_size` for a balance between speed and memory use.
- Set `total_ns` to your actual simulation time for accurate normalization.

## Troubleshooting

- **File Not Found:**  
  Check your file paths and names.
- **No Interactions Detected:**  
  Confirm that the ligand residue name matches the topology and that your trajectory actually contains interactions.

## License

GNU License

## Citation & Acknowledgements

- Built using [MDTraj](https://github.com/mdtraj/mdtraj).
- Inspired by standard protocols in computational structural biology.

---

**Author:**  
Arnab Mukherjee, PhD
arnabbiotech.gen@gmail.com
