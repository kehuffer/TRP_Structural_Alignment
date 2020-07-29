# Structural Alignment
Huffer, K.E., Aleksandrova, A.A., Jara-Oseguera, A., Forrest, L.R., & K.J. Swartz (2020). Global alignment and assessment of TRP channel transmembrane domain structures to explore functional mechanisms. eLife. 

## Requirements
Requires .xml file containing the following information for each structure:
- PDB ID 
- Chain order (counterclockwise from extracellular view)
- Residue range for each chain
- Subfamily (TRPV, TRPM, etc.)
- Method of structure determination (x-ray, cryo-EM)
- Imaging environment (detergent, nanodisc, amphipol)
- Ligand state (apo, activator, inhibitor)

Requires paths.txt file specifying locations of important directories and programs
- work_dir: the directory to which the structures, analysis, and output will be saved.  If not specified in paths.txt, the current working directory will be used.
- structs_info: the path to the required .xml file containing information about each structure
- frtmalign: the path to the frtmalign.exe executable
- hole: the path to the hole executable
- provided_struct: path to directory containing any pdbid.pdb files that are not found in OPM (see below)

Van der Waals radius file for HOLE, simple.rad

One reference structure for HOLE analysis, along with the xyz coordinates of a point within that structure's pore.

---

## Dependencies
[pyali, release 3.5](https://github.com/christang/pyali)

[HOLE](http://www.holeprogram.org/)

[Fr-TM-Align](http://cssb.biology.gatech.edu/skolnick/files/FrTMalign/index.html)

[MDAnalysis](https://www.mdanalysis.org/)

[SciPy](https://www.scipy.org/)

[Seaborn](https://seaborn.pydata.org/#)

[DSSP, release 3.0.0](https://github.com/cmbi/hssp)

---

## Options and Customization
Specifying which structures and domains to include in analysis
1) The struct_info191031.xml file contains TM domains of TRP channel structures up to a cutoff date of Oct. 31, 2019, and can be used to reproduce the analysis in the original paper. 
2) The struct_info200709.xml file contains TM domains of TRP channel structures up to a cutoff date of Jul. 7, 2020.
3) If you are interested in different structures or domains, you can create your own .xml file.

Options for obtaining the structure files:
Based on the PDB IDs provided in the .xml file
1) The program will first attempt to obtain the structure from the [Orientations of Proteins in Membranes Database] (https://opm.phar.umich.edu/)
2) If there is no structure in OPM, the program will then look in the provided_struct directory specified in the paths.txt file.  You may include in this directory a file named [pdbid].pdb for any structure you wish to use.  If you wish to have it embedded in the membrane like the OPM structures, you may first analyze it using the [PPM server] (https://opm.phar.umich.edu/ppm_server).
3) If there is no provided file with the pdbid in the provided_struct directory, the program will attempt to obtain the structure from the [RCSB Protein Data Bank] (https://www.rcsb.org/).
4) If there is no structure in the PDB, the program will ignore that pdbid and output a warning.

---

## Outputs
Analysis sub-folders containing each mobile structure aligned to one single stationary structure
- Fr-TM-Align output text files (filename: mobile_stationary.frtxt)
- Fr-TM-Align output rotation matrix file (filename: mobile_stationary.mat)
- Fr-TM-Align output superimposed coordinates in rasmol format (filename: mobile_stationary.sup; to view the superimposed structures by rasmol: `./rasmol - script TM.sup`)
- Fr-TM-Align parwise alignment Markx0 files (filename: mobile_stationary.aln)
- Fr-TM-Align pairwise alignment FASTA files (filename: mobile_stationary.fasta)
- mobile structure PDB file transformed using Fr-TM-Align transformation matrix so it aligns to the stationary structure (filename: mobile_stationary_full_align.pdb) (Note that even though this file includes all atoms/residues from the .pdb file, the Fr-TM-Align analysis was performed only on the regions indicated in the struct_info.xml file)
- Selected regions of mobile structure transformed using Fr-TM-Align transformation matrix so it aligns to the stationary structure (filename: mobile_stationary_tmem_align.pdb)

Analysis sub-folder containing each mobile structure aligned to hole_reference_struct specified in paths.txt
- HOLE output (filename: mobile_stationary_hole_out.txt)
- HOLE output (filename: mobile_stationary_hole_plot.pdf and mobile_stationary_hole_plot.png)
- HOLE output (filename: mobile_stationary_hole_plot.pdb)
- pore radii mapped to contributing residues (filename: mobile_stationary_radius.csv)
- dssp secondary structure assignment, FASTA format (filename: mobile_stationary_dssp.fa)
	
Parent analysis folder
- Summary of the data outputted by Fr-TM-Align for all pairs of structures (filename: frtmalign_output.csv) (contains pdbid of mobile structure, length of mobile structure, pdbid of stationary structure, length of stationary structure, length of alignment by Fr-TM-Align, RMSD of aligned structures, TM-score of aligned structures, and sequence identity of aligned structures)
- Separate files for each of the above metrics (alignment length, RMSD, TM-score, and sequence identity) for each pairwise alignment of structures.  PDB IDs in rows and columns represent mobile and stationary structures, respectively. (filenames: aligned_length.csv, RMSD.csv, TM_score.csv, and sequence_identity.csv)
- Clustermap as in Figure 2 of paper, with TM score shown for all pairs of mobile and stationary structures and clustered along the stationary axis to group most similar together.  The clustermap also contains labels for the subfamily, method, environment, and ligand states defined in struct_info.xml.  (filenames: clustermap.pdf, clustermap.png)
- Ordered list of PDB IDs after clustering along the stationary axis (filename: statinary_clustering_order.csv)
- Multiple sequence alignment FASTA files for each stationary/reference structure (filenames: stationary_full.ali and stationary_nogap.ali)
- full indicates that MSA includes all residues and gaps from all sequences
- nogap indicates that MSA omits any gaps or insertions in the reference structure, i.e. the MSA includes only residues or gaps in each sequence that align to a residue in the reference sequence
- Jalview pore radius annotation files (filenames: stationary_full_annot.jlv, stationary_nogap_annot.jlv)
- DSSP information outputted into pseudo-FASTA format that matches the multiple sequence alignment (filenames: stationary_full_dssp.fa, stationary_nogap_dssp.fa)
- For the hole_reference_struct specified in paths.txt, the MSA, radius, and dssp information is also saved in .csv files (filenames: stationary_full_msa.csv, stationary_full_radius.csv, stationary_full_dssp.csv, stationary_nogap_msa.csv, stationary_nogap_radius.csv, stationary_nogap_dssp.csv)

---

## Using and Visualizing the Outputs
Use PyMOL to open the sphpdb.pdb file containing spheres representing the pore determined with HOLE, and use the commands below to set the Van der Waals radius of each sphere equal to the pore radius value stored as the B-factor, then update the display to visualize the new sphere sizes:
```
alter sphpdb_name, vdw=b
rebuild
```

Multiple sequence alignments may be viewed with your choice of alignment program that accepts FASTA files (we used Jalview).  

To visualize the multiple sequence alignment with pore-contributing residues colored according to radius, use Jalview to open the multiple sequence alignment (pdbid_full.ali or pdbid_nogap.ali) and the corresponding annotations file (pdbid_full_annot.jlv or pdbid_nogap_annot.jlv) in the same project. 
To color the multiple sequence alignment by radius:
- Colour > By Annotation
- In the dialog box, select the tick-box next to "Per-sequence only," and select "radius" in the top dropdown menu. From there, you can customize the colors and thresholds.
- Please note that coloring based on Jalview annotations is typically normalized according to the min and max value of each individual sequence. To force all sequences to use the same color scale, we have changed the last two "radius" values of each sequence to a global min and global max, which you can see by scrolling to the end of the alignment.  These values are artifacts and do NOT represent pore contributions by those residues.

---

## Limitations

Fr-TM-Align is limited to 3000 total residues

Frequent changes to the OPM server may mean that this method of auto-downloading from OPM will become obsolete. If this happens, you may manually download and save structures into a provided_struct directory that you specify in paths.txt.  Alternatively, the code will attempt to aquire the structures from the PDB.

The HOLE program only accepts filenames with 70 or fewer characters.  We have used a method to automatically attempt to shorten long file names, which will output an error if it is unable to shorten the filename to a usable length.

The HOLE program requires a radius file specifying the Van der Waals radius for each atom, and the included simple2.rad file covers the atoms found in the TRP channel structures to date.  If structures in your analysis contain atoms not defined in the simple2.rad file, HOLE will fail and the HOLE output text file will contain the name of the atom whose radius could not be defined.  To fix this, update the simple2.rad file with the VdW radius of the atom.


## Authors
Katherine Huffer

Antoniya Aleksandrova
