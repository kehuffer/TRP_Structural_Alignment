import os
import sys
from frtmalign_2_msa import *

if len(sys.argv) < 2:
    raise SystemExit("Usage: python %s <path to paths.txt>" % sys.argv[0])
else:
    paths_file = sys.argv[1].strip()
# set working directory, output directories
# identify files containing necessary information, including paths to program files (frtmalign, HOLE, pyali), information for each structure
print('Info: Reading paths.txt.')
paths = paths_dic(paths_file)
if not os.path.isfile(paths['structs_info']):
    raise SystemExit('Error: File for structure info ' + paths['structs_info'] + ' does not exist.')

# establish structures to include in the analysis, as well as chain order and residues to be kept when cleaning files
# Note that paths['structs_info'] is an xml file containing information for each structure, including: pdbid, chain order, starting and ending residue numbers for transmembrane domain, and categories (subfamily, environment, method, and ligand)
print('Info: Reading xml file with structure info.')
struct_info_dic, struct_info_df = xml_parser(paths['structs_info'])

# acquire structures from OPM based on list of pdb codes and save to paths['pdb_dir']
# perform cleaning based on residue ranges specified in paths['structs_info'] xml file
# save cleaned structures to paths['clean_dir']
print('Info: Acquiring and preprocessing structures.')
for pdb_id, value in struct_info_dic.items():
    success = get_struct(pdb_id, paths['pdb_dir'], paths['provided_struct'])
    if success:
        num = strip_tm_chains(paths['clean_dir'], pdb_id, paths['pdb_dir'] + pdb_id + '.pdb', struct_info_dic[pdb_id]['tm_chains'])

# perform all-vs-all pairwise alignments using Fr-TM-Align on all structures in the paths['clean_dir']
# save Fr-TM-Align output, including aligned structures, transformation matrices, and text files, to paths['frtmalign_dir']
# create separate folders for each stationary structure within paths['frtmalign_dir']
# file names include the mobile PDB ID first and the stationary PDB ID second
print('Info: Running Fr-TM-Align to align all structures pairwise.')
batch_frtmalign(paths['clean_dir'], paths['frtmalign_dir'], paths['frtmalign'], paths['pdb_dir'], paths['clean_dir'])

# process Fr-TM-Align text files to extract data for all pairwise alignments (TM-score, RMSD, aligned length, sequence identity)
# save all outputs to paths['frtmalign_dir']
# output data to .csv files for further analysis
# create clustermap as in Figure 2, with TM scores clustered along the stationary axis and with color codes for qualitative categories (subfamily, environment, method, ligand)
# saves post-clustering list of pdbids to a .csv file for later use
print('Info: Processing, organizing, and saving Fr-TM-Align output.')
frtmalign_list = frtmalign_2_list(paths['frtmalign_dir'])
frtmalign_2_tables(frtmalign_list, struct_info_df, paths['frtmalign_dir'])

# create multiple sequence alignments, using each stationary structure as reference
# save alignment files in paths['frtmalign_dir']
# name alignment files based on reference (i.e. common stationary) structure
# full alignment includes all residues from all structures
# nogap alignment removes all gaps/insertions from reference structure
print('Info: Constructing multiple sequence alignments from pairwise structural alignments.')
batch_align_merger(paths['frtmalign_dir'], paths['frtmalign_dir']+'stationary_clustering_order.csv')

# perform HOLE pore radius analysis for each structure and determine minimum radius of each residue
# save HOLE output files in paths['frtmalign_dir']
print('Info: Running permeation pathway HOLE profile analysis on structures aligned to reference structure '+paths['hole_reference_struct'])
batch_hole(paths['frtmalign_dir'], struct_info_df, paths['hole'], paths['hole_reference_struct'], paths['vdw_radius_file'], paths['hole_reference_pore_point'])

# create Jalview annotation file for minimum radius of each residue in full and nogap multiple sequence alignments
# perform DSSP secondary structure analysis for each structure and determine secondary structure of each residue
# create pseudo-fasta file containing secondary structure of each residue in full and nogap multiple sequence alignments
# save output files in paths['frtmalign_dir']
print('Info: Running secondary structure dssp analysis on structures aligned to reference structure '+paths['hole_reference_struct']+' and making HOLE radius annotation files and dssp pseudo-FASTA files.')
batch_annotation(paths['frtmalign_dir'], paths['hole_reference_struct'], paths['norm_max_radius'], paths['clean_dir'])

print('Info: Running identity and similarity calculations on multiple sequence alignment aligned to reference structure '+paths['hole_reference_struct']+'.')
ident_sim_calc(paths['frtmalign_dir'], paths['hole_reference_struct'])

print('Info: Structural alignment is complete.')