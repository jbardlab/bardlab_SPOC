#!/usr/bin/python3

###############################################################################################################################################################
# This is the code for SPOC (Structure Prediction and Omics Classifier)for evaluating AlphaFold multimer predictions 
# @author: Ernst Schmid
# @date: January 2025 
#If you use this code or its output as part of a publication, please cite the AlphaFold2 multimer manuscript, the ColabFold manuscript and the SPOC manuscript.

# Updated on 04/25/2025 by Jared Bard to only run structure analysis
# This involved a number of changes. Foremost, I had to find code to calculate the pDockQ and pDockQv2 scores.
#    pDockQ: Bryant, Pozotti, and Eloffson. https://www.nature.com/articles/s41467-022-28865-w
#    pDockQ2: Zhu, Shenoy, Kundrotas, Elofsson. https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714
# Luckily, someone had already written a script to calculate these: https://github.com/DunbrackLab/IPSAE/.
# I asked Google Gemini 2.5 pro to help me incorporate the results of this script into the SPOC code (see call_ipsae_and_parse).
# There were also a number of parameters that the structure RF model used that SPOC did not, so i had to calculate those.
# Most of those were just variants of scores that were already calculated, so I implemented those largely inanalyze_complex_nobio.
# I was also aided by the descriptions of each parameter in Table S3 from the predictomes paper: 10.1016/j.molcel.2025.01.034.
###############################################################################################################################################################

import argparse, os, glob, re, random, math, lzma, gzip, csv, json, time, sys
import multiprocessing
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
from scipy.spatial.distance import cosine, euclidean

from joblib import load

#BLAST imports
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import tempfile
import subprocess

#dict for converting 3 letter amino acid code to 1 letter code
aa_3c_to_1c = {
    "ALA":'A',
    "CYS":'C',
    "ASP":'D',
    "GLU":'E',
    "PHE":'F',
    "GLY":'G',
    "HIS":'H',
    "ILE":'I',
    "LYS":'K',
    "LEU":'L',
    "MET":'M',
    "ASN":'N',
    "PRO":'P',
    "GLN":'Q',
    "ARG":'R',
    "SER":'S',
    "THR":'T',
    "VAL":'V',
    "TRP":'W',
    "TYR":'Y',
}
canonical_amino_acids = set(aa_3c_to_1c.values())

def call_ipsae_and_parse(ipsae_script_path, pdb_filepath, pae_filepath):
    """Calls ipsae.py and parses its output for pDockQ scores."""
    # --- Constants ---
    IPSAE_PAE_CUTOFF = "10"        # Default PAE cutoff for ipsae.py call
    IPSAE_DIST_CUTOFF = "10"       # Default distance cutoff for ipsae.py call
    pdockq = 0.0
    pdockq_v2 = 0.0
    ipsae_success = False

    # Determine the expected output file path from ipsae.py
    pdb_stem = pdb_filepath.replace(".pdb.gz", "").replace(".pdb.xz", "").replace(".pdb", "")
    # Construct the expected output filename based on ipsae.py logic
    pae_cutoff_str = IPSAE_PAE_CUTOFF.zfill(2) # e.g., "10" -> "10", "5" -> "05"
    dist_cutoff_str = IPSAE_DIST_CUTOFF.zfill(2)
    ipsae_output_file = f'{pdb_stem}_{pae_cutoff_str}_{dist_cutoff_str}.txt'

    # Construct the command with conda environment activation
    command = [
        '/bin/bash', '-c',
        f'source activate spoc_venv && python {ipsae_script_path} "{pae_filepath}" "{pdb_filepath}" {IPSAE_PAE_CUTOFF} {IPSAE_DIST_CUTOFF}'
    ]

    try:
        print(f"  Running ipsae.py for {os.path.basename(pdb_filepath)}...")
        # Run ipsae.py, capture output/errors, wait for completion
        result = subprocess.run(command, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            print(f"  ERROR: ipsae.py failed for {os.path.basename(pdb_filepath)}.")
            print(f"  Return Code: {result.returncode}")
            print(f"  Stderr: {result.stderr.strip()}")
            # Optionally print stdout too for debugging: print(f"  Stdout: {result.stdout.strip()}")
            # Don't attempt to parse output if script failed
            return pdockq, pdockq_v2, ipsae_success, ipsae_output_file

        # Check if the output file was created
        if not os.path.isfile(ipsae_output_file):
            print(f"  ERROR: ipsae.py ran but output file {ipsae_output_file} not found.")
            print(f"  Stdout: {result.stdout.strip()}") # Print stdout for clues
            return pdockq, pdockq_v2, ipsae_success, ipsae_output_file

        # Parse the output file
        with open(ipsae_output_file, 'r') as f:
            for line in f:
                if not line.strip() or line.startswith("#") or line.startswith("Chn1"):
                    continue # Skip comments, header, and empty lines

                parts = line.split()
                if len(parts) >= 13 and parts[4] == "max": # Identify the 'max' summary line
                    try:
                        # pDockQ is 11th value (index 10), pDockQ2 is 12th value (index 11)
                        # Indices based on the header in ipsae.py:
                        # Chn1 Chn2 PAE Dist Type ipSAE ipSAE_d0chn ipSAE_d0dom ipTM_af ipTM_d0chn pDockQ pDockQ2 LIS ...
                        # 0    1    2   3    4    5     6           7           8       9          10     11      12
                        pdockq = float(parts[10])
                        pdockq_v2 = float(parts[11])
                        ipsae_success = True
                        print(f"    Extracted pDockQ={pdockq:.4f}, pDockQ2={pdockq_v2:.4f}")
                        break # Found the max line, no need to read further
                    except (ValueError, IndexError) as parse_err:
                        print(f"  ERROR: Could not parse pDockQ/pDockQ2 from ipsae.py output line: {line.strip()}")
                        print(f"  Parse Error: {parse_err}")
                        # Continue trying other lines? No, max line should be unique or error out.
                        break

        if not ipsae_success:
            print(f"  ERROR: 'max' summary line not found or parsed in {ipsae_output_file}.")


    except subprocess.TimeoutExpired:
        print(f"  ERROR: ipsae.py timed out for {os.path.basename(pdb_filepath)}.")
    except Exception as e:
        print(f"  ERROR: Unexpected error running or parsing ipsae.py for {os.path.basename(pdb_filepath)}: {e}")

    # Clean up the generated output file? Optional.
    # if ipsae_success and os.path.exists(ipsae_output_file):
    #     try: os.remove(ipsae_output_file)
    #     except OSError as e: print(f"  Warning: Could not remove ipsae output file {ipsae_output_file}: {e}")

    return pdockq, pdockq_v2, ipsae_success, ipsae_output_file


def get_jaccard_index(list_a, list_b):
    """
    Calculate the Jaccard index between two lists.

    :param list_a: First list of binary values.
    :param list_b: Second list of binary values.
    :return: Jaccard index as a float.
    """
    # Ensure both lists are of the same length
    if len(list_a) != len(list_b):
        return 0
    
    intersection_count = 0
    union_count = 0
    
    for a, b in zip(list_a, list_b):
        if a == 1 and b == 1: intersection_count += 1
        if a == 1 or b == 1: union_count += 1
    
    # Prevent division by zero
    if union_count == 0: return 0
    
    jaccard_index = intersection_count / union_count
    return jaccard_index


def cosine_similarity(v1, v2):
    """
    Calculate the cosine similarity between two vectors.

    :param v1: First vector.
    :param v2: Second vector.
    :return: Cosine similarity as a float.
    """
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    cosine_sim = dot_product / (norm_v1 * norm_v2)
    return cosine_sim

def get_contact_type(a1, a2, d):

    if d < 1:
        return 'C'

    r1, a1 = a1.split(':')
    r2, a2 = a2.split(':')

    basic_atoms = ['NH2', 'NZ', 'ND1', 'NE', 'NH1']
    a1b = a1 in basic_atoms
    a2b = a2 in basic_atoms

    if r1 == 'H' and a1 == 'NE2':
        a1b = True

    if r2 == 'H' and a2 == 'NE2':
        a2b = True

    acidic_atoms = ['OE2', 'OD2', 'OXT']
    
    #serine, threonine, tyrosine,
    hbond_donors = ['OG','OG1','OH', 'OE2', 'OD2', 'NE1', 'ND2', 'NE2', 'NZ', 'NE', 'NH1', 'NH2', 'ND1', 'N', 'OXT']
    hbond_acceptors = ['OG','OG1','OH', 'OE1', 'OD1','OE2', 'OD2', 'O', 'NE1']
    
    if (a1b and a2 in acidic_atoms) or (a2b and a1 in acidic_atoms):
        if d <= 5:
            return 'S'
    if (a1 in hbond_donors and a2 in hbond_acceptors) or (a2 in hbond_donors and a1 in hbond_acceptors):
        if d <= 3:
            return 'H'
    if (a1 in acidic_atoms and a2 in acidic_atoms) or (a1 in basic_atoms and a2 in basic_atoms):
        if d <= 5:
            return 'R'
        
    return 'V'

def join_csv_files(files:list, output_name:str, sort_col:str = None, sort_ascending:bool = False, headers = None):
    """
        Join multiple CSV files into a single file.

        :param files (list): A list of file paths to CSV files to be joined.
        :param output_name (str): The name of the output file.
        :param sort_col (str, optional): The column header of the final CSV column by which to sort the rows by.
        :param sort_ascending (bool, optional): The sort direction to use when sorting the final output CSV.
        :param headers (list, optional): A list of column names for the output file. If not provided, the column names from the first input file are used.
    """
    if(len(files) < 1):
        return

    all_dfs = []
    for f in files:
        all_dfs.append(pd.read_csv(f))

    combo_df = pd.concat(all_dfs, ignore_index=True)

    if headers is not None:
        combo_df.columns = headers

    if sort_col:
        combo_df.sort_values(by=[sort_col], ascending=sort_ascending, inplace=True)
    combo_df.to_csv(output_name, index=None)


def distribute(lst:list, n_bins:int) -> list:
    """
        Returns a list containg n_bins number of lists that contains the items passed in with the lst argument

        :param lst: list that contains that items to be distributed across n bins
        :param n_bins: number of bins/lists across which to distribute the items of lst
    """ 
    if n_bins < 1:
       raise ValueError('The number of bins must be greater than 0')
    
    #cannot have empty bins so max number of bin is always less than or equal to list length
    n_bins = min(n_bins, len(lst))
    distributed_lists = []
    for i in range(0, n_bins):
        distributed_lists.append([])
    
    for i, item in enumerate(lst):
        distributed_lists[i%n_bins].append(item)

    return distributed_lists


def get_af_model_num(filename) -> int:
    """
        Returns the Alphafold model number from an input filestring as an int

        :param filename: string representing the filename from which to extract the model number
    """ 
    
    if "model_" not in filename: return 0

    model_num = int(re.findall(r'model_\d+', filename)[0].replace("model_", ''))
    return model_num


def get_seed_num(filename) -> int:
    """
        Returns the seed from an input filestring as an int

        :param filename: string representing the filename from which to extract the model number
    """ 
    
    if "seed_" not in filename: return 0

    seed = int(re.findall(r'seed_\d+', filename)[0].replace("seed_", ''))
    return seed



def get_filepaths_for_complex(path:str, complex_name:str, pattern:str = '*') -> list:
    """
        Helper methdof for returning a list of filepaths (strs) that match the specified GLOB pattern

        :param path: string representing the path/folder to search for complexes
        :param complex_name: string that represents the name of a complex
        :param pattern: string representing the pattern to use when searching for files belonging to complex. Ex: *.pdb, *.json, etc
    """

    glob_str = os.path.join(path, complex_name + pattern)
    return sorted(glob.glob(glob_str))


def get_data_from_json_file(json_filepath) -> list:
    """
        Returns a list of string values representing the pAE(predicated Aligned Error) values stored in the JSON output along with the PTM and IPTM values

        :param json_filepath: string representing the JSON filename from which to extract the PAE values
    """ 

    if not os.path.isfile(json_filepath):
        raise ValueError('Non existing PAE file was specified')

    scores_file = None 
    if(json_filepath.endswith('.xz')):
        scores_file = lzma.open(json_filepath, 'rt')
    elif(json_filepath.endswith('.gz')):
        scores_file = gzip.open(json_filepath,'rt')
    elif (json_filepath.endswith('.json')):
        scores_file = open(json_filepath,'rt')
    else:
        raise ValueError('pAE file with invalid extension cannot be analyzed. Only valid JSON files can be analyzed.')

    #read pae file in as text
    try:
        file_text = scores_file.read()
        pae_index = file_text.find('"pae":')
        ptm_index = file_text.find('"ptm":')
        iptm_index = file_text.find('"iptm":')
        scores_file.close()
    except:
        print("Could not parse json_filepath")
        raise ValueError('Could not parse JSON file')
    
    #Transform string representing 2d array into a 1d array of strings (each string is 1 pAE value). We save time by not unecessarily converting them to numbers before we use them.
    pae_data = file_text[pae_index + 6:file_text.find(']]', pae_index) + 2].replace('[','').replace(']','').split(',')

    if len(pae_data) != int(math.sqrt(len(pae_data)))**2:
        #all valid pAE files consist of an N x N matrice of scores
        raise ValueError('pAE values could not be parsed from files')
    
    ptm = float(file_text[ptm_index + 6:file_text.find(',', ptm_index)])
    iptm = float(file_text[iptm_index + 7:].replace('}', ''))
    return pae_data, ptm, iptm


def dist2(v1, v2) -> float:
    """
        Returns the square of the Euclian distance between 2 vectors carrying 3 values representing positions in the X,Y,Z axis

        :param v1: a vector containing 3 numeric values represening X, Y, and Z coordinates
        :param v2: a vector containing 3 numeric values represening X, Y, and Z coordinates
    """ 

    if len(v1) != 3 or len(v2) != 3:
        raise ValueError('3D coordinates require 3 values')

    return (v1[0] - v2[0])**2 + (v1[1] - v2[1])**2 + (v1[2] - v2[2])**2

def atom_from_pdb_line(atom_line:str) -> dict:
    """
        Parses a single line string in the standard PDB format and returns a list a dict that represents an atom with 3d coordinates and a type(element)

        :param atom_line: string representing a single line entry in a PDB file which contains information about an atom in a protein structure
    """ 
    coordinates = np.array([float(atom_line[30:38]), float(atom_line[38:46]), float(atom_line[46:54])])
    return {"type":atom_line[13:16].strip(),"xyz":coordinates,}

def get_closest_atoms(res1:dict, res2:dict):
    """
        Find the two closest atoms between two residues and returns the minimum distance as well as the closest atoms from each residue.

        :param res1: A dictionary representing the first residue. It should have a key 'atoms' that contains a list of dictionaries, where each dictionary represents an atom with keys 'xyz' (a list of three floats representing the x, y, and z coordinates) and 'name' (a string representing the name of the atom).
        :param res2: A dictionary representing the second residue. It should have a key 'atoms' that contains a list of dictionaries, where each dictionary represents an atom with keys 'xyz' (a list of three floats representing the x, y, and z coordinates) and 'name' (a string representing the name of the atom).
    """
    min_d2 = 1e6
    atoms = [None, None]
    for a1 in res1['atoms']:
        for a2 in res2['atoms']:
            d2 = dist2(a1["xyz"], a2["xyz"])
            if d2 < min_d2:
                min_d2 = d2
                atoms[0] = a1
                atoms[1] = a2
    return (min_d2, atoms)


def get_atom_contacts(res1, res2, d = 8):

    """
    Get atom contacts between two residues.

    :param res1: First residue.
    :param res2: Second residue.
    :param d: Maximum distance for contacts.
    :return: List of contacts and minimum distance.
    """
    contacts = []
    min_distance = 1e6
    min_d2 = d**2
    for a1 in res1['atoms']:
        for a2 in res2['atoms']:
            d2 = dist2(a1["xyz"], a2["xyz"])
            backbone = False
            if a1['type'] in ['C', 'CA', 'O', 'N'] and a2['type'] in ['C', 'CA', 'O', 'N']:
                backbone = True

            if d2 < min_d2:
                min_distance = min(min_distance, d2**0.5)
                t = get_contact_type(res1['type']+':'+a1['type'],res2['type']+':'+a2['type'],d2**0.5)
                contacts.append([backbone, t, d2**0.5])

    return contacts, min_distance


def get_atom(res, t='CA'):

    for a in res['atoms']:
        if a['type'] == t: return a
    return None


def get_res_contact_score(contact):
    """
    Calculate the contact score for a residue.

    :param contact: Contact information.
    :return: Contact score.
    """

    return len(contact['atom_contacts'])*0.5*(sum(contact['plddts']))/(1 + 0.5*sum(contact['paes']))


def get_lines_from_pdb_file(pdb_filepath:str) -> list:
    """
        Returns the contents of a protein databank file (PDB) file as a list of strings

        :param pdb_filepath: string representing the path of the PDB file to open and parse (can handle PDB files that have been compressed via GZIP or LZMA)
    """ 

    if not os.path.isfile(pdb_filepath):
        raise ValueError('Non existing PDB file was specified')

    pdb_file = None 
    if(pdb_filepath.endswith('.xz')):
        pdb_file = lzma.open(pdb_filepath, 'rt')
    elif(pdb_filepath.endswith('.gz')):
        pdb_file = gzip.open(pdb_filepath,'rt')
    elif(pdb_filepath.endswith('.pdb')):
        pdb_file = open(pdb_filepath,'rt')
    else:
        raise ValueError('Unable to parse a PDB file with invalid file extension')
    
    pdb_data = pdb_file.read()
    pdb_file.close()
    return pdb_data.splitlines()


def get_sequences(pdb_filepath:str)->dict:

    """
    Get sequences from a PDB file.

    :param pdb_filepath: PDB filename.
    :return: Dictionary of sequences indexed by chain.
    """
    sequences = {}
    last_chain = None
    for atom_line in get_lines_from_pdb_file(pdb_filepath):
        if atom_line[0:4] != 'ATOM':
            continue
        if(atom_line[13:16].strip() != 'N'):
            continue
        chain = atom_line[20:22].strip()
        if chain != last_chain:
            sequences[chain] = ''

        last_chain = chain

        aa_type = atom_line[17:20]
        if aa_type not in aa_3c_to_1c:
            return None
        
        aa_code = aa_3c_to_1c[aa_type]
        sequences[chain] += aa_code
    
    return sequences 

def get_data_from_structure(pdb_filepath:str, max_distance:float = 8, min_plddt:float = 70, within_chain = False):
    """
        Returns a dict that contains all amino acids between different chains that are in contact and meet the specified criteria

        :param pdb_filepath:string of the filepath to the PDB structure file to be parsed for contacts
        :param max_distance:the maximum allowed distance in Angstroms that the 2 residues must have in order to be considered in contact
        :param min_plddt:the minimum pLDDT(0-100) 2 residues must have to be considered in contact
    """ 

     #holds a more relaxed distance criteria for fast preliminary filtering  
    d2_n_cutoff = (max_distance + 20)**2

    d2_cutoff = max_distance**2
    last_chain = None 
    last_chain_2 = None
    abs_res_index = 0
    chain_index = -1
    chains = []

    residues = []

    #holds 3d coordinates of all amide nitrogens for all residues
    #organized as a 2d list with rows representing chains and columns are residues within the chain
    #Code assumes that the amino nitrogen is always the first atom of a residue encountered (is true in AF-M)
    N_coords = []
    pdb_res_sequence = ''
    plddts = []
    residue = None
    
    for atom_line in get_lines_from_pdb_file(pdb_filepath):
        if atom_line[0:4] != 'ATOM':
            continue
        
        atom_type = atom_line[13:16].strip()
        chain = atom_line[20:22].strip()
        aa_type = aa_3c_to_1c[atom_line[17:20]]
        is_nitrogen = atom_type == 'N'
        #in AlphaFold output PDB files, pLDDT values are stored in the "bfactor" column 
        bfactor = float(atom_line[60:66])
        if(is_nitrogen):
            #Keep track of what the absolute index of this residue in the file is 
            abs_res_index += 1
            if chain != last_chain_2: 
                if last_chain_2 is not None:
                    pdb_res_sequence += ':'
                last_chain_2 = chain
            
            pdb_res_sequence += aa_type
            plddts.append(bfactor)

        if bfactor < min_plddt:
            #No need to examine this atom since it does not meet the pLDDT cutoff, skip it
            continue

        atom = atom_from_pdb_line(atom_line)
        if is_nitrogen:

            #Every amino acid residue starts PDB entry with exactly one "N" atom, so when we see one we know we have just encountered a new residue
            if chain != last_chain:
                #Are we in a new chain? If so, increment chain index and create new list in "residues"
                chain_index += 1
                last_chain = chain
                N_coords.append([])
                residues.append([])
                chains.append(chain)

            residue = {"chain":chain, "atoms":[],'c_ix':int(atom_line[22:26]), "a_ix":abs_res_index, "type":aa_type, "plddt":bfactor}
            residues[chain_index].append(residue)

            #add nitrogen atom coordinates to coordinates list to allow for fast broad searching later
            N_coords[chain_index].append(atom['xyz'])

        if residue is not None:
            residue['atoms'].append(atom)
    
    contacts = []
    num_chains = len(chains)

    #loop through all the protein chains to find contacts between chains
    #strategy is to first look for residues in general proximity by just looking at the distance between their amide nitrogen
    for i in range(0, num_chains):
        chain_1_coords = N_coords[i]
        num_in_c1 = len(chain_1_coords)

        i2_start = i if within_chain else i + 1
        for i2 in range(i2_start, num_chains):

            chain_2_coords = N_coords[i2]
            num_in_c2 = len(chain_2_coords)

            #construct 2 3D numpy arrays to hold coordinates of all residue amide nitorgens
            c1_matrix = np.tile(chain_1_coords, (1, num_in_c2)).reshape(num_in_c1, num_in_c2, 3)
            c2_matrix = np.tile(chain_2_coords, (num_in_c1, 1)).reshape(num_in_c1, num_in_c2, 3)

            #calculate euclidian distance squared (faster) between all amide nitorgens of all residues
            d2s = np.sum((c1_matrix - c2_matrix)**2, axis=2)
            #get residue pairs where amide nitrogens are closer than the initial broad cutoff
            index_pairs = list(zip(*np.where(d2s < d2_n_cutoff)))

            #find closest atoms between residues that were found to be somewhat in proximity
            for c1_res_ix, c2_res_ix in index_pairs:
                
                r1 = residues[i][c1_res_ix]
                r2 = residues[i2][c2_res_ix]

                atom_contacts, min_d = get_atom_contacts(r1, r2, max_distance)
                if len(atom_contacts) > 0:

                    clashing = False
                    for ac in atom_contacts:
                        if ac[1] == 'C':
                            clashing = True
                            break

                    r1_ca = get_atom(r1)
                    r2_ca = get_atom(r2)

                    if r1_ca is None or r2_ca is None:
                        print("NO CA FOUND")
                        continue

                    #residues have atoms closer than specified cutoff, lets add them to the list
                    contacts.append({
                        'distance':min_d,
                        'atom_contacts':atom_contacts,
                        'clashing':clashing,
                        "aa1":{"chain":r1["chain"], "ca":r1_ca, "type":r1["type"], "c_ix":r1['c_ix'], "a_ix":r1['a_ix'], "plddt": r1["plddt"]},
                        "aa2":{"chain":r2["chain"], "ca":r2_ca, "type":r2["type"], "c_ix":r2['c_ix'], "a_ix":r2['a_ix'], "plddt": r2["plddt"]}
                    })
    return contacts, pdb_res_sequence, plddts


def get_data_nobio(pdb_filepath:str, pae_filepath:str, ipsae_script_path=str, max_distance:float=5,min_plddt:float=50, max_pae:float=15) -> dict:
    """
        Get contacts from a protein structure in PDB format that meet the specified distance and confidence criteria.

        :param pdb_filepath (str): The path to the PDB file.
        :param pae_filepath (str): The path to the predicted Alignment Error (pAE) file.
        :param max_distance (float): The maximum distance between two atoms for them to be considered in contact.
        :param min_plddt (float): The minimum PLDDT score required for a residue to be considered "well-modeled".
        :param max_pae (float): The maximum predicted Alignment Error allowed for a residue to be considered "well-modeled".
        :param pae_mode (str): The method to use for calculating predicted atomic error (pAE). Possible values are "avg" or "min".
    """

    # --- Call ipsae.py ---
    print("  Calling ipsae.py for pDockQ scores...")
    pdockq_score, pdockq_v2_score, ipsae_success, _ = call_ipsae_and_parse(ipsae_script_path,pdb_filepath, pae_filepath)
    print(f"pdockq_score: {pdockq_score}, pdockq_v2_score: {pdockq_v2_score}")
    # Note: We get the *aggregated* 'max' scores directly from ipsae.py output parser

    model_num = get_af_model_num(pdb_filepath)
    if model_num < 1 or model_num > 5:
        raise ValueError('There are only 5 Alphafold models, numbered 1 to 5. All PDB files and PAE files must have a valid AlphaFold model number to be analyzed.')

    if model_num != get_af_model_num(pae_filepath):
        raise ValueError(f"File mismatch, can only compare PDB and PAE files from same complex and the same AF2 model PDB:({pdb_filepath}) PAE:({pae_filepath})")

    try:
        pae_data, ptm, iptm = get_data_from_json_file(pae_filepath)
    except Exception as e:
        print(f"Exception occurred during JSON file parsing {e}")
        return None

    total_aa_length = int(math.sqrt(len(pae_data)))

    try:
        unfiltered_contacts, aa_sequence, unfiltered_plddts = get_data_from_structure(pdb_filepath, max_distance, 0, False)
        #first determine which residues are in physical contact(distance) and have a minimum pLDDT score (bfactor column)
        contacts, aa_sequence, plddts = get_data_from_structure(pdb_filepath, max_distance, min_plddt, False)
    except Exception as e:
        print(f"Exception occurred during PDB file parsing: {e}")
        return None

    unfiltered_plddts = np.array(unfiltered_plddts)
    plddt_min = np.min(unfiltered_plddts)
    plddt_max = np.max(unfiltered_plddts)
    plddt_mean = np.mean(unfiltered_plddts)
    
    clashing_residue_indices = {}
    for c in contacts:
        if c['clashing']:
            clashing_residue_indices[c['aa1']['a_ix']] = 1
            clashing_residue_indices[c['aa2']['a_ix']] = 1

    new_contacts = []
    total_atom_contacts = 0
    for c in contacts:
        if c['aa1']['a_ix'] in clashing_residue_indices or c['aa2']['a_ix'] in clashing_residue_indices:
            continue
        new_contacts.append(c) 
        total_atom_contacts += len(c['atom_contacts'])

    contacts = new_contacts

    unfiltered_if_paes = []
    for c in unfiltered_contacts:

        aas = [c['aa1'], c['aa2']]
        aa_indices = [aas[0]['a_ix'],aas[1]['a_ix']]
        pae_index_1 = total_aa_length*(aa_indices[0] - 1) + aa_indices[1] - 1
        pae_index_2 = total_aa_length*(aa_indices[1] - 1) + aa_indices[0] - 1

        if pae_index_1 >= len(pae_data) or pae_index_2 >= len(pae_data):
            raise ValueError(f"Something went wrong and we are attempting to access non-existant PAE values for PDB file: {pdb_filepath} from PAE file: {pae_filepath}")

        unfiltered_if_paes += [float(pae_data[pae_index_1]), float(pae_data[pae_index_2])]


    if len(unfiltered_if_paes) == 0: unfiltered_if_paes = [30]

    if len(contacts) < 1:
        return {'total_length':total_aa_length, 
                'contacts':{},
                'ptm':ptm, 
                'iptm':iptm, 
                'unfiltered_if_pae_mean':round(np.mean(np.array(unfiltered_if_paes)), 1),
                'plddt_min': plddt_min,
                'plddt_max': plddt_max,
                'plddt_mean':plddt_mean,
                'disorder_percent':round(100*((unfiltered_plddts < 50).sum()/total_aa_length), 0),
                'aa_sequence':aa_sequence,
                'pdockq': pdockq_score,
                'pdockq_v2': pdockq_v2_score}
    
    filtered_contacts = {}

    for c in contacts:

        aas = [c['aa1'], c['aa2']]
        aa_indices = [aas[0]['a_ix'],aas[1]['a_ix']]

        pae_values = [0, 0]
        pae_value = 0

        pae_index_1 = total_aa_length*(aa_indices[0] - 1) + aa_indices[1] - 1
        pae_index_2 = total_aa_length*(aa_indices[1] - 1) + aa_indices[0] - 1

        if pae_index_1 >= len(pae_data) or pae_index_2 >= len(pae_data):
            raise ValueError(f"Something went wrong and we are attempting to access non-existant PAE values for PDB file: {pdb_filepath} from PAE file: {pae_filepath}")

        #pae data contains string values, have to convert them to floats before using them for math calculations
        pae_values = [float(pae_data[pae_index_1]), float(pae_data[pae_index_2])]
        pae_value = min(pae_values[0],pae_values[1])
        
        if(pae_value > max_pae):
            #The PAE value of this residue pair is too high, skip it
            continue

        #Use the 2 chains IDS as a key
        chain_contact_id = aas[0]['chain'] + ":" + aas[1]['chain']
        if chain_contact_id not in filtered_contacts:
            filtered_contacts[chain_contact_id] = {}


        #Use the absolute indices of the two residues in the PDB file as the unique key for this pair/contact
        contact_id = str(aa_indices[0]) + '&' + str(aa_indices[1])
        filtered_contacts[chain_contact_id][contact_id] = {
            'chains':[aas[0]['chain'], aas[1]['chain']],
            'indices':[aas[0]['a_ix'], aas[1]['a_ix']],
            'inchain_indices':[aas[0]['c_ix'], aas[1]['c_ix']],
            'types':[aas[0]['type'], aas[1]['type']], 
            'cas':[aas[0]['ca'], aas[1]['ca']],
            'pae':pae_value,
            'paes':pae_values,
            'plddts':[aas[0]['plddt'], aas[1]['plddt']], 
            'model':model_num,
            'distance':c['distance'],
            'atom_contacts':c['atom_contacts']
        }

    return {'total_length':total_aa_length,
            'contacts':filtered_contacts, 
            'ptm':ptm, 'iptm':iptm,
            'unfiltered_if_pae_mean':round(np.mean(np.array(unfiltered_if_paes)), 1),
            'plddt_min': plddt_min,
            'plddt_max': plddt_max,
            'plddt_mean':plddt_mean,
            'disorder_percent':round(100*((unfiltered_plddts < 50).sum()/total_aa_length), 0),
            'aa_sequence':aa_sequence,
            'pdockq': pdockq_score,
            'pdockq_v2': pdockq_v2_score}

        
def calculate_interface_statistics_nobio(contacts:dict) -> dict:
    """
        Returns summary confidence statistics such as pAE and pLDDT values across all the contacts in an interface

        :param contacts:dict of contacts in an interface of the form {'chain1:chain2':{'1&400':{plddts:[75,70], paes:[10, 7]}, '4&600':{plddts:[68,77], paes:[8, 3]}}}
    """ 

    #plddts always range from 0 to 100
    plddt_sum = 0 
    plddt_min = 100
    plddt_max = 0
    plddt_avg = 0
   
    #paes always range from 0 to 30
    pae_avg = 0
    pae_sum = 0
    pae_min = 30
    pae_max = 0
    distance_avg = 0

    num_contacts = 0
    d_sum = 0
    clash_count = 0

    unique_residues = {}
    contact_scores = np.array([])

    num_significant_afm = 0
    plddt_diffs = []

    for interchain_id, interchain_contacts in contacts.items():
        for contact_id, contact in interchain_contacts.items():

            resix_1, resix_2 = contact_id.split('&')
            unique_residues[resix_1] = 1
            unique_residues[resix_2] = 1

            avg_plddt = np.mean(np.array(contact['plddts']))
            plddt_sum += avg_plddt
            plddt_max = max(plddt_max, avg_plddt)
            plddt_min = min(plddt_min, avg_plddt)

            d_sum += contact['distance']
            clash_count += 1 if contact['distance'] < 1 else 0

            pae_max = max(pae_max, contact['pae'])
            pae_min = min(pae_min, contact['pae'])
            pae_sum += contact['pae']

            score = get_res_contact_score(contact)
            contact_scores = np.append(contact_scores, score)
            num_contacts += 1
            plddt_diffs.append(abs(contact['plddts'][0] - contact['plddts'][1]))



    if len(contact_scores) < 1:
        contact_scores = np.append(contact_scores, 0)

    if num_contacts > 0:
        plddt_avg = round(plddt_sum/num_contacts, 1)
        pae_avg = round(pae_sum/num_contacts, 1)
        distance_avg = round(d_sum/num_contacts, 1)
        plddt_diff_avg = round(np.mean(plddt_diffs), 1)
    else:
        pae_min = 0
        plddt_min = 0
        distance_avg = 0
        plddt_diff_avg = 0
        plddt_avg = 0
        pae_avg = 0

    data = {'num_residue_contacts':num_contacts,
            'num_residues':len(unique_residues),
            'clash_count':clash_count,
            'plddt':[plddt_min, plddt_avg, plddt_max],
            'plddt_avg':plddt_avg,
            'plddt_diff_avg':plddt_diff_avg,
            'pae':[pae_min, pae_avg, pae_max],
            'distance_avg': distance_avg,
            'contact_score_avg':round(np.mean(contact_scores), 2),
            'contact_score_median':round(np.median(contact_scores), 2),
            'contact_score_max':round(np.max(contact_scores), 2),
            'contact_score_avg':round(np.mean(contact_scores), 2),
            'num_significant_afm_scores':num_significant_afm
    }

    return data

def summarize_interface_statistics_nobio(interfaces:dict) -> dict:
    """
        summarize_interface_statistics returns aggregate statistics over multiple interfaces across predictions from different models

        :param interfaces:dict of interfaces in the form 
            {1: {'chain1:chain2':{
                                '1&400':{'plddts':[75,70], 'paes':[10, 7]}, 
                                '4&600':{'plddts':[68,77], 'paes':[8, 3]}
                                }
                }, 
            4: {'chain1:chain2':{
                                '13&400':{'plddts':[77,91], 'paes':[5, 7]}, 
                                '49&600':{'plddts':[68,56], 'paes':[9, 3]}
                                }
                },     
            }
    """
    
    unique_contacts = {}
    max_num_models = 0
    num_models_run = 0

    for model_num, interchain_interfaces in interfaces.items():
        num_models_run += 1

        for interchain_str, contacts in interchain_interfaces.items():
            for contact_id, c in contacts.items():

                if contact_id not in unique_contacts:
                    unique_contacts[contact_id] = 1
                else:
                    unique_contacts[contact_id] += 1

                max_num_models = max(max_num_models, unique_contacts[contact_id])

    num_contacts = 0
    sum_num_models = 0
    num_contacts_with_max_n_models = 0
    for contact_id, observation_count in unique_contacts.items():

        num_contacts += 1
        sum_num_models += observation_count

        if observation_count == max_num_models:
            num_contacts_with_max_n_models += 1
    
    summary_stats = {
        'avg_n_models':round(sum_num_models/(num_contacts*num_models_run), 3) if num_contacts > 0 else 0,
        'max_n_models':round(max_num_models/num_models_run, 3),
        'num_contacts_with_max_n_models':num_contacts_with_max_n_models,
        'num_unique_contacts':num_contacts
    }
    return summary_stats

    
def analyze_complex_nobio(complex_name:str, complexes:dict, ipsae_script_path:str, single_complex_mode:bool=False):

    try:
        final_data = None
        start_time = time.time()

        interface_contacts = {}
        models = list(complexes[complex_name].keys())
        
        # Modified logic for single complex mode
        if single_complex_mode and len(models) == 1:
            # In single complex mode, use the same model 3 times
            single_model = models[0]
            models = [single_model, single_model, single_model]
            print(f'Single complex mode: Using model {single_model} three times for analysis')
        elif len(models) < 3: 
            print(f'ERROR: Only {len(models)} models were supplied for {complex_name}. We require at least 3 models per complex.')
            return complex_name, None
        elif len(models) > 3: 
            models = random.sample(models, 3)
            print(f"Found more than 3 models, randomly sampled 3 models {models}")

        # Get the actual model data (handling repeated models in single complex mode)
        model_data = []
        for model_num in models:
            if model_num in complexes[complex_name]:
                model_data.append((model_num, complexes[complex_name][model_num]))
            else:
                # This shouldn't happen but handle it gracefully
                print(f"Warning: Model {model_num} not found in complex data")
                return complex_name, None

        aa_sequences = get_sequences(model_data[0][1][0])  # Use first model's PDB
        if (len(aa_sequences) != 2):
            print(f'ERROR: Our pipeline currently only supports 2 chain complexes. Found {len(aa_sequences)} chains.')
            return complex_name, None

        interface_contacts = {}
        iptm_values = np.array([])
        pdockq_values = np.array([])
        pdockq_v2_values = np.array([])
        residue_contacts_across_predictions = np.array([])
        best_interface_stats = None

        # Process each model (may be the same model repeated in single complex mode)
        for idx, (model_num, (pdb_filepath, pae_filepath)) in enumerate(model_data):
            if single_complex_mode:
                print(f"Getting Data (iteration {idx+1}/3 of same model)")
            else:
                print(f"Getting Data for model {model_num}")
                
            data = get_data_nobio(pdb_filepath, pae_filepath, ipsae_script_path)
            contacts = data['contacts']
            
            # Use a unique key for interface_contacts even in single complex mode
            contact_key = f"{model_num}_{idx}" if single_complex_mode else model_num
            interface_contacts[contact_key] = contacts
            
            print("Calculating interface stats")
            if_stats = calculate_interface_statistics_nobio(contacts)
            
            residue_contacts_across_predictions = np.append(residue_contacts_across_predictions, if_stats['num_residue_contacts'])
            iptm_values = np.append(iptm_values, data['iptm'])
            pdockq_values = np.append(pdockq_values, data['pdockq'])
            pdockq_v2_values = np.append(pdockq_v2_values, data['pdockq_v2'])

            if_stats['ptm'] = data['ptm']
            if_stats['iptm'] = data['iptm']
            if_stats['unfiltered_if_pae_mean'] = data['unfiltered_if_pae_mean']
            if_stats['all_plddt_mean'] = data['plddt_mean']
            if_stats['all_plddt_min'] = data['plddt_min']
            if_stats['all_plddt_max'] = data['plddt_max']
            if_stats['all_disorder_percent'] =  data['disorder_percent']
            if_stats['if_residues_percent'] = round(if_stats['num_residues']/data['total_length'], 3) if data['total_length'] > 0 else 0
            if_stats['contacts_per_pae'] = round(if_stats['num_residue_contacts']/(if_stats['pae'][1] + 1), 3)
            if_stats['total_length'] = data['total_length']
            if_stats['pdockq'] = data['pdockq']

            if best_interface_stats is None:
                best_interface_stats = if_stats
                best_interface_stats['model_num'] = model_num
            else:
                if if_stats['contacts_per_pae'] > best_interface_stats['contacts_per_pae']:
                    best_interface_stats = if_stats
                    best_interface_stats['model_num'] = model_num

        print("Summarizing interface stats")
        stats = summarize_interface_statistics_nobio(interface_contacts)
        stats['mean_contacts_across_predictions'] = round(np.mean(residue_contacts_across_predictions), 0)
        stats['min_contacts_across_predictions'] = np.min(residue_contacts_across_predictions)
        stats['iptm_mean'] = round(np.mean(iptm_values), 3)
        stats['iptm_min'] = round(np.min(iptm_values), 3)
        stats['iptm_max'] = round(np.max(iptm_values), 3)
        stats['pdockq_e_max'] = round(np.max(pdockq_values), 4)
        stats['pdockq_e_v2_max'] = round(np.max(pdockq_v2_values), 4)
        stats['models_checked'] = models if not single_complex_mode else [models[0]]  # Show actual model(s) checked
        stats['single_complex_mode'] = single_complex_mode  # Add flag to output
        stats['best_num_residue_contacts'] = best_interface_stats['num_residue_contacts']
        stats['best_if_residues'] = best_interface_stats['num_residues']
        stats['best_if_residues_percent'] = best_interface_stats['if_residues_percent']
        stats['best_plddt_max'] = round(best_interface_stats['plddt'][2], 0)
        stats['best_plddt_avg'] = round(best_interface_stats['plddt_avg'], 0)
        stats['best_plddt_diff_mean'] = round(best_interface_stats['plddt_diff_avg'], 0)
        stats['best_pae_min'] = best_interface_stats['pae'][0]
        stats['best_unfiltered_if_pae_mean'] = best_interface_stats['pae'][1]
        stats['best_contact_score_max'] = round(best_interface_stats['contact_score_max'], 0)
        stats['best_contact_score_avg'] = round(best_interface_stats['contact_score_avg'], 0)
        stats['best_contact_score_median'] = round(best_interface_stats['contact_score_median'], 0)
        stats['best_num_significant_afm_scores'] = best_interface_stats['num_significant_afm_scores']
        stats['best_pdockq_e'] = round(best_interface_stats['pdockq'], 4)
        stats['best_all_disorder_percent'] = best_interface_stats['all_disorder_percent']
        stats['total_aa_length'] = best_interface_stats['total_length']
        stats['sequence'] = data['aa_sequence']

        final_data = stats

        end_time = time.time()
        print(f'Finished: {complex_name}. Took {round(end_time - start_time, 1)} seconds')

        return complex_name, final_data
    
    except Exception as e:
        print(f'ERROR: Exception {e} occurred for complex {complex_name}. Skipping. ')
        return complex_name, None


def analyze_complexes_nobio(complexes, ipsae_script_path:str, single_complex_mode:bool=False):
    all_data = {}
    errored_complexes = []

    print(f"Splitting analysis across {multiprocessing.cpu_count()} CPUs")
    
    # Create a pool of workers equal to the number of CPU cores
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        # Map process_complex function to each item in complexes dictionary
        results = pool.starmap(analyze_complex_nobio, 
                             [(name, complexes, ipsae_script_path, single_complex_mode) for name in complexes])

    # Collect results from all processes
    for complex_name, data in results:
        if data is None:
            errored_complexes.append(complex_name)
        else:
            all_data[complex_name] = data

    return all_data, errored_complexes


def main(folder_paths:list, name_filter:str, classifier, output_name:str, ipsae_script_path:str, single_complex_mode:bool=False):

    """
    Main function to run the analysis.

    :param folder_paths: List of folder paths containing data.
    :param name_filter: Filter for complex names.
    :param classifier: Random forest classifier model.
    :param output_name: Name for the output file. If not set will use a default name.
    :param ipsae_script_path: Path to the ipsae script.
    :param single_complex_mode: If True, allows analysis with a single model repeated 3 times.
    """

    #First find all the complexes that can be analyzed across all folders
    complexes = {}

    for folder_path in folder_paths:

        pdb_file_paths = glob.glob(os.path.join(folder_path, '*.pdb')) + glob.glob(os.path.join(folder_path, '*.pdb.??'))
        pae_file_paths = glob.glob(os.path.join(folder_path, '*.json')) + glob.glob(os.path.join(folder_path, '*.json.??'))

        for pdb_file_path in pdb_file_paths:

            if name_filter and name_filter not in pdb_file_path:
                continue

            pdb_filename = os.path.basename(pdb_file_path)
            if pdb_filename.split('.').pop() not in ['gz', 'xz', 'pdb']:
                continue

            model_num = get_af_model_num(pdb_file_path)
            if model_num < 1 or model_num > 5:
                continue

            seed_num = get_seed_num(pdb_file_path)

            complex_name = None

            if '_unrelaxed_' in pdb_filename:
                complex_name = pdb_filename.split('_unrelaxed_')[0]

            pae_filepath = None
            for f in pae_file_paths:
                if complex_name in f and f"model_{model_num}" in f and f"seed_{seed_num:03d}" in f: 
                    pae_filepath = f
                    break

            if pae_filepath is None:
                continue

            if pae_filepath.split('.').pop() not in ['gz', 'xz', 'json']:
                continue

            if complex_name not in complexes:
                complexes[complex_name] = {}
            
            complexes[complex_name][model_num] = [pdb_file_path, pae_filepath]

    # Check if we're in single complex mode and have appropriate data
    if single_complex_mode:
        for complex_name in list(complexes.keys()):
            if len(complexes[complex_name]) != 1:
                if len(complexes[complex_name]) > 1:
                    print(f"Warning: {complex_name} has {len(complexes[complex_name])} models in single complex mode. Using first model only.")
                    # Keep only the first model
                    first_model = min(complexes[complex_name].keys())
                    complexes[complex_name] = {first_model: complexes[complex_name][first_model]}
                # If 0 models, it will be caught in analyze_complex_nobio

    if len(complexes) > 0:

        print(f"Found {len(complexes)} complexes to analyze across {len(folder_paths)} folder(s).")
        if single_complex_mode:
            print("Running in SINGLE COMPLEX MODE - each complex will be analyzed using the same model 3 times")

        data, errored_complexes = analyze_complexes_nobio(complexes, ipsae_script_path, single_complex_mode)

        print(f"Errors occurred for {len(errored_complexes)} complexes.")
        if len(errored_complexes) > 0:
            if output_name: 
                if not output_name.endswith('.csv'):
                    output_name += '.csv'
                final_error_csv_path = output_name.replace('.csv', '_errored_complexes.csv')
            # Convert list to DataFrame
            df = pd.DataFrame(errored_complexes, columns=['complex_name'])
            # Write DataFrame to CSV
            df.to_csv(final_error_csv_path, index=False)  # Set index=False to not write row numbers
            print(f"Wrote out errored complexes to {final_error_csv_path}.")

        print("")
        if len(data) < 1:
            print("No complexes yielded any data, the script will terminate now.")
            return

        print("---------------------")
        for cname in data:
            cdata = data[cname]
            if type(cdata) == str: #only occurs when there is an error and the string error message has been assigned
                continue

            if cdata['num_unique_contacts'] < 5:
                #SPOC was only trained to work with complexes with at least 5 unique contacts
                print(f"Less than 5 contacts detected for {cname}. Assigning a SPOC_nobio score of 0.")
                cdata['SPOC_nobio_score'] = 0
                continue

            print(f"Running SPOC_nobio on {cname}")

            data_for_rf = {}
            for f in classifier.feature_names_in_:
                if f != 'single_complex_mode':  # Exclude this flag from RF features
                    data_for_rf[f] = cdata[f]

            data_for_rf = pd.DataFrame(data_for_rf, index=[0])
            SPOC_nobio_score = round(float(classifier.predict_proba(data_for_rf)[:, 1][0]), 3)
            cdata['SPOC_nobio_score'] = SPOC_nobio_score
            print(f"SPOC_nobio score = {SPOC_nobio_score}")


        final_datalist = []

        for cname in data:
            data[cname]['complex_name'] = cname
            final_datalist.append(data[cname])
        
        columns = list(final_datalist[0].keys())
        columns.remove('complex_name')
        columns.remove('SPOC_nobio_score')
        if 'single_complex_mode' in columns:
            columns.remove('single_complex_mode')
        columns = ['complex_name', 'SPOC_nobio_score', 'single_complex_mode'] + columns

        final_df = pd.DataFrame(final_datalist)
        final_df = final_df[columns]

        folder_name = folder_path.rstrip('/').split('/').pop()
        final_csv_path = folder_name + "_SPOC_nobio_analysis.csv"
        if output_name: 
            if not output_name.endswith('.csv'):
                output_name += '.csv'
            final_csv_path = output_name
  
        final_df.sort_values(by=['SPOC_nobio_score'], ascending=False, inplace=True)
        final_df.to_csv(final_csv_path, index=None)
        print("")
        print("---------------------")
        print(f"Wrote analysis data out to {final_csv_path}")

    else:
        print('ERROR: No complexes to analyze found')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run the SPOC_nobio classifier on data in the specified folder which should contain PDB and JSON PAE files.")
    # Define a positional argument for the folder path
    parser.add_argument(
        "folder_paths",
        type=str,
        nargs="+",
        help="Paths to folders containing the data (accepts glob patterns natively).",
    )
    parser.add_argument(
        "--name_filter",
        type=str,
        default=None,
        help="Analyze only complexes with names containing the specified string. I.e. 'MCM' will only analyze complexes with MCM in their name."
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Name for the final analysis output file."
    )
    parser.add_argument(
        "--rf_params",
        type=str,
        default="models/rf_afm_no_bio.joblib",
        help="file containing the random forest parameters"
    )
    parser.add_argument("--ipsae_script", type=str, required=True, help="Path to ipsae.py script.")
    
    # Add new argument for single complex mode
    parser.add_argument(
        "--single_complex",
        action="store_true",
        default=False,
        help="Enable single complex mode: analyze complexes with only one model by repeating it 3 times"
    )
    
    # Parse arguments
    args = parser.parse_args()
    if len(args.folder_paths) < 0:
        print("No folders to analyze were provided. Script is terminating.")
        exit()


    start_message = """
-------------------------------------------------------------------------------------------------------------------------
Thank you for using our modified SPOC_nobio analysis script.
Please note that SPOC_nobio is designed to work with AlphaFold2 multimer model data. 
-------------------------------------------------------------------------------------------------------------------------  
"""

    print(start_message)

    # Load SPOC_nobio random forest classifier mode
    script_path = os.path.realpath(__file__)
    SPOC_nobio_rf_path = os.path.normpath(args.rf_params)
    print(f"Loading SPOC_nobio Random Forest model from {SPOC_nobio_rf_path}")
    rf_classifier = load(SPOC_nobio_rf_path)
    print("Finished loading SPOC_nobio Random Forest model")
    
    if args.single_complex:
        print("*** SINGLE COMPLEX MODE ENABLED ***")
        print("Complexes with only one model will be analyzed by repeating the model 3 times")
    
    print("-----------------------------------")
    print("")

    ipsae_script_path = args.ipsae_script
    if not os.path.isfile(ipsae_script_path):
         script_dir = os.path.dirname(os.path.realpath(__file__))
         alt_path = os.path.join(script_dir, args.ipsae_script)
         if os.path.isfile(alt_path): ipsae_script_path = alt_path
         else: print(f"ERROR: ipsae.py script not found: '{args.ipsae_script}' or '{alt_path}'"); exit(1)

    main_start_time = time.time()
    main(args.folder_paths, args.name_filter, rf_classifier, args.output, ipsae_script_path, args.single_complex)
    main_end_time = time.time()
    print(f'SPOC_nobio Analysis complete. Took a total of {round(main_end_time - main_start_time, 1)} seconds')
    print("-----------------------------------")
    print("")
    print(f'If you use this code or its output as part of a publication, please cite the AlphaFold2 multimer manuscript, the ColabFold manuscript and the SPOC manuscript.')