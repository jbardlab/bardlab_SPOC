#!/usr/bin/python3

###############################################################################################################################################################
# This is the code for SPOC (Structure Prediction and Omics Classifier)for evaluating AlphaFold multimer predictions 
# @author: Ernst Schmid
# @date: January 2025 
#If you use this code or its output as part of a publication, please cite the AlphaFold2 multimer manuscript, the ColabFold manuscript and the SPOC manuscript.
###############################################################################################################################################################

import argparse, os, glob, re, random, math, lzma, gzip, csv, json, time
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



def initialize_databases():
    global uniprot_id_map, uniprot_id_to_entry_map, uniprot_seq_map, deep_loc2_vectors
    global af_missense, embeddings, depmap_data, crispr_ko_data, biogrid_data, co_expression_data

    script_path = os.path.realpath(__file__)
    db_path = os.path.dirname(script_path) + "/data"
    uniprot_id_map = {}
    uniprot_id_to_entry_map = {}
    uniprot_seq_map = {}

    print("Loading UniProt ID Map")
    with open(f'{db_path}/hs_proteome_uniprot_id_map.tsv', 'r') as file:
        for line in file:
            line = line.strip()
            fields = line.split('\t')
            if len(fields) != 3:
                continue
            uniprot_id, uniprot_entry_name, seq = fields
            uniprot_id_map[uniprot_entry_name] = uniprot_id
            uniprot_id_to_entry_map[uniprot_id] = uniprot_entry_name
            uniprot_seq_map[seq] = uniprot_id

    print("Loaded UniProt ID Map")

    print("Loading DeepLoc2 Localization Predictions")
    deep_loc2_vectors = {}
    with open(f'{db_path}/hs_proteome_deep_loc2.csv', 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        for row in csv_reader:
            entry_name = row[0].strip().replace('"', '')
            if entry_name not in uniprot_id_map:
                continue
            uniprot_id = uniprot_id_map[entry_name]
            deep_loc2_vectors[uniprot_id] = np.array([float(value) for value in row[1:]])

    print("Loaded DeepLoc2 Localization Predictions")

    print("Loading AlphaMissense data for human proteome")
    af_missense = {}
    with open(f'{db_path}/hs_proteome_afmissense_compressed.json', 'r') as file:
        af_missense = json.load(file)

    print("Loaded AlphaMissense data for human proteome")

    print("Loading T5 embedding vectors for human proteome")
    embeddings = {}
    with h5py.File(f'{db_path}/hs_proteome_embeddings.h5', "r") as file:
        for sequence_id, embedding in file.items():
            embeddings[sequence_id] = np.array(embedding)

    print("Loaded T5 embedding vectors for human proteome")

    print("Loading DepMap data for human proteome")
    depmap_data = {}
    with gzip.open(f'{db_path}/hs_proteome_depmap_vectors.gz', 'rt') as file:
        depmap_data = json.load(file)
        for id in depmap_data:
            depmap_data[id] = np.array(depmap_data[id])

    print("Loaded DepMap data for human proteome")

    print("Loading CRISPR KO data from BioORCS for human proteome")
    crispr_ko_data = {}
    with gzip.open(f'{db_path}/hs_CRISPR_DB_KO_binary_20240209.json.gz', 'rt') as file:
        crispr_ko_data = json.load(file)

    print("Loaded CRISPR KO data from BioORCS for human proteome")


    print("Loading BioGrid data for human proteome")
    biogrid_data = {}
    with open(f'{db_path}/hs_biogrid_pairs.json', 'rt') as file:
        biogrid_data = json.load(file)

    print("Loaded BioGrid data for human proteome")

    print("Loading mRNA coexpression data from COXPRESdb for human proteome")
    unique_ids = set()
    co_expression_data = {}
    with gzip.open(f'{db_path}/hs_coexpression_data.csv.gz', 'rt') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        for row in csv_reader:
            if len(row) < 3:
                continue

            e1, e2, score = row
            e1 = e1.strip().replace('"', '')
            e2 = e2.strip().replace('"', '')

            if e1 not in uniprot_id_map:continue
            if e2 not in uniprot_id_map:continue

            uid1 = uniprot_id_map[e1]
            uid2 = uniprot_id_map[e2]

            score = float(score.strip())
            unique_ids.add(uid1)
            unique_ids.add(uid2)
            sorted_uids = ":".join(sorted((uid1, uid2)))
            co_expression_data[sorted_uids] = score

    for uid in unique_ids:
        co_expression_data[f"{uid}:{uid}"] = 10.0

    del unique_ids
    print("Loaded mRNA coexpression data from COXPRESdb for human proteome")


def get_bio_data_for_uniprot_ids(seq_map:dict):

    """
    Retrieve biological data for given UniProt IDs.

    :param seq_map: Dictionary mapping sequence chains to UniProt IDs.
    :return: Dictionary containing various biological data.
    """

    global af_missense, embeddings, depmap_data, crispr_ko_data, biogrid_data, co_expression_data

    chains = list(seq_map.keys())
    data = {
        'crispr_jaccard':0,
        'crispr_shared_hit_count':0,
        'co_expression_score':0,
        't5_embedding_cosine_dist':0,
        't5_embedding_euclidian_dist':0,
        'depmap_cosine_dist':0,
        'depmap_euclidian_dist':1,
        'depmap_abs_diff':10,
        'colocalization_match_score':0,
        'biogrid_detect_count':0
    }

    if len(chains) != 2:
        return None
    
    ids = [seq_map[chains[0]]['uniprot_id'], seq_map[chains[1]]['uniprot_id']]
    ids.sort()

    if ids[0] in crispr_ko_data and ids[1] in crispr_ko_data:
        v1 = crispr_ko_data[ids[0]]
        v2 = crispr_ko_data[ids[1]]
        data['crispr_jaccard'] = get_jaccard_index(v1, v2)

        shared_hits = 0
        for hit1, hit2 in zip(v1,v2):
            if hit1 > 0 and hit2 > 0:
                shared_hits += 1
        data['crispr_shared_hit_count'] = shared_hits

    if ids[0] in embeddings and ids[1] in embeddings:
        v1 = embeddings[ids[0]]
        v2 = embeddings[ids[1]]
        data['t5_embedding_cosine_dist'] = cosine(v1, v2)
        data['t5_embedding_euclidian_dist'] = euclidean(v1, v2)


    if ids[0] in depmap_data and ids[1] in depmap_data:
        v1 = depmap_data[ids[0]]
        v2 = depmap_data[ids[1]]
        data['depmap_cosine_dist'] = cosine(v1, v2)
        data['depmap_euclidian_dist'] = euclidean(v1, v2)
        data['depmap_abs_diff'] = abs(np.mean(v2-v1))


    if ids[0] in deep_loc2_vectors and ids[1] in deep_loc2_vectors:
        v1 = deep_loc2_vectors[ids[0]]
        v2 = deep_loc2_vectors[ids[1]]
        data['colocalization_match_score'] = cosine_similarity(v1, v2)
    
    if ids[0] + ":" + ids[1] in co_expression_data:
          data['co_expression_score'] = co_expression_data[ids[0] + ":" + ids[1]]

    if ids[0] in biogrid_data and ids[1] in biogrid_data[ids[0]]:
        data['biogrid_detect_count'] = biogrid_data[ids[0]][ids[1]]

    return data


def get_af_missense_data(uniprot_id):

    """
    Retrieve AlphaMissense data for a given UniProt ID.

    :param uniprot_id: UniProt ID.
    :return: AlphaMissense data.
    """
    return af_missense.get(uniprot_id)



blastp_cache = {}
def get_best_hs_matches(sequence:str, evalue:float=None, max_hits:int=5):
    """
    Perform a BLAST search to find the best human sequence matches.

    :param sequence: Query sequence.
    :param evalue: E-value threshold for BLAST search.
    :param max_hits: Maximum number of hits to return.
    :return: List of BLAST result records.
    """

    if sequence in blastp_cache:
        return blastp_cache[sequence]

    if evalue is None:
        #Impute evalue to use based on query sequence length
        #Higher e values for shorter sequences as otherwise matches will be missed
        seq_len = len(sequence)
        if seq_len < 30:
            evalue = 10
        elif seq_len < 50:
            evalue = 1
        elif seq_len < 100:
            evalue = 0.01
    
    # Create a temporary file for the query sequence in FASTA format
    temp_fd, temp_path = tempfile.mkstemp(suffix=".fasta")
    with os.fdopen(temp_fd, 'w') as temp_file:
        seq_record = SeqRecord(Seq(sequence), id="query_seq", description="Query Sequence")
        SeqIO.write(seq_record, temp_file, "fasta")
    
    # Paths
    script_path = os.path.realpath(__file__)
    script_path = os.path.dirname(script_path)
    fasta_path = script_path + "/data/hs_proteome.fasta"
    blastp_path = script_path + "/tools/ncbi-blast-2.16.0+/bin/blastp"

    # Temporary output file for BLAST results
    output_fd, output_path = tempfile.mkstemp(suffix=".xml")
    os.close(output_fd)  # Close file descriptor as we only need the path

    print("Running PBLAST search")

    # Command-line blastp invocation
    blast_command = [
        blastp_path, "-query", temp_path, "-db", fasta_path, 
        "-evalue", str(evalue), "-outfmt", "5", "-out", output_path, 
        "-num_threads", str(os.cpu_count()),
        "-max_target_seqs", str(max_hits)
    ]

    # Run the command in the current environment
    subprocess.run(blast_command, env=os.environ, check=True)

    # Parse the BLAST XML output
    blast_result_records = []
    with open(output_path) as result_file:
        blast_records = NCBIXML.parse(result_file)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # Parsing the title to extract UniProt ID and Entry Name
                    title = alignment.title.split(' ', 1).pop()
                    parts = title.replace(' ', '').split('|')
                    uniprot_id = parts[0].strip() if len(parts) == 2 else "N/A"
                    entry_name = parts[1].strip().split(' ')[0] if len(parts) == 2 else "N/A"

                    residue_mapping = {}
                    query_index = hsp.query_start
                    subject_index = hsp.sbjct_start
                    for query_residue, match, subject_residue in zip(hsp.query, hsp.match, hsp.sbjct):
                        if query_residue != '-':
                            if subject_residue != '-':
                                residue_mapping[query_index] = subject_index
                                subject_index += 1
                            query_index += 1

                    if hsp.identities / len(sequence) < 0.90:
                        # less than 90% sequence similarity we ignore it
                        continue

                    blast_result_records.append({
                        "uniprot_id": uniprot_id,
                        "entry_name": entry_name,
                        "length": alignment.length,
                        "e_value": hsp.expect,
                        "query_start": hsp.query_start,
                        "query_end": hsp.query_end,
                        "subject_start": hsp.sbjct_start,
                        "subject_end": hsp.sbjct_end,
                        "identity": hsp.identities,
                        "alignment_length": hsp.align_length,
                        "residue_mapping": residue_mapping
                    })

    # Clean up temporary files
    os.remove(temp_path)
    os.remove(output_path)

    blastp_cache[sequence] = blast_result_records
    return blast_result_records


def get_uniprot_id_for_sequence(seq):

    if seq in uniprot_seq_map:
        uniprot_id = uniprot_seq_map[seq]
        res_map = {ix: ix for ix in range(1, len(seq) + 1)}
        return uniprot_id, res_map

    matches = get_best_hs_matches(seq)
    if len(matches) > 0:
        return matches[0]['uniprot_id'], matches[0]['residue_mapping']
            
    return None, None


def get_seq_to_proteome_map(sequences):

    uniprot_ids = {}
    for chain_id in sequences:
        aa_seq = sequences[chain_id]
        uniprot_id, res_map = get_uniprot_id_for_sequence(aa_seq)
        uniprot_ids[chain_id] = {'uniprot_id':uniprot_id, 'rmap':res_map}

    return uniprot_ids



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



def get_data(pdb_filepath:str, pae_filepath:str, map_data:dict, max_distance:float=5,min_plddt:float=50, max_pae:float=15) -> dict:
    """
        Get contacts from a protein structure in PDB format that meet the specified distance and confidence criteria.

        :param pdb_filepath (str): The path to the PDB file.
        :param pae_filepath (str): The path to the predicted Alignment Error (pAE) file.
        :param max_distance (float): The maximum distance between two atoms for them to be considered in contact.
        :param min_plddt (float): The minimum PLDDT score required for a residue to be considered "well-modeled".
        :param max_pae (float): The maximum predicted Alignment Error allowed for a residue to be considered "well-modeled".
        :param pae_mode (str): The method to use for calculating predicted atomic error (pAE). Possible values are "avg" or "min".
    """

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
                'aa_sequence':aa_sequence}
    
    filtered_contacts = {}

    afmiss_data = {}
    for chain_id, data in map_data.items():
        afmiss_data[chain_id] = get_af_missense_data(data['uniprot_id'])

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


        af_missense_scores = [] 
        for ix, aa in enumerate(aas):
            afm_values = afmiss_data[aa['chain']]
            
            if aa['c_ix'] not in map_data[aa['chain']]['rmap']:
                continue

            canonical_ix = map_data[aa['chain']]['rmap'][aa['c_ix']]
            res_ix = canonical_ix - 1
            if afm_values is None or res_ix > len(afm_values) - 1:
                continue
            
            af_missense_scores.append(afm_values[res_ix])

        if len(af_missense_scores) < 2:
            af_missense_scores = [0, 0]


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
            'atom_contacts':c['atom_contacts'],
            'af_missense':af_missense_scores
        }

    return {'total_length':total_aa_length,
            'contacts':filtered_contacts, 
            'ptm':ptm, 'iptm':iptm,
            'unfiltered_if_pae_mean':round(np.mean(np.array(unfiltered_if_paes)), 1),
            'plddt_min': plddt_min,
            'plddt_max': plddt_max,
            'plddt_mean':plddt_mean,
            'disorder_percent':round(100*((unfiltered_plddts < 50).sum()/total_aa_length), 0),
            'aa_sequence':aa_sequence}

        
def calculate_interface_statistics(contacts:dict) -> dict:
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

    for interchain_id, interchain_contacts in contacts.items():
        for contact_id, contact in interchain_contacts.items():

            resix_1, resix_2 = contact_id.split('&')
            unique_residues[resix_1] = 1
            unique_residues[resix_2] = 1

            avg_plddt = np.mean(np.array(contact['plddts']))
            plddt_sum += avg_plddt
            plddt_max = max(plddt_max, avg_plddt)
            plddt_min = min(plddt_min, avg_plddt)

            afm_score = 0.5*sum(contact['af_missense'])
            num_significant_afm += 1 if afm_score > 80 else 0

            d_sum += contact['distance']
            clash_count += 1 if contact['distance'] < 1 else 0

            pae_max = max(pae_max, contact['pae'])
            pae_min = min(pae_min, contact['pae'])
            pae_sum += contact['pae']

            score = get_res_contact_score(contact)
            contact_scores = np.append(contact_scores, score)
            num_contacts += 1


    if len(contact_scores) < 1:
        contact_scores = np.append(contact_scores, 0)

    if num_contacts > 0:
        plddt_avg = round(plddt_sum/num_contacts, 1)
        pae_avg = round(pae_sum/num_contacts, 1)
        distance_avg = round(d_sum/num_contacts, 1)
    else:
        pae_min = 0
        plddt_min = 0

    data = {'num_residue_contacts':num_contacts,
            'num_residues':len(unique_residues),
            'clash_count':clash_count,
            'plddt':[plddt_min, plddt_avg, plddt_max],
            'pae':[pae_min, pae_avg, pae_max],
            'distance_avg': distance_avg,
            'contact_score_avg':round(np.mean(contact_scores), 2),
            'contact_score_median':round(np.median(contact_scores), 2),
            'contact_score_max':round(np.max(contact_scores), 2),
            'num_significant_afm_scores':num_significant_afm
    }

    return data

def summarize_interface_statistics(interfaces:dict) -> dict:
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



def analyze_complex(complex_name:str, complexes:dict):

    try:
        final_data = None
        start_time = time.time()

        interface_contacts = {}
        models = list(complexes[complex_name].keys())
        if len(models) < 3: 
            print(f'ERROR: Only {len(models)} models were supplied for {complex_name}. We require at least 3 models per complex.')
            # all_data[complex_name] = f'ERROR: Only {len(models)} models were supplied for {complex_name}. We require at least 3 models per complex.'
            return complex_name, None
        if len(models) > 3: 
            models = random.sample(models, 3)
            print(f"Found more than 3 models, randomly sampled 3 models {models}")

        aa_sequences = get_sequences(complexes[complex_name][models[0]][0])
        if (len(aa_sequences) != 2):
            print(f'ERROR: Our pipeline currently only supports 2 chain complexes. Found {len(aa_sequences)} chains.')
            # all_data[complex_name] = 'ERROR: Our pipeline currently only supports 2 chain complexes.'
            return complex_name, None

        seq_map = get_seq_to_proteome_map(aa_sequences)
        uniprot_ids_list = []
        missing_uniprot_id_chains = []
        for chain_id in seq_map:
            if seq_map[chain_id]['uniprot_id'] is None:
                missing_uniprot_id_chains.append(chain_id)
            else:
                uniprot_ids_list.append(seq_map[chain_id]['uniprot_id'])

        if len(missing_uniprot_id_chains) > 0:
            print(f'ERROR: Could not map chains {",".join(missing_uniprot_id_chains)} to unique UniProt IDs in the human proteome (as found in UniProt in December 2024). We currently only support sequences that can be mapped to canonical human protein entries.')
            return complex_name, None

        bio_data = get_bio_data_for_uniprot_ids(seq_map)

        interface_contacts = {}
        iptm_values = np.array([])
        residue_contacts_across_predictions = np.array([])
        best_interface_stats = None

        for model_num in models:
            pdb_filepath, pae_filepath = complexes[complex_name][model_num]

            data = get_data(pdb_filepath, pae_filepath, seq_map)
            contacts = data['contacts']
            interface_contacts[model_num] = contacts
            if_stats = calculate_interface_statistics(contacts)
            
            residue_contacts_across_predictions = np.append(residue_contacts_across_predictions, if_stats['num_residue_contacts']);
            iptm_values = np.append(iptm_values, data['iptm'])

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

            if best_interface_stats is None:
                best_interface_stats = if_stats
                best_interface_stats['model_num'] = model_num
            else:

                if if_stats['contacts_per_pae'] > best_interface_stats['contacts_per_pae']:
                    best_interface_stats = if_stats
                    best_interface_stats['model_num'] = model_num


        stats = summarize_interface_statistics(interface_contacts)
        stats['mean_contacts_across_predictions'] = round(np.mean(residue_contacts_across_predictions), 0)
        stats['min_contacts_across_predictions'] = np.min(residue_contacts_across_predictions)
        stats['iptm_mean'] = round(np.mean(iptm_values), 3)
        stats['iptm_min'] = round(np.min(iptm_values), 3)
        stats['iptm_max'] = round(np.max(iptm_values), 3)
        stats['models_checked'] = models
        stats['uniprot_ids'] = ":".join(uniprot_ids_list)
        stats['sequence'] = data['aa_sequence']
        stats['best_num_residue_contacts'] = best_interface_stats['num_residue_contacts']
        stats['best_if_residues'] = best_interface_stats['num_residues']
        stats['best_plddt_max'] = round(best_interface_stats['plddt'][2], 0)
        stats['best_pae_min'] = best_interface_stats['pae'][0]
        stats['best_contact_score_max'] = round(best_interface_stats['contact_score_max'], 0)
        stats['best_num_significant_afm_scores'] = best_interface_stats['num_significant_afm_scores']

        for k in bio_data:
            stats[k] = round(bio_data[k], 3)

        final_data = stats

        end_time = time.time()
        print(f'Finished: {complex_name}. Took {round(end_time - start_time, 1)} seconds')

        return complex_name, final_data
    
    except Exception as e:
        print(f'ERROR: Exception {e} occurred for complex {complex_name}. Skipping. ')
        return complex_name, None


def analyze_complexes(complexes):
    all_data = {}
    errored_complexes = []

    print(f"Splitting analysis across {multiprocessing.cpu_count()} CPUs")
    
    # Create a pool of workers equal to the number of CPU cores
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        # Map process_complex function to each item in complexes dictionary
        results = pool.starmap(analyze_complex, [(name, complexes) for name in complexes])

    # Collect results from all processes
    for complex_name, data in results:
        if data is None:
            errored_complexes.append(complex_name)
        else:
            all_data[complex_name] = data

    return all_data, errored_complexes




def main(folder_paths:list, name_filter:str, classifier, output_name:str):

    """
    Main function to run the analysis.

    :param folder_paths: List of folder paths containing data.
    :param name_filter: Filter for complex names.
    :param classifier: Random forest classifier model.
    :param output_name: Name for the output file. If not set will use a default name. 
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

            complex_name = None

            if '_unrelaxed_' in pdb_filename:
                complex_name = pdb_filename.split('_unrelaxed_')[0]

            pae_filepath = None
            for f in pae_file_paths:
                if complex_name in f and f"model_{model_num}" in f: 
                    pae_filepath = f
                    break

            if pae_filepath.split('.').pop() not in ['gz', 'xz', 'json']:
                continue

            if complex_name not in complexes:
                complexes[complex_name] = {}
            
            complexes[complex_name][model_num] = [pdb_file_path, pae_filepath]


    if len(complexes) > 0:

        print(f"Found {len(complexes)} complexes to analyze across {len(folder_paths)} folder(s).")

        data, errored_complexes = analyze_complexes(complexes)

        print(f"Errors occurred for {len(errored_complexes)} complexes.")
        if len(errored_complexes) > 0:
            final_error_csv_path = folder_path.rstrip('/') + "_errored_complexes.csv"
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
                print(f"Less than 5 contacts detected for {cname}. Assigning a SPOC score of 0.")
                cdata['spoc_score'] = 0
                continue

            print(f"Running SPOC on {cname}")

            data_for_rf = {}
            for f in classifier.feature_names_in_:
                data_for_rf[f] = cdata[f]

            data_for_rf = pd.DataFrame(data_for_rf, index=[0])
            spoc_score = round(float(classifier.predict_proba(data_for_rf)[:, 1][0]), 3)
            cdata['spoc_score'] = spoc_score
            print(f"SPOC score = {spoc_score}")


        final_datalist = []

        for cname in data:
            data[cname]['complex_name'] = cname
            final_datalist.append(data[cname])
        
        columns = list(final_datalist[0].keys())
        columns.remove('complex_name')
        columns.remove('spoc_score')
        columns = ['complex_name', 'spoc_score'] + columns

        final_df = pd.DataFrame(final_datalist)
        final_df = final_df[columns]

        folder_name = folder_path.rstrip('/').split('/').pop()
        final_csv_path = folder_name + "_SPOC_analysis.csv"
        if output_name: 
            if not output_name.endswith('.csv'):
                output_name += '.csv'
            final_csv_path = output_name
  
        final_df.sort_values(by=['spoc_score'], ascending=False, inplace=True)
        final_df.to_csv(final_csv_path, index=None)
        print("")
        print("---------------------")
        print(f"Wrote analysis data out to {final_csv_path}")

    else:
        print('ERROR: No complexes to analyze found')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run the SPOC classifier on data in the specified folder which should contain PDB and JSON PAE files.")
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
    # Parse arguments
    args = parser.parse_args()
    if len(args.folder_paths) < 0:
        print("No folders to analyze were provided. Script is terminating.")
        exit()


    start_message = """
-------------------------------------------------------------------------------------------------------------------------
Thank you for using our SPOC analysis script.
Please note that SPOC is designed to work with AlphaFold2 multimer model data. 
If you are using a different model for your predictions, the SPOC results may not be accurate.

In addition, SPOC currently only works on HUMAN protein pairs predicted in at least 3 AF-M models. 
Complexes with less than 3 predictions or those that do not contain exactly 2 human protein chains will be skipped.
-------------------------------------------------------------------------------------------------------------------------  
"""

    print(start_message)

    print("Loading external databases")
    print("-----------------------------------")
    initialize_databases()
    print("-----------------------------------")
    print("Finished loading external databases")

    # Load SPOC random forest classifier model
    script_path = os.path.realpath(__file__)
    spoc_rf_path = os.path.dirname(script_path) + "/models/SPOC_rf_params.joblib"
    print(f"Loading SPOC Random Forest model from {spoc_rf_path}")
    rf_classifier = load(spoc_rf_path)
    print("Finished loading SPOC Random Forest model")
    print("-----------------------------------")
    print("")

    main_start_time = time.time()
    main(args.folder_paths, args.name_filter, rf_classifier, args.output)
    main_end_time = time.time()
    print(f'SPOC Analysis complete. Took a total of {round(main_end_time - main_start_time, 1)} seconds')
    print("-----------------------------------")
    print("")
    print(f'If you use this code or its output as part of a publication, please cite the AlphaFold2 multimer manuscript, the ColabFold manuscript and the SPOC manuscript.')
