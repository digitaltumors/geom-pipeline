###
#
#  Code to fetch the CIDs for the compounds in the GEOM dataset
#  and compute both the average and most probable WHIM descriptor
#  of each compound
#
#   Author: James Bannon
#   Github: @jbannon
####


import json
import os

# import flags
import os.path
import pickle
import time
from collections import defaultdict
from typing import Dict

import numpy as np
import pandas as pd
import pubchempy as pcp
import tqdm
from damply import dirs
from rdkit.Chem.rdMolDescriptors import CalcWHIM

RENAME_FIELDS = {
	'totalconfs': 'Total Conformations',
	'uniqueconfs': 'Unique Conformations',
	'temperature': 'Temperature',
	'ensembleneergy': 'Ensemble Energy',
	'ensemblentropy': 'Ensemble Entropy',
	'poplowestpct': 'Population Lowest Percentage',
	'ensemblefreeenergy': 'Ensemble Free Energy',
}


def fetch_molecule_names(smiles: str):
	"""
	Function to search for the PubChem CID of a molecule given in CID format
	"""
	try:
		search_results = pcp.get_compounds(identifier=smiles, namespace='smiles')
		cid = None if search_results[0].cid is None else search_results[0].cid
		inchi_key = (
			None if search_results[0].cid is None else search_results[0].inchikey
		)
	except:
		# hack to get around pubchem's throttle
		time.sleep(3.0)
		search_results = pcp.get_compounds(identifier=smiles, namespace='smiles')
		cid = None if search_results[0].cid is None else search_results[0].cid
		inchi_key = (
			None if search_results[0].cid is None else search_results[0].inchikey
		)
	return cid, inchi_key


def update_colData(colData: Dict, smiles: str, mol: Dict) -> int:

	cid, inchi = fetch_molecule_names(smiles)
	colData['PubChemCID'].append(cid)
	colData['SMILES'].append(smiles)
	colData['INCHI-KEY'].append(inchi)

	for field in RENAME_FIELDS:
		if field in mol:  # if this value has a measurement for the molecule
			val = mol[field]
		else:
			val = 'NA'

		colData[RENAME_FIELDS[field]].append(val)

	return cid


def main():

	colData = defaultdict(list)
	wavg_whims, avg_whims, hp_whims = {}, {}, {}

	processed_mols = []

	base_dir = dirs.RAWDATA / 'rdkit_folder'

	subset_names = (['drugs', 'qm9'],)
	file_names = (['summary_drugs.json', 'summary_qm9.json'],)

	for subset_name, file_name in zip(subset_names, file_names):
		with open(base_dir / 'summary_drugs.json', 'r') as f:
			molecule_info = json.load(f)

		for mol_smiles in (pbar := tqdm.tqdm(molecule_info.keys())):
			pbar.set_description(f'Working on {mol_smiles}')

			mol_metadata = molecule_info[mol_smiles]

			if not 'pickle_path' in mol_metadata.keys():
				continue

			pickle_file = base_dir / mol_metadata['pickle_path']

			if os.path.isfile(pickle_file):  # CAN DROP FOR FULL
				with open(pickle_file, 'rb') as input:
					mol_data = pickle.load(input)

				cid = update_colData(colData, mol_smiles, mol_data)
				conformations = mol_data['conformers']

				valid_idx = [
					i
					for i, conf in enumerate(conformations)
					if 'boltzmannweight' in conf
				]

				whims = np.array(
					[CalcWHIM(conformations[i]['rd_mol']) for i in valid_idx]
				)
				weights = np.array(
					[conformations[i]['boltzmannweight'] for i in valid_idx]
				)

				highest_prob_whim = whims[np.argmax(weights)]

				weighted_avg_whim = (weights.reshape(-1, 1) * whims).sum(0)

				avg_whim = np.mean(whims, axis=0)

				hp_whims[cid] = highest_prob_whim
				wavg_whims[cid] = weighted_avg_whim
				avg_whims[cid] = avg_whim

	# write ColData
	colData = pd.DataFrame(colData)
	colData.to_csv(dirs.RESULTS / 'colData.csv')

	# Write weighted average whims
	wavg_whims = pd.DataFrame(wavg_whims)
	wavg_whims = wavg_whims.to_csv(dirs.RESULTS / 'WHIMS_wavg.csv')

	# write average whims
	avg_whims = pd.DataFrame(avg_whims)
	avg_whims.to_csv(dirs.RESULTS / 'WHIMS_avg.csv')

	# wrrite highest prob whims
	hp_whims = pd.DataFrame(hp_whims)
	hp_whims.to_csv(dirs.RESULTS / 'WHIM_hp.csv')


if __name__ == '__main__':
	main()
