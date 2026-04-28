import csv
import json
import pickle
from contextlib import ExitStack
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import CalcWHIM

# Columns to keep (and rename)
SUMMARY_FIELD_MAP = {
	'Charge': 'Charge',
	'Total_Conformers': 'Total_Conformers',
	'Unique_Conformers': 'Unique_Conformers',
	'Temperature': 'Temperature',
	'Lowest_Energy': 'Lowest_Energy',
	'Ensemble_Energy': 'Ensemble_Energy',
	'Ensemble_Entropy': 'Ensemble_Entropy',
	'Ensemble_Free_Energy': 'Ensemble_Free_Energy',
	'Lowest_Energy_Population_Percent': 'Lowest_Energy_Population_Percent',
}


def read_manifest(path: Path) -> list[dict[str, str]]:
	with path.open(newline='', encoding='utf-8') as handle:
		return list(csv.DictReader(handle, delimiter='\t'))


def json_dumps(value: object) -> str:
	# handles "null" values with empty string
	if value in (None, ''):
		return ''
	return json.dumps(value, sort_keys=True)


def safe_mol_to_inchi(mol: Mol) -> str:
	try:
		return Chem.MolToInchi(mol)
	except Exception:
		return ''


def safe_mol_to_inchikey(mol: Mol) -> str:
	try:
		return inchi.MolToInchiKey(mol)
	except Exception:
		return ''


def load_pickle_payload(pickle_path: Path) -> dict[str, Any]:
	with pickle_path.open('rb') as handle:
		return pickle.load(handle)


def resolve_pickle_path(rdkit_root: Path, row: dict[str, str]) -> Path | None:
	pickle_path = row['Pickle_Path']
	if not pickle_path:
		return None
	return rdkit_root / pickle_path


# Get the WHIM dim dynamically just in case it varies
def infer_whim_dimension(
	manifest_rows: list[dict[str, str]], rdkit_root: Path
) -> tuple[int, dict[Path, dict[str, Any]]]:
	for row in manifest_rows:
		pickle_path = resolve_pickle_path(rdkit_root, row)
		if pickle_path is None or not pickle_path.is_file():
			continue
		payload = load_pickle_payload(pickle_path)
		for conformer in payload.get('conformers', []):
			mol = conformer.get('rd_mol')
			if mol is None or 'boltzmannweight' not in conformer:
				continue
			# also returns the cached payload for reuse!
			return len(CalcWHIM(mol)), {pickle_path: payload}
	message = (
		'Could not find any conformer with both rd_mol and boltzmannweight; '
		'cannot infer WHIM descriptor length.'
	)
	raise ValueError(message)


# qol utility
def blank_whim_row(compound_id: str, columns: list[str]) -> dict[str, str]:
	return {'Compound_ID': compound_id, **{column: '' for column in columns}}


# Init Snakemake variables
manifest_path = Path(snakemake.input.manifest_tsv)
rdkit_root = Path(snakemake.params.rdkit_root)
output_compound_master = Path(snakemake.output.compound_master)
output_conformer_metadata = Path(snakemake.output.conformer_metadata)
output_whim_hp = Path(snakemake.output.whim_hp)
output_whims_avg = Path(snakemake.output.whims_avg)
output_whims_wavg = Path(snakemake.output.whims_wavg)

manifest_rows = read_manifest(manifest_path)
if not manifest_rows:
	message = 'Manifest contains no rows.'
	raise ValueError(message)

whim_dimension, prefetched_payloads = infer_whim_dimension(manifest_rows, rdkit_root)
whim_columns = [f'WHIM_{index:03d}' for index in range(1, whim_dimension + 1)]

compound_columns = [
	'Compound_ID',
	'Source_Subset',
	'Source_SMILES',
	'Metadata_SMILES',
	'Canonical_SMILES',
	'Metadata_InChI',
	'Metadata_InChIKey',
	*SUMMARY_FIELD_MAP.values(),
	'Conformer_Record_Count',
	'Conformers_With_Boltzmann_Weight',
	'Conformers_Without_Boltzmann_Weight',
]

conformer_columns = [
	'Compound_ID',
	'Source_Subset',
	'Source_SMILES',
	'Metadata_InChIKey',
	'Conformer_Index',
	'Is_Valid_For_WHIM',
	'GEOM_ID',
	'GEOM_Set',
	'Degeneracy',
	'Total_Energy',
	'Relative_Energy',
	'Boltzmann_Weight',
	'Conformer_Weights_JSON',
	'RDKit_Atom_Count',
	'RDKit_Bond_Count',
	'Conformer_Record_JSON',
]

whim_output_columns = ['Compound_ID', *whim_columns]

output_compound_master.parent.mkdir(parents=True, exist_ok=True)
output_conformer_metadata.parent.mkdir(parents=True, exist_ok=True)
output_whim_hp.parent.mkdir(parents=True, exist_ok=True)
output_whims_avg.parent.mkdir(parents=True, exist_ok=True)
output_whims_wavg.parent.mkdir(parents=True, exist_ok=True)

metrics = {
	'manifest_rows': len(manifest_rows),
	'processed_compounds': 0,
	'skipped_missing_pickles': 0,
	'compounds_without_valid_whims': 0,
	'written_conformer_rows': 0,
}

seen_compound_ids = set()

# Use exit stack to handle multiple output files in parallel safely
with ExitStack() as stack:
	compound_handle = stack.enter_context(
		output_compound_master.open('w', newline='', encoding='utf-8')
	)
	conformer_handle = stack.enter_context(
		output_conformer_metadata.open('w', newline='', encoding='utf-8')
	)
	whim_hp_handle = stack.enter_context(
		output_whim_hp.open('w', newline='', encoding='utf-8')
	)
	whims_avg_handle = stack.enter_context(
		output_whims_avg.open('w', newline='', encoding='utf-8')
	)
	whims_wavg_handle = stack.enter_context(
		output_whims_wavg.open('w', newline='', encoding='utf-8')
	)

	compound_writer = csv.DictWriter(
		compound_handle, fieldnames=compound_columns, delimiter='\t'
	)
	conformer_writer = csv.DictWriter(
		conformer_handle, fieldnames=conformer_columns, delimiter='\t'
	)
	whim_hp_writer = csv.DictWriter(whim_hp_handle, fieldnames=whim_output_columns)
	whims_avg_writer = csv.DictWriter(whims_avg_handle, fieldnames=whim_output_columns)
	whims_wavg_writer = csv.DictWriter(
		whims_wavg_handle, fieldnames=whim_output_columns
	)

	compound_writer.writeheader()
	conformer_writer.writeheader()
	whim_hp_writer.writeheader()
	whims_avg_writer.writeheader()
	whims_wavg_writer.writeheader()

	for row in manifest_rows:
		pickle_path = resolve_pickle_path(rdkit_root, row)
		if pickle_path is None or not pickle_path.is_file():
			metrics['skipped_missing_pickles'] += 1
			continue

		payload = prefetched_payloads.pop(pickle_path, None)
		if payload is None:
			payload = load_pickle_payload(pickle_path)

		conformers = payload.get('conformers', [])
		representative_mol = next(
			(
				conf.get('rd_mol')
				for conf in conformers
				if conf.get('rd_mol') is not None
			),
			None,
		)

		canonical_smiles = (
			Chem.MolToSmiles(representative_mol, canonical=True)
			if representative_mol is not None
			else row['Source_SMILES']
		)
		metadata_smiles = canonical_smiles or row['Source_SMILES']
		metadata_inchi = (
			safe_mol_to_inchi(representative_mol) if representative_mol else ''
		)
		metadata_inchikey = (
			safe_mol_to_inchikey(representative_mol) if representative_mol else ''
		)

		compound_id = row['Source_SMILES']
		if compound_id in seen_compound_ids:
			message = (
				'Duplicate GEOM Source_SMILES cannot be used as Compound_ID: '
				f'{compound_id}'
			)
			raise ValueError(message)
		seen_compound_ids.add(compound_id)

		whims = []
		weights = []
		n_with_boltzmann = 0
		n_without_boltzmann = 0

		for conformer_index, conformer in enumerate(conformers):
			conformer_copy = {
				key: value for key, value in conformer.items() if key != 'rd_mol'
			}
			mol = conformer.get('rd_mol')

			# Check if the conformer has a Boltzmann weight; only those are valid for WHIM calculations
			has_boltzmann = 'boltzmannweight' in conformer
			is_valid_for_whim = mol is not None and has_boltzmann

			if has_boltzmann:
				n_with_boltzmann += 1
			else:
				n_without_boltzmann += 1

			if is_valid_for_whim:
				whim_values = np.asarray(CalcWHIM(mol), dtype=float)
				whims.append(whim_values)
				weights.append(float(conformer['boltzmannweight']))

			conformer_writer.writerow(
				{
					'Compound_ID': compound_id,
					'Source_Subset': row['Source_Subset'],
					'Source_SMILES': row['Source_SMILES'],
					'Metadata_InChIKey': metadata_inchikey,
					'Conformer_Index': conformer_index,
					'Is_Valid_For_WHIM': str(is_valid_for_whim).lower(),
					'GEOM_ID': conformer.get('geom_id', ''),
					'GEOM_Set': conformer.get('set', ''),
					'Degeneracy': conformer.get('degeneracy', ''),
					'Total_Energy': conformer.get('totalenergy', ''),
					'Relative_Energy': conformer.get('relativeenergy', ''),
					'Boltzmann_Weight': conformer.get('boltzmannweight', ''),
					'Conformer_Weights_JSON': json_dumps(
						conformer.get('conformerweights')
					),
					'RDKit_Atom_Count': mol.GetNumAtoms() if mol is not None else '',
					'RDKit_Bond_Count': mol.GetNumBonds() if mol is not None else '',
					'Conformer_Record_JSON': json_dumps(conformer_copy),
				}
			)
			metrics['written_conformer_rows'] += 1

		compound_writer.writerow(
			{
				'Compound_ID': compound_id,
				'Source_Subset': row['Source_Subset'],
				'Source_SMILES': row['Source_SMILES'],
				'Metadata_SMILES': metadata_smiles,
				'Canonical_SMILES': canonical_smiles,
				'Metadata_InChI': metadata_inchi,
				'Metadata_InChIKey': metadata_inchikey,
				**{
					output_field: row.get(input_field, '')
					for input_field, output_field in SUMMARY_FIELD_MAP.items()
				},
				'Conformer_Record_Count': len(conformers),
				'Conformers_With_Boltzmann_Weight': n_with_boltzmann,
				'Conformers_Without_Boltzmann_Weight': n_without_boltzmann,
			}
		)

		if not whims:
			blank_row = blank_whim_row(compound_id, whim_columns)
			whim_hp_writer.writerow(blank_row)
			whims_avg_writer.writerow(blank_row)
			whims_wavg_writer.writerow(blank_row)
			metrics['compounds_without_valid_whims'] += 1
			metrics['processed_compounds'] += 1
			continue

		whims_array = np.vstack(whims)
		weights_array = np.asarray(weights, dtype=float)
		highest_prob_whim = whims_array[np.argmax(weights_array)]
		avg_whim = np.mean(whims_array, axis=0)
		weighted_avg_whim = np.average(whims_array, axis=0, weights=weights_array)

		whim_hp_writer.writerow(
			{
				'Compound_ID': compound_id,
				**dict(zip(whim_columns, highest_prob_whim.tolist(), strict=True)),
			}
		)
		whims_avg_writer.writerow(
			{
				'Compound_ID': compound_id,
				**dict(zip(whim_columns, avg_whim.tolist(), strict=True)),
			}
		)
		whims_wavg_writer.writerow(
			{
				'Compound_ID': compound_id,
				**dict(zip(whim_columns, weighted_avg_whim.tolist(), strict=True)),
			}
		)
		metrics['processed_compounds'] += 1

print(
	'[process_geom] '
	f'manifest_rows={metrics["manifest_rows"]} '
	f'processed_compounds={metrics["processed_compounds"]} '
	f'skipped_missing_pickles={metrics["skipped_missing_pickles"]} '
	f'compounds_without_valid_whims={metrics["compounds_without_valid_whims"]} '
	f'written_conformer_rows={metrics["written_conformer_rows"]}',
	flush=True,
)
