import csv
import os
from pathlib import Path

import ijson

SUMMARY_FIELD_MAP = {
	'charge': 'Charge',
	'totalconfs': 'Total_Conformers',
	'uniqueconfs': 'Unique_Conformers',
	'temperature': 'Temperature',
	'lowestenergy': 'Lowest_Energy',
	'ensembleenergy': 'Ensemble_Energy',
	'ensembleentropy': 'Ensemble_Entropy',
	'ensemblefreeenergy': 'Ensemble_Free_Energy',
	'poplowestpct': 'Lowest_Energy_Population_Percent',
}


rdkit_root = Path(snakemake.params.rdkit_root)
max_compounds_per_subset = snakemake.params.max_compounds_per_subset
summary_paths = [Path(path) for path in snakemake.input.summary_jsons]
manifest_path = Path(snakemake.output.manifest_tsv)
summary_output_path = Path(snakemake.output.summary_tsv)

manifest_rows = []
summary_rows = []

for summary_path in summary_paths:
	subset = summary_path.stem.removeprefix('summary_')
	existing_pickles = set()
	if max_compounds_per_subset is None and (rdkit_root / subset).is_dir():
		subset_dir = rdkit_root / subset
		with os.scandir(subset_dir) as entries:
			for entry in entries:
				if entry.is_file() and entry.name.endswith('.pickle'):
					existing_pickles.add(f'{subset}/{entry.name}')

	rows = []
	scanned_rows = 0
	rows_with_pickle_path = 0
	rows_with_existing_pickle = 0
	summary_scan_complete = True

	with summary_path.open('rb') as handle:
		for source_index, (source_smiles, record) in enumerate(
			ijson.kvitems(handle, '', use_float=True)
		):
			scanned_rows += 1
			pickle_relpath = record.get('pickle_path', '')
			pickle_abspath = rdkit_root / pickle_relpath if pickle_relpath else None
			if max_compounds_per_subset is None:
				pickle_exists = pickle_relpath in existing_pickles
			else:
				pickle_exists = bool(pickle_abspath and pickle_abspath.is_file())

			if pickle_relpath:
				rows_with_pickle_path += 1
			if pickle_exists:
				rows_with_existing_pickle += 1

			row = {
				'Source_Subset': subset,
				'Source_Order': source_index,
				'Source_SMILES': source_smiles,
				'Pickle_Path': pickle_relpath,
			}
			for raw_field, output_field in SUMMARY_FIELD_MAP.items():
				row[output_field] = record.get(raw_field, '')

			if max_compounds_per_subset is None:
				rows.append(row)
				continue

			if pickle_exists and len(rows) < int(max_compounds_per_subset):
				rows.append(row)

			if len(rows) >= int(max_compounds_per_subset):
				summary_scan_complete = False
				break

	manifest_rows.extend(rows)
	summary_rows.append(
		{
			'Source_Subset': subset,
			'Summary_JSON_Path': str(summary_path),
			'Summary_Scan_Complete': str(summary_scan_complete).lower(),
			'Scanned_Rows': scanned_rows,
			'Rows_With_Pickle_Path': rows_with_pickle_path,
			'Rows_With_Existing_Pickle': rows_with_existing_pickle,
			'Rows_Missing_Pickle': rows_with_pickle_path - rows_with_existing_pickle,
			'Selected_Rows': len(rows),
			'Max_Compounds_Per_Subset': (
				'' if max_compounds_per_subset is None else max_compounds_per_subset
			),
		}
	)

manifest_columns = [
	'Source_Subset',
	'Source_Order',
	'Source_SMILES',
	'Pickle_Path',
	*SUMMARY_FIELD_MAP.values(),
]

summary_columns = [
	'Source_Subset',
	'Summary_JSON_Path',
	'Summary_Scan_Complete',
	'Scanned_Rows',
	'Rows_With_Pickle_Path',
	'Rows_With_Existing_Pickle',
	'Rows_Missing_Pickle',
	'Selected_Rows',
	'Max_Compounds_Per_Subset',
]

manifest_path.parent.mkdir(parents=True, exist_ok=True)
with manifest_path.open('w', newline='', encoding='utf-8') as handle:
	writer = csv.DictWriter(handle, fieldnames=manifest_columns, delimiter='\t')
	writer.writeheader()
	writer.writerows(manifest_rows)

with summary_output_path.open('w', newline='', encoding='utf-8') as handle:
	writer = csv.DictWriter(handle, fieldnames=summary_columns, delimiter='\t')
	writer.writeheader()
	writer.writerows(summary_rows)

print(
	'[build_manifest] '
	f'subsets={",".join(sorted({row["Source_Subset"] for row in manifest_rows}))} '
	f'rows={len(manifest_rows)} '
	f'manifest={manifest_path} '
	f'summary={summary_output_path}',
	flush=True,
)
