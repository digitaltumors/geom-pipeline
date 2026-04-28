import csv
from pathlib import Path

REQUIRED_COMPOUND_MASTER_COLUMNS = [
	'Compound_ID',
	'Source_Subset',
	'Source_SMILES',
	'Metadata_SMILES',
	'Canonical_SMILES',
	'Metadata_InChI',
	'Metadata_InChIKey',
	'Charge',
	'Total_Conformers',
	'Unique_Conformers',
	'Temperature',
	'Lowest_Energy',
	'Ensemble_Energy',
	'Ensemble_Entropy',
	'Ensemble_Free_Energy',
	'Lowest_Energy_Population_Percent',
	'Conformer_Record_Count',
	'Conformers_With_Boltzmann_Weight',
	'Conformers_Without_Boltzmann_Weight',
]

REQUIRED_ANNOTATIONDB_COLUMNS = [
	'inchikey',
	'cid',
	'name',
	'smiles',
]

COMPOUND_METADATA_COLUMNS = [
	*REQUIRED_COMPOUND_MASTER_COLUMNS,
	'AnnotationDB_CID',
	'AnnotationDB_Name',
	'AnnotationDB_SMILES',
	'In_AnnotationDB',
]

COLDATA_COLUMNS = [
	column for column in COMPOUND_METADATA_COLUMNS if column != 'Metadata_SMILES'
]


def read_tsv(path: Path) -> list[dict[str, str]]:
	with path.open(newline='', encoding='utf-8') as handle:
		return list(csv.DictReader(handle, delimiter='\t'))


def read_csv(path: Path) -> list[dict[str, str]]:
	with path.open(newline='', encoding='utf-8') as handle:
		return list(csv.DictReader(handle))


def validate_columns(
	rows: list[dict[str, str]], required_columns: list[str], label: str
) -> None:
	if not rows:
		message = f'{label} contains no rows'
		raise ValueError(message)

	missing = [column for column in required_columns if column not in rows[0]]
	if missing:
		message = f'{label} is missing required columns: {", ".join(missing)}'
		raise ValueError(message)


compound_master_path = Path(snakemake.input.compound_master)
annotationdb_path = Path(snakemake.input.annotationdb_cache)
compound_metadata_path = Path(snakemake.output.compound_metadata)
coldata_path = Path(snakemake.output.coldata)

compound_master = read_tsv(compound_master_path)
annotationdb_rows = read_csv(annotationdb_path)

validate_columns(
	compound_master, REQUIRED_COMPOUND_MASTER_COLUMNS, 'compound_master.tsv'
)
validate_columns(
	annotationdb_rows, REQUIRED_ANNOTATIONDB_COLUMNS, 'annotationdb_compounds.tsv'
)

annotationdb_by_inchikey = {}
for row in annotationdb_rows:
	inchikey = row['inchikey']
	if inchikey and inchikey not in annotationdb_by_inchikey:
		annotationdb_by_inchikey[inchikey] = row

compound_metadata = []
annotationdb_matches = 0
for row in compound_master:
	annotation_row = annotationdb_by_inchikey.get(row['Metadata_InChIKey'], {})
	in_annotationdb = annotation_row.get('cid', '') != ''
	if in_annotationdb:
		annotationdb_matches += 1

	output_row = {column: row[column] for column in REQUIRED_COMPOUND_MASTER_COLUMNS}
	output_row['AnnotationDB_CID'] = annotation_row.get('cid', '')
	output_row['AnnotationDB_Name'] = annotation_row.get('name', '')
	output_row['AnnotationDB_SMILES'] = annotation_row.get('smiles', '')
	output_row['In_AnnotationDB'] = str(in_annotationdb).lower()
	compound_metadata.append(output_row)

compound_metadata.sort(key=lambda row: row['Compound_ID'])

compound_metadata_path.parent.mkdir(parents=True, exist_ok=True)
with compound_metadata_path.open('w', newline='', encoding='utf-8') as handle:
	writer = csv.DictWriter(
		handle, fieldnames=COMPOUND_METADATA_COLUMNS, delimiter='\t'
	)
	writer.writeheader()
	writer.writerows(compound_metadata)

with coldata_path.open('w', newline='', encoding='utf-8') as handle:
	writer = csv.DictWriter(handle, fieldnames=COLDATA_COLUMNS)
	writer.writeheader()
	for row in compound_metadata:
		writer.writerow({column: row[column] for column in COLDATA_COLUMNS})

print(
	'[build_compound_metadata] '
	f'compound_rows={len(compound_metadata)} '
	f'annotationdb_matches={annotationdb_matches} '
	f'annotationdb_cache={annotationdb_path} '
	f'compound_metadata={compound_metadata_path} '
	f'coldata={coldata_path}',
	flush=True,
)
