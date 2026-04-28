# Data Sources

## GEOM Input Data

- **Name**: GEOM: energy-annotated molecular conformations
- **Primary format used here**: `rdkit_folder.tar.gz` downloaded into
  `data/rawdata/` and extracted locally
- **Upstream project**: [learningmatter-mit/geom](https://github.com/learningmatter-mit/geom)
- **Citation**: Axelrod S, Gómez-Bombarelli R. *GEOM, energy-annotated
  molecular conformations for property prediction and molecular generation*.
  Scientific Data 9, 185 (2022).
- **License**: CC0 1.0 Universal, according to the data drop README bundled with
  the Dataverse data bundle.

## Download Sources

- Harvard Dataverse dataset page:
  `https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/JNGTDF`
- RDKit archive:
  `https://dataverse.harvard.edu/api/access/datafile/4327252`
- MoleculeNet archive:
  `https://dataverse.harvard.edu/api/access/datafile/5858506`

These direct links remain in `config/pipeline.yaml` and are the default source
for raw data.

## Expected Local Layout

The workflow downloads into `data/rawdata/` and extracts the RDKit archive to:

```text
data/rawdata/rdkit_folder
```

Required files for the curated `drugs` and `qm9` path:

- `summary_drugs.json`
- `summary_qm9.json`
- `drugs/*.pickle`
- `qm9/*.pickle`

The optional `moleculenet` subset currently downloads the upstream
`molecule_net.tar.gz` archive into `data/rawdata/` but does not yet feed the
main curation tables.

## Large-File Considerations

The upstream GEOM files are large:

- `rdkit_folder.tar.gz` at about 50.1 GB
- `drugs_crude.msgpack.tar.gz` at about 42.7 GB
- `molecule_net.tar.gz` at about 2.2 GB
- `rdkit_folder/summary_drugs.json` at about 168 MB
- `rdkit_folder/summary_qm9.json` at about 48 MB

These sizes matter operationally:

- make sure `data/rawdata/` has sufficient free space before a real run
- all processed tables are written to `data/procdata/`

## AnnotationDB Enrichment

Compound enrichment uses the configurable AnnotationDB API base URL:

- config key: `annotationdb_api`
- default: `https://v2annotationdb.bhklab.ca`
- workflow endpoint: `<annotationdb_api>/compound/all`

The workflow caches the bulk response as
`data/rawdata/metadata/all_adb_compounds.csv` and joins it locally by
`Metadata_InChIKey`.

The cached table keeps only the normalized join and display fields:
`inchikey`, `cid`, `name`, and `smiles`. Processed metadata renames these to the
pipeline's internal `AnnotationDB_*` fields before MAE construction.

## Generated Data Products

The workflow generates:

- `data/procdata/metadata/geom_manifest.tsv`
- `data/procdata/metadata/geom_manifest_summary.tsv`
- `data/procdata/metadata/compound_master.tsv`
- `data/procdata/metadata/drug_metadata_raw.tsv`
- `data/procdata/metadata/colData_raw.csv`
- `data/procdata/conformers/conformer_metadata.tsv`
- `data/procdata/assays/WHIM_hp.csv`
- `data/procdata/assays/WHIMS_avg.csv`
- `data/procdata/assays/WHIMS_wavg.csv`

Final MAE-derived exports use public dot-style headers.
