# Developer Notes

## Environment Decisions

- The default runtime is pinned to Python 3.12.
- This is intentional: current conda `snakemake` builds on `osx-arm64` solve cleanly with Python 3.12.
- Core scientific packages and Snakemake are kept in Pixi conda dependencies.

## Workflow Decisions

- The root `Snakefile` delegates to `workflow/Snakefile`.
- Raw GEOM inputs are downloaded into `data/rawdata/` from Harvard Dataverse.
- The HDD-facing default curates only the `drugs` subset from the downloaded `rdkit_folder.tar.gz` archive.
- `qm9` and `moleculenet` remain config-valid side-product subsets, but they are not part of the default HDD-facing output.

## Identifier Decisions

- The compound-level MAE key is `GEOM.Source.SMILES`, using GEOM's own source SMILES key.
- The pipeline errors if `GEOM.Source.SMILES` is missing or duplicated instead of generating a replacement identifier.
- Public columns with HDD-shared names are intended to be directly joinable to the base HDD; source-specific fields use source-specific prefixes.

## Metadata Decisions

- RDKit derives `Metadata_SMILES`, `Metadata_InChI`, and `Metadata_InChIKey` locally from stored conformer molecules.
- AnnotationDB enrichment is intentionally minimal: PubChem CID, AnnotationDB name, AnnotationDB SMILES, and a match flag.
- The AnnotationDB `/compound/all` response is cached once at `data/rawdata/metadata/all_adb_compounds.csv` and downstream joins read the cache.
- Pickle paths are kept only as internal manifest fields needed to locate GEOM files during processing. They are dropped from compound metadata, conformer metadata, MAE `colData`, MAE metadata, and public table exports.
- Public-facing tables drop bulky bookkeeping and raw JSON payload columns.

## Output Decisions

- The three assays are WHIM summaries, not a top-three conformer selection: `WHIM_hp`, `WHIM_avg`, and `WHIM_wavg`.
- `conformer_metadata.tsv` preserves per-conformer GEOM metadata whether or not the conformer is valid for WHIM aggregation.
- Output headers are normalized to public dot-style names only at MAE/export boundaries; internal intermediates may retain source-style names for traceability.
