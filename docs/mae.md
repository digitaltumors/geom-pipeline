# MAE Structure

The final object is written to `data/results/geom_MultiAssayExperiment.rds`. The HDD-facing default is compound-level GEOM drugs.

| Component | Structure |
| --- | --- |
| Assays | `WHIM_hp`, `WHIM_avg`, and `WHIM_wavg`. |
| Assay columns | `GEOM.Source.SMILES`, the GEOM source molecule key. |
| Assay rows | WHIM descriptor features named `WHIM_001`, `WHIM_002`, and so on. |
| sampleMap | One row per assay/compound pair with `assay`, `primary`, and `colname`; both `primary` and `colname` are `GEOM.Source.SMILES`. |
| colData | Compound-level metadata. Each row is a GEOM drug compound and is keyed by `GEOM.Source.SMILES`. |
| rowData | Feature metadata with `Feature.ID`, `Feature.Index`, and `Assay.Name`. |

## metadata(mae)

| Object | Source | Purpose |
| --- | --- | --- |
| `Pipeline` | `config/pipeline.yaml` and Snakefile params | Named list containing `ID`, `Version`, and parsed run `Config`. |
| `Drug.Metadata` | `data/procdata/metadata/drug_metadata_raw.tsv` plus AnnotationDB cache | Harmonized compound metadata keyed by `GEOM.Source.SMILES`. |
