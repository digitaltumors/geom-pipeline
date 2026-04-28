# GEOM Pipeline

This pipeline curates GEOM molecular conformer data into WHIM descriptor assays, processed metadata, a `MultiAssayExperiment`, and an MAE-derived tabular archive. Sources include Harvard Dataverse GEOM archives and AnnotationDB `/compound/all` for PubChem CID enrichment.

Run from the repository root:

```bash
pixi run snakemake --cores <n>
```

Edit `config/pipeline.yaml` to change GEOM subsets, subset limits, source URLs, paths, and the AnnotationDB API base URL.
