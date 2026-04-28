# colData Columns

Final `colData` is built in `workflow/scripts/build_mae.R` from processed GEOM compound metadata and the cached AnnotationDB table. Columns are exported from the MAE to `data/results/geom_tables/colData.tsv`.

| Column | Type | Description | Computed from / origin |
| --- | --- | --- | --- |
| `GEOM.Source.SMILES` | character | Compound-level MAE primary key, using GEOM's source molecule key. | GEOM summary JSON molecule key copied through `Source_SMILES`; duplicates error before MAE construction. |
| `Pubchem.CID` | integer | PubChem compound identifier used to join back to the base HDD. | AnnotationDB `cid`, joined through RDKit-derived InChIKey. |
| `InChIKey` | character | Compound InChIKey used for metadata joining. | RDKit-derived `Metadata_InChIKey` from the representative GEOM conformer molecule. |
| `GEOM.RDKit.SMILES` | character | RDKit-derived SMILES retained for comparison with the source key. | `Metadata_SMILES`, falling back to source SMILES during processing if needed. |
| `In.AnnotationDB` | logical | Whether the compound has an AnnotationDB match or PubChem CID. | Processed `In_AnnotationDB` flag OR non-missing `Pubchem.CID`. |
| `AnnotationDB.Name` | character | AnnotationDB compound name. | AnnotationDB `name` from the cached `/compound/all` response. |
| `AnnotationDB.SMILES` | character | AnnotationDB SMILES string. | AnnotationDB `smiles` from the cached `/compound/all` response. |
| `GEOM.Source.Subset` | character | GEOM subset used for the record. | Manifest `Source_Subset`; default HDD-facing output is `drugs`. |
| `GEOM.Canonical.SMILES` | character | Canonical SMILES derived from the representative conformer molecule. | RDKit `Chem.MolToSmiles(..., canonical=TRUE)` during GEOM processing. |
| `GEOM.Charge` | numeric | Molecular charge reported by GEOM. | GEOM summary field `Charge`, numeric cast in MAE construction. |
| `GEOM.Total.Conformers` | integer | Total conformer count reported by GEOM. | GEOM summary field `Total_Conformers`, integer cast in MAE construction. |
| `GEOM.Unique.Conformers` | integer | Unique conformer count reported by GEOM. | GEOM summary field `Unique_Conformers`, integer cast in MAE construction. |
| `GEOM.Lowest.Energy` | numeric | Lowest conformer energy reported by GEOM. | GEOM summary field `Lowest_Energy`, numeric cast in MAE construction. |
| `GEOM.Ensemble.Energy` | numeric | Ensemble energy reported by GEOM. | GEOM summary field `Ensemble_Energy`, numeric cast in MAE construction. |
