# procdata

Processed GEOM intermediates live here. The pipeline assumes generated metadata under `procdata/metadata/`, conformer metadata under `procdata/conformers/`, and WHIM assay tables under `procdata/assays/`.

Expected usage is usually smaller than `rawdata` but can still reach several GB for the full `drugs` subset. Files here are ignored by git.
