from damply import dirs

configfile: "config/pipeline.yaml"

rule all:
    input: 
        dirs.RESULTS/ "colData.csv",
        dirs.RESULTS/ "WHIM_hp.csv",
        dirs.RESULTS/ "WHIMS_avg.csv",
        dirs.RESULTS/ "WHIMS_wavg.csv"



rule process_geom:
    input:
        dirs.RAWDATA/ "rdkit_folder"/ "summary_drugs.json",
        dirs.RAWDATA/ "rdkit_folder"/ "summary_qm9.json"  

    output: 
        dirs.RESULTS/ "colData.csv",
        dirs.RESULTS/ "WHIM_hp.csv",
        dirs.RESULTS/ "WHIMS_avg.csv",
        dirs.RESULTS/ "WHIMS_wavg.csv"
    script:
        process_geom.py
        

