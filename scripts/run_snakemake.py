import os
import pandas as pd
configfile: "/scratch/user/pangkx/repos/bardlab_colabfold/scripts/250805_N_G3BP1/config.yaml"

DATA_DIR = config["data_dir"]
SAMPLE_PAIRS = ["N_NP_938405.1"]
INPUT = config["input_dir"]
OUTPUT = config["output_dir"]
LOG_DIR = config["log_dir"]
DATABASE = config["database_dir"]
IMAGE = config["image_path"]
CUSTOM_SEARCH = config["custom_search"]
MERGE_A3M = config["merge_a3m_script"]
VMTOUCH_ENV = config["vmtouch_env"]
DATABASE1 = config["db1"]
DATABASE3 = config["db3"]
CUSTOM_DB = config["custom_db"]
DATE = config["date"]
TMP = config["TMPDIR"]
WEIGHTS = config["weights_path"]
SPOC_CONTAINER = config["spoc_container"]
SPOC_REPO = config["spoc_repo_dir"]

rule all:
    input:
        expand(OUTPUT+"/spoc_score_{sample_pairs}", sample_pairs=SAMPLE_PAIRS)

rule colabfold_batch:
    input:
        merged_a3m = INPUT + "/merge_a3m/{sample_pairs}.a3m"
    params:
        prefix = "{sample_pairs}",
        weights = WEIGHTS,
        tmp = TMP
    output:
        directory(OUTPUT + "/colabfold_batch/{sample_pairs}")
    singularity:
        IMAGE
    shell:
        """
        colabfold_batch \
            --data {params.weights} \
            --num-recycle 3 \
            --recycle-early-stop-tolerance 0.5 \
            --num-ensemble 1 \
            --num-seeds 5 \
            --random-seed 0 \
            --num-models 3 \
            --model-order "3,2,1" \
            --model-type alphafold2_multimer_v1 \
            --rank multimer \
            --stop-at-score 85.0 \
            --templates \
            --jobname-prefix {params.prefix} \
            "{input}" \
            "{output}"
        """

rule spoc_score:
    input:
        directory(OUTPUT+"/colabfold_batch/{sample_pairs}")
    params:
        spoc_container = SPOC_CONTAINER,
        spoc_repo = SPOC_REPO,
        tmp = TMP 
    output:
        directory(OUTPUT+"/spoc_score_{sample_pairs}")
    container:
        SPOC_CONTAINER
    shell:
        """
        mkdir {output}
        source activate spoc_venv && \
        cd {input} && \
        python {params.spoc_repo}/scripts/run_custom_nobio.py {input} --rf_params {params.spoc_repo}/models/rf_afm_no_bio.joblib --output {output}/spoc_nobio_output.csv --ipsae_script {params.spoc_repo}/scripts/ipsae.py
        """

rule spoc_score_individual:
    input:
        directory(OUTPUT+"/colabfold_batch/{sample_pairs}")
    params:
        spoc_container = SPOC_CONTAINER,
        spoc_repo = SPOC_REPO,
        tmp = TMP 
    output:
        directory(OUTPUT+"/spoc_score_{sample_pairs}")
    container:
        SPOC_CONTAINER
    shell:
        """
        mkdir {output}
        source activate spoc_venv && \
        cd {input} && \
        python {params.spoc_repo}/scripts/run_custom_nobio.py {input} --rf_params {params.spoc_repo}/models/rf_afm_no_bio.joblib --output {output}/spoc_nobio_output.csv --ipsae_script {params.spoc_repo}/scripts/ipsae.py
        """