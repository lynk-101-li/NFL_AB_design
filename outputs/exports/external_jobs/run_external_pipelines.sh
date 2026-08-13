#!/usr/bin/env bash
set -euo pipefail

# Generated command sheet for external structure and docking tools.
# Review config/external_pipelines.example.json before enabling jobs.

mkdir -p "real_runs/handoffs/nfl_design_20260813T1821587876520800"
# python3 ./scripts/prepare_real_model_jobs.py --request-dir ./outputs/exports/design_requests --target-manifest ./config/target_structure_manifest.json --profile smoke --job-scope canary --iggm-repo-dir ./third_party/IgGM --germinal-repo-dir ./third_party/germinal --germinal-af-params-dir <PATH_TO_AF_MULTIMER_PARAMS>

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/tfold_candidate_structure_prediction/selected_candidate_fv_chains"
# conda run -n tfold python ./third_party/tfold/projects/tfold_ab/predict.py --fasta outputs/exports/fasta/selected_candidate_fv_chains.fasta --output real_runs/results/nfl_design_20260813T1821587876520800/tfold_candidate_structure_prediction/selected_candidate_fv_chains

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/igfold_fv_modeling/selected_candidate_fv_chains"
# python3 -c <USE_THE_REVIEWED_IGFOLD_PYTHON_SNIPPET_IN_docs/real_model_installation.md>

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/immunebuilder_abodybuilder2_fv_modeling/selected_candidate_fv_chains"
# python3 -c <USE_THE_REVIEWED_ABODYBUILDER2_PYTHON_SNIPPET_IN_docs/real_model_installation.md>

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-e6e8d7351d0d_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-e6e8d7351d0d_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-e6e8d7351d0d_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-8df6f8cc4294_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-8df6f8cc4294_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-8df6f8cc4294_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-77576c381425_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-77576c381425_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-77576c381425_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-a9cf24db579e_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-a9cf24db579e_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-a9cf24db579e_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-c0d75e927e22_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-c0d75e927e22_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-c0d75e927e22_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-2ea27b071345_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-2ea27b071345_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-2ea27b071345_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-db51a6170580_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-db51a6170580_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-db51a6170580_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-870775f54c72_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-870775f54c72_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-870775f54c72_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-98a428e115e9_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-98a428e115e9_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-98a428e115e9_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-6d8910dd2d15_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-6d8910dd2d15_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-6d8910dd2d15_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-a78bf47c4ea2_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-a78bf47c4ea2_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-a78bf47c4ea2_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-9c287784eedf_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_complex_DN-9c287784eedf_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_complex_DN-9c287784eedf_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-377"
# python ./third_party/alphafold3/run_alphafold.py --json_path outputs/exports/af3_json/af3_sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-377.json --model_dir <PATH_TO_AF3_MODEL_DIR> --db_dir <PATH_TO_AF3_DATABASE_DIR> --output_dir real_runs/results/nfl_design_20260813T1821587876520800/alphafold3_complex/af3_sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-e6e8d7351d0d_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-e6e8d7351d0d_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-e6e8d7351d0d_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-8df6f8cc4294_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-8df6f8cc4294_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-8df6f8cc4294_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-77576c381425_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-77576c381425_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-77576c381425_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-a9cf24db579e_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-a9cf24db579e_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-a9cf24db579e_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-c0d75e927e22_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-c0d75e927e22_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-c0d75e927e22_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-2ea27b071345_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-2ea27b071345_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-2ea27b071345_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-db51a6170580_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-db51a6170580_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-db51a6170580_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-870775f54c72_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-870775f54c72_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-870775f54c72_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-98a428e115e9_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-98a428e115e9_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-98a428e115e9_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-6d8910dd2d15_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-6d8910dd2d15_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-6d8910dd2d15_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-a78bf47c4ea2_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-a78bf47c4ea2_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-a78bf47c4ea2_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-9c287784eedf_NEFL_280-377"
# conda run -n chai1 chai-lab fold outputs/exports/fasta/complex_DN-9c287784eedf_NEFL_280-377.fasta real_runs/results/nfl_design_20260813T1821587876520800/chai1_complex/complex_DN-9c287784eedf_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-e6e8d7351d0d_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-e6e8d7351d0d_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-e6e8d7351d0d_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-8df6f8cc4294_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-8df6f8cc4294_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-8df6f8cc4294_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-77576c381425_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-77576c381425_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-77576c381425_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-a9cf24db579e_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-a9cf24db579e_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-a9cf24db579e_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-c0d75e927e22_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-c0d75e927e22_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-c0d75e927e22_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-2ea27b071345_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-2ea27b071345_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-2ea27b071345_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-db51a6170580_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-db51a6170580_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-db51a6170580_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-870775f54c72_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-870775f54c72_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-870775f54c72_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-98a428e115e9_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-98a428e115e9_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-98a428e115e9_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-6d8910dd2d15_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-6d8910dd2d15_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-6d8910dd2d15_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-a78bf47c4ea2_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-a78bf47c4ea2_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-a78bf47c4ea2_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-9c287784eedf_NEFL_280-377"
# conda run -n boltz2 boltz predict <BOLTZ_YAML_DERIVED_FROM_outputs/exports/fasta/complex_DN-9c287784eedf_NEFL_280-377.fasta> --use_msa_server --out_dir real_runs/results/nfl_design_20260813T1821587876520800/boltz_complex/complex_DN-9c287784eedf_NEFL_280-377

mkdir -p "real_runs/results/nfl_design_20260813T1821587876520800/rosetta_interface_analysis/sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-377"
# <PATH_TO_ROSETTA_INTERFACE_ANALYZER> --input outputs/exports/fasta/sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-377.fasta --database <PATH_TO_ROSETTA_DATABASE> --outdir real_runs/results/nfl_design_20260813T1821587876520800/rosetta_interface_analysis/sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-377
