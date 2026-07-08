#!/usr/bin/env bash
set -euo pipefail

# Generated command sheet for external structure and docking tools.
# Review config/external_pipelines.example.json before enabling jobs.

mkdir -p "outputs/external_results/igfold_fv_modeling/validated_fv_chains"
# igfold --fasta outputs/exports/fasta/validated_fv_chains.fasta --outdir outputs/external_results/igfold_fv_modeling/validated_fv_chains

mkdir -p "outputs/external_results/abodybuilder3_fv_modeling/validated_fv_chains"
# abodybuilder3 --fasta outputs/exports/fasta/validated_fv_chains.fasta --output outputs/external_results/abodybuilder3_fv_modeling/validated_fv_chains

mkdir -p "outputs/external_results/alphafold3_complex/af3_complex_7-H11-D3-2-C7_NEFL_280-375"
# run_alphafold3 --json_path outputs/exports/af3_json/af3_complex_7-H11-D3-2-C7_NEFL_280-375.json --output_dir outputs/external_results/alphafold3_complex/af3_complex_7-H11-D3-2-C7_NEFL_280-375

mkdir -p "outputs/external_results/alphafold3_complex/af3_complex_15-C12-H6_NEFL_280-375"
# run_alphafold3 --json_path outputs/exports/af3_json/af3_complex_15-C12-H6_NEFL_280-375.json --output_dir outputs/external_results/alphafold3_complex/af3_complex_15-C12-H6_NEFL_280-375

mkdir -p "outputs/external_results/alphafold3_complex/af3_sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-375"
# run_alphafold3 --json_path outputs/exports/af3_json/af3_sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-375.json --output_dir outputs/external_results/alphafold3_complex/af3_sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-375

mkdir -p "outputs/external_results/chai1_complex/complex_7-H11-D3-2-C7_NEFL_280-375"
# chai-lab fold outputs/exports/fasta/complex_7-H11-D3-2-C7_NEFL_280-375.fasta --output-dir outputs/external_results/chai1_complex/complex_7-H11-D3-2-C7_NEFL_280-375

mkdir -p "outputs/external_results/chai1_complex/complex_15-C12-H6_NEFL_280-375"
# chai-lab fold outputs/exports/fasta/complex_15-C12-H6_NEFL_280-375.fasta --output-dir outputs/external_results/chai1_complex/complex_15-C12-H6_NEFL_280-375

mkdir -p "outputs/external_results/boltz_complex/complex_7-H11-D3-2-C7_NEFL_280-375"
# boltz predict outputs/exports/fasta/complex_7-H11-D3-2-C7_NEFL_280-375.fasta --out_dir outputs/external_results/boltz_complex/complex_7-H11-D3-2-C7_NEFL_280-375

mkdir -p "outputs/external_results/boltz_complex/complex_15-C12-H6_NEFL_280-375"
# boltz predict outputs/exports/fasta/complex_15-C12-H6_NEFL_280-375.fasta --out_dir outputs/external_results/boltz_complex/complex_15-C12-H6_NEFL_280-375

mkdir -p "outputs/external_results/rosetta_interface_analysis/sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-375"
# rosetta_interface_analyzer --input outputs/exports/fasta/sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-375.fasta --outdir outputs/external_results/rosetta_interface_analysis/sandwich_7-H11-D3-2-C7_15-C12-H6_NEFL_280-375
