#!/bin/bash
set -e
set -o pipefail

# USAGE: ./runPHACE.sh "id"

id=$1

# RUST PIPELINE

cargo run tolerance "${id}"
cargo run msa1 "${id}"
cargo run msa2 "${id}"

# ASR

iqtree2 -s MSA1/"${id}_MSA1.fasta" -te "${id}.treefile" -blfix -m Data/vals_MSA1.txt -asr --prefix "${id}_MSA1" --seed 20 -T AUTO
iqtree2 -s MSA2/"${id}_MSA2.fasta" -te "${id}.treefile" -blfix -m Data/vals_MSA2.txt -asr --prefix "${id}_MSA2" --seed 20 -T AUTO

mv "${id}_MSA1"* MSA1/
mv "${id}_MSA2"* MSA2/

# RUST PIPELINE CONTINUED

cargo run msa1matrix "${id}"
cargo run msa2matrix "${id}"
cargo run matrices "${id}"
cargo run phace "${id}"

# DIRECTORIES

mkdir -p Inputs/
mv "${id}_MaskedMSA.fasta" "${id}.treefile" "${id}.state" Inputs/

dir="Test_Cases/${id}"
mkdir -p "${dir}"
mv MSA1 MSA2 Matrices ToleranceScores Inputs "${id}_PHACE.csv" "${dir}/"
