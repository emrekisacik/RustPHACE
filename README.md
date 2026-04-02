# rustyPHACE

rustyPHACE is a high-performance, Rust-based implementation of the Phylogeny-Aware detection of molecular CoEvolution (PHACE) pipeline. It is a command-line tool designed to calculate co-evolution scores between amino acid positions using sequence alignments and phylogenetic trees.

**Note:** This repository contains the optimized Rust implementation of the original R-based version of PHACE, which can be found here: [https://github.com/CompGenomeLab/PHACE](https://github.com/CompGenomeLab/PHACE).

## Workflow

The entire workflow is automated via `runPHACE.sh`, running Rust modules and IQ-TREE. The steps are:

1. Tolerance Scoring (`tolerance`): Calculates position-specific tolerance scores.
2. MSA1 Generation (`msa1`): Generates the first modified MSA.
3. MSA2 Generation (`msa2`): Generates the second modified MSA based on gap profiles.
4. Ancestral State Reconstruction (ASR): Runs IQ-TREE to infer ancestral sequences for both generated MSAs.
5. Matrix Generation (`msa1matrix`, `msa2matrix`, `matrices`): Computes evolutionary changes, effective branch lengths, and transition matrices.
6. PHACE Scoring (`phace`): Calculates final pairwise co-evolution scores.
7. File Organization: Automatically moves output files into `Test_Cases/[id]/`.

## Prerequisites

To run RustPHACE, the following dependencies must be installed and available in your system's PATH:

* **Rust & Cargo** (v1.70 or newer)
* **IQ-TREE 2**
* **Bash environment** (Linux, macOS, or WSL for Windows)

## Installation

Clone the repository to your local machine. You must grant execution permissions to the main bash script using `chmod +x`.

```bash
git clone [https://github.com/CompGenomeLab/rustyPHACE.git](https://github.com/CompGenomeLab/rustyPHACE.git)
cd rustyPHACE
chmod +x runPHACE.sh
```

## Repository Structure

```text
PHACE/
├── src/                   # Source code directory
│   ├── lib.rs             # Library root
│   ├── main.rs            # Entry point
│   ├── tolerance.rs       # Calculate tolerance scores (ToleranceScore.R)
│   ├── msa1.rs            # Generate MSA1 (MSA1.R)
│   ├── msa2.rs            # Generate MSA2 (MSA2.R)
│   ├── msa1matrix.rs      # Process MSA1 results (Part1_MSA1.R)
│   ├── msa2matrix.rs      # Process MSA2 results (Part1_MSA2.R)
│   ├── matrices.rs        # Merge matrices (GetTotalChangeMatrix.R)
│   ├── phace.rs           # PHACE algorithm (PHACE.R)
│   ├── interphace.rs      # Optimized PHACE algorithm for interprotein analysis (NEW!)
├── Data/                  # Configuration files (used by IQTREE during ASR)
│   ├── vals_MSA1.txt      # MSA1 substitution model
│   └── vals_MSA2.txt      # MSA2 substitution model
├── Test_Cases/            # Directory to organize previous runs
├── Cargo.lock             # Locking dependency versions
├── Cargo.toml             # Project manifest and dependency configuration
└── runPHACE.sh            # Bash script to run the full pipeline (Rust + IQ-TREE)
```

## Usage

The pipeline requires a specific naming convention for your input files. Before running the pipeline, ensure the following files exist in your root directory, replacing `[id]` with your unique run identifier:

* `[id]_MaskedMSA.fasta` (The masked multiple sequence alignment)
* `[id].treefile` (The phylogenetic tree)
* `[id].state` (IQ-TREE state file)

### 1. Automated Execution (Recommended)

Execute the main script by simply providing your `[id]`:

```bash
./runPHACE.sh <id>
```

### 2. Manual Execution (Step-by-Step)

If you prefer not to use the automated bash script, or if the process is interrupted and you need to resume from a specific step, you can execute the modules manually.

**Step 1: Generate initial MSAs and Tolerance Scores**
```bash
cargo run tolerance <id>
cargo run msa1 <id>
cargo run msa2 <id>
```

**Step 2: Ancestral State Reconstruction (IQ-TREE)**
```bash
iqtree2 -s MSA1/<id>_MSA1.fasta -te <id>.treefile -blfix -m Data/vals_MSA1.txt -asr --prefix <id>_MSA1 -T AUTO
iqtree2 -s MSA2/<id>_MSA2.fasta -te <id>.treefile -blfix -m Data/vals_MSA2.txt -asr --prefix <id>_MSA2 -T AUTO

# Move generated files to their respective directories
mv <id>_MSA1* MSA1/
mv <id>_MSA2* MSA2/
```

**Step 3: Matrix Generation and Final Scoring**
```bash
cargo run msa1matrix <id>
cargo run msa2matrix <id>
cargo run matrices <id>
cargo run phace <id>
```

## Outputs

Upon completion (if using the automated script), a new directory will be created at `Test_Cases/[id]/`. This folder will contain:

* `Inputs/`: The original input files used for the run.
* `MSA1/` & `MSA2/`: Reconstructed MSA files from IQ-TREE.
* `Matrices/`: Intermediate numerical and transition matrices.
* `ToleranceScores/`: Calculated tolerance score data.
* `[id]_PHACE.csv`: The final output file containing the pairwise co-evolution scores.

## License

This project is licensed under the MIT License.