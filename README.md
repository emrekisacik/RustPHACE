# Rust Implementation of PHACE

The original R scripts have been converted to Rust for improved performance and modularity. Each step of the pipeline is triggered from a central entry point. Intermediate files are organized in directories during execution. The switch between Rust and IQTREE for ASR is automated using a Bash script. 

### Repository Structure

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
├── Data/                  # Configuration files (used by IQTREE)
│   ├── vals_MSA1.txt      # MSA1 substitution model
│   └── vals_MSA2.txt      # MSA2 substitution model
├── Test_Cases/            # Directory to organize the previous runs
├── Cargo.lock             # Locking dependency versions
├── Cargo.toml             # Project manifest and dependency configuration
└── runPHACE.sh            # Bash script to run the full pipeline (Rust + IQ-TREE)
```
### Pipeline Automation

Since ASR requires IQTREE, `runPHACE.sh` handles the switch between environments:
1. **Rust (Part 1):** Runs `tolerance.rs`, `msa1.rs`, and `msa2.rs`
2. **IQ-TREE (Terminal):** Executes IQ-TREE to perform ASR for both alignments.
3. **Rust (Part 2):** Runs `msa1matrix.rs`, `msa2matrix.rs`, `matrices.rs`, and `phace.rs`.
