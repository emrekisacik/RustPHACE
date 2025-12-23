use anyhow::{Context, Result};
use bio::io::fasta;
use std::fs;
use std::path::Path;

pub fn msa2(id: &str) -> Result<()> {

    // Define file paths dynamically based on the input ID
    let input_path = format!("{}_MaskedMSA.fasta", id);
    let output_path = format!("MSA2/{}_MSA2.fasta", id);

    // --- Step 1: Initialize Reader and Writer ---
    // Reads the masked multiple sequence alignment (MSA) file
    let reader = fasta::Reader::from_file(&input_path)
        .context(format!("Failed to read FASTA: {}", input_path))?;

    // Create the output directory "MSA2" if it doesn't already exist
    if let Some(parent) = Path::new(&output_path).parent() {
        fs::create_dir_all(parent)?;
    }

    // specific output file creation
    let file = fs::File::create(&output_path)
        .context(format!("Failed to create output file: {}", output_path))?;
    let mut writer = fasta::Writer::new(file);

    // --- Step 2: Process Sequences Streamwise ---
    // Iterates through the FASTA file record by record (memory efficient)
    for result in reader.records() {
        let record = result?;
        let original_seq = record.seq();

        // Transform the sequence based on gap presence
        // Rule:
        // - Gap ('-') -> 'G'
        // - Any other character (Amino Acid) -> 'C'
        let new_seq: Vec<u8> = original_seq.iter()
            .map(|&byte| {
                if byte == b'-' {
                    b'G' 
                } else {
                    b'C' 
                }
            })
            .collect();

        // Write the transformed sequence immediately to the output file
        writer.write(record.id(), record.desc(), &new_seq)?;
    }

    Ok(())
}