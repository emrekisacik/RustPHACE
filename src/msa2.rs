use anyhow::{Context, Result};
use bio::io::fasta;
use indicatif::{ProgressBar, ProgressStyle};
use std::fs;
use std::path::Path;

pub fn msa2(id: &str) -> Result<()> {
    let input_path = format!("{}_MaskedMSA.fasta", id);
    let output_path = format!("MSA2/{}_MSA2.fasta", id);

    let reader = fasta::Reader::from_file(&input_path)
        .context(format!("Failed to read FASTA: {}", input_path))?;

    if let Some(parent) = Path::new(&output_path).parent() {
        fs::create_dir_all(parent)?;
    }

    let file = fs::File::create(&output_path)
        .context(format!("Failed to create output file: {}", output_path))?;
    let mut writer = fasta::Writer::new(file);

    let records: Vec<_> = reader.records().collect::<Result<Vec<_>, _>>()?;
    
    let pb = ProgressBar::new(records.len() as u64);
    println!("\x1b[1;32mGenerating MSA2...\x1b[0m");
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.green.bold}] {pos}/{len} ({eta})").unwrap()
        .progress_chars("━╸ "));

    for record in records {
        pb.inc(1);
        let original_seq = record.seq();

        let new_seq: Vec<u8> = original_seq.iter()
            .map(|&byte| {
                if byte == b'-' {
                    b'G' 
                } else {
                    b'C' 
                }
            })
            .collect();

        writer.write(record.id(), record.desc(), &new_seq)?;
    }
    
    pb.finish();

    Ok(())
}