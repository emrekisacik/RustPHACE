use anyhow::{Context, Result};
use bio::io::fasta;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

const EPS: f64 = 1e-15;

pub fn msa1(id: &str) -> Result<()> {

    // Define file paths dynamically based on the input ID
    let csv_path = format!("ToleranceScores/{}.csv", id);
    let fasta_path = format!("{}_MaskedMSA.fasta", id);
    let output_path = format!("MSA1/{}_MSA1.fasta", id);

    // --- Step 1: Read and Transform Tolerance Scores ---
    let mut rdr = csv::ReaderBuilder::new()
        .from_path(&csv_path)
        .context(format!("Failed to read CSV: {}", csv_path))?;

    // Parse the header to map column indices to Amino Acid names
    let headers = rdr.headers()?.clone();
    let aa_names: Vec<String> = headers.iter().skip(1).take(20).map(|s| s.to_string()).collect();

    let mut scores_matrix: Vec<HashMap<String, f64>> = Vec::new();
    let log_eps = EPS.ln();

    // Iterate through each row (position) in the CSV
    for result in rdr.records() {
        let record = result?;
        let mut row_scores = HashMap::new();

        // Apply logarithmic transformation to normalize scores: 1 - log(val + eps) / log(eps)
        for (i, aa) in aa_names.iter().enumerate() {
            let val: f64 = record[i + 1].parse().unwrap_or(0.0);
            let transformed_val = 1.0 - (val + EPS).ln() / log_eps;
            row_scores.insert(aa.clone(), transformed_val);
        }
        scores_matrix.push(row_scores);
    }

    // --- Step 2: Read Input FASTA MSA ---
    // Note: The 'mut' keyword was removed here as per previous compiler warning fixes
    let reader = fasta::Reader::from_file(&fasta_path)
        .context(format!("Failed to read FASTA: {}", fasta_path))?;

    // Load all sequences into memory to allow in-place modification
    let mut records: Vec<(String, Vec<u8>)> = reader.records()
        .map(|r| {
            let rec = r.unwrap();
            (rec.id().to_string(), rec.seq().to_vec())
        })
        .collect();

    if records.is_empty() {
        return Ok(());
    }

    let seq_len = records[0].1.len();

    // --- Step 3: Process Each Column of the MSA ---
    for i in 0..seq_len {
        // Ensure we don't exceed the number of rows in our score matrix
        if i >= scores_matrix.len() {
            break; 
        }

        let sc = &scores_matrix[i];

        // Frequency count: Determine which Amino Acid is most common at this position (ignoring gaps)
        let mut counts: HashMap<u8, usize> = HashMap::new();
        for (_, seq) in &records {
            let aa = seq[i];
            if aa != b'-' {
                *counts.entry(aa).or_insert(0) += 1;
            }
        }

        // Skip processing if the column contains only gaps
        if counts.is_empty() {
            continue; 
        }

        // Identify the Dominant Amino Acid (highest frequency)
        let dom_byte = counts.into_iter()
            .max_by(|a, b| a.1.cmp(&b.1))
            .map(|(k, _)| k)
            .unwrap();
        
        let dom_str = String::from_utf8(vec![dom_byte]).unwrap();

        // Define Tolerance Threshold: Any AA with a score >= the dominant AA's score is "tolerated"
        let dom_score = sc.get(&dom_str).unwrap_or(&0.0);
        
        let mut tolerated_aas: Vec<u8> = Vec::new();
        tolerated_aas.push(dom_byte); // Ensure the dominant AA itself is included

        for (aa_name, score) in sc {
            if score > dom_score {
                if let Some(byte) = aa_name.bytes().next() {
                    tolerated_aas.push(byte);
                }
            }
        }

        // --- Step 4: Mutate Sequences Based on Tolerance ---
        // Rule:
        // - Gap ('-') -> Remains Gap
        // - Tolerated AA -> Mutate to 'C' (Conserved/Compatible)
        // - Non-Tolerated AA -> Mutate to 'A' (Aberrant/Different)
        for (_, seq) in &mut records {
            let aa = seq[i];
            if aa == b'-' {
                continue;
            } else if tolerated_aas.contains(&aa) {
                seq[i] = b'C';
            } else {
                seq[i] = b'A';
            }
        }
    }

    // --- Step 5: Write the Modified MSA to Output File ---
    if let Some(parent) = Path::new(&output_path).parent() {
        fs::create_dir_all(parent)?;
    }

    let file = fs::File::create(&output_path)
        .context(format!("Failed to create output file: {}", output_path))?;
    let mut writer = fasta::Writer::new(file);

    for (id, seq) in records {
        writer.write(&id, None, &seq)?;
    }

    Ok(())
}