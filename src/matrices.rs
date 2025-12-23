use std::error::Error;
use std::fs::{self, File};
use std::collections::HashMap;
use serde::Serialize;
use csv::QuoteStyle;

// --- Data Structures ---

#[derive(Serialize, Clone)]
pub struct BranchInfo {
    length: f64,
    name: String,
}

// --- Helper Functions ---

// Generic helper to read a CSV into a vector of string vectors (ignoring types for flexible processing)
fn read_csv_to_matrix(path: &str) -> Result<(Vec<String>, Vec<Vec<String>>), Box<dyn Error>> {
    let file = File::open(path).map_err(|e| format!("Failed to open input file '{}': {}", path, e))?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true) 
        .from_reader(file);

    let headers = rdr.headers()?.iter().map(|s| s.to_string()).collect();

    let mut matrix = Vec::new();
    for result in rdr.records() {
        let record = result?;
        let row: Vec<String> = record.iter().map(|s| s.to_string()).collect();
        matrix.push(row);
    }
    Ok((headers, matrix))
}

// Generic helper to write matrices to CSV with consistent quoting
fn write_matrix_to_csv<T: std::fmt::Display>(
    path: &str, 
    matrix: &[Vec<T>], 
    headers: &[String]
) -> Result<(), Box<dyn Error>> {
    
    let mut wtr = csv::WriterBuilder::new()
        .quote_style(QuoteStyle::Always) // Force quotes for compatibility
        .from_path(path)?;
    
    // Write headers
    wtr.write_record(headers)?;
    
    // Write data rows
    for row in matrix {
        let string_row: Vec<String> = row.iter().map(|x| x.to_string()).collect();
        wtr.write_record(&string_row)?;
    }
    wtr.flush()?;
    Ok(())
}

// --- MAIN FUNCTION ---

pub fn matrices(id: &str) -> Result<(), Box<dyn Error>> {
    
    fs::create_dir_all("Matrices")?;

    // 1. LOAD INPUT DATA
    // Read MSA1 (standard substitution) matrices
    let (_, msa1_aa_changes) = read_csv_to_matrix(&format!("MSA1/{}_MSA1_aa_changes.csv", id))?;
    let (_, msa1_all_changes) = read_csv_to_matrix(&format!("MSA1/{}_MSA1_all_changes.csv", id))?;
    let (_, msa1_aff_leaves) = read_csv_to_matrix(&format!("MSA1/{}_MSA1_aff_leaves.csv", id))?;
    let (_, msa1_aff_branches) = read_csv_to_matrix(&format!("MSA1/{}_MSA1_aff_branches.csv", id))?;

    // Read MSA2 (gap/indel) matrices
    let (_, msa2_all_changes) = read_csv_to_matrix(&format!("MSA2/{}_MSA2_all_changes.csv", id))?;
    let (_, msa2_aa_changes) = read_csv_to_matrix(&format!("MSA2/{}_MSA2_aa_changes.csv", id))?;
    let (_, msa2_aff_leaves) = read_csv_to_matrix(&format!("MSA2/{}_MSA2_aff_leaves.csv", id))?;
    let (_, msa2_aff_branches) = read_csv_to_matrix(&format!("MSA2/{}_MSA2_aff_branches.csv", id))?;

    // 2. PREPARE STRUCTURES
    let num_rows = msa1_aa_changes.len();
    let row_width = msa1_aa_changes[0].len();
    
    // Calculate total sequence positions (Columns - Label - BranchLength)
    let total_pos = row_width - 2;
    let output_headers: Vec<String> = (1..=total_pos).map(|i| format!("V{}", i)).collect();
    
    let mut branch_infos: Vec<BranchInfo> = Vec::with_capacity(num_rows);
    let mut branch_map: HashMap<String, usize> = HashMap::new();

    // Map branch names to row indices for quick lookup during merging
    for (idx, row) in msa1_aa_changes.iter().enumerate() {
        let name = row[0].clone();
        let length = row[row_width - 1].parse::<f64>().unwrap_or(0.0);
        
        branch_map.insert(name.clone(), idx);
        branch_infos.push(BranchInfo { length, name: name.clone() });
    }

    // Initialize result matrices with default values based on MSA1 size
    let mut res_all_changes = vec![vec![0.0; total_pos]; num_rows];
    let mut res_aa_changes = vec![vec![String::new(); total_pos]; num_rows];
    let mut res_aff_num = vec![vec![0.0; total_pos]; num_rows];
    let mut res_aff_br = vec![vec![String::new(); total_pos]; num_rows];

    // 3. MERGE LOGIC
    // Iterate through every position (column)
    for i in 0..total_pos {
        let data_col_idx = i + 1; // Skip 'Label' column

        // Step A: Populate defaults from MSA1
        for row_idx in 0..num_rows {
            let val_dif1_main = msa1_all_changes[row_idx][data_col_idx].parse::<f64>().unwrap_or(0.0);
            let mut val_aa1_main = msa1_aa_changes[row_idx][data_col_idx].clone();
            let val_num_eff1_main = msa1_aff_leaves[row_idx][data_col_idx].parse::<f64>().unwrap_or(0.0);
            let val_pos1_aff_br_main = msa1_aff_branches[row_idx][data_col_idx].clone();

            // Normalize specific non-informative states
            if val_aa1_main == "C-" || val_aa1_main == "A-" {
                val_aa1_main = "--".to_string();
            }

            res_all_changes[row_idx][i] = val_dif1_main;
            res_aa_changes[row_idx][i] = val_aa1_main;
            res_aff_num[row_idx][i] = val_num_eff1_main;
            res_aff_br[row_idx][i] = val_pos1_aff_br_main;
        }

        // Step B: Override with MSA2 data if a Gap event ("CG") is detected
        // Note: MSA2 tracks presence/absence (Gap vs Char), "CG" implies Change to Gap
        for r in 0..num_rows {
            let val_aa1_gap = &msa2_aa_changes[r][data_col_idx];
            
            if val_aa1_gap == "CG" {
                let branch_name_gap = &msa2_aa_changes[r][0];
                
                // Find corresponding row in MSA1 matrix using the branch map
                if let Some(&target_idx) = branch_map.get(branch_name_gap) {
                    let gap_dif1 = msa2_all_changes[r][data_col_idx].parse::<f64>().unwrap_or(0.0);
                    let gap_num_eff = msa2_aff_leaves[r][data_col_idx].parse::<f64>().unwrap_or(0.0);
                    let gap_aff_br = &msa2_aff_branches[r][data_col_idx];

                    // Overwrite MSA1 data with MSA2 gap data
                    res_all_changes[target_idx][i] = gap_dif1;
                    res_aff_num[target_idx][i] = gap_num_eff;
                    res_aff_br[target_idx][i] = gap_aff_br.clone();
                    // Mark as 'X-' to denote a gap event override
                    res_aa_changes[target_idx][i] = "X-".to_string();
                }
            }
        }
    }

    // 4. WRITE OUTPUTS
    // Save all merged matrices to the 'Matrices' directory
    let f1 = format!("Matrices/{}_mat_total_change.csv", id);
    write_matrix_to_csv(&f1, &res_all_changes, &output_headers)?;

    let f2 = format!("Matrices/{}_mat_aa_change.csv", id);
    write_matrix_to_csv(&f2, &res_aa_changes, &output_headers)?;

    let f3 = format!("Matrices/{}_mat_aff_branch_num.csv", id);
    write_matrix_to_csv(&f3, &res_aff_num, &output_headers)?;

    let f4 = format!("Matrices/{}_mat_aff_branches.csv", id);
    write_matrix_to_csv(&f4, &res_aff_br, &output_headers)?;

    // Write Branch Info (Lengths and Names)
    let f5 = format!("Matrices/{}_mat_info.csv", id);
    let mut wtr_info = csv::WriterBuilder::new()
        .quote_style(QuoteStyle::Always)
        .from_path(f5)?;
        
    wtr_info.write_record(&["V1", "V2"])?; 
    for info in branch_infos {
        wtr_info.write_record(&[info.length.to_string(), info.name])?;
    }
    wtr_info.flush()?;

    Ok(())
}