use bio::io::fasta;
use csv::ReaderBuilder;
use csv::Writer;
use ndarray::prelude::*;
use regex::Regex;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs;
use std::path::Path;


// ---------- Main ----------

pub fn phace(id: &str) -> Result<(), Box<dyn std::error::Error>> {
    // Define input/output paths (Matrices dir for inputs, MSA1/Root for others)
    let f_total = format!("Matrices/{}_mat_total_change.csv", id);
    let f_aa = format!("Matrices/{}_mat_aa_change.csv", id);
    let f_aff_num = format!("Matrices/{}_mat_aff_branch_num.csv", id);
    let f_aff_br = format!("Matrices/{}_mat_aff_branches.csv", id);
    let f_info = format!("Matrices/{}_mat_info.csv", id);
    let fasta_msa = format!("MSA1/{}_MSA1.fasta", id);
    let fasta_masked = format!("{}_MaskedMSA.fasta", id);
    let treefile = format!("MSA1/{}_MSA1.treefile", id);

    // Load matrices and branch info
    let mat_total_change = read_numeric_matrix(&f_total)?;
    let mat_aa_change = read_string_matrix(&f_aa)?;
    let mat_aff_branch_num = read_numeric_matrix(&f_aff_num)?;
    let mat_aff_branches = read_string_matrix(&f_aff_br)?;
    let (branch_lens, nodes) = read_info(&f_info)?;

    // Parse tree to get the correct order of tips
    let tip_order = extract_tips_from_newick(&treefile)?;
    if tip_order.is_empty() {
        return Err(format!("No tip labels parsed from treefile: {}", &treefile).into());
    }

    // Load FASTA sequences (standard and masked)
    let msa_map = read_fasta(&fasta_msa)?;
    let msa_org_map = read_fasta(&fasta_masked)?;

    // Align MSA data to match the tree's tip order
    let first_seq = msa_map.get(&tip_order[0]).ok_or("First tip not found in MSA fasta")?;
    let total_pos = first_seq.len();
    let mut msa: Vec<Vec<char>> = Vec::with_capacity(tip_order.len());
    let mut msa_org: Vec<Vec<char>> = Vec::with_capacity(tip_order.len());
    for name in &tip_order {
        let seq = msa_map.get(name).ok_or(format!("{} not found in MSA fasta", name))?;
        let seq_mask = msa_org_map.get(name).ok_or(format!("{} not found in Masked MSA fasta", name))?;
        msa.push(seq.clone());
        msa_org.push(seq_mask.clone());
    }

    let num_branch = branch_lens.len();
    let n_branches = mat_total_change.shape()[0];

    println!("Parsed tips: {}  positions: {}  branches(info): {}", tip_order.len(), total_pos, num_branch);
    println!("(mat_info nodes count = {})", nodes.len()); 

    // Calculate base weights for branches (inverse of mean change)
    // This normalizes fast-evolving vs slow-evolving branches
    let mut row_means: Vec<f64> = Vec::with_capacity(n_branches);
    for br in 0..n_branches {
        let row = mat_total_change.slice(s![br, ..]);
        let mean = row.mean().unwrap_or(0.0);
        row_means.push(mean);
    }
    let general_mean = if row_means.is_empty() { 1.0 } else { row_means.iter().sum::<f64>() / row_means.len() as f64 };
    for v in row_means.iter_mut() {
        if *v < general_mean { *v = general_mean; }
        *v = *v / general_mean;
    }
    let weight_branch_base: Vec<f64> = row_means.iter().map(|r| if *r == 0.0 { 1.0 } else { 1.0 / *r }).collect();

    let dash_re = Regex::new(r"-").unwrap();

    // Store results: (Pos1, Pos2, Score)
    let mut out: Vec<(usize, usize, f64)> = Vec::new();

    // Main Loop: Iterate over all unique pairs of positions
    for i1 in 0..(total_pos.saturating_sub(1)) {
        for i2 in (i1 + 1)..total_pos {
            let mut weight_branch = weight_branch_base.clone();

            // 1. Identify gapped sequences at either position
            let mut gap_names: HashSet<String> = HashSet::new();
            for (ridx, row) in msa.iter().enumerate() {
                if row[i1] == '-' { gap_names.insert(tip_order[ridx].clone()); }
                if row[i2] == '-' { gap_names.insert(tip_order[ridx].clone()); }
            }
            let _com_gap = gap_names.len();

            // 2. Check conservation in masked MSA
            let mut set1: HashSet<char> = HashSet::new();
            let mut set2: HashSet<char> = HashSet::new();
            for row in &msa_org {
                set1.insert(row[i1]);
                set2.insert(row[i2]);
            }
            let cons1 = set1.len();
            let cons2 = set2.len();

            // 3. Extract column vectors for current pair
            let dif1_col: Vec<f64> = mat_total_change.slice(s![.., i1]).to_vec();
            let dif2_col: Vec<f64> = mat_total_change.slice(s![.., i2]).to_vec();
            let num_eff1_col: Vec<f64> = mat_aff_branch_num.slice(s![.., i1]).to_vec();
            let num_eff2_col: Vec<f64> = mat_aff_branch_num.slice(s![.., i2]).to_vec();

            let pos1_aff_br_col: Vec<String> = mat_aff_branches.iter().map(|r| r[i1].clone()).collect();
            let pos2_aff_br_col: Vec<String> = mat_aff_branches.iter().map(|r| r[i2].clone()).collect();

            let dif_aa1_col: Vec<String> = mat_aa_change.iter().map(|r| r[i1].clone()).collect();
            let dif_aa2_col: Vec<String> = mat_aa_change.iter().map(|r| r[i2].clone()).collect();

            // 4. Handle Gap Events ("X-")
            // Find branches marked as "X-" (gap introduction)
            let w1_idx: Vec<usize> = dif_aa1_col.iter().enumerate().filter(|(_, v)| *v == "X-").map(|(i, _)| i).collect();
            let w2_idx: Vec<usize> = dif_aa2_col.iter().enumerate().filter(|(_, v)| *v == "X-").map(|(i, _)| i).collect();

            // Collect affected branches from the "X-" events
            let mut br1: Vec<String> = Vec::new();
            for &idx in &w1_idx {
                for piece in dash_re.split(&pos1_aff_br_col[idx]) {
                    if piece != "0" && !piece.is_empty() { br1.push(piece.to_string()); }
                }
            }
            let mut br2: Vec<String> = Vec::new();
            for &idx in &w2_idx {
                for piece in dash_re.split(&pos2_aff_br_col[idx]) {
                    if piece != "0" && !piece.is_empty() { br2.push(piece.to_string()); }
                }
            }
            
            // Adjust weights based on shared branches involved in gaps
            let mut uniq_br: HashSet<String> = HashSet::new();
            for b in br1.iter().chain(br2.iter()) { uniq_br.insert(b.clone()); }
            let ww = uniq_br.len();
            let com_ww = br1.iter().filter(|b| br2.contains(b)).count();
            let mut weight_inc1 = 1.0 - ((ww as f64 - com_ww as f64) / num_branch as f64);
            if weight_inc1 < 0.0 { weight_inc1 = 0.0; }

            let mut pos1_totalscore = dif1_col.clone();
            let mut pos2_totalscore = dif2_col.clone();

            // 5. Compute Coevolution Score
            let score3: f64;
            // Trivial conservation check (if one site is fully conserved, score is -1)
            if std::cmp::min(cons1, cons2) == 1 {
                score3 = -1.0;
            } else {
                // Determine locations of large differences
                let dif_tot: Vec<f64> = pos1_totalscore.iter().zip(pos2_totalscore.iter()).map(|(a,b)| a - b).collect();
                let loc1: Vec<usize> = dif_tot.iter().enumerate().filter(|(_, v)| **v >= 0.5).map(|(i, _)| i).collect();
                let loc2: Vec<usize> = dif_tot.iter().enumerate().filter(|(_, v)| **v <= -0.5).map(|(i, _)| i).collect();

                // Logic to zero out scores if change is driven by a single leaf observation (noise reduction)
                // Checks Pos1 (loc1)
                if !loc1.is_empty() {
                    let indices_to_check = if loc1.len() == 1 { vec![loc1[0]] } else { loc1.clone() };
                    for lx in indices_to_check {
                        let aff_lf = num_eff1_col[lx] as usize;
                        let aff_br_list: Vec<&str> = pos1_aff_br_col[lx].split('-').filter(|s| *s != "0" && !s.is_empty()).collect();
                        
                        if aff_lf == 1 && !aff_br_list.is_empty() {
                            // Find which tip corresponds to this branch
                            if let Some(row_idx) = tip_order.iter().position(|n| n == aff_br_list[0]) {
                                let alt = msa_org[row_idx][i1];
                                let obs_org = msa_org.iter().filter(|r| r[i1] == alt).count();
                                // If this state is observed only once, zero the score
                                if obs_org == 1 {
                                    pos1_totalscore[lx] = 0.0;
                                    pos2_totalscore[lx] = 0.0;
                                }
                            }
                        }
                    }
                }

                // Checks Pos2 (loc2) - Same logic as above
                if !loc2.is_empty() {
                    let indices_to_check = if loc2.len() == 1 { vec![loc2[0]] } else { loc2.clone() };
                    for lx in indices_to_check {
                        let aff_lf = num_eff2_col[lx] as usize;
                        let aff_br_list: Vec<&str> = pos2_aff_br_col[lx].split('-').filter(|s| *s != "0" && !s.is_empty()).collect();
                        
                        if aff_lf == 1 && !aff_br_list.is_empty() {
                             if let Some(row_idx) = tip_order.iter().position(|n| n == aff_br_list[0]) {
                                let alt = msa_org[row_idx][i2];
                                let obs_org = msa_org.iter().filter(|r| r[i2] == alt).count();
                                if obs_org == 1 {
                                    pos1_totalscore[lx] = 0.0;
                                    pos2_totalscore[lx] = 0.0;
                                }
                            }
                        }
                    }
                }

                // Zero out scores where only one side has a Gap event ("X-")
                let inc_gap1: Vec<usize> = (0..n_branches).filter(|&r| dif_aa1_col[r] == "X-" && dif_aa2_col[r] != "--" && dif_aa2_col[r] != "X-").collect();
                let inc_gap2: Vec<usize> = (0..n_branches).filter(|&r| dif_aa2_col[r] == "X-" && dif_aa1_col[r] != "--" && dif_aa1_col[r] != "X-").collect();
                for &idx in &inc_gap1 { pos1_totalscore[idx] = 0.0; }
                for &idx in &inc_gap2 { pos2_totalscore[idx] = 0.0; }

                // Final Scoring Validation
                let sum_p1: f64 = pos1_totalscore.iter().sum();
                let sum_p2: f64 = pos2_totalscore.iter().sum();
                let mx_p1: f64 = pos1_totalscore.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let mx_p2: f64 = pos2_totalscore.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

                // Ensure minimum signal exists before running WCCC
                if sum_p1 >= 1.0 && sum_p2 >= 1.0 && mx_p1 >= 0.5 && mx_p2 >= 0.5 {
                    // Update weights based on difference correlation
                    let df: Vec<f64> = pos1_totalscore.iter().zip(pos2_totalscore.iter()).map(|(a,b)| (a - b).abs()).collect();
                    let plc: Vec<usize> = df.iter().enumerate().filter(|(_, v)| **v >= 0.5).map(|(i, _)| i).collect();
                    for &p in &plc { weight_branch[p] = 1.0; }

                    // Construct combined correlation vector
                    let cc_vec: Vec<f64> = pos1_totalscore.iter().zip(pos2_totalscore.iter()).map(|(a,b)| a.max(*b)).collect();
                    let mut cc_nonzero = cc_vec.clone();
                    for c in cc_nonzero.iter_mut() { if *c == 0.0 { *c = 1.0; } }

                    // Final weight vector
                    let weight_branch2: Vec<f64> = weight_branch.iter().zip(cc_nonzero.iter()).map(|(w,c)| (w * c).sqrt()).collect();
                    
                    // Compute Weighted Concordance Correlation Coefficient (WCCC)
                    score3 = wccc(&pos1_totalscore, &pos2_totalscore, &weight_branch2) * weight_inc1;
                } else {
                    score3 = -1.0;
                }
            }

            println!("({},{})", i1 + 1, i2 + 1);
            out.push((i1 + 1, i2 + 1, score3));
        }
    }

    // Write final results to CSV
    let outname = format!("{}_PHACE.csv", id);
    let mut wtr = Writer::from_path(&outname)?;
    wtr.write_record(&["Pos1", "Pos2", "PHACE_Score"])?;
    for (p1, p2, score) in out {
        wtr.write_record(&[p1.to_string(), p2.to_string(), format!("{}", score)])?;
    }
    wtr.flush()?;

    Ok(())
}

// ---------- Helpers: reading matrices and FASTA ----------

fn read_numeric_matrix<P: AsRef<Path>>(path: P) -> Result<Array2<f64>, Box<dyn Error>> {
    let mut rdr = ReaderBuilder::new().has_headers(true).from_path(path)?;
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for result in rdr.records() {
        let rec = result?;
        let mut row: Vec<f64> = Vec::with_capacity(rec.len());
        for v in rec.iter() {
            let vtrim = v.trim().trim_matches('"');
            let val = if vtrim.is_empty() { 0.0 } else { vtrim.parse::<f64>().unwrap_or(0.0) };
            row.push(val);
        }
        rows.push(row);
    }
    if rows.is_empty() {
        return Ok(Array2::<f64>::zeros((0, 0)));
    }
    // Convert Vec<Vec<f64>> to Array2<f64>
    let nrow = rows.len();
    let ncol = rows[0].len();
    let mut arr = Array2::<f64>::zeros((nrow, ncol));
    for (i, r) in rows.into_iter().enumerate() {
        for (j, val) in r.into_iter().enumerate() {
            arr[[i, j]] = val;
        }
    }
    Ok(arr)
}

fn read_string_matrix<P: AsRef<Path>>(path: P) -> Result<Vec<Vec<String>>, Box<dyn Error>> {
    let mut rdr = ReaderBuilder::new().has_headers(true).from_path(path)?;
    let mut out: Vec<Vec<String>> = Vec::new();
    for result in rdr.records() {
        let rec = result?;
        out.push(rec.iter().map(|s| s.trim().trim_matches('"').to_string()).collect());
    }
    Ok(out)
}

#[derive(Debug, Deserialize)]
struct InfoRecord {
    #[serde(rename = "V1")]
    v1: String,
    #[serde(rename = "V2")]
    v2: String,
}
fn read_info<P: AsRef<Path>>(path: P) -> Result<(Vec<f64>, Vec<String>), Box<dyn Error>> {
    let mut rdr = ReaderBuilder::new().has_headers(true).from_path(path)?;
    let mut branch_lens: Vec<f64> = Vec::new();
    let mut nodes: Vec<String> = Vec::new();
    for result in rdr.deserialize() {
        let rec: InfoRecord = result?;
        let num = rec.v1.trim().trim_matches('"').parse::<f64>().unwrap_or(0.0);
        branch_lens.push(num);
        nodes.push(rec.v2.trim().trim_matches('"').to_string());
    }
    Ok((branch_lens, nodes))
}

fn read_fasta<P: AsRef<Path> + std::fmt::Debug>(path: P) -> Result<HashMap<String, Vec<char>>, Box<dyn Error>> {
    let reader = fasta::Reader::from_file(path)?;
    let mut map: HashMap<String, Vec<char>> = HashMap::new();
    for record in reader.records() {
        let rec = record?;
        let seq = rec.seq().iter().map(|&b| b as char).collect::<Vec<char>>();
        map.insert(rec.id().to_string(), seq);
    }
    Ok(map)
}

// Regex parsing for Newick format to extract tip labels in order
fn extract_tips_from_newick<P: AsRef<Path>>(path: P) -> Result<Vec<String>, Box<dyn Error>> {
    let s = fs::read_to_string(path)?;
    let re_single = Regex::new(r"[(,]\s*'([^']+)'\s*:").unwrap();
    let re_double = Regex::new(r#"[(,]\s*"([^"]+)"\s*:"#).unwrap();
    let re_unq = Regex::new(r"[(,]\s*([A-Za-z0-9_\.\-]+)\s*:").unwrap();

    let mut found: Vec<String> = Vec::new();
    // Capture different quoting styles for node names
    for cap in re_single.captures_iter(&s) {
        let lab = cap.get(1).unwrap().as_str().to_string();
        if !found.contains(&lab) { found.push(lab); }
    }
    for cap in re_double.captures_iter(&s) {
        let lab = cap.get(1).unwrap().as_str().to_string();
        if !found.contains(&lab) { found.push(lab); }
    }
    for cap in re_unq.captures_iter(&s) {
        let lab = cap.get(1).unwrap().as_str().to_string();
        if !found.contains(&lab) { found.push(lab); }
    }
    Ok(found)
}

// ---------- WCCC (Weighted Concordance Correlation Coefficient) ----------

fn weighted_mean(values: &[f64], weights: &[f64]) -> f64 {
    let sw: f64 = weights.iter().sum();
    if sw == 0.0 { 0.0 } else { values.iter().zip(weights.iter()).map(|(v,w)| v*w).sum::<f64>() / sw }
}

fn wccc(y: &[f64], x: &[f64], w: &[f64]) -> f64 {
    let wy = weighted_mean(y, w);
    let wx = weighted_mean(x, w);

    // Weighted Covariance
    let cov_w: f64 = y.iter().zip(x.iter()).zip(w.iter())
        .map(|((yi, xi), wi)| wi * (yi - wy) * (xi - wx)).sum();

    // Weighted Standard Deviations
    let sd_y_w = (y.iter().zip(w.iter()).map(|(yi, wi)| wi * (yi - wy).powi(2)).sum::<f64>()).sqrt();
    let sd_x_w = (x.iter().zip(w.iter()).map(|(xi, wi)| wi * (xi - wx).powi(2)).sum::<f64>()).sqrt();

    let denom = sd_y_w.powi(2) + sd_x_w.powi(2) + (wy - wx).powi(2);
    if denom == 0.0 { 0.0 } else { 2.0 * cov_w / denom }
}