use anyhow::{Context, Result, anyhow};
use bio::io::fasta;
use std::collections::{HashMap, HashSet};
use std::fs;

// STRUCTS

#[derive(Clone, Debug)]
pub struct TreeInfo {
    pub parent: Vec<usize>,
    pub node: Vec<usize>,
    pub label: Vec<String>,
    pub branch_length: Vec<f64>,
}

#[derive(Debug)]
pub struct DiffOutput {
    pub diff1: Vec<f64>,
    pub mat2: Vec<String>,
    pub eff_br: Vec<f64>,
    pub aff_br: Vec<String>,
    pub extra: (f64, String, f64, String),
}

#[derive(Debug, Clone)]
struct SimpleNode {
    label: String,
    len: f64,
    children: Vec<SimpleNode>,
}

// Normalize Floats
fn normalize(val: f64) -> f64 {
    // If value is exactly 0.0 (or -0.0), return plain 0.0 to avoid "-0.0" in output
    if val == 0.0 { 0.0 } else { val }
}

// MAIN EXECUTION

pub fn msa2matrix(id: &str) -> Result<()> {

    // Define input/output paths
    let file_fasta = format!("MSA2/{}_MSA2.fasta", id);
    let file_nwk = format!("MSA2/{}_MSA2.treefile", id);
    let file_rst = format!("MSA2/{}_MSA2.state", id);
    
    let out_all_changes = format!("MSA2/{}_MSA2_all_changes.csv", id);
    let out_aa_changes = format!("MSA2/{}_MSA2_aa_changes.csv", id);
    let out_aff_leaves = format!("MSA2/{}_MSA2_aff_leaves.csv", id);
    let out_aff_branches = format!("MSA2/{}_MSA2_aff_branches.csv", id);

    // 1. LOAD TREE
    use std::io::{stdout, Write};
    let _ = stdout().flush();

    // Parse Newick tree into flat arrays (parents, nodes, labels)
    let (tree_info, num_leaves, num_nodes, node_map) = process_newick_custom(&file_nwk)
        .context("Failed to parse Newick file")?;

    let num_branch = tree_info.parent.len(); 

    // 2. LOAD ANCESTRAL STATES
    let _ = stdout().flush();
    
    // Initialize CSV reader for the .state file (tab-delimited)
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .comment(Some(b'#'))
        .from_path(&file_rst)?;

    let mut x_matrix: Vec<Vec<String>> = Vec::new();
    let mut max_site_idx = 0;
    // Map internal NodeID -> Site Index -> Character State
    let mut internal_seqs: HashMap<usize, HashMap<usize, char>> = HashMap::new();

    // Iterate over ancestral state reconstruction records
    for (idx, result) in rdr.records().enumerate() {
        let record = result?;
        if record.is_empty() { continue; }

        let node_raw = record.get(0).unwrap_or("").replace("Node", "");
        let site_str = record.get(1).unwrap_or("0");
        let state_str = record.get(2).unwrap_or("-");
        
        let site: usize = site_str.parse().unwrap_or_else(|_| {
            eprintln!("Warning: Row {} has invalid site '{}'", idx + 2, site_str);
            0
        });

        if site > max_site_idx { max_site_idx = site; }

        let state_char = state_str.chars().next().unwrap_or('-');

        // Map the raw node label from .state file to our internal tree index
        let full_label = format!("Node{}", node_raw);
        if let Some(&node_id) = node_map.get(&full_label) {
            internal_seqs.entry(node_id).or_default().insert(site, state_char);
        }

        // Store raw prob data: Node, Site, State, p_val1, p_val2
        let mut row = Vec::new();
        row.push(node_raw);
        row.push(site_str.to_string());
        row.push(state_str.to_string());
        row.push(record.get(4).unwrap_or("0.0").to_string()); // p_C or similar
        row.push(record.get(5).unwrap_or("0.0").to_string()); // p_G or similar
        x_matrix.push(row);
    }

    // Reconstruct full sequences for internal nodes based on parsed map
    let total_pos = max_site_idx;
    let max_id = tree_info.label.len();
    let mut fasta_node = vec![vec!['-'; total_pos + 1]; max_id];

    for (node_id, site_map) in internal_seqs {
        if node_id < fasta_node.len() {
            for (site, aa) in site_map {
                if site > 0 && site <= total_pos {
                    fasta_node[node_id][site - 1] = aa; 
                }
            }
        }
    }

    // 3. LOAD FASTA (Tips)
    let fa_reader = fasta::Reader::from_file(&file_fasta)?;
    let mut temp_msa: HashMap<String, Vec<char>> = HashMap::new();
    
    for result in fa_reader.records() {
        let record = result?;
        let seq: Vec<char> = record.seq().iter().map(|&b| b as char).collect();
        temp_msa.insert(record.id().to_string(), seq);
    }

    // Align tip sequences to the tree leaf order
    let mut msa_upd: Vec<Vec<char>> = Vec::new();
    for i in 0..num_leaves {
        let label = &tree_info.label[i];
        if let Some(seq) = temp_msa.get(label) {
            msa_upd.push(seq.clone());
        } else {
             return Err(anyhow!("Error: Tip label '{}' in tree not found in FASTA.", label));
        }
    }

    // 4. TOPOLOGY LEVELS
    // Calculate processing order (leaves -> root or vice versa) for propagation
    let levels = calculate_levels_exact(&tree_info, num_leaves);

    // 5. INITIALIZE OUTPUTS
    let num_rows = num_branch + 1; 
    let mut data_all_changes = vec![vec!["0".to_string(); total_pos]; num_rows];
    let mut data_aa_changes = vec![vec!["0".to_string(); total_pos]; num_rows];
    let mut data_aff_leaves = vec![vec!["0".to_string(); total_pos]; num_rows];
    let mut data_aff_branches = vec![vec!["0".to_string(); total_pos]; num_rows];

    // 6. MAIN COMPUTATION LOOP
    // Iterate through every position in the alignment
    let ord = &tree_info.label[tree_info.label.len() - 1]; 
    let nodes_raxml_prev = node_map.clone(); 

    for ps in 0..total_pos {
        let diff_res = compute_difference(
            ps,
            ord,
            &levels,
            &tree_info,
            num_nodes,
            num_leaves,
            &x_matrix,
            &nodes_raxml_prev,
            &msa_upd,
            &fasta_node
        );

        // Fill Matrices with results for this position (normalizing floats)
        data_all_changes[0][ps] = normalize(diff_res.extra.0).to_string();
        data_aa_changes[0][ps] = diff_res.extra.1.clone();
        data_aff_leaves[0][ps] = normalize(diff_res.extra.2).to_string();
        data_aff_branches[0][ps] = diff_res.extra.3.clone();

        for i in 0..diff_res.diff1.len() {
            let row_idx = i + 1;
            if row_idx < num_rows {
                data_all_changes[row_idx][ps] = normalize(diff_res.diff1[i]).to_string();
                data_aa_changes[row_idx][ps] = diff_res.mat2[i].clone();
                data_aff_leaves[row_idx][ps] = normalize(diff_res.eff_br[i]).to_string();
                data_aff_branches[row_idx][ps] = diff_res.aff_br[i].clone();
            }
        }
    }

    // 7. WRITE CSV
    write_final_csv(&out_all_changes, &data_all_changes, &tree_info, num_leaves)?;
    write_final_csv(&out_aa_changes, &data_aa_changes, &tree_info, num_leaves)?;
    write_final_csv(&out_aff_leaves, &data_aff_leaves, &tree_info, num_leaves)?;
    write_final_csv(&out_aff_branches, &data_aff_branches, &tree_info, num_leaves)?;

    Ok(())
}

// HELPER FUNCTIONS

fn aa_to_num(aa: char) -> u8 {
    // Map C/G to numeric IDs for probability indexing
    match aa { 'C' => 1, 'G' => 2, _ => 21 }
}

fn write_final_csv(filename: &str, data: &Vec<Vec<String>>, tree: &TreeInfo, _num_leaves: usize) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .quote_style(csv::QuoteStyle::Always)
        .from_path(filename)?;
    
    let mut header = vec!["Label".to_string()];
    for i in 1..=data[0].len() { header.push(format!("Pos_{}", i)); }
    header.push("BranchLength".to_string());
    wtr.write_record(&header)?;

    // Closure to write a single row mapped to a tree node
    let mut write_row = |i: usize| -> Result<()> {
        let row_idx = i + 1; 
        if row_idx < data.len() {
            let child_node = tree.node[i];
            let label = if child_node < tree.label.len() {
                tree.label[child_node].clone()
            } else {
                format!("Node{}", child_node)
            };
            
            let mut row = vec![label];
            row.extend(data[row_idx].clone());
            // Sanitize branch length too
            row.push(normalize(tree.branch_length[i]).to_string());
            wtr.write_record(&row)?;
        }
        Ok(())
    };

    // Sort outputs by Node ID for consistency
    let mut indices: Vec<usize> = (0..tree.parent.len()).collect();
    indices.sort_by_key(|&i| tree.node[i]);

    for i in indices {
        write_row(i)?;
    }

    wtr.flush()?;
    Ok(())
}

pub fn compute_difference(
    ps: usize,
    ord: &str,
    levels: &Vec<Vec<String>>, 
    tree_new_info: &TreeInfo,
    _num_nodes: usize,
    num_leaves: usize,
    x: &Vec<Vec<String>>,
    nodes_raxml_prev: &HashMap<String, usize>,
    msa_upd: &Vec<Vec<char>>,
    fasta_node: &Vec<Vec<char>>,
) -> DiffOutput {

    let num_branch = tree_new_info.parent.len(); 
    let root_node_id = num_leaves; 

    // Filter x_matrix for the current position (ps)
    let ps_str = (ps + 1).to_string();
    let mut tt_rows: Vec<&Vec<String>> = x.iter().filter(|row| row[1] == ps_str).collect();

    tt_rows.sort_by(|a, b| {
        let a_n = a[0].parse::<f64>().unwrap_or(0.0);
        let b_n = b[0].parse::<f64>().unwrap_or(0.0);
        a_n.partial_cmp(&b_n).unwrap()
    });

    // Populate node probabilities map
    let mut node_probs: HashMap<usize, Vec<f64>> = HashMap::new();
    let ord_suffix = ord.replace("Node", "");
    let mut _root_prob = vec![0.0, 0.0];
    
    for row in &tt_rows {
        if let Ok(nid) = row[0].parse::<usize>() {
            // [MSA2] x_matrix cols 3 & 4 correspond to p_C and p_G probabilities
            let p1 = row[3].parse::<f64>().unwrap_or(0.0);
            let p2 = row[4].parse::<f64>().unwrap_or(0.0);
            node_probs.insert(nid, vec![p1, p2]);
            if row[0] == ord_suffix { _root_prob = vec![p1, p2]; }
        }
    }

    let max_id = tree_new_info.label.len();
    let mut all_prob = vec![vec![0.0; 2]; max_id];

    // Set leaf probabilities (Tips are hard 1.0 or 0.0)
    for i in 0..num_leaves {
        if ps < msa_upd[i].len() {
            let num = aa_to_num(msa_upd[i][ps]);
            if num == 1 { all_prob[i][0] = 1.0; } // C
            else if num == 2 { all_prob[i][1] = 1.0; } // G
        }
    }
    
    // Fill internal node probabilities from the map
    for (nid, probs) in node_probs {
        let label_key = format!("Node{}", nid);
        if let Some(&idx) = nodes_raxml_prev.get(&label_key) {
             if idx < all_prob.len() { all_prob[idx] = probs; }
        }
    }

    // Dynamic Root Logic: determine dominant state (C or G) at this site
    let mut counts = HashMap::new();
    for i in 0..num_leaves {
        if ps < msa_upd[i].len() {
            *counts.entry(msa_upd[i][ps]).or_insert(0) += 1;
        }
    }
    
    let c_count = *counts.get(&'C').unwrap_or(&0);
    let g_count = *counts.get(&'G').unwrap_or(&0);

    let (faa, fa) = if g_count > c_count {
        ('G', vec![0.0, 1.0])
    } else {
        ('C', vec![1.0, 0.0])
    };

    // Initialize vectors for difference calculations
    let mut diff = Vec::with_capacity(num_branch + 1);
    let mut mat2 = Vec::with_capacity(num_branch + 1);
    let mut num_of_ch = Vec::with_capacity(num_branch + 1);
    let mut aff_brs = Vec::with_capacity(num_branch + 1);

    // Dummy Root Calculation
    let root_right = root_node_id; 
    let d_vec_root = vec![all_prob[root_right][0] - fa[0], all_prob[root_right][1] - fa[1]];
    // Sanitize diff values
    let diff_root: f64 = normalize(d_vec_root.iter().filter(|&&v| v > 0.0).sum());
    let root_aa = if root_right < fasta_node.len() { fasta_node[root_right][ps] } else { '-' };
    let m1_root = faa;
    let m2_root = root_aa;
    
    if m1_root == m2_root {
        mat2.push("--".to_string());
        num_of_ch.push(0.0);
        aff_brs.push("0".to_string());
    } else {
        mat2.push(format!("{}{}", m1_root, m2_root));
        num_of_ch.push(1.0);
        aff_brs.push("Node0".to_string());
    }
    diff.push(diff_root);

    // Calculate Real Branch Differences (Parent vs Child)
    for i in 0..num_branch {
        let lf = tree_new_info.parent[i];
        let rg = tree_new_info.node[i];
        
        let d_vec = vec![all_prob[rg][0] - all_prob[lf][0], all_prob[rg][1] - all_prob[lf][1]];
        // Sum only positive differences to get net change
        let val_diff: f64 = normalize(d_vec.iter().filter(|&&v| v > 0.0).sum());
        diff.push(val_diff);

        // Determine Amino Acid / State change strings
        let aa_lf = if lf < num_leaves { msa_upd[lf][ps] } else { 
            if lf < fasta_node.len() { fasta_node[lf][ps] } else { '-' }
        };
        let aa_rg = if rg < num_leaves { msa_upd[rg][ps] } else {
            if rg < fasta_node.len() { fasta_node[rg][ps] } else { '-' }
        };

        if aa_lf == aa_rg {
            mat2.push("--".to_string());
            num_of_ch.push(0.0);
            aff_brs.push("0".to_string());
        } else {
            mat2.push(format!("{}{}", aa_lf, aa_rg));
            num_of_ch.push(1.0);
            let lbl = if rg < tree_new_info.label.len() { 
                tree_new_info.label[rg].clone() 
            } else { 
                format!("Node{}", rg) 
            };
            aff_brs.push(lbl);
        }
    }

    // Back Propagation Logic: Propagate changes up the tree to find effective branches
    let mut child_to_branch: HashMap<usize, usize> = HashMap::new();
    child_to_branch.insert(root_right, 0); 
    for i in 0..num_branch {
        child_to_branch.insert(tree_new_info.node[i], i + 1);
    }

    for level_grp in levels.iter() {
        for pair_str in level_grp {
            if pair_str == "0" { continue; }
            let parts: Vec<&str> = pair_str.split(',').collect();
            if parts.len() < 2 { continue; }
            let rg: usize = parts[0].parse().unwrap_or(0); // Parent
            let lf: usize = parts[1].parse().unwrap_or(0); // Child

            if let Some(&lc) = child_to_branch.get(&lf) {
                // If the current branch (lc) has no state change ("--"), accumulate info from below
                if mat2[lc] == "--" {
                    let keep = normalize(num_of_ch[lc] + 1.0);
                    
                    let child_label = if lf < tree_new_info.label.len() {
                        tree_new_info.label[lf].clone()
                    } else {
                        format!("Node{}", lf)
                    };
                    
                    let keep_lf = if aff_brs[lc] == "0" {
                         format!("{}-0", child_label)
                    } else {
                         format!("{}-{}", child_label, aff_brs[lc])
                    };

                    // Clear current branch info as it's passing through
                    num_of_ch[lc] = 0.0;
                    aff_brs[lc] = "0".to_string();

                    // Add to Parent's accumulator
                    if let Some(&prev) = child_to_branch.get(&rg) {
                        num_of_ch[prev] = normalize(num_of_ch[prev] + keep);
                        
                        let prev_str = aff_brs[prev].clone();
                        if prev_str == "0" {
                             aff_brs[prev] = format!("0-{}", keep_lf);
                        } else {
                             aff_brs[prev] = format!("{}-{}", prev_str, keep_lf); 
                        }
                    }
                }
            }
        }
    }

    let final_diff: Vec<f64> = diff[1..].iter().map(|&x| normalize(x)).collect();
    let final_eff: Vec<f64> = num_of_ch[1..].iter().map(|&x| normalize(x)).collect();

    DiffOutput {
        diff1: final_diff,
        mat2: mat2[1..].to_vec(),
        eff_br: final_eff,
        aff_br: aff_brs[1..].to_vec(),
        extra: (normalize(diff[0]), mat2[0].clone(), normalize(num_of_ch[0]), aff_brs[0].clone()),
    }
}

// CUSTOM NEWICK PARSER

fn process_newick_custom(path: &str) -> Result<(TreeInfo, usize, usize, HashMap<String, usize>)> {
    let content = fs::read_to_string(path)?;
    let root = parse_newick_string(&content).context("Parsing error")?;

    let mut tip_labels = Vec::new();
    collect_tips_custom(&root, &mut tip_labels);

    let num_leaves = tip_labels.len();
    let mut label_to_id = HashMap::new();
    for (i, lbl) in tip_labels.iter().enumerate() {
        label_to_id.insert(lbl.clone(), i);
    }

    let mut counter = num_leaves;
    let mut parents = Vec::new();
    let mut children = Vec::new();
    let mut lengths = Vec::new();
    let mut all_labels = vec![String::new(); num_leaves]; 
    for (i, l) in tip_labels.iter().enumerate() { all_labels[i] = l.clone(); }

    flatten_recursive_custom(&root, &mut counter, &mut label_to_id, &mut parents, &mut children, &mut lengths, &mut all_labels)?;

    let num_nodes = counter - num_leaves;
    let info = TreeInfo { parent: parents, node: children, label: all_labels, branch_length: lengths };
    Ok((info, num_leaves, num_nodes, label_to_id))
}

fn parse_newick_string(input: &str) -> Result<SimpleNode> {
    let clean = input.trim().trim_end_matches(';');
    
    // Recursive parser for nested brackets
    fn parse_subtree(s: &str) -> Result<SimpleNode> {
        let s = s.trim();
        let mut children = Vec::new();
        let mut remainder = s;

        if s.starts_with('(') {
            let mut balance = 0;
            let mut split_idx = 0;
            // Find matching closing parenthesis
            for (i, c) in s.chars().enumerate() {
                if c == '(' { balance += 1; }
                else if c == ')' { balance -= 1; }
                if balance == 0 {
                    split_idx = i;
                    break;
                }
            }
            
            let inner = &s[1..split_idx];
            remainder = &s[split_idx+1..]; 

            let mut start = 0;
            let mut depth = 0;
            // Split children by comma
            for (i, c) in inner.chars().enumerate() {
                if c == '(' { depth += 1; }
                else if c == ')' { depth -= 1; }
                else if c == ',' && depth == 0 {
                    children.push(parse_subtree(&inner[start..i])?);
                    start = i + 1;
                }
            }
            children.push(parse_subtree(&inner[start..])?);
        }

        // Parse Label and Branch Length (after :)
        let (label_part, len_part) = if let Some(idx) = remainder.find(':') {
            (&remainder[..idx], Some(&remainder[idx+1..]))
        } else {
            (remainder, None)
        };

        let label = label_part.trim().to_string();
        let len = if let Some(l) = len_part {
            l.parse::<f64>().unwrap_or(0.0)
        } else {
            0.0
        };

        Ok(SimpleNode { label, len, children })
    }

    parse_subtree(clean)
}

fn collect_tips_custom(node: &SimpleNode, tips: &mut Vec<String>) {
    if node.children.is_empty() {
        tips.push(node.label.clone());
    } else {
        for child in &node.children { collect_tips_custom(child, tips); }
    }
}

// Convert recursive tree structure into flat arrays
fn flatten_recursive_custom(
    node: &SimpleNode, 
    counter: &mut usize, 
    map: &mut HashMap<String, usize>, 
    p: &mut Vec<usize>, 
    c: &mut Vec<usize>, 
    l: &mut Vec<f64>, 
    lbls: &mut Vec<String>
) -> Result<usize> {
    let curr_id;
    if node.children.is_empty() {
        if let Some(&id) = map.get(&node.label) {
            curr_id = id;
        } else {
            return Err(anyhow!("Error: Tip label '{}' found in tree topology was not indexed in the tip collection phase.", node.label));
        }
    } else {
        curr_id = *counter;
        *counter += 1;
        
        let name = if node.label.is_empty() { format!("Node{}", curr_id) } else { node.label.clone() };
        map.insert(name.clone(), curr_id);
        
        if lbls.len() <= curr_id { lbls.resize(curr_id + 1, String::new()); }
        lbls[curr_id] = name;

        for child in &node.children {
            let child_id = flatten_recursive_custom(child, counter, map, p, c, l, lbls)?;
            p.push(curr_id);
            c.push(child_id);
            l.push(child.len);
        }
    }
    Ok(curr_id)
}

// Determine tree traversal levels for bottom-up/top-down operations
fn calculate_levels_exact(tree: &TreeInfo, _num_leaves: usize) -> Vec<Vec<String>> {
    let mut levels = Vec::new();
    let mut active_edges: HashSet<usize> = (0..tree.parent.len()).collect();
    
    loop {
        let mut remaining_parents: HashSet<usize> = HashSet::new();
        for &idx in &active_edges {
            remaining_parents.insert(tree.parent[idx]);
        }
        
        // Find tips in the current pruned tree (nodes that are children but not parents)
        let mut current_tips_indices = Vec::new();
        let mut sorted_active: Vec<usize> = active_edges.iter().cloned().collect();
        sorted_active.sort();

        for idx in sorted_active {
            let child = tree.node[idx];
            if !remaining_parents.contains(&child) {
                current_tips_indices.push(idx);
            }
        }
        
        if current_tips_indices.is_empty() {
            break;
        }

        let mut level_group = Vec::new();
        for &idx in &current_tips_indices {
            let p = tree.parent[idx];
            let c = tree.node[idx];
            level_group.push(format!("{},{}", p, c)); // PARENT, CHILD
            active_edges.remove(&idx);
        }
        levels.push(level_group);
    }
    
    levels
}