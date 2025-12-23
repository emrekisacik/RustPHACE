use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

// --- Configuration Constants ---
// Used for mapping characters to array indices
const AMINO_ACIDS: [&str; 20] = [
    "G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", 
    "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T"
];

// --- Data Structures ---

#[derive(Debug, Clone)]
struct TreeNode {
    label: String,
    parent_id: Option<usize>, // None for the root
    branch_length: f64,
    children: Vec<usize>,     // IDs of child nodes
}

struct ParsedTree {
    nodes: HashMap<usize, TreeNode>,
    leaves: Vec<usize>,   // IDs of leaf nodes
    internal: Vec<usize>, // IDs of internal nodes
    root_id: usize,       // Specific ID of the root
}

// --- Helper Functions ---

fn aa_to_num(aa: &str) -> usize {
    AMINO_ACIDS.iter().position(|&r| r == aa).unwrap_or(20)
}

// --- Parsers ---

/// Reads a FASTA file into a HashMap where Key = Header, Value = Vector of characters
fn read_fasta(path: &str) -> io::Result<HashMap<String, Vec<String>>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut msa = HashMap::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() { continue; }

        if trimmed.starts_with('>') {
            // If we hit a new header, save the previous sequence (if it exists)
            if !current_header.is_empty() {
                let seq_vec: Vec<String> = current_seq.chars().map(|c| c.to_string()).collect();
                msa.insert(current_header, seq_vec);
            }
            // Parse new header (remove '>')
            current_header = trimmed[1..].trim().to_string();
            current_seq = String::new();
        } else {
            // Append sequence lines (handles multiline FASTA)
            current_seq.push_str(trimmed);
        }
    }
    // Save the very last sequence in the file
    if !current_header.is_empty() {
        let seq_vec: Vec<String> = current_seq.chars().map(|c| c.to_string()).collect();
        msa.insert(current_header, seq_vec);
    }
    Ok(msa)
}

/// Parses a Newick format tree file manually to match the logic of R's 'ape' package
fn read_tree_structure(path: &str) -> io::Result<ParsedTree> {
    let file = File::open(path)?;
    let mut content = String::new();
    BufReader::new(file).read_line(&mut content)?;
    let content = content.trim().trim_end_matches(';');

    // 1. Tokenize the Newick string
    // Splits input into: '(', ')', ',', ':', and labels/lengths
    let mut tokens = Vec::new();
    let mut buffer = String::new();
    let chars: Vec<char> = content.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        match chars[i] {
            '(' | ')' | ',' | ':' => {
                if !buffer.is_empty() {
                    tokens.push(buffer.clone());
                    buffer.clear();
                }
                tokens.push(chars[i].to_string());
            }
            '\'' => {
                // Handle quoted labels
                buffer.push('\'');
                i += 1;
                while i < chars.len() && chars[i] != '\'' {
                    buffer.push(chars[i]);
                    i += 1;
                }
                buffer.push('\'');
            }
            _ => buffer.push(chars[i]),
        }
        i += 1;
    }
    if !buffer.is_empty() { tokens.push(buffer); }

    // Temporary structure for parsing before we re-index IDs
    struct TmpNode {
        lbl: String,
        len: f64,
        children: Vec<usize>,
        is_leaf: bool,
    }
    let mut tmp_nodes = Vec::new();
    
    // 2. Recursive Parsing Function
    // Consumes tokens and builds the TmpNode hierarchy
    fn parse_subtree(
        tokens: &[String], 
        pos: &mut usize, 
        nodes: &mut Vec<TmpNode>
    ) -> usize {
        let start_token = &tokens[*pos];
        let my_idx = nodes.len();
        nodes.push(TmpNode { lbl: String::new(), len: 0.0, children: vec![], is_leaf: false }); 

        if start_token == "(" {
            // Internal Node: Recurse for children
            *pos += 1; 
            loop {
                let child_idx = parse_subtree(tokens, pos, nodes);
                nodes[my_idx].children.push(child_idx);
                
                if *pos >= tokens.len() { break; }
                let t = &tokens[*pos];
                if t == "," {
                    *pos += 1; 
                } else if t == ")" {
                    *pos += 1; 
                    break;
                } else {
                    break;
                }
            }
            nodes[my_idx].is_leaf = false;
        } else {
            // Leaf Node
            nodes[my_idx].lbl = start_token.clone();
            nodes[my_idx].is_leaf = true;
            *pos += 1;
        }

        // Check for Label after closing parenthesis or leaf name
        if *pos < tokens.len() {
            let t = &tokens[*pos];
            if t != ":" && t != "," && t != ")" && t != ";" {
                nodes[my_idx].lbl = t.clone();
                *pos += 1;
            }
        }
        // Check for Branch Length
        if *pos < tokens.len() && tokens[*pos] == ":" {
            *pos += 1;
            if *pos < tokens.len() {
                nodes[my_idx].len = tokens[*pos].parse::<f64>().unwrap_or(0.0);
                *pos += 1;
            }
        }
        my_idx
    }

    let mut pos = 0;
    let root_tmp_index = parse_subtree(&tokens, &mut pos, &mut tmp_nodes);

    // 3. Re-indexing
    // Order: Leaves (1..N) -> Root (N+1) -> Other Internals (N+2..)
    fn traverse_collect(
        curr: usize, 
        nodes: &Vec<TmpNode>, 
        leaves: &mut Vec<usize>, 
        internals: &mut Vec<usize>
    ) {
        if nodes[curr].is_leaf {
            leaves.push(curr);
        } else {
            for &c in &nodes[curr].children {
                traverse_collect(c, nodes, leaves, internals);
            }
            internals.push(curr);
        }
    }

    let mut ordered_leaves = Vec::new();
    let mut ordered_internals_post = Vec::new(); 
    traverse_collect(root_tmp_index, &tmp_nodes, &mut ordered_leaves, &mut ordered_internals_post);

    let num_leaves = ordered_leaves.len();
    let mut final_nodes = HashMap::new();
    let mut new_id_map = HashMap::new(); 

    // Assign IDs: Leaves first
    for (i, &tmp_idx) in ordered_leaves.iter().enumerate() {
        new_id_map.insert(tmp_idx, i + 1);
    }
    // Assign Root ID
    let final_root_id = num_leaves + 1;
    new_id_map.insert(root_tmp_index, final_root_id);

    // Assign remaining Internal IDs
    let mut next_internal = num_leaves + 2;
    for &tmp_idx in &ordered_internals_post {
        if tmp_idx != root_tmp_index {
            new_id_map.insert(tmp_idx, next_internal);
            next_internal += 1;
        }
    }

    let mut leaves_final = Vec::new();
    let mut internals_final = Vec::new();

    // 4. Construct Final Tree
    for (tmp_idx, node) in tmp_nodes.iter().enumerate() {
        if !new_id_map.contains_key(&tmp_idx) { continue; }

        let my_id = new_id_map[&tmp_idx];
        let child_ids: Vec<usize> = node.children.iter()
            .map(|c| *new_id_map.get(c).unwrap())
            .collect();
        
        let final_node = TreeNode {
            label: node.lbl.clone(),
            parent_id: None, // Will be filled in the next step
            branch_length: node.len,
            children: child_ids,
        };
        
        final_nodes.insert(my_id, final_node);
        if node.is_leaf { leaves_final.push(my_id); } else { internals_final.push(my_id); }
    }

    // 5. Back-fill Parent IDs
    // Iterate all nodes, look at their children, and tell the children who their parent is
    let keys: Vec<usize> = final_nodes.keys().cloned().collect();
    for pid in keys {
        let children = final_nodes[&pid].children.clone();
        for cid in children {
            if let Some(child) = final_nodes.get_mut(&cid) {
                child.parent_id = Some(pid);
            }
        }
    }

    leaves_final.sort();
    internals_final.sort();

    Ok(ParsedTree {
        nodes: final_nodes,
        leaves: leaves_final,
        internal: internals_final,
        root_id: final_root_id,
    })
}

/// Reads the ancestral state reconstruction file (e.g., from IQ-TREE)
/// Returns a Map of (NodeName, Site) -> [Probabilities of 20 AAs]
fn read_state_file(path: &str) -> io::Result<(HashMap<(String, usize), Vec<f64>>, usize)> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Parse Header to find which columns correspond to which Amino Acids
    let header_line = lines.next().ok_or(io::Error::new(io::ErrorKind::InvalidData, "Empty state file"))??;
    let headers: Vec<String> = header_line.split('\t').map(|s| s.trim().to_string()).collect();

    let mut col_map = HashMap::new(); 
    let mut node_col_idx = 0;
    let mut site_col_idx = 1;

    for (i, h) in headers.iter().enumerate() {
        if h == "Node" { node_col_idx = i; }
        else if h == "Site" { site_col_idx = i; }
        else if h.starts_with("p_") {
            // e.g., "p_A", "p_C"
            let aa_str = h.trim_start_matches("p_");
            let aa_idx = aa_to_num(aa_str);
            if aa_idx < 20 {
                col_map.insert(i, aa_idx);
            }
        }
    }

    let mut data = HashMap::new();
    let mut max_site = 0;

    // Parse Data Rows
    for line in lines {
        let line = line?;
        if line.trim().is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        
        if cols.len() <= site_col_idx || cols.len() <= node_col_idx { continue; }

        let raw_node = cols[node_col_idx].trim();
        let node_label = raw_node.replace("Node", ""); 
        
        let site_str = cols[site_col_idx].trim();
        let site: usize = match site_str.parse() {
            Ok(s) => s,
            Err(_) => continue, 
        };
        if site > max_site { max_site = site; }

        let mut probs = vec![0.0; 20];
        for (&col_i, &aa_i) in &col_map {
            if col_i < cols.len() {
                let val: f64 = cols[col_i].trim().parse().unwrap_or(0.0);
                probs[aa_i] = val;
            }
        }
        data.insert((node_label, site), probs);
    }
    Ok((data, max_site))
}

// --- Logic ---

/// The Core Algorithm: Calculates tolerance scores for a specific position
/// Compares ancestral probabilities between Parent and Child nodes, weighted by distance
fn calculate_position_score(
    pos: usize,
    tree: &ParsedTree,
    state_data: &HashMap<(String, usize), Vec<f64>>,
    msa: &HashMap<String, Vec<String>>,
    weights: &HashMap<usize, f64>
) -> Vec<f64> {

    let num_leaves = tree.leaves.len();
    let num_nodes = tree.internal.len();
    let root_id = tree.root_id;

    // 1. Prepare Matrix Probabilities (Ancestral States for Internal Nodes)
    let mut matrix_prob: HashMap<usize, Vec<f64>> = HashMap::new();
    for &nid in &tree.internal {
        let label_raw = &tree.nodes[&nid].label;
        let label_clean = label_raw.replace("Node", "");
        let mut p = match state_data.get(&(label_clean, pos)) {
            Some(v) => v.clone(),
            None => vec![0.0; 20], 
        };
        // Apply threshold to reduce noise
        for val in &mut p { if *val < 0.01 { *val = 0.0; } }
        matrix_prob.insert(nid, p);
    }

    // 2. Prepare Leaf Vectors (One-hot encoding of observed MSA data)
    let mut prob_leaves: HashMap<usize, Vec<f64>> = HashMap::new();
    let mut gaps: HashSet<usize> = HashSet::new(); 
    for &lid in &tree.leaves {
        let label = &tree.nodes[&lid].label;
        let mut is_gap = true;
        let mut v = vec![0.0; 20];
        if let Some(seq) = msa.get(label) {
            // Convert 1-based pos to 0-based index
            if pos <= seq.len() && pos > 0 {
                let aa_char = &seq[pos - 1]; 
                let aa_idx = aa_to_num(aa_char);
                if aa_idx < 20 {
                    v[aa_idx] = 1.0; // Certainty for observed AA
                    is_gap = false;
                }
            }
        }
        if is_gap { gaps.insert(lid); }
        prob_leaves.insert(lid, v);
    }

    // 3. Calculation Loop
    // Iterates over every Amino Acid to calculate a score for that AA at this position
    let mut final_scores = vec![0.0; 20];
    let root_pr = matrix_prob.get(&root_id).cloned().unwrap_or(vec![0.0; 20]);

    for aa_i in 0..20 {
        let mut score = 0.0;

        // Part A: Internal Nodes Contribution
        if num_nodes > 1 {
            for &child_id in &tree.internal {
                if child_id == root_id { continue; }
                
                if let Some(parent_id) = tree.nodes[&child_id].parent_id {
                    // Compare Parent Prob vs Child Prob
                    if let Some(parent_probs) = matrix_prob.get(&parent_id) {
                        let child_probs = matrix_prob.get(&child_id).unwrap();
                        let pp = parent_probs[aa_i];
                        let cp = child_probs[aa_i];
                        

                        let mut diff = pp - cp;
                        if diff < 0.0 { diff = 0.0; }

                        // Weight by the Parent's distance factor
                        let w = weights.get(&parent_id).unwrap_or(&0.0);
                        score += w * diff;
                    }
                }
            }
        }

        // Part B: Root Contribution
        score += root_pr[aa_i];

        // Part C: Leaves Contribution
        let mut s1 = 0.0;
        for &lid in &tree.leaves {
            if gaps.contains(&lid) { continue; } // Skip gaps
            if let Some(parent_id) = tree.nodes[&lid].parent_id {
                if let Some(parent_probs) = matrix_prob.get(&parent_id) {
                    let leaf_val = prob_leaves[&lid][aa_i]; // 1.0 or 0.0
                    let parent_val = parent_probs[aa_i];
                    
                    let mut diff = leaf_val - parent_val;
                    if diff < 0.0 { diff = 0.0; }

                    let w = weights.get(&lid).unwrap_or(&0.0);
                    s1 += w * diff;
                }
            }
        }
        score += s1;

        // Part D: Normalize by total number of nodes
        final_scores[aa_i] = score / ((num_nodes + num_leaves) as f64);
    }

    final_scores
}

// --- Entry Point ---

pub fn tolerance(id: &str) -> Result<(), Box<dyn Error>> {
    // Define input/output filenames
    let file_fasta = format!("{}_MaskedMSA.fasta", id);
    let file_nwk = format!("{}.treefile", id);
    let file_rst = format!("{}.state", id);
    let output_path = format!("ToleranceScores/{}.csv", id);

    // Read Data
    let tree = read_tree_structure(&file_nwk)?; 
    let msa = read_fasta(&file_fasta)?;
    let (state_data, total_pos) = read_state_file(&file_rst)?;

    // --- Weights Calculation ---
    // Calculates a weight for each node based on its distance from the root
    // Uses a Gaussian-like kernel: exp(-dist^2 / mean_dist^2)
    let mut dists: HashMap<usize, f64> = HashMap::new();
    let mut stack = vec![(tree.root_id, 0.0)];
    
    // DFS to calculate cumulative branch lengths from root
    while let Some((curr, d)) = stack.pop() {
        dists.insert(curr, d);
        let node = &tree.nodes[&curr];
        for &child in &node.children {
            let child_len = tree.nodes[&child].branch_length;
            stack.push((child, d + child_len));
        }
    }

    let dist_values: Vec<f64> = dists.values().cloned().collect();
    if dist_values.is_empty() { return Err("Tree invalid.".into()); }
    
    let sum_dist: f64 = dist_values.iter().sum();
    let mean_dist = sum_dist / (dist_values.len() as f64);
    
    let mut weights = HashMap::new();
    for (id, d) in dists {
        let w = (- (d.powi(2)) / mean_dist.powi(2)).exp();
        weights.insert(id, w);
    }

    // --- Output Generation ---
    let _ = std::fs::create_dir_all("ToleranceScores");
    let file = File::create(&output_path)?;
    let mut writer = std::io::BufWriter::new(file);

    // Write Header
    write!(writer, "Pos/AA")?;
    for aa in AMINO_ACIDS.iter() { write!(writer, ",{}", aa)?; }
    writeln!(writer)?;

    // Process every site and write rows
    for ps in 1..=total_pos {
        let scores = calculate_position_score(ps, &tree, &state_data, &msa, &weights);
        write!(writer, "{}", ps)?;
        for s in scores { write!(writer, ",{}", s)?; }
        writeln!(writer)?;
    }
    
    Ok(())
}