use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

const AMINO_ACIDS: [&str; 20] = [
    "G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", 
    "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T"
];

fn canonicalize_id(label: &str) -> String {
    let mut s = label.replace("'", "").replace("\"", "");
    s = s.replace("Node", "").replace("Branch", "");
    let trimmed = s.trim();
    if let Ok(num) = trimmed.parse::<usize>() {
        return num.to_string();
    }
    trimmed.to_string()
}

#[derive(Debug, Clone)]
struct TreeNode {
    raw_label: String,
    clean_id: String,
    parent_id: Option<usize>,
    branch_length: f64,
    children: Vec<usize>,
}

struct ParsedTree {
    nodes: HashMap<usize, TreeNode>,
    leaves: Vec<usize>,
    internal: Vec<usize>,
    root_id: usize,
}

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
            if !current_header.is_empty() {
                let seq: Vec<String> = current_seq.chars().map(|c| c.to_string()).collect();
                msa.insert(current_header, seq);
            }
            current_header = trimmed[1..].trim().to_string();
            current_seq = String::new();
        } else {
            current_seq.push_str(trimmed);
        }
    }
    if !current_header.is_empty() {
        let seq: Vec<String> = current_seq.chars().map(|c| c.to_string()).collect();
        msa.insert(current_header, seq);
    }
    Ok(msa)
}

fn read_tree_structure(path: &str) -> io::Result<ParsedTree> {
    let file = File::open(path)?;
    let mut content = String::new();
    BufReader::new(file).read_line(&mut content)?;
    let content = content.trim().trim_end_matches(';');

    let mut tokens = Vec::new();
    let mut buffer = String::new();
    let chars: Vec<char> = content.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        match chars[i] {
            '(' | ')' | ',' | ':' => {
                let token = buffer.trim();
                if !token.is_empty() { tokens.push(token.to_string()); }
                buffer.clear();
                tokens.push(chars[i].to_string());
            }
            '\'' => {
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
    if !buffer.trim().is_empty() { tokens.push(buffer.trim().to_string()); }

    struct TmpNode { lbl: String, len: f64, children: Vec<usize>, is_leaf: bool }
    let mut tmp_nodes = Vec::new();

    fn parse_recur(tokens: &[String], pos: &mut usize, nodes: &mut Vec<TmpNode>) -> usize {
        let my_idx = nodes.len();
        nodes.push(TmpNode { lbl: String::new(), len: 0.0, children: vec![], is_leaf: false });
        
        if tokens[*pos] == "(" {
            *pos += 1;
            loop {
                let child = parse_recur(tokens, pos, nodes);
                nodes[my_idx].children.push(child);
                if *pos >= tokens.len() { break; }
                if tokens[*pos] == "," { *pos += 1; }
                else if tokens[*pos] == ")" { *pos += 1; break; }
                else { break; }
            }
            nodes[my_idx].is_leaf = false;
        } else {
            nodes[my_idx].lbl = tokens[*pos].clone();
            nodes[my_idx].is_leaf = true;
            *pos += 1;
        }

        if *pos < tokens.len() {
            let t = &tokens[*pos];
            if t != ":" && t != "," && t != ")" && t != ";" {
                nodes[my_idx].lbl = t.clone();
                *pos += 1;
            }
        }
        if *pos < tokens.len() && tokens[*pos] == ":" {
            *pos += 1;
            if *pos < tokens.len() {
                let len_str = tokens[*pos].replace(")", "").replace(",", ""); 
                nodes[my_idx].len = len_str.parse::<f64>().unwrap_or(0.0);
                *pos += 1;
            }
        }
        my_idx
    }

    let mut pos = 0;
    let root_tmp = parse_recur(&tokens, &mut pos, &mut tmp_nodes);

    fn collect_indices(curr: usize, nodes: &Vec<TmpNode>, leaves: &mut Vec<usize>, internals: &mut Vec<usize>) {
        if nodes[curr].is_leaf {
            leaves.push(curr);
        } else {
            for &c in &nodes[curr].children { collect_indices(c, nodes, leaves, internals); }
            internals.push(curr);
        }
    }
    let mut ordered_leaves = Vec::new();
    let mut ordered_internals = Vec::new();
    collect_indices(root_tmp, &tmp_nodes, &mut ordered_leaves, &mut ordered_internals);

    let num_leaves = ordered_leaves.len();
    let mut map_old_new = HashMap::new();
    let mut final_nodes = HashMap::new();
    let mut final_leaves = Vec::new();
    let mut final_internals = Vec::new();

    for (i, &old) in ordered_leaves.iter().enumerate() {
        let new_id = i + 1;
        map_old_new.insert(old, new_id);
        final_leaves.push(new_id);
    }
    let root_final_id = num_leaves + 1;
    map_old_new.insert(root_tmp, root_final_id);
    
    let mut next_internal = num_leaves + 2;
    for &old in &ordered_internals {
        if old != root_tmp {
            map_old_new.insert(old, next_internal);
            next_internal += 1;
        }
    }
    final_internals.push(root_final_id);
    for &old in &ordered_internals {
        if old != root_tmp { final_internals.push(map_old_new[&old]); }
    }
    final_internals.sort();

    for (old, new_id) in &map_old_new {
        let tmp = &tmp_nodes[*old];
        let children: Vec<usize> = tmp.children.iter().map(|c| map_old_new[c]).collect();
        let clean = canonicalize_id(&tmp.lbl);
        
        final_nodes.insert(*new_id, TreeNode {
            raw_label: tmp.lbl.clone(),
            clean_id: clean,
            parent_id: None,
            branch_length: tmp.len,
            children,
        });
    }

    let ids: Vec<usize> = final_nodes.keys().cloned().collect();
    for pid in ids {
        let children = final_nodes[&pid].children.clone();
        for cid in children {
            if let Some(child) = final_nodes.get_mut(&cid) {
                child.parent_id = Some(pid);
            }
        }
    }

    Ok(ParsedTree { nodes: final_nodes, leaves: final_leaves, internal: final_internals, root_id: root_final_id })
}

fn read_state_file(path: &str) -> io::Result<(HashMap<(String, usize), Vec<f64>>, usize)> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let mut header_line = String::new();
    
    while let Some(line_res) = lines.next() {
        let line = line_res?;
        let trimmed = line.trim();
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            header_line = line;
            break;
        }
    }

    if header_line.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "State file contains no valid header"));
    }

    let delimiter = if header_line.contains('\t') { '\t' } else { ' ' };
    let headers: Vec<String> = header_line.split(delimiter)
        .map(|s| s.trim().replace("\"", "").replace("'", "")) 
        .filter(|s| !s.is_empty())
        .collect();

    let mut col_map = HashMap::new(); 
    let mut node_col_idx = None;
    let mut site_col_idx = None;

    for (i, h) in headers.iter().enumerate() {
        let h_upper = h.to_uppercase();
        
        if h_upper == "NODE" { 
            node_col_idx = Some(i); 
        } else if h_upper == "SITE" { 
            site_col_idx = Some(i); 
        } else {
            let aa_clean = h_upper.replace("P_", "").trim().to_string();
            if let Some(pos) = AMINO_ACIDS.iter().position(|&a| a == aa_clean) {
                col_map.insert(i, pos);
            }
        }
    }

    if col_map.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "No AA columns found in header"));
    }
    
    let node_idx = node_col_idx.ok_or(io::Error::new(io::ErrorKind::InvalidData, "Missing 'Node' column"))?;
    let site_idx = site_col_idx.ok_or(io::Error::new(io::ErrorKind::InvalidData, "Missing 'Site' column"))?;

    let mut data = HashMap::new();
    let mut max_site = 0;

    for line in lines {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') { continue; }
        
        let cols: Vec<&str> = line.split(delimiter)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .collect();
        
        if cols.len() <= site_idx || cols.len() <= node_idx { continue; }

        let raw_node = cols[node_idx];
        let clean_node_id = canonicalize_id(raw_node); 
        
        let site_str = cols[site_idx];
        let site: usize = match site_str.parse() {
            Ok(s) => s,
            Err(_) => continue, 
        };
        if site > max_site { max_site = site; }

        let mut probs = vec![0.0; 20];
        for (&col_i, &aa_i) in &col_map {
            if col_i < cols.len() {
                probs[aa_i] = cols[col_i].parse().unwrap_or(0.0);
            }
        }
        data.insert((clean_node_id, site), probs);
    }
    
    Ok((data, max_site))
}

fn calculate_score(
    pos: usize, 
    tree: &ParsedTree, 
    state: &HashMap<(String, usize), Vec<f64>>,
    msa: &HashMap<String, Vec<String>>,
    weights: &HashMap<usize, f64>
) -> Vec<f64> {
    
    let mut node_probs: HashMap<usize, Vec<f64>> = HashMap::new();
    
    for &nid in &tree.internal {
        let clean = &tree.nodes[&nid].clean_id;
        let mut p = match state.get(&(clean.clone(), pos)) {
            Some(v) => v.clone(),
            None => vec![0.0; 20],
        };
        for x in &mut p { if *x < 0.01 { *x = 0.0; } }
        node_probs.insert(nid, p);
    }

    let mut leaf_probs: HashMap<usize, Vec<f64>> = HashMap::new();
    let mut gaps = HashSet::new();

    for &lid in &tree.leaves {
        let label = &tree.nodes[&lid].raw_label;
        let mut v = vec![0.0; 20];
        let mut is_gap = true;
        
        if let Some(seq) = msa.get(label) {
            if pos <= seq.len() && pos > 0 {
                let aa = seq[pos-1].to_uppercase();
                if let Some(idx) = AMINO_ACIDS.iter().position(|&a| a == aa) {
                    v[idx] = 1.0;
                    is_gap = false;
                }
            }
        }
        if is_gap { gaps.insert(lid); }
        leaf_probs.insert(lid, v);
    }

    let mut scores = vec![0.0; 20];
    let root_p = node_probs.get(&tree.root_id).cloned().unwrap_or(vec![0.0; 20]);

    for aa_i in 0..20 {
        let mut score = 0.0;

        for &child_id in &tree.internal {
            if child_id == tree.root_id { continue; }
            if let Some(pid) = tree.nodes[&child_id].parent_id {
                if let Some(pp) = node_probs.get(&pid) {
                    if let Some(cp) = node_probs.get(&child_id) {
                        let diff = (pp[aa_i] - cp[aa_i]).max(0.0);
                        let w = weights.get(&pid).unwrap_or(&0.0);
                        score += diff * w;
                    }
                }
            }
        }

        score += root_p[aa_i];

        let mut s_leaf = 0.0;
        for &lid in &tree.leaves {
            if gaps.contains(&lid) { continue; }
            if let Some(pid) = tree.nodes[&lid].parent_id {
                if let Some(pp) = node_probs.get(&pid) {
                    let diff = (leaf_probs[&lid][aa_i] - pp[aa_i]).max(0.0);
                    let w = weights.get(&lid).unwrap_or(&0.0);
                    s_leaf += diff * w;
                }
            }
        }
        score += s_leaf;
        scores[aa_i] = score / ((tree.internal.len() + tree.leaves.len()) as f64);
    }
    scores
}

pub fn tolerance(id: &str) -> Result<(), Box<dyn Error>> {
    let file_fasta = format!("{}_MaskedMSA.fasta", id);
    let file_nwk = format!("{}.treefile", id);
    let file_rst = format!("{}.state", id);
    let output_path = format!("ToleranceScores/{}.csv", id);

    let tree = read_tree_structure(&file_nwk)?;
    let msa = read_fasta(&file_fasta)?;
    let (state, total_pos) = read_state_file(&file_rst)?;
    
    let mut dists = HashMap::new();
    let mut stack = vec![(tree.root_id, 0.0)];
    
    while let Some((curr, d)) = stack.pop() {
        dists.insert(curr, d);
        for &child in &tree.nodes[&curr].children {
            let len = tree.nodes[&child].branch_length;
            stack.push((child, d + len));
        }
    }
    
    let sum_dist: f64 = dists.values().sum();
    let count_dist = dists.len() as f64;
    let mean_dist = sum_dist / count_dist;
    
    let mut weights = HashMap::new();
    for (node, d) in dists {
        let w = (- (d.powi(2)) / mean_dist.powi(2)).exp();
        weights.insert(node, w);
    }

    std::fs::create_dir_all("ToleranceScores")?;
    let file = File::create(output_path)?;
    let mut writer = std::io::BufWriter::new(file);

    write!(writer, "Pos/AA")?;
    for aa in AMINO_ACIDS { write!(writer, ",{}", aa)?; }
    writeln!(writer)?;

    for ps in 1..=total_pos {
        let scores = calculate_score(ps, &tree, &state, &msa, &weights);
        write!(writer, "{}", ps)?;
        for s in &scores { 
            if *s == 0.0 {
                write!(writer, ",0")?;
            } else {
                write!(writer, ",{:.18}", s)?; 
            }
        }
        writeln!(writer)?;
    }
    
    Ok(())
}