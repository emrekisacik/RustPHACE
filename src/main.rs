use std::env;
use std::process;

mod tolerance;
mod msa1;
mod msa2;
mod msa2matrix;
mod msa1matrix;
mod matrices;
mod phace;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: cargo run <mode> <id>");
        process::exit(1);
    }

    let mode = &args[1];
    let id = &args[2];

    match mode.as_str() {
        "tolerance" => {
             if let Err(e) = tolerance::tolerance(id) {
                eprintln!("Error: {}", e);
                process::exit(1);
            }
        },
        "msa1" => {
             if let Err(e) = msa1::msa1(id) {
                eprintln!("Error: {}", e);
                process::exit(1);
            }
        },
        "msa2" => {
             if let Err(e) = msa2::msa2(id) {
                eprintln!("Error: {}", e);
                process::exit(1);
            }
        },
        "msa1matrix" => {
             if let Err(e) = msa1matrix::msa1matrix(id) {
                eprintln!("Error: {}", e);
                process::exit(1);
            }
        },    
        "msa2matrix" => {
             if let Err(e) = msa2matrix::msa2matrix(id) {
                eprintln!("Error: {}", e);
                process::exit(1);
            }
        },
        "matrices" => {
            if let Err(e) = matrices::matrices(id) {
                eprintln!("Error: {}", e);
                process::exit(1);
            }
        },
        "phace" => {
            if let Err(e) = phace::phace(id) {
                eprintln!("Error: {}", e);
                process::exit(1);
            }
        },
        _ => {
            eprintln!("Invalid mode.");
            process::exit(1);
        }
    }
}