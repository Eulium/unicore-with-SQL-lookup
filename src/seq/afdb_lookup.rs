use crate::seq::fasta_io::write_fasta;
use crate::seq::create_gene_specific_fasta::read_db;
use crate::util::message::print_message as mprint;
use crate::util::message::println_message as mprintln;
use crate::util::command as cmd;
use crate::envs::error_handler as err;

use std::collections::HashMap;
use std::io::prelude::*;
use std::path::MAIN_SEPARATOR as SEP;
use std::process::Command as Cmd;
use duckdb::{params, Connection, Result};
use std::time::Instant;

use reqwest;

#[derive(Debug)]
struct AFDB_entry {
    name: String,
    seq: String,
    ss: String
}

fn download_table(path: &String) -> Result<(), Box<dyn std::error::Error>> {
    // create directory if not exists
    let mut md5_path = format!("{}{}md5", path, SEP);
    if std::path::Path::new(&md5_path).exists() {
        mprintln(&format!("Directory {} already exists.", md5_path), 4);
    } else {
        mprintln(&format!("Creating directory {}...", md5_path), 4);
        std::fs::create_dir_all(&md5_path)?;
    }

    // download the tables
    mprint(&"Downloading AFDB lookup tables (this may take a while)... 0.0%".to_string(), 3);
    for i in 0..256 {
        let hex = format!("{:02x}", i);
        let url = format!("https://unicore.steineggerlab.workers.dev/md5/{}.tsv.gz", hex);
        let file = format!("{}{}{}.tsv.gz", md5_path, SEP, hex);
        let mut resp = reqwest::blocking::get(&url)?;
        let mut file = std::fs::File::create(&file)?;
        std::io::copy(&mut resp, &mut file)?;
        mprint(&format!("\rDownloading AFDB lookup tables (this may take a while)... {:.1}%", (i as f64 + 1.0) / 2.56), 3);
    }
    mprintln(&"\rDownloading AFDB lookup tables (this may take a while)... 100.0% Done".to_string(), 3);

    //no need to decompress the tables, when using DuckDB
    mprint(&"Decompressing the tables... 0.0%".to_string(), 3);
    for i in 0..256 {
        let hex = format!("{:02x}", i);
        let file = format!("{}{}{}.tsv.gz", &md5_path, SEP, hex);
        cmd::run(Cmd::new("gzip").arg("-d").arg(&file));
        mprint(&format!("\rDecompressing the tables... {:.1}%", (i as f64 + 1.0) / 2.56), 3);
    }
    mprintln(&"\rDecompressing the tables... 100.0% Done".to_string(), 3);

    // create AFDB hive partitoned parquet files
    mprintln(&"\rCreating hive partitioned parquet file for faster lookup".to_string(), 3);
    let conn = Connection::open_in_memory()?;
    let files = format!("{}{}md5{}{}.tsv", path,SEP,SEP,'*');
    let AFDB_path = format!("{}{}AFDB.pq", path, SEP);
    let mut stmt = format!(
        r#"
        SET preserve_insertion_order TO false; 
        COPY
            (
            SELECT 
                hash,
                seq,
                hash[:2] as hex
            FROM
                read_csv('{files}', header = false, delim = '\t', names = ['hash', 'seq'])
            )
        TO
        '{AFDB_path}'
        (FORMAT 'parquet', PARTITION_BY hex, compression 'zstd');
        ;"#);
    conn.execute_batch(&stmt)?;

    Ok(())
}

pub fn run_afdb(fasta_data: &HashMap<String, String>, afdb_lookup: &String, converted_aa: &String, converted_ss: &String, combined_aa: &String) -> Result<(), Box<dyn std::error::Error>> {
    
    // Open DuckDB connection and set path to AFDB
    let conn = Connection::open_in_memory()?;
    let path = afdb_lookup.clone();
    let mut afdb_path = format!("{}{}AFDB.pq", path, SEP);
    let mut hive_path = format!("{}{}AFDB.pq{}hex=*{}data_*.parquet", path, SEP,SEP, SEP);

    // check if the directory is present
    if std::path::Path::new(&path).exists() && std::fs::File::open(&format!("{}{}hex=00{}data_0.parquet", path, SEP, SEP)).is_ok() { 
        afdb_path = path.clone();
        hive_path = format!("{}{}hex=*{}data_*.parquet", path, SEP, SEP);
    }
    if std::fs::File::open(&format!("{}{}hex=00{}data_0.parquet", afdb_path, SEP, SEP)).is_err() {
        mprintln(&"AFDB lookup tables not found.".to_string(), 0);
        mprint(&format!("Trying to download the tables to {} (~30GB). Continue? [y/n]: ", path), 0);
        let mut input = String::new();
        std::io::stdin().read_line(&mut input)?;
        if input.trim().to_lowercase() != "y" {
            mprintln(&"Download cancelled. Aborting the program.".to_string(), 0);
            std::process::exit(0);
        }
        // use path/md5/ to download and store tsv tables 
        // download_table will create parquet files in path/AFDB.pq
        download_table(&path)?;
    }

    // Load query sequences into DuckDB
    conn.execute_batch(r"
        CREATE TEMP TABLE query_seq 
            (
            name VARCHAR, 
            hash VARCHAR,
            seq VARCHAR, 
            hex VARCHAR
            );",
        )?;
    let mut app = conn.appender("query_seq")?;
    
    for (h, seq) in fasta_data {
        // add line feed to the end of the sequence
        let mut bytes = seq.clone().into_bytes(); bytes.push(10);
        let hash = format!("{:x}", md5::compute(bytes));
        let idx = usize::from_str_radix(&hash[..2], 16)?;
        let hex = format!("{:02x}", idx);
        app.append_row([h, &hash, seq, &hex])?;
    }
    app.flush()?;

    // Prep and run JOIN on first two hex digits and full hash
    // Join by hex allows for filter pushdown into parquet files
    // If no 3Di match is found, the AFDB entry will be 'missing'
    let query = format!(r#"
    SELECT 
        query_seq.name AS name,
        query_seq.seq AS seq,
        IF (AFDB.seq IS NULL , 'missing', AFDB.seq)  AS ss
    FROM
        query_seq
    LEFT JOIN
        read_parquet('{hive_path}') AS AFDB
    ON 
        AFDB.hex = query_seq.hex
    AND
        query_seq.hash = AFDB.hash
    ;"#
    ); 

    mprintln(&"\rLooking up the AFDB database...".to_string(), 3);
    let mut stmt = conn.prepare(&query)?;
    let (mut conv, mut pred) = (0, 0);
    let mut AFDB_entry_itr = stmt.query_map(params![], |row| 
            Ok(AFDB_entry {
                name: row.get(0)?,
                seq: row.get(1)?,
                ss: row.get(2)?,
            })
        )?;
    

    let mut converted_aa_data: HashMap<String, String> = HashMap::new();
    let mut converted_ss_data: HashMap<String, String> = HashMap::new();
    let mut combined_data: HashMap<String, String> = HashMap::new();

    // Load into appropriate hashmaps
    for afdb_entry in AFDB_entry_itr.flatten() {
        if afdb_entry.ss != "missing" {
            converted_aa_data.insert(afdb_entry.name.clone(), afdb_entry.seq.clone());
            converted_ss_data.insert(afdb_entry.name.clone(), afdb_entry.ss.clone());
            conv += 1;
        }
        else {
            combined_data.insert(afdb_entry.name.clone(), afdb_entry.seq.clone());
            pred += 1;
        }
    }

    mprintln(&"\rLooking up the AFDB database... Done".to_string(), 3);
    mprintln(&format!("{} sequences found from the lookup database", conv), 3);
    mprintln(&format!("{} sequences not found and will be predicted", pred), 3);

    write_fasta(&converted_aa, &converted_aa_data, true)?;
    write_fasta(&converted_ss, &converted_ss_data, true)?;
    write_fasta(&combined_aa, &combined_data, false)?;

    Ok(())
}

pub fn run_custom(fasta_data: &HashMap<String, String>, custom_lookup: &String, converted_aa: &String, converted_ss: &String, combined_aa: &String) -> Result<(), Box<dyn std::error::Error>> {
    // check if the directory is present
    let path = custom_lookup.clone();
    let ss_path = format!("{}_ss", path);
    if std::fs::File::open(&path).is_err() || std::fs::File::open(&ss_path).is_err() {
        err::error(err::ERR_GENERAL, Some("Custom lookup database does not exist or improperly formatted.".to_string()));
    }

    let mut converted_aa_data: HashMap<String, String> = HashMap::new();
    let mut converted_ss_data: HashMap<String, String> = HashMap::new();
    let mut combined_data: HashMap<String, String> = HashMap::new();

    // load table to memory
    mprintln(&"\nLoading the database...".to_string(), 3);
    let table_aa = read_db(&path);
    let table_ss = read_db(&ss_path);
    if table_aa.len() != table_ss.len() {
        err::error(err::ERR_GENERAL, Some("The custom lookup database is not properly formatted.".to_string()));
    }
    let mut table_map = HashMap::<String, String>::new();
    for (aa, ss) in table_aa.into_iter().zip(table_ss.into_iter()) {
        table_map.insert(aa, ss);
    }

    let cnt = fasta_data.len();
    let (mut conv, mut pred) = (0, 0);
    mprint(&"Looking up the custom database... 0.0%".to_string(), 3);
    for (i, (h, seq)) in fasta_data.iter().enumerate() {
        mprint(&format!("\rLooking up the custom database... {:.1}%", (i as f64 + 1.0) * 100.0 / (cnt as f64)), 3);
        match table_map.get(seq) {
            Some(converted_seq) => {
                converted_aa_data.insert(h.clone(), seq.clone());
                converted_ss_data.insert(h.clone(), converted_seq.clone());
                conv += 1;
            },
            None => {
                combined_data.insert(h.clone(), seq.clone());
                pred += 1;
            },
        }
    }
    mprintln(&"\rLooking up the custom database... 100.0% Done".to_string(), 3);
    mprintln(&format!("{} sequences found from the lookup database", conv), 3);
    mprintln(&format!("{} sequences not found and will be predicted", pred), 3);

    write_fasta(&converted_aa, &converted_aa_data, true)?;
    write_fasta(&converted_ss, &converted_ss_data, true)?;
    write_fasta(&combined_aa, &combined_data, false)?;

    Ok(())
}