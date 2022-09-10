#![feature(let_else)]

use bam;
use clap::Parser;
use std::{
    ffi::OsStr,
    io::{BufWriter, Write},
    path::PathBuf,
};

#[derive(Parser)]
#[clap(author, about)]
struct Cli {
    /// Input .bam file.
    #[clap(parse(from_os_str))]
    input: PathBuf,

    /// Output .seq file. Defaults to <input>.seq.
    #[clap(parse(from_os_str))]
    output: Option<PathBuf>,

    /// Disable trimming of soft clipped regions from the read.
    #[clap(long)]
    no_clip: bool,

    /// Only output reads at least this long.
    #[clap(long)]
    min_len: Option<usize>,
}

fn main() {
    let args = Cli::parse();

    assert_eq!(
        args.input.extension(),
        Some(OsStr::new("bam")),
        "Input file must have .bam extension!"
    );
    if let Some(output) = &args.output {
        assert_eq!(
            output.extension(),
            Some(OsStr::new("seq")),
            "Output file must have .seq extension!"
        );
    }

    let mut output = BufWriter::new(
        std::fs::File::options()
            .write(true)
            .create(true)
            .truncate(true)
            .open(args.output.unwrap_or(args.input.with_extension("seq")))
            .unwrap(),
    );

    let reader = bam::BamReader::from_path(args.input, 4).unwrap();
    let mut taken = 0;
    let mut skipped = 0;
    let mut total_len = 0;
    for record in reader {
        let record = match record {
            Ok(record) => record,
            Err(err) => {
                println!("Broken record: {err}");
                break;
            }
        };

        if !record.sequence().available() || record.cigar().is_empty() {
            skipped += 1;
            continue;
        }

        if let Some(l) = args.min_len {
            let mut len = record.query_len();
            if !args.no_clip {
                len -= record.cigar().soft_clipping(true);
                len -= record.cigar().soft_clipping(false);
            }
            if (len as usize) < l {
                skipped += 1;
                continue;
            }
        }

        let mut reference = Vec::new();
        let Ok(alignment_entries) = record.alignment_entries() else {
            skipped += 1;
            continue;
        };
        for entry in alignment_entries {
            if let Some((_, ref_nt)) = entry.ref_pos_nt() {
                reference.push(ref_nt);
            }
        }

        let read = record.sequence().to_vec();
        let mut read = &read[..];
        if !args.no_clip {
            let left_clip = record.cigar().soft_clipping(true) as usize;
            let right_clip = record.cigar().soft_clipping(false) as usize;
            read = &read[left_clip..read.len() - right_clip];
        }

        taken += 1;
        total_len += read.len();
        output.write_all(">".as_bytes()).unwrap();
        output.write_all(read).unwrap();
        output.write_all("\n".as_bytes()).unwrap();
        output.write_all("<".as_bytes()).unwrap();
        output.write_all(&reference).unwrap();
        output.write_all("\n".as_bytes()).unwrap();
    }
    println!("Records taken:   {taken:>8}");
    println!("Records skipped: {skipped:>8}");
    println!("Average length:  {:>8}", total_len / taken);
}
