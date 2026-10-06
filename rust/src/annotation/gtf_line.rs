//! One GTF line to the exon fields the annotation index reads.
//!
//! The grammar is noodles-gtf 0.37's (the parser this replaced), field by field, so
//! the same lines load: nine tab-separated columns (the ninth may hold tabs) after
//! trailing whitespace is trimmed; start and end positive integers; score `.` or a
//! float; strand `.`, `+` or `-`; frame `.` or 0–2; attributes as `key value`
//! entries, each value quoted (to the next `"`) or raw (to the next `;`), entries
//! split by an optional `;`, the first of a repeated key winning. An unquoted value
//! must be followed by a `;` unless it is empty at the end of the line (noodles
//! rejected such a line, and GTF2.2 ends every attribute with one). Four departures,
//! each a line noodles loaded wrongly:
//! - a whitespace run may separate a key from its value (noodles kept the run in a
//!   raw value, so `transcript_id  "T1"` gave the ID ` "T1"`, quotes included);
//! - raw values are trimmed;
//! - start > end is rejected (noodles kept a negative-length exon);
//! - a coordinate past `i32` is rejected (it wrapped negative in the index).
//!
//! Callers read the feature column first and check the rest only for the lines
//! they keep, so most of a GTF costs a scan to its third tab.

/// A line's columns after the feature (the caller has read the first three);
/// the attributes column runs to the end of the line.
pub(crate) struct Columns<'a> {
    start: &'a str,
    end: &'a str,
    score: &'a str,
    strand: &'a str,
    frame: &'a str,
    attributes: &'a str,
}

/// An exon line's fields: 0-based half-open coordinates, strand `+`, `-` or `.`.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct ExonLine<'a> {
    pub(crate) start: i32,
    pub(crate) end: i32,
    pub(crate) strand: char,
    pub(crate) transcript_id: Option<&'a str>,
}

/// Why a line was rejected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LineError {
    MissingColumns,
    InvalidStart,
    InvalidEnd,
    StartAfterEnd,
    CoordinateOverflow,
    InvalidScore,
    InvalidStrand,
    InvalidFrame,
    InvalidAttributes,
    /// Not a grammar error: the caller's rule for an exon without a usable ID.
    MissingTranscriptId,
    /// The line is not UTF-8 text.
    NotUtf8,
}

impl std::fmt::Display for LineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::MissingColumns => "fewer than nine tab-separated columns",
            Self::InvalidStart => "start is not a positive integer",
            Self::InvalidEnd => "end is not a positive integer",
            Self::StartAfterEnd => "start is after end",
            Self::CoordinateOverflow => "a coordinate exceeds 2,147,483,647",
            Self::InvalidScore => "score is neither '.' nor a number",
            Self::InvalidStrand => "strand is not '+', '-' or '.'",
            Self::InvalidFrame => "frame is not '.', 0, 1 or 2",
            Self::InvalidAttributes => "attributes are not `key value;` entries",
            Self::MissingTranscriptId => "transcript_id is missing or empty",
            Self::NotUtf8 => "the line is not UTF-8 text",
        })
    }
}

/// Split a line (no newline) into its columns.
pub(crate) fn columns(line: &str) -> Result<Columns<'_>, LineError> {
    let mut it = line.trim_end().splitn(9, '\t');
    let mut next = || it.next().ok_or(LineError::MissingColumns);
    for _ in 0..3 {
        next()?; // seqname, source, feature
    }
    Ok(Columns {
        start: next()?,
        end: next()?,
        score: next()?,
        strand: next()?,
        frame: next()?,
        attributes: next()?,
    })
}

impl<'a> Columns<'a> {
    /// The exon fields, after checking every column as noodles did.
    pub(crate) fn exon(&self) -> Result<ExonLine<'a>, LineError> {
        let start = position(self.start).ok_or(LineError::InvalidStart)?;
        let end = position(self.end).ok_or(LineError::InvalidEnd)?;
        if start > end {
            return Err(LineError::StartAfterEnd);
        }
        if end > i32::MAX as usize {
            return Err(LineError::CoordinateOverflow);
        }
        if self.score != "." && self.score.parse::<f32>().is_err() {
            return Err(LineError::InvalidScore);
        }
        let strand = match self.strand {
            "+" => '+',
            "-" => '-',
            "." => '.',
            _ => return Err(LineError::InvalidStrand),
        };
        if self.frame != "." && !matches!(self.frame.parse::<u8>(), Ok(0..=2)) {
            return Err(LineError::InvalidFrame);
        }
        let transcript_id = transcript_id(self.attributes)?;
        Ok(ExonLine { start: start as i32 - 1, end: end as i32, strand, transcript_id })
    }
}

/// A 1-based position: a positive integer (noodles' `Position`).
fn position(s: &str) -> Option<usize> {
    s.parse::<usize>().ok().filter(|&p| p > 0)
}

/// The first `transcript_id` in the attributes column, after checking every
/// entry parses. Scans bytes: every delimiter is ASCII, so each slice taken at
/// one is on a character boundary.
fn transcript_id(column: &str) -> Result<Option<&str>, LineError> {
    let b = column.as_bytes();
    let n = b.len();
    let ws = |c: u8| c == b' ' || c == b'\t';
    let find = |from: usize, want: u8| b[from..].iter().position(|&c| c == want).map(|p| from + p);
    if n == 0 {
        return Err(LineError::InvalidAttributes);
    }
    let mut tx = None;
    let mut i = 0;
    while i < n {
        let key_start = i;
        while i < n && !ws(b[i]) {
            i += 1;
        }
        if i == n {
            return Err(LineError::InvalidAttributes); // a key with no value
        }
        let key = &column[key_start..i];
        while i < n && ws(b[i]) {
            i += 1;
        }
        let value = if i < n && b[i] == b'"' {
            let close = find(i + 1, b'"').ok_or(LineError::InvalidAttributes)?;
            let v = &column[i + 1..close];
            i = close + 1;
            v
        } else {
            match find(i, b';') {
                Some(semi) => {
                    let v = column[i..semi].trim_end();
                    i = semi;
                    v
                }
                None if i == n => "",
                None => return Err(LineError::InvalidAttributes),
            }
        };
        let rest = column[i..].trim_start();
        i = n - rest.strip_prefix(';').unwrap_or(rest).trim_start().len();
        if key == "transcript_id" && tx.is_none() {
            tx = Some(value);
        }
    }
    Ok(tx)
}

#[cfg(test)]
mod tests {
    use super::*;

    const A: &str = r#"gene_id "G1"; transcript_id "T1";"#;

    fn line(start: &str, end: &str, score: &str, strand: &str, frame: &str, attrs: &str) -> String {
        format!("1\tsrc\texon\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attrs}")
    }

    fn exon(l: &str) -> Result<ExonLine<'_>, LineError> {
        columns(l)?.exon()
    }

    fn ok(tx: &'static str) -> ExonLine<'static> {
        ExonLine { start: 99, end: 200, strand: '+', transcript_id: Some(tx) }
    }

    #[test]
    fn noodles_lines_load_the_same_fields() {
        // Each of these loaded under noodles with these fields (diffed in the
        // GTF-loading measurement); the parser must agree.
        for (name, attrs) in [
            ("plain", A.to_string()),
            ("no final semicolon", r#"gene_id "G1"; transcript_id "T1""#.to_string()),
            ("unquoted value with its semicolon", r#"gene_id "G1"; transcript_id "T1"; exon_number 1;"#.to_string()),
            ("transcript first", r#"transcript_id "T1"; gene_id "G1";"#.to_string()),
            ("no space after ;", r#"gene_id "G1";transcript_id "T1";"#.to_string()),
            ("look-alike key", r#"orig_transcript_id "X9"; gene_id "G1"; transcript_id "T1";"#.to_string()),
            ("unquoted values", "gene_id G1; transcript_id T1;".to_string()),
            ("quoted semicolon", r#"note "a; b"; gene_id "G1"; transcript_id "T1";"#.to_string()),
            ("unquoted number", r#"gene_id "G1"; transcript_id "T1"; exon_number 1;"#.to_string()),
            ("repeated key", r#"gene_id "G1"; transcript_id "T1"; tag "basic"; tag "CCDS";"#.to_string()),
            ("first of a repeat wins", r#"gene_id "G1"; transcript_id "T1"; transcript_id "T2";"#.to_string()),
            ("utf-8 value", r#"gene_id "G1"; transcript_id "T1"; note "é";"#.to_string()),
        ] {
            assert_eq!(exon(&line("100", "200", ".", "+", ".", &attrs)), Ok(ok("T1")), "{name}");
        }
        for (name, l) in [
            ("crlf", format!("{}\r", line("100", "200", ".", "+", ".", A))),
            ("trailing tab", format!("{}\t", line("100", "200", ".", "+", ".", A))),
            ("trailing spaces", format!("{}  ", line("100", "200", ".", "+", ".", A))),
            ("score and frame", line("100", "200", "0.5", "+", "2", A)),
            ("float forms", line("100", "200", "1e5", "+", "0", A)),
        ] {
            assert_eq!(exon(&l), Ok(ok("T1")), "{name}");
        }
        let unstranded = line("100", "200", ".", ".", ".", A);
        assert_eq!(exon(&unstranded).unwrap().strand, '.');
        let minus = line("1", "1", ".", "-", ".", A);
        let minus = exon(&minus).unwrap();
        assert_eq!((minus.start, minus.end, minus.strand), (0, 1, '-'));
    }

    #[test]
    fn noodles_rejections_are_kept() {
        for (name, l, want) in [
            ("start 0", line("0", "200", ".", "+", ".", A), LineError::InvalidStart),
            ("negative start", line("-5", "200", ".", "+", ".", A), LineError::InvalidStart),
            ("non-numeric start", line("1x0", "200", ".", "+", ".", A), LineError::InvalidStart),
            ("non-numeric end", line("100", "2o0", ".", "+", ".", A), LineError::InvalidEnd),
            ("bad score", line("100", "200", "abc", "+", ".", A), LineError::InvalidScore),
            ("strand ?", line("100", "200", ".", "?", ".", A), LineError::InvalidStrand),
            ("strand x", line("100", "200", ".", "x", ".", A), LineError::InvalidStrand),
            ("empty strand", line("100", "200", ".", "", ".", A), LineError::InvalidStrand),
            ("frame 7", line("100", "200", ".", "+", "7", A), LineError::InvalidFrame),
            ("empty frame", line("100", "200", ".", "+", "", A), LineError::InvalidFrame),
            ("escaped quote", line("100", "200", ".", "+", ".", r#"gene_id "G1"; transcript_id "T\"1";"#), LineError::InvalidAttributes),
            ("unclosed quote", line("100", "200", ".", "+", ".", r#"gene_id "G1; transcript_id "T1";"#), LineError::InvalidAttributes),
            ("key without value", line("100", "200", ".", "+", ".", r#"gene_id "G1"; lonely"#), LineError::InvalidAttributes),
            ("unquoted last value, no ;", line("100", "200", ".", "+", ".", r#"gene_id "G1"; transcript_id "T1"; exon_number 1"#), LineError::InvalidAttributes),
            ("no semicolons", line("100", "200", ".", "+", ".", "transcript_id T1 gene_id G1"), LineError::InvalidAttributes),
            ("trailing comment", line("100", "200", ".", "+", ".", r#"gene_id "G1"; transcript_id "T1"; # note"#), LineError::InvalidAttributes),
            ("extra columns", line("100", "200", ".", "+", ".", "gene_id \"G1\"; transcript_id \"T1\";\textra\tcols"), LineError::InvalidAttributes),
            ("empty attributes", "1\tsrc\texon\t100\t200\t.\t+\t.\t".to_string(), LineError::MissingColumns),
        ] {
            assert_eq!(exon(&l).map(|_| ()), Err(want), "{name}");
        }
        assert_eq!(columns("1\tsrc\texon\t100\t200\t.\t+\t.").err(), Some(LineError::MissingColumns));
        assert_eq!(columns("1 src exon 100 200 . + . x").err(), Some(LineError::MissingColumns));
    }

    #[test]
    fn whitespace_runs_separate_key_and_value() {
        for attrs in [
            r#"gene_id  "G1";  transcript_id  "T1";"#,
            "gene_id\t\"G1\"; transcript_id \t \"T1\";",
            "gene_id G1 ; transcript_id T1 ;",
        ] {
            assert_eq!(exon(&line("100", "200", ".", "+", ".", attrs)), Ok(ok("T1")), "{attrs}");
        }
    }

    #[test]
    fn lines_noodles_loaded_wrongly_are_rejected() {
        assert_eq!(exon(&line("300", "200", ".", "+", ".", A)).map(|_| ()), Err(LineError::StartAfterEnd));
        assert_eq!(
            exon(&line("100", "3000000000", ".", "+", ".", A)).map(|_| ()),
            Err(LineError::CoordinateOverflow)
        );
        let max = line("100", "2147483647", ".", "+", ".", A);
        assert_eq!(exon(&max).unwrap().end, i32::MAX);
    }

    #[test]
    fn an_empty_unquoted_value_may_end_the_column() {
        // noodles accepted it (the line trim makes it unreachable through `columns`).
        assert_eq!(transcript_id(r#"transcript_id "T1"; tag "#), Ok(Some("T1")));
        assert_eq!(transcript_id(r#"transcript_id "T1"; tag x"#), Err(LineError::InvalidAttributes));
    }

    #[test]
    fn missing_ids_are_none() {
        let l = line("100", "200", ".", "+", ".", r#"gene_id "G1"; orig_transcript_id "X9";"#);
        assert_eq!(exon(&l).unwrap().transcript_id, None);
        let l = line("100", "200", ".", "+", ".", r#"transcript_id "";"#);
        assert_eq!(exon(&l).unwrap().transcript_id, Some(""));
    }
}
