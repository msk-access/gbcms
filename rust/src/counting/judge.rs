//! The haplotype judge (`--read-judge haplotype`, #304): one read judgment for every
//! small-variant row. A read is scored whole against the row's local haplotypes (every
//! context of the other listed rows within read reach, crossed with every state of the
//! row's site) and counts ALT or REF only by one confident base; otherwise it is depth
//! only. Design: `docs/proposals/haplotype-judge.md`.

/// How far a read reaches each side of a row's REF span: rows within it can share a read
/// with the row (2x250 and 2x300 runs included), so they shape its haplotypes.
pub(crate) const READ_REACH: i64 = super::carrier::WHOLE_READ_REACH;

/// Which read judgment counts a row's reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ReadJudge {
    /// The size-specific rules (RJ-1 to RJ-27): the default.
    Rules,
    /// The haplotype judge.
    Haplotype,
}

/// Parse the `read_judge` token strictly, as `alignment_backend` is: a typo must fail
/// loudly rather than count with the other judgment.
pub(crate) fn parse_read_judge(token: &str) -> Result<ReadJudge, String> {
    match token {
        "rules" => Ok(ReadJudge::Rules),
        "haplotype" => Ok(ReadJudge::Haplotype),
        other => Err(format!("unknown read_judge {other:?}; expected \"rules\" or \"haplotype\"")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_read_judge_token_is_parsed_strictly() {
        assert_eq!(parse_read_judge("rules").unwrap(), ReadJudge::Rules);
        assert_eq!(parse_read_judge("haplotype").unwrap(), ReadJudge::Haplotype);
        for bad in ["Haplotype", "haplotypes", "", "rule"] {
            assert!(parse_read_judge(bad).is_err(), "{bad:?} must be rejected");
        }
    }
}
