//! Records the commit a build came from, as `GBCMS_BUILD_COMMIT` (read by
//! `build_commit()`), so two builds of one version are told apart in output
//! provenance (`#gbcms v6.6.0.dev0 (9c371263)`).
//!
//! Taken from the `GBCMS_BUILD_COMMIT` environment variable when set (Docker and
//! the release wheels build from a copy without `.git`), else from git. Empty
//! when neither is available (an sdist built outside a checkout).
use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn main() {
    println!("cargo:rerun-if-env-changed=GBCMS_BUILD_COMMIT");
    let commit = std::env::var("GBCMS_BUILD_COMMIT")
        .ok()
        .filter(|c| !c.trim().is_empty())
        .or_else(|| git(&["rev-parse", "HEAD"]))
        .unwrap_or_default()
        .trim()
        .to_ascii_lowercase();
    // Only a hex object name is a commit; keep the usual 8-character short form.
    let commit = if commit.len() >= 7 && commit.bytes().all(|b| b.is_ascii_hexdigit()) {
        commit[..commit.len().min(8)].to_string()
    } else {
        String::new()
    };
    // Rebuild when the checkout moves to another commit.
    if let Some(dir) = git(&["rev-parse", "--absolute-git-dir"]) {
        println!("cargo:rerun-if-changed={dir}/HEAD");
        println!("cargo:rerun-if-changed={dir}/refs/heads");
        println!("cargo:rerun-if-changed={dir}/packed-refs");
    }
    println!("cargo:rustc-env=GBCMS_BUILD_COMMIT={commit}");
}
