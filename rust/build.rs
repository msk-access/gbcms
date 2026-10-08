//! Records the commit a build came from, as `GBCMS_BUILD_COMMIT` (read by
//! `build_commit()`), so two builds of one version are told apart in output
//! provenance (`#gbcms v6.6.0.dev0 (9c371263)`).
//!
//! Taken from the `GBCMS_BUILD_COMMIT` environment variable when it holds a
//! commit (Docker and the release wheels build from a copy without `.git`), else
//! from git, when this source is tracked by the repository git finds (an sdist
//! unpacked inside another repository is not). Empty otherwise. The commit is the
//! checkout's HEAD when the extension was built: a dirty tree, or Python edits
//! after the last Rust build (an editable install), are not reflected.
use std::path::Path;
use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

/// The usual 8-character short form of a hex object name, or None.
fn short_commit(raw: &str) -> Option<String> {
    let c = raw.trim().to_ascii_lowercase();
    (c.len() >= 7 && c.bytes().all(|b| b.is_ascii_hexdigit())).then(|| c[..c.len().min(8)].to_string())
}

fn main() {
    println!("cargo:rerun-if-env-changed=GBCMS_BUILD_COMMIT");
    let from_env = std::env::var("GBCMS_BUILD_COMMIT").ok().and_then(|c| short_commit(&c));
    // git only when it tracks this very file (not a repository the source sits in).
    let tracked = git(&["ls-files", "--error-unmatch", "build.rs"]).is_some();
    let commit = from_env
        .or_else(|| tracked.then(|| git(&["rev-parse", "HEAD"])).flatten().and_then(|c| short_commit(&c)))
        .unwrap_or_default();
    // Rebuild when the checkout moves to another commit. Watch only paths that
    // exist (cargo treats a missing one as always changed): HEAD in the git dir;
    // branch refs in the common dir (a worktree's git dir holds none).
    if tracked {
        let mut watch = Vec::new();
        if let Some(dir) = git(&["rev-parse", "--absolute-git-dir"]) {
            watch.push(format!("{dir}/HEAD"));
        }
        if let Some(common) = git(&["rev-parse", "--path-format=absolute", "--git-common-dir"]) {
            watch.push(format!("{common}/refs/heads"));
            watch.push(format!("{common}/packed-refs"));
        }
        for path in watch.iter().filter(|p| Path::new(p).exists()) {
            println!("cargo:rerun-if-changed={path}");
        }
    }
    println!("cargo:rustc-env=GBCMS_BUILD_COMMIT={commit}");
}
