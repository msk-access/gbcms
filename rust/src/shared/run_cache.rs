//! A small once-per-run cache for the inputs that every counting call of a run reads again:
//! the GTF annotation index and the RNA editing database.
//!
//! One run counts several BAMs one after another, and the MNP rescue counts each BAM a second
//! time; each of those calls used to parse the same GTF and the same editing database from
//! disk. A parse is a pure function of the name it is given, the file's bytes and, for the GTF,
//! the set of chromosomes asked for, so a result kept under a key that holds all of them is the
//! same result. The name counts as given, because the loaders read it: the editing database
//! picks gzip from a `.gz` ending, and the GTF's warnings print the name. So the key holds the
//! name as given, the file it resolves to (links followed, a relative name made absolute), on
//! Unix its device and inode, its length and its modification time. Two names for one file are
//! two entries; the same name read from another working folder, a file replaced at the same
//! path, and a file changed between two calls are each read again; a file that cannot be
//! inspected is never cached (the loader then reports its own error as before).
//!
//! The cache holds a few entries and drops the oldest when it is full, so a long Python
//! session that reads many different files does not keep them all.

use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::SystemTime;

/// What a cached parse depends on: the name, the file and, for the GTF, the chromosomes asked
/// for.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MemoKey {
    /// The name as the loader is given it: the loader reads it (the gzip choice, the text of its
    /// warnings), so another name for the same file is another parse.
    name: String,
    /// The file's resolved path: the same relative name in another working folder is another
    /// file.
    path: PathBuf,
    /// The file's device and inode: a file replaced at the same path (a rename over it) is
    /// another file, even with the same length and modification time.
    #[cfg(unix)]
    file_id: (u64, u64),
    len: u64,
    modified: SystemTime,
    /// The sorted chromosome names the parse was filtered to (empty when it has no filter).
    scope: Vec<String>,
}

impl MemoKey {
    /// The key for the file named `name` filtered to `scope`, or None when the file cannot be
    /// inspected.
    pub fn for_file<'a>(name: &str, scope: impl IntoIterator<Item = &'a String>) -> Option<Self> {
        let path = std::fs::canonicalize(name).ok()?;
        let meta = std::fs::metadata(&path).ok()?;
        #[cfg(unix)]
        let file_id = {
            use std::os::unix::fs::MetadataExt;
            (meta.dev(), meta.ino())
        };
        let mut scope: Vec<String> = scope.into_iter().cloned().collect();
        scope.sort_unstable();
        Some(Self {
            name: name.to_string(),
            path,
            #[cfg(unix)]
            file_id,
            len: meta.len(),
            modified: meta.modified().ok()?,
            scope,
        })
    }
}

/// How many parses one cache keeps.
const MEMO_CAPACITY: usize = 4;

pub struct Memo<V> {
    slots: Mutex<Vec<(MemoKey, Arc<V>)>>,
}

impl<V> Memo<V> {
    pub const fn new() -> Self {
        Self { slots: Mutex::new(Vec::new()) }
    }

    /// The value kept under `key`, or the one `build` makes (kept, unless `key` is None or
    /// `build` fails). The second element is true when the value was already kept. The
    /// parse runs outside the lock.
    pub fn get_or_try_build<E>(
        &self,
        key: Option<MemoKey>,
        build: impl FnOnce() -> Result<V, E>,
    ) -> Result<(Arc<V>, bool), E> {
        let Some(key) = key else {
            return Ok((Arc::new(build()?), false));
        };
        {
            let slots = self.slots.lock().unwrap_or_else(|e| e.into_inner());
            if let Some((_, v)) = slots.iter().find(|(k, _)| *k == key) {
                return Ok((Arc::clone(v), true));
            }
        }
        let value = Arc::new(build()?);
        let mut slots = self.slots.lock().unwrap_or_else(|e| e.into_inner());
        if let Some((_, v)) = slots.iter().find(|(k, _)| *k == key) {
            return Ok((Arc::clone(v), true)); // another call finished first: share its value
        }
        if slots.len() >= MEMO_CAPACITY {
            slots.remove(0);
        }
        slots.push((key, Arc::clone(&value)));
        Ok((value, false))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temp_file(name: &str, body: &str) -> String {
        let p = std::env::temp_dir().join(format!("gbcms_run_cache_{}_{}", std::process::id(), name));
        std::fs::File::create(&p).unwrap().write_all(body.as_bytes()).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn builds_once_per_key_and_rebuilds_for_another_scope_or_a_changed_file() {
        let path = temp_file("a.txt", "one");
        let memo: Memo<u32> = Memo::new();
        let calls = std::cell::Cell::new(0u32);
        let build = || -> Result<u32, ()> {
            calls.set(calls.get() + 1);
            Ok(calls.get())
        };
        let scope_a = vec!["1".to_string(), "2".to_string()];
        let scope_a_shuffled = vec!["2".to_string(), "1".to_string()];
        let scope_b = vec!["1".to_string()];

        let (v1, hit1) = memo.get_or_try_build(MemoKey::for_file(&path, &scope_a), build).unwrap();
        let (v2, hit2) = memo.get_or_try_build(MemoKey::for_file(&path, &scope_a_shuffled), build).unwrap();
        assert!(!hit1 && hit2, "the second call, with the same chromosomes in another order, reuses the parse");
        assert!(Arc::ptr_eq(&v1, &v2));
        assert_eq!(calls.get(), 1);

        let (_, hit3) = memo.get_or_try_build(MemoKey::for_file(&path, &scope_b), build).unwrap();
        assert!(!hit3, "another chromosome set is another parse");
        assert_eq!(calls.get(), 2);

        // The file changes (its length): read again.
        std::fs::write(&path, "one, longer").unwrap();
        let (_, hit4) = memo.get_or_try_build(MemoKey::for_file(&path, &scope_a), build).unwrap();
        assert!(!hit4, "a changed file is read again");
        assert_eq!(calls.get(), 3);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn a_failed_build_is_not_kept_and_a_missing_file_is_never_cached() {
        let path = temp_file("b.txt", "x");
        let memo: Memo<u32> = Memo::new();
        let key = MemoKey::for_file(&path, std::iter::empty());
        assert!(memo.get_or_try_build(key.clone(), || Err::<u32, &str>("boom")).is_err());
        let (v, hit) = memo.get_or_try_build(key, || Ok::<u32, &str>(7)).unwrap();
        assert_eq!((*v, hit), (7, false));
        std::fs::remove_file(&path).ok();
        assert!(MemoKey::for_file(&path, std::iter::empty()).is_none());
    }

    #[test]
    fn a_name_is_keyed_as_given_and_by_the_file_it_resolves_to() {
        // Two folders, a file of the same name, length and modification time in each.
        let base = std::env::temp_dir().join(format!("gbcms_run_cache_{}_names", std::process::id()));
        let (a, b) = (base.join("a"), base.join("b"));
        std::fs::create_dir_all(&a).unwrap();
        std::fs::create_dir_all(&b).unwrap();
        let when = SystemTime::UNIX_EPOCH + std::time::Duration::from_secs(1_700_000_000);
        for (dir, body) in [(&a, "one"), (&b, "two")] {
            let f = std::fs::File::create(dir.join("x.gtf")).unwrap();
            (&f).write_all(body.as_bytes()).unwrap();
            f.set_modified(when).unwrap();
        }
        let key = |p: &std::path::Path| MemoKey::for_file(p.to_str().unwrap(), std::iter::empty()).unwrap();
        assert_ne!(key(&a.join("x.gtf")), key(&b.join("x.gtf")), "the same name in another folder is another file");
        // The loaders read the name (the gzip choice, the warnings), so two names for one file are
        // two entries: a reuse gives what a fresh load under that name gives.
        assert_ne!(key(&a.join("x.gtf")), key(&b.join("..").join("a").join("x.gtf")), "one file, two spellings");
        assert_eq!(key(&a.join("x.gtf")), key(&a.join("x.gtf")), "one file, one spelling");
        #[cfg(unix)]
        {
            let link = base.join("link.gtf");
            std::os::unix::fs::symlink(a.join("x.gtf"), &link).unwrap();
            assert_ne!(key(&link), key(&a.join("x.gtf")), "a link and its target are two names");
            // A file renamed over the path, with the same length and modification time, is another file.
            let before = key(&a.join("x.gtf"));
            let f = std::fs::File::create(base.join("new.gtf")).unwrap();
            (&f).write_all(b"one").unwrap();
            f.set_modified(when).unwrap();
            std::fs::rename(base.join("new.gtf"), a.join("x.gtf")).unwrap();
            assert_ne!(key(&a.join("x.gtf")), before, "a file replaced at the same path is read again");
        }
        std::fs::remove_dir_all(&base).ok();
    }

    #[test]
    fn keeps_a_few_entries_and_drops_the_oldest() {
        let path = temp_file("c.txt", "y");
        let memo: Memo<usize> = Memo::new();
        for i in 0..(MEMO_CAPACITY + 1) {
            let scope = vec![format!("c{i}")];
            memo.get_or_try_build(MemoKey::for_file(&path, &scope), || Ok::<usize, ()>(i)).unwrap();
        }
        let first = vec!["c0".to_string()];
        let (_, hit) = memo.get_or_try_build(MemoKey::for_file(&path, &first), || Ok::<usize, ()>(0)).unwrap();
        assert!(!hit, "the oldest entry was dropped");
        std::fs::remove_file(&path).ok();
    }
}
