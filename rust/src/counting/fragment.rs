//! Re-export of fragment utilities from `shared::fragment`.
//!
//! The canonical implementation lives in `crate::shared::fragment`.
//! This re-export maintains backward compatibility for existing
//! `super::fragment::*` imports within the counting module.

pub use crate::shared::fragment::*;
