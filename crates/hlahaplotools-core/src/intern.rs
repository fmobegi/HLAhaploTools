use std::collections::HashMap;

/// Integer handle for an allele name such as `A*02:01` or `DRB1*04:01:01G`.
///
/// IDs are handed out in order of first appearance, so the same input file
/// always produces the same IDs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AlleleId(pub u32);

/// Two-way mapping between allele names and [`AlleleId`]s.
///
/// Each allele also remembers which locus it belongs to, as an index into
/// [`crate::Dataset::loci`].
#[derive(Debug, Default, Clone)]
pub struct AlleleTable {
    names: Vec<String>,
    loci: Vec<u16>,
    index: HashMap<String, AlleleId>,
}

impl AlleleTable {
    pub fn new() -> Self {
        Self::default()
    }

    /// Return the ID for `name`, adding it if it is new.
    pub fn intern(&mut self, locus: u16, name: &str) -> AlleleId {
        if let Some(&id) = self.index.get(name) {
            return id;
        }
        let id = AlleleId(u32::try_from(self.names.len()).expect("more than u32::MAX alleles"));
        self.names.push(name.to_owned());
        self.loci.push(locus);
        self.index.insert(name.to_owned(), id);
        id
    }

    pub fn get(&self, name: &str) -> Option<AlleleId> {
        self.index.get(name).copied()
    }

    pub fn name(&self, id: AlleleId) -> &str {
        &self.names[id.0 as usize]
    }

    pub fn locus(&self, id: AlleleId) -> u16 {
        self.loci[id.0 as usize]
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interning_is_stable_and_deduplicated() {
        let mut t = AlleleTable::new();
        let a = t.intern(0, "A*01:01");
        let b = t.intern(0, "A*02:01");
        assert_eq!(t.intern(0, "A*01:01"), a);
        assert_eq!((a, b), (AlleleId(0), AlleleId(1)));
        assert_eq!(t.name(b), "A*02:01");
        assert_eq!(t.get("A*02:01"), Some(b));
        assert_eq!(t.get("A*03:01"), None);
        assert_eq!(t.len(), 2);
    }
}
