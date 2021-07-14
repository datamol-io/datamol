**Added:**

* More tests for the `dm.similarity` modules + check against RDKit equivalent methods.
* `dm.same_mol(mol1, mol2)` to check whether 2 molecules are the same based on their InChiKey.

**Changed:**

* use `scipy` in `dm.similarity.pdist()`.
* Raise an error when a molecule is invalid in `dm.similarity.pdist/cdist`.

**Deprecated:**

* `dm.similarity.pdist()` nows returns only the dist matrix without the `valid_idx` vector.

**Removed:**

* <news item>

**Fixed:**

* A bug returning an inconsistent dist matrix with `dm.similarity.pdist()`.

**Security:**

* <news item>
