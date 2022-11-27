**Added:**

* Support for `max_num_mols` when reading SDF files. Useful when files are large and debugging code.
* Support for returning the molecules that did not featurize well. Useful when we need to know which one failed
* Support for more compression formats when reading SDF files using `fssep.open(..., compression="infer")`
