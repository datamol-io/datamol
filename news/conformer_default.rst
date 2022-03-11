**Added:**

* Add `dm.conformers.keep_conformers` in order to only keep one or multiple conformers from a molecules.

**Changed:**

* Change the conformer generation arguments to use `useRandomCoords=True` by default.
* Start using explicit `Optional` instead of implicit `Optional` for typing.
* Start using relative imports instead of absolute ones.
* When conformers are not minimized, sort them by energy (can be turned to False).

**Deprecated:**

* <news item>

**Removed:**

* Remove `fallback_to_random_coords` argument from `generate_conformers`.

**Fixed:**

* <news item>

**Security:**

* <news item>
