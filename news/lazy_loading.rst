**Added:**

* <news item>

**Changed:**

* All the datamol modules and objects are now lazy loaded. It means that loading now happens on-demand. Preliminary tests suggest the datamol import time decreases by 20-fold (from 1s to 50ms on a regular Ubuntu laptop) without affecting the subsequent calls to the modules and objects. This is a major improvement for the datamol usability. This new behaviour is enabled by default but can be disabled by setting the environment variable `DATAMOL_DISABLE_LAZY_LOADING` to `1`.
* Move the fs module to its dedicated section in the docs. Fix #160.

**Deprecated:**

* <news item>

**Removed:**

* Remove unused, broken and uncovered `datamol.fragment.assemble_fragment_iter()` function.

**Fixed:**

* <news item>

**Security:**

* <news item>
