**Added:**

* Support for selfies<2.0.0 in tests

**Changed:**

* Behaviour of all *inchi* functions to return None with a warning instead of silently returning an empty string
* Order of str evaluation on convertion function. `isinstance(str)` is now evaluated before `is None`

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Bug in unique_id making this evaluation falling back on 'd41d8cd98f00b204e9800998ecf8427e' on unsupported inputs. Instead None is returned now

**Security:**

* <news item>
