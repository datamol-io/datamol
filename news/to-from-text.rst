**Added:**

* Add `to_smi` and `from_smi` in the IO module.
* Support filelike object in io module.

**Changed:**

* Switch tests of the IO module to regular functions.

**Deprecated:**

* In the IO module, use `urlpath` instead of `file_uri` to follow `fsspec` conventions.

**Removed:**

* <news item>

**Fixed:**

* Fix bug in read_excel where sheet_name wasnt being used.

**Security:**

* <news item>
