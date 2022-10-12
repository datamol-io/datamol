**Added:**

* <news item>

**Changed:**

* Bump upstream GH actions versions.
* `dm.fs.copy_dir` now uses the internal fsspec `copy` when the two source and destination fs are the same. It makes the copy much faster.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Use `os.PathLike` to recognize a broader range of string-based path inputs in the `dm.fs` module. It prevents file objects such as `py._path.local.LocalPath` not being recognized as path.

**Security:**

* <news item>
