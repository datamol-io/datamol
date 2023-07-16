# Usage

## How to use

Datamol has been designed to be used with a single import:

```python
import datamol as dm
```

All `datamol` functions are available under `dm`.

## Lazy loading

datamol uses lazy loading to dynamically expose all its API without imposing a long import time during `import datamol as dm`. In case of trouble you can always disable lazy loading by setting the environment variable `DATAMOL_DISABLE_LAZY_LOADING` to `1`. Please report any issue [on the datamol repo](https://github.com/datamol-io/datamol/issues).
