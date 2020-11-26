from typing import Any

hascEntropy: int

def PyInfoEntropy(results: Any): ...
def PyInfoGain(varMat: Any): ...

InfoEntropy: Any
InfoGain: Any
InfoEntropy = PyInfoEntropy
InfoGain = PyInfoGain
