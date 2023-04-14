import subprocess
import time
import numpy as np

from tqdm.auto import tqdm

n = 5
duration_list = []
for i in tqdm(range(n)):
    start = time.time()
    subprocess.call(["python", "-c", "import datamol"])
    duration = time.time() - start
    duration_list.append(duration)

duration_average = np.mean(duration_list)
duration_std = np.std(duration_list)

print(f"duration: {duration_average:.2f}±{duration_std:.2f}s")
