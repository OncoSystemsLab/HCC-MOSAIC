import sys
import numpy as np
import pandas as pd
import scipy
import sklearn
import torch

print("Python:", sys.executable, sys.version)
for module in [np, pd, scipy, sklearn, torch]:
    print(module.__name__, module.__version__)
assert np.allclose(torch.from_numpy(np.array([1., 2.], dtype=np.float32)).numpy(), [1, 2])
torch.set_num_threads(8)
torch.use_deterministic_algorithms(True)
print("Device: CPU; NumPy/Torch bridge: PASS")
