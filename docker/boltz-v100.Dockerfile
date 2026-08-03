FROM ovo-boltz:latest

# PyTorch 2.11/cu130 wheels omit sm_70 kernels. Boltz 2.2.1 requires
# torch>=2.2, so pin the official cu121 build that retains V100 support.
RUN python3 -m pip install --no-cache-dir --force-reinstall \
    torch==2.5.1 \
    --index-url https://download.pytorch.org/whl/cu121

RUN python3 - <<'PY'
import importlib.metadata as metadata
import torch

assert torch.__version__.startswith("2.5.1"), torch.__version__
assert torch.version.cuda == "12.1", torch.version.cuda
assert metadata.version("boltz") == "2.2.1"
print(
    {
        "torch": torch.__version__,
        "torch_cuda": torch.version.cuda,
        "boltz": metadata.version("boltz"),
    }
)
PY
