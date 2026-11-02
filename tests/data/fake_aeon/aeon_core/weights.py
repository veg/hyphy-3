# Fake aeon_core.weights. FAKE_AEON_API=new exposes only HF_HUB_CACHE (aeon-core 0.2.0);
# anything else exposes only CACHE_DIR (aeon-core 0.1.x).
import os
from pathlib import Path

DEFAULT_VARIANT = "general"
HF_REPO_ID = "datamonkey/hyphaeon"

if os.environ.get("FAKE_AEON_API") == "new":
    HF_HUB_CACHE = Path("/nonexistent/huggingface/hub")
else:
    CACHE_DIR = Path("/nonexistent/aeon")


def get_variant_filename(v):
    return f"model.{v}.pt"


def load_weights(weights, map_location=None):
    return {}
