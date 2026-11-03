# Stub: reports every weights file as cached (this module's own path, so no data file is needed);
# config files are reported missing, so the helper falls back to an empty config. Ignores cache_dir.
def try_to_load_from_cache(repo_id, filename, cache_dir=None):
    return None if filename.endswith(".json") else __file__
def hf_hub_download(*a, **k):
    raise RuntimeError("stub: no network")
