# Stub: reports every file as cached (a stub JSON file), ignoring cache_dir.
def try_to_load_from_cache(repo_id, filename, cache_dir=None):
    import os
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "config.json")
def hf_hub_download(*a, **k):
    raise RuntimeError("stub: no network")
