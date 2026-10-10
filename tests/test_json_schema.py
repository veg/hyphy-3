import os
import sys
import json
import pytest

try:
    import jsonschema
except ImportError:
    pytest.skip("jsonschema not installed", allow_module_level=True)

# Add python binding path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "python")))
import _hyphy3 as hp

def test_fel_modern_and_legacy_json_schema():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    aln_path = os.path.join(repo_root, "..", "hyphy", "hyphy-next", "benchmarks", "data", "cd2.fna")
    tree_path = os.path.join(repo_root, "..", "hyphy", "hyphy-next", "benchmarks", "data", "cd2_fitted.nwk")
    schema_path = os.path.join(repo_root, "schemas", "v3", "fel.v3.schema.json")

    assert os.path.exists(aln_path), f"Missing alignment: {aln_path}"
    assert os.path.exists(tree_path), f"Missing tree: {tree_path}"
    assert os.path.exists(schema_path), f"Missing schema: {schema_path}"

    with open(schema_path) as sf:
        schema = json.load(sf)

    aln = hp.Alignment.from_fasta(aln_path)
    tree = hp.Tree.from_newick_file(tree_path)

    fel = hp.FELAnalyzer(tree, aln)
    fel.run(False, False)

    # 1. Test Modern JSON (v3.0)
    modern_json_str = fel.to_modern_json()
    modern = json.loads(modern_json_str)

    # Validate against formal JSON Schema
    jsonschema.validate(instance=modern, schema=schema)

    # Check key versioning and structure
    assert modern["$schema"] == "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/fel.v3.schema.json"
    assert modern["schema_version"] == "3.0.0"
    assert modern["analysis"]["id"] == "fel"
    assert modern["software"]["name"] == "hyphy"
    assert modern["software"]["version"] == "3.0.0"
    assert len(modern["dataset"]["taxa"]) == 10
    assert modern["dataset"]["codon_sites"] == 187

    # Check columnar site results
    site_data = modern["site_results"]["data"]
    assert len(site_data["site"]) == 187
    assert len(site_data["alpha"]) == 187
    assert len(site_data["beta"]) == 187
    assert len(site_data["lrt"]) == 187
    assert len(site_data["p_value"]) == 187
    assert len(site_data["classification"]) == 187

    # 2. Test Legacy JSON Backward Compatibility
    legacy_json_str = fel.to_legacy_json()
    legacy = json.loads(legacy_json_str)

    assert "MLE" in legacy
    assert "headers" in legacy["MLE"]
    assert "content" in legacy["MLE"]
    assert "0" in legacy["MLE"]["content"]
    assert len(legacy["MLE"]["content"]["0"]) == 187
    assert "branch attributes" in legacy
    assert "fits" in legacy
    assert "input" in legacy

    # 3. Test Parity between Modern Columnar and Legacy Array-of-Arrays
    for i in range(187):
        leg_row = legacy["MLE"]["content"]["0"][i]
        assert abs(site_data["alpha"][i] - leg_row[0]) < 1e-9
        assert abs(site_data["beta"][i] - leg_row[1]) < 1e-9
        assert abs(site_data["alpha_null"][i] - leg_row[2]) < 1e-9
        assert abs(site_data["lrt"][i] - leg_row[3]) < 1e-9
        assert abs(site_data["p_value"][i] - leg_row[4]) < 1e-9

    print("All Modern & Legacy JSON Schema tests passed successfully!")

if __name__ == "__main__":
    test_fel_modern_and_legacy_json_schema()
