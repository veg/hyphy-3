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

    print("FEL Modern & Legacy JSON Schema tests passed successfully!")

def test_meme_modern_and_legacy_json_schema():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    aln_path = os.path.join(repo_root, "..", "hyphy", "hyphy-next", "benchmarks", "data", "cd2.fna")
    tree_path = os.path.join(repo_root, "..", "hyphy", "hyphy-next", "benchmarks", "data", "cd2_fitted.nwk")
    schema_path = os.path.join(repo_root, "schemas", "v3", "meme.v3.schema.json")

    assert os.path.exists(schema_path), f"Missing schema: {schema_path}"
    with open(schema_path) as sf:
        schema = json.load(sf)

    aln = hp.Alignment.from_fasta(aln_path)
    tree = hp.Tree.from_newick_file(tree_path)

    meme = hp.MEMEAnalyzer(tree, aln, hp.MG94Parameters())
    meme.run(False, False)

    # 1. Modern JSON
    modern = json.loads(meme.to_modern_json())
    jsonschema.validate(instance=modern, schema=schema)
    assert modern["schema_version"] == "3.0.0"
    assert modern["analysis"]["id"] == "meme"
    assert len(modern["site_results"]["data"]["site"]) == 187
    assert len(modern["site_results"]["data"]["beta_plus"]) == 187

    # 2. Legacy JSON
    legacy = json.loads(meme.to_legacy_json())
    assert "MLE" in legacy
    assert "0" in legacy["MLE"]["content"]
    assert len(legacy["MLE"]["content"]["0"]) == 187

    # 3. Parity
    for i in range(187):
        assert abs(modern["site_results"]["data"]["alpha"][i] - legacy["MLE"]["content"]["0"][i][0]) < 1e-9
        assert abs(modern["site_results"]["data"]["beta1"][i] - legacy["MLE"]["content"]["0"][i][1]) < 1e-9
        assert abs(modern["site_results"]["data"]["p1"][i] - legacy["MLE"]["content"]["0"][i][2]) < 1e-9
        assert abs(modern["site_results"]["data"]["beta_plus"][i] - legacy["MLE"]["content"]["0"][i][3]) < 1e-9
        assert abs(modern["site_results"]["data"]["p_plus"][i] - legacy["MLE"]["content"]["0"][i][4]) < 1e-9
        assert abs(modern["site_results"]["data"]["lrt"][i] - legacy["MLE"]["content"]["0"][i][5]) < 1e-9
        assert abs(modern["site_results"]["data"]["p_value"][i] - legacy["MLE"]["content"]["0"][i][6]) < 1e-9

    print("MEME Modern & Legacy JSON Schema tests passed successfully!")

def test_busted_modern_and_legacy_json_schema():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    aln_path = os.path.join(repo_root, "..", "hyphy", "hyphy-next", "benchmarks", "data", "cd2.fna")
    tree_path = os.path.join(repo_root, "..", "hyphy", "hyphy-next", "benchmarks", "data", "cd2.nwk")
    schema_path = os.path.join(repo_root, "schemas", "v3", "busted.v3.schema.json")

    assert os.path.exists(schema_path), f"Missing schema: {schema_path}"
    with open(schema_path) as sf:
        schema = json.load(sf)

    aln = hp.Alignment.from_fasta(aln_path)
    tree = hp.Tree.from_newick_file(tree_path)

    busted = hp.BUSTEDAnalyzer.create_and_fit(tree, aln)
    res = busted.run()

    # 1. Modern JSON
    modern = json.loads(busted.to_modern_json(res))
    jsonschema.validate(instance=modern, schema=schema)
    assert modern["schema_version"] == "3.0.0"
    assert modern["analysis"]["id"] == "busted"
    assert len(modern["site_results"]["data"]["site"]) == 187
    assert len(modern["site_results"]["data"]["evidence_ratio"]) == 187

    # 2. Legacy JSON
    legacy = json.loads(busted.to_legacy_json(res))
    assert "test results" in legacy
    assert "Evidence Ratios" in legacy
    assert "fits" in legacy

    # 3. Parity
    assert abs(modern["statistical_tests"]["gene_level"]["statistic_value"] - legacy["test results"]["LRT"]) < 1e-9
    assert abs(modern["statistical_tests"]["gene_level"]["p_value"] - legacy["test results"]["p-value"]) < 1e-9
    for i in range(187):
        assert abs(modern["site_results"]["data"]["evidence_ratio"][i] - legacy["Evidence Ratios"]["optimized null"][0][i]) < 1e-9
        assert abs(modern["site_results"]["data"]["unconstrained_log_l"][i] - legacy["Site Log Likelihood"]["unconstrained"][0][i]) < 1e-9
        assert abs(modern["site_results"]["data"]["constrained_log_l"][i] - legacy["Site Log Likelihood"]["constrained"][0][i]) < 1e-9

    print("BUSTED Modern & Legacy JSON Schema tests passed successfully!")

def test_absrel_modern_and_legacy_json_schema():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    aln_path = os.path.join(repo_root, "..", "hyphy", "hyphy-next", "benchmarks", "data", "cd2.fna")
    tree_path = os.path.join(repo_root, "..", "hyphy", "hyphy-next", "benchmarks", "data", "cd2.nwk")
    schema_path = os.path.join(repo_root, "schemas", "v3", "absrel.v3.schema.json")

    assert os.path.exists(schema_path), f"Missing schema: {schema_path}"
    with open(schema_path) as sf:
        schema = json.load(sf)

    aln = hp.Alignment.from_fasta(aln_path)
    tree = hp.Tree.from_newick_file(tree_path)

    settings = hp.ABSRELSettings()
    settings.max_rate_classes = 1
    settings.p_threshold = 0.05

    absrel = hp.ABSRELAnalyzer.create(tree, aln, settings)
    res = absrel.run()

    # 1. Modern JSON
    modern = json.loads(res.to_modern_json(tree, aln))
    jsonschema.validate(instance=modern, schema=schema)
    assert modern["schema_version"] == "3.0.0"
    assert modern["analysis"]["id"] == "absrel"
    assert len(modern["branch_results"]["data"]["branch"]) == len(res.tested_branches)

    # 2. Legacy JSON
    legacy = json.loads(res.to_legacy_json(tree, aln))
    assert "test results" in legacy
    assert "fits" in legacy
    assert "branch attributes" in legacy

    # 3. Parity
    assert abs(modern["model_fits"]["baseline_mg94"]["log_likelihood"] - legacy["fits"]["Baseline MG94xREV"]["Log Likelihood"]) < 1e-9
    assert abs(modern["model_fits"]["full_adaptive"]["log_likelihood"] - legacy["fits"]["Full adaptive model"]["Log Likelihood"]) < 1e-9

    print("aBSREL Modern & Legacy JSON Schema tests passed successfully!")

def test_relax_modern_and_legacy_json_schema():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    candidates = [
        os.path.join(repo_root, "..", "hyphy", "tests", "data", "Fig4E.nex"),
        os.path.join(repo_root, "tests", "data", "Fig4E.nex"),
    ]
    aln_path = next((p for p in candidates if os.path.exists(p)), candidates[0])
    schema_path = os.path.join(repo_root, "schemas", "v3", "relax.v3.schema.json")

    assert os.path.exists(schema_path), f"Missing schema: {schema_path}"
    with open(schema_path) as sf:
        schema = json.load(sf)

    aln = hp.Alignment.load(aln_path)
    tree = hp.Tree.from_newick(aln.embedded_tree_newick)

    settings = hp.RELAXSettings()
    settings.refine_branch_lengths = False

    relax = hp.RELAXAnalyzer.create(tree, aln, settings)
    res = relax.run()

    # 1. Modern JSON
    modern = json.loads(res.to_modern_json(tree, aln))
    jsonschema.validate(instance=modern, schema=schema)
    assert modern["schema_version"] == "3.0.0"
    assert modern["analysis"]["id"] == "relax"
    assert "hypothesis_test" in modern["statistical_tests"]
    assert "k_parameter" in modern["statistical_tests"]["hypothesis_test"]

    # 2. Legacy JSON
    legacy = json.loads(res.to_legacy_json(tree, aln))
    assert "test results" in legacy
    assert "fits" in legacy
    assert "branch attributes" in legacy

    # 3. Parity
    assert abs(modern["statistical_tests"]["hypothesis_test"]["statistic_value"] - legacy["test results"]["LRT"]) < 1e-9
    assert abs(modern["statistical_tests"]["hypothesis_test"]["p_value"] - legacy["test results"]["p-value"]) < 1e-9
    assert abs(modern["statistical_tests"]["hypothesis_test"]["k_parameter"] - legacy["test results"]["relaxation or intensification parameter"]) < 1e-9

    print("RELAX Modern & Legacy JSON Schema tests passed successfully!")

if __name__ == "__main__":
    test_fel_modern_and_legacy_json_schema()
    test_meme_modern_and_legacy_json_schema()
    test_busted_modern_and_legacy_json_schema()
    test_absrel_modern_and_legacy_json_schema()
    test_relax_modern_and_legacy_json_schema()
