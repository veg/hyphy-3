RequireVersion ("2.5.0");

LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/genetic_code.bf");
LoadFunctionLibrary("libv3/models/DNA/GTR.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");

alignment_file = "/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2_reduced.fna";
tree_file      = "/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2_reduced.nwk";

DataSet ds = ReadDataFile (alignment_file);
code_info = alignments.LoadGeneticCode ("Universal");
DataSetFilter dsf_codon = CreateFilter (ds, 3, "", "", "TAA,TAG,TGA");

fscanf (tree_file, "String", tree_string);

model_mg94 = model.generic.DefineModel("models.codon.MG_REV.ModelDescription", "mg94", {"0": "terms.global", "1": code_info[terms.code]}, "dsf_codon", null);
ExecuteCommands ("UseModel (" + model_mg94[terms.id] + ");");
Tree tree_codon = tree_string;

LikelihoodFunction lf_mg94 = (dsf_codon, tree_codon);

// Fix MG94 parameters
mg94.omega = 0.5;
mg94.theta_AC = 1.0;
mg94.theta_AT = 1.0;
mg94.theta_CG = 1.0;
mg94.theta_CT = 2.0;
mg94.theta_GT = 1.0;

LFCompute (lf_mg94, LF_START_COMPUTE);
LFCompute (lf_mg94, lf_val_mg94);
LFCompute (lf_mg94, LF_DONE_COMPUTE);

fprintf (stdout, "HyPhy 2.5 Total MG94 LogL = ", lf_val_mg94, "\n");

// Site-by-site log likelihoods
ConstructCategoryMatrix (site_l, lf_mg94, COMPLETE);
fprintf (stdout, "HyPhy 2.5 Site Log-Likelihoods:\n", site_l, "\n");

// Branch lengths
bls = BranchLength (tree_codon, -1);
fprintf (stdout, "HyPhy 2.5 Branch lengths:\n", bls, "\n");

// Branch names
b_names = BranchName (tree_codon, -1);
fprintf (stdout, "HyPhy 2.5 Branch names:\n", b_names, "\n");

// Codon frequencies
efv = model_mg94[terms.efv];
fprintf (stdout, "HyPhy 2.5 First 5 Codon freqs:\n", efv[0], ", ", efv[1], ", ", efv[2], ", ", efv[3], ", ", efv[4], "\n");
