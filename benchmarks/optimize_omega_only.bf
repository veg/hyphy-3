RequireVersion ("2.5.0");
LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
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

// Fix GTR parameters
mg94.theta_AC := 1.0;
mg94.theta_AT := 1.0;
mg94.theta_CG := 1.0;
mg94.theta_CT := 2.0;
mg94.theta_GT := 1.0;

// Fix all branch lengths to their initial values from the Newick tree
tree_branches = BranchName (tree_codon, -1);
for (i = 0; i < Columns (tree_branches); i += 1) {
    b_name = tree_branches[i];
    ExecuteCommands ("tree_codon." + b_name + ".t := " + Eval("tree_codon." + b_name + ".t") + ";");
}

// Ensure omega is the ONLY free parameter
mg94.omega = 1.0;

Optimize (res, lf_mg94);

fprintf (stdout, "HyPhy 2.5 Optimized omega = ", mg94.omega, "\n");
fprintf (stdout, "HyPhy 2.5 Optimized LogL  = ", res[1][0], "\n");
