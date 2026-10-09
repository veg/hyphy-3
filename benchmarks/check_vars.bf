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

fprintf (stdout, "Rate matrix key: ", model_mg94[terms.model.rate_matrix], "\n");
ExecuteCommands ("q_name = model_mg94[terms.model.rate_matrix]; Q_mat = ^q_name;");
fprintf (stdout, "Q_mat dimensions: ", Rows(Q_mat), " x ", Columns(Q_mat), "\n");
fprintf (stdout, "Q_mat[0][0]: ", Q_mat[0][0], "\n");
fprintf (stdout, "Q_mat[0][1]: ", Q_mat[0][1], "\n");
fprintf (stdout, "Q_mat[0][2]: ", Q_mat[0][2], "\n");
fprintf (stdout, "Q_mat[0][3]: ", Q_mat[0][3], "\n");
