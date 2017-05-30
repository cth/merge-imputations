#!/home/fng514/bin/julia 
# Christian Theil Have, 2017.
# Merging of imputed data. 
# 
# Merges a selected set of individuals from a pair of imputations on the same panel.
# If the same individual is present in more both imputations, then the script will
# select the best-possible imputed genotypes for the individuals from either imputation 
# on a per genotype basis. INFO scores as well as AN (allele number) and (AC) allele 
# count is recalculated and the output VCF file is written with the same fields.
# 
# It is not as generic as I would like it to be: 
# Currently, designed to work with the format produced by Sangers imputation server.
# Adaptation will be need for other formats.

using ArgParse
function main(args)
	s = ArgParseSettings("Example for merge.jl: " *  # description 
		"flags, options help, " *
		"required arguments.")

	@add_arg_table s begin
		"--vcfs"
		nargs = '*'
		help = "A list of VCF files to merge"
		"--merge-mode"
		help = "One of best,mean,weightedmean. Default is best."
		default = "best"
		"--out"
	end

	parsed_args = parse_args(s)
	println("Parsed args:")
	for (key,val) in parsed_args
		println("  $key  =>  $(repr(val))")
	end

	@assert length(parsed_args["vcfs"]) > 0


	# Merge method:
	if parsed_args["merge-mode"] == "best"
		mergefun = best_gp
	elseif parsed_args["merge-mode"] == "mean"
		mergefun = mean_gp
	elseif parsed_args["merge-mode"] == "weightedmean"
		mergefun = confidence_weighted_mean_gp
	elseif parsed_args["merge-mode"] == "bestinfo"
		mergefun = best_info	
	else
		throw("Invalid merge-mode : $(parsed_args["merge-mode"])")
	end	

	merged_imputed(parsed_args["vcfs"], parsed_args["out"], mergefun)
end

main(ARGS)
