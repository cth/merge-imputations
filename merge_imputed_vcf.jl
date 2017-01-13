#!/home/fng514/bin/julia
#$ -S /home/fng514/bin/julia
#$ -cwd 
# Christian Theil Have, 2016.

using ArgParse

function process_info(old_info1, old_info2, calculated_info_score)
	old_info1
end

@everywhere function process_genotypes( fields1, fields2, keep, dict1, dict2)
	gpidx = first(find(x -> x=="GP", split(fields1[9],':')))

	maxgp(x) = foldl(max,map(x->parse(Float64,x), split(split(x,':')[gpidx],',')))

	best_genotypes = Array{String,1}()

	for indv in keep
		max_gp_value = 0.0
		genotype_field = nothing
		try 
			tmp_max_gp = maxgp(fields1[dict1[indv]])
			if max_gp_value < maxgp(fields1[dict1[indv]])
				max_gp_value = tmp_max_gp 
				genotype_field = fields1[dict1[indv]]
			end
		end

		try 
			tmp_max_gp = maxgp(fields2[dict2[indv]])
			if max_gp_value < maxgp(fields2[dict2[indv]])
				max_gp_value = tmp_max_gp 
				genotype_field = fields2[dict2[indv]]
			end
		end

		if genotype_field != nothing
			push!( best_genotypes, genotype_field)
		end
	end

	@assert length(keep) == length(best_genotypes)

	best_genotypes
end

# FIXME: Stub for now 
@everywhere calculate_info_score(genotypes) = 1.0

# We know that the files have exactly the same number of header lines
@everywhere function process_line(line1, line2, keep=nothing, dict1=nothing, dict2=nothing)
	if (ismatch(r"^#CHROM", line1) && ismatch(r"^#CHROM", line2))
		fields1 = split(chomp(line1))
		fields2 = split(chomp(line2))

		if (keep == nothing)
			keep = unique(vcat(fields1[10:end],fields2[10:end]))
		else
			keep = intersect( keep, unique(vcat(fields1[10:end],fields2[10:end])))
		end

		write(STDERR, "Keeping $(length(keep)) individuals\n") 

		# reverse lookup positions
		#dict1 = Dict(fields[i] => i for i in fields1)
		#dict1 = Dict(fields[i] => i for i in fields2)
 		dict1 = [ fields1[i] => i for i in 1:length(fields1) ]
		dict2 = [ fields2[i] => i for i in 1:length(fields2) ]

		return (join(foldl(vcat,[fields1[1:9],keep,['\n']]),'\t'),keep,dict1,dict2)

	elseif (ismatch(r"^#", line1) && ismatch(r"^#", line2))
		return (line1,keep,dict1,dict2)
	else
		fields1 = split(line1)
		fields2 = split(line2)

		# Make sure that we are processing the same genotype
		for i in 1:7
			@assert fields1[i] == fields2[i]
		end
		# INFO fields fields1[8] != fields2[8] 
		@assert fields1[9] == fields2[9]

		new_genotypes = process_genotypes(fields1, fields2, keep, dict1, dict2)
		new_info = process_info(fields1[8], fields2[8], calculate_info_score(new_genotypes))

		return (join(foldl(vcat,[fields1[1:7],[new_info],new_genotypes,['\n']]),'\t'),keep,dict1,dict2)
	end
end

# FIXME: We may need to do this with gzip readers instead
function main(args)
	s = ArgParseSettings("Example 2 for merge.jl: " *  # description
		"flags, options help, " *
		"required arguments.")

	@add_arg_table s begin
		"vcf1"
		"vcf2"
		"--keep" 
		help = "A file with a list of ids to keep. One id per line."
	end

	parsed_args = parse_args(s) # the result is a Dict{String,Any}
	println("Parsed args:")
	for (key,val) in parsed_args
		println("  $key  =>  $(repr(val))")
	end

	println(dump(parsed_args))

	if parsed_args["keep"] != nothing
		keep = map(chomp, readstrings(parsed_args["keep"]))
	else
		keep = nothing
	end

	vcf1 = open(parsed_args["vcf1"])
	vcf2 = open(parsed_args["vcf2"])

	lines1 = eachline(vcf1)
	lines2 = eachline(vcf2)

	state1 = start(lines1)
	state2 = start(lines2)

	keep_positions1 = nothing
	keep_positions2 = nothing

	while !done(lines1, nothing) && !done(lines2,nothing)
		(current_line1, state1) = next(lines1, state1)
		(current_line2, state2) = next(lines2, state2)

		(processed_line, keep, keep_positions1, keep_positions2) =  process_line(current_line1,current_line2, keep, keep_positions1,keep_positions2) 

		print(processed_line)
	end

	close(vcf1)
	close(vcf2)
end

main(ARGS)

