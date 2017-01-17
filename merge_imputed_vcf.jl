#!/home/fng514/bin/julia
#$ -S /home/fng514/bin/julia
#$ -cwd 
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
using BGZFStreams
using Memoize

######## Calculation of INFO scores ########
tiny_constant = 10.0^-30
check_probs(pr) = pr
probs(pr::Tuple{Float64,Float64,Float64}) = check_probs([pr...])
probs(pr::Array{Float64,1}) = check_probs(pr)
probs(pr::Array{Tuple{Float64,Float64,Float64},1}) = normalize([foldr((x,y)->(x[1]+y[1],x[2]+y[2], x[3]+y[3]), (.0,.0,.0), pr)...],1)
freqA(g) = 0.5*probs(g)[2] + probs(g)[1]
expected_dosage_variance(g,gprobs) = 4*gprobs[1] + gprobs[2] - (2*freqA(gprobs))^2
observed_dosage_variance(g,gprobs) = sum(map(x->(2*probs(x)[1] + probs(x)[2])^2,g) - (2*freqA(gprobs))^2)  / length(g)
function rsquared_hat(likelihoods::Array{Tuple{Float64,Float64,Float64},1}) 
	gprobs = probs(likelihoods)
	expected = expected_dosage_variance(likelihoods,gprobs)
	observed = observed_dosage_variance(likelihoods,gprobs)
	if expected == observed == 0.0
		1.0
	else
		observed / (expected+tiny_constant)
	end
end

@everywhere function process_info(old_info1, old_info2,  allele_number, allele_count,calc_info_score)
	field1 = split(old_info1,';')[1] # We expect first field is RefPanelAF
	string("$(field1);AN=$(allele_number);AC=$(allele_count);INFO=$(round(calc_info_score,3))")
end

# We know that the files has exactly the same number of header lines
@everywhere function process_line(line1, line2, keep=nothing, dict1=nothing, dict2=nothing)
	if (ismatch(r"^#CHROM", line1) && ismatch(r"^#CHROM", line2))
		fields1 = split(chomp(line1))
		fields2 = split(chomp(line2))

		if (keep == nothing)
			keep = unique(vcat(fields1[10:end],fields2[10:end]))
		else
			keep = intersect( keep, unique(vcat(fields1[10:end],fields2[10:end])))
			#println(keep)
		end
		#println(keep)

		write(STDERR, "Keeping $(length(keep)) individuals\n") 

		# reverse lookup positions
		dict1 = Dict(fields1[i] => i for i in 1:length(fields1))
		dict2 = Dict(fields2[i] => i for i in 1:length(fields2))

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

		#(new_genotypes,new_info) = process_genotypes(fields1, fields2, keep, dict1, dict2, best_genotypes, genotype_likelihoods)

		gpidx = first(find(x -> x=="GP", split(fields1[9],':')))
		gtidx = first(find(x -> x=="GT", split(fields1[9],':')))

		best_genotypes = Array{SubString{String},1}(length(keep))
		genotype_likelihoods = Array{Tuple{Float64,Float64,Float64},1}(length(keep))

		write(STDERR,"processing genotype $(fields1[3])\n")

		# FIXME: This function is depends on phased genotypes
		function count_alleles(x) 
			tmp=split(x,':')[gtidx]
			if tmp == "0|0"
				0
			elseif tmp == "0|1" || tmp == "1|0"
				1
			elseif tmp == "1|1"
				2
			else
				throw("Unknown allele $(tmp)")
			end
		end

		allele_number = length(keep) * 2
		allele_count = 0
		allele_incr = 0

		indv_processed = 0 
		single_genotype_likelihood1=[.0,.0,.0]
		single_genotype_likelihood2=[.0,.0,.0]

		@time for i in 1:length(keep)
			indv = keep[i]
			tmp_max_gp = 0.0
			max_gp_value = 0.0

			if haskey(dict1,indv) 
				likelihood_fields = split(split(fields1[dict1[indv]],':')[gpidx],',')
				println(likelihood_fields)
				for i in 1:3
					single_genotype_likelihood1[i]=parse(Float64,likelihood_fields[i])
					if single_genotype_likelihood1[i] > tmp_max_gp
						tmp_max_gp = single_genotype_likelihood1[i] 
					end
				end
				if max_gp_value < tmp_max_gp
					max_gp_value = tmp_max_gp 
					best_genotypes[i] = fields1[dict1[indv]]
					genotype_likelihoods[i] = (single_genotype_likelihood1...)
				end
			elseif haskey(dict2,indv)	
				likelihood_fields = split(split(x,':')[gpidx],',')
				for i in 1:3
					single_genotype_likelihood2[i]=parse(Float64,likelihood_fields[i])
					if single_genotype_likelihood2[i] > tmp_max_gp
						tmp_max_gp = single_genotype_likelihood2[i] 
					end
				end
				if max_gp_value < tmp_max_gp
					max_gp_value = tmp_max_gp 
					best_genotypes[i] = fields1[dict1[indv]]
					genotype_likelihoods[i] = (single_genotype_likelihood2...)
				end
			end

			allele_count = allele_count + count_alleles(best_genotypes[i])
			indv_processed =  indv_processed + 1
		end

		new_info = process_info(fields1[8], fields2[8], allele_number, allele_count, rsquared_hat(genotype_likelihoods))

		#join(new_genotypes,'\t')

		return (join(foldl(vcat,[fields1[1:7],[new_info],best_genotypes,["\n"]]),'\t'),keep,dict1,dict2)
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
		"--out"
	end

	parsed_args = parse_args(s) # the result is a Dict{String,Any}
	println("Parsed args:")
	for (key,val) in parsed_args
		println("  $key  =>  $(repr(val))")
	end

	if parsed_args["keep"] != nothing
		keep = Array{String,1}()
		open(parsed_args["keep"]) do keepfile
			for line in eachline(keepfile)
				push!(keep,chomp(line))
			end
		end
	else
		keep = nothing
	end

	#vcf1 = open(parsed_args["vcf1"])
	#vcf2 = open(parsed_args["vcf2"])

	vcf1=BGZFStream(parsed_args["vcf1"])
	vcf2=BGZFStream(parsed_args["vcf2"])

	#out = open(parsed_args["out"],"w")
	out = BGZFStream(parsed_args["out"],"w")

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

		write(out,processed_line)
	end

	close(vcf1)
	close(vcf2)
	close(out)
end

main(ARGS)

# TEST: Merging the a vcf with itself should yield itself 
