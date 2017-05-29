#!/home/fng514/bin/julia #$ -S /home/fng514/bin/julia #$ -cwd # Christian Theil Have, 2017.
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

typealias GenotypeProbabilities Array{Tuple{Float64,Float64,Float64},1}
typealias GenotypeProbability Tuple{Float64,Float64,Float64}

@everywhere probs(pr::Array{Tuple{Float64,Float64,Float64},1}) = normalize([foldr((x,y)->(x[1]+y[1],x[2]+y[2], x[3]+y[3]), (.0,.0,.0), pr)...],1)
@everywhere freqA(g) = 0.5*probs(g)[2] + probs(g)[1]

@everywhere function r²(ps::Array{Tuple{Float64,Float64,Float64},1})
	n = length(ps)

	esum = e2sum = fsum = 0.0

	for l in ps
		norm = sum(l) 
		l = [ i/norm for i in l ]
		esum += l[2] + 2l[3]
		e2sum += (l[2] + 2l[3])^2
		fsum += l[2] + 4l[3]
	end

	θ = esum / 2n

	(1.0 > freqA(ps) > 0.0) ?  abs(1 - (fsum - e2sum) / (2n * θ * (1.0 - θ))) : 1.0
end 

@everywhere likelihoods_to_dosage(p_aa, p_ab, p_bb) = 2p_bb + p_ab
@everywhere likelihoods_to_dosage{T<:AbstractFloat}(triple::Tuple{T,T,T}) = likelihoods_to_dosage(triple[1],triple[2],triple[3])

@everywhere function max_index(p)
	max_index = 0
	max_probability = 0.0
	for i in 1:length(p)
		if p[i] > max_probability
			max_index = i
        		max_probability = p[i]
        	end
	end
	max_index
end

@everywhere max_prob(p) = p[max_index(p)]

putf(f) = @sprintf("%0.3f",f)

alt_allele_count(p::Tuple{Float64,Float64,Float64}) = max_index(p)-1
alt_allele_count(p::Array{Tuple{Float64,Float64,Float64},1}) = sum(map(alt_allele_count,p))

GT(p::Tuple{Float64,Float64,Float64}) = ("0/0","0/1","1/1")[max_index(p)]
DS(p::Tuple{Float64,Float64,Float64}) = putf(likelihoods_to_dosage(p[1],p[2],p[3]))
GP(p::Tuple{Float64,Float64,Float64}) = join(map(putf,p),',')

INFO(l)=string("AC=",length(l)*2,";AN=",alt_allele_count(l),";R2=", r²(l))

@everywhere best_gp(ps) = ps[max_index(map(max_prob,ps))]

@everywhere mean_gp(ps) = tuple(probs(ps)...)

@everywhere function confidence_weighted_mean_gp(ps)
	cs = map(p -> foldl(max,p) - foldl(min,p), ps)
	tuple(probs([ (cs[i]*ps[i][1], cs[i]*ps[i][2], cs[i]*ps[i][3]) for i in 1:length(ps) ])...)
end

@everywhere best_info(ps) = () # Only a Placeholder

function merge_snp(fields, mergefun) 
	@assert length(fields) > 0

	println(join(fields[1][1:4],' '))

	gpidx = [ first(find(x -> x=="GP"||x=="GL", split(f[9],':'))) for f in fields ] 

	# Create an array Array{Tuple{Float64,Float64,Float64},1} for each individual
	best_likelihoods = []
	let snp_probs(i,indv) = (map(x->parse(Float64,x),split(split(fields[i][indv],':')[gpidx[i]],','))...)
		# Not very multi-dispatchy, but ...
		if mergefun == best_info
			# This needs access to others genotypes as well
			# likelihoods is an array  of Array{Tuple{Float64,Float64,Float64},1} for each method
			likelihoods = [ [ snp_probs(idx,indv) for indv in 10:length(first(fields)) ] for idx in 1:length(fields) ]
			# Select likehoods of method with highest overall INFO score
			best_likelihoods = likelihoods[max_index(map(r², likelihoods))]
		else	
			# Create an array Array{Tuple{Float64,Float64,Float64},1} for each individual:
			best_likelihoods = @parallel vcat for indv in 10:length(fields[1]) 
				# Create an array Array{Tuple{Float64,Float64,Float64},1} for each method and collapse with "mergefun"
				mergefun( filter(isna,[ snp_probs(i,indv) for i in 1:length(fields) ]) )
			end
		end
	end

	join(foldl(vcat,
		[fields[1][1:7],
		INFO(best_likelihoods),
		["GT:DS:GP"],
		map(p->join([GT(p),DS(p),GP(p)],':'), best_likelihoods),
		['\n']]),'\t')
end

vcf_header="""##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]">
##FORMAT=<ID=GP,Number=3,Type=Float,Description="Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 ">
##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy">
"""


# Two vcf lines are considered to represent the same variant if 
# chromosome AND position AND ref+alt allele are the same (but order of alt/ref may be opposite)
function vcfline_same_variant(f1,f2)
	f1[1] == f2[1] && f1[2] == f2[2] && sort(f1[3:4]) == sort(f2[3:4])
end

# Greater-than operator for vcf lines (split into fields) 
function vcfline_before(f1,f2)
	if f1[1] < f2[1] # Chromosome smaller
		true
	elseif f1[1] == f2[1] && f1[2] < f2[2] # Position smaller 
		true
	elseif f1[1] == f2[1] && f1[2] == f2[2] && string(sort(f1[3:4])) < string(sort(f2[3:4])) # allele ordering
		true
	else
		false
	end	
end

revidx(lines) = [ Dict(f[i] => i for i in 10:length(f)) for f in in [ split(line) for line in lines ] ]

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

	vcfs = [ BGZFStream(vcf) for vcf in parsed_args["vcfs"] ]
	out = BGZFStream(parsed_args["out"],"w")

	lines = []
	states = []
	current_lines = []
	for vcf in vcfs
		lines_vcf = eachline(vcf)
		push!(lines,lines_vcf)
		push!(states,start(lines_vcf))
		push!(current_lines,"##")
	end

	# Read past header lines
	for i in 1:length(lines)
		while ismatch(r"^##", current_lines[i])
			current_lines[i], states[i] = next(lines[i], states[i])
		end
	end

	@assert (all(map(line->ismatch(r"^#CHROM", line),current_lines)))

	# Build sample indices
	sample_idx=revidx(current_lines)


	# TODO: 
	# And idea for a strategy that requires minimal changes to merge_snp function
	# would be to "reconstruct" full lines with all samples in a paricular order
	# where some samples may be NA. 
	# 2) As an optimization we may only have to do this when not all files contain
	# all samples in the same order...

	samples_in_order = sort(foldl(union,map(keys,sample_idx)))

	write(out,vcf_header)
	write(out,current_lines[1])

	for i in 1:length(lines)
		current_lines[i], states[i] = next(lines[i], states[i])
	end
	current_fields = [ split(line) for line in current_lines ] 

	while all(map(x->!done(x,nothing), lines))
		# Find smallest position
		smallest =  1
		for i in 1:length(lines)
			if vcfline_before(current_fields[i], current_fields[smallest]) 
				smallest=i
			end
		end

		# Select only lines with smallest position for merging in this step and read next line from those files
		process_lines = []
		for i in 1:length(lines)
			if vcfline_same_variant(current_fields[i], current_fields[smallest])
				current_lines[i], states[i] = next(lines[i], states[i])
				current_fields[i] = split(current_lines[i])
				push!((i,process_lines), current_fields[i])
			end
		end

		# I am going to need the index for each process_lines (To reverse map individuals
		# I need to know hwich file they are from

		# TODO: 
		# And idea for a strategy that requires minimal changes to merge_snp function
		# would be to "reconstruct" full lines with all samples in a paricular order
		# where some samples may be NA. 

		proclines = []
		for (idx,line) in process_lines
			aline=[]
			sidx = sample_idx[idx] 
			for i in 1:length(procline)
				if < 10 # Copy over non-genotype fields
					push!(aline, line[i])
				elseif haskey(samples[idx],samples_in_order[i-10])
					push!(aline,line[sample_idx[idx][samples_in_order[i-10]]])
				else
					push!(aline,NA)	
				end
			end
			push!(proclines,aline)
		end

		write(out,merge_snp(proclines,mergefun))
	end

	# Close open files
	for vcf in vcfs 
		close(vcf)
	end
	close(out)
end

main(ARGS)
