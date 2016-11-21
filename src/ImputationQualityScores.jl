module ImputationQualityScores
	export dosage_to_likelihoods,likelihoods_to_dosage, rsquared_hat, proper_info

	using Distributions
	
	"The Dosage represents the predicted dosage of the non reference `bb` allele given the data available, 
	it will always have a value between 0 and 2.
	The formula is Dosage = Pr(Het|Data) + 2*Pr(Alt|Data)"
	function likelihoods_to_dosage(p_aa, p_ab, p_bb)
		p_bb * 2 + p_ab
	end

	likelihoods_to_dosage{T<:AbstractFloat}(triple::Tuple{T,T,T}) = likelihoods_to_dosage(triple[1],triple[2],triple[3])

	# Note that there is a one-to-infinite correspondence between doages and genotype likelihoods.
	# This is just one of many possible conversions
	function dosage_to_likelihoods(dosage)
		β=Beta(1,2)
		(pdf(β,dosage)/2, pdf(β,abs(dosage-1))/2, pdf(β,abs(dosage-2))/2)
	end

	function minor_allele_frequency(dosages)
		μ=mean(dosages) 
		(μ > 1.0) ? 2-μ : μ
	end

	function p_hwe(likelihoods)
	end

	# http://www.pypedia.com/index.php/imputation_r_square_hat
	"Calculates the (estimated) fraction of variance in unobserved 0/1/2 genotype explained by the the individual mean genotypes."
	function rsquared_hat(dosages::Array{Float64,1}) 
		p = mean(dosages) / 2
		return var(dosages) / 2 * p * (1-p)
	end

	# Adapted from Python code from 
	function rsquared_hat(likelihoods::Array{Tuple{Float64,Float64,Float64},1})
		rsquared_hat(likelihoods_to_dosage.(likelihoods))
	end
	
	# Adapted from http://www.pypedia.com/index.php/imputation_r_square_hat
	function proper_info(likelihoods)
		n = length(likelihoods)
		meanAB = sum([l[2] for l in likelihoods])
		meanBB = sum([l[3] for l in likelihoods])
		sumX2 = sum([x[2] + (4.0 *x[3]) for x in likelihoods])
		sumXbar2 = sum([(x[2] + (2.0 * x[3]))^2.0 for x in likelihoods])
		meanX = (meanAB+2.*meanBB)/n
		rSqHat = (sumXbar2/n-(meanX^2))/(sumX2/n-(meanX^2))
	end

	# Helli

end # module
