using ImputationQualityScores
using Base.Test

@test dosage_to_likelihoods(0.) == (1.,0.,0.) 
@test dosage_to_likelihoods(1.) == (0.,1.,0.) 
@test dosage_to_likelihoods(2.) == (0.,0.,1.) 

@test likelihoods_to_dosage(1.,0.,0.) == (0.)  
@test likelihoods_to_dosage(0.,1.,0.) == (1.)  
@test likelihoods_to_dosage(1.,0.,1.) == (2.)  

for i in 0:0.01:2
	@test i ≈ likelihoods_to_dosage(dosage_to_likelihoods(i))
end

@test proper_info([(0.3333333, 0.3333333, 0.3333333), (1.0, 0.0, 0.0)])≈0.428571367347

