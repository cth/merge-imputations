### merge-imputations - Merging of imputed data. 

MergeImp
Mimp.jl
MagicImp.jl

MergImp

ImputedGenotypes

Merges a selected set of individuals from a pair of imputations on the same panel.
If the same individual is present in both imputations, then the script will
select the best-possible imputed genotypes for the individuals from either imputation 
on a per genotype basis. INFO scores as well as AN (allele number) and (AC) allele 
count are recalculated and the output VCF file is written with the same fields.
 
It is not as generic as I would like it to be. Currently, designed to work with the
format produced by Sangers imputation server. Adaptation will be needed for other formats.

TODO list:

- Support multiple (sanger, michigan) imputation formats
- support missing genotypes
- support genotypes where ref/alt alleles are swapped. 
- As a Julia package as well as standalone

Guide: 

- If you have no or very few overlapping samples, choose best
- If you partially overlapping samples choose weighted mean
- If you have completely overlapping sampling you should choose mean


