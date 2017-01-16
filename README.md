### merge-imputations - Merging of imputed data. 

Merges a selected set of individuals from a pair of imputations on the same panel.
If the same individual is present in both imputations, then the script will
select the best-possible imputed genotypes for the individuals from either imputation 
on a per genotype basis. INFO scores as well as AN (allele number) and (AC) allele 
count are recalculated and the output VCF file is written with the same fields.
 
It is not as generic as I would like it to be. Currently, designed to work with the
format produced by Sangers imputation server. Adaptation will be needed for other formats.
