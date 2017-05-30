# Eventually, I wnat it to be a module, but for now we will just create some test cases
#using MergImp
#using Base.Test

bases = [ "A", "G", "C", "T" ]
function random_base_pair()
	pair=[ bases[rand(1:4)] for i in 1:2 ]
	while first(pair)==last(pair)
		pair=[ bases[rand(1:4)] for i in 1:2 ]
	end
	return pair
end

random_genotype() = join(map(x->@sprintf("%0.3f",x),normalize([rand(),rand(),rand()],1)),":")

generate_snp(position) = vcat(position,random_base_pair() )
generate_individual(i) = "indv$i"
generate_n(func,n) = [ func(i) for i in 1:n ]


vcf_header = "##fileformat=VCFv4.1\n##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"...\">"
vcf_header_line(individuals) = join(vcat(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"],individuals),'\t')
vcf_line(snp,individuals) = join(vcat("42",snp,join(snp,'_'),".",".","GP",[random_genotype() for i in individuals]),'\t')

generate_imputed_vcf(snps, individuals) = join(vcat(vcf_header, vcf_header_line(individuals), [ vcf_line(snp,individuals) for snp in snps ]),'\n')



