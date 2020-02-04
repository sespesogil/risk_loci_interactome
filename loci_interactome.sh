
module load bedtools/2.24.0
module load circos
module load homer/4.10

SNPs=   # bedfile of SNPs/risk loci of interest
Hic=  #interaction frequencies from HiC-Pro
REF=  # to genome fastq
samples=("cluster1" "cluster2" "cluster3" "cluster4" "cluster5" "cluster6" "cluster7" "cluster8" "cluster9" "cluster10" "cluster11")

awk -F'\t' '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$7}' $Hic > HiC_targets.bed

length=${#samples[@]}
echo ${length}

for (( i=0; i<${length}; i++ ))
do
        varname=${samples[${i}]}

        bedtools intersect -a $SNPs -b $varname.bed -wa > $SNPs.$varname.bed

# should include intersection of SNPs only in cluster1 


bedtools intersect -a $SNPs.$varname.bed -b $Hic -wo > temp_a  #FORWARD

awk -F'\t' '{print $7"\t"$8"\t"$9"\t"$4"\t"$5"\t"$6"\t"$10}' temp_a > temp_a2

bedtools intersect -a $varname.bed  -b temp_a2 -wo > temp_a3

# results= SNPS per cluster [1-3] / anchor [4-6] / target [7-9]/ freq [10]

bedtools intersect -a $SNPs.$varname.bed -b HiC_targets.bed -wo > temp_b  #REVERSE

awk -F'\t' '{print $7"\t"$8"\t"$9"\t"$4"\t"$5"\t"$6"\t"$10}' temp_b > temp_b2

bedtools intersect -a $varname.bed  -b temp_b2 -wo > temp_b3

awk -F'\t' '{print $7"\t"$8"\t"$9"\t"$10}' temp_b3 | cat -n | awk -F' ' '{print "interaction"$1"\t"$2"\t"$3"\t"$4"\t""thickness="$5}' > $SNPs.$varname.connectome_inside.anchor1.final.bed
awk -F'\t' '{print $4"\t"$5"\t"$6"\t"$10}' temp_b3 | cat -n | awk -F' ' '{print "interaction"$1"\t"$2"\t"$3"\t"$4"\t""thickness="$5}' > $SNPs.$varname.connectome_inside.target1.final.bed
awk -F'\t' '{print $7"\t"$8"\t"$9"\t"$10}' temp_a3 | cat -n | awk -F' ' '{print "interaction"$1"\t"$2"\t"$3"\t"$4"\t""thickness="$5}' > $SNPs.$varname.connectome_inside.anchor2.final.bed
awk -F'\t' '{print $4"\t"$5"\t"$6"\t"$10}' temp_a3 | cat -n | awk -F' ' '{print "interaction"$1"\t"$2"\t"$3"\t"$4"\t""thickness="$5}' > $SNPs.$varname.connectome_inside.target2.final.bed

cat $SNPs.$varname.connectome_inside.target1.final.bed $SNPs.$varname.connectome_inside.anchor1.final.bed | sort > $SNPs.$varname.connectome_inside1.interactions.bed

cat $SNPs.$varname.connectome_inside.target2.final.bed $SNPs.$varname.connectome_inside.anchor2.final.bed | sort > $SNPs.$varname.connectome_inside2.interactions.bed

awk -F"\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' $SNPs.$varname.connectome_inside1.interactions.bed > $SNPs.$varname.connectome_inside.interactions1.forinteraction.bed
awk -F"\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' $SNPs.$varname.connectome_inside2.interactions.bed > $SNPs.$varname.connectome_inside.interactions2.forinteraction.bed


# mkdir motif_analysis
# findMotifsGenome.pl motif_discovery_peaks.$varname.bed $REF motif_analysis/motif_analysis.$varname
# creating the karyotype

sort -k1,1 -k2,2n $varname.bed > $varname.sorted.bed

awk -v OFS='\t' '$1 in min{min[$1]=$2<min[$1]?$2:min[$1]; max[$1]=$3>max[$3]?$3:max[$1] ; next } {min[$1]=$2;max[$1]=$3}; END{for(x in min)print x, min[x], max[x]}' $varname.sorted.bed > $varname.uniq.bed

awk -F"\t" '{print "chr - "" "$1"\t"$1"\t"$2"\t"$3"\t"$1}' $varname.uniq.bed > $varname.karyotype.bed

awk -F"\t" '{print $1}' $varname.bed | uniq > chromosomes.bed

cut -f1 chromosomes.bed | paste -s  > t_chromosomes

sed -i.bak $'s/\t/;/g' t_chromosomes

pipe_to_conf=$(head t_chromosomes)

cp SNPs.config SNPs.$varname.1.conf
cp SNPs.config SNPs.$varname.2.conf

sed -i "s/chromosomes =/chromosomes = $(head t_chromosomes)/" SNPs.$varname.1.conf
sed -i "s/file =/file = $(echo $SNPs.$varname.connectome_inside1.interactions.bed)/" SNPs.$varname.1.conf
sed -i "s/karyotype =/karyotype = $(echo $varname.karyotype.bed)/" SNPs.$varname.1.conf

circos -conf SNPs.$varname.1.conf


sed -i "s/chromosomes =/chromosomes = $(head t_chromosomes)/" SNPs.$varname.2.conf
sed -i "s/file =/file = $(echo $SNPs.$varname.connectome_inside2.interactions.bed)/" SNPs.$varname.2.conf
sed -i "s/karyotype =/karyotype = $(echo $varname.karyotype.bed)/" SNPs.$varname.2.conf

circos -conf SNPs.$varname.2.conf


done
