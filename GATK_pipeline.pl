# usr/bin/perl -w
use strict;
# usage perl GATK_ssr.pl working_directory(no end/)
# change your sleep time
# if you rerun the DBimport step, remove all files produced before
# check log file each step, grep "Runtime" |wc -l should get you the number of your file

# indexed genome
my $ref = "$ARGV[0]/IRGSP-1.0/IRGSP-1.0_genome.fasta";
# make sure that you add at least 400bp upstream or downstream for target region
# the reduce.list contain non-overlapping regions, merged in R beforehand, GeneRegion is exact the position of the gene, region extends 400bp beyond your gene for accuracy purpose in SNP/InDel call
my $region = "$ARGV[0]/IRGSP-1.0/LOC.list";
my $GeneRegion ="$ARGV[0]/Harkamal_GeneLen/Gene.list"; 
# Input Bam files with Input.txt file
my $Bam = "$ARGV[0]/Bam_SortIndex";

open(OUT1,">$ARGV[0]/Code/VariantCall_I.sh");
open(OUT2,">$ARGV[0]/Variant_call/cohort.sample_map");
open(OUT3,">$ARGV[0]/Acc_genome/RemoveStar_VI.sh");
open(OUT4,">$ARGV[0]/Acc_genome/AccGenome_VII.sh");
open(OUT5,">$ARGV[0]/Code/PersoEpcr.sh");

# 1. Creat regional reference genome
open(OUT6,">$ARGV[0]/Code/Extract_genome.sh");
open(IN1,"$GeneRegion");
my @position;
while(<IN1>){
        chomp;
	push @position,$_;
   	print OUT6 "samtools faidx $ref $_ >$ARGV[0]/IRGSP-1.0/ExtractG/$_\n";

}
close(OUT6);
close(IN1);

# 2. Call variant
my $tcontrol = 1;
open(IN,"$Bam/Input.txt"); 
while(<IN>){
	chomp;
	my @array = split(/\./,$_);
	my $file_name = $array[0];
	## 2.1 linux code for variant call
	# regional call
        print OUT1 "nohup gatk HaplotypeCaller -R $ref -I $Bam/$_ -O $ARGV[0]/Variant_call/$file_name.g.vcf -L $region -ERC GVCF --max-alternate-alleles 25 >$ARGV[0]/Variant_call/Log/$file_name.log 2>&1 &\n";
	if(($tcontrol>14)&&($tcontrol%15==0)){
		print OUT1 "sleep 120\n";
	}
	## whole genome call
	# print OUT1 "nohup gatk HaplotypeCaller -R $ref -I $Bam/$_ -O $ARGV[0]/Variant_call/$file_name.gvcf -ERC GVCF --max-alternate-alleles 25 >$ARGV[0]/Variant_call/Log/$file_name.log 2>&1 &\n";
	## 2.2 sample file for DBimport
	print OUT2 "$file_name\t$ARGV[0]/Variant_call/$file_name.g.vcf\n";
	## 2.3 remove rows containing star; if we have alternative read >5 or alternative read > reference read, kept; gunzip produced NStar.vcf; index vcf file
	print OUT3 "nohup ```awk '{if(NR<=54){print \$0}else{if(\$5~/,|\\*/){next}else{split(\$10,a,\":\");num=split(a[2],b,\",\");if((num==1)||(b[2]>b[1])||(b[2]>5)){print \$0}}}}' $file_name.vcf > $file_name.NStar.vcf;bgzip $file_name.NStar.vcf;gatk IndexFeatureFile -F $file_name.NStar.vcf.gz >Log/$file_name.log 2>&1``` &\n";
	if(($tcontrol>14)&&($tcontrol%15==0)){
                print OUT3 "sleep 30\n";
        }
	## 2.4 produce consensus genome, if repeat, should remove the file produced before
	print OUT4 "nohup ```";
	foreach my $pos(@position){
		print OUT4 "bcftools consensus -f $ARGV[0]/IRGSP-1.0/ExtractG/$pos $file_name.NStar.vcf.gz >> $file_name.fa;";
	}
	print OUT4 "``` &\n";
	if(($tcontrol>14)&&($tcontrol%15==0)){
                print OUT4 "sleep 20\n";
        }
	# 2.5 personal epcr
	print OUT5 "nohup ```famap -tN -b $ARGV[0]/e_pcr/$file_name.famap $ARGV[0]/Acc_genome/$file_name.fa;fahash -b $ARGV[0]/e_pcr/$file_name.hash -w 10 -f3 $ARGV[0]/e_pcr/$file_name.famap;re-PCR -S $ARGV[0]/e_pcr/$file_name.hash -n 3 -g 1 -o $ARGV[0]/e_pcr/$file_name.rePCR $ARGV[0]/query_prime_16167.txt;grep -v \"^#\" $ARGV[0]/e_pcr/$file_name.rePCR >$ARGV[0]/e_pcr/$file_name.FPCR``` &\n";
	if(($tcontrol>14)&&($tcontrol%15==0)){
                print OUT5 "sleep 30\n";
        }
        $tcontrol = $tcontrol+1;

}
close(IN);
close(OUT1);
close(OUT2);
close(OUT3);
close(OUT4);
close(OUT5);

## 3. DBimport and genotype call
my @count=("01","02","03","04","05","06","07","08","09","10","11","12");
open(OUT7,">$ARGV[0]/Code/DBimport_II.sh");
open(OUT8,">$ARGV[0]/Code/GenoCall_III.sh");

foreach my $count (@count){
	print OUT7 "gatk GenomicsDBImport --sample-name-map $ARGV[0]/Variant_call/cohort.sample_map --reader-threads 10 --genomicsdb-workspace-path $ARGV[0]/DBimport/chr$count"."_directory --intervals chr$count >$ARGV[0]/DBimport/Log/chr$count.log 2>&1\n";
	print OUT8 "nohup gatk GenotypeGVCFs -R $ARGV[0]/IRGSP-1.0/IRGSP-1.0_genome.fasta -V gendb://$ARGV[0]/DBimport/chr$count"."_directory -O $ARGV[0]/Geno_call/chr$count"."_output.vcf >$ARGV[0]/Geno_call/Log/chr$count.log 2>&1 &\n";

}

close(OUT7);
close(OUT8);

## IV. merge produced vcf files into one
# IV. picard GatherVcfs I=chr01_output.vcf  I=chr02_output.vcf I=chr03_output.vcf I=chr04_output.vcf I=chr05_output.vcf I=chr06_output.vcf I=chr07_output.vcf I=chr08_output.vcf I=chr09_output.vcf I=chr10_output.vcf I=chr11_output.vcf I=chr12_output.vcf  O=merged.vcf

## V. split merged vcf by samples, run the following code in a sh file
open(SPLIT,">split.sample.sh");
print SPLIT "for sample in `bcftools query -l $ARGV[0]/Geno_call/merged.vcf`\n";
print SPLIT "do\n";
print SPLIT "bcftools view -Ov -a -U -s \$sample -o $ARGV[0]/Acc_genome/\"\${sample}\".vcf $ARGV[0]/Geno_call/merged.vcf\n";
print SPLIT "done\n";
close(SPLIT);



