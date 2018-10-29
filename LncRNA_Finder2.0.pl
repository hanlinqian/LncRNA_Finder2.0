#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my $ncbi_blast="blastn"; # could be modified according to users' computational environment
my $ncbi_blast_makeblastdb="makeblastdb"; # could be modified according to users' computational environment
my $cpc="CPC2.py";#could be modified according to users' computational environment
my $bowtie_program="bowtie2"; #could be modified according to users' computational environment
my $bowtie_build_program="bowtie2-build"; #could be modified according to users' computational environment
#----------------------------main interface that the program starts from--------------------------------------
&main(); 
sub main(){
	my $USAGE = qq(
Welcome to use LncRNA_Finder2.0 to identify long noncoding RNAs only according to the nucleotide attributes of input sequences!

USAGE:
   perl LncRNA_Finder2.0.pl -i <transcript.fasta> -k <housekeeping.fasta> -o <output prefix> [-t <# of thread>] [-r <minimum lncRNA length>] [-f <maximum ORF length>] [-e <E-value of alignment>]
Options:
   -i <transcript.fasta>
   -k <housekeeping.fasta>
   -o <output prefix>
   -h help
   -t <int> number of thread for the computation || default=4
   -r <int> minimum lncRNA length || default=200
   -f <int> maximum potential ORF length of lncRNAs || default=100
   -e E-value of the alignment against protein database || default=1.0e-3
   
Requirement:
   Need to install BLAST+ standalone package, bowtie2 and cpc2 in your local environment correctly and set the running path of external programs at the begining of the pipeline
	);
	my ($inputfile,$hkfile,$srnafile,$profile,$outfile); #required parameters
	my ($mintslength,$maxorflength,$evalue,$thread);#optional parameters
	#
	my $parastr = &GetOptions("i=s{1}"=>\$inputfile,
							  "k=s{1}"=>\$hkfile,
							  "o=s{1}"=>\$outfile,
							  "t=i{0,1}"=>\$thread,
							  "r=i{0,1}"=>\$mintslength,
							  "f=i{0,1}"=>\$maxorflength,
							  "e=s{0,1}"=>\$evalue
	);
	GetOptions ("help|?");

	unless ($parastr && defined($outfile)) {
	   print $USAGE."\n";
	   exit 1;
	}else{
		print "----------------------------------------------------------\n\n";
		unless(-e $inputfile){
			print "input transcript sequence file-$inputfile does not exist\n";
			die("Pipeline halted\n");
		}else{
			print "input transcript sequence file: $inputfile\n";
		}		
		unless(-e $hkfile){
			print "input housekeeping sequence file-$hkfile does not exist\n";
			die("Pipeline halted\n");
		}else{
			print "input housekeeping sequence file: $hkfile\n";
		}		
		#check the optional parameters
		unless(defined($thread)){
			$thread=4;
		}
		unless(defined($mintslength)){
			$mintslength=200;
		}
		unless(defined($maxorflength)){
			$maxorflength=100;
		}
		unless(defined($evalue)){
			$evalue=1.0e-3;
		}
		print "----------------------------------------------------------\n\n";
		print "LncRNA_Finder starts to work with the following parameters:\n";
		print "number of thread for computation:\t$thread\n";
		print "minimum lncRNA length:\t$mintslength\n";
		print "maximum potential ORF length of lncRNAs:\t$maxorflength\n";
		print "----------------------------------------------------------\n\n";
	}
		
	my %oriseqhash=getOriSeqinfo($inputfile);
	print "processing transcript length/ORF filter...\n";
	print "----------------------------------------------------------\n\n";
	#transcript length/ORF filter
	my $tslorf=$outfile."-RNAlen_ORFlen_info.txt";
	my $lenfilterleftseqfile=$outfile."-lenfilterleft.fa";
	open LENOUT,">$tslorf" || die "Cannot create result file:$!";
	open LEFTOUT,">$lenfilterleftseqfile" || die "Cannot create result file:$!";
	print LENOUT "seqname\ttranscriptlen\tORF1\tORF2\tORF3\trevORF1\trevORF2\trevORF3\n";
	while (my ($tkey,$tvalue)=each(%oriseqhash)){
		my $toutstr=$tkey;
		$toutstr.="\t".length($tvalue);
		my $checkstr=checkprotein(translate_frame($tvalue, 1));
		my $checkstr2=checkprotein(translate_frame($tvalue, 2));
		my $checkstr3=checkprotein(translate_frame($tvalue, 3));
		$toutstr.="\t".length($checkstr)."\t".length($checkstr2)."\t".length($checkstr3);
		my $revcom = revcom($tvalue);
		my $checkstr4=checkprotein(translate_frame($revcom, 1));
		my $checkstr5=checkprotein(translate_frame($revcom, 2));
		my $checkstr6=checkprotein(translate_frame($revcom, 3));
		$toutstr.="\t".length($checkstr4)."\t".length($checkstr5)."\t".length($checkstr6);
		print LENOUT $toutstr."\n";
		if(length($checkstr)>$maxorflength || length($checkstr2)>$maxorflength || length($checkstr3)>$maxorflength || length($checkstr4)>$maxorflength || length($checkstr5)>$maxorflength || length($checkstr6)>$maxorflength || length($tvalue)<$mintslength){
			$oriseqhash{$tkey}="";
		}else{
			print LEFTOUT ">".$tkey."\n".$tvalue."\n";
		}
	}
	close(LEFTOUT);
	close(LENOUT);
	print "transcript length/ORF filter completed\n";
	print "----------------------------------------------------------\n\n";
	print "processing CPC filter...\n";
	print "----------------------------------------------------------\n\n";
	#cpc potential protein score calculation
	my $cpc_res_file=$outfile."-cpc_res_table.txt";
	my $cpccmd="$cpc -i $lenfilterleftseqfile -o $cpc_res_file";
	print $cpccmd."\n";
	if(system($cpccmd)!=0){
		die("CPC program has not been installed properly to be invoked by lncRNA_Finder Pipeline\nPipeline\tHalted\n");
	}
	#parse the cpc result
	open CPCIN,"$cpc_res_file" || die "Cannot open $cpc_res_file:$!";
	my $cpcstr;
	while($cpcstr=<CPCIN>){
		my @barr=split(/\t/,trim($cpcstr));
		if(trim($barr[7]) eq "coding"){
			$oriseqhash{trim($barr[0])}="";
		}
	}
	close(CPCIN);
	my $CPCleftseqfile=$outfile."-CPC_left.fa";
	open LEFTOUT,">$CPCleftseqfile" || die "Cannot create result file:$!";
	while (my ($tkey,$tvalue)=each(%oriseqhash)){
		if($tvalue ne ""){
			print LEFTOUT ">".$tkey."\n".$tvalue."\n";
		}
	}
	close(LEFTOUT);
	print "CPC filter completed\n";
	print "----------------------------------------------------------\n\n";
	print "processing housekeeping RNAs filter...\n";
	print "----------------------------------------------------------\n\n";
	#align the transcript to housekeeping RNAs
	my $HKformatdbcmd=$ncbi_blast_makeblastdb." -in $hkfile -dbtype nucl -parse_seqids -out hk_db";
	print $HKformatdbcmd."\n";
	if(system($HKformatdbcmd)!=0){
		die("NCBI BLAST package has not been installed properly to be invoked by lncRNA_Finder Pipeline\nPipeline\tHalted\n");
	}
	my $hk_blast_file=$outfile."-housekeeping_blast_res.txt";
	`$ncbi_blast -db hk_db -query $CPCleftseqfile -num_threads $thread -evalue $evalue -out $hk_blast_file -outfmt "6 qseqid sseqid qcovs"`;
	#print $HKblastcmd."\n";
	#parse the protein DB blast result
	open BLASTIN,"$hk_blast_file" || die "Cannot open $hk_blast_file:$!";
	my $hkblaststr;
	while($hkblaststr=<BLASTIN>){
		my @barr=split(/\t/,trim($hkblaststr));
		if ($barr[2] >= 98){
			$oriseqhash{trim($barr[0])}="";
		}
	}
	close(BLASTIN);
	my $HKdbleftseqfile=$outfile."-putativelncRNA.fa";
	open LEFTOUT,">$HKdbleftseqfile" || die "Cannot create result file:$!";
	while (my ($tkey,$tvalue)=each(%oriseqhash)){
		if($tvalue ne ""){
			print LEFTOUT ">".$tkey."\n".$tvalue."\n";
		}
	}
	close(LEFTOUT);
	print "housekeeping RNAs filter completed\n";
	print "----------------------------------------------------------\n\n";
	print "----------------------------------------------------------\n\n";
	print "putative lncRNA classification done!!!\n";	
	print "Congratulations! LncRNA_Finder has completed all the filtering analyses!!!\n";
	print "Please check out the following result files:\n";
	print "Putative lncRNA file:\t$HKdbleftseqfile\n";	
	print "----------------------------------------------------------\n\n";
}
sub getOriSeqinfo(){
	my $orifile=shift;
	my %resarr;
	open ORI,$orifile or die "Cannot open $orifile:$!";
	my $tstr;
	while($tstr=<ORI>){
		if(trim($tstr) ne "" && $tstr=~/>/){
			my @tarr=split(/\s+/,trim($tstr));
			my $tseqstr="";
			my $nstr="";
			while ($nstr=<ORI>){
				if($nstr=~/>/){
					my $tname=substr(trim($tarr[0]),1);
					$resarr{$tname}=$tseqstr;
					seek(ORI,-length($nstr),1);
					last;
				}else{
					$tseqstr.=trim($nstr);
				}
			}
		}
	}
	close(ORI);
	return %resarr;
}
sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }#else{

#            print STDERR "Bad codon \"$codon\"!!\n";
#            exit;
#    }
}

sub dna2peptide {

    my($dna) = @_;

    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aa( substr($dna,$i,3) );
    }

    return $protein;
}

sub translate_frame {

    my($seq, $start, $end) = @_;

    my $protein;

    # To make the subroutine easier to use, you won't need to specify
    #  the end point-it will just go to the end of the sequence
    #  by default.
    unless($end) {
        $end = length($seq);
    }

    # Finally, calculate and return the translation
    return dna2peptide ( substr ( $seq, $start - 1, $end -$start + 1) );
}

# revcom 
#
# A subroutine to compute the reverse complement of DNA sequence

sub revcom {

    my($dna) = @_;

    # First reverse the sequence
    my $revcom = reverse $dna;

    # Next, complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C
    $revcom =~ tr/ACGTacgt/TGCAtgca/;

    return $revcom;
}

sub checkprotein(){
	my $tstr=shift;
	my $resstr="";
	my @tarr=split(/\*/,trim($tstr));
	my $maxlen=0;
	my $maxstr="NA";
	my $tpro="";
	pop(@tarr);
	foreach $tpro(@tarr){
		if(trim($tpro) ne ""){
			if(index($tpro,"M")!=-1){
				if($maxlen<length(trim($tpro))-index(trim($tpro),"M")){
					$maxlen=length(trim($tpro))-index(trim($tpro),"M");
					$maxstr=substr(trim($tpro),index(trim($tpro),"M"));
				}
			}
		}
	}
	$resstr=$maxstr;
	return $resstr;
}

sub trim(){
	my $string=shift;
	$string=~s/^\s+//;
	$string=~s/\s+$//;
	return $string;
}