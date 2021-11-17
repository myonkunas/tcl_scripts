# To render the protein with hydrophobicity in b column
#Amino acid scale: Hydrophobicity scale based on free energy of transfer (kcal/mole).
#Author(s): Guy H.R.
#Reference: Biophys J. 47:61-70(1985).

#array set bvalue {
#ALA  0.100  
#ARG  1.910  
#ASN  0.480  
#ASP  0.780  
#CYS -1.420  
#GLN  0.950  
#GLU  0.830  
#GLY  0.330  
#HIS -0.500  
#ILE -1.130  
#LEU -1.180  
#LYS  1.400  
#MET -1.590  
#PHE -2.120  
#PRO  0.730  
#SER  0.520  
#THR  0.070  
#TRP -0.510  
#TYR -0.210  
#VAL -1.270  
#
#}

#Scaled hydrophobicity from 
#http://psyche.uthct.edu/shaun/SBlack/aagrease.html
#from Dr. Shaun D. Black (University of Texas Health Center at Tyler)
array set bvalue {
ALA 0.616
CYS 0.680
ASP 0.028
GLU 0.043
PHE 1.00
GLY 0.501
HIS 0.165
ILE 0.943
LYS 0.283
LEU 0.943
MET 0.738
ASN 0.236
PRO 0.711
GLN 0.251
ARG 0.000
SER 0.359
THE 0.450
VAL 0.825
TRP 0.878
TYR 0.880

}
# Suppose change the current top molecule
set all [atomselect top all]
$all set beta 0

foreach residuename [array names bvalue ] {
	set res [atomselect top "resname $residuename"]
	$res set beta $bvalue($residuename)
	}



#Hydrophobicity Scale Values
#Amino Acid;  Engleman-Steitz; Hopp-Woods; Kyte-Doolittle; Janin; Chothia; Eisenberg-Weiss
#PHE -3.7 -2.5 2.8 0.5 0.0 0.61
#MET -3.4 -1.3 1.9 0.4 -0.24 0.26
#ILE -3.1 -1.8 4.5 0.7 0.24 0.73
#LEU -2.8 -1.8 3.8 0.5 -0.12 0.53
#VAL -2.6 -1.5 4.2 0.6 0.09 0.54
#CYS -2.0 -1.0 2.5 0.9 0.0 0.04
#TRP -1.9 -3.4 -0.9 0.3 -0.59 0.37
#ALA -1.6 -0.5 1.8 0.3 -0.29 0.25
#THR -1.2 -0.4 -0.7 -0.2 -0.71 -0.18
#GLY -1.0 0.0 -0.4 0.3 -0.34 0.16
#SER -0.6 0.3 -0.8 -0.1 -0.75 -0.26
#PRO 0.2 0.0 -1.6 -0.3 -0.9 -0.07
#TYR 0.7 -2.3 -1.3 -0.4 -1.02 0.02
#HIS 3.0 -0.5 -3.2 -0.1 -9.94 -0.40
#GLN 4.1 0.2 -3.5 -0.7 -1.53 -0.69
#ASN 4.8 0.2 -3.5 -0.5 -1.18 -0.64
#GLU 8.2 3.0 -3.5 -0.7 -0.90 -0.62
#LYS 8.8 3.0 -3.9 -1.8 -2.05 -1.1
#ASP 9.2 3.0 -3.5 -0.6 -1.02 -0.72
#ARG 12.3 3.0 -4.5 -1.4 -2.71 -1.8
#THRESHOLD VALUES
#Hydrophobic -1.4 -0.75 0.70 0.10 -0.47 0.10
#Hydrophilc 1.85 1.65 -2.4 -0.45 -0.98 -0.51


#Amino acid scale: Hydropathicity.
#Author(s): Kyte J., Doolittle R.F.
#
#Reference: J. Mol. Biol. 157:105-132(1982).
#
#Amino acid scale values:
#
#Ala:  1.800  
#Arg: -4.500  
#Asn: -3.500  
#Asp: -3.500  
#Cys:  2.500  
#Gln: -3.500  
#Glu: -3.500  
#Gly: -0.400  
#His: -3.200  
#Ile:  4.500  
#Leu:  3.800  
#Lys: -3.900  
#Met:  1.900  
#Phe:  2.800  
#Pro: -1.600  
#Ser: -0.800  
#Thr: -0.700  
#Trp: -0.900  
#Tyr: -1.300  
#Val:  4.200  


#Scaled Hydrophobicity of the Physiological L-alpha-Amino Acids
#from Dr. Shaun D. Black (University of Texas Health Center at Tyler)
#Parameters for the Unmodified Physiological L-alpha-Amino Acids
#Amino Acid 3-Letter Code 1-Letter Code Molecular Weight? Hydrophobicity?
#Alanine Ala A 89.09 0.616
#Cysteine Cys C 121.16 0.680
#Aspartate Asp D 133.10 0.028
#Glutamate Glu E 147.13 0.043
#Phenylalanine Phe F 165.19 1.00
#Glycine Gly G 75.07 0.501
#Histidine His H 155.16 0.165
#Isoleucine Ile I 131.18 0.943
#Lysine Lys K 146.19 0.283
#Leucine Leu L 131.18 0.943
#Methionine Met M 149.21 0.738
#Asparagine Asn N 132.12 0.236
#Proline Pro P 115.13 0.711
#Glutamine Gln Q 146.15 0.251
#Arginine Arg R 174.20 0.000
#Serine Ser S 105.09 0.359
#Threonine The T 119.12 0.450
#Valine Val V 117.15 0.825
#Tryptophan Trp W 204.23 0.878
#Tyrosine Tyr Y 181.19 0.880
#
#  ?   The molecular weights given are those of the neutral, free amino acids; residue weights can be obtained by subtraction of one equivalent of water (18 g/mol).
#
#  ?   The hydrophobicities given are the "Scaled" values from computational log(P) determinations by the "Small Fragment Approach" (see, "Development of Hydrophobicity Parameters to Analyze Proteins Which Bear Post- or Cotranslational Modifications" Black, S.D. and Mould, D.R. (1991) Anal. Biochem. 193, 72-82). The equation used to scale raw log(P) values to the scaled values given is as follows: Scaled Parameters = (Raw Parameters + 2.061)/4.484 .
#
#Trend of Hydrophobicity Parameters for the Physiological L-alpha-Amino Acids
#Phe > Leu = Ile > Tyr = Trp > Val > Met > Pro > Cys > Ala > Gly >
#Thr > Ser > Lys > Gln > Asn > His > Glu > Asp > Arg
