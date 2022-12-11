# nucleosome_dna

The Matlab function seq2nucleosome_grad_descent.m finds  
an optimal DNA configuration on a nucleosome
for a given sequence S, by minimising the cgDNA+ model energy
subject to constraints on the phosphate positions.

Input:

S    - any chosen DNA sequence of length 147.
       For example: S = randseq(147);

Output:

wopt - optimal configuration on a nucleosome for the sequence S;
U    - potential (free) energy, needed to deform DNA sequence S
       into the configuration wopt. 

Example of usage:

S = randseq(147);
[U, wopt] = seq2nucleosome_grad_descent(S);

If you find this code useful, please cite:

Giniūnaitė, R. and Petkevičiūtė-Gerlach, D. (2022) Predicting the configuration and energy of 
DNA in a nucleosome by coarse-grain modelling. Phys. Chem. Chem. Phys., 2022, 26124–26133.

--------------------------------------------------------------------

The cgDNA+ (cgNA+) model, also used inside the function seq2nucleosome_grad_descent.m, constructs an equilibrium 
shape coordinate vector w and a stiffness matrix K for a double stranded DNA of a given sequence S
and a model parameter set params. 

The potential (free) energy U for any configuration x of a sequence S then reads:
U = (w-x)' * K * (w-x) (here w and x are column vectors and ' is a vector transpose). 

Input:

S      - any chosen DNA sequence;
          For example: S = 'ATGAAGC';
params - a model parameter set;    
          For example: params = load('/Users/daiva/dna/mc/Di-hemi-meth-hmeth.mat');

Output:

w    - equilibrium configuration coordinate vector for the sequence S;
K    - stiffness matrix for the sequence S. 

Example of usage:

params = load('/Users/daiva/dna/mc/Di-hemi-meth-hmeth.mat');
S = 'ATGAAGC';
[w, K] = constructSeqParms(S, params);

The sequence S can include a methylated and a hydroxymethylated cytosine in a CpG dinucleotide step.
A symmetrically methylated CpG dinucleotide (when cytosines are methylated on both strands) is referred to as MN, e.g, S = ATTMNAC;
a symmetrically hydroxymethylated CpG dinucleotide: HK, e.g, S = ATTHKAC;
an asymetrically methylated CpG (methylation only on one strand): MG or CN;
an asymetrically hydroxymethylated CpG (methylation only on one strand): HG or CK.


For the cgDNA+ (cgNA+) model, please cite:

Sharma, R. cgNA+: A sequence-dependent coarse-grain model of double-stranded nucleic acids. PhD thesis EPFL (2022) # 9792.

Patelli, A. S. A sequence-dependent coarse-grain model of B-DNA with explicit description of bases and phosphate 
groups parametrised from large scale Molecular Dynamics simulations. PhD thesis EPFL (2019) # 9552.

Petkevičiūtė, D., Pasi, M., Gonzalez, O., and Maddocks, J. (2014) cgDNA: a software package for the prediction of 
sequence-dependent coarse-grain free energies of B-form DNA. Nucleic Acids Res., 42(20), e153–e153.
