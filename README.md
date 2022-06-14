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

Rasa Giniunaite and Daiva Petkeviciute-Gerlach.
Predicting the sequence-dependent configuration and
energy of DNA in a nucleosome by coarse-grain modelling.
Submitted (2022).

For the cgDNA+ model, please cite:

Patelli, A. S., 2019. A sequence-dependent coarse-grain model of B-DNA with explicit description of bases and phosphate
groups parametrised from large scale Molecular Dynamics simulations. Ph.D. thesis, EPFL.

Petkeviciute, D., M. Pasi, O. Gonzalez, and J. Maddocks, 2014. cgDNA: a software package for the prediction of
sequence-dependent coarse-grain free energies of B-form DNA. Nucleic Acids Res. 42:e153–e153.
