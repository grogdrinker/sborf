#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  
#  Copyright 2019 Gabriele Orlando <orlando.gabriele89@gmail.com>
#  

from sborf.src import parse,utils
from sborf.optimize import optimize
import os,time,sys
from sborf.predict import run_prediction
from sborf.explain import run_explaination
from sborf.src.writeExplaination import visualize_codons,write_codons_visualization

def parse_input(sequence):
    """Resolve the ``sequence`` argument into a {name: seq} dict.

    Accepts, in order of priority:
      * a path to a fasta file on disk;
      * fasta content copy-pasted directly on the command line (contains '>');
      * a single raw sequence (any whitespace/newlines are stripped).
    """
    if os.path.exists(sequence):
        sequences = parse.leggifasta(sequence)
    elif ">" in sequence:
        sequences = parse.leggifasta_string(sequence)
    else:
        sequences = {"inputSequence": "".join(sequence.split())}
    if not sequences:
        print("### ERROR! ###")
        print("No sequence was found in the input. Paste one sequence, or a FASTA in which a line starting with '>' precedes each sequence.")
        sys.exit(0)
    return sequences

STOP_CODONS = ("TAA", "TAG", "TGA")

def modify_sequences(dna):
    """Normalise cDNA sequences: no whitespace, upper case, U read as T, and the stop codon
    removed if the sequence ends with one (each sequence of a FASTA is handled on its own)."""
    for i in dna.keys():
        s = "".join(dna[i].split()).upper().replace("U", "T")
        if len(s) % 3 == 0 and s[-3:] in STOP_CODONS:
            s = s[:-3]
        dna[i] = s
    return dna

def find_dna_problems(sequences, per_sequence_limit=200):
    """Every problem of sequences already normalised by modify_sequences, as
    (sequence id, problem, position) tuples; the position is empty for whole-sequence problems."""
    allowedLetters = parse.allowed_codons
    problems = []
    for k in sequences:
        s = sequences[k]
        if len(s) % 3 != 0:
            problems.append((k, "the length is "+str(len(s))+" nucleotides, which is not a multiple of 3 (cDNA is made of triplets)", ""))
            continue
        codons = [s[i:i + 3] for i in range(0, len(s), 3)]
        if not codons:
            problems.append((k, "no codons (the sequence is empty, or it is only a stop codon)", ""))
            continue
        found = 0
        for position, codon in enumerate(codons, 1):
            if codon in STOP_CODONS:
                problem = "stop codon "+codon+" inside the sequence"
            elif codon not in allowedLetters:
                problem = "invalid triplet '"+codon+"' (a cDNA can only contain the letters A, C, G and T)"
            else:
                continue
            found += 1
            if found <= per_sequence_limit:
                problems.append((k, problem, "codon "+str(position)+" (nucleotides "+str(3*position-2)+"-"+str(3*position)+")"))
        if found > per_sequence_limit:
            problems.append((k, "and "+str(found-per_sequence_limit)+" more problems of the same kind (only the first "+str(per_sequence_limit)+" are listed)", ""))
    return problems

def quality_check_DNASeq(sequences, error_report=None):
    """Check sequences already normalised by modify_sequences and stop the run if anything is wrong.
    One problem: a specific message. Several: a summary message, with the full list (sequence id,
    problem, position) written to the tab-separated file error_report if a path is given, otherwise printed."""
    problems = find_dna_problems(sequences)
    if not problems:
        return
    print("### ERROR! ###")
    if len(problems) == 1:
        k, problem, where = problems[0]
        print("Sequence '"+k+"': "+problem+(", at "+where if where else "")+".")
    elif error_report is not None:
        with open(error_report, "w") as f:
            f.write("Sequence\tProblem\tPosition\n")
            for k, problem, where in problems:
                f.write(k+"\t"+problem+"\t"+where+"\n")
        names = sorted({p[0] for p in problems})
        print("Multiple errors in "+("multiple sequences" if len(names) > 1 else "sequence '"+names[0]+"'")+": look at the log report for more information.")
    else:
        for k, problem, where in problems[:50]:
            print("Sequence '"+k+"': "+problem+(", at "+where if where else "")+".")
        if len(problems) > 50:
            print("... and "+str(len(problems)-50)+" more problems.")
    sys.exit(0)

def quality_check_aaSeq(sequences):
    allowedLetters = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
    for k in sequences:
        upperseq = sequences[k].upper()
        if len(upperseq) >= 30 and set(upperseq) <= set("ACGTU"):
            print("### ERROR! ###")
            print("The input "+k+" looks like a DNA sequence (it only contains the letters A, C, G, T). Optimize expects an amino-acid sequence: use predict or explain for DNA, or translate the sequence into a protein first.")
            sys.exit(0)
        for l in upperseq:
            if not l in allowedLetters:
                print("### ERROR! ###")
                print("Non-amino acid letter "+l+" found in sequence"+k)
                sys.exit(0)

def run_optimization(sequence, organism, outfile,verbose=True,iterations=100,mating_parents=200,pop_size=1500,target_sol=1.0,num_optimized_seqs=1,device="cpu",cancel_check=None):
    print("Running optimization")
    old_time = time.time()
    sequences = parse_input(sequence)
    sequences = {k: v.upper() for k, v in sequences.items()}

    quality_check_aaSeq(sequences)

    final_diz= {}
    for sname in sequences.keys():

        s = optimize(target_seq = sequences[sname], organism=organism, verbose=verbose,iterations=iterations,num_parents_mating=mating_parents,pop_size=pop_size,TARGET_SOL=target_sol,num_optimized_seqs=num_optimized_seqs,device=device,cancel_check=cancel_check)
        if num_optimized_seqs!=1:
            for k in range(len(s)):
                final_diz[sname+"_"+str(k)] = s[k]
        else:
            final_diz[sname] = s[0]

    if outfile is None:
        print("###########")
        print("# SBORF RESULTS #")
        for k in final_diz.keys():
            print(">"+k+"\n"+ final_diz[k]+"\n")
        print("###################")
    else:
        f=open(outfile,"w")
        for k in final_diz.keys():
            f.write(">"+k+"\n"+ final_diz[k]+"\n")
        f.close()

    print("Done in ",round(time.time()-old_time,4),"seconds. The monkeys are listening")

def prediction(sequence,organism,outfile=None,device="cpu",error_report=None):
    print("Running prediction")
    old_time = time.time()
    sequences = parse_input(sequence)
    sequences = modify_sequences(sequences)
    quality_check_DNASeq(sequences, error_report)
    pred = run_prediction(sequences,organism,printPreds=False,device=device)
    if outfile is None:
        print("###########")
        print("# SBORF RESULTS #")
        for k in pred.keys():
            print(k, round(pred[k],3))
        print("###################")
    else:
        f=open(outfile,"w")
        f.write("Name\tScore\n")
        for k in pred.keys():
            f.write(k+ "\t"+str(round(pred[k],3))+"\n")
        f.close()
    print("Done in ",round(time.time()-old_time,4),"seconds. The monkeys are listening")

def explaination(sequence,organism,outfolder=None,device="cpu",error_report=None):
    print("Running explaination")
    print("### DISCLAIMER ### the explain feature is still in development. Results may change and should be interpreted with caution.")
    old_time = time.time()
    sequences = parse_input(sequence)
    sequences = modify_sequences(sequences)
    quality_check_DNASeq(sequences, error_report)

    pred = run_explaination(sequences,organism,device=device)

    if outfolder is None:
        print("###########")
        print("# SBORF RESULTS #")
        for k in pred.keys():
            print(">"+k)
            visualize_codons([sequences[k][i:i + 3] for i in range(0, len(sequences[k]), 3)], pred[k])
        print("###################")
    else:
        os.makedirs(outfolder, exist_ok=True)
        for k in pred.keys():
            outfile = outfolder+"/"+k+".html"
            write_codons_visualization([sequences[k][i:i + 3] for i in range(0, len(sequences[k]), 3)], pred[k],outfile=outfile)
    print("Done in ",round(time.time()-old_time,4),"seconds. The monkeys are listening")

def main():
    import argparse,sys,torch
    import textwrap
    parser = argparse.ArgumentParser(
        prog='Sborf',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
             if you have problems or you bugs,
             mail orlando.gabriele89@gmail.com.
             
             The monkeys are listening
             '''))


    parser.add_argument('--command',"-c", default="predict",help="The action to perform. predict evaluates with sborf cDNA sequence(s), optimize runs an optimization of amino acid sequence(s) to obtain an encoding with interaction probability defined by target_optimization",choices=['predict', "optimize","explain"])
    parser.add_argument('sequence', help='If predict, the input is either a cDNA sequence or a fasta file with cDNA sequences, if optimize, the input is either a amino acid sequence or a fasta file with amino acid sequences')
    parser.add_argument('--iterations',"-i", help='number of iterations for the genetic optimization. Ignored if command is predict',default=100,type=int)
    parser.add_argument('--mating_parents', help='number of parents mating for the genetic optimization. Ignored if command is predict',default=200,type=int)
    parser.add_argument('--pop_size', help='size of a generation for the genetic optimization. Ignored if command is predict',default=1500,type=int )
    parser.add_argument('--target_optimization',"-t", help='target fitness for the genetic optimization. 1 means the sequence is optimized TO MAXIMIZE total expression, 0 means the sequence is optimized TO MINIMIZE PROTEIN abundance. REMEMBER MAXIMISING PROTEIN ABUNDANCE DOES NOT ENSURE MAXIMAL FUNCTIONAL PROTEIN. Ignored if command is predict',default=1 ,type=int)
    parser.add_argument('--num_optimized_seqs',"-n", help='number of output optimized sequences per input sequence. Ignored if command is predict',default=1,type=int )
    parser.add_argument('--outfile',"-o", help='the output file. if not provided, it prints on screen. If command is explain, it a folder is expected',default=None )
    parser.add_argument('--silent',"-v", action='store_true',help='does not print text while optimizing. ### CAREFULL ### you will not be able to check if the optimization converges or not')
    parser.add_argument("--organism", choices=["coli", "cerevisiae", "musculus"], default="cerevisiae",help="select the organism used for the optimization" )
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu",help="device where to run sborf" )

    args = parser.parse_args()

    if args.command == "predict":
        prediction(args.sequence,outfile=args.outfile,organism=args.organism,device=args.device)
    elif args.command == "optimize":
        run_optimization(args.sequence,outfile=args.outfile,verbose = not args.silent,iterations=args.iterations,mating_parents=args.mating_parents,pop_size=args.pop_size,target_sol=args.target_optimization,num_optimized_seqs=args.num_optimized_seqs,organism=args.organism,device=args.device)
    elif args.command == "explain":
        explaination(args.sequence,outfolder=args.outfile,organism=args.organism,device=args.device)

if __name__ == '__main__':

    main()


