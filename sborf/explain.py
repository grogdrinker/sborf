from sborf.src import NNWrappers
from sborf.src import parse
from sborf.src.writeExplaination import visualize_codons

def run_explaination(seqs, organism, printPreds=False,device="cpu"):
    wrapper = NNWrappers.NNwrapper(device=device)
    wrapper.load(organism=organism)
    if type(seqs) is list:

        seqs = {i:seqs[i] for i in range(len(seqs))}

    elif type(seqs) is str:
        seqs = parse.leggifasta(seqs)

    assert type(seqs) is dict

    names = []
    x = []
    for prot in seqs.keys():
        names += [prot]
        x += [[seqs[prot][i:i + 3] for i in range(0, len(seqs[prot]), 3)]]

    preds = wrapper.explain(x)
    diz={}
    for k in range(len(names)):
        diz[names[k]] = preds[k]

    if printPreds:
        print("###########")
        print("# RESULTS #")
        for k in diz.keys():
            visualize_codons([seqs[k][i:i+3] for i in range(0, len(seqs[k]), 3)], diz[k])
        print("###########")

    return diz