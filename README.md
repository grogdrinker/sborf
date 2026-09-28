# sborf

SBORF is a neural network that reads the codon usage of a coding DNA sequence (cDNA)
and estimates how strongly it will be expressed as protein. With it you can:

- **predict** the expression level of cDNA sequences (a score between 0 and 1);
- **optimize** the DNA encoding of a protein, to obtain a sequence that reaches a target expression level;
- **explain** a prediction, i.e. see how much every single codon contributes (a value between -1 and 1 per codon).

The models are organism-specific. Available organisms: `cerevisiae` (*S. cerevisiae*, the default),
`coli` (*E. coli*) and `musculus` (*M. musculus*).

## Installation

SBORF needs Python 3.9 or newer and the libraries `torch`, `numpy` and `scikit-learn`. It runs on a CPU;
a GPU is used if you ask for it (`--device cuda`).

We suggest a separate conda environment:

```sh
conda create -n sborf python=3.12
conda activate sborf
```

Install PyTorch following the instructions for your machine on https://pytorch.org/get-started/locally/
(optional: `pip` installs a default build of PyTorch if you skip this step), then install SBORF from this folder:

```sh
pip install .
```

(run it in the folder that contains this `README.md`, or give the path to it, e.g. `pip install ./sborf_bin`).
For a development install that follows your edits use `pip install -e .`.

Check that it works:

```sh
sborf --help
```

## Input

- **predict** and **explain** take **cDNA**: the coding sequence, in triplets, letters A, C, G, T (U is read as T,
  lower case is fine).
- **optimize** takes a **protein**: amino-acid letters.
- A single sequence can be given directly, or a FASTA file with any number of sequences (`>` header lines).
- A stop codon (TAA, TAG or TGA) at the end of a cDNA sequence is removed automatically.

If something is wrong with the input, the run stops with a message that names the sequence and the position,
for example: `Sequence 'x': stop codon TAA inside the sequence, at codon 10 (nucleotides 28-30).`

## Command line

```
sborf [-c {predict,optimize,explain}] SEQUENCE_OR_FILE [options]
```

`SEQUENCE_OR_FILE` is a sequence typed on the command line, or the path of a FASTA file.

| Option | Meaning | Default |
|---|---|---|
| `-c`, `--command` | `predict`, `optimize` or `explain` | `predict` |
| `--organism` | `cerevisiae`, `coli` or `musculus` | `cerevisiae` |
| `-o`, `--outfile` | output file (a folder for `explain`); printed on screen if not given | none |
| `--device` | `cpu` or `cuda` | `cuda` if available, otherwise `cpu` |
| `-t`, `--target_optimization` | optimize: `1` maximizes the expression, `0` minimizes it | `1` |
| `-n`, `--num_optimized_seqs` | optimize: number of designs returned per protein | `1` |
| `-i`, `--iterations` | optimize: iterations of the genetic algorithm | `100` |
| `--pop_size` | optimize: population size of the genetic algorithm | `1500` |
| `--mating_parents` | optimize: parents mating in each generation | `200` |
| `-v`, `--silent` | optimize: do not print the progress | off |

## A working example

The folder `examples/` contains two small inputs, all natural *S. cerevisiae* sequences:

- `example_cdna.fasta`: the cDNA of two genes, `YOR012W` and `YHR059W`;
- `example_protein.fasta`: the protein `YDR045C` (110 amino acids).

Run the commands from the folder that contains this `README.md`.

### 1. Predict

```sh
sborf examples/example_cdna.fasta --organism cerevisiae
```

```
Running prediction
###########
# SBORF RESULTS #
YOR012W 0.399
YHR059W 0.599
###################
Done in  0.093 seconds. The monkeys are listening
```

Every sequence gets a score between 0 (low expression) and 1 (high expression). Here `YHR059W` is predicted
to be expressed more than `YOR012W`. To save the scores in a tab-separated file:

```sh
sborf examples/example_cdna.fasta --organism cerevisiae -o predictions.tsv
```

`predictions.tsv` contains:

```
Name	Score
YOR012W	0.399
YHR059W	0.599
```

### 2. Explain

```sh
sborf examples/example_cdna.fasta -c explain --organism cerevisiae -o explanations
```

This writes one HTML file per sequence in the folder `explanations/` (created if needed):
`explanations/YOR012W.html` and `explanations/YHR059W.html`. Open them in a browser: each is a table
with the position, the codon and its score, colored from red (the codon lowers the predicted expression) to
green (it raises it). The scores are between -1 and 1, and 0 (white) is the median contribution of a codon
in that organism, so colors can be compared between sequences. The explain feature is still in development.

### 3. Optimize

```sh
sborf examples/example_protein.fasta -c optimize --organism cerevisiae -n 3 -i 20 --pop_size 200 --mating_parents 50 -o optimized.fasta
```

This designs three cDNA sequences that encode `YDR045C` and are predicted to be highly expressed
(the natural DNA of this protein scores 0.243). It takes a few seconds here because the genetic algorithm is
run with reduced settings (20 iterations, population of 200); the defaults (100 iterations, population of
1500) explore more and take longer. `optimized.fasta` is a normal FASTA file:

```
>YDR045C_0
ATGTTATCATTTTGCCCATCGTGCAACAACATGTTATTAATCACTTCGGGTGACTCCGGTGTTTACACATTAGCATGTCGTTCC...
>YDR045C_1
...
>YDR045C_2
...
```

The designs are sorted from the best to the worst. The optimization is stochastic, so the sequences
are different at every run. To check them, predict their expression like in step 1:

```sh
sborf optimized.fasta --organism cerevisiae
```

## Using SBORF from Python

The same steps in a Python script (run it from the folder that contains this `README.md`):

```python
from sborf.src import parse
from sborf.sborf_standalone import modify_sequences, quality_check_DNASeq
from sborf.predict import run_prediction
from sborf.explain import run_explaination
from sborf.optimize import optimize

# predict: {sequence id: cDNA} -> {sequence id: score between 0 and 1}
sequences = modify_sequences(parse.leggifasta("examples/example_cdna.fasta"))
quality_check_DNASeq(sequences)   # stops with an explanatory message if a sequence is not valid cDNA
scores = run_prediction(sequences, organism="cerevisiae")
for name, score in scores.items():
    print(name, round(score, 3))

# explain: {sequence id: one value between -1 and 1 for every codon}
per_codon = run_explaination(sequences, organism="cerevisiae")
for name, values in per_codon.items():
    print(name, len(values), "codons, min %.2f, max %.2f" % (min(values), max(values)))

# optimize: protein -> list of cDNA sequences, best first
protein = parse.leggifasta("examples/example_protein.fasta")["YDR045C"]
designs = optimize(protein, organism="cerevisiae", iterations=20, pop_size=200,
                   num_parents_mating=50, num_optimized_seqs=2)
scores = run_prediction(dict(enumerate(designs)), organism="cerevisiae")
for i, dna in enumerate(designs):
    print("design", i, len(dna), "nucleotides, predicted score", round(scores[i], 3))
```

Its output (the designs and their scores change at every run):

```
YOR012W 0.399
YHR059W 0.599
YOR012W 137 codons, min -1.00, max 0.34
YHR059W 130 codons, min -1.00, max 1.00
design 0 330 nucleotides, predicted score 0.889
design 1 330 nucleotides, predicted score 0.876
```

`modify_sequences` prepares the sequences (upper case, removes a stop codon at the end) and
`quality_check_DNASeq` validates them. `run_prediction` and `run_explaination` expect sequences prepared
in this way: cDNA in upper case, a length that is a multiple of 3, and no stop codon.

## Contents of the package

```
sborf/
  sborf_standalone.py   the `sborf` command
  predict.py            run_prediction
  explain.py            run_explaination
  optimize.py           optimize
  src/                  neural network, parsing and output code
  models/               trained models of the three organisms
```

## License

Apache License 2.0, see `LICENSE`.

## Help

For bug reports, feature requests and technical questions please contact gabriele.orlando@umontpellier.fr
