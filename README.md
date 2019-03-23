phylomito
===========

A python script for mitochondrial supermatrix phylogenomics.

```
usage: phylomito.py [-h] -i [INPATH] [-o [OUTPATH]]
                    [-e [EXTENSION [EXTENSION ...]]] [-b [BOOTSTRAP]]
                    [-p [PROTEIN]] [-g [GENE_TREE]] [-d [DLOOP]]

```

This program is licensed under GPLv3.

## Quick start:

* Save your genebank files in a folder (for example, ./genebank/) and create a folder for the output (./output/). 
* Make sure your genebank files have the extension '.gbk' or '.gb'.
* Run the command:
```
python phylomito.py -i ./genebank/ -o ./outpath/
```
* Your results will be in the ./outpath/ folder. The final tree file will be named `all_nuc.phy_phyml_tree.txt` by default.

## Requisites:

You need to install [PhyML](http://www.atgc-montpellier.fr/phyml/binaries.php), [Clustal Omega](http://www.clustal.org/omega/), [python 3](https://www.python.org/downloads/) and the [Biopython library](http://biopython.org/wiki/Download) to run this program.

This program was tested on a Linux machine.

## How does it work:

This program finds all genebank files (mitogenomes) in a folder and saves, in a multifasta file, each gene that is present in all mitogenomes. 
These files are aligned with CLUSTALW and the alignment is concatenated in a single file (`all_nuc.aln`, by default) that contains all aligned genes from all mitogenomes. 
Phyml uses this file to generate a Maximum Likelihood tree.
 
## Advanced features:

* You can generate an amino acid alignment and phylogeny using the -p (or --protein) flag. The default is nucleotidic alignment and phylogeny.
* Running the program with the -g (or --gene_tree) flag will generate a tree for every gene. 
* Using -g along with -d (or --dloop) will generate a tree of the DLOOP region. 
The supermatrix tree will also include this region. Do not use the -d flag with the -p flag, as it will translate the DLOOP region, generating gibberish.
* Default number of bootstrap resamples is 100. You can change this with the -b (or --bootstrap) flag. 
Changing this will affect how long it takes to run the phylogeny.

## Common errors:

The most common problem is bad formatted genebank files. The error will look like this:

```
mitochondria1.gb
mitochondria2.gb
mitochondria3.gb
COX_1 is not a known gene. Replace the CDS gene id with one of the following:
ND1  ND2  COX1  COX2  ATP8  ATP6  ND3  ND4L  ND4  ND5  CYTB  ND6  COX3
Traceback (most recent call last):
  File "/path/phylomito/phylomito.py", line 315, in <module>
    main(args)
  File "/path/phylomito/phylomito.py", line 70, in main
    split_seqs(inpath, outpath, protein, extension, dloop)
  File "/path/phylomito/phylomito.py", line 140, in split_seqs
    gene_key = gene_dict[header]
KeyError: 'COX_1'
```

Where `mitochondria3.gb` is the file where the error was found.
In this case the solution is to edit your file and replace all occurrences of `COX_1` with `COX1`.

