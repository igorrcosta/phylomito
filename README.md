phylomito
===========

A python script for mitochondrial supermatrix phylogenomics.
```
usage: phylomito.py [-h] -i [INPATH] [-o [OUTPATH]]
                    [-e [EXTENSION [EXTENSION ...]]] [-b [BOOTSTRAP]]
                    [-p [PROTEIN]] [-g [GENE_TREE]] [-d [DLOOP]]

Phylomito is a simple pipeline to automatize mitochondrial super-matrix
phylogenomic, using clustaw, phyML and mrbayes.

optional arguments:
  -h, --help            show this help message and exit
  -i [INPATH], --inpath [INPATH]
                        Path to the folder with genbank sequences. (default:
                        None)
  -o [OUTPATH], --outpath [OUTPATH]
                        Path were the alignments and phylogenetic tree will be
                        saved. (default: /home/igor/phylomito/)
  -e [EXTENSION [EXTENSION ...]], --extension [EXTENSION [EXTENSION ...]]
                        Extension for the genbank files. (default: ['.gbk',
                        '.gb'])
  -b [BOOTSTRAP], --bootstrap [BOOTSTRAP]
                        Number of bootstrap repetitions on PhyML. (default:
                        101)
  -p [PROTEIN], --protein [PROTEIN]
                        Set this flag for protein sequences alignment and
                        phylogeny. (default: False)
  -g [GENE_TREE], --gene_tree [GENE_TREE]
                        Set this flag if you want to make a tree for every
                        gene. (default: False)
  -d [DLOOP], --dloop [DLOOP]
                        Flag to include DLOOP region in the alignment.
                        (default: False)
```

Quickstart:
===========


* Save your genebank files in a folder (for example, ./genebank/) and create a folder for the output (./output/). 
* Make sure your genebank files have the extension '.gbk' or '.gb'.
* Run the command:
```
python phylomito.py -i ./genebank/ -o ./outpath/
```

The most commom problem during your run is a genebank file with bad format. The error will look like this:

```
mitochondria1.gb
mitochondria2.gb
mitochondria3.gb
COX_1 is not a known gene. Replace the CDS gene id with one of the following:
ND1  ND2  COX1  COX2  ATP8  ATP6  ND3  ND4L  ND4  ND5  CYTB  ND6  COX3
Traceback (most recent call last):
  File "/home/igor/phylomito/phylomito.py", line 315, in <module>
    main(args)
  File "/home/igor/phylomito/phylomito.py", line 70, in main
    split_seqs(inpath, outpath, protein, extension, dloop)
  File "/home/igor/phylomito/phylomito.py", line 140, in split_seqs
    gene_key = gene_dict[header]
KeyError: 'COX_1'
```

Where *mitochondria3.gb* is the file where the error was found.

