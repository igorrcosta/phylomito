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
                        100)
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

##Quickstart:

* Save your genebank files in a folder (for example, ./genebank/) and create a folder for the output (./output/). 
* Make sure your genebank files have the extension '.gbk' or '.gb'.
* Run the command:
```
python phylomito.py -i ./genebank/ -o ./outpath/
```
* Your results will be in the ./outpath/ folder. The final tree file will be named `all_nuc.phy_phyml_tree.txt` by default.

##Requisites:

You need to install [PhyML](http://www.atgc-montpellier.fr/phyml/binaries.php), [CLUSTALW](http://www.clustal.org/download/current/), [python 2.7](https://www.python.org/downloads/release/python-2710/) and the [Biopython library](http://biopython.org/wiki/Download) to run this program.

This program was tested on a linux machine.

##How does it work:

This program finds all genebank files (mitogenomes) in a folder and saves, in a multifasta file, each gene that is present in all mitogenomes. These files are aligned with CLUSTALW and the alignment is concatenated in a single file (`all_aa.aln`, by default) that contains all aligned genes from all mitogenomes. The file is  Phyml uses this file to generate a Maximum Likelihood tree.
 
##Advanced features:

* You can generate an amino acid alignment and phylogeny using the -p (or --protein) flag. The default is nucleotidic alignment and phylogeny.
* Running the program with the -g (or --gene_tree) flag will generate a tree for every gene. Using this along with -d (or --dloop) to generate a tree of the DLOOP region. The supermatrix tree will also include this region. Do not use the -d flag with the -p flag, as it will translate the DLOOP region, generating giberish.
* Default number of bootstap resamples is 100. You can change this with the -b (or --bootstrap) flag. Changing this will affect how long it takes to run the phylogeny.

##Commom errors:

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

Where `mitochondria3.gb` is the file where the error was found.
In this case the solution is to edit your file and replace all occurences of `COX_1` with `COX1`.

