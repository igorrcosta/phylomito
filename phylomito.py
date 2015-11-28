#!/usr/bin/env python
# -*- coding: utf-8 -*-
# phylomito.py

'''Mitochondrial phylogeny using the supermatrix method.'''

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

import os
import shlex
import argparse
from pprint import pprint
from subprocess import Popen
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna


genes = [('ND1', 'NAD1'), ('ND2', 'NAD2'), ('COX1', 'CO1'), ('COX2', 'CO2'), ('ATP8', 'ATPase 8'), ('ATP6', 'ATPase 6'), ('ND3', 'NAD3'), ('ND4L', 'NAD4L'), ('ND4', 'NAD4'), ('ND5', 'NAD5'), ('CYTB', 'Cyt B', 'COB'), ('ND6', 'NAD6'), ('COX3', 'CO3'), ('DLOOP',)]
known_genes = ['ND1', 'ND2', 'COX1', 'COX2', 'ATP8', 'ATP6', 'ND3', 'ND4L', 'ND4', 'ND5', 'CYTB', 'ND6', 'COX3', 'DLOOP']

gene_dict = {}
for n, gs in enumerate(genes):
    for g in gs:
        gene_dict[g] = known_genes[n]

def argument_parser(hlp = False):
    '''phylomito.py -i /path/to/genbank/ -p -o /path/to/output/
    Output default: current working directory.'''

    default_out = os.getcwd() + '/'
    parser = argparse.ArgumentParser(description = 'Phylomito is a simple pipeline to automatize mitochondrial super-matrix phylogenomic, using clustaw, phyML and mrbayes.',\
                                     argument_default = None, fromfile_prefix_chars = '@')
    parser.add_argument('-i', '--inpath', nargs = '?', type = str, required = True,\
                        dest = 'inpath', help = 'Path to the folder with genbank sequences. (default: %(default)s)')
    parser.add_argument('-o', '--outpath', nargs = '?', type = str, default = default_out,\
                        dest = 'outpath', help = 'Path were the alignments and phylogenetic tree will be saved. (default: %(default)s)')
    parser.add_argument('-l', '--log', nargs = '?', type = str, default = 'phyml.log',\
                        dest = 'log', help = 'Log file. (default: %(default)s)')
    parser.add_argument('-e', '--extension', nargs = '*', type = str, default = ['.gbk', '.gb'],\
                        dest = 'extension', help = 'Extension for the genbank files. (default: %(default)s)')
    parser.add_argument('-b', '--bootstrap', nargs = '?', type = int, default = 100 ,\
                        dest = 'bootstrap', help = 'Number of bootstrap repetitions on PhyML. (default: %(default)s)')
    parser.add_argument('-p', '--protein', nargs = '?', const = True, default = False,\
                        dest = 'protein', help = 'Set this flag for protein sequences alignment and phylogeny. (default: %(default)s)')
    parser.add_argument('-g', '--gene_tree', nargs = '?', const = True, default = False,\
                        dest = 'gene_tree', help = 'Set this flag if you want to make a tree for every gene. (default: %(default)s)')
    parser.add_argument('-d', '--dloop', nargs = '?', const = True, default = False,\
                        dest = 'dloop', help = 'Flag to include DLOOP region in the alignment. (default: %(default)s)')
    parser.add_argument('-t', '--code_table', nargs = '?', type = int, default = 2,\
                        dest = 'code_table', help = 'Genetic table to be used for translation. Only matters if --protein flag is used. Default is Vertebrate Mitochondrial Code. See all tables at http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. (default: %(default)s)')
    if hlp:
        args = parser.parse_args(['-h'])
    else:
        args = parser.parse_args().__dict__
    return args

def main(args):
    #Args processing.
    protein = args['protein']
    inpath = args['inpath']
    code_table = args['code_table']
    log_file = args['log']
    skip_phyml = False
    if not inpath.endswith('/'):
        inpath += '/'
    outpath = args['outpath']
    if not outpath.endswith('/'):
        outpath += '/'
    extension = args['extension']
    dloop = args['dloop']
    bootstrap = str(args['bootstrap'])
    if code_table not in [2,3,4,5,9,13,14,16,21,22,23,24]: #All mitochondrial codes
        print 'WARNING: Genetic code n.{} is not Mitochondial!'.format('code_table')
    if not dloop:
        known_genes.remove('DLOOP')
    try:
        a = open(outpath + log_file, 'w')
        a.close()
    except:
        print 'Was not able to open {}. Check your permissions.'.format(outpath+log_file)
        raise
    #seqfile will be the file with the concatenated alignment.
    if protein:
        seqfile = outpath + 'all_aa'
    	command = 'phyml -d aa -m JTT -b ' + bootstrap + ' -v 0.0 -c 4 -a 4 -f m -i ' + outpath + 'all_aa.phy'
    else:
        seqfile = outpath + 'all_nuc'
    	command = 'phyml -m GTR -b ' + bootstrap + ' -v 0.0 -c 4 -a 4 -f m -i ' + outpath + 'all_nuc.phy'
    sp_code = split_seqs(inpath, outpath, protein, extension, code_table, dloop) #Save each gene in a fasta file
    run_clustalw(outpath, protein) #Align all fasta files
    join_seqs(outpath, protein) #Concatenate all alignments
    fastatophy(seqfile + '.aln', seqfile + '.phy') #Save alignments in phylip format for PhyML.
    #fastatophy(seqfile + '.aln', seqfile + '.nex', 'fasta', 'nexus', protein=protein) #Save alignemnts in nexus format, for MrBayes. (not implemented yet)
    print 'Running command:', command
    with open(outpath + log_file, 'a') as log:
        a = Popen(shlex.split(command), stdout=log, stderr=log)
        a.wait()
    if args['gene_tree']:
        aln_genes = [f for f in os.listdir(outpath) if ('.aln' in f and 'all' not in f and '.phy' not in f)]
        for gene_file in aln_genes:
            tree_file = gene_tree(outpath + gene_file, protein, bootstrap, outpath + log_file)
            final_tree = tree_file.replace('_phyml_tree.txt', '_final_tree.txt')
            tree_code(tree_file, sp_code, final_tree)
                
def split_seqs(inpath, outpath, protein, extensions, table, dloop = False):
    'if protein, translates to mitochondrial protein'
    seq_dic = {gene:[] for gene in known_genes}
    mitos = []
    for e in extensions:
        mits = [mit for mit in os.listdir(inpath) if mit.endswith(e)]
        mitos += mits
    if len(mitos) < 2:
        print 'Less than 2 files found. Check your extension and inpath flags!'
        return 0
    size = 0
    sp_dict = {}
    for n, f in enumerate(mitos):
        print 'Reading genebank file:', f #genebank file
        true_sp = ''
        try:
            i = SeqIO.read(inpath + f, 'genbank')
        except ValueError:
            print 'File', f, 'was not recognized. Check formating and genebank header.'
            raise
        sp = str(n)
        for seq in i.features:
            if seq.type == 'source':
                true_sp = '_'.join(seq.qualifiers['organism'][0].split())
                sp_dict[sp] = true_sp
            if seq.type == 'CDS':
                s = i[seq.location.start:seq.location.end].seq
                if seq.strand == -1:
                    s = s.reverse_complement()
                if protein:
                    if 'translation' in seq.qualifiers:
                        s = Seq(seq.qualifiers['translation'][0], IUPAC.protein)
                    else:
                        s = s.translate(table=table)
		try:
                    header = seq.qualifiers['gene'][0].upper()
		except:
		    header = seq.qualifiers['product'][0].upper()
                rec = SeqRecord(s, description = '', id = sp + '_' + header)
                try:
                    gene_key = gene_dict[header]
		    seq_dic[gene_key].append(rec)
                except:
                    print header + ' is not a known gene. Replace the CDS gene id with one of the following:'
                    for g in known_genes:
                        print g + ' ',
                    raise
                size += len(s)
            if (seq.type.lower() == 'd-loop' or seq.type.lower() == 'dloop') and not protein and dloop:
                s = i[seq.location.start:seq.location.end].seq
                if seq.strand == -1:
                    s = s.reverse_complement()
                header = 'DLOOP'
                rec = SeqRecord(s, description = '', id = sp + '_' + header)
                seq_dic[header].append(rec)
                size += len(s)
        if not true_sp:
            print 'File', f, 'has no source feature!' 
        size = 0

    for i in seq_dic:
        try:
            assert len(seq_dic[i]) >= len(mitos)
        except:
            print 'Warning: {0}. This gene is not present in all genbank files.({1}/{2})'.format(i, len(seq_dic[i]), len(mitos))
            print 'Gene removed.'
            continue
        a = open(outpath + i + '.fasta', 'w')
        a.close()
        SeqIO.write(seq_dic[i], outpath + i + '.fasta', 'fasta')
    with open('species_code.txt', 'w') as sp_file:
        for k in sorted(sp_dict.keys()):
            sp_file.write(k + ' ' + sp_dict[k] + '\n')
    return sp_dict
    
def run_clustalw(outpath, protein = False):

    for f in os.listdir(outpath):
        if f.endswith('.fasta'):
	    fp = outpath + f
            if not protein:
                command = 'clustalw2 -INFILE=' + fp +\
                          ' -ALIGN -OUTPUT=FASTA -OUTFILE=' + outpath + f.split('.')[0] + '_nuc.aln'
            else:
                command = 'clustalw2 -INFILE=' + fp +\
                          ' -ALIGN -TYPE=PROTEIN -OUTPUT=FASTA -OUTFILE=' + outpath + f.split('.')[0] + '_aa.aln'
            if not protein:
                command2 = 'clustalw -INFILE=' + fp +\
                          ' -ALIGN -OUTPUT=FASTA -OUTFILE=' + outpath + f.split('.')[0] + '_nuc.aln'
            else:
                command2 = 'clustalw -INFILE=' + fp +\
                          ' -ALIGN -TYPE=PROTEIN -OUTPUT=FASTA -OUTFILE=' + outpath + f.split('.')[0] + '_aa.aln'
            with open(outpath+'log_clustalw.txt', 'a') as log:
                log.write(fp + ' ' + command + '\n')
                try:
                    print 'Running command:', command 
                    a = Popen(shlex.split(command), stdout=log, stderr=log)
                    a.wait()
                except:
                    a = Popen(shlex.split(command2), stdout=log, stderr=log)
                    a.wait()

def join_seqs(path, protein = False):
    if protein:
        end = '_aa.aln'
    else:
        end = '_nuc.aln'
    sp_dic = {} #Species dict, each mitogenome is one species.
    for f in os.listdir(path):
        if f.endswith(end) and 'all' not in f:
            for seq in SeqIO.parse(path + f, 'fasta'):
                sp = seq.description.split('_')[0]
                if sp in sp_dic.keys():
                    sp_dic[sp].seq = sp_dic[sp].seq + seq.seq
                else:
                    sp_dic[sp] = SeqRecord(seq = Seq(str(seq.seq)), id = sp, description = '') #sp_dic = {species1:str(gene1)+str(gene2), species2: str(gene1)+str(gene2), ...}
    a = open(path + 'all' + end, 'w')
    a.close()
    SeqIO.write(sp_dic.values(), path + 'all' + end, 'fasta')

def fastatophy(infile, outfile, format_in = 'fasta', format_out = 'phylip', protein = True):
    seq_records = []
    with open(infile, 'r') as handle:
        i = SeqIO.parse(handle, format_in)
        for seq in i:
            if format_out == 'nexus':
                if protein:
                    seq.seq.alphabet = IUPAC.protein
                else:
                    seq.seq.alphabet = IUPAC.unambiguous_dna
            seq_records.append(seq)
    with open(outfile, 'wb') as out:
        try:
            SeqIO.write(seq_records, out, format_out)
        except:
            print 'Could not open file:', infile, '| Check your permissions.'
            raise

def gene_tree(aln_file, protein, bootstrap, log_file):
    
    if protein:
        command = 'nohup phyml -d aa -m JTT -b ' + bootstrap + ' -v 0.0 -c 4 -a 4 -f m -i '
    else:
        command = 'nohup phyml -m GTR -b ' + bootstrap + ' -v 0.0 -c 4 -a 4 -f m -i '
    phy_file = ''.join(aln_file.split('.')[:-1]) + '.phy'
    fastatophy(aln_file, phy_file)
    command = command + phy_file
    print 'Running command:', command
    with open(log_file, 'a') as log:
        a = Popen(shlex.split(command), stdout=log, stderr=log)
        a.wait()
    tree_file = phy_file + '_phyml_tree.txt'
    return tree_file

def remove_gap(aligned_fasta, outfile):
    with open(outfile, 'w') as o:
        for line in open(aligned_fasta, 'r'):
            if '>' not in line:
                o.write(line.replace('-', ''))
            else:
                o.write(line)

def file_handler(aln_file, nuc_file, outfile, alphabet = 'Vertebrate Mitochondrial'):
    out_rec = []
    for aln in SeqIO.parse(aln_file, 'fasta'):
        for nuc in SeqIO.parse(nuc_file, 'fasta'):
            if nuc.id == aln.id:
                out = bt(aln, nuc, alphabet)
                out_seq = Seq(out, generic_dna)
                out_rec.append(SeqRecord(out_seq, id = nuc.id.split('_')[0], description = nuc.description.split('_')[0]))
    SeqIO.write(out_rec, outfile, 'fasta')


def bt(aln, nuc, alphabet):
    'Back translate a nucleotidic sequence based on an amino acid (gapped) sequence'
    prot_seq = aln.seq
    nucl_seq = nuc.seq
    gaps = 0
    bt_seq = ''
    if len(nucl_seq)%3 != 0:
        print 'Nucleotide sequence is not divisible by 3, removing excess nucleotides.'
        nucl_seq = nucl_seq[:-(len(nucl_seq)%3)]
    if len(nucl_seq)/3 < len(str(prot_seq).replace('-', '')):
        print len(nucl_seq)/3, str(prot_seq).replace('-', '')
        raise ValueError('Nucleotide sequence is smaller than protein sequence times 3!')
    for n, aa in enumerate(list(prot_seq)):
        if aa == '-':
            bt_seq += '---'
            gaps += 1
        else:
            pos = n * 3 - gaps * 3
            codon = nucl_seq[pos:pos+3]
            translated = codon.translate(table=alphabet)
            if aa != str(translated):
                print 'Translation error!'
                print 'aminoacid/position:', aa, n
                print 'codon/translated', codon, translated
                print nuc.id, aln.id
                return 0
            else:
                bt_seq += str(codon)
    return bt_seq

def run_bt_mito(path):
    genes = ['ND1', 'ND2', 'COX1', 'COX2', 'ATP8', 'ATP6', 'ND4L', 'ND4', 'ND5', 'CYTB', 'ND6', 'COX3'] #sem o ND3!
    pairs = [(gene + '_aa.aln', gene + '.fasta', gene + '_codon.aln') for gene in genes]
    for p in pairs:
        file_handler(p[0], p[1], p[2])

def tree_code(tree_file, species_code, out_file):
    t = open(tree_file, 'r') 
    tree = t.read()[:-1]
    t.close()
    #print tree
    for k in sorted(species_code.keys())[::-1]:
        tree = tree.replace('(' + str(k), '(' + species_code[k])
        tree = tree.replace(',' + str(k), ',' + species_code[k])
    with open(out_file, 'w') as out:
        out.write(tree)    

if __name__ == '__main__':
    args = argument_parser()
    main(args)
