#!/usr/bin/env python
# -*- coding: utf-8 -*-
# prtolog_phylo.py

'''Orthologs phylogeny supermatrix method.'''
from __future__ import print_function
from __future__ import division

from builtins import zip
from builtins import str
from builtins import range
from past.utils import old_div
__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

import os
import sys
import shlex
import argparse
from subprocess import Popen
from copy import deepcopy
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna



def main(seq_file, genomes):
    #Args processing.
    outpath = os.getcwd() +'/'
    log_file = 'ortholog.log'
    protein = True
    bootstrap = '100'
    try:
        a = open(outpath + log_file, 'w')
        a.close()
    except:
        print('Was not able to open {}. Check your permissions.'.format(outpath + log_file))
        raise
    seq_names = parse_seqs(seq_file, genomes)
    file_names = separate_seqs(seq_names, genomes) #Save each gene in a fasta file
    run_clustalw(outpath, protein) #Align all fasta files
    with open('file_names', 'w') as fn:
        for f in file_names:
            fn.write(f + '\n')
    for f in file_names:
        fastatophy(f + '_aa.aln', f + '.phy') #Save alignments in phylip format for PhyML.
        command = 'phyml -d aa -m JTT -b ' + bootstrap + ' -v 0.0 -c 4 -a 4 -f m -i ' + f + '.phy'
        print('Running command:', command)
        with open(outpath + log_file, 'a') as log:
            a = Popen(shlex.split(command), stdout=log, stderr=log)
            a.wait()
    joined_file = join_seqs(outpath, protein) #Concatenate all alignments
    #joined_file = 'all_aa.aln'
    fastatophy(joined_file, joined_file + '.phy') #Save alignments in phylip format for PhyML.
    command = 'phyml -d aa -m JTT -b ' + bootstrap + ' -v 0.0 -c 4 -a 4 -f m -i ' + joined_file + '.phy'
    with open(outpath + log_file, 'a') as log:
        a = Popen(shlex.split(command), stdout=log, stderr=log)
        a.wait()
                
def parse_seqs(seq_file, genomes):
    seq_names = []
    with open(seq_file, 'r') as ortho_file:
        for l in ortho_file:
            seq_ids = [s.split('(')[0] for s in l.split(':')[1].split()]
            assert len(seq_ids) == len(genomes)
            seq_names.append(seq_ids)
    return seq_names

def separate_seqs(seq_names, genomes):
    file_names = []
    for s in seq_names:
        with open(s[0] + '.fasta', 'w') as fasta_file:
            for seq_name, genome in zip(s, genomes):
                seq = find(seq_name, genome)
                #print seq_name, seq[0]
                #assert seq_name == seq[0]
                header = '>'+genome.split('/')[-1]
                fasta_file.write(header + '\n' + seq[1] + '\n')
            file_names.append(s[0])
    return file_names

def find(seq_name, genome):
    '''return the first sequence that contains seq_name on it's fasta ID'''
    seq = ''
    seq_id = ''
    with open(genome, 'r') as f:
        for seq_line in f:
            if seq_line[0] == '>':
                if seq and seq_name in seq_id:
    	            return seq_id, seq
                seq_id = seq_line[1:-1]
                seq = ''
            else:
                seq += seq_line[:-1]
        if seq_name in seq_id:
            return seq_id, seq
            
def run_clustalw(outpath, protein=False):

    for f in os.listdir(outpath):
        if f.endswith('.fasta'):
            fp = outpath + f
            print(fp)
            if not protein:
                command = 'clustalw2 -INFILE=' + fp +\
                          ' -ALIGN -OUTPUT=FASTA -OUTFILE=' + outpath + '.'.join(f.split('.')[:-1]) + '_nuc.aln'
            else:
                command = 'clustalw2 -INFILE=' + fp +\
                          ' -ALIGN -TYPE=PROTEIN -OUTPUT=FASTA -OUTFILE=' + outpath + '.'.join(f.split('.')[:-1]) + '_aa.aln'
            if not protein:
                command2 = 'clustalw -INFILE=' + fp +\
                          ' -ALIGN -OUTPUT=FASTA -OUTFILE=' + outpath + '.'.join(f.split('.')[:-1]) + '_nuc.aln'
            else:
                command2 = 'clustalw -INFILE=' + fp +\
                          ' -ALIGN -TYPE=PROTEIN -OUTPUT=FASTA -OUTFILE=' + outpath + '.'.join(f.split('.')[:-1]) + '_aa.aln'
            with open(outpath+'log_clustalw.txt', 'a') as log:
                log.write(fp + ' ' + command + '\n')
                try:
                    print('Running command:', command) 
                    a = Popen(shlex.split(command), stdout=log, stderr=log)
                    a.wait()
                except:
                    a = Popen(shlex.split(command2), stdout=log, stderr=log)
                    a.wait()

def join_seqs(path, protein=False):
    if protein:
        end = '_aa.aln'
    else:
        end = '_nuc.aln'
    sp_dic = {} #Species dict, each mitogenome is one species.
    for f in os.listdir(path):
        if f.endswith(end) and 'all' not in f:
            for seq in SeqIO.parse(path + f, 'fasta'):
                sp = seq.description[:4]
                if sp in list(sp_dic.keys()):
                    sp_dic[sp].seq = sp_dic[sp].seq + seq.seq
                else:
                    sp_dic[sp] = SeqRecord(seq = Seq(str(seq.seq)), id = sp, description = '') #sp_dic = {species1:str(gene1)+str(gene2), species2: str(gene1)+str(gene2), ...}
    outfile = path + 'all' + end
    a = open(outfile, 'w')
    a.close()
    SeqIO.write(list(sp_dic.values()), outfile, 'fasta')
    return outfile

def fastatophy(infile, outfile, format_in='fasta', format_out='phylip', protein=True):
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
            print('Could not open file:', infile, '| Check your permissions.')
            raise

def gene_tree(aln_file, protein, bootstrap, log_file):
    
    if protein:
        command = 'nohup phyml -d aa -m JTT -b ' + bootstrap + ' -v 0.0 -c 4 -a 4 -f m -i '
    else:
        command = 'nohup phyml -m GTR -b ' + bootstrap + ' -v 0.0 -c 4 -a 4 -f m -i '
    phy_file = aln_file[:-4] + '.phy'
    fastatophy(aln_file, phy_file)
    command = command + phy_file
    print('Running command:', command)
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

def file_handler(aln_file, nuc_file, outfile, alphabet='Vertebrate Mitochondrial'):
    out_rec = []
    for aln in SeqIO.parse(aln_file, 'fasta'):
        for nuc in SeqIO.parse(nuc_file, 'fasta'):
            if nuc.id == aln.id:
                out = bt(aln, nuc, alphabet)
                out_seq = Seq(out, generic_dna)
                out_rec.append(SeqRecord(out_seq, id=nuc.id.split('_')[0], description=nuc.description.split('_')[0]))
    SeqIO.write(out_rec, outfile, 'fasta')


def bt(aln, nuc, alphabet):
    'Back translate a nucleotidic sequence based on an amino acid (gapped) sequence'
    prot_seq = aln.seq
    nucl_seq = nuc.seq
    gaps = 0
    bt_seq = ''
    if len(nucl_seq)%3 != 0:
        print('Nucleotide sequence is not divisible by 3, removing excess nucleotides.')
        nucl_seq = nucl_seq[:-(len(nucl_seq)%3)]
    if old_div(len(nucl_seq),3) < len(str(prot_seq).replace('-', '')):
        print(old_div(len(nucl_seq),3), str(prot_seq).replace('-', ''))
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
                print('Translation error!')
                print('aminoacid/position:', aa, n)
                print('codon/translated', codon, translated)
                print(nuc.id, aln.id)
                return 0
            else:
                bt_seq += str(codon)
    return bt_seq

def run_bt_mito(path):
    genes = deepcopy(KNOWN_GENES)
    gene.remove('ND3')
    try:
        gene.remove('DLOOP')
    except ValueError:
        pass
    pairs = [(gene + '_aa.aln', gene + '.fasta', gene + '_codon.aln') for gene in genes]
    for p in pairs:
        file_handler(p[0], p[1], p[2])

def tree_code(tree_file, species_list, out_file):
    t = open(tree_file, 'r') 
    tree = t.read()[:-1]
    t.close()
    #print tree
    for k in reversed(list(range(len(species_list)))):
        tree = tree.replace('(' + str(k), '(' + species_list[k])
        tree = tree.replace(',' + str(k), ',' + species_list[k])
    with open(out_file, 'w') as out:
        out.write(tree)    

if __name__ == '__main__':
    #seq file is at /data2/Marcela/marcela.fasta.an/4genes.4taxa.no.gallo 
    main(sys.argv[1], sys.argv[2:]) #seq_file, genomes
    #args = argument_parser()
    #main(args)
