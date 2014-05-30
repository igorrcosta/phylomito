#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Filogenia mitocondrial usando super-matriz'''

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

'''todo:
1) Fazer a filogenia com mitocondria completa e por alinhamento de codons.
2) Adicionar modeltest
3) Suporte para skipped bases no ND3
4) Adicionar codigo de especies.
5) Arvore só do DLOOP
6) Implementar NJ'''

genes = ['ND1', 'ND2', 'COX1', 'COX2', 'ATP8', 'ATP6', 'ND3', 'ND4L', 'ND4', 'ND5', 'CYTB', 'ND6', 'COX3', 'DLOOP']

def argument_parser(hlp = False):
    '''filomito.py -i /path/to/genbank/ -p -o /path/to/output/
    Output default: current working directory.'''

    default_out = os.getcwd() + '/'
    parser = argparse.ArgumentParser(description = 'Filomito is a simple pipeline to automatize mitochondrial super-matrix phylogenomic, using clustaw, phyML and mrbayes.',\
                                     argument_default = None, fromfile_prefix_chars = '@')
    parser.add_argument('-i', '--inpath', nargs = '?', type = str, required = True,\
                        dest = 'inpath', help = 'Path to the folder with genbank sequences. (default: %(default)s)')
    parser.add_argument('-o', '--outpath', nargs = '?', type = str, default = default_out,\
                        dest = 'outpath', help = 'Path were the alignments and phylogenetic tree will be saved. (default: %(default)s)')
    parser.add_argument('-e', '--extension', nargs = '?', type = str, default = '.gbk',\
                        dest = 'extension', help = 'Extension for the genbank files. (default: %(default)s)')
    parser.add_argument('-p', '--protein', nargs = '?', const = True, default = False,\
                        dest = 'protein', help = 'Set this flag for protein sequences alignment and phylogeny. (default: %(default)s)')
    parser.add_argument('-g', '--gene_tree', nargs = '?', const = True, default = False,\
                        dest = 'gene_tree', help = 'Set this flag if you want to make a tree for every gene. (default: %(default)s)')
    parser.add_argument('-d', '--dloop', nargs = '?', const = True, default = False,\
                        dest = 'dloop', help = 'Flag to include DLOOP region in the alignment. (default: %(default)s)')
    if hlp:
        args = parser.parse_args(['-h'])
    else:
        args = parser.parse_args().__dict__
    return args

def main(args):
    protein = args['protein']
    inpath = args['inpath']
    outpath = args['outpath']
    extension = args['extension']
    dloop = args['dloop']
    if not dloop:
        genes.remove('DLOOP')
    separador(inpath, outpath, protein, extension, dloop)
    run_clustalw(outpath, protein)
    junta(outpath, protein)
    if protein:
        fastatophy(outpath + 'all_aa.aln', outpath + 'all_aa.phy')
        fastatophy(outpath + 'all_aa.aln', outpath + 'all_aa.nex', 'fasta', 'nexus')
    else:
        fastatophy(outpath + 'all_nuc.aln', outpath + 'all_nuc.phy')
        fastatophy(outpath + 'all_nuc.aln', outpath + 'all_nuc.nex', 'fasta', 'nexus', protein = False)
    command_nuc = 'nohup phyml -m GTR -b 100 -v 0.0 -c 4 -a 4 -f m -i ' + outpath + 'all_nuc.phy'
    command_aa = 'nohup phyml -d aa -m JTT -b 100 -v 0.0 -c 4 -a 4 -f m -i '+ outpath + 'all_aa.phy'
    with open(outpath + 'log_phyml.txt', 'a') as log:
        if protein:
            log.write(command_aa + '\n')
            a = Popen(shlex.split(command_aa), stdout=log, stderr=log)
            a.wait()
        else:
            log.write(command_nuc + '\n')
            a = Popen(shlex.split(command_nuc), stdout=log, stderr=log)
            a.wait()
    if args['gene_tree']:
        gene_tree(outpath)

def separador(inpath, outpath, protein = False, extension = '.gbk', dloop = False):
    'If protein, translates with certebrate mitochondrial code'
    seq_dic = {gene:[] for gene in genes}
    mitos = [mit for mit in os.listdir(inpath) if mit.endswith(extension)]
    size = 0
    spec_dict = {}
    for n, f in enumerate(mitos):
        true_spec = ''
        i = SeqIO.read(inpath + f, 'genbank')
        for seq in i.features:
            if seq.type == 'source':
                true_spec = '_'.join(seq.qualifiers['organism'][0].split())
                spec = str(n)
                spec_dict[spec] = true_spec
            if seq.type == 'CDS':
                s = i[seq.location.start:seq.location.end].seq
                if seq.strand == -1:
                    s = s.reverse_complement()
                if protein:
                    s = s.translate(table="Vertebrate Mitochondrial")
                header = seq.qualifiers['gene'][0]
                rec = SeqRecord(s, description = '', id = spec + '_' + header)
                seq_dic[header].append(rec)
                size += len(s)
            if seq.type == 'misc_feature' and not protein and dloop:
                s = i[seq.location.start:seq.location.end].seq
                if seq.strand == -1:
                    s = s.reverse_complement()
                header = 'DLOOP'
                rec = SeqRecord(s, description = '', id = spec + '_' + header)
                seq_dic[header].append(rec)
                size += len(s)
        if not true_spec:
            print 'File', f, 'has no source feature!'
        size = 0

    for i in seq_dic:
        try:
            assert len(seq_dic[i]) == len(mitos)
        except:
            print 'Warning: {0}. This gene is not present in all genbank files.({1}/{2})'.format(i, len(seq_dic[i]), len(mitos))
            print 'Gene removed.'
            continue
        a = open(outpath + i + '.fasta', 'w')
        a.close()
        SeqIO.write(seq_dic[i], outpath + i + '.fasta', 'fasta')
    pprint(spec_dict)

def run_clustalw(path, protein = False):

    for f in os.listdir(path):
        if f.endswith('.fasta'):
            os.chdir(path)
            fp = path + f
            if not protein:
                command = 'clustalw -INFILE=' + fp +\
                          ' -ALIGN -OUTPUT=FASTA -OUTFILE=' + fp.split('.')[0] + '_nuc.aln'
            else:
                command = 'clustalw -INFILE=' + fp +\
                          ' -ALIGN -TYPE=PROTEIN -OUTPUT=FASTA -OUTFILE=' + fp.split('.')[0] + '_aa.aln'
            with open('log.txt', 'a') as log:
                log.write(fp + ' ' + command)
                a = Popen(shlex.split(command), stdout=log, stderr=log)
                a.wait()

def junta(path, protein = False):
    if protein:
        end = '_aa.aln'
    else:
        end = '_nuc.aln'
    spec_dic = {}
    for f in os.listdir(path):
        if f.endswith(end):
            for seq in SeqIO.parse(path + f, 'fasta'):
                spec = seq.description.split('_')[0]
                if spec in spec_dic.keys():
                    spec_dic[spec].seq = spec_dic[spec].seq + seq.seq
                else:
                    spec_dic[spec] = SeqRecord(seq = Seq(str(seq.seq)), id = spec, description = '') #spec_dic = {especie1:str(gene1)+str(gene2), especie2: str(gene1)+str(gene2), ...}
    if protein:
        a = open(path + 'all_aa.aln', 'w')
        a.close()
        SeqIO.write(spec_dic.values(), path + 'all_aa.aln', 'fasta')
    else:
        a = open(path + 'all_nuc.aln', 'w')
        a.close()
        SeqIO.write(spec_dic.values(), path + 'all_nuc.aln', 'fasta')

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
            print infile
            raise

def gene_tree(path):
    aln_genes = [f for f in os.listdir(path) if ('.aln' in f and 'all' not in f and '.phy' not in f)]
    for f in aln_genes:
        command = 'nohup phyml -d aa -m JTT -b 100 -v 0.0 -c 4 -a 4 -f m -i '
        fastatophy(path + f, path + ''.join(f.split('.')[:-1]) + '.phy')
        command = command + path + ''.join(f.split('.')[:-1]) + '.phy'
        with open(path + 'log_phyml.txt', 'a') as log:
            log.write(command + '\n')
            a = Popen(shlex.split(command), stdout=log, stderr=log)
            a.wait()

def tira_gap(aligned_fasta, outfile):
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
                if out == 0:
                    print nuc.id, aln.id
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
        print 'seq de nucleotideos não é multipla de três, editando sequencia'
        nucl_seq = nucl_seq[:-(len(nucl_seq)%3)]
    if len(nucl_seq)/3 < len(str(prot_seq).replace('-', '')):
        print len(nucl_seq)/3, str(prot_seq).replace('-', '')
        raise ValueError('seq de nucleotídeos/3 é menor que a seq de proteína!')
    for n, aa in enumerate(list(prot_seq)):
        if aa == '-':
            bt_seq += '---'
            gaps += 1
        else:
            pos = n * 3 - gaps * 3
            codon = nucl_seq[pos:pos+3]
            translated = codon.translate(table=alphabet)
            if aa != str(translated):
                print 'aminoacid/position:', aa, n
                print 'codon/translated', codon, translated
                return 0
            else:
                bt_seq += str(codon)
    return bt_seq

def run_bt_mito(path):
    genes = ['ND1', 'ND2', 'COX1', 'COX2', 'ATP8', 'ATP6', 'ND4L', 'ND4', 'ND5', 'CYTB', 'ND6', 'COX3'] #sem o ND3!
    pairs = [(gene + '_aa.aln', gene + '.fasta', gene + '_codon.aln') for gene in genes]
    for p in pairs:
        file_handler(p[0], p[1], p[2])


if __name__ == '__main__':
    args = argument_parser()
    main(args)
