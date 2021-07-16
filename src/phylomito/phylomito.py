#!/usr/bin/env python
# -*- coding: utf-8 -*-
# phylomito.py

'''Mitochondrial phylogeny using the supermatrix method.'''

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'


import os
import shlex
import argparse
from builtins import str
from builtins import range
from subprocess import Popen
from copy import deepcopy
from Bio import SeqIO#, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


GENES = [('ND1', 'NAD1'), ('ND2', 'NAD2'), ('COX1', 'CO1'), ('COX2', 'CO2'),
         ('ATP8', 'ATPase 8'), ('ATP6', 'ATPase 6'), ('ND3', 'NAD3'),
         ('ND4L', 'NAD4L'), ('ND4', 'NAD4'), ('ND5', 'NAD5'),
         ('CYTB', 'Cyt B', 'COB'), ('ND6', 'NAD6'), ('COX3', 'CO3'), ('DLOOP',)]
KNOWN_GENES = ['ND1', 'ND2', 'COX1', 'COX2', 'ATP8', 'ATP6', 'ND3', 'ND4L',
               'ND4', 'ND5', 'CYTB', 'ND6', 'COX3', 'DLOOP']
GENE_DICT = {}
for n, gs in enumerate(GENES):
    for g in gs:
        GENE_DICT[g] = KNOWN_GENES[n]

def argument_parser(hlp=False):
    '''phylomito.py -i /path/to/genbank/ -p -o /path/to/output/
    Output default: current working directory.'''

    default_out = os.getcwd() + '/'
    parser = argparse.ArgumentParser(description='Phylomito is a simple pipeline\
                                     to automatize mitochondrial super-matrix\
                                     phylogenomic, using Clustal and PhyML.',\
                                     argument_default=None, fromfile_prefix_chars='@')
    parser.add_argument('-i', '--inpath', nargs='?', type=str, required=True,\
                        dest='inpath', help='Path to the folder with genbank\
                        sequences. (default: %(default)s)')
    parser.add_argument('-o', '--outpath', nargs='?', type=str, default=default_out,\
                        dest='outpath', help='Path were the alignments and phylogenetic\
                        tree will be saved. (default: %(default)s)')
    parser.add_argument('-l', '--log', nargs='?', type=str, default='phyml.log',\
                        dest='log', help='Log file. (default: %(default)s)')
    parser.add_argument('-e', '--extension', nargs='*', type=str, default=['.gbk', '.gb'],\
                        dest='extension', help='Extension for the genbank files.\
                        (default: %(default)s)')
    parser.add_argument('-b', '--bootstrap', nargs='?', type=int, default=100,\
                        dest='bootstrap', help='Number of bootstrap repetitions\
                        on PhyML. (default: %(default)s)')
    parser.add_argument('-n', '--nucleotide', nargs='?', const=True, default=False,\
                        dest='nucleotide', help='Set this flag for nucleotide sequences\
                        alignment and phylogeny. Default use protein sequences.')
    parser.add_argument('-s', '--skip_phyml', nargs='?', const=True, default=False,\
                        dest='skip_phyml', help='Set this flag to skip the phylogenetic\
                        tree reconstruction and stop the analysis after the alignment.\
                        (default: %(default)s).')
    parser.add_argument('-g', '--gene_tree', nargs='?', const=True, default=False,\
                        dest='gene_tree', help='Set this flag if you want to make a\
                        tree for every gene. (default: %(default)s)')
    parser.add_argument('-m', '--model', nargs='?', type=str, default='JTT',\
                        dest='model', help='Substitution model to use with phyml.\
                        default to GTR with nucleotide sequences and JTT with proteins. (default: %(default)s)')
    parser.add_argument('-d', '--dloop', nargs='?', const=True, default=False,\
                        dest='dloop', help='Flag to include DLOOP region in the\
                        alignment. (default: %(default)s)')
    parser.add_argument('-t', '--code_table', nargs='?', type=int, default=2,\
                        dest='code_table', help='Genetic table to be used for\
                        translation. Only matters if --protein flag is used.\
                        Default is Vertebrate Mitochondrial Code. See all tables\
                        at http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.\
                        (default: %(default)s)')
    if hlp:
        args = parser.parse_args(['-h'])
    else:
        args = parser.parse_args().__dict__
    return args

def main():
    args = argument_parser()
    #Args processing.
    protein = not args['nucleotide']
    if not protein and args['model'] == 'JTT':
        model = 'GTR'
    else:
        model = args['model']
    inpath = args['inpath']
    code_table = args['code_table']
    log_file = args['log']
    skip_phyml = args['skip_phyml']
    if not inpath.endswith('/'):
        inpath += '/'
    outpath = args['outpath']
    if not outpath.endswith('/'):
        outpath += '/'
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    extension = args['extension']
    dloop = args['dloop']
    bootstrap = str(args['bootstrap'])
    #All mitochondrial codes
    if code_table not in [2, 3, 4, 5, 9, 13, 14, 16, 21, 22, 23, 24]:
        print('WARNING: Genetic code n.{} is not Mitochondial!'.format('code_table'))
    #if not dloop:
    #    KNOWN_GENES.remove('DLOOP')
    try:
        a = open(outpath + log_file, 'w')
        a.close()
    except:
        print('Was not able to open {}. Check your permissions.'.format(outpath + log_file))
        raise
    #seqfile will be the file with the concatenated alignment.
    if protein:
        seqfile = outpath + 'all_aa'
        command = 'phyml -d aa -m JTT -b ' + bootstrap + ' -v 0.0 -c 4 -a 4 -f m -i '\
                  + outpath + 'all_aa.phy'
    else:
        seqfile = outpath + 'all_nuc'
        command = 'phyml -m GTR -b ' + bootstrap + ' -v 0.0 -c 4 -a 4 -f m -i '\
                  + outpath + 'all_nuc.phy'
    #Save each gene in a fasta file
    sp_list = split_seqs(inpath, outpath, protein, extension, code_table, dloop) 
    run_clustalw(outpath, protein) #Align all fasta files
    join_seqs(outpath, protein) #Concatenate all alignments
    if not skip_phyml:
        fastatophy(seqfile + '.aln.fasta', seqfile + '.phy') #Save alignments in phylip format for PhyML.
        #Save alignemnts in nexus format, for MrBayes. (not implemented yet)
        #fastatophy(seqfile + '.aln.fasta', seqfile + '.nex', 'fasta', 'nexus', protein=protein) 
        print('Running command:', command)
        with open(outpath + log_file, 'a') as log:
            a = Popen(shlex.split(command), stdout=log, stderr=log)
            a.wait()
        tree_code(seqfile + '.phy_phyml_tree.txt', sp_list, seqfile + '.phy_final_tree.txt')
        if args['gene_tree']:
            aln_genes = [f for f in os.listdir(outpath) if ('.aln.fasta' in f and 'all' not in f and '.phy' not in f)]
            for gene_file in aln_genes:
                print(outpath+gene_file)
                tree_file = gene_tree(outpath + gene_file, protein, bootstrap, outpath + log_file)
                final_tree = tree_file.replace('_phyml_tree.txt', '_final_tree.txt')
                tree_code(tree_file, sp_list, final_tree)
                
def find_files(inpath, extensions):
    fastatophy(seqfile + '.aln.fasta', seqfile + '.phy') #Save alignments in phylip format for PhyML.
    #fastatophy(seqfile + '.aln.fasta', seqfile + '.nex', 'fasta', 'nexus', protein=protein) #Save alignemnts in nexus format, for MrBayes. (not implemented yet)
    print('Running command:', command)
    with open(outpath + log_file, 'a') as log:
        a = Popen(shlex.split(command), stdout=log, stderr=log)
        a.wait()
    tree_code(seqfile + '.phy_phyml_tree.txt', sp_list, seqfile + '.phy_final_tree.txt')
    if args['gene_tree']:
        aln_genes = [f for f in os.listdir(outpath) if ('.aln.fasta' in f and 'all' not in f and '.phy' not in f)]
        for gene_file in aln_genes:
            print(outpath+gene_file)
            tree_file = gene_tree(outpath + gene_file, protein, bootstrap, outpath + log_file)
            final_tree = tree_file.replace('_phyml_tree.txt', '_final_tree.txt')
            tree_code(tree_file, sp_list, final_tree)
                
def split_seqs(inpath, outpath, protein, extensions, table, dloop=False):
    'if protein, translates to mitochondrial protein'
    genes = deepcopy(KNOWN_GENES)
    if not dloop and 'DLOOP' in genes:
        genes.remove('DLOOP')
    seq_dic = {gene:[] for gene in genes}
    mitos = []
    for e in extensions:
        mits = [mit for mit in os.listdir(inpath) if mit.endswith(e)]
        mitos += mits
    if len(mitos) < 2:
        print('Less than 2 files found. Check your extension and inpath flags!')
        return 0
    mitos.sort()
    size = 0
    sp_list = []
    present_genes = {gene:False for gene in genes}
    mitos.sort()
    for n, f in enumerate(mitos):
        print('Reading genebank file:', f) #genebank file
        true_sp = ''
        try:
            genebank_data = SeqIO.read(inpath + f, 'genbank')
        except ValueError:
            print('File', f, 'was not recognized. Check formating and genebank header.')
            raise
        sp = str(n) #species code
        for seq in genebank_data.features:
            if seq.type == 'source':
                true_sp = '_'.join(seq.qualifiers['organism'][0].split())
                sp_list.append(true_sp)
            if seq.type == 'CDS':
                s = genebank_data[seq.location.start:seq.location.end].seq
                if seq.strand == -1:
                    s = s.reverse_complement()
                if protein:
                    if 'translation' in seq.qualifiers:
                        s = Seq(seq.qualifiers['translation'][0])
                    else:
                        s = s.translate(table=table)
                try:
                    header = seq.qualifiers['gene'][0].upper()
                except:
                    header = seq.qualifiers['product'][0].upper()
                rec = SeqRecord(s, description='', id=sp + '_' + header)
                try:
                    gene_key = GENE_DICT[header]
                    seq_dic[gene_key].append(rec)
                    present_genes[gene_key] = True
                except:
                    print(header + ' is not a known gene. Replace the CDS gene id\
                                    with one of the following:')
                    for gene_name in KNOWN_GENES:
                        print(gene_name + ' ', end=' ')
                    raise
                size += len(s) #size of the multiple alignment
            if seq.type == 'misc_feature' and dloop:
                if 'control region' in list(seq.qualifiers.values())[0]:
                    s = genebank_data[seq.location.start:seq.location.end].seq
                    if seq.strand == -1:
                        s = s.reverse_complement()
                    header = 'DLOOP'
                    rec = SeqRecord(s, description='', id=sp + '_' + header)
                    seq_dic[header].append(rec)
                    size += len(s)
            if (seq.type.lower() == 'd-loop' or seq.type.lower() == 'dloop') and not protein and dloop:
                s = genebank_data[seq.location.start:seq.location.end].seq
                if seq.strand == -1:
                    s = s.reverse_complement()
                header = 'DLOOP'
                rec = SeqRecord(s, description='', id=sp + '_' + header)
                seq_dic[header].append(rec)
                size += len(s)
        if not true_sp:
            print('File', f, 'has no source feature!') 
        size = 0
    
        for gene in present_genes:
            if not present_genes[gene]:
                print('File ', f, ' is missing gene ', gene)
    for gene in seq_dic:
        try:
            assert len(seq_dic[gene]) >= len(mitos)
        except AssertionError:
            print('Warning: {0}. This gene is not present in all genbank files.({1}/{2})'\
                  .format(gene, len(seq_dic[gene]), len(mitos)))
            print('Gene removed.')
            continue
        a = open(outpath + gene + '.fasta', 'w')
        a.close()
        SeqIO.write(seq_dic[gene], outpath + gene + '.fasta', 'fasta')
    with open(outpath + 'species_code.txt', 'w') as sp_file:
        for n, sp in enumerate(sp_list):
            sp_file.write(str(n) + ' ' + sp + '\n')
    return sp_list
    
def run_clustalw(outpath, protein=False):

    for f in os.listdir(outpath):
        if f.endswith('.fasta'):
            file_path = outpath + f
            if not protein:
                command = 'clustalo -INFILE=' + file_path +\
                          ' -ALIGN -OUTPUT=FASTA -OUTFILE='\
                          + outpath + f.split('.')[0] + '_nuc.aln.fasta'
            else:
                command = 'clustalo -INFILE=' + file_path +\
                          ' -ALIGN -TYPE=PROTEIN -OUTPUT=FASTA -OUTFILE='\
                          + outpath + f.split('.')[0] + '_aa.aln.fasta'
            if not protein:
                command2 = 'clustalw -INFILE=' + file_path +\
                          ' -ALIGN -OUTPUT=FASTA -OUTFILE='\
                          + outpath + f.split('.')[0] + '_nuc.aln.fasta'
            else:
                command2 = 'clustalw -INFILE=' + file_path +\
                          ' -ALIGN -TYPE=PROTEIN -OUTPUT=FASTA -OUTFILE='\
                          + outpath + f.split('.')[0] + '_aa.aln.fasta'
            with open(outpath+'log_clustalw.txt', 'a') as log:
                log.write(file_path + ' ' + command + '\n')
                try:
                    print('Running command:', command) 
                    a = Popen(shlex.split(command), stdout=log, stderr=log)
                    a.wait()
                except:
                    a = Popen(shlex.split(command2), stdout=log, stderr=log)
                    a.wait()

def join_seqs(path, protein=False):
    if protein:
        end = '_aa.aln.fasta'
    else:
        end = '_nuc.aln.fasta'
    sp_dic = {} #Species dict, each mitogenome is one species.
    for f in os.listdir(path):
        if f.endswith(end) and 'all' not in f:
            for seq in SeqIO.parse(path + f, 'fasta'):
                sp = seq.description.split('_')[0]
                if sp in list(sp_dic.keys()):
                    sp_dic[sp].seq = sp_dic[sp].seq + seq.seq
                else:
                    #sp_dic = {species1:str(gene1)+str(gene2), species2: str(gene1)+str(gene2),}
                    sp_dic[sp] = SeqRecord(seq=Seq(str(seq.seq)), id=sp, description='') 
    clear_file = open(path + 'all' + end, 'w')
    clear_file.close()
    SeqIO.write(list(sp_dic.values()), path + 'all' + end, 'fasta')

def fastatophy(infile, outfile, format_in='fasta', format_out='phylip', protein=True):
    seq_records = []
    with open(infile, 'r') as handle:
        in_parsed = SeqIO.parse(handle, format_in)
        for seq in in_parsed:
            seq_records.append(seq)
    with open(outfile, 'w') as out:
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
        run = Popen(shlex.split(command), stdout=log, stderr=log)
        run.wait()
    tree_file = phy_file + '_phyml_tree.txt'
    return tree_file

def remove_gap(aligned_fasta, outfile):
    with open(outfile, 'w') as out:
        for line in open(aligned_fasta, 'r'):
            if '>' not in line:
                out.write(line.replace('-', ''))
            else:
                out.write(line)

def file_handler(aln_file, nuc_file, outfile, alphabet='Vertebrate Mitochondrial'):
    out_rec = []
    for aln in SeqIO.parse(aln_file, 'fasta'):
        for nuc in SeqIO.parse(nuc_file, 'fasta'):
            if nuc.id == aln.id:
                out = back_translate(aln, nuc, alphabet)
                out_seq = Seq(out)
                out_rec.append(SeqRecord(out_seq, id=nuc.id.split('_')[0],\
                               description=nuc.description.split('_')[0]))
    SeqIO.write(out_rec, outfile, 'fasta')


def back_translate(aln, nuc, alphabet):
    'Back translate a nucleotidic sequence based on an amino acid (gapped) sequence'
    prot_seq = aln.seq
    nucl_seq = nuc.seq
    gaps = 0
    bt_seq = ''
    if len(nucl_seq)%3 != 0:
        print('Nucleotide sequence is not divisible by 3, removing excess nucleotides.')
        nucl_seq = nucl_seq[:-(len(nucl_seq)%3)]
    if len(nucl_seq) // 3 < len(str(prot_seq).replace('-', '')):
        print(len(nucl_seq) // 3, str(prot_seq).replace('-', ''))
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
            bt_seq += str(codon)
    return bt_seq

def run_bt_mito():
    genes = deepcopy(KNOWN_GENES)
    genes.remove('ND3')
    try:
        genes.remove('DLOOP')
    except ValueError:
        pass
    pairs = [(gene + '_aa.aln.fasta', gene + '.fasta', gene + '_codon.aln.fasta') for gene in genes]
    for pair in pairs:
        file_handler(pair[0], pair[1], pair[2])

def tree_code(tree_file, species_list, out_file):
    tree_temp = open(tree_file, 'r') 
    tree = tree_temp.read()[:-1]
    tree_temp.close()
    #print tree
    for key in reversed(list(range(len(species_list)))):
        tree = tree.replace('(' + str(key), '(' + species_list[key])
        tree = tree.replace(',' + str(key), ',' + species_list[key])
    with open(out_file, 'w') as out:
        out.write(tree)    

if __name__ == '__main__':
    main()
