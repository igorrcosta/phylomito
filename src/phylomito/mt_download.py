import argparse
import pickle
import gzip
from os import getcwd, listdir, remove, mkdir
from os.path import join as os_join, isfile, isdir, abspath, dirname
from shutil import copy, copyfileobj
from time import sleep
from urllib.request import urlopen
from urllib.error import HTTPError
from Bio import Entrez, SeqIO
from anytree import PreOrderIter

def argument_parser():

    default_out = getcwd() + '/'
    data_path = os_join(abspath(dirname(__file__)), 'data/')
    parser = argparse.ArgumentParser(description = 'Downloads mitochondrial genomes.',\
                                     argument_default = None, fromfile_prefix_chars = '@')
    parser.add_argument('-t', '--taxonomy', nargs = '?', type = str, required = True,\
                        dest = 'taxa', help = '''Taxonomy id or scientific name of the clade of interest.
                         The program will download all refseq mitochondrial genomes from this clade.''')
    parser.add_argument('-o', '--outpath', nargs = '?', type = str, default = default_out,\
                        dest = 'outpath', help = 'Folder were the mitogenomes\' genbank files will be saved. (default: %(default)s)')
    parser.add_argument('-a', '--download_all', action='store_true', dest = 'download_all',\
                         help = '''Download and process all mitochondria from refseq and save them on the mitogenomes (-m) path.
                         This may take a few minutes (<10).''')
    parser.add_argument('-m', '--mitogenomes_path', nargs = '?', type = str, default = os_join(data_path, 'refseq/'),\
                        dest = 'mitogenomes_path', help = 'Path to the folder with refseq mitogenomes.')
    parser.add_argument('-d', '--dictionary_path', nargs = '?', type = str, default = os_join(data_path, 'taxid_to_filenames'),\
                        dest = 'dictionary_path', help = 'Path to the taxid_to_filenames dictionary. (default: %(default)s)')
    parser.add_argument('-p', '--phylogeny_path', nargs = '?', type = str, default = os_join(data_path, 'mitochondria_tree'),\
                        dest = 'phylogeny_path', help = 'Path to the mitochondrial taxonomy tree dump. (default: %(default)s)')
    args = parser.parse_args().__dict__
    return args

def main():
    args = argument_parser()
    #print(args)
    assert isfile(args['phylogeny_path'])
    assert isfile(args['dictionary_path'])
    assert isdir(args['outpath'])
    if not isdir(args['mitogenomes_path']):
        mkdir(args['mitogenomes_path'])
    if args['download_all']:
        download_and_parse(args)
        return
    with open(args['phylogeny_path'], 'rb') as tree_dump:
        mito_tree = pickle.load(tree_dump) 
    with open(args['dictionary_path'], 'rb') as dict_dump:
        taxid_file_names = pickle.load(dict_dump) 
    if not args['taxa'].isnumeric():
        print(f'Searching for organism {args["taxa"]} on NCBI Taxonomy')
        taxid = get_taxid(args['taxa'])
        if taxid:
            args['taxa'] = taxid
            print(f'Found taxonomy id {args["taxa"]}')
        else:
            print(f'Organism not found')
            return
    if args['taxa'] not in mito_tree:
        print('This taxa has no mitochondrial genomes on refseq')
        return
    target_mitogenomes = iterate_tree(mito_tree, args['taxa'], taxid_file_names)
    downloaded_mitogenomes = listdir(args['mitogenomes_path'])
    for m in target_mitogenomes:
        if m not in downloaded_mitogenomes:
            gid = '.'.join(m.split('.')[:-1])
            download_genbank(gid, args['outpath'])
        else:
            copy(args['mitogenomes_path'] + m, args['outpath'] + m)

def download_and_parse(args):
    max_attempts = 5
    attempts = 0
    interval = 10
    base_url = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/'
    files = ['mitochondrion.1.genomic.gbff.gz', 'mitochondrion.2.genomic.gbff.gz']
    for f in files:
        while attempts < max_attempts:
            sleep(interval)
            try:
                response = urlopen(base_url + f, timeout = 5)
                content = response.read()
                with open(f, 'wb') as outfile:
                    outfile.write(content)
                    break
            except HTTPError as e:
                attempts += 1
                print(type(e))
    gb_files = []
    for f in files:
        gb_file = '.'.join(f.split('.')[:-1])
        gb_files.append(gb_file)
        with gzip.open(f, 'rb') as infile, open(gb_file, 'wb') as outfile:
            for l in infile:
                outfile.write(l)
    for f in gb_files:
        all_rec = SeqIO.index(f, 'genbank')
        for rec in all_rec:
            file_name = rec.split('.')[0] + '.gb'
            with open(os_join(args['mitogenomes_path'], file_name), 'wb') as outfile:
                outfile.write(all_rec.get_raw(rec))
    for f in files + gb_files:
        remove(f)


def download_genbank(gid, outpath):
    Entrez.email='igor.bioinfo@gmail.com'
    handle = Entrez.efetch(db='nuccore', id=gid, rettype='gb', retmode='text')
    with open(os_join(outpath, gid + '.gb'), 'w') as local_file:
        local_file.write(handle.read())
    handle.close()

def iterate_tree(tree, taxid, taxid_to_file):
    target_mitogenomes = set()
    for node in PreOrderIter(tree[taxid]):
        if node.name in taxid_to_file:
            files = taxid_to_file[node.name]
            for f in files:
                target_mitogenomes.add(f)
    return target_mitogenomes

def get_taxid(term):
    Entrez.email='igor.bioinfo@gmail.com'
    handle = Entrez.esearch('taxonomy', term=term)#rettype='uilist', retmode='text', id=o)
    result = Entrez.read(handle)
    handle.close()
    if len(result['IdList']) > 1:
        print(f'more than one taxid: {result["IdList"]}')
    if not result['IdList']:
        print('not found')
    return result['IdList'][0]

if __name__ == '__main__':
    main()
