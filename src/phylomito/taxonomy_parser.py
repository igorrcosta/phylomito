import re
import pickle
import gzip
from urllib.request import urlopen
from os import listdir
from collections import defaultdict

from Bio import SeqIO
from anytree import Node

#https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER = 206
#https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.genomic.gbff.gz
#https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.genomic.gbff.gz
#https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

def parse_names(names_path):
    taxid_name = {}
    small_taxid_name = {}
    duplicates = {}
    auth_reg = re.compile('\([A-Z][a-z]*, \d{4}\)$')
    synonym = {}
    taxid_dict = {}
    assert auth_reg.search('Achatina fulica (Ferussac, 1821)'), 'Regex not working'
    with open(names_path) as infile:
        for l in infile:
            taxid = l.split('|')[0].strip()
            name = l.split('|')[1].strip()
            if name == 'environmental samples':
                continue
            status = l.split('|')[3].strip()
            if status in ['scientific name', 'synonym', 'equivalent name', 'includes']:
                if auth_reg.search(name):
                    #print('auth in name:', name)
                    name = '('.join(name.split('(')[:-1]).strip()
                if name not in taxid_name:
                    taxid_name[name] = taxid
                    if len(name.split()) == 1:
                        taxid_dict[taxid] = name
                    if len(name.split()) > 2:
                        small_name = ' '.join(name.split()[:-1])
                        if small_name not in taxid_name:
                            small_taxid_name[small_name] = taxid
                else:
                    if taxid != taxid_name[name]:
                        duplicates[name] = taxid
                    #print(l.strip())
                    #print(taxid_name[name])
                if status == 'synonym' and len(name.split()) == 1 and name[0].isupper():
                    synonym[name] = taxid
    for k in small_taxid_name:
        if k not in taxid_name:
            taxid_name[k] = small_taxid_name[k]
    assert 'Achatina fulica' in taxid_name, 'not in name'
    return taxid_name, duplicates, synonym, taxid_dict

def parse_refseq(refseq):
    auth_reg = re.compile('\([A-Z][a-z]*, \d{4}\)$')
    taxid_name, duplicates, synonym, taxid_dict = parse_names('taxonomy/names.dmp')
    organism_taxid = {} 
    typo = {'Pseudosiderastrea tayami':'Pseudosiderastrea tayamai', 'Phylloscopus inornata':'Phylloscopus inornatus'}
    for mito_file in refseq:
        for record in SeqIO.parse(mito_file, 'genbank'):
            organism = record.annotations['organism']
            if organism in typo:
                organism = typo[organism]
            # Removes tags like (Smith, 1999) from the end of the organism name
            if auth_reg.search(organism):
                print(f'removing auth {organism}') 
                organism = '('.join(organism.split('(')[:-1]).strip()
            if organism not in taxid_name:
                # Found a synonym genera 
                if organism.split()[0] in synonym:
                    new_genera = synonym[organism.split()[0]]
                    organism = ' '.join([taxid_dict[new_genera]] + organism.split()[1:])
                else:
                    if organism.split()[0] not in taxid_name:
                        raise ValueError(f'{organism} not in name dict \n {record.annotations}')
                    else:
                        # Found only genera
                        organism = organism.split()[0]
            if organism in duplicates:
                print(f'{organism} in duplicates')
            if organism in typo: # Can be in typo after changing genera
                organism = typo[organism]
            if record.annotations['organism'] not in organism_taxid:
                organism_taxid[record.annotations['organism']] = taxid_name[organism] 
    with open('./data/organism_taxid', 'wb') as pickle_file:
        pickle.dump(organism_taxid, pickle_file)

def split_refseq(refseq):
    # Spliting as binary data so biopython doesn't clip the organism's name.
    for f in refseq:
        all_rec = SeqIO.index(f, 'genbank')
        for rec in all_rec:
            file_name = rec.split('.')[0] + '.gb'
            with open('refseq/' + file_name, 'wb') as outfile:
                outfile.write(all_rec.get_raw(rec))
        
def parse_taxonomy(taxonomy_dump, all_mitochondria_dump):
    with open(all_mitochondria_dump, 'rb') as pickle_dump:
        mitochondria_taxid = pickle.load(pickle_dump)
    all_mitochondria_taxid = mitochondria_taxid.values()
    
    taxonomy = {}
    fields = ['tax_id', 'parent', 'rank', 'embl', 'division', 'inherited',
              'genetic_code', 'inherited_gc', 'mitochondrial_gc', 'inherited_mgc',
              'genbank_hidden', 'subtree_hidden']
    with open(taxonomy_dump) as infile:
        for l in infile:
            parsed = {field:n.strip() for field, n in zip(fields, l.split('|'))}
            taxonomy[parsed['tax_id']] = parsed['parent']
    mito_tree = {}
    parents = set()
    for mitochondria in all_mitochondria_taxid:
        if mitochondria not in mito_tree:
            mito_tree[taxonomy[mitochondria]] = Node(taxonomy[mitochondria])
            mito_tree[mitochondria] = Node(mitochondria, parent=mito_tree[taxonomy[mitochondria]])
            parents.add(taxonomy[mitochondria])
    while parents:
        parent = parents.pop()
        if parent == taxonomy[parent]: #root:
            continue
        if taxonomy[parent] not in mito_tree:
            mito_tree[taxonomy[parent]] = Node(taxonomy[parent])
        mito_tree[parent].parent = mito_tree[taxonomy[parent]]
        parents.add(taxonomy[parent])
    with open('./data/mito_tree', 'wb') as outfile:
        pickle.dump(mito_tree, outfile)

def make_file_name_dict(organism_taxid_dump):
    file_name_dict = defaultdict(set)
    with open(organism_taxid_dump, 'rb') as pickle_dump:
        organism_taxid = pickle.load(pickle_dump)
    print(f'number of species: {len(organism_taxid)}')
    refseq_folder = 'refseq/'
    for f in listdir(refseq_folder):
        record = SeqIO.read(refseq_folder + f, 'genbank')
        organism = record.annotations['organism']
        try:
            file_name_dict[organism_taxid[organism]].add(f)
        except:
            print(organism, f)
            raise
    with open('./data/taxid_to_filenames', 'wb') as out:
        pickle.dump(file_name_dict, out)

def test_tree(pickle_dump):
    with open(pickle_dump, 'rb') as dmp_file:
        mito_tree = pickle.load(dmp_file)
    for k in mito_tree:
        if k != '1':
            assert mito_tree[k].parent is not None, k
        else:
            assert mito_tree[k].parent is None

def download_refseq(files):
    max_attempts = 5
    attempts = 0
    interval = 10
    base_url = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/'
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

if __name__ == '__main__':
    refseq = ['mitochondria/mitochondrion.1.genomic.gbff', 'mitochondria/mitochondrion.2.genomic.gbff']
    download_refseq(refseq)
    parse_refseq(refseq)
    infile = './taxonomy/nodes.dmp'
    pickle_dump = './data/organism_taxid'
    parse_taxonomy(infile, pickle_dump)
    mito_tree_dump = './data/mito_tree'
    test_tree(mito_tree_dump)
    split_refseq(refseq)
    pickle_dump = './data/organism_taxid'
    make_file_name_dict(pickle_dump)
    
