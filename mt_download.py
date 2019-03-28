from __future__ import print_function
from future import standard_library

standard_library.install_aliases()

import argparse
from os import getcwd
from os.path import join
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError
from Bio import Entrez, SeqIO

def argument_parser():

    default_out = getcwd() + '/'
    parser = argparse.ArgumentParser(description = 'Downloads mitochondrial genomes.',\
                                     argument_default = None, fromfile_prefix_chars = '@')
    parser.add_argument('-t', '--taxonomy', nargs = '?', type = str, required = True,\
                        dest = 'taxa', help = 'Taxonomy id.')
    parser.add_argument('-o', '--outpath', nargs = '?', type = str, default = default_out,\
                        dest = 'outpath', help = 'Path were the mitogenomes\' genbank files will be saved. (default: %(default)s)')
    args = parser.parse_args().__dict__
    return args

def download_mt(taxa, outpath):
  '''Downloads all mitochondria complete genomes from the organisms in "taxa".
  Saves them in genbank format inside "outpath" folder.'''
  
  #https://www.ncbi.nlm.nih.gov/genomes/OrganelleResource.cgi?taxid=9223&cmd=download1  
  data = dict(taxid=taxa, cmd='download1')
  try:
    accesion_ids = urlopen(Request("https://www.ncbi.nlm.nih.gov/genomes/OrganelleResource.cgi", urlencode(data)))
  except HTTPError as e:
    print(e.read())
  Entrez.email = 'igorrcosta@hotmail.com'
  accesion_list = accesion_ids.read().split()
  n = len(accesion_list)
  if n:
      print(n, 'mitocondrial genomes found.')
  else:
      print('No mitocondrial genomes found in this taxon id.') 
  for i in accesion_list:
      # Use this? https://www.ncbi.nlm.nih.gov/genome/browse#!/organelles/
      print('Downloading genome id:', i)
      handle = Entrez.efetch(db='nuccore', id=i, rettype='gb', retmode='text')
      try:
        with open(join(outpath, i + '.gbk'), 'w') as local_file:
            SeqIO.write(SeqIO.read(handle, "genbank"), local_file, "genbank")
            handle.close()
      except:
        print('Cannot open file: Check if you have writing permissions.')
        raise
if __name__ == '__main__':
    args = argument_parser()
    download_mt(args['taxa'], args['outpath'])
