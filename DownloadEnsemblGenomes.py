#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

# This script is used to download genomes from Ensembl Bacteria (http://bacteria.ensembl.org/)

import ijson
import ftplib
from urllib import urlopen

def parse_args():
	parser = argparse.ArgumentParser(description='''
A script for accessing the current release of Ensembl Bacteria and downloading
complete genomes.  Will download nucleotide, amino acid, and CDS information, as
well as a metadata table.
	''')
	parser.add_argument('outdir', type=str,help='directory to download genomes to')
	return parser.parse_args()

def parse_json():
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	ens.cwd('pub/bacteria/current')

	# for line in ens.retrlines('LIST'):
	# 	print line,
		# print line.rstrip().split()[-1]

	for j in ens.pwd().split("/"):
		if j.startswith("release"):
			print "Current release of EnsemblBacteria:", j

	items = ijson.items(urlopen("ftp://ftp.ensemblgenomes.org/{}/species_metadata_EnsemblBacteria.json".format(ens.pwd())),'item')
	count = 0
	finished_genomes = {}
	o = open("test.txt",'w')
	for js in items:
		count += 1

		thisline = []
		finished_genomes[js["assembly_id"]] = {}
		thisline.append(js["assembly_id"])
		for feat in ['assembly_level','base_count','name', 'strain', 'dbname','taxonomy_id']:
			finished_genomes[js["assembly_id"]][feat] = js[feat]
			thisline.append(js[feat])
		finished_genomes[js["assembly_id"]]['contigs'] = len(js["sequences"])
		thisline.append(str(len(js["sequences"])))
		finished_genomes[js["assembly_id"]]['ngenes'] = js["annotations"]["nProteinCoding"]
		thisline.append(js["annotations"]["nProteinCoding"])

		o.write("\t".join([str(x) for x in thisline])+"\n")

		if count % 100 == 0:
			print count, "JSON records parsed."

	print len(finished_genomes), "of these are finished genomes."


	return finished_genomes

def main():
	finished_genomes = parse_json()



if __name__ == '__main__':
	main()
