#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

# This script is used to download genomes from Ensembl Bacteria (http://bacteria.ensembl.org/)

import ijson
import ftplib
import argparse, os, errno, sys
from urllib import urlopen

def parse_args():
	parser = argparse.ArgumentParser(description='''
A script for accessing the current release of Ensembl Bacteria and downloading
complete genomes.  Will download nucleotide, amino acid, and CDS information, as
well as a metadata table.
	''')
	parser.add_argument('outdir', type=str,help='directory to download genomes to')
	return parser.parse_args()

def parse_json(outdir):
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	ens.cwd('pub/bacteria/current')

	# for line in ens.retrlines('LIST'):
	# 	print line,
		# print line.rstrip().split()[-1]

	for j in ens.pwd().split("/"):
		if j.startswith("release"):
			print "Current release of EnsemblBacteria:", j
			o = open(os.path.join(outdir,"{}.txt".format(j)),'w')
			break

	fields = ["assembly_id",'assembly_level','base_count','name', 'strain', 'dbname','species','taxonomy_id','contigs','protein_coding_genes']
	o.write("\t".join(fields)+"\n")

	items = ijson.items(urlopen("ftp://ftp.ensemblgenomes.org/{}/species_metadata_EnsemblBacteria.json".format(ens.pwd())),'item')
	count = 0
	finished_genomes = {}

	for js in items:
		count += 1

		thisline = []
		thisline.append(js["assembly_id"])
		for feat in ['assembly_level','base_count','name', 'strain', 'dbname','species','taxonomy_id']:
			thisline.append(js[feat])
		thisline.append(str(len(js["sequences"])))
		thisline.append(js["annotations"]["nProteinCoding"])

		if js["assembly_level"] == "chromosome":
			finished_genomes[js["assembly_id"]] = {}
			for feat in ['assembly_level','base_count','name', 'strain', 'dbname','species','taxonomy_id']:
				finished_genomes[js["assembly_id"]][feat] = js[feat]
			finished_genomes[js["assembly_id"]]['contigs'] = len(js["sequences"])
			finished_genomes[js["assembly_id"]]['ngenes'] = js["annotations"]["nProteinCoding"]

		o.write("\t".join([str(x) for x in thisline])+"\n")

		if count % 100 == 0:
			print count, "JSON records parsed."

		#### remove! only for testing purposes
		if count > 500:
			break

	print len(finished_genomes), "of these are finished genomes."

	ens.close()
	return finished_genomes

def setupdirs(outdir):
	try:
		os.makedirs(outdir)
	except OSError as exc:
		if exc.errno == errno.EEXIST:
			print "Database folder exists:", outdir
			print "Exiting to prevent overwriting..."
			sys.exit()

	for f in ["nuc","prot","dna","gff"]:
		try:
			os.makedirs(os.path.join(os.path.join(outdir,f)))
		except OSError:
			print "Subfolder exists:", outdir

	print outdir
	return

def extract_genomes(fg, outdir):
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	for f in fg:
		print f
		print fg[f]["name"]
		print fg[f]["dbname"]
		ens.cwd("/pub/bacteria/current/fasta/{}/{}/dna".format("_".join(fg[f]["dbname"].split("_")[0:3]),fg[f]["species"]))
		print ens.pwd()
		print ens.nlst()
		for filepath in ens.nlst():
			if filepath.endswith(".dna.toplevel.fa.gz"):
				o = open(os.path.join(outdir,"dna",fg[f]["species"]+".fa.gz"),'wb')
				ens.retrbinary("RETR " + filepath, o.write)


def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)
	setupdirs(outdir)
	finished_genomes = parse_json(outdir)
	extract_genomes(finished_genomes, outdir)


if __name__ == '__main__':
	main()
