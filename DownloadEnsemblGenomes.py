#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

# This script is used to download genomes from Ensembl Bacteria (http://bacteria.ensembl.org/)

import ijson
import ftplib
import psycopg2
import datetime
import argparse, os, errno, sys, subprocess
from urllib import urlopen

def parse_args():
	parser = argparse.ArgumentParser(description='''
A script for accessing the current release of Ensembl Bacteria and downloading
complete genomes.  Will download nucleotide, amino acid, and CDS information, as
well as a metadata table for SQL.
	''')
	parser.add_argument('outdir', type=str,help='directory to download genomes to and name of psql library')
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

	p = open(os.path.join(outdir, "downloaded_genomes.txt"),'w')
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
			p.write("\t".join([str(x) for x in thisline])+"\n")

		o.write("\t".join([str(x) for x in thisline])+"\n")

		if count % 100 == 0:
			print count, "JSON records parsed."

	print len(finished_genomes), "of these are finished genomes."
	o.close()
	p.close()

	ens.close()
	return finished_genomes,j

def setupdirs(outdir):
	try:
		os.makedirs(outdir)
	except OSError as exc:
		if exc.errno == errno.EEXIST:
			print "Database folder exists:", outdir
			print "Exiting to prevent overwriting..."
			sys.exit()

	for f in ["cds","pep","dna","ncrna", "gff3", "genbank"]:
		try:
			os.makedirs(os.path.join(os.path.join(outdir,f)))
		except OSError:
			print "Subfolder exists:", outdir

	print outdir
	return

def get_fasta_files(fg, outdir):
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	print "Downloading fasta files..."
	count = 0
	for f in fg:
		ens.cwd("/pub/bacteria/current/fasta/{}/{}/dna".format("_".join(fg[f]["dbname"].split("_")[0:3]),fg[f]["species"]))
		for filepath in ens.nlst():
			if filepath.endswith(".dna.toplevel.fa.gz"):
				download_and_unzip(ens,filepath,os.path.join(outdir,"dna",fg[f]["species"]+".dna.fa.gz"))
		ens.cwd("../pep")
		for filepath in ens.nlst():
			if filepath.endswith(".pep.all.fa.gz"):
				download_and_unzip(ens,filepath,os.path.join(outdir,"pep",fg[f]["species"]+".pep.fa.gz"))
		ens.cwd("../cds")
		for filepath in ens.nlst():
			if filepath.endswith(".cds.all.fa.gz"):
				download_and_unzip(ens,filepath,os.path.join(outdir,"cds",fg[f]["species"]+".cds.fa.gz"))
		ens.cwd("../ncrna")
		for filepath in ens.nlst():
			if filepath.endswith(".ncrna.fa.gz"):
				download_and_unzip(ens,filepath,os.path.join(outdir,"ncrna",fg[f]["species"]+".ncrna.fa.gz"))
		count += 1
		if count % 10 == 0:
			print count, "files downloaded."
	ens.close()

def download_and_unzip(ftp,f,outfile):
	o = open(outfile,'wb')
	ftp.retrbinary("RETR " + f, o.write)
	o.close()
	cmds = ["gunzip",outfile]
	proc = subprocess.Popen(cmds)
	proc.wait()
	return

def get_gff_files(fg, outdir):
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	print "Downloading GFF3 files..."
	count = 0
	for f in fg:
		ens.cwd("/pub/bacteria/current/gff3/{}/{}".format("_".join(fg[f]["dbname"].split("_")[0:3]),fg[f]["species"]))
		for filepath in ens.nlst():
			fields = filepath.split(".")
			if ".".join(fields[3:]) == ("gff3.gz"):
				download_and_unzip(ens,filepath,os.path.join(outdir,"gff3",fg[f]["species"]+".gff3.gz"))
		count += 1
		if count % 10 == 0:
			print count, "files downloaded."
	ens.close()

def setup_psqldb():
	con = psycopg2.connect(user='ryan', dbname="postgres", host='localhost', password='')
	con.set_isolation_level(0)
	cur = con.cursor()
	cur.execute("CREATE DATABASE genomedb")
	cur.close()
	con.close()

	con = psycopg2.connect(user='ryan', dbname="genomedb", host='localhost', password='')
	cur = con.cursor()
	cur.execute("""CREATE TABLE genome_metadata (id serial PRIMARY KEY, assembly_id varchar, source varchar, base_count int,
		species varchar, taxonomy_id int, contigs int, protein_coding_genes int, date_added date, date_modified date);""")
	con.commit()
	cur.close()
	con.close()
	return

def populate_psqldb(fg, EV):
	print "Updating SQL database..."
	con = psycopg2.connect(user='ryan', dbname="genomedb", host='localhost', password='')
	cur = con.cursor()
	for f in fg:
		sql = """INSERT INTO genome_metadata (assembly_id, base_count, species, taxonomy_id, contigs,
			protein_coding_genes, source, date_added, date_modified) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s);"""
		vals = [f]
		for key in ["base_count", "species", "taxonomy_id", "contigs","ngenes"]:
			vals.append(fg[f][key])
		vals.append("ensembl-"+EV)
		vals.append(datetime.datetime.now())
		vals.append(datetime.datetime.now())
		cur.execute(sql, vals)
	con.commit()
	cur.close()
	con.close()

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)

	# this step will fail if the genomeDB psql database already exists - this is intentional
	setup_psqldb()

	setupdirs(outdir)
	finished_genomes, ENSEMBL_VERSION = parse_json(outdir)
	get_fasta_files(finished_genomes, outdir)
	get_gff_files(finished_genomes, outdir)
	populate_psqldb(finished_genomes,ENSEMBL_VERSION)

if __name__ == '__main__':
	main()
