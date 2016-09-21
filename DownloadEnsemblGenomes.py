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
	parser.add_argument('outdir', type=str,help='directory to download genomes')
	parser.add_argument('sql_user',type=str, help='user for sql')
	parser.add_argument('sql_db',type=str,help='name of sql database')
	parser.add_argument('sql_host', type=str, help='host for sql: "localhost" for MacOSX psql install, "172.18.0.71" for bugaboo')
	parser.add_argument('--names', type=str, help='comma-separated list of keywords to find in name of strain to download draft genomes. i.e. "pseudomonas,salmonella,syringae,K12". Only complete/exact matches will be downloaded and spelling counts!')
	return parser.parse_args()

def parse_json(outdir, assemblies, names):
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
			if js["assembly_id"] not in assemblies:
				finished_genomes[js["assembly_id"]] = {}
				for feat in ['assembly_level','base_count','name', 'strain', 'dbname','species','taxonomy_id']:
					finished_genomes[js["assembly_id"]][feat] = js[feat]
				finished_genomes[js["assembly_id"]]['contigs'] = len(js["sequences"])
				finished_genomes[js["assembly_id"]]['ngenes'] = js["annotations"]["nProteinCoding"]
		else:
			fields = js['species'].split("_")
			for n in names:
				if n in fields:
					if js["assembly_id"] not in assemblies:
						finished_genomes[js["assembly_id"]] = {}
						for feat in ['assembly_level','base_count','name', 'strain', 'dbname','species','taxonomy_id']:
							finished_genomes[js["assembly_id"]][feat] = js[feat]
						finished_genomes[js["assembly_id"]]['contigs'] = len(js["sequences"])
						finished_genomes[js["assembly_id"]]['ngenes'] = js["annotations"]["nProteinCoding"]

		o.write("\t".join([str(x) for x in thisline])+"\n")

		if count % 500 == 0:
			print count, "JSON records parsed."

	print count, "total JSON records parsed."
	print len(assemblies), "found in", os.path.basename(outdir)+"."
	print len(finished_genomes), "remaining to download."
	o.close()
	ens.close()
	return finished_genomes,j

def setupdirs(outdir):
	try:
		os.makedirs(outdir)
	except OSError as exc:
		if exc.errno == errno.EEXIST:
			print "Database folder exists:", outdir

	for f in ["cds","pep","dna","ncrna", "gff3", "genbank"]:
		try:
			os.makedirs(os.path.join(os.path.join(outdir,f)))
		except OSError:
			print "Subfolder exists:", outdir

	print outdir
	return

def get_files(fg, outdir, EV, sql_user, sql_db, sql_host):
	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()
	print "Downloading genome files..."
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
		ens.cwd("/pub/bacteria/current/gff3/{}/{}".format("_".join(fg[f]["dbname"].split("_")[0:3]),fg[f]["species"]))
		for filepath in ens.nlst():
			fields = filepath.split(".")
			if ".".join(fields[3:]) == ("gff3.gz"):
				download_and_unzip(ens,filepath,os.path.join(outdir,"gff3",fg[f]["species"]+".gff3.gz"))
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
		count += 1
		print fg[f]["species"], "processed.", count, "files downloaded."

	ens.close()
	cur.close()
	con.close()
	return

def download_and_unzip(ftp,f,outfile):
	o = open(outfile,'wb')
	ftp.retrbinary("RETR " + f, o.write)
	o.close()
	cmds = ["gunzip",outfile]
	proc = subprocess.Popen(cmds)
	proc.wait()
	return

def query_sql(sql_user,sql_db,sql_host):
	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()
	cur.execute("SELECT * FROM genome_metadata")
	records = cur.fetchall()
	assemblies = []
	for r in records:
		assemblies.append(r[1])

	cur.close()
	con.close()
	return assemblies

def dump_flat_file(outdir,sql_user,sql_db,sql_host):
	print "Dumping SQL database to genome_metadata.txt..."
	o = open(os.path.join(outdir, "genome_metadata.txt"),'w')
	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()

	cur.execute("SELECT * FROM genome_metadata")
	o.write("\t".join([desc[0] for desc in cur.description])+"\n")
	records = cur.fetchall()
	for r in records:
		o.write("\t".join([str(x) for x in r])+"\n")

	cur.close()
	con.close()
	o.close()
	return

def main():
	args = parse_args()

	outdir = os.path.abspath(args.outdir)
	sql_user = args.sql_user
	sql_db = args.sql_db
	sql_host = args.sql_host
	if args.names == None:
		names = []
	else:
		names = args.names.split(",")
	for n in names:
		print n

	setupdirs(outdir)

	assemblies = query_sql(sql_user,sql_db,sql_host)

	finished_genomes, ENSEMBL_VERSION = parse_json(outdir,assemblies,names)
	get_files(finished_genomes, outdir, ENSEMBL_VERSION, sql_user, sql_db, sql_host)
	dump_flat_file(outdir,sql_user,sql_db,sql_host)

if __name__ == '__main__':
	main()
