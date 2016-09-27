#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os, psycopg2, shutil
import sys, datetime
from Bio import SeqIO, Entrez

def parse_args():
	parser = argparse.ArgumentParser(description='''
This script will add a local genome annotated with Prokka to the genomedb.""
	''')
	parser.add_argument('prokka', type=str,help='relative path to prokka directory')
	parser.add_argument('outdir', type=str,help='directory to download genomes')
	parser.add_argument('species_id', type=str,help='species_id for strain identification (Must be unique in genomedb!)')
	parser.add_argument('tax_id', type=str, help='taxonomy id for NCBI lookup. (Pseudomonas = "286")')
	parser.add_argument('sql_user',type=str, help='user for sql')
	parser.add_argument('sql_db',type=str,help='name of sql database')
	parser.add_argument('sql_host', type=str, help='host for sql: "localhost" for MacOSX psql install, "172.18.0.71" for bugaboo')
	return parser.parse_args()

def check_unique(species_id,sql_user,sql_db,sql_host):
	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()
	cur.execute("SELECT * FROM genome_metadata")
	records = cur.fetchall()
	strains = []
	for r in records:
		strains.append(r[4])

	if species_id not in strains:
		print "Species ID is unique! Moving on..."
	else:
		print "Species ID is not unique. Select a new ID."
		print "Exiting script..."
		sys.exit()

	cur.close()
	con.close()
	return


def copy_files(outdir,prokka,species_id):
	CDS = []
	RNA = []
	stats = {}
	print "Copying files..."
	for f in os.listdir(prokka):
		if f.endswith(".faa"):
			shutil.copy(os.path.join(prokka,f),os.path.join(outdir,"pep","{}.pep.fa".format(species_id)))
		elif f.endswith(".fna"):
			basecount = 0
			contigcount = 0
			for seq in SeqIO.parse(open(os.path.join(prokka,f),'r'),'fasta'):
				basecount += len(str(seq.seq))
				contigcount += 1
			stats['basecount'] = basecount
			stats['contigcount'] = contigcount
			shutil.copy(os.path.join(prokka,f),os.path.join(outdir,"dna","{}.dna.fa".format(species_id)))
		elif f.endswith(".gff"):
			for line in open(os.path.join(prokka,f)):
				if line.startswith("#"):
					if line.rstrip().endswith("FASTA"):
						break
					else:
						pass
				else:
					vals = line.rstrip().split("\t")
					geneID = vals[8].split(";")[0].split("=")[1]
					if vals[2] == "CDS":
						CDS.append(geneID)
					elif vals[2].endswith("RNA"):
						RNA.append(geneID)
					else:
						pass
			shutil.copy(os.path.join(prokka,f),os.path.join(outdir,"gff3","{}.gff3".format(species_id)))
		else:
			pass
	for f in os.listdir(prokka):
		if f.endswith(".ffn"):
			count = 0
			CDSfile = open(os.path.join(outdir,"cds","{}.cds.fa".format(species_id)),'w')
			RNAfile = open(os.path.join(outdir,"ncrna","{}.ncrna.fa".format(species_id)),'w')
			for seq	in SeqIO.parse(open(os.path.join(prokka,f),'r'),'fasta'):
				if seq.id in CDS:
					SeqIO.write(seq,CDSfile,'fasta')
					count += 1
				elif seq.id in RNA:
					SeqIO.write(seq,RNAfile,'fasta')
				else:
					pass
			stats['ngenes'] = count
			CDSfile.close()
			RNAfile.close()
		else:
			pass
	return stats

def update_SQLtables(species_id,tax_id,sql_user,sql_db,sql_host,stats):
	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()
	sql = """INSERT INTO genome_metadata (assembly_id, base_count, species, taxonomy_id, contigs,
		protein_coding_genes, source, date_added, date_modified) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s);"""
	vals = [species_id+"_v1",stats['basecount'],species_id,tax_id,stats['contigcount'],stats['ngenes'],"prokka_in_house", datetime.datetime.now(),datetime.datetime.now()]
	cur.execute(sql, vals)

	sql = """INSERT INTO taxonomy (species_id, taxonomy_id, superkingdom, phylum, class,
			__order__, family, genus, species, date_added, date_modified) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);"""
	Entrez.email = "schmelnyk@gmail.com"
	tax_vals = [species_id,tax_id,None,None,None,None,None,None,None,datetime.datetime.now(),datetime.datetime.now()]
	records = Entrez.parse(Entrez.efetch(db="taxonomy",id=tax_id,retmode="xml"))
	for r in records:
		for l in r["LineageEx"]:
			if l["Rank"] == "superkingdom":
				tax_vals[2] = l["ScientificName"]
			elif l["Rank"] == "phylum":
				tax_vals[3] = l["ScientificName"]
			elif l["Rank"] == "class":
				tax_vals[4] = l["ScientificName"]
			elif l["Rank"] == "order":
				tax_vals[5] = l["ScientificName"]
			elif l["Rank"] == "family":
				tax_vals[6] = l["ScientificName"]
			elif l["Rank"] == "genus":
				tax_vals[7] = l["ScientificName"]
			elif l["Rank"] == "species":
				tax_vals[8] = l["ScientificName"]
			else:
				pass
	cur.execute(sql, tax_vals)
	con.commit()

	return

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
	o.close()

	print "Dumping SQL database to taxonomy.txt..."
	o = open(os.path.join(outdir, "taxonomy.txt"),'w')
	cur.execute("SELECT * FROM taxonomy")
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
	prokka = os.path.abspath(args.prokka)
	outdir = os.path.abspath(args.outdir)
	species_id = args.species_id
	tax_id = args.tax_id
	sql_user = args.sql_user
	sql_db = args.sql_db
	sql_host = args.sql_host

	check_unique(species_id,sql_user,sql_db,sql_host)
	stats = copy_files(outdir,prokka,species_id)
	update_SQLtables(species_id,tax_id,sql_user,sql_db,sql_host,stats)
	dump_flat_file(outdir,sql_user,sql_db,sql_host)

if __name__ == '__main__':
	main()
