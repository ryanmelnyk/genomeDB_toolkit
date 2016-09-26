#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, psycopg2, os
from Bio import Entrez
import datetime


def parse_args():
	parser = argparse.ArgumentParser(description='''
Populate the taxonomy psql table with taxonomy information extracted from NCBI Taxonomy.
	''')
	parser.add_argument('outdir', type=str,help='directory to download genomes')
	parser.add_argument('sql_user',type=str, help='user for sql')
	parser.add_argument('sql_db',type=str,help='name of sql database')
	parser.add_argument('sql_host', type=str, help='host for sql: "localhost" for MacOSX psql install, "172.18.0.71" for bugaboo')
	return parser.parse_args()

def get_genomedb_data(strains,sql_user,sql_db,sql_host):

	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()
	cur.execute("SELECT * FROM genome_metadata")
	records = cur.fetchall()
	taxdata = {}
	for r in records:
		if r[4] not in strains:
			taxdata[r[4]] = r[5]

	cur.close()
	con.close()
	return taxdata

def fetch_tax(taxdata,sql_user,sql_db,sql_host):
	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()
	Entrez.email = "schmelnyk@gmail.com"
	sql = """INSERT INTO taxonomy (species_id, taxonomy_id, superkingdom, phylum, class,
			__order__, family, genus, species, date_added, date_modified) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);"""
	count = len(taxdata)
	print "Extracting", str(count), "taxonomy records from Entrez-NCBI..."
	for t in taxdata:
		values = [t,taxdata[t],None,None,None,None,None,None,None,datetime.datetime.now(),datetime.datetime.now()]
		records = Entrez.parse(Entrez.efetch(db="taxonomy",id=str(taxdata[t]),retmode="xml"))
		for r in records:
			for l in r["LineageEx"]:
				if l["Rank"] == "superkingdom":
					values[2] = l["ScientificName"]
				elif l["Rank"] == "phylum":
					values[3] = l["ScientificName"]
				elif l["Rank"] == "class":
					values[4] = l["ScientificName"]
				elif l["Rank"] == "order":
					values[5] = l["ScientificName"]
				elif l["Rank"] == "family":
					values[6] = l["ScientificName"]
				elif l["Rank"] == "genus":
					values[7] = l["ScientificName"]
				elif l["Rank"] == "species":
					values[8] = l["ScientificName"]
				else:
					pass
		count -= 1
		if count % 100 == 0:
			print count, "records remaining..."

		cur.execute(sql,values)
		con.commit()
	cur.close()
	con.close()
	print "Done!"

	return

def dump_flat_file(outdir,sql_user,sql_db,sql_host):
	print "Dumping SQL database to taxonomy.txt..."
	o = open(os.path.join(outdir, "taxonomy.txt"),'w')
	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()

	cur.execute("SELECT * FROM taxonomy")
	o.write("\t".join([desc[0] for desc in cur.description])+"\n")
	records = cur.fetchall()
	for r in records:
		o.write("\t".join([str(x) for x in r])+"\n")

	cur.close()
	con.close()
	o.close()
	return

def query_sql(sql_user,sql_db,sql_host):
	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()
	cur.execute("SELECT * FROM taxonomy")
	records = cur.fetchall()
	strains = []
	for r in records:
		strains.append(r[1])

	cur.close()
	con.close()
	return strains

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)
	sql_user = args.sql_user
	sql_db = args.sql_db
	sql_host = args.sql_host

	strains = query_sql(sql_user,sql_db,sql_host)
	taxdata = get_genomedb_data(strains, sql_user,sql_db,sql_host)
	fetch_tax(taxdata,sql_user,sql_db,sql_host)
	dump_flat_file(outdir,sql_user,sql_db,sql_host)

if __name__ == '__main__':
	main()
