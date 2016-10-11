#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, psycopg2, os
from Bio import Entrez
import datetime


def parse_args():
	parser = argparse.ArgumentParser(description='''
Change taxonomy in the PSQL database and redump taxonomy.txt
	''')
	parser.add_argument('outdir', type=str,help='directory to download genomes')
	parser.add_argument('sql_user',type=str, help='user for sql')
	parser.add_argument('sql_db',type=str,help='name of sql database')
	parser.add_argument('sql_host', type=str, help='host for sql: "localhost" for MacOSX psql install, "172.18.0.71" for bugaboo')
	parser.add_argument('tax_changes',type=str, help='path to file containing species_id and new taxonomy id')
	return parser.parse_args()

def query_sql(sql_user,sql_db,sql_host,tax):
	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()

	for t in tax:
		print "Fetching {}...".format(t)
		cur.execute("SELECT * FROM taxonomy WHERE species_id = %s",[t])
		for old_vals in cur.fetchall():
			print old_vals
			for r in Entrez.parse(Entrez.efetch(db="taxonomy",id=str(tax[t]),retmode="xml")):
				changed = False
				for l in r["LineageEx"]:
					if l["Rank"] == "superkingdom":
						if l["ScientificName"] != old_vals[3]:
							cur.execute("UPDATE taxonomy SET superkingdom = %s WHERE id = %s;",[l["ScientificName"],old_vals[0]])
							changed = True
					elif l["Rank"] == "phylum":
						if l["ScientificName"] != old_vals[4]:
							cur.execute("UPDATE taxonomy SET phylum = %s WHERE id = %s;",[l["ScientificName"],old_vals[0]])
							changed = True
					elif l["Rank"] == "class":
						if l["ScientificName"] != old_vals[5]:
							cur.execute("UPDATE taxonomy SET class = %s WHERE id = %s;",[l["ScientificName"],old_vals[0]])
							changed = True
					elif l["Rank"] == "order":
						if l["ScientificName"] != old_vals[6]:
							cur.execute("UPDATE taxonomy SET __order__ = %s WHERE id = %s;",[l["ScientificName"],old_vals[0]])
							changed = True
					elif l["Rank"] == "family":
						if l["ScientificName"] != old_vals[7]:
							cur.execute("UPDATE taxonomy SET family = %s WHERE id = %s;",[l["ScientificName"],old_vals[0]])
							changed = True
					elif l["Rank"] == "genus":
						if l["ScientificName"] != old_vals[8]:
							cur.execute("UPDATE taxonomy SET genus = %s WHERE id = %s;",[l["ScientificName"],old_vals[0]])
							changed = True
					elif l["Rank"] == "species":
						if l["ScientificName"] != old_vals[9]:
							cur.execute("UPDATE taxonomy SET species = %s WHERE id = %s;",[l["ScientificName"],old_vals[0]])
							changed = True
					else:
						pass
				if changed:
					cur.execute("UPDATE taxonomy SET date_modified = %s WHERE id = %s;",[datetime.datetime.now(),old_vals[0]])
	con.commit()
	cur.close()
	con.close()
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

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)
	sql_user = args.sql_user
	sql_db = args.sql_db
	sql_host = args.sql_host
	tax = {vals[0]:vals[1] for vals in [line.rstrip().split("\t") for line in open(os.path.abspath(args.tax_changes))]}

	query_sql(sql_user,sql_db,sql_host,tax)
	dump_flat_file(outdir,sql_user,sql_db,sql_host)

if __name__ == '__main__':
	main()
