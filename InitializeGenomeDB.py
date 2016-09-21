#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, psycopg2

def parse_args():
	parser = argparse.ArgumentParser(description='''
Sets up the genomedb PostgreSQL database and the genome_metadata table.
	''')
	return parser.parse_args()

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

def main():
	setup_psqldb()
	pass

if __name__ == '__main__':
	main()
