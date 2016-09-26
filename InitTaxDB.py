#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, psycopg2

def parse_args():
	parser = argparse.ArgumentParser(description='''
Sets up the taxonomy table.
	''')
	parser.add_argument('sql_user',type=str, help='user for sql')
	parser.add_argument('sql_db',type=str,help='name of sql database')
	parser.add_argument('sql_host', type=str, help='host for sql: "localhost" for MacOSX psql install, "172.18.0.71" for bugaboo')
	return parser.parse_args()

def setup_psqldb(sql_user,sql_db,sql_host):
	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()
	cur.execute("""ALTER TABLE genome_metadata ADD CONSTRAINT unique_species_name UNIQUE (species);""")
	con.commit()
	cur.execute("""CREATE TABLE taxonomy (id serial PRIMARY KEY, species_id varchar references genome_metadata(species), taxonomy_id int, superkingdom varchar, phylum varchar, class varchar, __order__ varchar, family varchar, genus varchar, species varchar, date_added date, date_modified date)""")
	cur.execute("""ALTER TABLE taxonomy ADD CONSTRAINT unique_species_id UNIQUE (species_id);""")
	con.commit()
	cur.close()
	con.close()
	return

def main():
	args = parse_args()
	sql_user = args.sql_user
	sql_db = args.sql_db
	sql_host = args.sql_host
	setup_psqldb(sql_user,sql_db,sql_host)

if __name__ == '__main__':
	main()
