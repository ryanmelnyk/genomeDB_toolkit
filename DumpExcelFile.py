#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os, psycopg2
import xlsxwriter

def parse_args():
	parser = argparse.ArgumentParser(description='''
Dump taxonomy and genome metadata to a single excel file to add to ownCloud.
	''')
	parser.add_argument('outdir', type=str,help='directory to download genomes')
	parser.add_argument('sql_user',type=str, help='user for sql')
	parser.add_argument('sql_db',type=str,help='name of sql database')
	parser.add_argument('sql_host', type=str, help='host for sql: "localhost" for MacOSX psql install, "172.18.0.71" for bugaboo')
	return parser.parse_args()

def get_sql_data(sql_user,sql_db,sql_host):
	dat = {}

	con = psycopg2.connect(user=sql_user, dbname=sql_db, host=sql_host, password='')
	cur = con.cursor()
	cur.execute("SELECT * FROM genome_metadata")
	records = cur.fetchall()
	for r in records:
		dat[r[4]] = [r[1],r[2],r[3],r[5],r[6],r[7]]

	cur.execute("SELECT * FROM taxonomy")
	records = cur.fetchall()
	for r in records:
		for i in range(3,10):
			dat[r[1]].append(r[i])


	cur.close()
	con.close()
	return dat

def write_xl(outdir, dat):
	workbook = xlsxwriter.Workbook(os.path.join(outdir, 'genomedb.xlsx'))
	worksheet = workbook.add_worksheet()

	row = 0

	for d in dat:
		col = 0
		worksheet.write(row, col, d)
		col += 1
		for s in dat[d]:
			worksheet.write(row, col, s)
			col += 1
		row += 1
	return

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)
	sql_user = args.sql_user
	sql_db = args.sql_db
	sql_host = args.sql_host

	dat = get_sql_data(sql_user,sql_db,sql_host)
	write_xl(outdir,dat)

if __name__ == '__main__':
	main()
