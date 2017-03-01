#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os
import ftplib
import subprocess

def parse_args():
	parser = argparse.ArgumentParser(description='''
Given a list of strain names, accesses the Ensembl FTP server to download genbank files.
	''')
	parser.add_argument('strains', type=str,help='path to text file with strain names')
	parser.add_argument('genomedb',type=str,help='path to genomedb folder')
	return parser.parse_args()

def download(strains, genomedb):
	ens = ftplib.FTP('ftp.ensemblgenomes.org')
	ens.login()

	for line in open(os.path.join(genomedb,"release-32.txt")):
		vals = line.rstrip().split("\t")
		if vals[6] in strains:
			wd = 'pub/release-32/bacteria/genbank/{}_collection/{}'.format("_".join(vals[5].split("_")[0:2]),vals[6])
			print wd
			ens.cwd(wd)

			for filepath in ens.nlst():
				if filepath.endswith(".dat.gz"):
					o = open(strains[0]+".gbk.gz",'wb')
					ens.retrbinary("RETR " + filepath, o.write)
					o.close()
					cmds = ["gunzip",strains[0]+".gbk.gz"]
					proc = subprocess.Popen(cmds)
					proc.wait()
		else:
			pass

	ens.close()

def main():
	args = parse_args()
	strains = [line.rstrip() for line in open(os.path.abspath(args.strains),'r')]
	print strains
	genomedb = os.path.abspath(args.genomedb)


	download(strains,genomedb)

if __name__ == '__main__':
	main()
