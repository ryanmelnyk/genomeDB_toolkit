#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os, shutil

def parse_args():
	parser = argparse.ArgumentParser(description='''
Organize ensembl files for running phylophlan.
	''')
	parser.add_argument('path_to_genomedb', type=str,help='relative path to genomedb set up by DownloadEnsemblGenomes.py')
	parser.add_argument('path_to_phylophlan', type=str,help='relative path to phylophlan directory')
	return parser.parse_args()

def main():
	args = parse_args()
	gdb = os.path.abspath(args.path_to_genomedb)
	phylo = os.path.abspath(args.path_to_phylophlan)

	os.mkdir(os.path.join(phylo,"input","genomedb"))

	count = 0
	for f in os.listdir(os.path.join(gdb,"pep")):
		count += 1
		shutil.copy(os.path.join(gdb,"pep",f),os.path.join(phylo,"input","genomedb",f.split(".")[0]+".faa"))
		if count % 200 == 0:
			print count, "files copied..."
	pass

if __name__ == '__main__':
	main()
