import re
import sys
import argparse

def parse_args():
	parser = argparse.ArgumentParser(description='De_dup vcf.')
	parser.add_argument('-i', dest="input", type = str, required = True)
	parser.add_argument('-o', dest="output", type = str, required = True)

	return parser.parse_args()

def main():
	args = parse_args()

	vcffile = args.input
	fwrite = args.output
	open(fwrite,'w').write(''.join(set(open(vcffile).readlines())))

if __name__=="__main__":
#	args = parse_args()
	main()

	