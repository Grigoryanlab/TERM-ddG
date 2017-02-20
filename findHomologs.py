import sys, os
sys.path.insert(1, '/home/ifsdata/grigoryanlab/library/FanPythonMods/') # needed to find General if FanPythonMods is not already in the user's path
from Analyze import createHomoProfile

def main():
	fasta_file, chain, output_file = getArgs()
	temp_file = '%s.tmp' % output_file
	createHomoProfile(fasta_file, temp_file)
	output = parseHomoProfileOutput(temp_file, chain)
	try:
		with open(output_file, 'w') as f:
			f.write(output)
	except:
		error('Error writing to %s' % output_file)
	os.remove(temp_file)

def error(message):
	print(message)
	exit(1)

def getArgs():
	args = sys.argv
	num_args = len(args)
	if num_args != 4:
		error('Usage: python findHomologs.py FASTA_FILE CHAIN OUTPUT_FILE')
	fasta_file, chain, output_file = args[1:4]
	return fasta_file, chain, output_file

def parseHomoProfileOutput(temp_file, chain):
	try:
		with open(temp_file, 'r') as f:
			lines = map(lambda line: line.strip(), f.readlines())
			chain_lines = filter(lambda line: chainInLine(line)==chain, lines)
			last_chain_line = chain_lines[-1]
			fields = last_chain_line.split(' ')
			pdbs = ' '.join(fields[1:])
			return pdbs
	except:
		error('Error opening the temp file %s' % temp_file)

def chainInLine(line):
	fields = line.split(' ')
	if len(fields) >= 1:
		name_info = fields[0]
		infos = name_info.split('|')
		if len(infos) >= 1:
			pdb = infos[0]
			pdb__chain = pdb.split(':')
			if len(pdb__chain) >= 2:
				chain = pdb__chain[1]
				return chain
			else:
				return None
		else:
			return None
	else:
		return None

if __name__ == '__main__':
	main()
