import sys, os
from Analyze import createHomoProfile

def main():
	fasta_file, pdb_id, chain, output_file = getArgs()
	temp_file = '%s.tmp' % output_file
	createHomoProfile(fasta_file, temp_file)
	output = parseHomoProfileOutput(temp_file, pdb_id, chain)
	try:
		with open(output_file, 'w') as f:
			f.write(output)
	except:
		error('Error writing to %s' % output_file)
	finally:
		os.remove(temp_file)

def error(message):
	print(message)
	exit(1)

def getArgs():
	args = sys.argv
	num_args = len(args)
	if num_args != 5:
		error('Usage: python findHomologs.py FASTA_FILE PDB_ID CHAIN OUTPUT_FILE')
	fasta_file, pdb_id, chain, output_file = args[1:5]
	return fasta_file, pdb_id, chain, output_file

def parseHomoProfileOutput(temp_file, pdb_id, chain):
	try:
		with open(temp_file, 'r') as f:
			lines = map(lambda line: line.strip(), f.readlines())
			chain_lines = filter(lambda line: chainInLine(line)==chain, lines)
			last_chain_line = chain_lines[-1]
			fields = last_chain_line.split(' ')
			pdb = '%s_%s' % (pdb_id, chain)
			homologous_pdbs = ' '.join(fields[1:])
			pdbs = '%s %s' % (pdb, homologous_pdbs)
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
