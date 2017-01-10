import sys

AAs = {0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F', 5: 'G', 6: 'H', 7: 'I', 8: 'K', 9: 'L', 10: 'M', 11: 'N', 12: 'P', 13: 'Q', 14: 'R', 15: 'S', 16: 'T', 17: 'V', 18: 'W', 19: 'Y'}

class Residue:
	def __init__(self, chain, position, amino_acid):
		self.chain = chain
		self.position = position
		self.amino_acid = amino_acid

def main():
	pdb_id, residues_str, mutation_list_file = get_args()
	residues = parseResiduesStr(residues_str)
	mutation_list = createMutationList(pdb_id, residues)
	storeMutationList(mutation_list_file, mutation_list)

def get_args():
	args = sys.argv
	num_args = len(args)
	if num_args != 4:
		error('Usage: python generateMutationList.py PDB_ID RESIDUES MUTATION_LIST_FILE\nThe residue list should be separated by commas and each residue should be in the form CHAIN-POSITION-AMINO_ACID')
	pdb_id, residues, mutation_list_file = args[1:4]
	return pdb_id, residues, mutation_list_file

def error(message):
	print(message)
	exit(1)

def parseResiduesStr(residues_str):
	residue_strs = residues_str.split(',')
	residues = map(parseResidueStr, residue_strs)
	return residues

def parseResidueStr(residue_str):
	fields = residue_str.split('-')
	num_fields = len(fields)
	if num_fields != 3:
		error('Each residue must contain three fields separated by hyphens: CHAIN-POSITION-AMINO_ACID')
	chain, position, amino_acid = fields[0:3]
	return Residue(chain, position, amino_acid)

def createMutationList(pdb_id, residues):
	mutation_list = reduce(lambda list, residue: addResidueMutations(pdb_id, list, residue), residues, [])
	return mutation_list

def addResidueMutations(pdb_id, mutation_list, residue):
	mutation_lines = map(lambda i: createMutationLine(pdb_id, mutateResidue(residue, i), residue.amino_acid), range(20))
	new_mutation_list = mutation_list + mutation_lines
	return new_mutation_list

def mutateResidue(residue, mutated_aa_index):
	return Residue(residue.chain, residue.position, AAs[mutated_aa_index])

def createMutationLine(pdb_id, residue, native_aa):
	return '%s\t%s\t%s\t%s\t%s\t0' % (pdb_id, residue.chain, native_aa, residue.position, residue.amino_acid)

def storeMutationList(mutation_list_file, mutation_list):
	mutated_list_contents = '\n'.join(mutation_list)
	mutated_list_fh = open(mutation_list_file, 'w')
	mutated_list_fh.write(mutated_list_contents)
	mutated_list_fh.close()

main()
