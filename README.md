# TERM-ddG
Method for computing stability changes in proteins upon mutations, published in https://doi.org/10.1371/journal.pone.0178272

Authors: Fan Zheng, Gevorg Grigoryan

This package contains the source code for a method of predicting the effect of amino-acid point mutations on protein thermostability, reported in our paper:

**Sequence statistics of tertiary structural motifs reflect protein stability**, F. Zheng, G. Grigoryan, PLoS ONE, 12(5): e0178272, 2017

## Prerequisites

This package is written in Python and MATLAB and works under a UNIX/Linux-like environment (should be possible to make it work in MacOS, but we have not tried). We tested our code with Python 2.7 and MATLAB 2015b. It requires the following to work:

### python libraries:
  * SciPy (https://www.scipy.org/)
  * ProDy (http://prody.csb.pitt.edu/downloads/)

### third-party programs:
* USEARCH8.0 -- a program for rapidly finding similarities within sets of sequences (http://www.drive5.com/usearch/download.html)
* BLAST -- NCBI blast-2.2.26
* Matlab 2015b or higher (Octave should also work, but we have not tested this extensively)

### other tools from our lab:
* MASTER (Method of Accerlerated Search for Tertiary Ensemble Representatives); download from http://www.grigoryanlab.org/master/.

> For our paper, we have parsed the Protein Data Bank (PDB) and created a database that enables the MASTER program to search the tertiary motifs in the queried structure. The original database can be downloaded from the Grigoryan lab website via `rsync -varz arteni.cs.dartmouth.edu::masterDB-ddG/ /local/path`. One can also follow our protocol at https://vimeo.com/120274509 to customize the database, such as limiting search in a subset of the PDB, or simply update the database with the latest PDB data.

> Whichever database you end up using, use the script `create_seqdb.py` included with this package to create a sequence database, by running:
`python create_seqdb.py <list> <db_file>`, where `<list>` is a text file in which each line specifies the path to each PDB file in the database. and then set `PATH_seqdb` in `General.py` to `<db_file>`.

* confind -- a program that identifies mutually-influencing pairs of positions in proteins; download from http://www.grigoryanlab.org/confind/

> confind comes with an appropriately pre-formatted Dubrack 2010 banckbone-dependent rotamer library, which can be found under `rotlibs` in your downloaded confind folder

Once the required programs and data sources are installed, their paths should be updated in the "Interface" section of `General.py`. In particular, set the following variables:
  * `PATH_blast` -- path to NCBI blast `bin` directory (we tested with `blast-2.2.26`)
  * `PATH_confind` -- path to confind binary
  * `PATH_rotLib` -- path to rotamer library directory that comes with ConFind
  * `PATH_master` -- path to directory with MASTER binaries (`master` and `createPDS`)
  * `PATH_usearch` -- path to USEARCH binary
  * `PATH_seqdb` -- path to the sequence database that corresponds to the structural database used for prediction

Because the calculation involves a large number of MASTER searches and other procedures, the code is set up to submit "jobs." By default, submitting a job just means running it locally and waiting until it finishes (i.e., running in serial), but the code also provisions the use of a Sun Grid Engine (SGE) computing cluster (i.e., where jobs are submitted to an SGE cluster and the code waits for them to finish). It is easy to generalize this to an arbitrary cluster by changing `Cluster.py`. Specifically, only functions `submit` and `checkJobRun` need to change. We recommend that you designate a new value for the member variable `type` of class `jobOnCluster` that corresponds to your specific cluster, and define additional cases in functions `submit` and `checkJobRun` to deal with that cluster.


## Get started

An example of program output can be found in `data_demo`. The following instructions will show how these results can be generated from scratch.

1. Prepare a PDB file for the protein to be predicted. For example, `1EY0.pdb`.

2. Prepare a list of point mutations. For example, `1EY0.s350.tab`. The file is a tab-delimited file, and it should conform the following format:  
  `PDB_id Chain_id Wild_type_aa Position Mutant_aa experimental_ddG`  
  and the first column should agree with the name of the PDB file. For example, if the elements in the first column are `1EY0` then the PDB file should be named `1EY0.pdb` (no underscores).

3. Run the following command:  
 `python mutationListIteration.py --l 1EY0.s350.tab --homof S2648.homo`

Once the database has been prepared, users should make sure the `--db` argument in `mutationListIteration.py` point to the path of the database.

> We recommend to have a file to keep track of the PDB structures that have homologous relationships to the query protein. These PDB structures will be excluded from statistics of structure motif search. Otherwise, we believe the results will be biased. Please see the details in our paper. That is the purpose of the `--homof` flag. In `1EY0.homo` the PDB IDs following `1EY0_A` are excluded. Such results can be created from NCBI `blastpgb` program. The list here includes all PDB chains that are similar to `1EY0_A` at a cutoff e=1.0, which is consistent to our paper. This step requires a sequence database corresponding to the set of structures used in MASTER search. The sequence database can be downloaded at: (**glab**). After downloaded, the `PATH_seqdb` variable in the `General.py` needs to be changed. 

> If a customized database is used, the sequence database can be generated by the script `create_seqdb.py`. The first input argument of this script is a text file, where each line is a path to each .pdb file in the database. The `extractPDB` code in the MASTER package can convert .pds files to .pdb files. 


4. Upon finishing the previous command, run:  
 `python submitMatlab.py --l 1EY0.s350.tab --o 1ey0.dat`  
  and `1ey0.dat` will contain the results. A positive score indicates that the mutation is predicted as destabilizing.

5. We also provide the option to score all amino acid at a position. For that, simply run:  
 `python submitMatlab.py --l 1EY0.s350.tab --o 1ey0.dat.scan --scan`  
 and `1ey0.dat.scan` with show the score of 20 amino acids, with the wild-type amino acid specified in column 3 assigned as zeroes, and the other amino acids are relative to the wild-type amino acid.
> The step 5 is independent of step 4. However, if step 4 has been carried out, adding a flag `--oo` to the step 5 command will ultilize the existing results and give output much faster.

## Citations

If you use this tool in your research, please cite the following publications:

* Zheng, Fan and Gevorg Grigoryan. "Sequence statistics of tertiary structural motifs reflect protein stability", F. Zheng, G. Grigoryan, PLoS ONE, 12.5 (2017): e0178272.

* Zheng, Fan, Jian Zhang, and Gevorg Grigoryan. "Tertiary structural propensities reveal fundamental sequence/structure relationships." Structure 23.5 (2015): 961-971.

We also appreciate the authors of the following works that are utilized in this tool:

* Zhou, Jianfu, and Gevorg Grigoryan. "Rapid search for tertiary fragments reveals protein sequence–structure relationships." Protein Science 24.4 (2015): 508-524.

* Bakan, Ahmet, et al. "Evol and ProDy for bridging protein sequence evolution and structural dynamics." Bioinformatics 30.18 (2014): 2681-2683.

* Edgar, Robert C. "Search and clustering orders of magnitude faster than BLAST." Bioinformatics 26.19 (2010): 2460-2461.


