# Butterknife_RIP

## 1. createcounttable.R

### Explanation
Create the input data for `ribodiff`.

### Call
Rscript createcounttable.R -i data_path -s singal_name -b background_name -c conbditions (comma separated, e.g.: Wildtype,Mutant) -r number_of_replicates -o output_path

### Input
Count data from `htseq-count`.

## 2. Ribodiff

### Explanation
Ribodiff allows to call differentially translated regions across two different conditions, incorporating Control data.

### Install
- Ribodiff can be found [here](https://github.com/ratschlab/RiboDiff). All credits go for the developers; published in [Zhong et al.](https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtw585)
- Install ribodiff via bioconda in an virtual environment (ribodiff)
- activate env
- Clone Ribodiff from git
- Go into the scripts folder where TE.py is

### Call
python2 TE.py -e conditionsheet -c counttable -o outputfile

### Input
Data from `createcounttable.R`

## 3. Riborex

### Explanation
Riborex allows to call differentially translated regions across two different conditions, incorporating Control data.

### Install
- Riborex was published in [Li et al.](https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtx047)
- Install all packages stated in the header of `riborex.R`.

### Call
Rscript riborex -i data_path -s singal_name -b background_name -c conbditions (comma separated, e.g.: Wildtype,Mutant) -r number_of_replicates -o output_path

### Input
Count data from `htseq-count`.

## 4. ribodiff_riborex_compare.py

### Explanation
Script filters results of ribodiff and riborex based on p-value threshold (<0.05) and compare the two lists of ribodiff and riborex to find regions metnioned in both lists.

### Call
python ribodiff_riborex_compare.py

### Input
Results of `riborex.R` and `ribodiff`.

## 5. annotate_gene_lists.sh

### Explanation
Very small scripts to annotate the results of `ribodiff_riborex_compare.py`.

### Call
annotate_gene_lists.sh regions_file annoation_file output_file

### Input
Result of `ribodiff_riborex_compare.py`
