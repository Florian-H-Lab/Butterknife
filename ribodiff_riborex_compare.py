
def find_diff_genes(results, t_pval, t_pvaladj, sep, indpval, indpvaladj):
    diffpval_genes = list()
    diffpvaladj_genes = list()
    for line in results:
        entries = line.split(sep)

        pval = entries[indpval]
        pvaladj = entries[indpvaladj]

        if(pval != "nan"):
            pval = float(pval)
            if(pval < t_pval):
                diffpval_genes.append(entries[0])

        if(pvaladj != "nan"):
            pvaladj = float(pvaladj)
            if(pvaladj < t_pvaladj):
                diffpvaladj_genes.append(entries[0])

    total = [diffpval_genes, diffpvaladj_genes]
    return (total)


def write_list_to_file(li, filename):
    file = open(filename, "w")
    for ele in li:
        file.write(ele + "\n")
    file.close()

control = "IgG"

data_path = "/home/flow/Documents/Bruno/"

ribodiff_results = open(data_path + "ribodiff/Med12_" + control + "_results.txt")
riborex_results = open(data_path + "riborex/Med12_" + control + "_results.txt")

header = ribodiff_results.readline()
header = riborex_results.readline()

t_pval = .05
t_pvaladj = .1

ribodiff_genes = find_diff_genes(ribodiff_results, t_pval, t_pvaladj, "\t", 2, 3)
ribodiff_results.close()

riborex_genes = find_diff_genes(riborex_results, t_pval, t_pvaladj, ",", 5, 6)
riborex_results.close()

ribodiff_diffpval_genes = ribodiff_genes[0]
ribodiff_diffpvaladj_genes = ribodiff_genes[1]

riborex_diffpval_genes = riborex_genes[0]
riborex_diffpvaladj_genes = riborex_genes[1]

write_list_to_file(ribodiff_diffpval_genes, data_path + "ribodiff/Med12_" + control + "_diffpval_genes.txt")
write_list_to_file(ribodiff_diffpvaladj_genes, data_path + "ribodiff/Med12_" + control + "_diffpvaladj_genes.txt")
write_list_to_file(riborex_diffpval_genes, data_path + "riborex/Med12_" + control + "_diffpval_genes.txt")
write_list_to_file(riborex_diffpvaladj_genes, data_path + "riborex/Med12_" + control + "_diffpvaladj_genes.txt")


compare_file = open(data_path + "compare_ribodiff_riborex/Med12_" + control + "_diffpval_genes.txt", "w")

for ele in ribodiff_diffpval_genes:
  if ele in riborex_diffpval_genes:
      compare_file.write(ele + "\n")
compare_file.close()



