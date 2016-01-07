## 1. Input and Output

Input will need the following source data ("xxx" means "real" or "simu" corresponding to real dataset or simulated dataset):

```
../data_xxx
../data_xxx/batch_individuals.txt
../data_xxx/batch_samples.txt
../data_xxx/genotype
../data_xxx/genotype/chrX
../data_xxx/genotype/chrX/SNP_dosage_IndividualID.txt
../data_xxx/genotype/chrX/SNP_info.txt
../data_xxx/list_individuals.txt
../data_xxx/expression.txt
../data_xxx/gene_tss.txt
../data_xxx/gene_xymt.txt
../data_xxx/list_samples_train.txt
```

it will also needs the following to initialize all the parameters (these files are prepared by other routines):

```
../result_init
../result_init/para_init_train_batch_batch_hidden.txt
../result_init/para_init_train_batch_hidden_gene.txt
../result_init/para_init_train_cis.txt
../result_init/para_init_train_snp_cellenv.txt
../result_init/para_init_train_cellenv_gene.txt
```

and the output will generate the following data (the model, or the set of coefficients):

```
../result
../result/etissue_list.txt
../result/para_batch_batch_hidden.txt
../result/para_batch_hidden_gene.txt
../result/para_snp_cellenv.txt
../result/para_cellenv_gene
../result/para_cellenv_gene/etissueX.txt
../result/para_cis_gene
../result/para_cis_gene/etissueX.txt
```

The following dir's are for debugging the program:

```
../temp_data (associated code: "opt_subroutine.cpp", to check some intermediate results, such as the derivatives, and the value before and after the logistic activation function)
../result_inter (associated code: "opt_para_save.cpp", called in "optimization.cpp" to save the results after each iteration across all tissues; not used in practice, as we will always quit the program when something is wrong in the learning)
```
