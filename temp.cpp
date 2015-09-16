



	//========================================================================
	// two step: forward propagation (get the function values); backward propagation (get the parameter derivatives)
	//========================================================================
	//========================================================================
	// step#1: ... (cis-; cell env; batch)
	//========================================================================
	//========================================================================
	// ****************************** [part1] cis- *********************************
	// for cis-, two issues:
	// 1. if this is a XYMT gene, jump;
	// 2. we use (gene_cis_index[gene].second - gene_cis_index[gene].first + 1) as the length of the cis- parameter array
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
		if ( got != gene_xymt_rep.end() )
		{
			(*expr_con_pointer)[i] = 0;
		}
		else
		{
			(*expr_con_pointer)[i] = 0;
			int chr = gene_tss[gene].chr;
			int num = gene_cis_index[gene].second - gene_cis_index[gene].first + 1;
			for(int k=0; k<num; k++)
			{
				int pos = gene_cis_index[gene].first + k;
				float dosage = (*dosage_list_pointer)[chr-1][pos];  // dosage at k position
				(*expr_con_pointer)[i] += dosage * para_cis_gene[etissue_index][i][k];
			}
		}
	}

	// ********************* [part2] cell env relevant parameters *********************
	// from snp to cell env variables
	for(int i=0; i<num_cellenv; i++)
	{
		(*cellenv_con_pointer)[i] = 0;
		long count = 0;
		for(int j=0; j<22; j++)  // across all the chromosomes
		{
			int chr = j+1;
			for(long k=0; k<snp_name_list[j].size(); k++)
			{
				float dosage = (*dosage_list_pointer)[j][k];
				(*cellenv_con_pointer)[i] += dosage * para_snp_cellenv[i][count];
				count ++;
			}
		}
	}
	//$$$$$$$$$$$ perform the activation function here (logistic or something else) $$$$$$$$$$$$
	for(int i=0; i<num_cellenv; i++)
	{
		(*cellenv_con_pointer)[i] = 1 / ( 1 + exp( - (*cellenv_con_pointer)[i] ));
	}
	// from cell env variables to genes
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		for(int j=0; j<num_cellenv; j++)
		{
			(*expr_con_pointer)[i] += para_cellenv_gene[etissue_index][i][j] * (*cellenv_con_pointer)[j];
		}
	}

	// ********************* [part3] linear or non-linear batches *********************
	// from original batch to hidden batch
	for(int i=0; i<num_batch_hidden; i++)
	{
		(*batch_hidden_con_pointer)[i] = 0;
		for(int j=0; j<num_batch; j++)
		{
			// DEBUG
			// print out the item
			//cout << (*batch_hidden_con_pointer)[i] << endl;
			//cout << (*batch_list_pointer)[j] << endl;
			//cout << para_batch_batch_hidden[i][j] << endl;

			(*batch_hidden_con_pointer)[i] += (*batch_list_pointer)[j] * para_batch_batch_hidden[i][j];
		}
	}
	//$$$$$$$$$$$ perform the activation function here (logistic or something else) $$$$$$$$$$$$
	for(int i=0; i<num_batch_hidden; i++)
	{
		(*batch_hidden_con_pointer)[i] = 1 / ( 1 + exp( - (*batch_hidden_con_pointer)[i] ));
	}
	// from hidden batch to genes
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		for(int j=0; j<num_batch_hidden; j++)
		{
			(*expr_con_pointer)[i] += para_batch_hidden_gene[i][j] * (*batch_hidden_con_pointer)[j];
		}
	}



	//========================================================================
	//========================================================================
	// step#2: ... (cis-;  cell env; batch)
	//========================================================================
	//========================================================================
	// *********************** [part1] cis- ************************
	// pseudo: (expected rpkm - real rpkm) * genotype
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
		if ( got != gene_xymt_rep.end() )
		{
			continue;
		}
		else
		{
			int chr = gene_tss[gene].chr;
			int num = gene_cis_index[gene].second - gene_cis_index[gene].first + 1;
			for(int k=0; k<num; k++)
			{
				int pos = gene_cis_index[gene].first + k;
				float dosage = (*dosage_list_pointer)[chr-1][pos];  // dosage at k position
				(*para_dev_cis_gene_pointer)[etissue_index][i][k] += ((*expr_con_pointer)[i] - (*expr_list_pointer)[i]) * dosage;
			}
		}
	}

	// ***************** [part2] cell env relevant parameters *****************
	// from cell env to genes
	// pseudo: (expected rpkm - real rpkm) * cell_env
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		for(int j=0; j<num_cellenv; j++)
		{
			(*para_dev_cellenv_gene_pointer)[etissue_index][i][j] += ((*expr_con_pointer)[i] - (*expr_list_pointer)[i]) * (*cellenv_con_pointer)[j];
		}
	}
	// from snp to cell env
	// pseudo: [ \sum w3 * (expected rpkm - real rpkm) ] * g'(w2 * x1) * x1
	for(int i=0; i<num_cellenv; i++)
	{
		float temp = 0;
		for(int t=0; t<num_gene; t++)
		{
			temp += para_cellenv_gene[etissue_index][t][i] * ((*expr_con_pointer)[t] - (*expr_list_pointer)[t]);
		}
		temp *= (*cellenv_con_pointer)[i] * ( 1 - (*cellenv_con_pointer)[i] );
		long count = 0;
		for(int j=0; j<22; j++)  // across all the chromosomes
		{
			int chr = j+1;
			for(long k=0; k<snp_name_list[j].size(); k++)
			{
				float dosage = (*dosage_list_pointer)[j][k];
				(*para_dev_snp_cellenv_pointer)[i][count] += temp * dosage;
				count ++;
			}
		}
	}

	// ********************* [part3] linear or non-linear batches *********************
	// from hidden batch to genes
	// pseudo: (expected rpkm - real rpkm) * hidden batch var
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		for(int j=0; j<num_batch_hidden; j++)
		{
			(*para_dev_batch_hidden_gene_pointer)[i][j] += ((*expr_con_pointer)[i] - (*expr_list_pointer)[i]) * (*batch_hidden_con_pointer)[j];
		}
	}
	// from original batch to hidden batch
	// pseudo: [ \sum w5 * (expected rpkm - real rpkm) ] * g'(w4 * x2) * x2
	for(int i=0; i<num_batch_hidden; i++)
	{
		float temp = 0;
		for(int t=0; t<num_gene; t++)
		{
			temp += para_batch_hidden_gene[t][i] * ((*expr_con_pointer)[t] - (*expr_list_pointer)[t]);
		}
		temp *= (*batch_hidden_con_pointer)[i] * ( 1 - (*batch_hidden_con_pointer)[i] );
		for(int j=0; j<num_batch; j++)
		{
			float batch_value = (*batch_list_pointer)[j];
			(*para_dev_batch_batch_hidden_pointer)[i][j] += temp * batch_value;
		}
	}


