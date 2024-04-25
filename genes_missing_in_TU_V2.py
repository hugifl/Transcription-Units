import pandas as pd
import re

TU_df = pd.read_csv('/cluster/home/hugifl/TranscriptionalUnits/TU_start_end_added_promoters_V2.tsv', sep='\t')
gene_df = pd.read_csv('/cluster/home/hugifl/spacer_coverage_input/ECOCYC_genes.txt', sep='\t')

gene_df.dropna(inplace=True)


def find_unrepresented_genes(TU_df, gene_df):
    
    # First, create a set of unique genes from TU_df
    unique_genes_set = set()
    for genes in TU_df['Genes']:
        unique_genes_set.update(genes.split(','))

    # Check for genes in gene_df not represented in TU_df
    unrepresented_genes = []
    for index, row in gene_df.iterrows():
        gene = row['Gene_Name']
        direction = row['Direction']
        if gene not in unique_genes_set:
            unrepresented_genes.append((gene, direction))

    # Create a DataFrame for unrepresented genes including their direction
    unrepresented_genes_df = pd.DataFrame(unrepresented_genes, columns=['Unrepresented_Genes', 'Direction'])

    # Initialize new columns
    
    unrepresented_genes_df['number_of_overlapping_TUs'] = 0
    unrepresented_genes_df['overlapping_TUs'] = ''
    unrepresented_genes_df['overlapping_TUs_directions'] = ''
    unrepresented_genes_df['number_of_spanning_TUs'] = 0
    unrepresented_genes_df['spanning_TUs'] = ''
    unrepresented_genes_df['spanning_TUs_directions'] = ''
    unrepresented_genes_df['overlaps a gene in the TU'] = 0
    unrepresented_genes_df['TU contains only genes that are not in gene_df'] = 0

    # Check for overlap with TUs
    counter = 0
    for index, row in unrepresented_genes_df.iterrows():
        counter += 1
        print(f"Checking gene {counter} of {unrepresented_genes_df.shape[0]}")
        gene_name = row['Unrepresented_Genes']
        gene_left = gene_df[gene_df['Gene_Name'] == gene_name]['Left'].values[0]
        gene_right = gene_df[gene_df['Gene_Name'] == gene_name]['Right'].values[0]

        overlapping_TUs = []
        overlapping_TUs_directions = []
        completely_in_TU = []
        completely_in_TU_directions = []
        for _, tu_row in TU_df.iterrows():
            TU_direction = tu_row['Direction']
            if not (gene_right < tu_row['TU_start'] or gene_left > tu_row['TU_end']):
                overlapping_TUs.append(tu_row['Name'])
                overlapping_TUs_directions.append(TU_direction)
            if gene_left >= tu_row['TU_start'] and gene_right <= tu_row['TU_end']:
                completely_in_TU.append(tu_row['Name']) 
                completely_in_TU_directions.append(TU_direction)
                genes_in_tu = tu_row['Genes'].split(',')
                for gene_in_tu in genes_in_tu:
                    gene_in_tu_data = gene_df[gene_df['Gene_Name'] == gene_in_tu]
                    if not gene_in_tu_data.empty:
                        gene_in_tu_left = gene_in_tu_data['Left'].values[0]
                        gene_in_tu_right = gene_in_tu_data['Right'].values[0]

                        if not (gene_right < gene_in_tu_left or gene_left > gene_in_tu_right):
                            unrepresented_genes_df.at[index, 'overlaps a gene in the TU'] = 1
                            break
                    else:
                        unrepresented_genes_df.at[index, 'TU contains only genes that are not in gene_df'] = 1


        unrepresented_genes_df.at[index, 'number_of_overlapping_TUs'] = len(overlapping_TUs)
        unrepresented_genes_df.at[index, 'overlapping_TUs'] = ', '.join(overlapping_TUs)
        unrepresented_genes_df.at[index, 'number_of_spanning_TUs'] = len(completely_in_TU)
        unrepresented_genes_df.at[index, 'spanning_TUs'] = ', '.join(completely_in_TU)
        unrepresented_genes_df.at[index, 'overlapping_TUs_directions'] = ', '.join(overlapping_TUs_directions)
        unrepresented_genes_df.at[index, 'spanning_TUs_directions'] = ', '.join(completely_in_TU_directions)

        

    return unrepresented_genes_df



def add_new_TUs(TU_df, unrepresented_genes_df):
    # Iterate through the unrepresented_genes_df
    for index, row in unrepresented_genes_df.iterrows():
        gene_name = row['Unrepresented_Genes']
        gene_direction = row['Direction']
        overlapping_TUs_count = row['number_of_overlapping_TUs']
        overlapping_TUs_directions = row['overlapping_TUs_directions'].split(', ')
        TU_to_add = False

        # Check the criteria for adding the gene as a new TU
        if overlapping_TUs_count == 0:
            TU_to_add = True
        elif overlapping_TUs_count == 1 and gene_direction != overlapping_TUs_directions[0]:
            TU_to_add = True
        elif overlapping_TUs_count > 1 and not all(dir == gene_direction for dir in overlapping_TUs_directions):
            TU_to_add = True

        # Add the gene as a new TU if it meets the criteria
        if TU_to_add:
            new_TU = pd.Series({
                'Name': gene_name,
                'Direction': gene_direction,
                'TU_start': gene_df[gene_df['Gene_Name'] == gene_name]['Left'].values[0],
                'TU_end': gene_df[gene_df['Gene_Name'] == gene_name]['Right'].values[0],
                'Genes': gene_name
            })
            TU_df = TU_df.append(new_TU, ignore_index=True)

    return TU_df

def add_to_existing_TUs(TU_df, unrepresented_genes_df, gene_df):
    # Iterate through each unrepresented gene
    for index, row in unrepresented_genes_df.iterrows():
        gene_name = row['Unrepresented_Genes']
        gene_direction = row['Direction']
        overlaps_gene_in_TU = row['overlaps a gene in the TU']
        TU_only_contains_unrepresented_genes = row['TU contains only genes that are not in gene_df']
        spanning_TUs = row['spanning_TUs'].split(', ')

        # Check if the gene meets the criteria for adding to an existing TU
        for TU_name in spanning_TUs:
            TU_row = TU_df[TU_df['Name'] == TU_name]
            # Check if TU direction matches gene direction
            if not TU_row.empty and TU_row['Direction'].values[0] == gene_direction:
                # Add the gene to this TU
                current_genes = TU_row['Genes'].values[0]
                updated_genes = current_genes + ',' + gene_name
                TU_df.loc[TU_df['Name'] == TU_name, 'Genes'] = updated_genes

    return TU_df

unrepresented_genes_df = find_unrepresented_genes(TU_df, gene_df)
print(unrepresented_genes_df.shape[0])

unrepresented_genes_df.to_csv('/cluster/home/hugifl/TranscriptionalUnits/genes_missing_in_TUs.tsv', sep='\t', index=False)

TU_df_updated = add_new_TUs(TU_df, unrepresented_genes_df)

TU_df_updated = add_to_existing_TUs(TU_df_updated, unrepresented_genes_df, gene_df)

TU_df_updated.to_csv('/cluster/home/hugifl/TranscriptionalUnits/TU_start_end_added_promoters_updated_V2.tsv', sep='\t', index=False)