import pandas as pd

# Read the TSV file into a DataFrame
df = pd.read_csv('/cluster/home/hugifl/TranscriptionalUnits/TUs_with_coordinates_genes.txt', sep='\t')
print("unique" + str(df.shape[0]))
# Assuming df is your DataFrame
df = df.drop_duplicates()
print("unique" + str(df.shape[0]))


promoter_df = pd.read_csv('/cluster/home/hugifl/ECOCYC_promoters_with_direction.txt', sep='\t')
promoter_df.dropna(inplace=True)

def rename_duplicates_with_numbering(df):
    # Get the first column name
    first_column = df.columns[0]

    # Create a dictionary to keep track of counts
    counts = dict()

    # Iterate through each element in the first column and count occurrences
    for value in df[first_column]:
        counts[value] = counts.get(value, 0) + 1

    # Process only those values that have duplicates
    for value, count in counts.items():
        if count > 1:
            current_count = 1
            for i in df.index[df[first_column] == value]:
                df.at[i, first_column] = f"{value}_{current_count}"
                current_count += 1

    return df

def extract_tu_coordinates(coord_str):
    # Check if coord_str is a string and not NaN
    if isinstance(coord_str, str) and 'NIL' not in coord_str:
        # Remove square brackets and replace arrow symbols with a standard delimiter
        coords = coord_str.replace('[', '').replace(']', '').replace(' &rarr; ', '-').replace(' &larr; ', '-')

        # Extract the coordinates
        _, coords = coords.split(' ', 1)
        start, end = coords.split('-')

        return int(start.replace(',', '')), int(end.replace(',', ''))
    return None, None


# Function to extract promoter coordinates
def extract_promoter_coordinates(coord_str):
    # Check if coord_str is a string and not NaN
    if isinstance(coord_str, str) and 'NIL' not in coord_str:
        # Remove square brackets and commas, then extract promoter coordinate
        promoter_coord = coord_str.strip('[]').replace(',', '')
        return int(promoter_coord.split(' ')[1])
    return None

# Function to search for matching promoters
def find_matching_promoters(tu_name, promoter_df):
    # Ensure tu_name is a string
    tu_name_str = str(tu_name)
    tu_name_first_gene_str = tu_name_str.split('-')[0]
    matching_promoters = promoter_df[promoter_df['Promoter_Name'].str.contains(tu_name_first_gene_str, regex=False, na=False)]
    return matching_promoters

def extract_TU_genes(gene_str):
    gene_str_out = gene_str.replace(" // ", ",")
    return gene_str_out 



# Initialize the new dataframe
TU_start_end = pd.DataFrame(columns=['Name', 'Direction', 'TU_start', 'TU_end'])

# Process each row in the original dataframe
for index, row in df.iterrows():
    TU_start, TU_end = extract_tu_coordinates(row['Sequence - coordinates of DNA region'])
    direction = row['Direction']
    genes = extract_TU_genes(row['Genes of transcription unit'])
    promoter_info = row['Transcription units - promoter of operon']
    promoter_coord = None
    if not pd.isna(promoter_info):
        promoter_coord = extract_promoter_coordinates(row['Sequence - coordinates of DNA region of promoter'])

    if promoter_coord is None or 'NIL' in promoter_info:
        # Find matching promoters and add rows for each
        matching_promoters = find_matching_promoters(row['Transcription-Units'], promoter_df)

        if not matching_promoters.empty:
            tu_count = 1
            valid_promoters = 0
            for _, promoter_row in matching_promoters.iterrows():
                use_promoter = False
                promoter_pos = promoter_row['Absolute_Plus_1_Position']
                if direction == '+':
                    if promoter_row['Direction'] == direction and abs(promoter_pos - TU_start) < 500:
                        TU_start = promoter_pos
                        valid_promoters += 1
                        use_promoter = True
                else:
                    if promoter_row['Direction'] == direction and abs(promoter_pos - TU_end) < 500:
                        TU_end = promoter_pos
                        valid_promoters += 1
                        use_promoter = True

                if use_promoter:
                    TU_start_end = TU_start_end.append({'Name': f"{row['Transcription-Units']}_{tu_count}",
                                                        'Direction': direction,
                                                        'TU_start': TU_start,
                                                        'TU_end': TU_end,
                                                        'Genes': genes}, ignore_index=True)
                    tu_count += 1

        if matching_promoters.empty or valid_promoters == 0:
            # No matching promoters found or promoter info is NIL, use original TU coordinates
            TU_start_end = TU_start_end.append({'Name': row['Transcription-Units'],
                                                'Direction': direction,
                                                'TU_start': TU_start,
                                                'TU_end': TU_end,
                                                'Genes': genes}, ignore_index=True)            
    else:
        # Use the promoter from the original row
        if direction == '+':
            TU_start = promoter_coord
        else:
            TU_end = promoter_coord
        TU_start_end = TU_start_end.append({'Name': row['Transcription-Units'],
                                            'Direction': direction,
                                            'TU_start': TU_start,
                                            'TU_end': TU_end,
                                            'Genes': genes}, ignore_index=True)


TU_start_end = rename_duplicates_with_numbering(TU_start_end)


print(TU_start_end.head())
TU_start_end.to_csv('/cluster/home/hugifl/TranscriptionalUnits/TU_start_end_added_promoters.tsv', sep='\t', index=False)
