import pandas as pd


df = pd.read_csv('/cluster/home/hugifl/TU_analysis/overlapping_TUs.tsv', sep='\t')

# Remove duplicate rows
df = df.drop_duplicates()

# Custom function to check if two strings are almost identical
def almost_identical(str1, str2):
    if str1 == str2:  # Check for exact match
        return True
    if len(str1) != len(str2):
        return False
    diff_count = sum(1 for a, b in zip(str1, str2) if a != b)
    return diff_count == 1

# Filter the DataFrame
df = df[~df.apply(lambda row: almost_identical(row['TU_Name'], row['Overlapping_TU']), axis=1)]
df.to_csv('/cluster/home/hugifl/TU_analysis/overlapping_TUs_reduced.tsv', sep='\t', index=False)

import pandas as pd

df.reset_index(drop=True, inplace=True)



# Create concatenated strings for comparison
df['Concatenated'] = df['TU_Name'] + df['Overlapping_TU']
# Identify rows to keep
rows_to_keep = set()
rows_to_not_check = set()
for i in range(len(df)):
    print(i)
    if i in rows_to_keep:
        continue
    if i in rows_to_not_check:
        continue
    for j in range(i + 1, len(df)):
        if almost_identical(df.iloc[i]['Concatenated'], df.iloc[j]['Concatenated']):
            if i not in rows_to_keep:
                rows_to_keep.add(i)
                rows_to_not_check.add(j)
            else:
                rows_to_not_check.add(j)
    rows_to_keep.add(i)

# Keep only the marked rows
df_reduced = df.loc[rows_to_keep]

print(df_reduced)

df_reduced.to_csv('/cluster/home/hugifl/TU_analysis/overlapping_TUs_reduced_more.tsv', sep='\t', index=False)