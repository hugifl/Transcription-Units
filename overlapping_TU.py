import pandas as pd

df = pd.read_csv('/cluster/home/hugifl/TU_analysis/TU_start_end_added_promoters_new.tsv', sep='\t')


# Function to check if two TUs overlap
def is_overlap(start1, end1, start2, end2):
    return max(start1, start2) <= min(end1, end2)

overlap_info = []
for i in range(len(df)):
    print(i)
    for j in range(i + 1, len(df)): 
        print(j)
        if is_overlap(df.iloc[i]["TU_start"], df.iloc[i]["TU_end"], df.iloc[j]["TU_start"], df.iloc[j]["TU_end"]):
            overlap_info.append({
                "TU_Name": df.iloc[i]["Name"],
                "TU_Direction": df.iloc[i]["Direction"],
                "Overlapping_TU": df.iloc[j]["Name"],
                "Overlapping_Direction": df.iloc[j]["Direction"]
            })

# Create a DataFrame from the overlap information
overlap_df = pd.DataFrame(overlap_info)
overlap_df.to_csv('/cluster/home/hugifl/TU_analysis/overlapping_TUs.tsv', sep='\t', index=False)

print("Transcriptional Units that overlap with others:", overlap_df.head())
