import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read the TSV file into a DataFrame
df = pd.read_csv('/cluster/home/hugifl/TranscriptionalUnits/TU_start_end_added_promoters_updated.tsv', sep='\t')
outdir = '/cluster/home/hugifl/TranscriptionalUnits/'
# Print the first five rows of the DataFrame
print(df.head())
# add length column
df['TU_length'] = abs(df['TU_start'] - df['TU_end'])

# print max length
print("max length: " + str(df['TU_length'].max()))
print("min length: " + str(df['TU_length'].min()))
print("mean length: " + str(df['TU_length'].mean()))
print("median length: " + str(df['TU_length'].median()))
print("std length: " + str(df['TU_length'].std()))

# sort by length
df_sorted = df.sort_values(by=['TU_length'], ascending=[True])
print("minum 10 length: ", df_sorted.head(10))



# Plotting with matplotlib
plt.figure(figsize=(10, 6))
plt.hist(df['TU_length'], bins=50, alpha=0.7)
plt.title('Distribution of Transcriptional Units Length')
plt.xlabel('TU Length')
plt.ylabel('Counts')
plt.savefig(outdir + 'TU_length_distribution.png')
# Plotting with seaborn for a smoother density plot
#plt.figure(figsize=(10, 6))
#sns.histplot(df['TU_length'], kde=True, bins=50)
#plt.title('Distribution of Transcriptional Units Length')
#plt.xlabel('TU Length')
#plt.ylabel('Density')
#plt.show()