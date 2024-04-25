import pandas as pd

def create_gff3_from_tsv(tsv_file, gff3_file):
    df = pd.read_csv(tsv_file, sep='\t')

    with open(gff3_file, 'w') as gff3:
        gff3.write("##gff-version 3\n")
        gff3.write(f"#!genome-build U00096.3\n")

        for _, row in df.iterrows():
            seqid = "U00096.3" 
            source = "EcoCyc" 
            type_ = "TU" 
            start = row['TU_start']
            end = row['TU_end']
            score = "."  
            strand = row['Direction']
            phase = "."  
            attributes = f"ID={row['Name']};Name={row['Name']};Genes={row['Genes']}"

            gff3_line = f"{seqid}\t{source}\t{type_}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n"
            gff3.write(gff3_line)

def create_gff3_from_tsv_2(tsv_file, gff3_file):
    df = pd.read_csv(tsv_file, sep='\t')

    with open(gff3_file, 'w') as gff3:
        gff3.write("##gff-version 3\n")
        gff3.write(f"#!genome-build NC_000913.3\n")

        for _, row in df.iterrows():
            seqid = "NC_000913.3" 
            source = "EcoCyc" 
            type_ = "TU" 
            start = row['TU_start']
            end = row['TU_end']
            score = "."  
            strand = row['Direction']
            phase = "."  
            attributes = f"ID={row['Name']};Name={row['Name']};Genes={row['Genes']}"

            gff3_line = f"{seqid}\t{source}\t{type_}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n"
            gff3.write(gff3_line)



# Usage
tsv_filename = '/cluster/home/hugifl/TranscriptionalUnits/TU_start_end_added_promoters_updated_V2.tsv'  # Replace with your TSV file path
gff3_filename = '/cluster/home/hugifl/TranscriptionalUnits/transcriptional_units_V2_U00096.3.gff3'  # Output file name
gff3_filename_2 = '/cluster/home/hugifl/TranscriptionalUnits/transcriptional_units_V2_NC_000913.3.gff3'  # Output file name

create_gff3_from_tsv(tsv_filename, gff3_filename)
create_gff3_from_tsv_2(tsv_filename, gff3_filename_2)