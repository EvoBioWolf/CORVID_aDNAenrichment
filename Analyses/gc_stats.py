import os
import pandas as pd

#AT/GC
input_directory = "/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/00_baitscomparison/gc"  # Update with your actual directory path
output_file = "GC_AT_dropout_summary.txt"

# Initialize a list to store the results
data = []

# Loop through each file in the directory
for filename in os.listdir(input_directory):
    if filename.endswith("_summary_metrics.txt"):
        sample_name = "_".join(filename.split("_")[:-2])
        with open(os.path.join(input_directory, filename), 'r') as file:
            for line in file:
                # Look for the summary line that contains GC and AT Dropout metrics
                if line.startswith("## METRICS CLASS"):
                    next(file)  
                    summary_line = next(file).strip().split("\t")
                    gc_dropout = summary_line[6]
                    at_dropout = summary_line[5]
                    data.append([sample_name, gc_dropout, at_dropout])
                    break

# Open the output file and write the data
data.sort(key=lambda x: x[0])
with open(output_file, 'w') as outfile:
    outfile.write("Sample\tGC_Dropout\tAT_Dropout\n")
    for row in data:
        outfile.write("\t".join(row) + "\n")

print(f"GC and AT Dropout summary saved to {output_file}")

#for mybaits GC coverage
output_file = "GC_coverage_summary_mybaits.txt"

# Initialize a dictionary to store dataframes for each sample
coverage_data = {}

# Loop through each file in the directory
for filename in os.listdir(input_directory):
    if filename.endswith("_mybaits_gc_bias_metrics.txt"):
        sample_name = "_".join(filename.split("_")[:-3])
        file_path = os.path.join(input_directory, filename)
        df = pd.read_csv(file_path, sep="\t", comment='#')  # Assuming tab-separated values and skipping comment lines
        
        # Select relevant columns: GC content and NORMALIZED_COVERAGE, rename for merging
        if 'GC' in df.columns and 'NORMALIZED_COVERAGE' in df.columns:
            df = df[['GC', 'NORMALIZED_COVERAGE']].rename(columns={'NORMALIZED_COVERAGE': sample_name})
            # Add the dataframe to the dictionary
            coverage_data[sample_name] = df
        else:
            print(f"Warning: File {filename} does not contain expected columns. Skipping.")

# Merge all dataframes on 'GC_CONTENT' column
merged_df = pd.DataFrame()
for sample, df in coverage_data.items():
    if merged_df.empty:
        merged_df = df
    else:
        merged_df = merged_df.merge(df, on='GC', how='outer')

sorted_columns = ['GC'] + sorted([col for col in merged_df.columns if col != 'GC'])
merged_df = merged_df[sorted_columns]

# Save the merged DataFrame to a text file with tab-separated values
merged_df.to_csv(output_file, sep='\t', index=False)

print(f"Normalized coverage values merged and saved to {output_file}")

#for mybaits GC coverage
output_file = "GC_coverage_summary_twist.txt"

# Initialize a dictionary to store dataframes for each sample
coverage_data = {}

# Loop through each file in the directory
for filename in os.listdir(input_directory):
    if filename.endswith("_TE_gc_bias_metrics.txt"):
        sample_name = "_".join(filename.split("_")[:-3])
        file_path = os.path.join(input_directory, filename)
        df = pd.read_csv(file_path, sep="\t", comment='#')  # Assuming tab-separated values and skipping comment lines
        
        # Select relevant columns: GC content and NORMALIZED_COVERAGE, rename for merging
        if 'GC' in df.columns and 'NORMALIZED_COVERAGE' in df.columns:
            df = df[['GC', 'NORMALIZED_COVERAGE']].rename(columns={'NORMALIZED_COVERAGE': sample_name})
            # Add the dataframe to the dictionary
            coverage_data[sample_name] = df
        else:
            print(f"Warning: File {filename} does not contain expected columns. Skipping.")

# Merge all dataframes on 'GC_CONTENT' column
merged_df = pd.DataFrame()
for sample, df in coverage_data.items():
    if merged_df.empty:
        merged_df = df
    else:
        merged_df = merged_df.merge(df, on='GC', how='outer')

sorted_columns = ['GC'] + sorted([col for col in merged_df.columns if col != 'GC'])
merged_df = merged_df[sorted_columns]

# Save the merged DataFrame to a text file with tab-separated values
merged_df.to_csv(output_file, sep='\t', index=False)

print(f"Normalized coverage values merged and saved to {output_file}")
