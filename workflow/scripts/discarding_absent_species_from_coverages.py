import pandas as pd
import sys
import os

dir_presence = os.path.dirname(sys.argv[2])
completeness = sys.argv[3]
kingdom = sys.argv[4]

coverage_df = pd.read_csv(sys.argv[1], sep='\t')
coverage_df.set_index("Contig", inplace=True)
coverage_df.rename(columns=lambda x: x.replace("_filtered_mapping_for_coverage_"+completeness+"_"+kingdom+" Trimmed Mean", ""), inplace=True)
coverage_df.rename(columns=lambda x: x.replace("_filtered_mapping_for_coverage_"+completeness+"_"+kingdom+" Read Count", ""), inplace=True)

for sample in coverage_df.columns:
    depth_df = pd.read_csv(dir_presence+"/metapresence_for_"+sample+"_from_"+completeness+"_"+kingdom+"_abundances.tsv", sep='\t')
    for mag, coverage in coverage_df[sample].items():
        if mag not in depth_df['genome'].values:
            coverage_df.at[mag, sample] = 0

coverage_df.to_csv(sys.argv[5], index=True, sep ='\t')
