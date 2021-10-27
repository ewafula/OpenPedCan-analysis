# import required module
from pathlib import Path
import pandas as pd
 
# assign directory
snp_array_dir = 'hg19-files'
 
# iterate over snp files in the direcotry 
# create a merged MYCN df all tumor files
mycn_snps_dfs = []
snp_array_files = Path(snp_array_dir).glob('*tumor_*_hg19.txt')
for snp_array_file in snp_array_files:
    print(snp_array_file)
    df = pd.read_csv(snp_array_file, sep="\t", skiprows=10, dtype={"Position": "int64"})
    mycn_df = df[(df.Chr == "2") & (df.Position >= 16080686) & (df.Position <= 16087129)]
    print(mycn_df)
    mycn_snps_dfs.append(mycn_df)
    mycn_snps_df = pd.concat(mycn_snps_dfs, sort=False, ignore_index=True)
print(mycn_snps_df)
mycn_snps_df["Sample ID"] = mycn_snps_df["Sample ID"].apply(lambda x: x.split("_")[0])
mycn_snps_df.to_csv("mycn_tumor_focal_SNParray.tsv", sep="\t", index=False, encoding="utf-8")