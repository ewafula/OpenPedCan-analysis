# import required module
import pandas as pd
 
# read segment arrays
# paired SNParray
paired_seg_df = pd.read_csv("nexus_tumor_paired_824_geo.seg", sep="\t", dtype={"Start.Position": "int64", "End.Position": "int64"})
#paired_seg_df = paired_seg_df[(paired_seg_df["Chromosome"] == "2") & (paired_seg_df["Start.Position"] >= 12200001) & (paired_seg_df["End.Position"] <= 16700000)]
paired_mycn_df = paired_seg_df[(paired_seg_df["Chromosome"] == "2") & (16080686 >= paired_seg_df["Start.Position"]) & (16080686 <= paired_seg_df["End.Position"]) & (16087129 >= paired_seg_df["Start.Position"]) & (16087129 <= paired_seg_df["End.Position"])]
paired_mycn_df["Sample"] = paired_mycn_df["Sample"].apply(lambda x: x.split("-")[0])
# print(paired_mycn_df)

# single SNParray
single_seg_df = pd.read_csv("nexus_tumor_single_924_geo.seg", sep="\t", dtype={"Start.Position": "int64", "End.Position": "int64"})
#single_seg_df = single_seg_df[(single_seg_df["Chromosome"] == "2") & (single_seg_df["Start.Position"] >= 12200001) & (single_seg_df["End.Position"] <= 16700000)]
single_mycn_df = single_seg_df[(single_seg_df["Chromosome"] == "2") & (16080686 >= single_seg_df["Start.Position"]) & (16080686 <= single_seg_df["End.Position"]) & (16087129 >= single_seg_df["Start.Position"]) & (16087129 <= single_seg_df["End.Position"])]
# print(single_mycn_df)


metadata_df = pd.read_csv("metadata.txt", sep="\t", dtype=str)
metadata_df.columns = metadata_df.columns.str.replace(': ','_')
for row in metadata_df.itertuples():
    paired_mycn_df.loc[paired_mycn_df["Sample"] == row.Sample, "Sample"] = row.characteristics_usi
    single_mycn_df.loc[single_mycn_df["Sample"] == row.Sample, "Sample"] = row.characteristics_usi

print(paired_mycn_df)
print(single_mycn_df)
paired_mycn_df.to_csv("mycn_tumor_paired_segment_SNParray.tsv", sep="\t", index=False, encoding="utf-8")
single_mycn_df.to_csv("mycn_tumor_single_segment_SNParray.tsv", sep="\t", index=False, encoding="utf-8")