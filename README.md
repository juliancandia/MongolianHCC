# MongolianHCC
Analysis pipeline of RNA-Seq and WES data from a Mongolian hepatocellular carcinoma cohort.

NOTE: The Whole Exome Sequencing-based analysis relies on the DATA/ORIGINAL/mutect2_merged.maf file, which (at 218 MB) was too large to upload to a GitHub folder. This file was split into two files, which were subsequently compressed. To rebuild the original file needed for the analysis, proceed to the DATA/ORIGINAL folder and perform the following command-line operations:

gunzip mutect2_merged_TOP.maf.gz mutect2_merged_BOTTOM.maf.gz

cat mutect2_merged_TOP.maf mutect2_merged_BOTTOM.maf > mutect2_merged.maf
