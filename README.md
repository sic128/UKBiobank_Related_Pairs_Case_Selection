# UKBiobank_Related_Pairs_Case_Selection v1
Given a phenotype list of individuals in UKBB, this script generates an unrelated subset that maximizes case individuals, based on a designated pihat (or kinship) value.

## Getting started

#### Python dependencies
- argsparse
- pandas
- numpy
- time
- sys
- os

#### Data
Some Data files will be required to run UKBiobank_Related_Pairs_Case_Selection. Files need to be downloaded from UK Biobank and modified to ensure the correct format.
The format is described below

1. Kinship File
   - Three columns: "IID1 IID2 Kinship"
   - Space Delimited, no header
2. Phenotype File
   - Three columns: "FID IID Phenotype"
   - Space Delimited, no header
3. Sample File
   - Two Columns: "FID IID"
   - Space Delimited, no header

## Typical analysis
For a typical analysis on UK Biobank (UKBB), use kinship file provided by UKBB as input. Generate sample file used for your analysis, usually a subset of UKBB subjects. Using provided phenotype files from UKBB, generate a phenotype file. Determine the pihat or kinship value to use as a threshold for relatedness (See flags for more info). Run Case selection script.


## Flags for analysis
### Inputs

##### `--pheno` 
Required
Phenotype file (Three column "FID IID Phenotype" ): Should only contain 3 class of values in the second column: Case Indicator, Control Indicator, and NA. Case-Control Indicator Coding can be user specified, but missing must be set to "NA".


        # Example Phenotype File:
        ##############
        #FID1 IID1 0
        #FID2 IID2 0
        #FID3 IID3 NA
        #FID4 IID4 1
        #FID5 IID5 0
        ##############


##### `--case_value`
####Required
Which Value in phenotype file second column represents cases. Specify indicator code (Number/Alpha/Text) to select cases (Alternatively can be used to select controls)


##### `--pihat`
####Optional
Pihat threshold (Optional Flag, can use --kinship_threshold instead)
Note: Anyone pair of relationship below the chosen pihat threshold will be considered unrelated. Pairs that are above the chosen pihat threshold are considered related.

        # For reference:
        # Pihat   Relationship
        # 1       Identical twins
        # 0.5     First Degree Relative (Siblings, Parents, Children)
        # 0.25    Second Degree Relatives (grandparents, grandchildren, aunts, uncles, nephews, nieces or half-siblings)
        # 0.125   Third Degree Relatives (first-cousins, great-grandparents or great grandchildren)
        # 0.0625  First-Cousins once removed
        # 0.03125 Fourth Degree Relatives (Second-Cousins)


##### `--kinship_threshold`
Optional
Kinship threshold (Optional Flag, can use --pihat instead)


##### `--kinship`
Required
UKBB kinship file. Three columns: "FID IID Kinship". Space delimited and no header.

        # Example Kinship File:
        ##############
        #IID1 IID2 0.5
        #IID2 IID3 0.25
        #IID1 IID3 0.15
        #IID5 IID6 0.03
        #IID7 IID8 0.4
        ##############


##### `--samples`
Required
UKBB full list of samples we are considering. Two columns: "FID IID", Space delimited, no headers.

         #Example UKBB sample File:
         ####################
         #FID1 IID1
         #FID2 IID2
         #FID3 IID3
         ####################


##### `--output`
Required
Output file name

### Sample command
```
python3 Related_Pair_Selection_UKBB.py \
--pheno pheno.txt \
--case_value 1 \
--pihat 0.0884  \
--kinship kinship.txt \
--samples sample.txt \
--output output
```

## Running the software on dbGaP cohorts
TODO

## Contact
Steven Cao(email: yic127@ucsd.edu) developed and maintains this software.
