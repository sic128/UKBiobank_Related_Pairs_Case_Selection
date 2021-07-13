#!/usr/bin/env python3
################################################################################################
#
#
#       Name : Related_Pair_Selection_UKBB.py 
#       
#       Description : Given a phenotype list of individuals in UKBB, this 
#                     script generates an unrelated subset that maximizes case individuals, based on
#                     a designated pihat value. 
#       
#       
#       Usage : python3 /data/nrnb03/REF/SOFTWARE/PYTHONCode/Related_Pair_Selection_UKBB.py --pheno [pheno file] --case_value [Case Value]
#                       --pihat [pihat threshold] --kinship [kinship table] --samples [samples list] --output [output file]
#       
#
#       Example : python3 /data/nrnb03/REF/SOFTWARE/PYTHONCode/Related_Pair_Selection_UKBB.py --pheno pheno.txt --case_value 1 
#                --pihat 0.125 --kinship kinship.txt --samples all.txt --output test
#       
#       By: rsalem and scao
#       
#       Date: 7/2/2021 
#       
#       
##################################################################################################
# READ NOTES:
# Note 1: Input files are assumed to be space delimited with no headers. 
# Note 2: FID and IID should be the same.
# Note 3: Assumes that the kinship file contains all relationships amongst your samples.
# Note 4: If samples are missing phenotype, imputing NA for the sample
# Note 5: Pihat and kinship are not the same. UKBB by default gives a file with kinship values which we use as the input relatedness file for simplicity,
#         but the --pihat flag requires a pihat threshold. 2*kinship=pihat. See below for pihat values and its corresponding relationship.
# Note 6: Outputs a file with one column containing IIDs of the unrelated subset that maximizes case individuals.

# Argument --pheno: Phenotype file (Two column "IID Phenotype" ): Should only contain 3 class of values in the second column: 
#                   Case Indicator, Control Indicator, and NA. Case-Control Indicator Coding can be user specified, but missing must be set to "NA".  
    
    
        # Example Phenotype File:
        ##############
        #FID1 IID1 0
        #FID2 IID2 0
        #FID3 IID3 NA
        #FID4 IID4 1
        #FID5 IID5 0
        ##############

# Argument --case_value: Which Value in phenotype file second column represents cases
#                        Specify indicator code (Number/Alpha/Text) to select cases (Alternatively can be used to select controls)

# Argument --pihat: Pihat threshold 
#                   Note: Anyone pair of relationship below the chosen pihat threshold will be considered unrelated
#                   Pairs that are above the chosen pihat threshold are considered related.
    
        # For reference:
        # Pihat   Relationship
        # 1       Identical twins
        # 0.5     First Degree Relative (Siblings, Parents, Children)
        # 0.25    Second Degree Relatives (grandparents, grandchildren, aunts, uncles, nephews, nieces or half-siblings)
        # 0.125   Third Degree Relatives (first-cousins, great-grandparents or great grandchildren) 
        # 0.0625  First-Cousins once removed
        # 0.03125 Fourth Degree Relatives (Second-Cousins)

# Arg --output: Output name

# Arg --kinship: UKBB kinship file: Three columns: "FID IID Kinship"

# Arg --samples: UKBB full list of samples we are considering. Two columns: "FID IID"

         #Example UKBB sample File:
         ####################
         #FID1 IID1
         #FID2 IID2
         #FID3 IID3
         #################### 


# Import libraries
import argparse
import pandas as pd
import numpy as np
import time
import sys
import os


# Name: check_Pheno 
# Description: Check if phenotype file second column contains only 3 types of values (Case value, Control Value and NA) 
#
# Parameters (1) : 1. pheno (str): phenotype file path 
#
# Output : None 

def check_Pheno(pheno):

    data = pd.read_csv(pheno,header=None,names=["FID","IID","Status"],delimiter=" ")  
    values = data.Status.unique() 

    if (len(values)>3) or (len(values)<1):
        print("Phenotype File Second Column incorrect. Make sure that there are no headers and contains only 3 possible elements: Case Code, Control Code, NA.")   
        print("Exiting")
        sys.exit()

# Name: get_Related 
# Description: Outputs a dataframe with two columns, these samples are related to each other according to the 
#              pihat threshold 
#
# Parameters (3) : 1. all_samples (pd data frame) : dataframe with two columns: FID,IID. All samples
#                     we are considering. 
#                  2. kinship_file (str) : path to kinship file input
#                  3. kinship_threshold (float) : kinship threshold to consider relatedness 
#
# Output : Array : First element of array is score. Second element is how many relatives it has.

# Outputs a dataframe with two columns, these samples are related to each other according to the pihat threshold
def get_Related(all_samples,kinship_file,kinship_threshold):

    # Load kinship file
    kinship = pd.read_csv(kinship_file,header=None,names=["IID1","IID2","Kinship"],delimiter=" ")
    # Filter kinship file so that pairs that are remaining are all in sample file
    kinship2 = kinship[kinship.IID1.isin(all_samples.IID)]
    kinship3 = kinship2[kinship2.IID2.isin(all_samples.IID)]
    # Filter out by kinship threshold
    kinship4 = kinship3[kinship3["Kinship"] >= kinship_threshold ]
    # Filter the Kinship column out
    kinship5 = kinship4[["IID1", "IID2"]]

    return kinship5

# Name: get_Strict_Unrelated 
# Description: Outputs a dataframe with one column, whom are individuals who are not related to anyone 
#              because they don't meet the pihat threshold with other samples 
#
# Parameters (2) : 1. related_people (pd data frame) : dataframe created with two columns: IID1, IID2
#                     of related individuals 
#                  2. all_samples (pd data frame) : dataframe with two columns: FID,IID. All samples 
#                     we are considering. 
#
# Output : pd data frame with two columns: FID,IID 

def get_Strict_Unrelated(related_people,all_samples):

    # Load sample file
    unrelated = all_samples[~all_samples.IID.isin(related_people.IID1)] 
    unrelated2 = unrelated[~unrelated.IID.isin(related_people.IID2)] 
    
    return unrelated2

# Name: get_Pheno 
# Description: Filter out phenotype file and make sure that all samples we are running the script on has 
#              phenotypes. If sample doesn't have phenotype, impute NA.
#
# Parameters (2) : 1. pheno_file (str): input phenotype file path 
#                  2. related people (pd data frame) : data frame created with two columns: IID1,IID2
#                     These are related pairs of people.  
#
# Output : pd data frame with two columns (IID Status) 

def get_Pheno(pheno_file,related_people):
   
    # Load pheno file
    pheno = pd.read_csv(pheno_file,header=None,names=["FID","IID","Status"],delimiter=" ",dtype=str)  
    # Extract column of related people and turn into an array
    first = related_people['IID1'].to_numpy()
    second = related_people['IID2'].to_numpy()
    merged = np.concatenate((first, second), axis=None)
    unique_merged=np.unique(merged)
    # Convert array to dataframe
    list_related = pd.DataFrame(unique_merged, columns = ['IID'],dtype=str)
    # Merge
    joined = pd.merge(list_related,pheno,on='IID',how='left')
    joined_filter = joined[["IID","Status"]]
    
    return joined_filter

# Name: dict_Pheno 
# Description: Exporting phenotype file into a dictionary 
#              
#
# Parameters (2) : 1. phenotype (pd dataframe) : phenotype dataframe generated beforehand 
#                  2. case_value (string) : user-inputted value designating case
#
# Output : Dictionary : key is IID, value is phenotype status 

def dict_Pheno(phenotype,case_value):

    pheno_dict={} 
    # Processing Phenotype Dataframe, one row at a time. 
    for index, row in phenotype.iterrows():
        if (row['Status'] == case_value):
            result = int(1)
        elif (pd.isnull(row['Status'])):
            result = "NA"
        else:
            result = int(0)
        # Add to dictionary
        pheno_dict[row['IID']] = result

    return pheno_dict

# Name: calculate_score 
# Description: Calculates a Score for a case individual. For each
#              other case individual this person is related to, score +1.
#
# Parameters (3) : 1. ID (string) : IID of person 
#                  2. phenotype_dict (dictionary) : Dictionary Containing phenotype information 
#                  3. relatedness_df (pandas dataframe): Dataframe containing our related pairs 
#
# Output : Array : First element of array is score. Second element is how many relatives it has.

def calculate_score(ID,phenotype_dict,relatedness_df):
    #array to output
    output_array=[]
    
    #search up ID in the relatedness df and add to the list_related everyone this person is related to
    ID=str(ID)
    x=relatedness_df[relatedness_df.IID1 == ID]
    l1 = x["IID2"].tolist()
    x=relatedness_df[relatedness_df.IID2 == ID]
    l2 = x["IID1"].tolist()
    list_related=l1+l2
    
    score=0
    # Go through this list and calculate score
    for i in range(0,len(list_related)):
        case=phenotype_dict[str(list_related[i])]
        if (case==1):
            score=score+1
    
    # Prepare output array
    output_array.append(score)
    output_array.append(len(list_related))
    
    return (output_array)
    

# Name: calculate_score_controls
# Description: Calculates a Score for a control individual. For each
#              other control individual this person is related to, score +1.
#
# Parameters (3) : 1. ID (string) : IID of person 
#                  2. phenotype_dict (dictionary) : Dictionary Containing phenotype information 
#                  3. relatedness_df (pandas dataframe): Dataframe containing our related pairs 
#
# Output : Array : First element of array is score. Second element is how many relatives it has.

def calculate_score_controls(ID,phenotype_dict,relatedness_df):
    #array to output
    output_array=[]
    
    #search up ID in the relatedness df and add to the list_related everyone this person is related to
    ID=str(ID)
    x=relatedness_df[relatedness_df.IID1 == ID]
    l1 = x["IID2"].tolist()
    x=relatedness_df[relatedness_df.IID2 == ID]
    l2 = x["IID1"].tolist()
    list_related=l1+l2
    
    score=0
    # Go through this list and calculate score
    for i in range(0,len(list_related)):
        case=phenotype_dict[str(list_related[i])]
        if (case==0):
            score=score+1

    # Prepare output array
    output_array.append(score)
    output_array.append(len(list_related))
    
    return(output_array)
    

# Name: calculate_score_nas
# Description: Calculates a Score for an NA individual. For each
#              other na individual this person is related to, score +1.
#
# Parameters (3) : 1. ID (string) : IID of person 
#                  2. phenotype_dict (dictionary) : Dictionary Containing phenotype information 
#                  3. relatedness_df (pandas dataframe): Dataframe containing our related pairs 
#
# Output : Array : First element of array is score. Second element is how many relatives it has.

def calculate_score_nas(ID,phenotype_dict,relatedness_df):
    #array to output
    output_array=[]
    
    #search up ID in the relatedness df and add to the list_related everyone this person is related to
    ID=str(ID)
    x=relatedness_df[relatedness_df.IID1 == ID]
    l1 = x["IID2"].tolist()
    x=relatedness_df[relatedness_df.IID2 == ID]
    l2 = x["IID1"].tolist()
    list_related=l1+l2
    
    score=0
    # Go through this list and calculate score
    for i in range(0,len(list_related)):
        case=phenotype_dict[str(list_related[i])]
        if (case=="NA"):
            score=score+1

    # Prepare output array
    output_array.append(score)
    output_array.append(len(list_related))
    
    return(output_array)

def main(args):

    ### Output input paramters
    print("Starting Related Pair Case Prioritization for UKBB samples...")
    print(" ")
    print("~~~~~~~~~~~ Input Parameters ~~~~~~~~~~~~~~~~~")
    print(f'Pheno File: {args.pheno}')
    print(f'Kinship File: {args.kinship}')
    print(f'Samples File: {args.samples}')
    print(f'Case Value: {args.case_value}')
    print(f'Pihat Threshold: {args.pihat}')
    print(f'Output File: {args.output}')
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(" ")
    
    # Make kinship variable
    kinship = args.pihat/2

    # Check if phenotype file second column contains the correct values
    check_Pheno(args.pheno)

    # Load Samples and print total number of samples
    all_samples = pd.read_csv(args.samples,header=None,names=["FID","IID"],delimiter=" ")
    print(f'Total number of samples: {len(all_samples.index)}')
    
    # Generate a Table of the related People. Output a pandas dataframe with two columns IID1 and IID2 in which individuals are related 
    # to each other 
    related_people = get_Related(all_samples,args.kinship,kinship)
    print(f'Number of related pairs: {len(related_people.index)}')
    # Make sure both columns are strings
    related = related_people.applymap(str)
    

    # Generate a Table of people who are not related to anyone based on pihat value
    unrelated_strict = get_Strict_Unrelated(related_people,all_samples)
    print(f'Number of people who are not related to anyone else: {len(unrelated_strict)}')

    # Filter out phenotype file and make sure that all samples we are running the script on has phenotypes
    # If sample doesn't have phenotype, impute NA 
    final_pheno = get_Pheno(args.pheno,related_people)
    print(f'Number of people who have relatedness: {len(final_pheno.index)}')
 
    # Start Time Output
    localtime = time.asctime( time.localtime(time.time()) )
    print(" ")
    print("Starting Selection Step")
    print(f'Start Time: {localtime}')

    # Make a phenotype dictionary where key is ID and maps to phenotype status
    phenotype = dict_Pheno(final_pheno,args.case_value)

    # Establish list of cases, this is the list I use, and when this list is empty I'm done
    list_cases=[]
    for i, k in enumerate(phenotype):
        if (phenotype[k]==1):
            list_cases.append(k)
            
    # Establish list of controls.
    list_controls=[]
    for i, k in enumerate(phenotype):
        if (phenotype[k]==0):
            list_controls.append(k)

    # Establish list of NA's.
    list_nas=[]
    for i, k in enumerate(phenotype):
        if (phenotype[k]=="NA"):
            list_nas.append(k)


    # Print out Original List of cases, controls and NAs.
    print("")
    print(f'Original Number of Cases: {len(list_cases)}')
    print("")
    print(f'Original Number of Controls: {len(list_controls)}')
    print("")
    print(f'Original Number of NAs: {len(list_nas)}')
    print("")

    new_score=dict()
    # Go through list of cases, calculate score for the cases and include all cases that have a score of 0. Update
    # related file by removing these pairs. 
    for i in range(0,len(list_cases)):
        arr=calculate_score(list_cases[i],phenotype,related)
        x=int(arr[0])
        new_score[str(list_cases[i])]=x


    # Remove all scores of 0 and update related list (remove row where you see the person you removed). For list cases, remove all
    # scores of 0 and look at the people it was related to. Remove them if they are on the list cases list.
    # Add array is a list of people who we will eventually add to the list to use

    add_array=[]

    for i, k in enumerate(new_score):
        if (new_score[k]==0):
            add_array.append(k)
            list_cases.remove(k)
            
            ID=str(k)
            x=related[related.IID1 == ID]
            l1 = x["IID2"].tolist()
            x=related[related.IID2 == ID]
            l2 = x["IID1"].tolist()
            list_related=l1+l2
            
            #Remove those who are on related to the selection from list_controls or list_nas
            for i in range(0,len(list_related)):
                if list_related[i] in list_controls:
                    list_controls.remove(list_related[i])
                if list_related[i] in list_nas:
                    list_nas.remove(list_related[i])
            index = related[ related['IID1'] == k ].index
            related.drop(index , inplace=True)
            index = related[ related['IID2'] == k ].index
            related.drop(index , inplace=True)

    # List of global variables:

    # new_score=dict of scores of cases at the moment
    # add_array = list of people who can be used (output list)
    # phenotype = phenotype status 
    # related = related pairs status (This will be changed every iteration)
    # list_cases = people who are cases, calculated before. Will keep removing people. Algorithm finishes when everyone is removed from this list.
    # list_controls= people who are controls
    # list_nas = people who are nas

    # Use a while loop to loop through list_cases, until list_cases is empty.
    # Use the Scoring dictionary. Choose someone with the lowest score and add him to add array and remove him, updating related list,
    # list_cases. Update Scores. Repeat until noone is left in List Cases

    new_score=dict()
    num_relative=dict()

    # First, Calculate Scores:
    for i in range(0,len(list_cases)):
        arr=calculate_score(list_cases[i],phenotype,related)
        x=int(arr[0])
        new_score[str(list_cases[i])]=x
        num_relative[str(list_cases[i])]=int(arr[1])

    while(len(list_cases)>0):
        
        max_value=1000000000
        max_k=""
        max_rel=100000000000
        for i in range (0,len(list_cases)):
            k=list_cases[i]
            if (new_score[k]<=max_value and num_relative[k]<max_rel):
                max_value=new_score[k]
                max_k=k
                max_rel=num_relative[k]

        #Deal with eliminating max_k's pair
        ID=str(max_k)
        x=related[related.IID1 == ID]
        l1 = x["IID2"].tolist()
        x=related[related.IID2 == ID]
        l2 = x["IID1"].tolist()
        list_related=l1+l2

        for i in range(0,len(list_related)):
            
            #Remove those who are on related to the selection from both lists    
            if (phenotype[list_related[i]]==1):
                if list_related[i] in list_cases:
                    
                    # We need to access everyone in list_related and change their score
                    # They have 1 less point because one less case pair. Also have 1 less relative
                    # This way we don't have to recalculate the score for everyone!
                    new_score[str(list_related[i])]=new_score[str(list_related[i])]-1
                    num_relative[str(list_related[i])]=num_relative[str(list_related[i])]-1
                    list_cases.remove(list_related[i])
                    
            if (phenotype[list_related[i]]==0):
                if list_related[i] in list_controls:
                    list_controls.remove(list_related[i])
            if (phenotype[list_related[i]]=="NA"):
                 if list_related[i] in list_nas:
                    list_nas.remove(list_related[i])

        # eliminate max_k
        add_array.append(max_k)
        list_cases.remove(str(max_k))
        index = related[ related['IID1'] == max_k ].index
        related.drop(index , inplace=True)
        index = related[ related['IID2'] == max_k ].index
        related.drop(index , inplace=True)

    print(f'Number of Cases Selected: {len(add_array)}')

    # Deal with Controls 

    # Use a while loop to loop through list_controls, until list_controls is empty.
    # Use the Scoring dictionary. Choose someone with the lowest score and add him to second add array and remove him, updating related list,
    # list_controls. Update Scores. Repeat until noone is left in List Controls

    second_add_array=[]
    new_score=dict()
    num_relative=dict()
        
    # Go through list of cases, calculate score for the cases.Update related file by removing these pairs. 
    for i in range(0,len(list_controls)):
        arr=calculate_score_controls(list_controls[i],phenotype,related)
        x=int(arr[0])
        new_score[str(list_controls[i])]=x
        num_relative[str(list_controls[i])]=int(arr[1])

    while(len(list_controls)>0):
        
        max_value=1000000000000
        max_k=""
        max_rel=100000000000
        for i in range (0,len(list_controls)):
            k=list_controls[i]
            if (new_score[k]<=max_value and num_relative[k]<max_rel):
                max_value=new_score[k]
                max_k=k
                max_rel=num_relative[k]
                
        #Deal with eliminating max_k's pair
        ID=str(max_k)
        x=related[related.IID1 == ID]
        l1 = x["IID2"].tolist()
        x=related[related.IID2 == ID]
        l2 = x["IID1"].tolist()
        list_related=l1+l2
        
        for i in range(0,len(list_related)):
        
            #Remove those who are on related to the selection from both lists
            if (phenotype[list_related[i]]==1):
                if list_related[i] in list_cases:
                    list_cases.remove(list_related[i])
                    print("should never be ran")
            if (phenotype[list_related[i]]==0):
                if list_related[i] in list_controls:
                    list_controls.remove(list_related[i])
                    # We need to access everyone in list_related and change their score
                    # They have 1 less point because one less case pair. Also have 1 less relative
                    # This way we don't have to recalculate the score for everyone! 
                    new_score[str(list_related[i])]=new_score[str(list_related[i])]-1
                    num_relative[str(list_related[i])]=num_relative[str(list_related[i])]-1
            if (phenotype[list_related[i]]=="NA"):
                 if list_related[i] in list_nas:
                    list_nas.remove(list_related[i])

        # eliminate max_k
        second_add_array.append(max_k)
        list_controls.remove(str(max_k))
        index = related[ related['IID1'] == max_k ].index
        related.drop(index , inplace=True)
        index = related[ related['IID2'] == max_k ].index
        related.drop(index , inplace=True)

    print("")
    print(f'Number of Controls Selected: {len(second_add_array)}')
    print("")

    # Deal with NAs 

    # Use a while loop to loop through list_nas, until list_nas is empty.
    # Use the Scoring dictionary. Choose someone with the lowest score and add him to the third add array and remove him, updating related list,
    # list_nas. Update Scores. Repeat until noone is left in List Nas

    third_add_array=[]
    new_score=dict()
    num_relative=dict()
        
    # Go through list of cases, calculate score for the cases.Update related file by removing these pairs. 
    for i in range(0,len(list_nas)):
        arr=calculate_score_nas(list_nas[i],phenotype,related)
        x=int(arr[0])
        new_score[str(list_nas[i])]=x
        num_relative[str(list_nas[i])]=int(arr[1])

    while(len(list_nas)>0):
        
        max_value=1000000000000
        max_k=""
        max_rel=100000000000
        for i in range (0,len(list_nas)):
            k=list_nas[i]
            if (new_score[k]<=max_value and num_relative[k]<max_rel):
                max_value=new_score[k]
                max_k=k
                max_rel=num_relative[k]

        #Deal with eliminating max_k's pair
        ID=str(max_k)
        x=related[related.IID1 == ID]
        l1 = x["IID2"].tolist()
        x=related[related.IID2 == ID]
        l2 = x["IID1"].tolist()
        list_related=l1+l2
        
        for i in range(0,len(list_related)):
            
            #Remove those who are on related to the selection from both lists
            if (phenotype[list_related[i]]==1):
                if list_related[i] in list_cases:
                    list_cases.remove(list_related[i])
                    print("should never be ran")
            if (phenotype[list_related[i]]==0):
                if list_related[i] in list_controls:
                    list_controls.remove(list_related[i])
                    print("should never be ran")
            if (phenotype[list_related[i]]=="NA"):
                 if list_related[i] in list_nas:
                    list_nas.remove(list_related[i])
                    # We need to access everyone in list_related and change their score
                    # They have 1 less point because one less case pair. Also have 1 less relative
                    # This way we don't have to recalculate the score for everyone! 
                    new_score[str(list_related[i])]=new_score[str(list_related[i])]-1
                    num_relative[str(list_related[i])]=num_relative[str(list_related[i])]-1

        # eliminate max_k
        third_add_array.append(max_k)
        list_nas.remove(str(max_k))
        index = related[ related['IID1'] == max_k ].index
        related.drop(index , inplace=True)
        index = related[ related['IID2'] == max_k ].index
        related.drop(index , inplace=True)

    print(f'Number of NAs Selected: {len(third_add_array)}')
    print("")

    # Merge Arrays and Output:
    final_array=add_array+second_add_array+third_add_array

    # Output list:

    with open(args.output, 'w') as f:
        for i in range(0,len(final_array)):
            f.write(final_array[i])
            f.write('\n')    
        for index, row in unrelated_strict.iterrows():
            f.write(str(row['IID']))
            f.write('\n')    

    # Final Output List
    print(f'Final Output File N (Unrelateds chosen by Python Script merged with strict unrelateds): {len(final_array)+len(unrelated_strict.index)}')
    print(" ")

    # End Time Output:
    localtime = time.asctime( time.localtime(time.time()) )
    print(f'End Time: {localtime}')
 

# Working with arg parser
parser = argparse.ArgumentParser()
parser.add_argument('--pheno', help=' Phenotype file (FID IID Phenotype)', type=str,required=True)
parser.add_argument('--case_value', help='Which case indicator/value are we prioritizing', type=str,required=True)
parser.add_argument('--pihat', help='Pihat threshold', type=float,required=True)
parser.add_argument('--output', help='Output Name', type=str,required=True)
parser.add_argument('--kinship', help='Kinship File (FID IID Kinship)', type=str,required=True)
parser.add_argument('--samples', help='List of all samples (FID IID)', type=str,required=True)

if __name__ == '__main__':
	args = parser.parse_args()
	main(args)


