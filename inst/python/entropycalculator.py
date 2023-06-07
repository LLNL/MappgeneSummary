#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 14:30:26 2023

@author: phoebeadamyan
"""

import argparse
import pandas as pd
import math 

#%%
#DEFINE FUNCTIONS

def filter_input_dataframe(df):
    
    #isolate columns i.e.keep the ones i want  
    df = df.loc[:,['POS','REF','ALT','AF']]
    #remove all strings with more than one base: ALT and REF
    df = df[(df['REF'].str.len() == 1)&(df['ALT'].str.len() == 1)]
    #drop duplicates
    df = df.drop_duplicates(subset=['POS', 'ALT'])
    
    return df


# new function creating a dictionary of alt and ref
def convert_dataframe_to_dictionary(DF):
    POS_dict = {}

    for index,row in DF.iterrows():
        P=row['POS']
    
    #if position already exists, adding new row 
        if P in POS_dict:
            POS_dict[P]['F'][row['ALT']]=row['AF']
            
        else:
            POS_dict[P] = {'REF':row['REF'], 
                              'F':{row['ALT']:row['AF']}}
        
    for key,value in POS_dict.items():
        SUM_ALT = sum(value['F'].values())
        REF = 1 - SUM_ALT
    
        POS_dict[key]['F'][value['REF']] = REF
        
    return POS_dict

    
    
# Entropy function
def S_pos(D):  
    
    Pos_sum = 0
    
    for i,AF in D.items():
        if AF > 0:
            Pos_sum += AF*math.log(AF)
            
    return(-Pos_sum)

#Average 
def S_avg(S_N, N=29903):
    S_sum = 0
    
    for S in S_N:
        S_sum += S
    
    return(S_sum/N)

#%%
#MAIN

#takes whatever is put in the command line and converts it to string
def main():
    parser = argparse.ArgumentParser(description='Process data frame')
    parser.add_argument('--dataframe', type=str, help='path to data frame file')
    parser.add_argument('--outputdata', type=str, help='data transferred to csv file')

    args = parser.parse_args()

    df = pd.read_table(args.dataframe)

    df = filter_input_dataframe(df)
    POS_dict = convert_dataframe_to_dictionary(df)

    F = {}

    for P,v in POS_dict.items():
        F[P] = S_pos(v['F'])
        
    df = pd.DataFrame(F.items(),columns=['POS','Entropy'])
    
    df.to_csv(args.outputdata, index=False)

#%%

if __name__ == '__main__':
    main()
    
