#!/usr/bin/env python

'''A script which returns the mutual information between the predictions of a
    model and a test data set.'''

from __future__ import division
#Our standard Modules
import argparse
import numpy as np
import scipy as sp
import sys
import pandas as pd
#Our miscellaneous functions
#This module will allow us to easily tally the letter counts at a particular position

import mpathic.utils as utils
import mpathic.EstimateMutualInfoforMImax as EstimateMutualInfoforMImax
import mpathic.qc as qc
import mpathic.numerics as numerics
from mpathic import SortSeqError
import mpathic.io as io
import matplotlib.pyplot as plt
import scipy


def main(
        data_df,model_df,
        start=0,end=None,err=False,coarse_graining_level=0,
        rsquared=False,return_freg=False):

    #determine whether you are working with RNA, DNA, or protein.
    #this also should determine modeltype (MAT, NBR, PAIR).
    dicttype, modeltype = qc.get_model_type(model_df)

    #get column header for the sequence column.
    seq_cols = qc.get_cols_from_df(data_df,'seqs')
    if not len(seq_cols)==1:
        raise SortSeqError('Dataframe has multiple seq cols: %s'%str(seq_cols))

    #create dictionary that goes from, for example, nucleotide to number and
    #visa versa.
    seq_dict,inv_dict = utils.choose_dict(dicttype,modeltype=modeltype)

    #set name of sequences column based on type of sequence
    type_name_dict = {'dna':'seq','rna':'seq_rna','protein':'seq_pro'}
    seq_col_name = type_name_dict[dicttype]

    if not end:
        seqL = len(data_df[seq_col_name][0]) - start
    else:
        seqL = end-start
    #throw out wrong length sequences.
    #Cut the sequences based on start and end, and then check if it makes sense 
    if (start != 0 or end):
        data_df.loc[:,seq_col_name] = \
            data_df.loc[:,seq_col_name].str.slice(start,end)
        right_length = data_df.loc[:,seq_col_name].apply(len) == (seqL)
        if not right_length.all():
            sys.stderr.write('''Not all sequences are the same length! 
                       Throwing out incorrect sequences!''')
            data_df = data_df.loc[right_length,:] 
        data_df = data_df.reset_index(drop=True)
        
        if modeltype =='MAT':
            if seqL != len(model_df.loc[:,'pos']):
                raise SortSeqError('model length does not match dataset length')
        elif modeltype =='NBR':
            if seqL != len(model_df.loc[:,'pos'])+1:
                raise SortSeqError('model length does not match dataset length')
        elif modeltype == 'PAIR':
            if int(scipy.misc.comb(seqL,2)) != len(model_df.loc[:,'pos']):
                raise SortSeqError('model length does not match dataset length')

    #get column names of the counts columns (excluding total counts 'ct')
    col_headers = utils.get_column_headers(data_df)
    if 'ct' not in data_df.columns:
                data_df['ct'] = data_df[col_headers].sum(axis=1)

    #remove empty rows.
    data_df = data_df[data_df.ct != 0]

    #determine sequence length.        
    


    #make a numpy array out of the model data frame
    model_df_headers = ['val_' + str(inv_dict[i]) for i in range(len(seq_dict))]
    value = np.array(model_df[model_df_headers])  

    #now we evaluate the expression of each sequence according to the model.
    #first convert to matrix representation of sequences
    seq_mat,wtrow = numerics.dataset2mutarray(data_df.copy(),modeltype)
    temp_df = data_df.copy()

    #evaluate energy of each sequence
    temp_df['val'] = numerics.eval_modelmatrix_on_mutarray(value,seq_mat,wtrow) 

    #sort based on value
    temp_sorted = temp_df.sort_values(by='val')
    temp_sorted.reset_index(inplace=True,drop=True)

    #freg is a regularized plot which show how sequences are distributed
    #in energy space.
    if return_freg:
        fig,ax = plt.subplots()   
        MI,freg = EstimateMutualInfoforMImax.alt4(temp_sorted,coarse_graining_level=coarse_graining_level,return_freg=return_freg)
        plt.imshow(freg,interpolation='nearest',aspect='auto')
        
        plt.savefig(return_freg)
    else:
        MI = EstimateMutualInfoforMImax.alt4(temp_sorted,coarse_graining_level=coarse_graining_level,return_freg=return_freg)

    #if we want to calculate error then use bootstrapping.
    if not err:
        Std = np.NaN
    else:
        data_df_for_sub = data_df.copy()
        sub_MI = np.zeros(15)
        for i in range(15):
            sub_df = data_df_for_sub.sample(int(len(data_df_for_sub.index)/2))
            sub_df.reset_index(inplace=True,drop=True)
            sub_MI[i],sub_std = main(
                sub_df,model_df,err=False)
        Std = np.std(sub_MI)/np.sqrt(2)

    #we can return linfoot corrolation (rsquared) or return MI.
    if rsquared:
        return (1-2**(-2*MI)),(1-2**(-2*Std))
    else:
        return MI,Std
     
def wrapper(args):
    
    data_df = io.load_dataset(args.dataset)    	    
    # Take input from standard input or through the -i flag.
    if args.model:
        model_df = io.load_model(args.model)
    else:
        model_df = io.load_model(sys.stdin)

    MI,Std = main(
        data_df,model_df,start=args.start,
        end=args.end,err=args.err,coarse_graining_level = args.coarse_graining_level,
        rsquared=args.rsquared,return_freg=args.return_freg)
    
    #format output
    output_df = pd.DataFrame([MI],columns=['info'])
    
    #if you calculated error add column to your data frame
    if args.err:
        output_df = pd.concat([output_df,pd.Series(Std,name='info_err')],axis=1)
  
    if args.out:
        outloc = open(args.out,'w')
    else:
        outloc = sys.stdout
    
    #set output option, this will remove column length restriction
    pd.set_option('max_colwidth',int(1e8))
    #write to file.
    output_df.to_string(
        outloc, index=False,col_space=10,float_format=utils.format_string)

# Connects argparse to wrapper
def add_subparser(subparsers):
    p = subparsers.add_parser('predictiveinfo')
    p.add_argument('-rs','--rsquared',action='store_true',help='return effective r squared')
    p.add_argument('-ds','--dataset')
    p.add_argument(
        '--err',action='store_true',help='''Flag to use if you want to
        calculate error''')
    p.add_argument(
        '-s','--start',type=int,default=0,help ='''Position to start your 
        analyzed region''')
    p.add_argument(
        '-e','--end',type=int,default = None, 
        help='''Position to end your analyzed region''')
    p.add_argument(
        '-fr','--return_freg',type=str,
        help='''return regularized plot and save it to this file name''')
    p.add_argument(
        '-cg','--coarse_graining_level',default=0,type=float,help='''coarse graining
        level to use for mutual information calculation, higher values will
        speed up computation''')
    p.add_argument(
        '-m', '--model', default=None,help='''Model file, otherwise input
        through the standard input.''')            
    p.add_argument('-o', '--out', default=None)
    p.set_defaults(func=wrapper)
