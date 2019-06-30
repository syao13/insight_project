# data analysis and wrangling
import pandas as pd
import numpy as np
import sys
from scipy.stats import pearsonr

# import data
data = pd.read_csv("sanger1018_brainarray_ensemblgene_rma.txt", sep='\t')

def calc_correlaton(i):
    for i in range(i, data.shape[1] - 1):
        for j in range(i+1, data.shape[1]):
            rsquared = pearsonr(data.loc[i,:].tolist()[1:], data.loc[j,:].tolist()[1:])[0]
            if rsquared >= 0.5 or rsquared <= -0.5:
                print(",".join([str(i) for i in [i, j, data.iloc[i,0], data.iloc[j,0], rsquared]]))
                        
if __name__ == '__main__':
    calc_correlaton(int(sys.argv[1]))




