# -*- coding: utf-8 -*-

import gc
import pandas as pd
import numpy as np
import re
import scipy.stats as stats
import matplotlib.pyplot as plt


data=pd.read_csv('data/Raw Data_GeneSpring.txt',delimiter='\t')

pattern=r'GSM\d*_(M|F)_(Ns|Sm|SM).*'

# determining the number of column of the form 'GSM..'
count=0
for col in data.columns:
  if re.match(pattern,col):
    count+=1

gc.collect()

N=np.zeros((count,4))
D=np.zeros((count,4))
stat_dict = {
    'row_index': [],
    'f_stat': [],
    'p_value': []
}

for row_index in range(len(data.index)): #data.loc[row_index] is the 'row_index'th row
  col_index=0
  X=np.array(data.loc[row_index][1:-3])
  for col in data.loc[row_index].index: #Creating D,N,X matrices
    if re.search(pattern, col):
      gender=re.search(pattern, col).group(1)
      smoker=re.search(pattern, col).group(2)
      if gender=='M':
        N[col_index][0]=1
        if smoker=='Ns':
          N[col_index][3]=1
          D[col_index][1]=1
        elif smoker=='SM' or smoker=='Sm':
          N[col_index][2]=1
          D[col_index][0]=1
      elif gender=='F':
        N[col_index][1]=1
        if smoker=='Ns':
          N[col_index][3]=1
          D[col_index][3]=1
        elif smoker=='SM' or smoker=='Sm':
          N[col_index][2]=1
          D[col_index][2]=1
      
      col_index+=1   

  rank_D = np.linalg.matrix_rank(D)  
  rank_N=  np.linalg.matrix_rank(N)
  r=(1/(rank_D-rank_N))/(1/(count-rank_D))
  N_=N@(np.linalg.pinv(((N.T)@N)))@(N.T)
  D_=D@(np.linalg.pinv(((D.T)@D)))@(D.T)
  row_N,col_N=N_.shape
  num=(X.T)@(np.eye(count)-N_)@X
  row_D,col_D=D_.shape
  den=(X.T)@(np.eye(count)-D_)@X 
  df1 = rank_D - rank_N
  df2 = count - rank_D
  if den != 0:
    stat_dict['f_stat'].append(r*((num/den)-1))
    stat_dict['row_index'].append(row_index)
    stat_dict['p_value'].append(1 - stats.f.cdf(stat_dict['f_stat'][-1], df1, df2))
  else:
    continue
  if row_index%300==0:  
    gc.collect()

stat_dict['p_value']

# Number of rows with p-value <= 5%
gene=[]
for i in range(len(stat_dict['row_index'])):
  if stat_dict['p_value'][i]<=0.05:
    gene.append(i)
fiveper_pvalue=len(gene)

gc.collect()
def plot_me(nobins):
    plt.hist(np.array(stat_dict['p_value'],dtype=object),edgecolor='black',bins=nobins)
    plt.title('Histogram of p-values')
    plt.xlabel('p-values')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()
    plt.savefig("histogram")

gc.collect()
plot_me(10)
gc.collect()
plot_me(20)
gc.collect()
plot_me(100)
gc.collect()
plot_me(1000)
gc.collect()
plot_me(10000)
gc.collect()
plot_me(len(stat_dict['row_index']))

