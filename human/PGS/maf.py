# Title : Plotting MAF for PGS variants of interest
# Author : Dr. Alice M. Godden



import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
from matplotlib import rcParams

# bold fonts
rcParams['font.weight']='bold'
rcParams['axes.labelweight']='bold'
rcParams['axes.titleweight']='bold'
rcParams['xtick.labelsize']=14
rcParams['ytick.labelsize']=14

xls=pd.ExcelFile('panukbiobank_matches.xlsx')
su=pd.read_excel(xls,'minor allele frequencies_SU')
mc=pd.read_excel(xls,'minor allele frequencies_MC')

def extract(df):
    maf_c=[]; maf_o=[]
    for s in df['info-MAF'].astype(str):
        m=re.search(r'Maf_C=([0-9\\.]+)', s)
        n=re.search(r'Maf_O=([0-9\\.]+)', s)
        maf_c.append(float(m.group(1)) if m else np.nan)
        maf_o.append(float(n.group(1)) if n else np.nan)
    df['Maf_C']=maf_c; df['Maf_O']=maf_o
    df=df.dropna(subset=['Maf_C','Maf_O'])
    return df

su=extract(su)
mc=extract(mc)

def bubble_plot(df,name):
    genes=df['nearest_genes']
    x=np.arange(len(df))
    size=(df['Maf_C']+df['Maf_O'])/2*500  # bubble size
    plt.figure(figsize=(12,6))
    plt.scatter(x, df['Maf_C'], s=size, color='blue', alpha=0.6, label='Maf_Centre')
    plt.scatter(x, df['Maf_O'], s=size, color='red', alpha=0.6, label='Maf_Outer')
    plt.xticks(x, genes, rotation=45, fontweight='bold')
    plt.ylabel('MAF', fontweight='bold')
    # Fix y‑axis 0–1
    plt.ylim(0, 1)
    plt.title(f'{name} MAF score for variants', fontweight='bold')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'bubble_{name}.png')

bubble_plot(su,'SU')
bubble_plot(mc,'MC')

