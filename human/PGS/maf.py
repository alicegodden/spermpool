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


# histogram


import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
from matplotlib import rcParams

# ---------------------------------------------------------
#  FONT / STYLE
# ---------------------------------------------------------
rcParams['font.weight'] = 'bold'
rcParams['axes.labelweight'] = 'bold'
rcParams['axes.titleweight'] = 'bold'
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14

# ---------------------------------------------------------
#  LOAD DATA
# ---------------------------------------------------------
xls = pd.ExcelFile('panukbiobank_matches.xlsx')
su = pd.read_excel(xls, 'minor allele frequencies_SU')
mc = pd.read_excel(xls, 'minor allele frequencies_MC')

# ---------------------------------------------------------
#  EXTRACT Maf_C / Maf_O FROM 'info-MAF'
# ---------------------------------------------------------
def extract(df):
    maf_c, maf_o = [], []
    for s in df['info-MAF'].astype(str):
        m = re.search(r'Maf_C=([0-9\.]+)', s)
        n = re.search(r'Maf_O=([0-9\.]+)', s)
        maf_c.append(float(m.group(1)) if m else np.nan)
        maf_o.append(float(n.group(1)) if n else np.nan)
    df = df.copy()
    df['Maf_C'] = maf_c
    df['Maf_O'] = maf_o
    # keep rows where both values are present
    df = df.dropna(subset=['Maf_C', 'Maf_O']).reset_index(drop=True)
    return df

su = extract(su)
mc = extract(mc)

# ---------------------------------------------------------
#  GROUPED BAR PLOT: Maf_C vs Maf_O per gene
# ---------------------------------------------------------
def bar_per_gene_grouped(df, name, sort_by=None, ascending=False):
    """
    Makes a grouped bar chart for one dataset (SU or MC):
      - x-axis: genes
      - bars: Maf_C (blue) and Maf_O (red)
    Args:
        df: dataframe with columns ['nearest_genes', 'Maf_C', 'Maf_O']
        name: label for title and filename ('SU' or 'MC')
        sort_by: optional, one of ['Maf_C','Maf_O','delta'] to sort bars
        ascending: sort order
    """
    d = df.copy()

    # Optional sorting
    if sort_by == 'delta':
        d['delta'] = d['Maf_C'] - d['Maf_O']
        d = d.sort_values('delta', ascending=ascending)
    elif sort_by in ('Maf_C', 'Maf_O'):
        d = d.sort_values(sort_by, ascending=ascending)

    d = d.reset_index(drop=True)
    genes = d['nearest_genes'].astype(str)
    x = np.arange(len(d))
    width = 0.42

    # Auto-size width for long gene lists
    fig_width = max(12, len(d) * 0.8)
    plt.figure(figsize=(fig_width, 6))

    # Bars
    plt.bar(x - width/2, d['Maf_C'], width, color='blue', alpha=0.8, edgecolor='black', label='Maf_C')
    plt.bar(x + width/2, d['Maf_O'], width, color='red',  alpha=0.8, edgecolor='black', label='Maf_O')

    # Axes & formatting
    plt.xticks(x, genes, rotation=45, ha='right', fontweight='bold')
    plt.ylabel('MAF', fontweight='bold')
    plt.ylim(0, 1)
    plt.title(f'{name}: MAF per gene (Maf_C vs Maf_O)', fontweight='bold')
    plt.legend()
    plt.tight_layout()

    # Save & show
    out = f'bar_{name}_mafc_difference_grouped.png'
    plt.savefig(out, dpi=300)
    plt.show()
    plt.close()

# ---------------------------------------------------------
#  PLOTS
# ---------------------------------------------------------
# Simple, unsorted:
#bar_per_gene_grouped(su, 'SU')
#bar_per_gene_grouped(mc, 'MC')

# If you prefer sorted by Maf_C descending, uncomment:
#bar_per_gene_grouped(su, 'SU', sort_by='Maf_C', ascending=False)
#bar_per_gene_grouped(mc, 'MC', sort_by='Maf_C', ascending=False)

# Or sort by the difference (Maf_C - Maf_O):
bar_per_gene_grouped(su, 'SU', sort_by='delta', ascending=False)
bar_per_gene_groupe
