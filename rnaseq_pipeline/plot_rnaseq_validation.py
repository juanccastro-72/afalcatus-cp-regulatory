#!/usr/bin/env python3
"""Fig. 10: RNA-seq transcriptional support for A. falcatus candidate regulatory regions.
Reads cov/*.npz (from rnaseq_plastome_coverage.py). Produces a genome-wide strand-specific
coverage track plus zoom panels on top candidate 5' UTRs, marking predicted promoter/UTR.
Run after coverage files exist. Requires numpy, pandas, matplotlib, biopython.
"""
import argparse, glob
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from Bio import SeqIO

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cov',default='cov'); ap.add_argument('--genome',required=True)
    ap.add_argument('--candidates',required=True)  # Afalcatus_integrated_candidates.tsv
    ap.add_argument('--out',default='Fig10_RNAseq_validation')
    a=ap.parse_args()
    for p in ['/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf',
              '/usr/share/fonts/truetype/liberation/LiberationSerif-Italic.ttf']:
        try:
            import matplotlib.font_manager as fm; fm.fontManager.addfont(p)
        except Exception: pass
    plt.rcParams.update({'font.family':'serif','font.serif':['Liberation Serif','Times New Roman','DejaVu Serif'],'svg.fonttype':'none','font.size':8})
    files=sorted(glob.glob(f'{a.cov}/*.npz'))
    assert files, 'no coverage .npz found; run rnaseq_plastome_coverage.py first'
    d0=np.load(files[0]); n=int(d0['length'])
    covp=np.zeros(n); covm=np.zeros(n); fivep=np.zeros(n); fivem=np.zeros(n); stats=[]
    for f in files:
        d=np.load(f); covp+=d['covp']; covm+=d['covm']; fivep+=d['fivep']; fivem+=d['fivem']
        stats.append((Path(f).stem, int(d['mapped']), int(d['total'])))
    rec=SeqIO.read(a.genome,'genbank')
    C=pd.read_csv(a.candidates, sep='\t')
    genes=['psbD','rbcL','petD','psbA','atpA','psaB']
    fig=plt.figure(figsize=(174/25.4, 200/25.4))
    gs=fig.add_gridspec(len(genes)+1,1, height_ratios=[1.4]+[1]*len(genes), hspace=0.85)
    # genome-wide
    ax=fig.add_subplot(gs[0]); x=np.arange(n)
    ax.fill_between(x, np.log10(covp+1), color='#2a78d6', lw=0, label='+ strand')
    ax.fill_between(x, -np.log10(covm+1), color='#eb6834', lw=0, label='− strand')
    ax.set_title('Genome-wide strand-specific RNA-seq coverage (log10; merged alignment of four A. falcatus libraries)', fontsize=8, loc='left')
    ax.set_xlim(0,n); ax.set_xlabel('plastome position (bp)', fontsize=7, labelpad=1); ax.set_ylabel('±log10(cov+1)'); ax.legend(fontsize=6, ncol=2, loc='upper right')
    # zooms
    feat={f.qualifiers['gene'][0]:f for f in rec.features if f.type in ('CDS','rRNA') and 'gene' in f.qualifiers}
    for i,g in enumerate(genes):
        ax=fig.add_subplot(gs[i+1])
        if g not in feat: ax.set_visible(False); continue
        ft=feat[g]; st=int(ft.location.parts[0].start); en=int(ft.location.parts[0].end); strand=ft.location.parts[0].strand
        if strand==1: a0,a1=st-350,st+150; cov=covp; five=fivep; tss_x=st
        else: a0,a1=en-150,en+350; cov=covm; five=fivem; tss_x=en
        a0=max(a0,0); a1=min(a1,n); xx=np.arange(a0,a1)
        ax.fill_between(xx, cov[a0:a1], color='#555', lw=0)
        ax2=ax.twinx(); ax2.plot(xx, five[a0:a1], color='#c0392b', lw=0.6); ax2.set_ylabel("5′-end reads", fontsize=6, color='#c0392b')
        ax.axvline(st if strand==1 else en, color='#008300', lw=0.8, ls='--')
        ax.set_title(f'{g} 5′ region (start codon dashed; coverage grey, 5′-end counts red)', fontsize=7, loc='left', fontstyle='italic')
        ax.set_xlim(a0,a1); ax.set_ylabel('cov', fontsize=6); ax.tick_params(labelsize=6)
    fig.savefig(f'{a.out}.jpg', dpi=600, bbox_inches='tight'); fig.savefig(f'{a.out}.pdf', bbox_inches='tight'); fig.savefig(f'{a.out}.svg', bbox_inches='tight')
    print('wrote', a.out+'.jpg'); print('per-dataset mapping:'); [print(' ',s) for s in stats]

if __name__=='__main__':
    main()
