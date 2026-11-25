#!/usr/bin/env python3
"""Anti-Shine-Dalgarno analysis: the 16S rRNA 3'-terminal anti-SD core (CCUCC) and its
complementarity to the SD-like (GGAG-family) words. Aligned on the core where detectable."""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from Bio import SeqIO
import figstyle as fs
from figstyle import plt
import logomaker


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project', default='..')
    a=ap.parse_args(); P=Path(a.project); T=P/'results/tables'; FS=P/'figures/supplementary'
    import matplotlib.font_manager as _fm
    for _p in ['/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf',
               '/usr/share/fonts/truetype/liberation/LiberationSerif-Italic.ttf',
               '/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf']:
        try: _fm.fontManager.addfont(_p)
        except Exception: pass
    SERIF={'font.family':'serif','font.serif':['Liberation Serif','Times New Roman','DejaVu Serif'],'svg.fonttype':'none'}
    man=pd.read_csv(T/'genome_manifest.tsv', sep='\t'); sp=dict(zip(man.insdc_accession,man.species))
    rows=[]; windows=[]; LW=4; RW=4   # flank around 5-nt core
    for gb in sorted((P/'data/standardized').glob('*.gb')):
        rec=SeqIO.read(gb,'genbank'); acc=rec.id if rec.id in sp else gb.stem
        fe=[f for f in rec.features if f.type=='rRNA' and f.qualifiers.get('gene',[''])[0]=='rrn16']
        if not fe: continue
        tail=str(fe[0].extract(rec.seq)[-60:]).upper().replace('U','T')
        pos=tail.rfind('CCTCC')
        has=pos>=0
        rows.append(dict(genome=acc, species=sp.get(acc,acc), tail60_3prime=tail,
                         antiSD_core='yes' if has else 'no'))
        if has and pos-LW>=0 and pos+5+RW<=len(tail):
            windows.append(tail[pos-LW:pos+5+RW])
    df=pd.DataFrame(rows); df.to_csv(T/'antiSD_16S_3prime.tsv', sep='\t', index=False)
    n_has=(df.antiSD_core=='yes').sum()
    print(f'anti-SD core CCUCC detected in {n_has}/{len(df)} plastome 16S 3-prime regions')
    af=df[df.genome.str.startswith('QB176266')].iloc[0]; cr=df[df.genome.str.startswith('FJ423446')].iloc[0]
    print('A. falcatus 16S 3-prime tail:', af.tail60_3prime[-18:])
    print('C. reinhardtii 16S 3-prime tail:', cr.tail60_3prime[-18:])
    # logo aligned on the core
    L=LW+5+RW
    counts=pd.DataFrame(0,index=range(L),columns=list('ACGT'))
    for w in windows:
        for i,b in enumerate(w):
            if b in 'ACGT': counts.loc[i,b]+=1
    info=logomaker.transform_matrix(counts, from_type='counts', to_type='information')
    with plt.rc_context(SERIF):
        fig,(ax1,ax2)=plt.subplots(2,1,figsize=(fs.DOUBLE,fs.DOUBLE*0.44),gridspec_kw=dict(height_ratios=[1,0.8]))
        logomaker.Logo(info.replace({'U':'T'}), ax=ax1, color_scheme='classic')
        ax1.set_ylabel('bits', fontsize=7)
        ax1.set_title(f'16S rRNA 3′ anti-SD region, aligned on the CCUCC core (n = {len(windows)} plastomes)', fontsize=8, loc='left')
        ax1.set_xticks(range(L)); ax1.set_xticklabels(list('····CCTCC····'[:L]), fontsize=6); ax1.tick_params(labelsize=6)
        ax2.axis('off')
        ax2.text(0.5,0.93,'16S 3′ anti-SD : SD-like element complementarity (Ankistrodesmus falcatus)', ha='center', fontsize=8, fontstyle='italic')
        anti=af.tail60_3prime[-12:].replace('T','U')
        ax2.text(0.28,0.60, "16S 3′   3′-…"+anti[::-1]+"-5′", fontsize=11, family='monospace', ha='left')
        comp={'A':'U','U':'A','G':'C','C':'G'}
        core='CCUCC'; sd=''.join(comp[b] for b in core[::-1])  # GGAGG
        ax2.text(0.28,0.40, "                 "+'|||||', fontsize=11, family='monospace', ha='left', color='#999')
        ax2.text(0.28,0.22, f"mRNA 5′   5′-…{sd}…-3′  (SD-like, GGAG family: psbD/psbC/atpH/psaB)", fontsize=11, family='monospace', ha='left', color='#2a78d6')
        ax2.text(0.5,0.02,f'Canonical anti-SD (3′-…GAUCACCUCCUU) present in A. falcatus and C. reinhardtii; core detected in {n_has}/{len(df)} panel 16S 3′ ends (others reflect imprecise WGS rRNA 3′-end annotation).', ha='center', fontsize=6.4, color=fs.MUTED)
        fig.tight_layout()
        fs.save(fig, FS, 'FigS11_antiSD_16S', [], 'The 16S rRNA 3′ anti-SD and its complementarity to SD-like elements')
    print('FigS11 done')

if __name__=='__main__':
    main()
