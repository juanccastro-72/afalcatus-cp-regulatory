#!/usr/bin/env python3
"""Ancestral-state reconstruction (Fitch parsimony) of the TATAATAT −10 octamer
promoter signature across the plastome ML tree, to test where the signature was lost."""
import argparse, re
from pathlib import Path
import numpy as np, pandas as pd
from Bio import Phylo
import figstyle as fs
from figstyle import plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from figures import FOCAL, tree_layout


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'; FM = P/'figures/main'
    import matplotlib.font_manager as _fm
    for _p in ['/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf',
               '/usr/share/fonts/truetype/liberation/LiberationSerif-Italic.ttf',
               '/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf']:
        try: _fm.fontManager.addfont(_p)
        except Exception: pass
    SERIF = {'font.family':'serif','font.serif':['Liberation Serif','Times New Roman','DejaVu Serif'],'svg.fonttype':'none'}
    S = pd.read_csv(T/'genome_level_signal_summary.tsv', sep='\t')
    man = pd.read_csv(T/'genome_manifest.tsv', sep='\t')
    sp = dict(zip(man.insdc_accession, man.species)); grp = dict(zip(S.genome, S.group))
    ACC_VER = {a2.split('.')[0]: a2 if '.' in a2 else a2+'.1' for a2 in man.insdc_accession}
    enr = {r.genome: int((r.oct_P_enriched < 0.05) and (r.oct_OE > 1)) for r in S.itertuples()}
    oe = dict(zip(S.genome, S.oct_OE))
    tree = Phylo.read(P/'results/phylogeny/plastome_ML.treefile','newick')
    accmap = {x.split('.')[0]: x for x in man.insdc_accession}
    tip_acc = lambda name: accmap[name.rsplit('_',1)[1]]
    # ---- Fitch parsimony ----
    up = {}
    def post(cl):
        if cl.is_terminal():
            up[cl] = {enr[tip_acc(cl.name)]}; return up[cl]
        sets = [post(k) for k in cl.clades]
        inter = set.intersection(*sets)
        up[cl] = inter if inter else set.union(*sets); return up[cl]
    post(tree.root)
    state = {}
    def pre(cl, par):
        if par is None:
            state[cl] = 1 if up[cl]=={1} else (0 if up[cl]=={0} else 1)  # root: default to enriched if ambiguous (ancestor)
        else:
            state[cl] = state[par] if state[par] in up[cl] else min(up[cl])
        for k in cl.clades: pre(k, cl)
    pre(tree.root, None)
    # count transitions
    losses, gains = [], []
    def walk(cl):
        for k in cl.clades:
            if state[k] != state[cl]:
                (losses if state[cl]==1 and state[k]==0 else gains).append(k)
            walk(k)
    walk(tree.root)
    # report
    rows=[]
    for cl in tree.find_clades():
        if cl.is_terminal():
            acc=tip_acc(cl.name)
            rows.append(dict(node=sp[acc], accession=ACC_VER[acc.split('.')[0]], group=grp[acc],
                             octamer_OE=round(oe[acc],2), enriched=enr[acc], reconstructed=state[cl]))
    pd.DataFrame(rows).to_csv(T/'octamer_ancestral_states.tsv', sep='\t', index=False)
    # ---- internal-node reconstructed states (Pt 8: internal-node ancestral reconstructions) ----
    node_id = {}
    for i, cl in enumerate(tree.get_nonterminals()):
        node_id[cl] = f'N{i:02d}'
    irows = []
    for cl in tree.get_nonterminals():
        tips_here = [sp[tip_acc(t.name)] for t in cl.get_terminals()]
        irows.append(dict(internal_node=node_id[cl],
                          n_descendant_tips=len(tips_here),
                          reconstructed_state=('enriched' if state[cl]==1 else 'not_enriched'),
                          fitch_uppass_set='/'.join(str(x) for x in sorted(up[cl])),
                          is_root=(cl is tree.root),
                          example_descendants='; '.join(tips_here[:5]) + ('...' if len(tips_here)>5 else '')))
    pd.DataFrame(irows).to_csv(T/'octamer_ancestral_internal_nodes.tsv', sep='\t', index=False)
    # ---- sensitivity recoding (Pt 5: do not equate non-significant with absence) ----
    srows = []
    for g in S.genome:
        oeg = oe[g]; pv = float(S.loc[S.genome==g,'oct_P_enriched'].iloc[0])
        if pv < 0.05 and oeg > 1: code = 'enriched (P<0.05, O/E>1)'
        elif oeg >= 1: code = 'not significantly enriched (O/E>=1, P>=0.05)'
        else: code = 'below expectation (O/E<1)'
        srows.append(dict(genome=g, species=sp.get(g,g), octamer_OE=round(oeg,2),
                          P_enriched=round(pv,4), detection_state=code))
    pd.DataFrame(srows).to_csv(T/'octamer_detection_states_sensitivity.tsv', sep='\t', index=False)
    def desc_tips(cl):
        return [sp[tip_acc(t.name)] for t in cl.get_terminals()]
    summary=[f"Tips coded as octamer-enriched: {sum(enr.values())}/{len(enr)}",
             f"Inferred losses of the octamer signature: {len(losses)}",
             f"Inferred (re)gains: {len(gains)}",
             f"Root state (chlorophycean ancestor): {'enriched' if state[tree.root]==1 else 'not enriched'}"]
    for k,ls in (('LOSS',losses),('GAIN',gains)):
        for cl in ls:
            tx=desc_tips(cl); summary.append(f"{k} on branch to clade ({len(tx)} tips): {', '.join(tx[:6])}{'...' if len(tx)>6 else ''}")
    (T/'octamer_ancestral_summary.txt').write_text('\n'.join(summary))
    print('\n'.join(summary))
    # ---- figure ----
    tips, ty, td = tree_layout(tree)
    with plt.rc_context(SERIF):
        fig = plt.figure(figsize=(fs.SINGLE*2.0, fs.DOUBLE*0.82))
        ax = fig.add_axes([0.0,0.06,0.62,0.92])
        ENR_C, NO_C = '#2a78d6', '#c0c0c0'
        for c in tree.find_clades():
            for k in c.clades:
                cc = ENR_C if state[k]==1 else NO_C
                ax.plot([td[c],td[c]],[ty[c],ty[k]], color=cc, lw=1.4, solid_capstyle='round')
                ax.plot([td[c],td[k]],[ty[k],ty[k]], color=cc, lw=1.4, solid_capstyle='round')
        xmax=max(td.values())
        for cl in losses:
            ax.plot(td[cl], ty[cl], marker='X', ms=6, color='#c0392b', mec='white', mew=0.6, zorder=6)
        for t in tips:
            acc=tip_acc(t.name); fcl=(acc==FOCAL)
            ax.plot(td[t]+0.004*xmax, ty[t], 'o', ms=4.2, mec=ENR_C if enr[acc] else NO_C,
                    mfc=(ENR_C if enr[acc] else 'white'), mew=0.9)
            ax.text(td[t]+0.016*xmax, ty[t], f'{sp[acc]}', va='center', fontsize=5.0,
                    fontstyle='italic', fontweight='bold' if fcl else 'normal',
                    color=fs.FOCAL if fcl else (fs.GROUP_COLOR[grp[acc]]))
        ax.set_ylim(len(tips)+0.6,-0.6); ax.set_xlim(0,xmax*1.95); ax.axis('off')
        ax.plot([0.01,0.11],[len(tips)]*2,color=fs.INK,lw=0.8); ax.text(0.06,len(tips)+0.4,'0.1 subst./site',fontsize=5,ha='center',va='top')
        leg=[Line2D([],[],marker='o',ls='',mfc=ENR_C,mec=ENR_C,label='Octamer-enriched (tip)'),
             Line2D([],[],marker='o',ls='',mfc='white',mec=NO_C,label='Not enriched (tip)'),
             Line2D([],[],color=ENR_C,lw=2,label='Branch reconstructed: enriched'),
             Line2D([],[],color=NO_C,lw=2,label='Branch reconstructed: not enriched'),
             Line2D([],[],marker='X',ls='',color='#c0392b',mec='white',label='Inferred loss of octamer signature')]
        ax.legend(handles=leg, loc='lower right', fontsize=5.6, frameon=False, bbox_to_anchor=(1.98,0.0))
        fs.save(fig, FM, 'Fig9_octamer_ancestral_reconstruction', [],
                'Fitch-parsimony ancestral reconstruction of the TATAATAT −10 octamer promoter signature')
    print('Fig9 done')

if __name__=='__main__':
    main()
