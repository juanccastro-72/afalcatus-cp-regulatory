"""Excel and publication-resolution figures for the plastid candidate inventory."""
from pathlib import Path
import textwrap
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

COLORS={'PEP_paired_boxes':'#285E8E','algal_TATAATAT':'#25948A','isolated_TATAAT':'#98BCB4',
        'SD_like':'#D69839','downstream_hairpin':'#8D609A'}
LABELS={'PEP_paired_boxes':'PEP paired boxes','algal_TATAATAT':'Algal TATAATAT',
        'isolated_TATAAT':'Exact TATAAT','SD_like':'SD-like','downstream_hairpin':"3′ hairpins"}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,
    'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,
    'svg.fonttype':'none','svg.hashsalt':'chloroplast-inventory-20261023',
    'savefig.facecolor':'white','axes.titleweight':'bold'})


def add_summary(data):
    meta=dict(zip(data['Metadata'].field,data['Metadata'].value))
    rows=[]
    def add(section,metric,value,note=''):
        rows.append({'section':section,'metric':metric,'value':value,'interpretation':note})
    for key in ['accession','organism','genome_length_bp','GC_percent','n_product_features','n_distinct_gene_names','n_boundary_contexts','n_mixed_strand_products','IR_reverse_complement_identity_percent']:
        add('Genome',key,meta.get(key,''))
    for typ,g in data['Products'].groupby('feature_type'):add('Annotations',typ+'_copies',len(g))
    inv=data['Inventory_associations'];sites=data['Unique_sites']
    for cat in COLORS:
        g=inv[inv.element_type==cat]
        add('Candidates',cat+'_associations',len(g),'Candidate-to-context associations, not independent functional elements.')
        add('Candidates',cat+'_unique_sites',int((sites.element_type==cat).sum()),'Unique within this motif/structure class and strand; classes may overlap.')
        add('Coverage',cat+'_product_features',g.product_id.nunique())
    hp=data['Hairpins'];pr=data['Promoters']
    add('Candidates','hairpins_nonoverlapping_per_context',int(hp.selected_nonredundant.sum()) if len(hp) else 0)
    add('Candidates','hairpins_sequence_prioritized',int((hp.priority=='sequence_prioritized').sum()),'Heuristic, not calibrated confidence; overlapping alternatives included.')
    add('Candidates','promoters_sequence_prioritized',int((pr.priority=='sequence_prioritized').sum()),'Heuristic, not calibrated confidence.')
    add('Neighborhoods','multicontext_neighborhoods',int((data['Neighborhoods'].n_contexts>1).sum()),'Candidate physical groups, NOT established operons.')
    add('Neighborhoods','singletons',int((data['Neighborhoods'].n_contexts==1).sum()),'Not called operons.')
    for test,g in data['Null_tests'].groupby('test'):
        add('Shuffled controls',test+'_contexts_tested',len(g))
        add('Shuffled controls',test+'_q_lt_0_05',int((g.BH_q_within_test_family<.05).sum()),'Approximate conditional anchor-level screening; not per-hit function or a uniformly sampled null.')
        add('Shuffled controls',test+'_minimum_q',g.BH_q_within_test_family.min())
    tr=data['Translation_QC']
    add('Quality control','CDS_translation_matches',int(tr.matches_annotation.sum()) if len(tr) else 0)
    add('Quality control','coordinate_round_trip_pass',bool(data['Coordinate_QC'].sequence_round_trip_pass.all()) if len(data['Coordinate_QC']) else True)
    add('Limitations','functional_validation','None','No experimental regulatory data supplied. Full inference limitations are in README and methods.')
    data['Summary']=pd.DataFrame(rows)
    return data


def save(fig, folder, name, caption, rows):
    fig.savefig(folder/f'{name}.jpg',dpi=600,pil_kwargs={'quality':96,'subsampling':0})
    # Matplotlib colorbars otherwise rasterize gradients automatically.
    for artist in fig.findobj():
        if hasattr(artist,'set_rasterized'):artist.set_rasterized(False)
    fig.savefig(folder/f'{name}.svg',format='svg')
    rows.append({'figure':name,'JPEG':f'figures/{name}.jpg','SVG':f'figures/{name}.svg',
                 'caption':caption,'SVG_text':'outlined paths; no font dependency',
                 'JPEG_dpi':600,'width_inches':fig.get_figwidth(),'height_inches':fig.get_figheight()})
    plt.close(fig)


def create_figures(data,record,folder,temperature):
    folder=Path(folder);rows=[];inv=data['Inventory_associations'];sites=data['Unique_sites']
    n=len(record);regions=data['Genome_regions'];a=data['Boundary_contexts']
    fig,axes=plt.subplots(3,1,figsize=(11,7),sharex=True,gridspec_kw={'height_ratios':[1,1.4,1.5]})
    regcolors=['#BCD5E8','#A7CEC3','#E8D6B9','#A7CEC3']
    for i,r in regions.iterrows():
        axes[0].broken_barh([(r.start_1b/1000,r.length_bp/1000)],(0,1),facecolors=regcolors[i%4])
        lab='IR' if r.feature_type=='repeat_region' else ('LSC' if 'large' in r.description else 'SSC')
        axes[0].text((r.start_1b+r.length_bp/2)/1000,.5,f'{lab}\n{r.length_bp:,} bp',ha='center',va='center',fontsize=8)
    axes[0].set(yticks=[],ylim=(0,1),title=f'{record.annotations.get("organism","")} chloroplast | {record.id} | {n:,} bp')
    for _,r in a.iterrows():
        y=.25 if r.strand==1 else -.65
        axes[1].broken_barh([((r.context_start_1b-1)/1000,(r.context_end_1b-r.context_start_1b+1)/1000)],(y,.35),facecolors='#3F6574' if r.feature_type=='CDS' else '#A6A793')
    axes[1].set(yticks=[-.48,.42],yticklabels=['− strand','+ strand'],ylim=(-.85,.85),title='Annotated product contexts (trans-spliced exon contexts shown separately)')
    cats=list(COLORS)
    for i,cat in enumerate(cats):
        sub=sites[sites.element_type==cat]
        x=[int(p.split(',')[0])/1000 for p in sub.position_path_1b]
        axes[2].scatter(x,[i]*len(x),marker='|',s=100,c=COLORS[cat],linewidths=1)
    axes[2].set(yticks=range(len(cats)),yticklabels=[LABELS[c] for c in cats],xlabel='Genomic coordinate (kb)',xlim=(0,n/1000),title='Category-specific unique candidate sites (not validated elements)')
    fig.tight_layout();save(fig,folder,'01_genome_candidate_map','Linearized circular genome; annotations are context spans. Site ticks show first base in transcription orientation; overlapping motif classes are not independent.',rows)

    fig,axes=plt.subplots(1,2,figsize=(10,4.7))
    counts=inv.element_type.value_counts();uniq=sites.element_type.value_counts()
    x=np.arange(len(cats))
    axes[0].barh(x-.18,[counts.get(c,0) for c in cats],height=.34,color=[COLORS[c] for c in cats],label='Context associations')
    axes[0].barh(x+.18,[uniq.get(c,0) for c in cats],height=.34,facecolor='white',edgecolor=[COLORS[c] for c in cats],linewidth=1.5,label='Unique sites per class')
    axes[0].set(yticks=x,yticklabels=[LABELS[c] for c in cats],xlabel='Count',title='Inventory size');axes[0].invert_yaxis();axes[0].legend(fontsize=7)
    covered=[]
    for c in cats:covered.append(inv[inv.element_type==c].product_id.nunique())
    axes[1].barh(x,covered,color=[COLORS[c] for c in cats])
    axes[1].set(yticks=x,yticklabels=[LABELS[c] for c in cats],xlabel='Annotated products with ≥1 candidate',title='Product coverage');axes[1].invert_yaxis()
    for i,v in enumerate(covered):axes[1].text(v+.3,i,str(v),va='center')
    fig.suptitle('Sequence-based candidates: overlapping classes and shared context windows',fontsize=12,fontweight='bold')
    fig.tight_layout();save(fig,folder,'02_inventory_counts','Hairpin counts include overlapping maximal alternatives; isolated TATAAT may occur inside octamers or paired promoters. Counts must not be summed as independent regulatory elements.',rows)

    pr=data['Promoters'];pairs=pr[pr.element_type=='PEP_paired_boxes'];sd=data['SD_like'];hp=data['Hairpins']
    fig,axes=plt.subplots(2,2,figsize=(10,7))
    axes[0,0].hist(-pairs.relative_end_to_boundary,bins=np.arange(0,321,20),color=COLORS['PEP_paired_boxes'],edgecolor='white')
    axes[0,0].set(xlabel='−10 box last base: bp upstream of annotated boundary',ylabel='Associations',title='Paired promoter location (not TSS distance)')
    axes[0,1].hist(pairs.spacer_bp,bins=np.arange(14.5,20.5,1),color=COLORS['PEP_paired_boxes'],edgecolor='white')
    axes[0,1].set(xlabel='Bases between −35 and −10 boxes',ylabel='Associations',title='Promoter-pair spacer')
    if len(sd):axes[1,0].hist(sd.spacer_to_start_codon_bp,bins=np.arange(2.5,16.5,1),color=COLORS['SD_like'],edgecolor='white')
    axes[1,0].set(xlabel='Bases after SD-like motif, before CDS start',ylabel='Associations',title='Corrected SD-like spacing')
    if len(hp):axes[1,1].hist(hp.distance_after_boundary_bp,bins=np.arange(.5,121.5,10),color=COLORS['downstream_hairpin'],edgecolor='white')
    axes[1,1].set(xlabel='Hairpin first base: bp after annotated boundary',ylabel='Associations',title='Downstream structure location')
    fig.tight_layout();save(fig,folder,'03_candidate_spacing','Distances reference annotated product/exon boundaries, not measured RNA ends or TSSs.',rows)

    fig,axes=plt.subplots(1,2,figsize=(10,4.5))
    if len(hp):
        scatter=axes[0].scatter(hp.stem_pairs,hp.intended_energy_kcal_mol,c=hp.tail_T_fraction,cmap='viridis',s=25,alpha=.65,edgecolors='none',vmin=0,vmax=1)
        fig.colorbar(scatter,ax=axes[0],label='T fraction in available next 8 bases')
        axes[0].axhline(-5,color='#555555',ls='--',lw=.8)
        axes[1].scatter(hp.intended_energy_kcal_mol,hp.MFE_kcal_mol,c=COLORS['downstream_hairpin'],s=20,alpha=.5)
        lo=min(hp.intended_energy_kcal_mol.min(),hp.MFE_kcal_mol.min())-1
        axes[1].plot([lo,1],[lo,1],ls='--',lw=.8,c='#555555')
    axes[0].set(xlabel='Perfect Watson–Crick stem pairs',ylabel='Intended-structure free energy (kcal/mol)',title=f'Downstream hairpins at {temperature:g} °C')
    axes[1].set(xlabel='Intended-structure free energy (kcal/mol)',ylabel='Unconstrained MFE (kcal/mol)',title='RNA folding supports structure, not function')
    fig.tight_layout();save(fig,folder,'04_hairpin_folding','ViennaRNA Turner 2004 energy model. Perfect stems 4–12 bp, loops 3–10 nt, intended energy ≤−2 kcal/mol. A hairpin can mediate RNA maturation/stability rather than termination.',rows)

    nt=data['Null_tests'];tests=list(nt.test.unique())
    fig,axes=plt.subplots(1,2,figsize=(10,4.6))
    for i,test in enumerate(tests):
        vals=nt[nt.test==test]
        axes[0].scatter(np.full(len(vals),i),-np.log10(vals.empirical_P),s=15,alpha=.4,c=list(COLORS.values())[i%5])
        axes[1].scatter(np.full(len(vals),i),-np.log10(vals.BH_q_within_test_family),s=15,alpha=.4,c=list(COLORS.values())[i%5])
    shortlabs=['PEP max\nscore','Octamer\ncount','TATAAT\ncount','SD max\nlength','Hairpin\nmin energy']
    for ax in axes:
        ax.set(xticks=range(len(tests)),xticklabels=shortlabs[:len(tests)])
        ax.axhline(-np.log10(.05),ls='--',lw=.8,c='#555555');ax.set_ylabel('−log10(value)')
    axes[0].set_title('Conditional empirical P-values');axes[1].set_title('BH q-values within test family')
    fig.suptitle('Dinucleotide-preserving shuffled-window controls (per context, not per hit)',fontweight='bold',fontsize=11)
    fig.tight_layout();save(fig,folder,'05_shuffle_controls','Randomized Euler shuffles preserve dinucleotide counts and ends but are not uniformly sampled. P=(1+extreme)/(B+1); BH correction within each test family. Shared windows, repeated IR copies and overlapping motif classes are dependent.',rows)

    g=data['Gene_summary'];columns=['PEP_pairs','algal_octamers','SD_like_sites','nonredundant_hairpins']
    # All products are included; two pages avoid unreadable 110-row plots.
    for page,start in enumerate(range(0,len(g),55),1):
        sub=g.iloc[start:start+55]
        fig,ax=plt.subplots(figsize=(8,max(6,len(sub)*.19+1.7)))
        matrix=sub[columns].to_numpy(dtype=float)
        # A vector QuadMesh avoids raster <image> elements in SVG/CorelDRAW.
        image=ax.pcolormesh(np.arange(5)-.5,np.arange(len(sub)+1)-.5,
                            np.log1p(matrix),cmap='Blues',shading='flat',rasterized=False)
        ax.set_ylim(len(sub)-.5,-.5)
        ax.set(yticks=np.arange(len(sub)),yticklabels=[f'{r.gene} [{r.product_id}]' for _,r in sub.iterrows()],
               xticks=np.arange(4),xticklabels=['PEP pairs','TATAATAT','SD-like','3′ hairpins\n(nonoverlapping)'],
               title=f'All annotated products: candidate associations (page {page})')
        ax.tick_params(axis='y',labelsize=6.8)
        for y in range(len(sub)):
            for x in range(4):ax.text(x,y,str(int(matrix[y,x])),ha='center',va='center',fontsize=6.5,color='white' if np.log1p(matrix[y,x])>image.norm.vmax*.6 else '#1D3445')
        fig.colorbar(image,ax=ax,fraction=.03,pad=.04,label='ln(1 + association count)')
        fig.tight_layout();save(fig,folder,f'06_product_coverage_page_{page}','Every annotated product copy is included and identified separately. Mixed-strand trans-spliced contexts are aggregated by product; duplicated IR genes are not collapsed.',rows)

    # Draw selected structures with actual ViennaRNA layout, not a generic cartoon.
    if len(hp):
        import RNA
        best=hp[hp.selected_nonredundant].sort_values('intended_energy_kcal_mol').drop_duplicates('product_id').head(6)
        fig,axs=plt.subplots(2,3,figsize=(11,7.2))
        for ax,(_,r) in zip(axs.flat,best.iterrows()):
            sequence=r.sequence_5to3.replace('T','U');structure=r.MFE_structure
            coordinates=RNA.get_xy_coordinates(structure)
            xy=np.array([[coordinates.get(i).X,coordinates.get(i).Y] for i in range(len(sequence))])
            ax.plot(xy[:,0],xy[:,1],c='#B7BDC0',lw=1,zorder=1)
            pt=RNA.ptable(structure)
            for i in range(len(sequence)):
                if pt[i+1]>i+1:
                    j=pt[i+1]-1;ax.plot(xy[[i,j],0],xy[[i,j],1],c='#8D609A',lw=1.1,zorder=1)
            for i,base in enumerate(sequence):
                ax.text(*xy[i],base,ha='center',va='center',fontsize=6.2,color='#173E50',bbox={'boxstyle':'circle,pad=.12','fc':'white','ec':'none'},zorder=2)
            ax.set_title(f'{r.gene} [{r.product_id}]\nMFE {r.MFE_kcal_mol:.1f} kcal/mol; +{r.distance_after_boundary_bp} bp',fontsize=9)
            ax.set_aspect('equal');ax.axis('off');ax.margins(.15)
        for ax in axs.flat[len(best):]:ax.axis('off')
        fig.suptitle('Examples of strongest predicted downstream RNA structures',fontweight='bold')
        fig.tight_layout();save(fig,folder,'07_top_RNA_structures','Unconstrained MFE structures of six strongest nonoverlapping candidates from distinct product features; selected by intended-stem energy. These are not validated terminators.',rows)
    data['Figure_manifest']=pd.DataFrame(rows)
    data['Figure_manifest'].to_csv(folder.parent/'results'/'Figure_manifest.tsv',sep='\t',index=False)


def create_workbook(data,root,path):
    readme=[
        ('Purpose','Exploratory sequence inventory, NOT a comprehensive experimentally validated regulatory atlas.'),
        ('Input','Annotated NCBI GenBank genome QB176266.1; raw record and SHA256 included. No user sequencing or transcriptomics data were supplied.'),
        ('Main results','Gene_summary covers every annotated product feature, including zero-hit products. Inventory_associations records candidate/context associations. Unique_sites deduplicates identical class/strand/position paths.'),
        ('Nonindependent classes','An isolated TATAAT can be inside an octamer or paired-box promoter. Multiple genes can share a hit. Overlapping alternative hairpins are retained; selected_nonredundant supplies a nonoverlapping per-context subset. Class-specific unique sites are not independent elements.'),
        ('Coordinates','All exported coordinates are 1-based inclusive except BED (0-based half-open). position_path_1b is exact genomic positions in 5prime-to-3prime orientation. For wraps, min/max coordinates are bounding values only; use segments/path.'),
        ('Boundaries are not RNA ends','Upstream is 300 bp before product/exon annotation; downstream is 120 bp after it. No TSS/UTR boundaries were measured. Full windows are searched even across other coding features; overlap_CDS_bp explicitly flags them.'),
        ('Trans-splicing','Mixed-strand psaA and rpl32 are separate exon contexts. SD search only at first annotated CDS exon. Other cis-spliced features use biological first/last parts. GenBank part ordering is assumed, not experimentally proven.'),
        ('Promoters','Paired PEP-like boxes: TTGACA and TATAAT, <=1 mismatch each, 15–19 intervening bases. Also exact algal TATAATAT and exact isolated TATAAT; -35 motifs alone are diagnostic only. Thresholds are heuristic, not trained/calibrated on this species.'),
        ('NEP correction','Original five-base NEP regex was not YRTA. Correct YRTA=[CT][AG]TA and GAA words are retained only in Raw_motif_diagnostics; neither establishes a NEP promoter, particularly in a chlorophyte.'),
        ('SD-like','Exact AGGAGG/GGAGG/AGGAG/GAGG/AGGA/GGAG, 3–15 bases AFTER motif and BEFORE start codon. Contained shorter motifs are removed. Absence does not imply failed translation in algae.'),
        ('Hairpins','Antiparallel reverse-complement stems 4–12 bp, loops 3–10 nt, no GU/bulges; intended-structure energy <=−2 kcal/mol. ViennaRNA unconstrained MFE, intended energy, pair retention and T-rich tail context provided. Structures may affect RNA processing/stability rather than termination.'),
        ('Temperature','25 C is an explicit modeling assumption, not a measured growth temperature. Turner2004 RNA parameters; do not compare energies from different model settings as if identical.'),
        ('Priorities','sequence_prioritized is a heuristic sequence/context label, NOT high confidence or biological validation. PEP: <=1 total mismatch, 16–18 spacer, no CDS overlap. Hairpin: >=6 pairs, <=−5 kcal/mol, >=75% intended pairs retained, no CDS overlap.'),
        ('Null model','1,000 randomized Euler dinucleotide-preserving upstream/SD shuffles and 200 downstream shuffles per context. Ends fixed; sampling not uniform. Conditional diagnostic P-values, not validated promoter probabilities.'),
        ('Statistics','P=(1+extreme)/(B+1). PEP and SD use best score, motifs use count, hairpins use strongest intended energy. BH adjustment separately within each test family. Values merged to hits are ANCHOR-level tests, not per-hit tests. q=.05 does not establish function.'),
        ('Operon limitation','Neighborhoods: physically adjacent same-strand contexts with -50..200 bp gap; opposite-strand neighbors block links; circular last/first handled. Singleton contexts are not operons. No cotranscription evidence supplied.'),
        ('Not assayed','No validated TSSs, promoters, true terminators, PPR/protein-binding sites, riboswitches, RNA thermometers, RNA editing, regulatory sRNAs, or nuclear-encoded regulators are inferred by this workflow. Imperfect stems/GU/bulges, distant promoters >300 bp and context >120 bp remain blind spots.'),
        ('Validation needed','Strand-specific/differential RNA-seq or 5prime-RACE for TSSs, 3prime-RACE/Term-seq and run-on assays for RNA ends vs termination, RT-PCR/long-read RNA for operons/trans-splicing, reporter mutagenesis and ribosome profiling for function.'),
        ('CorelDRAW','SVG 1.1 vector exports have outlined text (font independent) and no embedded raster images. Structurally checked, but not opened in proprietary CorelDRAW here. JPEGs are 600 dpi with stored resolution metadata.'),
        ('Reproducibility','Improved analysis/reporting/test scripts, frozen dependencies, inputs, methods, audit, manifests, TSVs, BED/GFF and FASTA included in ZIP. No software was installed on your own computer; installation was in this workspace.')]
    read=pd.DataFrame(readme,columns=['topic','details'])
    extra={}
    for name,filename in [('Code_audit','code_audit.tsv'),('References','references.tsv'),('Unit_tests','unit_test_results.tsv'),('Deliverable_QC','deliverable_qc.tsv')]:
        f=Path(root)/'logs'/filename
        if f.exists():extra[name]=pd.read_csv(f,sep='\t').fillna('')
    # Self-describing data dictionary for every exported tabular column.
    known={
        'association_id':'Unique candidate-to-context association ID; not a unique genomic site.',
        'product_id':'Unique GenBank product-feature ID; duplicated gene copies remain separate.',
        'anchor_id':'Unique oriented product/exon boundary context; not a measured transcription unit.',
        'genomic_segments_1b':'Semicolon-separated exact 1-based inclusive runs in transcription order.',
        'position_path_1b':'Comma-separated genomic positions, 1-based, in 5prime-to-3prime sequence order.',
        'anchor_empirical_P':'Conditional shuffled-window P-value for the entire anchor statistic; NOT per-hit.',
        'anchor_BH_q':'BH adjusted P within the corresponding anchor-level test family.',
        'BH_q_within_test_family':'BH adjusted P across tested anchors within one statistic family.',
        'sequence_5to3':'DNA alphabet sequence in the candidate transcript orientation.',
        'MFE_kcal_mol':'Unconstrained RNA minimum free energy at the configured temperature; kcal/mol.',
        'intended_energy_kcal_mol':'Energy of the enumerated exact stem-loop structure; kcal/mol.',
        'selected_nonredundant':'Greedy strongest-energy nonoverlapping subset within the same anchor.',
        'overlap_CDS_bp':'Number of hit bases overlapping any annotated CDS exon, either strand.',
        'overlap_other_products':'Other annotated product IDs with overlapping exon bases.',
        'priority':'Heuristic sequence prioritization only; not a calibrated confidence score.',
        'genomic_start_1b':'Lowest occupied genomic coordinate, 1-based inclusive; wrap hits require segments/path.',
        'genomic_end_1b':'Highest occupied genomic coordinate, 1-based inclusive; wrap hits require segments/path.'}
    dic=[]
    for sheet,df in data.items():
        for col in df.columns:
            definition=known.get(col,'')
            if not definition:
                if col.endswith('_bp'):definition='Length or distance in base pairs; see sheet and methods for reference boundary.'
                elif 'strand'==col:definition='+1 forward, -1 reverse; empty for mixed/unspecified.'
                else:definition=col.replace('_',' ')+'; interpreted in the context of this worksheet.'
            dic.append({'sheet':sheet,'column':col,'dtype':str(df[col].dtype),'definition':definition})
    frames={'README':read}
    if 'Summary' in data:frames['Summary']=data['Summary']
    frames.update(data);frames.update(extra);frames['Data_dictionary']=pd.DataFrame(dic)
    with pd.ExcelWriter(path,engine='xlsxwriter',engine_kwargs={'options':{'strings_to_formulas':False,'strings_to_urls':False}}) as writer:
        book=writer.book
        head=book.add_format({'bold':True,'font_color':'white','bg_color':'#244B60','text_wrap':True,'valign':'top'})
        text=book.add_format({'valign':'top'});wrap=book.add_format({'text_wrap':True,'valign':'top'})
        number=book.add_format({'num_format':'0.0000','valign':'top'})
        for name,df in frames.items():
            sheet=name[:31]
            if len(df)>1048574:raise ValueError(f'Excel row limit exceeded for {sheet}; split sheets required')
            if len(df.columns)==0:df=pd.DataFrame({'status':['No rows in this category']})
            df=df.replace([np.inf,-np.inf],np.nan)
            df.to_excel(writer,sheet_name=sheet,index=False)
            ws=writer.sheets[sheet];ws.freeze_panes(1,0);ws.set_zoom(90)
            ws.set_row(0,32,head)
            for j,c in enumerate(df.columns):
                sample=df[c].head(200).astype(str).map(len)
                width=min(46,max(12,len(c)+2,float(sample.quantile(.8))+2 if len(sample) else 12))
                ws.set_column(j,j,width,number if pd.api.types.is_float_dtype(df[c]) else text)
                ws.write(0,j,c,head)
            if len(df):ws.autofilter(0,0,len(df),len(df.columns)-1)
            if name=='README':
                ws.set_column(0,0,25,wrap);ws.set_column(1,1,115,wrap)
                for i,(_,r) in enumerate(df.iterrows(),1):ws.set_row(i,max(34,14*int(np.ceil(len(r.details)/105))))
            if name=='Null_tests' and len(df):
                col=df.columns.get_loc('BH_q_within_test_family')
                ws.conditional_format(1,col,len(df),col,{'type':'cell','criteria':'<','value':.05,'format':book.add_format({'bg_color':'#D8ECDC'})})
            if name=='Gene_summary':
                for c in ['PEP_pairs','algal_octamers','SD_like_sites','nonredundant_hairpins']:
                    j=df.columns.get_loc(c);ws.conditional_format(1,j,len(df),j,{'type':'data_bar','bar_color':'#8CBDB5'})
        overview=book.add_worksheet('Figure_previews');overview.hide_gridlines(2)
        overview.set_column(0,18,11)
        row=0
        for _,r in data.get('Figure_manifest',pd.DataFrame()).iterrows():
            overview.write(row,0,r.figure,head)
            overview.insert_image(row+1,0,str(Path(root)/r.JPEG),{'x_scale':.12,'y_scale':.12})
            row+=int(r.height_inches*600*.12/20)+5
