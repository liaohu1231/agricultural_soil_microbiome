#!/usr/bin/env python3
"""
Easyfig‑like plot for contigs containing target genes.
- Pure Matplotlib implementation, no dna_features_viewer dependency.
- Unannotated genes: arrows only, no labels.
- Contig names on left, similarity bands between adjacent contigs.
- For long contigs (>50kb) zoom to target gene ± flank.
- Output PDF with editable text.
- FIXED: single-polygon arrow, no seam, proper aspect ratio.
"""

import argparse
import os
import subprocess
import tempfile
import random
import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from Bio import SeqIO

# ------------------- 1. Parse Prodigal .faa header -------------------
def parse_faa_header(faa_file):
    genes = []
    for record in SeqIO.parse(faa_file, "fasta"):
        desc = record.description
        if desc.startswith('>'):
            desc = desc[1:]
        parts = [p.strip() for p in desc.split('#')]
        if len(parts) < 4:
            continue
        gene_id_full = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        strand_code = int(parts[3])
        strand = '+' if strand_code == 1 else '-'
        if '_' in gene_id_full:
            contig = '_'.join(gene_id_full.split('_')[:-1])
        else:
            contig = gene_id_full
        genes.append({
            'gene_id': record.id,
            'contig': contig,
            'start': start,
            'end': end,
            'strand': strand,
            'seq': str(record.seq)
        })
    return pd.DataFrame(genes)

# ------------------- 2. Load KEGG annotation -------------------
def load_kegg_annotation(kegg_file):
    try:
        df = pd.read_csv(kegg_file, sep='\t', header=0)
    except:
        df = pd.read_csv(kegg_file, sep=',', header=0)

    gene_col = None
    for col in ['gene name', '#query', 'Gene ID', 'gene_id', 'query']:
        if col in df.columns:
            gene_col = col
            break
    if gene_col is None:
        gene_col = df.columns[0]

    ko_col = None
    for col in ['K Number', 'ko', 'KEGG_ko', 'KO', '#ko', 'k_number']:
        if col in df.columns:
            ko_col = col
            break
    if ko_col is None:
        for col in df.columns:
            if 'ko' in col.lower():
                ko_col = col
                break
    gene_to_ko = {}
    if ko_col:
        gene_to_ko = dict(zip(df[gene_col], df[ko_col]))

    name_col = None
    for col in ['Gene name', 'Description', 'Product', 'definition']:
        if col in df.columns:
            name_col = col
            break
    if name_col is None:
        name_col = gene_col
    gene_to_name = dict(zip(df[gene_col], df[name_col]))

    return gene_to_ko, gene_to_name

# ------------------- 3. Determine contig order and target genes -------------------
def get_contig_order_and_target_genes(gene_list_file, gene_to_name, pos_df):
    with open(gene_list_file) as f:
        target_names = [line.strip() for line in f if line.strip()]
    print(f"Target genes: {target_names}")

    gene_to_name_lower = {gid: str(name).lower() for gid, name in gene_to_name.items()}
    target_lower = [n.lower() for n in target_names]

    contig_order = []
    seen_contigs = set()
    contig_target_genes = {}

    for t in target_lower:
        for gid, name_lower in gene_to_name_lower.items():
            if t in name_lower:
                contig = pos_df.loc[pos_df['gene_id'] == gid, 'contig'].values
                if len(contig) == 0:
                    continue
                c = contig[0]
                contig_target_genes.setdefault(c, []).append(gid)
                if c not in seen_contigs:
                    contig_order.append(c)
                    seen_contigs.add(c)
                break
    print(f"Contig order: {contig_order}")
    return contig_order, contig_target_genes

# ------------------- 4. Extract contig sequences -------------------
def extract_contig_sequences(fasta_file, contig_names):
    records = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}
    return {name: records[name] for name in contig_names if name in records}

# ------------------- 5. Run minimap2 -------------------
def run_minimap2_all_vs_all(seqs_dict, out_paf, identity_threshold=70):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        for name, seq in seqs_dict.items():
            f.write(f'>{name}\n{seq}\n')
        temp_fasta = f.name
    cmd = f"minimap2 -x asm5 --paf-no-hit -c {temp_fasta} {temp_fasta} > {out_paf}"
    subprocess.run(cmd, shell=True, check=True)
    os.unlink(temp_fasta)

    hits = []
    with open(out_paf) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 12:
                continue
            qseq = parts[0]
            qstart = int(parts[2])
            qend = int(parts[3])
            sseq = parts[5]
            sstart = int(parts[7])
            send = int(parts[8])
            identity = 100.0
            for field in parts[12:]:
                if field.startswith('id:f:'):
                    identity = float(field.split(':')[2]) * 100
                    break
            if identity >= identity_threshold:
                hits.append({
                    'qseq': qseq, 'qstart': qstart, 'qend': qend,
                    'sseq': sseq, 'sstart': sstart, 'send': send,
                    'identity': identity
                })
    return hits

# ------------------- 6. Single-polygon arrow drawing (no seam) -------------------
def draw_gene_arrow(ax, start, end, y, strand, color, label=None, label_size=8):
    """
    Draw a single gene arrow as one polygon (rectangle + triangle seamlessly).
    Arrowhead length is dynamically adjusted.
    """
    width = end - start
    if width <= 0:
        return
    height = 0.6
    y_bottom = y - height/2
    y_top = y + height/2

    # Arrowhead length: between 10 and 150 bp, or 15% of width (whichever smaller)
    arrow_len = max(10, min(150, width * 0.15))

    if strand == '+':
        # Right-pointing arrow
        rect_width = width - arrow_len
        if rect_width < 0:
            rect_width = 0
        # Vertices: bottom-left, bottom-right (end of rectangle), arrow tip, top-right, top-left
        vertices = [
            (start, y_bottom),                               # bottom-left of rectangle
            (start + rect_width, y_bottom),                  # bottom-right of rectangle
            (start + rect_width + arrow_len, y),             # arrow tip
            (start + rect_width, y_top),                     # top-right of rectangle
            (start, y_top)                                   # top-left of rectangle
        ]
        poly = Polygon(vertices, closed=True, facecolor=color, edgecolor='black', linewidth=0.5)
        ax.add_patch(poly)
        # Label positioned at the center of the rectangle part (if any)
        if label and rect_width > 0:
            ax.text(start + rect_width/2, y, label, ha='center', va='center', fontsize=label_size)
        elif label and rect_width <= 0:
            ax.text(start + width/2, y, label, ha='center', va='center', fontsize=label_size)
    else:
        # Left-pointing arrow
        rect_width = width - arrow_len
        if rect_width < 0:
            rect_width = 0
        vertices = [
            (start + arrow_len + rect_width, y_bottom),      # bottom-right of rectangle
            (start + arrow_len, y_bottom),                   # bottom-left of rectangle (arrow base)
            (start, y),                                      # arrow tip
            (start + arrow_len, y_top),                      # top-left of rectangle
            (start + arrow_len + rect_width, y_top)          # top-right of rectangle
        ]
        poly = Polygon(vertices, closed=True, facecolor=color, edgecolor='black', linewidth=0.5)
        ax.add_patch(poly)
        if label and rect_width > 0:
            ax.text(start + arrow_len + rect_width/2, y, label, ha='center', va='center', fontsize=label_size)
        elif label and rect_width <= 0:
            ax.text(start + width/2, y, label, ha='center', va='center', fontsize=label_size)

def plot_single_contig(contig, group_df, ko_color_map, ax, display_range=None):
    group_df = group_df.sort_values('start')
    contig_len = group_df['end'].max()
    if display_range is not None:
        start_lim, end_lim = display_range
        group_df = group_df[(group_df['start'] >= start_lim) & (group_df['end'] <= end_lim)].copy()
        if group_df.empty:
            ax.set_xlim(start_lim, end_lim)
            ax.set_ylim(0, 1)
            return contig_len
    else:
        start_lim, end_lim = 1, contig_len

    # Add small margin to avoid clipping arrowheads
    margin = 100
    ax.set_xlim(start_lim - margin, end_lim + margin)

    for _, row in group_df.iterrows():
        ko = row['ko']
        color = ko_color_map.get(ko, '#D3D3D3')
        label = ko if ko != 'Unknown' else None
        strand = row['strand']
        draw_gene_arrow(ax, row['start'], row['end'], y=0.5, strand=strand,
                        color=color, label=label, label_size=8)

    ax.set_ylim(0, 1)
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(axis='y', left=False, labelleft=False)
    return contig_len

# ------------------- 7. Add curved similarity band -------------------
def add_curved_bands(ax1, ax2, hit, fig, cmap=plt.cm.RdYlGn, linewidth=2, alpha=0.7):
    x1_mid = (hit['qstart'] + hit['qend']) / 2
    x2_mid = (hit['sstart'] + hit['send']) / 2
    y1 = 0.5
    y2 = 0.5
    trans1 = ax1.transData
    trans2 = ax2.transData
    p1 = trans1.transform((x1_mid, y1))
    p2 = trans2.transform((x2_mid, y2))
    color = cmap(hit['identity'] / 100)
    con = ConnectionPatch(
        xyA=p1, xyB=p2, coordsA='figure points', coordsB='figure points',
        arrowstyle='-', color=color, linewidth=linewidth, alpha=alpha,
        connectionstyle="arc3,rad=0.3"
    )
    fig.add_artist(con)

# ------------------- Main -------------------
def main():
    parser = argparse.ArgumentParser(description="Easyfig‑like plot for contigs with target genes (pure matplotlib)")
    parser.add_argument("--fasta", required=True, help="Contig nucleotide sequences (FASTA)")
    parser.add_argument("--faa", required=True, help="Prodigal .faa file")
    parser.add_argument("--kegg", required=True, help="KEGG annotation table (TSV)")
    parser.add_argument("--gene_list", required=True, help="Text file with target gene names (one per line)")
    parser.add_argument("--outdir", default="./easyfig_output", help="Output directory")
    parser.add_argument("--identity_threshold", type=float, default=70, help="Minimum identity %% for bands")
    #parser.add_argument("--dpi", type=int, default=300, help="Output DPI")
    parser.add_argument("--flank", type=int, default=20000, help="Flanking region (bp) for long contigs (>50kb)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Parse data
    pos_df = parse_faa_header(args.faa)
    gene_to_ko, gene_to_name = load_kegg_annotation(args.kegg)
    pos_df['ko'] = pos_df['gene_id'].map(gene_to_ko).fillna('Unknown')
    pos_df['gene_name'] = pos_df['gene_id'].map(gene_to_name).fillna(pos_df['gene_id'])

    # 2. Get contig order
    contig_order, contig_target_genes = get_contig_order_and_target_genes(
        args.gene_list, gene_to_name, pos_df
    )
    if not contig_order:
        print("No target genes found. Exiting.")
        return

    # 3. Extract sequences
    seqs = extract_contig_sequences(args.fasta, contig_order)
    missing = set(contig_order) - set(seqs.keys())
    if missing:
        print(f"Warning: contigs missing in FASTA: {missing}")

    # 4. Run minimap2
    paf_file = os.path.join(args.outdir, "alignments.paf")
    hits = run_minimap2_all_vs_all(seqs, paf_file, args.identity_threshold)
    print(f"Found {len(hits)} alignment hits (identity >= {args.identity_threshold}%)")

    # 5. Prepare KO colour map
    all_kos = set(pos_df[pos_df['contig'].isin(contig_order)]['ko']) - {'Unknown'}
    random.seed(42)
    ko_color_map = {ko: "#{:06x}".format(random.randint(0, 0xFFFFFF)) for ko in all_kos}
    ko_color_map['Unknown'] = '#D3D3D3'

    # 6. Create subplots
    n = len(contig_order)
    fig_height = max(1.5, 0.8 * n)   # at least 3 inches
    fig, axes = plt.subplots(n, 1, figsize=(15, fig_height))
    if n == 1:
        axes = [axes]
    plt.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.05, hspace=1)
#    plt.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.05)

    contig_axes = {}
    for i, contig in enumerate(contig_order):
        group = pos_df[pos_df['contig'] == contig].copy()
        contig_len = group['end'].max()
        display_range = None
        if contig_len > 15000:
            target_ids = contig_target_genes.get(contig, [])
            if target_ids:
                target_df = group[group['gene_id'].isin(target_ids)]
                if not target_df.empty:
                    min_start = target_df['start'].min()
                    max_end = target_df['end'].max()
                    display_start = max(1, min_start - args.flank)
                    display_end = min(contig_len, max_end + args.flank)
                    display_range = (display_start, display_end)
                    print(f"Contig {contig} length {contig_len} bp > 50kb, showing [{display_start}, {display_end}]")
        ax = axes[i]
        plot_single_contig(contig, group, ko_color_map, ax, display_range)
        # Contig name on left
        ax.text(-0.05, 0.5, contig, transform=ax.transAxes, ha='right', va='center', fontsize=10)
        contig_axes[contig] = ax

    # 7. Add similarity bands (adjacent only)
    fig.canvas.draw()
    contig_index = {c: idx for idx, c in enumerate(contig_order)}
    for hit in hits:
        q = hit['qseq']
        s = hit['sseq']
        if q == s:
            continue
        if q not in contig_axes or s not in contig_axes:
            continue
        if abs(contig_index[q] - contig_index[s]) != 1:
            continue
        ax_q = contig_axes[q]
        ax_s = contig_axes[s]
        add_curved_bands(ax_q, ax_s, hit, fig, cmap=plt.cm.RdYlGn, linewidth=2, alpha=0.7)

    # 8. Save PDF
    out_pdf = os.path.join(args.outdir, "easyfig_comparison.pdf")
    plt.savefig(out_pdf, bbox_inches="tight")
    print(f"Saved figure to {out_pdf}")

if __name__ == "__main__":
    main()
