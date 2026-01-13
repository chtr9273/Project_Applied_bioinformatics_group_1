#!/usr/bin/env python3
"""
Plot overlap between two BAM files
Usage: python plot_overlap.py genome1_name genome2_name unique1 unique2 overlap count1 count2 outdir
"""

import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Read arguments
genome1 = sys.argv[1]
genome2 = sys.argv[2]
unique1 = int(sys.argv[3])
unique2 = int(sys.argv[4])
overlap = int(sys.argv[5])
count1 = int(sys.argv[6])
count2 = int(sys.argv[7])
outdir = sys.argv[8]

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# ===== Venn Diagram =====
v = venn2(subsets=(unique1, unique2, overlap),
          set_labels=(genome1, genome2),
          ax=ax1)

# Colors
v.get_patch_by_id('10').set_color('#3498db')
v.get_patch_by_id('01').set_color('#e74c3c')
v.get_patch_by_id('11').set_color('#9b59b6')

# Labels with percentages 
pct1 = (unique1 / count1 * 100) if count1 > 0 else 0
pct2 = (unique2 / count2 * 100) if count2 > 0 else 0
pct_overlap1 = (overlap / count1 * 100) if count1 > 0 else 0
pct_overlap2 = (overlap / count2 * 100) if count2 > 0 else 0

if v.get_label_by_id('10'):
    v.get_label_by_id('10').set_text(f'{unique1:,}\n({pct1:.1f}%)')
    v.get_label_by_id('10').set_fontsize(14)  

if v.get_label_by_id('01'):
    v.get_label_by_id('01').set_text(f'{unique2:,}\n({pct2:.1f}%)')
    v.get_label_by_id('01').set_fontsize(14)  

if v.get_label_by_id('11'):
    v.get_label_by_id('11').set_text(f'{overlap:,}\n({pct_overlap1:.1f}% / {pct_overlap2:.1f}%)')
    v.get_label_by_id('11').set_fontsize(14)  

# Set labels
for text in v.set_labels:
    text.set_fontsize(20)  # Species names

ax1.set_title('Read Overlap Between Genomes', fontsize=20, fontweight='bold', pad=20)  
ax1.text(0.5, -0.15, f'Total reads: {genome1}={count1:,} | {genome2}={count2:,}',
         ha='center', transform=ax1.transAxes, fontsize=14, style='italic')  

# ===== Bar Chart =====
categories = [f'{genome1}\nUnique', 'Overlapping', f'{genome2}\nUnique']
values = [unique1, overlap, unique2]
colors = ['#3498db', '#9b59b6', '#e74c3c']

bars = ax2.bar(categories, values, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)

# Add value labels on bars 
for bar in bars:
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
            f'{int(height):,}',
            ha='center', va='bottom', fontsize=16, fontweight='bold')  

ax2.set_ylabel('Number of Reads', fontsize=18, fontweight='bold')  
ax2.set_title('Read Distribution', fontsize=20, fontweight='bold', pad=20) 
ax2.grid(axis='y', alpha=0.3, linestyle='--')
ax2.set_ylim(0, max(values) * 1.15)

# X-axis labels 
ax2.tick_params(axis='x', labelsize=13)  # Category labels 
ax2.tick_params(axis='y', labelsize=14)  # Y-axis numbers

plt.suptitle(f'BAM Read Comparison: {genome1} vs {genome2}',
             fontsize=22, fontweight='bold', y=0.98)  

plt.tight_layout(rect=[0, 0.03, 1, 0.96])
plt.savefig(f'{outdir}/overlap_plot.png', dpi=300, bbox_inches='tight')
print(f"\n  Saved: {outdir}/overlap_plot.png")