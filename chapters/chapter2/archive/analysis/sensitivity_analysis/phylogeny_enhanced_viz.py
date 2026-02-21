#!/usr/bin/env python3
"""
Phylogeny-Enhanced HOG Visualization
Creates side-by-side UMAP + phylogenetic tree figures for selected protein families

Focuses on 2 important human proteins:
1. P25685 (DNAJB1) - Heat Shock Protein 40 family (HOG 801468)
2. P11717 (IGF2R) - Insulin-like Growth Factor 2 Receptor family (HOG 792940)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, ConnectionPatch
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from collections import defaultdict
import re
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE_PATH = Path(__file__).parent.parent.parent
DATA_PATH = BASE_PATH / "datasets" / "cafa3_merged"
ANALYSIS_PATH = BASE_PATH / "analysis"
OUTPUT_PATH = Path(__file__).parent

# Taxonomy ID to species name mapping (common model organisms)
TAXONOMY_NAMES = {
    9606: "Human",
    10090: "Mouse", 
    10116: "Rat",
    7955: "Zebrafish",
    7227: "Drosophila",
    6239: "C. elegans",
    559292: "S. cerevisiae",
    284812: "S. pombe",
    3702: "A. thaliana",
    44689: "D. discoideum",
    7668: "S. purpuratus",
    36329: "P. falciparum",
    237561: "C. albicans",
    5664: "L. major",
}

# Color palette for different hierarchical depths
DEPTH_COLORS = [
    '#1f77b4',  # Root - blue
    '#ff7f0e',  # Level 1 - orange
    '#2ca02c',  # Level 2 - green
    '#d62728',  # Level 3 - red
    '#9467bd',  # Level 4 - purple
    '#8c564b',  # Level 5 - brown
    '#e377c2',  # Level 6 - pink
    '#7f7f7f',  # Level 7 - gray
    '#bcbd22',  # Level 8+ - olive
]

# Species colors for UMAP
SPECIES_COLORS = {
    9606: '#e41a1c',    # Human - red
    10090: '#377eb8',   # Mouse - blue
    10116: '#4daf4a',   # Rat - green  
    7955: '#984ea3',    # Zebrafish - purple
    7227: '#ff7f00',    # Drosophila - orange
    6239: '#ffff33',    # C. elegans - yellow
    559292: '#a65628',  # S. cerevisiae - brown
    284812: '#f781bf',  # S. pombe - pink
    3702: '#999999',    # A. thaliana - gray
    44689: '#66c2a5',   # D. discoideum - teal
}


def load_data():
    """Load all required data files"""
    print("Loading data...")
    
    # Load HOG cache
    hog_cache = pd.read_csv(DATA_PATH / "hog_cache.csv")
    print(f"  HOG cache: {len(hog_cache)} entries")
    
    # Load UMAP coordinates for HOG analysis
    umap_coords = pd.read_csv(ANALYSIS_PATH / "hog_analysis" / "hog_umap_coordinates.csv")
    print(f"  UMAP coordinates: {len(umap_coords)} points")
    
    return hog_cache, umap_coords


def parse_hog_hierarchy(hog_id):
    """
    Parse HOG ID into hierarchy levels
    Example: HOG:E0801468.10ivs.5885b -> [801468, '10ivs', '5885b']
    """
    if pd.isna(hog_id) or hog_id == '':
        return []
    
    # Remove HOG:E0 prefix
    clean_id = re.sub(r'^HOG:E0?', '', str(hog_id))
    
    # Split by dots
    parts = clean_id.split('.')
    
    # First part is root HOG (numeric)
    hierarchy = []
    if parts:
        try:
            hierarchy.append(int(parts[0]))
        except ValueError:
            hierarchy.append(parts[0])
        
        # Rest are child HOG suffixes
        hierarchy.extend(parts[1:])
    
    return hierarchy


def build_hog_tree(hog_cache, roothog_id):
    """
    Build a tree structure for a given root HOG
    Returns a nested dict representing the hierarchy
    """
    subset = hog_cache[hog_cache['roothog_id'] == roothog_id].copy()
    
    tree = {
        'id': roothog_id,
        'name': f'HOG:{roothog_id}',
        'proteins': [],
        'children': {},
        'depth': 0
    }
    
    for _, row in subset.iterrows():
        hierarchy = parse_hog_hierarchy(row['hog_id'])
        if not hierarchy:
            continue
            
        # Navigate/create path in tree
        current = tree
        for i, level in enumerate(hierarchy):
            if i == 0:
                # Root level - add protein if this is the root HOG
                if len(hierarchy) == 1:
                    current['proteins'].append({
                        'entry_id': row['EntryID'],
                        'oma_id': row.get('oma_id', ''),
                        'species': str(row.get('oma_id', ''))[:5] if pd.notna(row.get('oma_id')) else 'UNKN'
                    })
            else:
                # Child level
                level_str = str(level)
                if level_str not in current['children']:
                    current['children'][level_str] = {
                        'id': level_str,
                        'name': f'.{level_str}',
                        'proteins': [],
                        'children': {},
                        'depth': i
                    }
                current = current['children'][level_str]
                
                # If this is the last level, add the protein here
                if i == len(hierarchy) - 1:
                    current['proteins'].append({
                        'entry_id': row['EntryID'],
                        'oma_id': row.get('oma_id', ''),
                        'species': str(row.get('oma_id', ''))[:5] if pd.notna(row.get('oma_id')) else 'UNKN'
                    })
    
    return tree


def get_major_child_hogs(tree, max_children=8, min_proteins=3):
    """
    Get the most populated child HOGs for visualization
    Returns list of (path, node, depth, protein_count) tuples
    """
    results = []
    
    def traverse(node, path=""):
        # Count proteins in this subtree
        protein_count = len(node['proteins'])
        for child in node['children'].values():
            protein_count += count_proteins(child)
        
        if protein_count >= min_proteins:
            results.append((path, node, node['depth'], protein_count))
        
        for child_id, child_node in node['children'].items():
            traverse(child_node, path + '.' + child_id if path else child_id)
    
    traverse(tree)
    
    # Sort by protein count and take top ones
    results.sort(key=lambda x: x[3], reverse=True)
    return results[:max_children]


def count_proteins(node):
    """Count total proteins in a subtree"""
    total = len(node['proteins'])
    for child in node['children'].values():
        total += count_proteins(child)
    return total


def draw_phylogeny_tree(ax, tree, hog_cache, umap_coords, roothog_id, title, target_protein=None):
    """
    Draw a simplified phylogenetic tree showing HOG hierarchy
    """
    ax.set_xlim(-0.5, 10)
    ax.set_ylim(-0.5, 10)
    ax.axis('off')
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    
    # Get major branches
    major_branches = get_major_child_hogs(tree, max_children=10, min_proteins=5)
    
    if not major_branches:
        ax.text(5, 5, "No significant\nsub-HOGs found", ha='center', va='center', fontsize=10)
        return {}
    
    # Draw tree structure
    root_x, root_y = 0.5, 5
    
    # Draw root node
    root_circle = plt.Circle((root_x, root_y), 0.3, color=DEPTH_COLORS[0], ec='black', lw=2, zorder=5)
    ax.add_patch(root_circle)
    ax.text(root_x, root_y + 0.5, f'Root\n{roothog_id}', ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    # Position branches
    n_branches = min(len(major_branches), 8)
    y_positions = np.linspace(9, 1, n_branches)
    
    branch_colors = {}  # Map branch path to color for UMAP matching
    
    for i, (path, node, depth, pcount) in enumerate(major_branches[:n_branches]):
        y = y_positions[i]
        x = 2.5 + min(depth, 5) * 1.2
        
        # Color based on depth
        color = DEPTH_COLORS[min(depth, len(DEPTH_COLORS)-1)]
        branch_colors[path] = color
        
        # Draw connection line
        ax.plot([root_x + 0.3, 1.5], [root_y, root_y], 'k-', lw=1.5, alpha=0.6)
        ax.plot([1.5, 1.5], [root_y, y], 'k-', lw=1.5, alpha=0.6)
        ax.plot([1.5, x - 0.25], [y, y], 'k-', lw=1.5, alpha=0.6)
        
        # Draw node
        node_circle = plt.Circle((x, y), 0.25, color=color, ec='black', lw=1.5, zorder=5)
        ax.add_patch(node_circle)
        
        # Label with abbreviated path and protein count
        short_path = path if len(path) < 15 else path[:12] + '...'
        label = f'{short_path}\n({pcount} proteins)'
        ax.text(x + 0.4, y, label, ha='left', va='center', fontsize=7)
        
        # Highlight target protein branch if specified
        if target_protein and any(p['entry_id'] == target_protein for p in node['proteins']):
            ax.add_patch(plt.Circle((x, y), 0.35, fill=False, ec='gold', lw=3, zorder=4))
            ax.text(x, y - 0.5, f'← {target_protein}', ha='center', va='top', fontsize=7, 
                   color='gold', fontweight='bold')
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=DEPTH_COLORS[0], 
                  markersize=10, label='Root HOG'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=DEPTH_COLORS[1], 
                  markersize=10, label='Child Level 1'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=DEPTH_COLORS[2], 
                  markersize=10, label='Child Level 2+'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=7)
    
    return branch_colors


def draw_umap_panel(ax, umap_coords, hog_cache, roothog_id, branch_colors, title):
    """
    Draw UMAP panel with proteins colored by their HOG branch
    """
    # Filter to this root HOG
    hog_umap = umap_coords[umap_coords['roothog_id'] == roothog_id].copy()
    
    if len(hog_umap) == 0:
        ax.text(0.5, 0.5, "No UMAP data\nfor this HOG", ha='center', va='center', 
               transform=ax.transAxes, fontsize=12)
        ax.set_title(title, fontsize=12, fontweight='bold')
        return
    
    # Merge with hog_cache to get hog_id
    hog_data = hog_cache[hog_cache['roothog_id'] == roothog_id][['EntryID', 'hog_id']]
    hog_umap = hog_umap.merge(hog_data, on='EntryID', how='left')
    
    # Assign colors based on HOG branch
    colors = []
    for _, row in hog_umap.iterrows():
        hierarchy = parse_hog_hierarchy(row['hog_id'])
        if len(hierarchy) <= 1:
            colors.append(DEPTH_COLORS[0])  # Root
        else:
            # Find matching branch
            path = '.'.join(str(h) for h in hierarchy[1:])
            color_found = False
            for branch_path, color in branch_colors.items():
                if path.startswith(branch_path) or branch_path.startswith(path):
                    colors.append(color)
                    color_found = True
                    break
            if not color_found:
                colors.append('#cccccc')  # Gray for unmatched
    
    # Also color by species
    species_colors = []
    for tax in hog_umap['taxonomyID']:
        species_colors.append(SPECIES_COLORS.get(tax, '#cccccc'))
    
    # Create scatter plot
    scatter = ax.scatter(hog_umap['umap_x'], hog_umap['umap_y'], 
                        c=colors, s=20, alpha=0.7, edgecolors='white', linewidths=0.3)
    
    # Highlight human proteins
    human_mask = hog_umap['taxonomyID'] == 9606
    if human_mask.any():
        ax.scatter(hog_umap.loc[human_mask, 'umap_x'], 
                  hog_umap.loc[human_mask, 'umap_y'],
                  c='none', s=50, edgecolors='red', linewidths=1.5, 
                  marker='o', label='Human', zorder=10)
    
    ax.set_xlabel('UMAP 1', fontsize=10)
    ax.set_ylabel('UMAP 2', fontsize=10)
    ax.set_title(title, fontsize=12, fontweight='bold')
    
    # Add species legend
    legend_handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', 
                  markeredgecolor='red', markersize=8, label='Human'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#377eb8', 
                  markersize=8, label='Mouse'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff7f00', 
                  markersize=8, label='Drosophila'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#a65628', 
                  markersize=8, label='Yeast'),
    ]
    ax.legend(handles=legend_handles, loc='upper right', fontsize=7, title='Species')
    
    return hog_umap


def create_combined_figure(hog_cache, umap_coords, roothog_id, protein_name, protein_id, family_name):
    """
    Create a combined figure with UMAP + Phylogeny tree for a protein family
    """
    fig = plt.figure(figsize=(16, 8))
    
    # Create grid: [UMAP panel | Phylogeny tree panel]
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.3)
    
    ax_umap = fig.add_subplot(gs[0])
    ax_tree = fig.add_subplot(gs[1])
    
    # Build tree
    tree = build_hog_tree(hog_cache, roothog_id)
    
    # Draw phylogeny tree first to get branch colors
    branch_colors = draw_phylogeny_tree(
        ax_tree, tree, hog_cache, umap_coords, roothog_id,
        f'HOG Hierarchy for {family_name}',
        target_protein=protein_id
    )
    
    # Draw UMAP panel
    draw_umap_panel(
        ax_umap, umap_coords, hog_cache, roothog_id, branch_colors,
        f'UMAP Embedding: {protein_name} ({protein_id})\n{family_name} - HOG {roothog_id}'
    )
    
    # Add overall title
    fig.suptitle(f'Protein Family Analysis: {family_name}', fontsize=14, fontweight='bold', y=1.02)
    
    # Add annotation
    fig.text(0.5, -0.02, 
             f'Human reference protein: {protein_name} ({protein_id}) | '
             f'Total proteins in family: {count_proteins(tree)} across {len(set(hog_cache[hog_cache["roothog_id"]==roothog_id]["oma_id"].str[:5]))} species',
             ha='center', va='top', fontsize=9, style='italic')
    
    return fig


def create_dual_family_figure(hog_cache, umap_coords):
    """
    Create a figure showing both protein families
    """
    fig = plt.figure(figsize=(18, 14))
    
    # Create 2x2 grid
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1], height_ratios=[1, 1], 
                         wspace=0.25, hspace=0.35)
    
    # Family 1: DNAJB1 / HSP40 (HOG 801468)
    ax1_umap = fig.add_subplot(gs[0, 0])
    ax1_tree = fig.add_subplot(gs[0, 1])
    
    tree1 = build_hog_tree(hog_cache, 801468)
    branch_colors1 = draw_phylogeny_tree(
        ax1_tree, tree1, hog_cache, umap_coords, 801468,
        'HOG Hierarchy: Heat Shock Protein 40 Family',
        target_protein='P25685'
    )
    draw_umap_panel(
        ax1_umap, umap_coords, hog_cache, 801468, branch_colors1,
        'UMAP: DNAJB1 (P25685) - HSP40 Co-chaperone Family\nHOG 801468 | 2230 proteins | 77 taxa'
    )
    
    # Family 2: IGF2R (HOG 792940)
    ax2_umap = fig.add_subplot(gs[1, 0])
    ax2_tree = fig.add_subplot(gs[1, 1])
    
    tree2 = build_hog_tree(hog_cache, 792940)
    branch_colors2 = draw_phylogeny_tree(
        ax2_tree, tree2, hog_cache, umap_coords, 792940,
        'HOG Hierarchy: IGF2 Receptor Family',
        target_protein='P11717'
    )
    draw_umap_panel(
        ax2_umap, umap_coords, hog_cache, 792940, branch_colors2,
        'UMAP: IGF2R (P11717) - Mannose 6-Phosphate Receptor Family\nHOG 792940 | 504 proteins | 21 taxa'
    )
    
    # Add overall title
    fig.suptitle('Phylogeny-Enhanced HOG Analysis\nHuman Protein Families with Hierarchical Orthologous Groups', 
                fontsize=14, fontweight='bold', y=0.98)
    
    # Add legend/explanation
    fig.text(0.5, 0.01, 
             'Left panels: UMAP embedding colored by HOG hierarchy depth | '
             'Right panels: Simplified phylogenetic tree of major HOG branches\n'
             'Red circles highlight human proteins | Node colors match between UMAP and tree',
             ha='center', va='bottom', fontsize=9, style='italic', 
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    return fig


def main():
    """Main execution"""
    print("=" * 60)
    print("Phylogeny-Enhanced HOG Visualization")
    print("=" * 60)
    
    # Load data
    hog_cache, umap_coords = load_data()
    
    # Create dual family figure
    print("\nCreating combined dual-family figure...")
    fig_dual = create_dual_family_figure(hog_cache, umap_coords)
    
    output_file = OUTPUT_PATH / "phylogeny_enhanced_hog.png"
    fig_dual.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_file}")
    
    # Create individual figures for each family
    print("\nCreating individual family figures...")
    
    # HSP40 family
    fig_hsp = create_combined_figure(
        hog_cache, umap_coords, 801468,
        'DNAJB1', 'P25685', 'Heat Shock Protein 40 (HSP40) Co-chaperone Family'
    )
    output_hsp = OUTPUT_PATH / "phylogeny_hsp40_family.png"
    fig_hsp.savefig(output_hsp, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_hsp}")
    
    # IGF2R family
    fig_igf = create_combined_figure(
        hog_cache, umap_coords, 792940,
        'IGF2R', 'P11717', 'Mannose 6-Phosphate / IGF2 Receptor Family'
    )
    output_igf = OUTPUT_PATH / "phylogeny_igf2r_family.png"
    fig_igf.savefig(output_igf, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_igf}")
    
    plt.close('all')
    
    # Print summary
    print("\n" + "=" * 60)
    print("Summary of Generated Visualizations:")
    print("=" * 60)
    print(f"1. {output_file.name}")
    print("   - Combined figure with both protein families")
    print("   - Shows UMAP + phylogeny tree for each")
    print(f"\n2. {output_hsp.name}")
    print("   - DNAJB1 (P25685) - HSP40 family")
    print("   - 2230 proteins across 77 taxa")
    print(f"\n3. {output_igf.name}")
    print("   - IGF2R (P11717) - M6P/IGF2 receptor family")
    print("   - 504 proteins across 21 taxa")
    print("=" * 60)


if __name__ == "__main__":
    main()
