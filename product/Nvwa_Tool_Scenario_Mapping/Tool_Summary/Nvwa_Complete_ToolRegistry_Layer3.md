# Nvwa Bio — Complete Tool Registry
*April 2026 | All features | Based on validated scenarios*

---

## Global Setup

### `set_figure_params`
```python
sc.settings.set_figure_params(dpi=150, facecolor='white')
# Do NOT set figsize globally
```

### `data_ingestion`
```python
# Ensure adata.X is logcounts
adata.X = adata.layers['logcounts']

# Fix UMAP coordinates
adata.obsm['X_umap'] = np.array(adata.obsm['X_umap'])

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
```

---

## Condition Color Standard

Colors are assigned based on user-defined reference and treatment order, not hardcoded.

```python
def get_condition_colors(conditions, reference):
    # reference condition always gets blue, others get red/green/purple in order
    color_palette = ['#2E86C1', '#E74C3C', '#27AE60', '#8E44AD']
    ordered = [reference] + [c for c in conditions if c != reference]
    return dict(zip(ordered, color_palette[:len(ordered)]))
```

| Order | Role | Color | Hex |
|--|--|--|--|
| Reference / Control | — | Blue | `#2E86C1` |
| Treatment / Experimental | — | Red | `#E74C3C` |
| 3rd condition | — | Green | `#27AE60` |
| 4th condition | — | Purple | `#8E44AD` |

---

## Data Operation Tools

### `subset_data` — single value
```python
adata_sub = adata[adata.obs[col] == value]
```

### `subset_data` — multiple values
```python
adata_sub = adata[adata.obs[col].isin(values)]
```

### `create_combined_column`
```python
adata.obs['celltype_condition'] = (
    adata.obs[col1].astype(str) + '_' +
    adata.obs[col2].astype(str)
)
```

### `create_group_column`
- **Parameters:** `adata, condition_col, group_map`

```python
adata.obs['de_group'] = 'other'
for group_name, conditions in group_map.items():
    mask = adata.obs[condition_col].isin(conditions)
    adata.obs.loc[mask, 'de_group'] = group_name
```

---

## UMAP Tools

### `umap_plot`
- **Parameters:** `adata, CELLTYPE_COL`

```python
sc.pl.umap(
    adata,
    color=[CELLTYPE_COL],
    title='UMAP colored by cell type',
    show=False,
    frameon=False
)
plt.savefig('fig_umap.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

### `umap_plot_split`
- **Parameters:** `adata, CELLTYPE_COL, CONDITION_COL, reference`
- **Figsize:** `(5 * n_conditions, 5)`

```python
conditions = adata.obs[CONDITION_COL].unique()
n = len(conditions)
fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))

x_min, x_max = adata.obsm['X_umap'][:, 0].min(), adata.obsm['X_umap'][:, 0].max()
y_min, y_max = adata.obsm['X_umap'][:, 1].min(), adata.obsm['X_umap'][:, 1].max()

for ax, c in zip(axes, conditions):
    sc.pl.umap(
        adata[adata.obs[CONDITION_COL] == c],
        color=[CELLTYPE_COL],
        ax=ax,
        show=False,
        title=f'{c}',
        frameon=False
    )
    ax.set_xlim(x_min - 0.5, x_max + 0.5)
    ax.set_ylim(y_min - 0.5, y_max + 0.5)

plt.tight_layout()
plt.savefig('fig_umap_split.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

---

## QC Tools

### `qc_violin_plot`
- **Parameters:** `adata, CELLTYPE_COL, qc_metrics=['nFeature_RNA', 'nCount_RNA', 'pct_counts_mt']`

```python
sc.pl.violin(
    adata,
    keys=qc_metrics,
    groupby=CELLTYPE_COL,
    rotation=45,
    show=False,
    multi_panel=True
)
for ax in plt.gcf().get_axes():
    ax.set_xticklabels(ax.get_xticklabels(), ha='right')
plt.savefig('fig_qc_violin.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

### `qc_violin_plot_split`
- **Parameters:** `adata, CELLTYPE_COL, CONDITION_COL, qc_metrics=['nFeature_RNA', 'nCount_RNA', 'pct_counts_mt']`

```python
conditions = adata.obs[CONDITION_COL].unique()

for c in conditions:
    sc.pl.violin(
        adata[adata.obs[CONDITION_COL] == c],
        keys=qc_metrics,
        groupby=CELLTYPE_COL,
        rotation=45,
        show=False,
        multi_panel=True
    )
    fig = plt.gcf()
    for ax in fig.get_axes():
        ax.set_xticklabels(ax.get_xticklabels(), ha='right')
    fig.suptitle(c, fontsize=14, fontweight='bold', y=1.02)
    plt.savefig(f'fig_qc_violin_{c}.png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.show()
```

### `qc_summary_table`
- **Parameters:** `adata, qc_metrics=['nFeature_RNA', 'nCount_RNA', 'pct_counts_mt'], pct_mt_threshold=10, n_feature_threshold=200`

```python
summary = adata.obs[qc_metrics].agg(['min', 'max', 'mean', 'median']).T
summary.columns = ['Min', 'Max', 'Mean', 'Median']
summary = summary.round(2)

flags = []
for metric in qc_metrics:
    flag = ''
    if metric == 'pct_counts_mt' and summary.loc[metric, 'Median'] > pct_mt_threshold:
        flag = f'⚠️ Elevated mitochondrial fraction (median > {pct_mt_threshold}%). May indicate cell stress.'
    elif metric == 'nFeature_RNA' and summary.loc[metric, 'Median'] < n_feature_threshold:
        flag = f'⚠️ Low gene detection (median < {n_feature_threshold}). Consider reviewing filtering thresholds.'
    flags.append(flag)

summary['Flag'] = flags
print(summary)
```

---

## Feature Plot Tools

### `feature_plot`
- **Parameters:** `gene, adata`
- **Colormap:** `Reds`

```python
sc.pl.umap(
    adata,
    color=[gene],
    color_map='Reds',
    title=f'{gene}',
    show=False,
    frameon=False
)
plt.savefig(f'fig_feature_plot_{gene}.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

### `feature_plot_split`
- **Parameters:** `gene, adata, CONDITION_COL`
- **Figsize:** `(5 * n_conditions, 5)`
- **Colormap:** `Reds`
- **Internal:** `get_expression_range`, `get_umap_range`

```python
# get_expression_range
gene_idx = adata.var_names.get_loc(gene)
vmin = adata.X[:, gene_idx].min()
vmax = adata.X[:, gene_idx].max()

# get_umap_range
x_min, x_max = adata.obsm['X_umap'][:, 0].min(), adata.obsm['X_umap'][:, 0].max()
y_min, y_max = adata.obsm['X_umap'][:, 1].min(), adata.obsm['X_umap'][:, 1].max()

conditions = adata.obs[CONDITION_COL].unique()
n = len(conditions)
fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))

for ax, c in zip(axes, conditions):
    sc.pl.umap(
        adata[adata.obs[CONDITION_COL] == c],
        color=[gene],
        color_map='Reds',
        vmin=vmin,
        vmax=vmax,
        ax=ax,
        show=False,
        title=f'{gene} — {c}'
    )
    ax.set_xlim(x_min - 0.5, x_max + 0.5)
    ax.set_ylim(y_min - 0.5, y_max + 0.5)

plt.tight_layout()
plt.savefig(f'fig_feature_plot_split_{gene}.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

---

## Violin Plot Tools

### `violin_plot`
- **Parameters:** `gene, adata, groupby`
- **Figsize:** `(8, 4)`

```python
fig, ax = plt.subplots(figsize=(8, 4))
sc.pl.violin(
    adata,
    keys=[gene],
    groupby=groupby,
    rotation=45,
    ax=ax,
    show=False
)
ax.set_title(f'{gene} across {groupby}')
ax.set_xlabel('')
ax.set_xticklabels(ax.get_xticklabels(), ha='right')
plt.savefig(f'fig_violin_{gene}.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

### `violin_plot_split`
- **Parameters:** `gene, adata, CELLTYPE_COL, CONDITION_COL`
- **Figsize:** `(15, 6)`

```python
df = pd.DataFrame({
    'expression': adata[:, gene].X.toarray().flatten(),
    'celltype': adata.obs[CELLTYPE_COL].values,
    'condition': adata.obs[CONDITION_COL].values
})

fig, ax = plt.subplots(figsize=(15, 6))
sns.violinplot(
    data=df,
    x='celltype',
    y='expression',
    hue='condition',
    ax=ax,
    inner='box',
    scale='width',
    split=False
)
ax.set_ylim(bottom=0)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
ax.set_xlabel('')
ax.set_title(f'{gene} by cell type and condition')
plt.tight_layout()
plt.savefig(f'fig_violin_split_{gene}.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

---

## Dot Plot Tools

### `dot_plot`
- **Parameters:** `genes, adata, groupby`

```python
sc.pl.dotplot(
    adata,
    var_names=genes,
    groupby=groupby,
    show=False
)
plt.savefig('fig_dotplot.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

### `dot_plot_top_markers`
- **Parameters:** `adata, markers_df, CELLTYPE_COL, n_genes=5`

```python
celltype_order = adata.obs[CELLTYPE_COL].cat.categories.tolist()

top_genes_ordered = []
for ct in celltype_order:
    genes = (
        markers_df[markers_df['group'] == ct]
        .nlargest(n_genes, 'logfoldchanges')['names']
        .tolist()
    )
    top_genes_ordered.extend(genes)

sc.pl.dotplot(
    adata,
    var_names=top_genes_ordered,
    groupby=CELLTYPE_COL,
    show=False
)
plt.savefig('fig_dotplot_markers.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

### `color_yaxis_by_condition` (internal)
- **Parameters:** `ax, CONDITION_COL, adata, reference`

```python
condition_colors = get_condition_colors(adata.obs[CONDITION_COL].unique(), reference)
mainax = ax.get_axes()['mainplot_ax']

for label in mainax.get_yticklabels():
    text = label.get_text()
    for cond, color in condition_colors.items():
        if text.endswith('_' + cond):
            label.set_color(color)
            break
```

---

## Heatmap Tools

### `heatmap_plot`
- **Parameters:** `genes, adata, groupby`
- **Colormap:** `viridis`

```python
sc.pl.heatmap(
    adata,
    var_names=genes,
    groupby=groupby,
    show_gene_labels=True,
    cmap='viridis',
    swap_axes=True,
    show=False
)
for ax in plt.gcf().get_axes():
    ax.tick_params(axis='x', labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
plt.savefig('fig_heatmap.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

### `heatmap_plot_split`
- **Parameters:** `genes, adata, CELLTYPE_COL, CONDITION_COL, reference`
- **Colormap:** `viridis`
- **Internal:** `create_combined_column`, `color_xaxis_by_condition`

```python
if 'celltype_condition' not in adata.obs.columns:
    adata.obs['celltype_condition'] = (
        adata.obs[CELLTYPE_COL].astype(str) + '_' +
        adata.obs[CONDITION_COL].astype(str)
    )

sc.pl.heatmap(
    adata,
    var_names=genes,
    groupby='celltype_condition',
    show_gene_labels=True,
    cmap='viridis',
    swap_axes=True,
    show=False
)

fig = plt.gcf()
condition_colors = get_condition_colors(adata.obs[CONDITION_COL].unique(), reference)

groupby_ax = None
for ax in fig.get_axes():
    labels = [l.get_text() for l in ax.get_xticklabels()]
    if any('_' in l for l in labels):
        groupby_ax = ax
        break

if groupby_ax:
    for label in groupby_ax.get_xticklabels():
        text = label.get_text()
        for cond, color in condition_colors.items():
            if text.endswith('_' + cond):
                label.set_color(color)
                break
        label.set_fontsize(6)

for ax in fig.get_axes():
    if ax != groupby_ax:
        ax.tick_params(axis='y', labelsize=6)

plt.savefig('fig_heatmap_split.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

---

## Find Markers Tools

### `find_all_markers`
- **Parameters:** `adata, CELLTYPE_COL, method='wilcoxon', pvals_adj_threshold=0.05, logfc_threshold=0, n_top_genes=10`

```python
sc.tl.rank_genes_groups(
    adata,
    groupby=CELLTYPE_COL,
    method=method,
    pts=True
)

markers_df = sc.get.rank_genes_groups_df(adata, group=None)
markers_filtered = markers_df[
    (markers_df['pvals_adj'] < pvals_adj_threshold) &
    (markers_df['logfoldchanges'] > logfc_threshold)
]

for celltype, group in markers_filtered.groupby('group'):
    top_genes = group.nlargest(n_top_genes, 'logfoldchanges')['names'].tolist()
    print(f"{celltype}: {', '.join(top_genes)}")

markers_filtered.to_csv('markers_all_celltypes.csv', index=False)
print('\nFull results saved to markers_all_celltypes.csv')
```

### `find_markers`
- **Parameters:** `adata, CELLTYPE_COL, celltype, method='wilcoxon', pvals_adj_threshold=0.05, logfc_threshold=0, n_top_genes=10`

```python
sc.tl.rank_genes_groups(
    adata,
    groupby=CELLTYPE_COL,
    groups=[celltype],
    reference='rest',
    method=method,
    pts=True
)

markers_df = sc.get.rank_genes_groups_df(adata, group=celltype)
markers_filtered = markers_df[
    (markers_df['pvals_adj'] < pvals_adj_threshold) &
    (markers_df['logfoldchanges'] > logfc_threshold)
]

top_genes = markers_filtered.nlargest(n_top_genes, 'logfoldchanges')['names'].tolist()
print(f"{celltype}: {', '.join(top_genes)}")

markers_filtered.to_csv(f'markers_{celltype}.csv', index=False)
print(f'\nFull results saved to markers_{celltype}.csv')
```

---

## DE Analysis Tools

### `run_de`
- **Parameters:** `adata, groupby_col, group1, group2, method='wilcoxon', pvals_adj_threshold=0.05, n_top_genes=10`

```python
sc.tl.rank_genes_groups(
    adata,
    groupby=groupby_col,
    groups=[group1],
    reference=group2,
    method=method,
    pts=True
)

de_df = sc.get.rank_genes_groups_df(adata, group=group1)
de_filtered = de_df[de_df['pvals_adj'] < pvals_adj_threshold]

n_total = len(de_filtered)
n_up_g1 = len(de_filtered[de_filtered['logfoldchanges'] > 0])
n_up_g2 = len(de_filtered[de_filtered['logfoldchanges'] < 0])

print(f"DE Summary: {group1} vs {group2}")
print(f"Total significantly changed genes (adj_p < {pvals_adj_threshold}): {n_total}")
print(f"Upregulated in {group1}: {n_up_g1}")
print(f"Upregulated in {group2}: {n_up_g2}")

print(f"\nTop {n_top_genes} upregulated in {group1}:")
top_g1 = de_filtered[de_filtered['logfoldchanges'] > 0].nlargest(n_top_genes, 'logfoldchanges')
print(top_g1[['names', 'logfoldchanges', 'pvals_adj']].to_string(index=False))

print(f"\nTop {n_top_genes} upregulated in {group2}:")
top_g2 = de_filtered[de_filtered['logfoldchanges'] < 0].nsmallest(n_top_genes, 'logfoldchanges')
print(top_g2[['names', 'logfoldchanges', 'pvals_adj']].to_string(index=False))

de_filtered.to_csv(f'de_{group1}_vs_{group2}.csv', index=False)
print(f'\nFull results saved to de_{group1}_vs_{group2}.csv')
print(f'\nNote: Results based on {method} test (single-cell level). For publication, consider pseudo-bulk analysis with DESeq2.')
```

### `volcano_plot`
- **Parameters:** `de_df, group1, group2, pvals_adj_threshold=0.05, lfc_threshold=0.25, lfc_clip=6, n_label=10`

```python
from adjustText import adjust_text

de_results = de_df.copy()
de_results['neg_log10_pval'] = -np.log10(de_results['pvals_adj'].clip(lower=1e-300))
de_results['logfoldchanges_clipped'] = de_results['logfoldchanges'].clip(-lfc_clip, lfc_clip)

de_results['color'] = 'lightgrey'
de_results.loc[
    (de_results['pvals_adj'] < pvals_adj_threshold) & (de_results['logfoldchanges'] > lfc_threshold), 'color'
] = '#E74C3C'
de_results.loc[
    (de_results['pvals_adj'] < pvals_adj_threshold) & (de_results['logfoldchanges'] < -lfc_threshold), 'color'
] = '#2E86C1'

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(
    de_results['logfoldchanges_clipped'],
    de_results['neg_log10_pval'],
    c=de_results['color'],
    alpha=0.5,
    s=10
)
ax.axhline(-np.log10(pvals_adj_threshold), color='black', linestyle='--', linewidth=0.8)
ax.axvline(lfc_threshold, color='grey', linestyle='--', linewidth=0.8)
ax.axvline(-lfc_threshold, color='grey', linestyle='--', linewidth=0.8)

top_up = de_results[de_results['color'] == '#E74C3C'].nlargest(n_label, 'neg_log10_pval')
top_down = de_results[de_results['color'] == '#2E86C1'].nlargest(n_label, 'neg_log10_pval')
top_genes = pd.concat([top_up, top_down])

texts = []
for _, row in top_genes.iterrows():
    texts.append(ax.text(
        row['logfoldchanges_clipped'],
        row['neg_log10_pval'],
        row['names'],
        fontsize=7
    ))
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))

ax.set_xlabel('Log2 Fold Change')
ax.set_ylabel('-log10(adj p-value)')
ax.set_title(f'Volcano Plot: {group1} vs {group2}')
plt.tight_layout()
plt.savefig(f'fig_volcano_{group1}_vs_{group2}.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

---

## Cell Type Composition Tools

### `composition_barplot`
- **Parameters:** `adata, CELLTYPE_COL, CONDITION_COL, mode, reference, celltypes=None`
- **Modes:** `'count'`, `'proportion_by_celltype'`, `'proportion_by_condition'`

```python
if celltypes:
    adata = adata[adata.obs[CELLTYPE_COL].isin(celltypes)].copy()

count_df = adata.obs.groupby([CELLTYPE_COL, CONDITION_COL]).size().reset_index(name='count')
celltypes_list = adata.obs[CELLTYPE_COL].cat.categories.tolist()
conditions = list(adata.obs[CONDITION_COL].unique())
condition_colors = get_condition_colors(conditions, reference)

if mode == 'count':
    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(celltypes_list))
    width = 0.35
    for i, cond in enumerate(conditions):
        counts = [
            count_df[(count_df[CELLTYPE_COL] == ct) & (count_df[CONDITION_COL] == cond)]['count'].values[0]
            if len(count_df[(count_df[CELLTYPE_COL] == ct) & (count_df[CONDITION_COL] == cond)]) > 0
            else 0 for ct in celltypes_list
        ]
        ax.bar(x + i * width, counts, width, label=cond, color=condition_colors[cond], alpha=0.8)
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(celltypes_list, rotation=45, ha='right')
    ax.set_ylabel('Cell Count')
    ax.set_title('Cell number by cell type and condition')
    ax.legend()
    ax.set_ylim(0, ax.get_ylim()[1] * 1.1)

elif mode == 'proportion_by_celltype':
    count_df['proportion'] = count_df.groupby(CELLTYPE_COL)['count'].transform(lambda x: x / x.sum())
    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(celltypes_list))
    width = 0.35
    for i, cond in enumerate(conditions):
        proportions = [
            count_df[(count_df[CELLTYPE_COL] == ct) & (count_df[CONDITION_COL] == cond)]['proportion'].values[0]
            if len(count_df[(count_df[CELLTYPE_COL] == ct) & (count_df[CONDITION_COL] == cond)]) > 0
            else 0 for ct in celltypes_list
        ]
        ax.bar(x + i * width, proportions, width, label=cond, color=condition_colors[cond], alpha=0.8)
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(celltypes_list, rotation=45, ha='right')
    ax.set_ylabel('Proportion')
    ax.set_title('Proportion of each cell type across conditions')
    ax.legend()
    ax.set_ylim(0, 1.1)

elif mode == 'proportion_by_condition':
    count_df['proportion'] = count_df.groupby(CONDITION_COL)['count'].transform(lambda x: x / x.sum())
    celltype_colors = dict(zip(celltypes_list, sc.pl.palettes.default_20[:len(celltypes_list)]))
    n_celltypes = len(celltypes_list)
    fig_height = max(4, n_celltypes * 0.3)
    fig, ax = plt.subplots(figsize=(4, fig_height))
    x = np.arange(len(conditions))
    bottom = np.zeros(len(conditions))
    for ct in celltypes_list:
        proportions = [
            count_df[(count_df[CELLTYPE_COL] == ct) & (count_df[CONDITION_COL] == cond)]['proportion'].values[0]
            if len(count_df[(count_df[CELLTYPE_COL] == ct) & (count_df[CONDITION_COL] == cond)]) > 0
            else 0 for cond in conditions
        ]
        ax.bar(x, proportions, bottom=bottom, label=ct, color=celltype_colors[ct], alpha=0.8, width=0.5)
        bottom += np.array(proportions)
    ax.set_xticks(x)
    ax.set_xticklabels(conditions)
    ax.set_xlim(-0.5, len(conditions) - 0.5)
    ax.set_ylabel('Proportion')
    ax.set_title('Cell type composition per condition')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.set_ylim(0, 1.05)

plt.tight_layout()
plt.savefig(f'fig_composition_{mode}.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.show()
```

---

## Pending Tools

| Tool | Reason |
|--|--|
| `dot_plot_split` | Seurat-style split dot plot with condition-colored dots — requires custom matplotlib implementation |

---
