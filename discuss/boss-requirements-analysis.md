# Boss's Testing Requirements Analysis

## Document Overview
The boss provided a comprehensive testing framework with 4 levels of complexity for scRNA-seq visualization.

## Testing Levels

### 🟢 Level 1: Basic Visualization
**Purpose**: Check if Agent knows basic Seurat plotting functions (DimPlot, VlnPlot, FeaturePlot)

1. **Global Overview (UMAP)**
   - Prompt: "Please visualize the cell clusters using a UMAP plot."
   - Expected: `DimPlot(pbmc, reduction = "umap")`
   - **Our equivalent**: `umap_plot` tool with cluster coloring

2. **Single Gene Expression (Feature Plot)**
   - Prompt: "Show me the expression of the gene 'MS4A1' on the UMAP."
   - Expected: `FeaturePlot(pbmc, features = "MS4A1")`
   - **Our equivalent**: `feature_plot` tool

3. **Violin Plot**
   - Prompt: "Generate a violin plot for the genes 'CD3E' and 'CD8A' across all clusters."
   - Expected: `VlnPlot(pbmc, features = c("CD3E", "CD8A"))`
   - **Our equivalent**: `violin_plot` tool (currently single gene only - NEEDS ENHANCEMENT)

### 🟡 Level 2: Multi-Feature & Comparison
**Purpose**: Check if Agent can handle lists and complex chart forms

1. **Dot Plot**
   - Prompt: "Create a dot plot to show the expression of these markers: MS4A1, GNLY, CD3E, CD14, FCER1A, FCGR3A, LYZ, PPBP, CD8A."
   - Expected: `DotPlot(pbmc, features = c(...))`
   - **Our equivalent**: `dotplot` tool ✅

2. **Heatmap**
   - Prompt: "Generate a heatmap for the top 10 marker genes of each cluster."
   - Simplified: "Draw a heatmap for these genes: MS4A1, CD79A, CD3E, CD8A, GNLY, NKG7."
   - Expected: `DoHeatmap(pbmc, features = ...)`
   - **Our equivalent**: `heatmap_plot` tool ✅

3. **Correlation Scatter Plot**
   - Prompt: "Plot the correlation between 'CD3E' and 'CD8A' expression."
   - Expected: `FeatureScatter(pbmc, feature1 = "CD3E", feature2 = "CD8A")`
   - **Our equivalent**: ❌ MISSING - Need to add scatter plot tool

### 🔴 Level 3: Customization
**Purpose**: Check if Agent can modify parameters - the essence of "conversational analysis"

1. **Split By Group**
   - Prompt: "Can you show the UMAP again, but split the view by cluster? I want to see each cluster in a separate panel."
   - Expected: `DimPlot(pbmc, split.by = "seurat_clusters")`
   - **Our equivalent**: ❌ MISSING - Need to add split parameter to umap_plot

2. **Styling (Colors & Labels)**
   - Prompt: "Draw the UMAP plot again, but this time label the clusters directly on the plot and remove the legend."
   - Expected: `DimPlot(pbmc, label = TRUE) + NoLegend()`
   - **Our equivalent**: ❌ MISSING - Need to add label/legend parameters

3. **Side-by-Side Plots (Patchwork)**
   - Prompt: "Please show me the UMAP plot and the Violin plot for 'CD14' side by side."
   - Expected: `library(patchwork); plot1 + plot2`
   - **Our equivalent**: ❌ MISSING - Currently shows plots sequentially, not side-by-side

### 🔥 Level 4: Reasoning + Plotting
**Purpose**: Check if Agent understands biological context and can convert natural language to code logic

1. **Find & Visualize**
   - Prompt: "Which cluster has the highest expression of 'Lysozyme' (LYZ)? Please highlight that specific cluster on a UMAP."
   - Logic: Agent needs to calculate AverageExpression, find cluster ID, then use `DimPlot(..., cells.highlight = ...)`
   - **Our equivalent**: ❌ MISSING - Need reasoning + highlighting capability

2. **Cell Type Annotation**
   - Prompt: "Based on the expression of MS4A1, which cluster represents B cells? Rename that cluster to 'B cells' and show the UMAP with the new label."
   - Logic: `RenameIdents(pbmc, 'Old_ID' = 'B cells')` -> `DimPlot()`
   - **Our equivalent**: ❌ MISSING - Need cluster renaming capability

## Current System Gaps

### Critical Issues
1. **Hardcoded immune markers** - System assumes all datasets have CD3E, MS4A1, etc.
2. **No tissue-specific marker detection** - Can't adapt to brain, liver, kidney datasets
3. **Limited customization** - Can't modify plot parameters (labels, colors, split views)
4. **No reasoning tools** - Can't calculate which cluster has highest expression
5. **No annotation tools** - Can't rename clusters based on marker expression

### Missing Tools
1. `scatter_plot` - For gene-gene correlation
2. `calculate_average_expression` - For reasoning about cluster characteristics
3. `highlight_cells` - For emphasizing specific clusters
4. `rename_clusters` - For cell type annotation
5. Enhanced `violin_plot` - Support multiple genes
6. Enhanced `umap_plot` - Support split.by, label, legend parameters

### System Prompt Issues
1. Hardcoded marker genes don't generalize to non-immune datasets
2. No mechanism to detect tissue type and suggest appropriate markers
3. No guidance for handling Ensembl IDs vs gene symbols

## Recommended Implementation Priority

### P0 (Critical - Fixes Current Bugs)
1. **Dynamic marker gene detection** - Remove hardcoded immune markers
2. **Gene ID format handling** - Detect and handle Ensembl IDs
3. **Tissue-agnostic interpretation** - Don't assume immune cells

### P1 (High - Enables Level 3 Testing)
1. Add plot customization parameters (label, split.by, colors)
2. Support multiple genes in violin_plot
3. Add scatter_plot tool

### P2 (Medium - Enables Level 4 Testing)
1. Add calculate_average_expression tool
2. Add highlight_cells capability
3. Add rename_clusters tool
4. Enhance agent reasoning for "find highest expression" queries

### P3 (Nice to Have)
1. Side-by-side plot layout
2. Advanced styling options
3. Export to publication-ready formats

## Next Steps

1. **Immediate**: Fix the hardcoded marker gene issue
2. **Short-term**: Add missing Level 1-2 capabilities
3. **Medium-term**: Add Level 3 customization features
4. **Long-term**: Add Level 4 reasoning capabilities
