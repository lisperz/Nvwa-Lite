

| 🚀  Nvwa MVP 1.0 Quick Start Guide  —  Beta Pilot *Welcome — we're thrilled to have you as a pilot user. This guide covers everything you need to get started.* |
| :---- |

# **📂  1\. Getting Started**

* **Access:** Use the unique link provided to log in to your private workspace.

* **Data format:** We currently support processed .h5ad (AnnData) files — post-analysis output from Scanpy or Seurat pipelines. Nvwa interprets your results; it does not re-run raw data processing or FASTQ pipelines.

* **Privacy:** Your data runs in an isolated, encrypted AWS container. All uploaded data is permanently and automatically deleted 32 hours after your session ends.

| 📋  File Requirements Your .h5ad file should already contain completed clustering, UMAP embeddings, and ideally cell type annotations. Nvwa is designed for the last mile of your analysis — interpreting and visualizing results you have already generated. |
| :---- |

# **🛠️  2\. What Nvwa Can Do For You (V1.0)**

Nvwa is optimized for single-dataset exploration and biological interpretation:

* **Automated QC:** Instantly check mitochondrial content, gene distributions, and cell quality metrics with interpretive feedback.

* **Clustering & UMAP:** Visualize cell populations and cluster structure from your pre-computed embeddings.

* **Marker Gene Discovery:** Identify and tabulate the top distinguishing genes for any cluster or cell type.

* **Differential Expression Analysis:** Compare gene expression across clusters using three modes: all clusters vs. all, a single cluster vs. the rest, or pairwise between two specific clusters. Each result includes biological interpretation.

* **Gene Expression Visualization:** Generate high-quality Violin, Dot, Feature, and Heatmap plots on demand for any gene or gene set.

* **Biological Interpretation:** After each analysis, Nvwa provides biological interpretation — you can ask follow-up questions to go deeper.

* **Export:** Download a comprehensive bundle containing figures, DE result tables, and analysis logs.

# **💡  3\. Try These Prompts**

Copy and paste these examples to see Nvwa in action:

| Initial Audit | *Run standard QC checks and show me the UMAP clusters.* |
| :---- | :---- |

| Quality Check | *Are there any clusters with unusually high mitochondrial content? Show me the QC plots.* |
| :---- | :---- |

| Marker Search | *Find the top 10 marker genes for each cell type and provide the results in a table.* |
| :---- | :---- |

| Gene Visualization | *Show me a DotPlot of CD3D, CD8A, and CD14 expression across all clusters.* |
| :---- | :---- |

| Biological Interpretation | *Based on the markers, what cell types do you think Cluster 0 and Cluster 1 represent, and what is their likely biological significance?* |
| :---- | :---- |

| Differential Expression | *Compare Cluster 2 vs Cluster 5 — what are the top differentially expressed genes and what does it suggest biologically?* |
| :---- | :---- |

# **⚠️  4\. Important Constraints**

| Out of Scope (V1.0) Multi-sample integration or batch correction Trajectory inference / RNA velocity Raw FASTQ or unprocessed count matrix ingestion Spatial transcriptomics (coming post-V1.0) |
| :---- |

* **The Gatekeeper:** If Nvwa issues a WARN or FAIL message, please read the feedback carefully. It is designed to prevent over-interpretation of low-quality data — this is a feature, not a bug.

* **Concurrency:** To ensure stability during beta, each user is limited to one active run at a time.

* **Reproducibility:** Nvwa is designed for consistent outputs. Temperature is set to 0 to minimize variability.

# **📩  5\. Your Feedback Is Our Roadmap**

Nvwa is built by scientists, for scientists. As a pilot user, your experience directly shapes the product. Please tell us immediately if:

* **Something breaks:** or produces an unexpected result

* **A biological interpretation seems off:** your domain expertise is the ground truth

* **You wish Nvwa could do X:** your feature requests become our next sprint

| 📬  Reach Yalu Directly Email: yalu@nvwabio.com *Response time: within 24 hours. For urgent issues during active sessions, contact directly via the channel shared with you at onboarding.* |
| :---- |

***Enjoy exploring your data.***

*—  The Nvwa Team*