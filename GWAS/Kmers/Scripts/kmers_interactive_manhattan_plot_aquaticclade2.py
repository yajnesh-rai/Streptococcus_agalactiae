import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import plotly.io as pio


# Load GWAS_kmers CSV file with blast annotations
script_dir = os.path.dirname(os.path.abspath(__file__))  # Folder where script is located
data_folder = script_dir  # Assuming CSV is in the same folder
csv_file = os.path.join(data_folder, "kmer_results_fdr_with_blast_gene_aquaticclade2.csv")

df = pd.read_csv(csv_file)

df = df[df['genomic_pos'] >= 0].copy()
df['genomic_pos'] = pd.to_numeric(df['genomic_pos'])
df['neg_log10_p'] = -np.log10(df['lrt-pval-FDR'])


# Color based on host group association
df['strain_type'] = df['beta'].apply(lambda b: "Clinical" if b > 0 else "Aquatic")
point_colors = {"Clinical": "red", "Aquatic": "blue"}

# Plotly figure
fig = go.Figure()

# Add scatter points
for stype, color in point_colors.items():
    subset = df[df['strain_type']==stype]
    fig.add_trace(go.Scatter(
        x=subset['genomic_pos'],
        y=subset['neg_log10_p'],
        mode='markers',
        marker=dict(color=color, size=7),
        name=stype,
        hovertemplate=
            "Genomic pos: %{x}<br>" +
            "-log10 p-FDR: %{y}<br>" +
            "Gene: %{customdata[0]}<br>" +
            "Product: %{customdata[1]}<br>" +
            "Beta: %{customdata[2]}",
        customdata=subset[['gene','product','beta']].values
    ))

# Highlight genomic regions 
regions = {
    "Genomic islands": {"color":"cyan", "coords":[(1990501,2000871)]},
    "Prophages": {"color":"yellow", "coords":[(620951, 677171)]},
    "Insertion sequence IS1381": {"color":"green", "coords":[(155553, 156412),
   (985461, 986320), (1139127, 1139986), (1505615, 1506474), (1903804, 1904663),
   (1907788, 1908647)]},
    "group-II-introns": {"color":"violet", "coords":[(96917, 98776),
   (101976, 103832), (521954, 523810),
   (536871, 538729), (574205, 576063),
(753679, 755536), (923145, 925009),
(1512993, 1514851), (1565879, 1567735),
(1990436, 1992292)]},
    "Insertion sequence IS1563": {"color":"pink", "coords":[   (268372, 269862), (506862, 508353),  (570226, 571674), (726296, 727771), (764489, 765949),   (1669437, 1670894), (1823321, 1824817)]}
}

for rname, rdata in regions.items():
    for start, end in rdata["coords"]:
        fig.add_vrect(x0=start, x1=end, fillcolor=rdata["color"], opacity=0.3, line_width=0)

# Add region legend
for rname, rdata in regions.items():
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(size=15, color=rdata["color"]),
        showlegend=True,
        name=rname
    ))

# Significance threshold line
fig.add_hline(y=-np.log10(0.01), line_dash="dash", line_color="grey", annotation_text="p=0.01", annotation_position="top left")

# Layout
fig.update_layout(
    title="Interactive Manhattan plot of GWAS kmers",
    xaxis_title="Genomic position - ref; aquaticclade2_NC_018646",
    yaxis_title="-log10(LRT p-val-FDR)",
    legend_title_text="Strain / Region",
    template="plotly_white"
)

# Show plot
pio.renderers.default = "browser"
fig.show()