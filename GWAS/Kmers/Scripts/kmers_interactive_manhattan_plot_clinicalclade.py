import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import plotly.io as pio


# Load GWAS_kmers CSV file with blast annotations
script_dir = os.path.dirname(os.path.abspath(__file__))  # Folder where script is located
data_folder = script_dir  # Assuming CSV is in the same folder
csv_file = os.path.join(data_folder, "kmer_results_fdr_with_blast_gene_clinicalclade.csv")

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
    "Genomic islands": {"color":"cyan", "coords":[(2071253, 2086164),(2042672, 2050461),(1292358, 1305209),(704649, 720913)]},
    "Prophages": {"color":"yellow", "coords":[(544482, 586079),(654859, 707017)]},
    "Insertion sequence IS1381": {"color":"green", "coords":[(152078, 152937),
    (484003, 484862), (1030995, 1031854),
(1184872, 1185731), (1959307, 1960166),
(1995017, 1995876)]},
    "Insertion sequence IS1563": {"color":"pink", "coords":[(265841, 267331), (507187, 508677), (604247, 605695), (892696, 894147), (1875572, 1877041)]}
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
    xaxis_title="Genomic position - ref; clinicalclade_NC_007432",
    yaxis_title="-log10(LRT p-val-FDR)",
    legend_title_text="Strain / Region",
    template="plotly_white"
)

# Show plot
pio.renderers.default = "browser"
fig.show()