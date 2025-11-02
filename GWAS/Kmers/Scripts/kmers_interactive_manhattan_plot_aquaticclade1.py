import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import plotly.io as pio


# Load GWAS_kmers CSV file with blast annotations
script_dir = os.path.dirname(os.path.abspath(__file__))  # Folder where script is located
data_folder = script_dir  # Assuming CSV is in the same folder
csv_file = os.path.join(data_folder, "kmer_results_fdr_with_blast_gene_aquaticclade1.csv")

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
    "Genomic islands": {"color":"cyan", "coords":[(305453,380005),(1854245,1868187),(1121926,1137580),(1891278,1901648)]},
    "Prophages": {"color":"yellow", "coords":[(1139258,1191428)]},
    "Insertion sequence IS1381": {"color":"green", "coords":[(198618,199477),(648772,649631),(800965,801824),(1658518,1659377),
                                                (1983501,1984360),(1987484,1988343)]},
    "group-II-introns": {"color":"violet", "coords":[(136093,137949),(190243,192099),(863386,865242),(1032849,1034705),
                                                  (1093854,1095710),(1237365,1239221),(1287769,1289625),(1711098,1712954),
                                                  (1887582,1889438),(1899857,1901713)]},
    "Insertion sequence IS1563": {"color":"pink", "coords":[(32934,34391),(1022436,1023896),(1049748,1051199),(1062176,1063651),
                                              (1241752,1243200),(1304777,1306267),(1494574,1496034),(1545073,1546563),
                                              (2068596,2070092)]}
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
    xaxis_title="Genomic position - ref; aquaticclade1_NZ_CP016501",
    yaxis_title="-log10(LRT p-val-FDR)",
    legend_title_text="Strain / Region",
    template="plotly_white"
)

# Show plot
pio.renderers.default = "browser"
fig.show()