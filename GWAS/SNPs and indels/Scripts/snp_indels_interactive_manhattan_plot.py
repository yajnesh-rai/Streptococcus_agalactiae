import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from Bio import SeqIO
import os

# Set data folder path (relative)
data_folder = os.path.join(os.path.dirname(__file__), "../data")


# Load GWAS SNP results
snp_file = os.path.join(data_folder, "pyseer_results_coresnp_fdr_annotation.txt")
snp_gwas = pd.read_csv(snp_file, sep="\t")
snp_gwas['lrt-pval-FDR'] = pd.to_numeric(snp_gwas['lrt-pval-FDR'], errors='coerce')
snp_gwas['position'] = snp_gwas['variant'].str.rsplit('_', n=3, expand=True)[1].astype(float)
snp_gwas['minus_log10_p'] = -np.log10(snp_gwas['lrt-pval-FDR'])
snp_gwas['type'] = 'SNP'

# Load GWAS Indel results
indel_file = os.path.join(data_folder, "pyseer_results_indels_fdr_annotation.txt")
indel_gwas = pd.read_csv(indel_file, sep="\t")
indel_gwas['lrt-pval-FDR'] = pd.to_numeric(indel_gwas['lrt-pval-FDR'], errors='coerce')
indel_gwas['position'] = indel_gwas['variant'].str.rsplit('_', n=3, expand=True)[1].astype(float)
indel_gwas['minus_log10_p'] = -np.log10(indel_gwas['lrt-pval-FDR'])
indel_gwas['type'] = 'Indel'

# Combine SNPs and Indels
gwas = pd.concat([snp_gwas, indel_gwas], ignore_index=True)

# Load GenBank file of reference used for snp/indel calling and map gene/product
gb_file = os.path.join(data_folder, "ref_NZ_CP012419.gb")
genes = []

for record in SeqIO.parse(gb_file, "genbank"):
    for feature in record.features:
        if feature.type == "CDS":
            start = int(feature.location.start)
            end = int(feature.location.end)
            gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
            product = feature.qualifiers.get("product", ["unknown_product"])[0]
            genes.append({
                "chromosome": record.id,
                "start": start,
                "end": end,
                "gene": gene_name,
                "product": product
            })

genes_df = pd.DataFrame(genes)

# Map variants to genes
def map_gene_product(pos):
    hits = genes_df[(genes_df["start"] <= pos) & (genes_df["end"] >= pos)]
    if not hits.empty:
        return hits.iloc[0]["gene"], hits.iloc[0]["product"]
    else:
        return "intergenic", "intergenic_region"

mapped = gwas["position"].apply(map_gene_product)
gwas["gene"] = mapped.apply(lambda x: x[0])
gwas["product"] = mapped.apply(lambda x: x[1])

# Assign colors based on annotation
def color_by_annotation(row):
    ann = str(row.get("annotation", "")).lower()
    t = row["type"]
    if t == "SNP":
        if "stopgained" in ann and "aquatic" in ann:
            return "blue"
        elif "stopgained" in ann and "clinical" in ann:
            return "red"
        elif "nonsynonymous" in ann:
            return "purple"
        elif "synonymous" in ann:
            return "pink"
        else:
            return "grey"
    else:
        if "pseudo" in ann and "aquatic" in ann:
            return "blue"
        elif "pseudo" in ann and "clinical" in ann:
            return "red"
        elif "variant" in ann:
            return "pink"
        else:
            return "grey"

gwas["color"] = gwas.apply(color_by_annotation, axis=1)

# Create interactive Manhattan plot with Plotly
fig = px.scatter(
    gwas,
    x="position",
    y="minus_log10_p",
    color="color",
    symbol="type",
    symbol_map={"SNP": "circle", "Indel": "star"},
    color_discrete_map="identity",
    custom_data=["variant", "annotation", "gene", "product", "minus_log10_p", "type"],
    labels={
        "position": "Genomic position (bp)",
        "minus_log10_p": "-log10(FDR-corrected p-value)"
    },
    title="Interactive Manhattan Plot (SNPs + Indels with Gene/Product Mapping)"
)

# Hover template
fig.update_traces(
    hovertemplate=(
        "<b>Variant:</b> %{customdata[0]}<br>"
        "<b>Annotation:</b> %{customdata[1]}<br>"
        "<b>Gene:</b> %{customdata[2]}<br>"
        "<b>Product:</b> %{customdata[3]}<br>"
        "<b>-log10(FDR p):</b> %{customdata[4]:.3f}<br>"
        "<b>Type:</b> %{customdata[5]}<br>"
        "<b>Position:</b> %{x:.0f} bp"
    )
)

# Significance line
sig_threshold = 0.01
sig_y = -np.log10(sig_threshold)
fig.add_shape(
    type="line",
    x0=gwas['position'].min(),
    x1=gwas['position'].max(),
    y0=sig_y,
    y1=sig_y,
    line=dict(color="black", dash="dash", width=1),
    layer="above"
)
fig.add_annotation(
    x=gwas['position'].max(),
    y=sig_y + 0.5,
    text="Sig_threshold FDR-p-val = 0.01",
    showarrow=False,
    font=dict(size=10, color="black"),
    xanchor="right"
)

# Legends
legend_items = {
    "Stopgained (aquatic)": "blue",
    "Stopgained (clinical)": "red",
    "Nonsynonymous SNP": "purple",
    "Synonymous SNP/ Variant product from indels": "pink",
    "intergenic / non-significant": "grey"
}
for name, col in legend_items.items():
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode="markers",
        marker=dict(size=10, color=col, symbol="circle"),
        name=name,
        showlegend=True
    ))

# Layout settings
fig.update_traces(marker=dict(size=8, line=dict(width=0)))
fig.update_layout(
    template="simple_white",
    width=950,
    height=420,
    font=dict(size=11),
    title_font=dict(size=14),
    margin=dict(l=60, r=40, t=60, b=50),
    legend_title_text="Annotation"
)

# Show interactive plot
pio.renderers.default = "browser"
fig.show()
