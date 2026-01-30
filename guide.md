# Extending Your Rabies Proteomics Project

A roadmap for advancing your bioinformatics skills with practical exercises and visualization practice.

## Part 1: New Computational Analyses

### 1.1 Secondary Structure Prediction

**Concept**: Predict if amino acids form alpha-helices, beta-sheets, or loops

**Tools**: Use the `PSIPRED` algorithm or simpler `Chou-Fasman` prediction

**Implementation Difficulty**: Intermediate

```python
# Install: pip install propy
from propy import PyPro

def predict_secondary_structure(sequence):
    """Predict secondary structure using Chou-Fasman"""
    # Returns % helix, % sheet, % coil
    pass
```

**What you'll learn**: 
- Dihedral angles and Ramachandran plots
- How local amino acid properties determine structure
- Pattern recognition in sequences

**Visualization opportunity**: Stacked bar charts or circular plots showing structure composition

---

### 1.2 Transmembrane Domain Prediction

**Concept**: Predict which parts of the protein cross through membranes

**Why relevant for rabies**: The G protein has transmembrane domains

**Tools**: `TMHMM`, `Phobius`, or simple hydropathy analysis

```python
def find_transmembrane_regions(sequence, window=9):
    """
    Use Kyte-Doolittle hydropathy scale to identify 
    transmembrane helices (hydrophobic stretches)
    """
    # Returns list of (start, end) coordinates
    pass
```

**What you'll learn**:
- Hydropathy scales and their biochemical basis
- Membrane protein topology
- Signal peptide recognition

**Visualization opportunity**: Linear protein diagrams with colored domains

---

### 1.3 Disorder Prediction

**Concept**: Identify which regions are "floppy" (intrinsically disordered)

**Why it matters**: Disordered regions are often functionally important but not visible in crystal structures

**Tools**: `IUPred`, `PONDR`, or simple approaches using charge/hydropathy ratio

```python
from biopython_extensions import disorder

def predict_disorder(sequence):
    """
    Regions high in charged/polar and low in hydrophobic 
    residues are often disordered
    """
    return disorder.iupred_like_score(sequence)
```

**What you'll learn**:
- Intrinsically disordered proteins (IDPs)
- Post-translational modification sites
- Why wet lab structure determination fails for flexible regions

**Visualization opportunity**: Line plots of disorder probability along sequence

---

### 1.4 Epitope Prediction

**Concept**: Find linear B-cell and T-cell epitopes (vaccine targets)

**Why relevant**: G protein is target for rabies antibodies

**Tools**: `Bepipred`, `BepiPred`, or hydrophilicity-based approaches

```python
def predict_linear_epitopes(sequence):
    """
    Epitopes tend to be on surface (hydrophilic, flexible)
    Use Emini surface probability + flexibility
    """
    pass

def predict_mhc_binding(sequence, mhc_allele='HLA-A*02:01'):
    """
    MHC-peptide binding prediction
    (requires external tool or pre-trained ML model)
    """
    pass
```

**What you'll learn**:
- Immunology basics
- Surface probability calculation
- MHC binding motifs
- Why certain regions are immunogenic

**Visualization opportunity**: Circular heatmaps showing immunogenic regions

---

### 1.5 Signal Peptide and Cleavage Site Detection

**Concept**: Predict where proteins are cleaved and where signal peptides are removed

**Why it matters**: Rabies proteins are post-translationally modified

```python
def predict_signal_peptide(sequence):
    """
    Signal peptides: hydrophobic, at N-terminus, 
    followed by cleavage site
    """
    pass

def predict_cleavage_sites(sequence, protease='furin'):
    """
    Furin cleavage sites: RXXK/R pattern in glycoprotein
    """
    pass
```

**Visualization opportunity**: Schematic diagrams of protein domains and cleavage products

---

### 1.6 Codon Usage and Expression Optimization

**Concept**: Analyze nucleotide sequences (if available) for codon bias

**Why it matters**: Heterologous expression depends on host codon preferences

```python
def codon_usage_analysis(dna_sequence, organism='human'):
    """
    Calculate codon frequency per amino acid
    Compare to organism-specific preferences
    """
    pass

def calculate_cai(dna_sequence, reference_codons):
    """
    Codon Adaptation Index: how well does sequence 
    match host codon preferences?
    Higher CAI = better expression potential
    """
    pass
```

**What you'll learn**:
- The wobble hypothesis
- How organisms have codon preferences
- Synthetic biology optimization

---

### 1.7 Phylogenetic Analysis

**Concept**: Compare your proteins to orthologs from other virus strains and related viruses

**Why it matters**: Understanding evolution and function conservation

```python
def build_phylogenetic_tree(protein_sequences_dict):
    """
    Align sequences
    Calculate evolutionary distances
    Build UPGMA or neighbor-joining tree
    """
    pass
```

**Tools**: `Bio.Phylo`, `ete3`, `EBI's JalView`

**What you'll learn**:
- Evolutionary distances
- Tree building algorithms
- Molecular evolution

**Visualization opportunity**: Dendrogram trees, evolutionary distance heatmaps

---

### 1.8 Domain and Motif Finding

**Concept**: Identify known protein domains (Pfam, InterPro)

**Why it matters**: Domains are functional units; understanding them explains protein function

```python
def search_pfam_domains(sequence, protein_id):
    """
    Query InterProScan or Pfam database
    Returns: domain positions, E-values, descriptions
    """
    pass
```

**What you'll learn**:
- Protein domain databases
- Hidden Markov Models
- Functional annotation

**Visualization opportunity**: Domain architecture diagrams (like on UniProt)

---

### 1.9 Charge Distribution Analysis

**Concept**: Map positive/negative charges along the protein

**Why it matters**: Charge distribution affects protein-protein interactions, RNA binding

```python
def calculate_charge_distribution(sequence, window=7):
    """
    Moving window calculation of net charge
    pH-dependent (pKa of His, Lys, Arg, Asp, Glu)
    """
    return [(i, window_charge) for i, window_charge in enumerate(...)]
```

**Visualization opportunity**: Line plots, heatmaps, 3D visualizations

---

### 1.10 Homology Modeling and Structure Analysis

**Concept**: Build 3D structure predictions

**Advanced concept**: Use AlphaFold2 or Swiss-Model

```python
def predict_structure_alphafold(sequence):
    """
    Use ColabFold API or local AlphaFold2
    Returns PDB file with predicted structure
    """
    pass

def download_pdb_structure(pdb_id):
    """
    Download experimental structures for comparison
    """
    pass
```

**Tools**: `py3Dmol`, `biopython.PDB`, `AlphaFold2`

**What you'll learn**:
- Machine learning for structure prediction
- PDB file format
- Structure quality metrics (pLDDT, PAE)

**Visualization opportunity**: 3D molecular visualizations!

---

## Part 2: Visualization Library Practice

### Level 1: Static Plots (Matplotlib/Seaborn)

#### Exercise 1.1: Multi-panel protein diagrams

```python
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_protein_architecture(protein_name, domains, length):
    """
    Create UniProt-style linear protein diagram
    
    Shows:
    - Protein length as horizontal bar
    - Colored domains at correct positions
    - Key features (signal peptides, cleavage sites, etc.)
    """
    fig, ax = plt.subplots(figsize=(14, 3))
    
    # Protein backbone
    ax.add_patch(patches.Rectangle((0, 0), length, 1, 
                                   facecolor='lightgray'))
    
    # Add domains with different colors
    colors = ['red', 'blue', 'green', 'orange', 'purple']
    for i, (name, start, end) in enumerate(domains):
        ax.add_patch(patches.Rectangle((start, 0), end-start, 1,
                                       facecolor=colors[i % len(colors)],
                                       alpha=0.7,
                                       edgecolor='black'))
        ax.text((start+end)/2, 0.5, name, ha='center', va='center',
               fontsize=9, fontweight='bold')
    
    ax.set_xlim(0, length)
    ax.set_ylim(-0.5, 1.5)
    ax.set_aspect('equal')
    ax.set_title(f'{protein_name} Architecture')
    plt.tight_layout()
    return fig
```

**Skills**: Patch objects, positioning, layering

---

#### Exercise 1.2: Heatmap variations

```python
def plot_advanced_heatmaps(properties_df):
    """Multiple heatmap styles practice"""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Standard heatmap
    sns.heatmap(properties_df.corr(), ax=axes[0,0], 
                cmap='coolwarm', center=0, annot=True)
    
    # 2. Clustered heatmap (requires scipy)
    from scipy.cluster import hierarchy
    from scipy.spatial.distance import pdist
    g = sns.clustermap(properties_df.corr(), cmap='viridis')
    
    # 3. Annotated heatmap with text colors
    sns.heatmap(properties_df.set_index('Protein ID'),
                ax=axes[1,0], cmap='RdYlGn', 
                annot=True, fmt='.2f',
                cbar_kws={'label': 'Value'})
    
    # 4. Diverging heatmap for fold changes
    # (simulated comparison data)
    fold_change = (properties_df.set_index('Protein ID').T / 
                   properties_df.set_index('Protein ID').T.mean())
    sns.heatmap(fold_change, ax=axes[1,1], 
                cmap='PiYG', center=1, vmin=0.5, vmax=1.5)
    
    plt.tight_layout()
```

**Skills**: Color mapping, clustering, annotations, diverging colormaps

---

#### Exercise 1.3: Multi-scale plots

```python
def plot_residue_properties_along_sequence(sequence, properties_array):
    """
    Plots properties calculated per residue position
    Good for: hydropathy, disorder, epitope scores, etc.
    """
    fig, axes = plt.subplots(4, 1, figsize=(16, 8), 
                             sharex=True)
    
    # Property 1: Smooth line plot
    axes[0].plot(properties_array, linewidth=2, color='steelblue')
    axes[0].fill_between(range(len(properties_array)), 
                         properties_array, alpha=0.3)
    axes[0].set_ylabel('Property Value')
    
    # Property 2: Bar plot
    colors = ['red' if x > 0.5 else 'blue' 
              for x in properties_array]
    axes[1].bar(range(len(properties_array)), 
                properties_array, color=colors)
    axes[1].set_ylabel('Score')
    
    # Property 3: Scatter with size/color
    axes[2].scatter(range(len(properties_array)), 
                    properties_array,
                    s=properties_array*100,  # size varies
                    c=properties_array,       # color varies
                    cmap='viridis', alpha=0.6)
    axes[2].set_ylabel('Value')
    
    # Property 4: Position-specific annotation
    # Highlight regions of interest
    axes[3].imshow(np.array([properties_array]), 
                   aspect='auto', cmap='RdYlGn')
    axes[3].set_ylabel('Protein')
    axes[3].set_xlabel('Residue Position')
    
    plt.tight_layout()
```

**Skills**: Subplots, fill_between, scatter styling, imshow

---

### Level 2: Interactive Plots (Plotly)

#### Exercise 2.1: Interactive scatter and 3D plots

```python
import plotly.express as px
import plotly.graph_objects as go

def interactive_protein_comparison(properties_df):
    """
    3D scatter: MW vs pI vs Instability
    Hover for details, click legend to toggle
    """
    fig = px.scatter_3d(
        properties_df,
        x='Molecular Weight (kDa)',
        y='Isoelectric Point',
        z='Instability Index',
        color='Protein ID',
        size='Length (aa)',
        hover_data=['Protein Name', 'GRAVY'],
        title='Rabies Proteins in 3D Property Space',
        labels={'Molecular Weight (kDa)': 'MW (kDa)',
                'Isoelectric Point': 'pI',
                'Instability Index': 'II'},
    )
    fig.update_traces(marker=dict(size=8))
    return fig

def interactive_sequence_heatmap(matrix_data):
    """
    Hoverable heatmap: hover to see exact values
    Click to zoom, double-click to reset
    """
    fig = go.Figure(data=go.Heatmap(
        z=matrix_data,
        hoverongaps=False,
        hovertemplate='Position: %{x}<br>Position: %{y}<br>Identity: %{z:.3f}<extra></extra>',
        colorscale='Viridis'
    ))
    fig.update_layout(title='Pairwise Sequence Identity (Hoverable)')
    return fig
```

**Skills**: 3D plots, hover templates, interactivity, sizing

---

#### Exercise 2.2: Subplots and synchronized charts

```python
from plotly.subplots import make_subplots

def dashboard_multiple_perspectives(properties_df, aa_composition_df):
    """Create a dashboard with multiple synchronized plots"""
    
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'scatter'}, {'type': 'bar'}],
               [{'type': 'scatter'}, {'type': 'box'}]],
        subplot_titles=('MW vs pI', 'Protein Length',
                       'Aromaticity vs GRAVY', 'II Distribution')
    )
    
    # Scatter plot
    fig.add_trace(
        go.Scatter(x=properties_df['Molecular Weight (kDa)'],
                   y=properties_df['Isoelectric Point'],
                   mode='markers+text',
                   text=properties_df['Protein ID'],
                   textposition='top center',
                   marker=dict(size=10)),
        row=1, col=1
    )
    
    # Bar chart
    fig.add_trace(
        go.Bar(x=properties_df['Protein ID'],
               y=properties_df['Length (aa)'],
               marker_color='indianred'),
        row=1, col=2
    )
    
    # Another scatter
    fig.add_trace(
        go.Scatter(x=properties_df['Aromaticity'],
                   y=properties_df['GRAVY'],
                   mode='markers',
                   marker=dict(size=12, color=properties_df['Instability Index'],
                              colorscale='Viridis', showscale=True)),
        row=2, col=1
    )
    
    # Box plots
    for protein in properties_df['Protein ID']:
        aa_data = aa_composition_df[aa_composition_df['Protein'] == protein].iloc[0]
        fig.add_trace(
            go.Box(y=list(aa_data[1:]), name=protein),
            row=2, col=2
        )
    
    fig.update_xaxes(title_text='MW (kDa)', row=1, col=1)
    fig.update_yaxes(title_text='pI', row=1, col=1)
    
    return fig
```

**Skills**: Subplots, multiple trace types, color scales, text annotations

---

#### Exercise 2.3: Animated plots

```python
def animated_sequence_alignment_progress(alignment_file):
    """
    Animate the alignment process step by step
    """
    import plotly.graph_objects as go
    
    # Load alignment and create frames
    frames = []
    for i in range(10, len(alignment), 10):
        # Each frame shows alignment up to position i
        frame_data = [go.Scatter(y=distances_at_position_i)]
        frames.append(go.Frame(data=frame_data, name=str(i)))
    
    fig = go.Figure(
        data=[frames[0].data[0]],
        frames=frames
    )
    
    fig.update_layout(
        updatemenus=[dict(type='buttons',
                        buttons=[
                            dict(label='Play', method='animate',
                                 args=[None, {'frame': {'duration': 500}}]),
                            dict(label='Pause', method='animate',
                                 args=[[None], {'frame': {'duration': 0}}])
                        ])]
    )
    return fig
```

**Skills**: Animation, frames, playback controls

---

### Level 3: Advanced Visualizations

#### Exercise 3.1: Network diagrams

```python
import plotly.graph_objects as go
import networkx as nx

def protein_interaction_network():
    """
    Show which proteins are similar to each other
    Node size = protein length
    Edge weight = sequence identity
    """
    
    # Create network
    G = nx.Graph()
    
    # Add nodes (proteins)
    for idx, row in properties_df.iterrows():
        G.add_node(row['Protein ID'], 
                  size=row['Length (aa)'] / 10)
    
    # Add edges (sequence similarity)
    for i, prot1 in enumerate(proteins):
        for j, prot2 in enumerate(proteins):
            if i < j:
                similarity = identity_matrix.loc[prot1, prot2]
                if similarity > 0.2:  # Only if similar enough
                    G.add_edge(prot1, prot2, weight=similarity)
    
    # Layout
    pos = nx.spring_layout(G, k=2, iterations=50)
    
    # Create Plotly figure
    edge_x, edge_y = [], []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
    
    edge_trace = go.Scatter(x=edge_x, y=edge_y,
                           line=dict(width=0.5, color='gray'),
                           mode='lines')
    
    node_x, node_y = [], []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
    
    node_trace = go.Scatter(x=node_x, y=node_y,
                           mode='markers+text',
                           text=list(G.nodes()),
                           textposition='top center',
                           marker=dict(size=20, color='lightblue'))
    
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(title='Rabies Protein Similarity Network',
                     showlegend=False)
    return fig
```

**Skills**: Network graphs, spring layouts, network analysis

---

#### Exercise 3.2: Parallel coordinates

```python
def parallel_coordinates_proteins(properties_df):
    """
    Compare multiple properties simultaneously
    Hover to see individual protein paths
    """
    fig = go.Figure(data=
        go.Parcoords(
            line=dict(color=properties_df['Instability Index'],
                     colorscale='Viridis',
                     showscale=True),
            dimensions=[
                dict(label='MW (kDa)',
                     values=properties_df['Molecular Weight (kDa)']),
                dict(label='pI',
                    values=properties_df['Isoelectric Point']),
                dict(label='Aromaticity',
                    values=properties_df['Aromaticity']),
                dict(label='Length',
                    values=properties_df['Length (aa)']),
            ]
        )
    )
    fig.update_layout(title='Protein Properties: Parallel Coordinates')
    return fig
```

**Skills**: Multidimensional data visualization, brushing

---

#### Exercise 3.3: Sankey diagrams

```python
def protein_classification_sankey():
    """
    Show flow from proteins → structural class → function
    """
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color='black', width=0.5),
            label=['NP', 'P', 'M', 'G', 'L',  # Source
                   'Structural', 'Enzymatic',  # Class
                   'RNA binding', 'Catalysis', 'Membrane', 'Assembly'],  # Function
        ),
        link=dict(
            source=[0, 1, 1, 2, 3, 4,  # indices of source
                    5, 5, 6, 6, 6, 5],
            target=[5, 5, 6, 5, 5, 6,  # indices of target
                    8, 9, 10, 8, 11, 9],
            value=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        )
    )])
    return fig
```

**Skills**: Flow diagrams, categorical relationships

---

#### Exercise 3.4: Sunburst and treemap

```python
def protein_taxonomy_sunburst(data):
    """
    Hierarchical sunburst: Organism → Strain → Protein → Feature
    Click inner ring to zoom in
    """
    fig = go.Figure(go.Sunburst(
        labels=['Rabies', 'Strain A', 'Strain B', 'NP', 'P', 'G',
                'Feature1', 'Feature2'],
        parents=['', 'Rabies', 'Rabies', 'Strain A', 'Strain A', 'Strain B',
                 'NP', 'NP'],
        values=[100, 50, 50, 20, 30, 50, 10, 10],
        marker=dict(colorscale='RdYlBu')
    ))
    return fig

def protein_features_treemap(properties_df):
    """
    Treemap: size = MW, color = stability
    Good for hierarchical part-to-whole relationships
    """
    fig = go.Figure(go.Treemap(
        labels=properties_df['Protein ID'],
        parents=[''] * len(properties_df),
        values=properties_df['Molecular Weight (kDa)'],
        marker=dict(colors=properties_df['Instability Index'],
                   colorscale='RdYlGn_r',
                   cmid=30),
        text=properties_df['Protein Name'],
        textposition='middle center'
    ))
    return fig
```

**Skills**: Hierarchical visualization, interactive drilling

---

### Level 4: Scientific Visualization Libraries

#### Exercise 4.1: PyMOL (3D molecular structures)

```python
import pymol
from pymol import cmd

def visualize_protein_structure(pdb_file):
    """
    Load and visualize protein structure
    Color by property, create beautiful renders
    """
    pymol.finish_launching()
    
    cmd.load(pdb_file, 'protein')
    
    # Color by B-factor (flexibility)
    cmd.spectrum('b', 'blue_white_red', 'protein', 
                minimum=20, maximum=60)
    
    # Show secondary structure
    cmd.dss()
    
    # Show cartoon representation
    cmd.hide('everything')
    cmd.show('cartoon', 'protein')
    cmd.show('sticks', 'protein and ss s')  # Beta sheets in sticks
    
    # Zoom and center
    cmd.center()
    cmd.zoom()
    
    # Save high-quality image
    cmd.ray(2400, 2400)
    cmd.png('protein_structure.png', dpi=300)
```

**Skills**: 3D visualization, structural biology

---

#### Exercise 4.2: py3Dmol (Interactive 3D in Jupyter)

```python
import py3Dmol

def interactive_3d_structure(pdb_data):
    """
    Embed interactive 3D structure viewer in Jupyter notebook
    """
    viewer = py3Dmol.view(query='pdbid:1RZF')  # Rabies virus G protein
    viewer.setStyle({}, {'cartoon': {'color': 'spectrum'}})
    viewer.zoomTo()
    return viewer._make_html()
```

---

#### Exercise 4.3: Plotly Mesh3D (Custom 3D visualization)

```python
def plot_3d_charge_distribution(sequence, window=5):
    """
    3D visualization of charge distribution
    X,Y = 2D sequence layout, Z = charge
    """
    
    # Calculate local charges
    charges = calculate_charge_distribution(sequence, window)
    
    # Arrange in 2D lattice (like a helix unrolled)
    angles = np.arange(len(charges)) * 2 * np.pi / 3.6  # Helix pitch
    x = np.cos(angles)
    y = np.sin(angles)
    z = charges
    
    fig = go.Figure(data=[go.Mesh3d(
        x=x, y=y, z=z,
        opacity=0.7,
        color=z,
        colorscale='RdBu',
        showscale=True
    )])
    
    fig.update_layout(
        title='Charge Distribution Along Protein',
        scene=dict(xaxis_title='Angle', yaxis_title='Angle', zaxis_title='Charge')
    )
    return fig
```

---

### Level 5: Create an Interactive Dashboard

#### Exercise 5.1: Streamlit dashboard

```python
# rabies_dashboard.py
import streamlit as st
import plotly.express as px
from rabies_proteomics_project import RabiesVirusAnalyzer

st.set_page_config(layout='wide', page_title='Rabies Proteomics Explorer')

st.title('🦇 Rabies Virus Proteomics Explorer')

# Sidebar controls
st.sidebar.header('Analysis Options')
selected_proteins = st.sidebar.multiselect(
    'Select proteins to analyze:',
    ['NP', 'P', 'M', 'G', 'L'],
    default=['NP', 'P', 'M', 'G', 'L']
)

# Load data
@st.cache_data
def load_analysis():
    analyzer = RabiesVirusAnalyzer()
    analyzer.create_example_sequences()
    return analyzer.analyze_biophysical_properties()

properties_df = load_analysis()
filtered_df = properties_df[properties_df['Protein ID'].isin(selected_proteins)]

# Main dashboard
col1, col2 = st.columns(2)

with col1:
    st.subheader('Molecular Weight Comparison')
    fig = px.bar(filtered_df, x='Protein ID', y='Molecular Weight (kDa)',
                 color='Protein ID')
    st.plotly_chart(fig, use_container_width=True)

with col2:
    st.subheader('Stability Analysis')
    fig = px.scatter(filtered_df, x='Protein ID', y='Instability Index',
                     size='Length (aa)', color='Instability Index',
                     color_continuous_scale='RdYlGn_r')
    st.plotly_chart(fig, use_container_width=True)

# Detailed table
st.subheader('Detailed Properties')
st.dataframe(filtered_df, use_container_width=True)

# Download results
csv = filtered_df.to_csv(index=False)
st.download_button(
    label='Download data as CSV',
    data=csv,
    file_name='rabies_proteins.csv',
    mime='text/csv'
)
```

**Run with**: `streamlit run rabies_dashboard.py`

---

#### Exercise 5.2: Dash app (more control than Streamlit)

```python
# app.py
import dash
from dash import dcc, html, Input, Output
import plotly.express as px

app = dash.Dash(__name__)

# Load data
properties_df = load_analysis()

app.layout = html.Div([
    html.H1('Rabies Virus Proteomics Dashboard'),
    
    html.Div([
        dcc.Dropdown(
            id='protein-selector',
            options=[{'label': p, 'value': p} 
                    for p in properties_df['Protein ID']],
            value='G',
            style={'width': '100%'}
        ),
    ], style={'width': '48%', 'display': 'inline-block'}),
    
    dcc.Graph(id='property-graph'),
    dcc.Graph(id='composition-graph'),
])

@app.callback(
    Output('property-graph', 'figure'),
    Input('protein-selector', 'value')
)
def update_property_graph(selected_protein):
    filtered = properties_df[properties_df['Protein ID'] == selected_protein]
    return px.bar(filtered, x=['MW', 'pI', 'Aromaticity'],
                 title=f'{selected_protein} Properties')

if __name__ == '__main__':
    app.run_server(debug=True)
```

---

## Part 3: Challenge Projects

### Challenge 1: Comparative Proteomics
Fetch and analyze proteins from multiple rabies strains
- Compare properties across strains
- Identify conserved vs. variable regions
- Create comparison heatmaps

### Challenge 2: Structure-Function Analysis
1. Fetch PDB structures (if available)
2. Load structures with BioPython.PDB or py3Dmol
3. Color by property (B-factor, surface accessibility)
4. Overlay sequence-based predictions on 3D structure

### Challenge 3: Vaccine Design Workflow
1. Predict epitopes across G protein
2. Check MHC binding predictions
3. Visualize all potential epitopes
4. Create "epitope map" showing immunogenic regions

### Challenge 4: Protein Expression Prediction
1. Add codon usage analysis
2. Predict secondary structure
3. Predict disorder regions
4. Create "expressibility score" combining factors

### Challenge 5: Interactive Scientific Paper
Create a self-contained HTML report with:
- Interactive figures
- Embedded text explanations
- Ability to load user's own sequences
- Export results

---

## Recommended Learning Path

### Week 1: Computational Extensions
- [ ] Implement secondary structure prediction (1.1)
- [ ] Add transmembrane prediction (1.2)
- [ ] Learn disorder scoring (1.3)

### Week 2-3: Matplotlib/Seaborn Mastery
- [ ] Exercise 1.1 (protein diagrams)
- [ ] Exercise 1.2 (heatmap variations)
- [ ] Exercise 1.3 (multi-scale plots)

### Week 3-4: Plotly Skills
- [ ] Exercise 2.1 (3D scatter, interactive heatmap)
- [ ] Exercise 2.2 (dashboards)
- [ ] Exercise 2.3 (animations)

### Week 4-5: Advanced Visualization
- [ ] Exercise 3.1-3.4 (networks, parallel coordinates, etc.)
- [ ] Try py3Dmol (4.2)

### Week 5-6: Build Something
- [ ] Create Streamlit dashboard (5.1)
- [ ] Complete Challenge 1 or 2
- [ ] Deploy online (Streamlit Cloud is free)

---

## Resources

### Visualization Libraries Documentation
- **Matplotlib**: https://matplotlib.org/stable/tutorials/index
- **Seaborn**: https://seaborn.pydata.org/tutorial.html
- **Plotly**: https://plotly.com/python/
- **Streamlit**: https://docs.streamlit.io/
- **py3Dmol**: https://3dmol.csb.pitt.edu/

### Bioinformatics Tutorials
- **BioPython**: https://biopython.org/wiki/Documentation
- **Protein Structure**: https://www.coursera.org/learn/protein-structure
- **Epitope Prediction**: https://www.iedb.org/

### Datasets for Practice
- **UniProt**: https://www.uniprot.org/ (protein sequences & annotations)
- **PDB**: https://www.rcsb.org/ (3D structures)
- **NCBI**: https://www.ncbi.nlm.nih.gov/ (all sequences)
- **Virus database**: https://ictv.global/ (viral proteins)

### Communities & Help
- **BioPython Mailing List**: https://biopython.org/wiki/Mailing_lists
- **Stack Overflow**: Tag with `biopython`, `plotly`, `streamlit`
- **Reddit**: r/bioinformatics, r/dataviz

---

## Tips for Success

1. **Start small**: Pick one computation and visualize it well before moving to next
2. **Use real data**: Work with actual rabies sequences early
3. **Iterate on visualization**: Spend time making plots publication-quality
4. **Build incrementally**: Add one feature at a time
5. **Document everything**: Comment your code; future you will thank you
6. **Share your work**: Blog about what you learn; it reinforces understanding
7. **Look at published papers**: See how real scientists visualize their data
8. **Combine skills**: Use multiple libraries in one project (Plotly + Streamlit, etc.)

Good luck! This is an exciting path in bioinformatics. 🧬📊
