import networkx as nx
import matplotlib.pyplot as plt

amino_acid_properties = {
    'A': {'hydrophobicity': 0.61, 'alpha_CH_chemical_shifts': 4.35, 'positive_charge': 0, 'negative_charge': 0, 'polarity': -0.06, 'melting_point': 297, 'hydration_potential': 1.94, 'rf_rank': 9.9, 'isoelectric_point': 6.0, 'buriability': 13.4, 'linker_index': 0.0166, 'nneig_index': 50.76, 'sweig_index': -0.414, 'prift_index': -0.96, 'prils_index':-0.26},
    'C': {'hydrophobicity': 1.07, 'alpha_CH_chemical_shifts': 4.65, 'positive_charge': 0, 'negative_charge': 0, 'polarity': 1.36, 'melting_point': 178, 'hydration_potential': -1.24, 'rf_rank': 2.8, 'isoelectric_point': 5.05, 'buriability': 22.6, 'linker_index': 0.5724, 'nneig_index': 58.74, 'sweig_index': 0.162, 'prift_index': 4.54, 'prils_index':0.83},
    'D': {'hydrophobicity': 0.46, 'alpha_CH_chemical_shifts': 4.76, 'positive_charge': 0, 'negative_charge': 1, 'polarity': -0.80, 'melting_point': 270, 'hydration_potential': -10.95, 'rf_rank': 2.8, 'isoelectric_point': 2.77, 'buriability': 8.2, 'linker_index': -0.1278, 'nneig_index': 43.17, 'sweig_index': -1.310, 'prift_index': -5.68, 'prils_index':-1.30},
    'E': {'hydrophobicity': 0.47, 'alpha_CH_chemical_shifts': 4.29, 'positive_charge': 0, 'negative_charge': 1, 'polarity': -0.77, 'melting_point': 249, 'hydration_potential': -10.20, 'rf_rank': 3.2, 'isoelectric_point': 3.22, 'buriability': 7.3, 'linker_index': -0.1794, 'nneig_index': 43.48, 'sweig_index': -1.218, 'prift_index': -3.86, 'prils_index':-0.73},
    'F': {'hydrophobicity': 2.02, 'alpha_CH_chemical_shifts': 4.66, 'positive_charge': 0, 'negative_charge': 1, 'polarity': 1.27, 'melting_point': 270, 'hydration_potential': -10.95, 'rf_rank': 18.8, 'isoelectric_point': 5.48, 'buriability': 23.9, 'linker_index': -0.1278, 'nneig_index': 43.17, 'sweig_index': -1.310, 'prift_index': -5.68, 'prils_index':-1.30},
    'G': {'hydrophobicity': 0.07, 'alpha_CH_chemical_shifts': 3.97, 'positive_charge': 0, 'negative_charge': 0, 'polarity': -0.41, 'melting_point': 290, 'hydration_potential': 2.39, 'rf_rank': 5.6, 'isoelectric_point': 5.97, 'buriability': 7.0, 'linker_index': -0.0442, 'nneig_index': 50.27, 'sweig_index': -0.684, 'prift_index': -1.28, 'prils_index':-0.40},
    'H': {'hydrophobicity': 0.61, 'alpha_CH_chemical_shifts': 4.63, 'positive_charge': 1, 'negative_charge': 0, 'polarity': 0.49, 'melting_point': 277, 'hydration_potential': -10.27, 'rf_rank': 8.2, 'isoelectric_point': 7.59, 'buriability': 11.3, 'linker_index': 0.1643, 'nneig_index': 49.33, 'sweig_index': -0.630, 'prift_index': -0.62, 'prils_index':-0.18},
    'I': {'hydrophobicity': 2.22, 'alpha_CH_chemical_shifts': 3.95, 'positive_charge': 0, 'negative_charge': 0, 'polarity': 1.31, 'melting_point': 284, 'hydration_potential': 2.15, 'rf_rank': 17.1, 'isoelectric_point': 6.02, 'buriability': 20.3, 'linker_index': 0.2758, 'nneig_index': 57.30, 'sweig_index': 1.237, 'prift_index': 5.54, 'prils_index':1.10},
    'K': {'hydrophobicity': 1.15, 'alpha_CH_chemical_shifts': 4.36, 'positive_charge': 1, 'negative_charge': 0, 'polarity': -1.18, 'melting_point': 337, 'hydration_potential': -19.92, 'rf_rank': 4.6, 'isoelectric_point': 10.76, 'buriability': 8.5, 'linker_index': -0.0762, 'nneig_index': 48.66, 'sweig_index': -0.584, 'prift_index': 0.75, 'prils_index':0.08},
    'L': {'hydrophobicity': 1.53, 'alpha_CH_chemical_shifts': 4.17, 'positive_charge': 0, 'negative_charge': 0, 'polarity': 1.21, 'melting_point': 297, 'hydration_potential': 1.94, 'rf_rank': 9.9, 'isoelectric_point': 6.0, 'buriability': 13.4, 'linker_index': 0.0166, 'nneig_index': 50.76, 'sweig_index': -0.414, 'prift_index': -0.96, 'prils_index':-0.26},
    'M': {'hydrophobicity': 1.18, 'alpha_CH_chemical_shifts': 4.52, 'positive_charge': 0, 'negative_charge': 0, 'polarity': 1.27, 'melting_point': 283, 'hydration_potential': -9.68, 'rf_rank': 14.9, 'isoelectric_point': 5.74, 'buriability': 15.7, 'linker_index': -0.0786, 'nneig_index': 45.80, 'sweig_index': -0.916, 'prift_index': -1.94, 'prils_index':-0.46},
    'N': {'hydrophobicity': 0.06, 'alpha_CH_chemical_shifts': 4.75, 'positive_charge': 0, 'negative_charge': 0, 'polarity': -0.48, 'melting_point': 236, 'hydration_potential': -9.68, 'rf_rank': 5.4, 'isoelectric_point': 5.41, 'buriability': 7.6, 'linker_index': -0.0786, 'nneig_index': 45.80, 'sweig_index': -0.916, 'prift_index': -1.94, 'prils_index':-0.46},
    'P': {'hydrophobicity': 0.00, 'alpha_CH_chemical_shifts': 4.44, 'positive_charge': 0, 'negative_charge': 0, 'polarity': 0, 'melting_point': 178, 'hydration_potential': -1.24, 'rf_rank': 2.8, 'isoelectric_point': 5.05, 'buriability': 22.6, 'linker_index': 0.5724, 'nneig_index': 58.74, 'sweig_index': 0.162, 'prift_index': 4.54, 'prils_index':0.83},
    'Q': {'hydrophobicity': 0.00, 'alpha_CH_chemical_shifts': 4.37, 'positive_charge': 0, 'negative_charge': 0, 'polarity': -0.73, 'melting_point': 185, 'hydration_potential': -9.38, 'rf_rank': 9.0, 'isoelectric_point': 5.65, 'buriability': 8.5, 'linker_index': -0.1051, 'nneig_index': 46.09, 'sweig_index': -0.905, 'prift_index': -5.30, 'prils_index':-0.83},
    'R': {'hydrophobicity': 0.60, 'alpha_CH_chemical_shifts': 4.38, 'positive_charge': 1, 'negative_charge': 1, 'polarity': -0.84, 'melting_point': 238, 'hydration_potential': -19.92, 'rf_rank': 4.6, 'isoelectric_point': 10.76, 'buriability': 8.5, 'linker_index': -0.0762, 'nneig_index': 48.66, 'sweig_index': -0.584, 'prift_index': 0.75, 'prils_index':0.08},
    'S': {'hydrophobicity': 0.05, 'alpha_CH_chemical_shifts': 4.50, 'positive_charge': 0, 'negative_charge': 0, 'polarity': -0.5, 'melting_point': 185, 'hydration_potential': -9.38, 'rf_rank': 9.0, 'isoelectric_point': 5.65, 'buriability': 8.2, 'linker_index': -0.1051, 'nneig_index': 46.09, 'sweig_index': -0.905, 'prift_index': -5.30, 'prils_index':-0.83},
    'T': {'hydrophobicity': 0.05, 'alpha_CH_chemical_shifts': 4.35, 'positive_charge': 0, 'negative_charge': 1, 'polarity': -0.27, 'melting_point': 249, 'hydration_potential': -10.20, 'rf_rank': 3.2, 'isoelectric_point': 5.66, 'buriability': 10.3, 'linker_index': -0.1794, 'nneig_index': 43.48, 'sweig_index': -1.218, 'prift_index': -3.86, 'prils_index':-0.73},
    'V': {'hydrophobicity': 1.32, 'alpha_CH_chemical_shifts': 3.95, 'positive_charge': 0, 'negative_charge': 0, 'polarity': 1.09, 'melting_point': 284, 'hydration_potential': 2.15, 'rf_rank': 17.1, 'isoelectric_point': 5.96, 'buriability': 19.5, 'linker_index': 0.2758, 'nneig_index': 57.30, 'sweig_index': 1.237, 'prift_index': 5.54, 'prils_index':1.10},
    'W': {'hydrophobicity': 2.65, 'alpha_CH_chemical_shifts': 4.70, 'positive_charge': 0, 'negative_charge': 0, 'polarity': 0.88, 'melting_point': 290, 'hydration_potential': -5.88, 'rf_rank': 17.1, 'isoelectric_point': 5.89, 'buriability': 24.5, 'linker_index': -0.0442, 'nneig_index': 50.27, 'sweig_index': -0.684, 'prift_index': -1.28, 'prils_index':-0.40},
    'Y': {'hydrophobicity': 1.88, 'alpha_CH_chemical_shifts': 4.60, 'positive_charge': 0, 'negative_charge': 0, 'polarity': 0.33, 'melting_point': 277, 'hydration_potential': -10.27, 'rf_rank': 15.0, 'isoelectric_point': 5.66, 'buriability': 19.5, 'linker_index': 0.1643, 'nneig_index': 49.33, 'sweig_index': -0.630, 'prift_index': -0.62, 'prils_index':-0.18}
}

# Function to create the graph for a given peptide sequence
def create_graph_with_properties(peptide_sequence, amino_acid_groups):
    G = nx.Graph()
    previous_group = None  # Track the last group for connecting group nodes

    for i, aa in enumerate(peptide_sequence):
        aa_group = amino_acid_groups.get(aa, None)
        aa_properties = amino_acid_properties.get(aa, None)

        if aa_group and aa_properties:
            # Create a unique node name for each amino acid occurrence
            amino_acid_node = f"{aa}_{i}" # Ensures that the amino acid has both a group and defined properties

            # Add the amino acid node with its physicochemical properties
            G.add_node(amino_acid_node, label=aa, type='amino_acid', **aa_properties)

            # Add group node if it doesn't exist, initialize with empty properties
            if not G.has_node(aa_group):
                G.add_node(aa_group, label=aa_group, type='group', hydrophobicity=0, alpha_CH_chemical_shifts=0, positive_charge=0, negative_charge=0, polarity=0, melting_point=0, hydration_potential=0, rf_rank=0, isoelectric_point=0, buriability=0, linker_index=0, nneig_index=0, sweig_index=0, prift_index=0, prils_index=0, count=0)
                # The group node is initialized with zero values for all physicochemical properties and a count of 0 (to track the number of amino acids in this group).

            # Add an edge between the amino acid and its group
            G.add_edge(amino_acid_node, aa_group)

            # Update the group node with the amino acid properties
            group_data = G.nodes[aa_group]
            group_data['hydrophobicity'] += aa_properties['hydrophobicity']
            group_data['alpha_CH_chemical_shifts'] += aa_properties['alpha_CH_chemical_shifts']
            group_data['positive_charge'] += aa_properties['positive_charge']
            group_data['negative_charge'] += aa_properties['negative_charge']
            group_data['polarity'] += aa_properties['polarity']
            group_data['melting_point'] += aa_properties['melting_point']
            group_data['hydration_potential'] += aa_properties['hydration_potential']
            group_data['rf_rank'] += aa_properties['rf_rank']
            group_data['isoelectric_point'] += aa_properties['isoelectric_point']
            group_data['buriability'] += aa_properties['buriability']
            group_data['linker_index'] += aa_properties['linker_index']
            group_data['nneig_index'] += aa_properties['nneig_index']
            group_data['sweig_index'] += aa_properties['sweig_index']
            group_data['prift_index'] += aa_properties['prift_index']
            group_data['prils_index'] += aa_properties['prils_index']
            group_data['count'] += 1  # Track the number of amino acids in the group

            # Connect groups if they are adjacent
            if previous_group and previous_group != aa_group:
                G.add_edge(previous_group, aa_group)

            # Update previous group
            previous_group = aa_group

    # Calculate average properties for group nodes
    for node in G.nodes(data=True):
        if node[1]['type'] == 'group' and node[1]['count'] > 0:
            node[1]['hydrophobicity'] /= node[1]['count']
            node[1]['alpha_CH_chemical_shifts'] /= node[1]['count']
            node[1]['positive_charge'] /= node[1]['count']
            node[1]['negative_charge'] /= node[1]['count']
            node[1]['polarity'] /= node[1]['count']
            node[1]['melting_point'] /= node[1]['count']
            node[1]['hydration_potential'] /= node[1]['count']
            node[1]['rf_rank'] /= node[1]['count']
            node[1]['isoelectric_point'] /= node[1]['count']
            node[1]['buriability'] /= node[1]['count']
            node[1]['linker_index'] /= node[1]['count']
            node[1]['nneig_index'] /= node[1]['count']
            node[1]['sweig_index'] /= node[1]['count']
            node[1]['prift_index'] /= node[1]['count']
            node[1]['prils_index'] /= node[1]['count']

    return G

    # node[1] gives you access to the dictionary containing the properties of the node.
    # This allows you to modify or read specific properties of the node, such as hydrophobicity, alpha_CH_chemical_shifts.
    # node[1]['type'] == 'group' checks if the node is a group node (as opposed to an amino acid node).
    # node[1]['count'] > 0 ensures that the group node has amino acids associated with it (the count is greater than zero).

# Function to visualize the graph
def visualize_graph(G, title="Peptide Graph with Properties", save_path=None):
    pos = nx.spring_layout(G)
    node_labels = nx.get_node_attributes(G, 'label')
    
    # Differentiate node colors between amino acids and groups
    node_colors = ['lightgreen' if G.nodes[node].get('type') == 'group' else 'lightblue' for node in G.nodes()]
    
    # Draw the graph
    nx.draw(G, pos, with_labels=True, labels=node_labels, node_color=node_colors, node_size=500, font_size=10, font_color="black")
    plt.title(title)
    plt.savefig(save_path)
    # plt.show()
