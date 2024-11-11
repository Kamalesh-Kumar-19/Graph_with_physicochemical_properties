import torch
from torch_geometric.data import Data
import networkx as nx

def convert_to_pytorch_geometric(G):
    # Create a mapping from node names (strings) to integer indices
    node_mapping = {node: i for i, node in enumerate(G.nodes())}
    
    # Map the edges to integer indices
    edge_index = torch.tensor([[node_mapping[u], node_mapping[v]] for u, v in G.edges()], dtype=torch.long).t().contiguous()

    # Extract node features: physicochemical properties (15 features now)
    node_features = []
    for node, data in G.nodes(data=True):
        if 'hydrophobicity' in data:
            # Use the physicochemical properties as node features
            node_features.append([
                data['hydrophobicity'],       
                data['alpha_CH_chemical_shifts'],   
                data['positive_charge'],   
                data['negative_charge'],        
                data['polarity'],
                data['melting_point'],
                data['hydration_potential'],
                data['rf_rank'],
                data['isoelectric_point'],
                data['buriability'],
                data['linker_index'],
                data['nneig_index'],
                data['sweig_index'],
                data['prift_index'],
                data['prils_index']         
            ])
        else:
            # If node doesn't have properties (for group nodes), use a default zero vector for all 15 features
            node_features.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    # Convert node features to torch tensor
    node_features = torch.tensor(node_features, dtype=torch.float)

    # Create PyTorch Geometric Data object
    data = Data(x=node_features, edge_index=edge_index)

    return data
