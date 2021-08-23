def parse_pae_file(pae_json_file):
    import json, numpy

    with open(pae_json_file, 'rt') as f:
        data = json.load(f)[0]
    
    r1, r2, d = data['residue1'],data['residue2'],data['distance']

    size = max(r1)

    matrix = numpy.empty((size,size))

    matrix.ravel()[:] = d

    return matrix

def domains_from_pae_matrix(pae_matrix, pae_power=1, pae_cutoff=5, graph_resolution=1):
    '''
    Takes a predicted aligned error (PAE) matrix representing the predicted error in distances between each 
    pair of residues in a model, and uses a graph-based community clustering algorithm to partition the model
    into approximately rigid groups.

    Arguments:

        * pae_matrix: a (n_residues x n_residues) numpy array. Diagonal elements should be set to some non-zero
          value to avoid divide-by-zero warnings
        * pae_power (optional, default=1): each edge in the graph will be weighted proportional to (1/pae**pae_power)
        * pae_cutoff (optional, default=5): graph edges will only be created for residue pairs with pae<pae_cutoff
        * graph_resolution (optional, default=1): regulates how aggressively the clustering algorithm is. Smaller values
          lead to larger clusters. Value should be larger than zero, and values larger than 5 are unlikely to be useful.

    Returns: a series of lists, where each list contains the indices of residues belonging to one cluster.
    '''
    import networkx as nx
    import numpy
    weights = 1/pae_matrix**pae_power
    # Within the NetworkX greedy_modularity_communities() method the weight is used to define a term
    # that becomes one element of a dict key. Therefore the use of floating point is inadvisable (can 
    # lead to random `KeyError`s). So, we convert the weights to integers, first multiplying by a 
    # sufficiently large number to make sure everything is left of the decimal place.
    # Note that as of 19 August 2021 `KeyError`s can still happen - this has been reported to the
    # NetworkX developers (https://github.com/networkx/networkx/issues/4992) and possible fixes are being 
    # explored.
    weights = (weights * 1e6).astype(numpy.int)

    g = nx.Graph()
    size = weights.shape[0]
    g.add_nodes_from(range(size))
    for i in range(size):
        for j in range(i+1, size):
            pae = pae_matrix[i,j]
            if pae < pae_cutoff:
                weight = weights[i,j]
                if weight > 0:
                    g.add_edge(i, j, weight=weight)

    from networkx.algorithms import community

    clusters = community.greedy_modularity_communities(g, weight='weight', resolution=graph_resolution)
    return clusters

_defaults = {
    'output_file':  'clusters.csv',
    'pae_power':    1.0,
    'pae_cutoff':   5.0,
    'resolution':   1.0
}

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Extract pseudo-rigid domains from an AlphaFold PAE matrix.')
    parser.add_argument('pae_file', type=str, help="Name of the PAE JSON file.")
    parser.add_argument('--output_file', type=str, default=_defaults['output_file'], help=f'Name of output file (comma-delimited text format. Default: {_defaults["output_file"]}')
    parser.add_argument('--pae_power', type=float, default=_defaults['pae_power'], help=f'Graph edges will be weighted as 1/pae**pae_power. Default: {_defaults["pae_power"]}')
    parser.add_argument('--pae_cutoff', type=float, default=_defaults['pae_cutoff'], help=f'Graph edges will only be created for residue pairs with pae<pae_cutoff. Default: {_defaults["pae_cutoff"]}')
    parser.add_argument('--resolution', type=float, default=_defaults['resolution'], help=f'Higher values lead to stricter (i.e. smaller) clusters. Default: {_defaults["resolution"]}')
    args = parser.parse_args()
    pae = parse_pae_file(args.pae_file)
    clusters = domains_from_pae_matrix(pae, pae_power=args.pae_power, pae_cutoff=args.pae_cutoff, graph_resolution=args.resolution)
    max_len = max([len(c) for c in clusters])
    clusters = [list(c) + ['']*(max_len-len(c)) for c in clusters]
    output_file = args.output_file
    with open(output_file, 'wt') as outfile:
        for c in clusters:
            outfile.write(','.join([str(e) for e in c])+'\n')
    print(f'Wrote {len(clusters)} clusters to {output_file}. Biggest cluster contains {max_len} residues.')
    
