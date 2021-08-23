# pae_to_domains
Graph-based community clustering approach to extract protein domains from a predicted aligned error matrix

## Overview
Using a predicted aligned error matrix corresponding to an AlphaFold2 model (e.g. as downloaded from https://alphafold.ebi.ac.uk/), returns a series of lists of residue indices, where each list corresponds to a set of residues clustering together into a pseudo-rigid domain.
