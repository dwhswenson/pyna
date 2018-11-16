import itertools
# hydrogen bond motif has lists h-bonds with donor, then acceptor
# the possible hydrogens are found on-the-fly
wc_motif = {
    frozenset({'DG', 'DC'}): [('DG-N1', 'DC-N3'),
                              ('DG-N2', 'DC-O2'),
                              ('DC-N4', 'DG-O6')],
    frozenset({'DA', 'DT'}): [('DT-N3', 'DA-N1'),
                              ('DA-N6', 'DT-O4')]
}

hg_motif = {
    frozenset({'DA', 'DT'}): [()]
}

base_pairing_partner = {'DG': 'DC',
                        'DC': 'DG',
                        'DA': 'DT',
                        'DT': 'DA'}

def get_atom_from_string(string, name_to_residue):
    (res_str, atom_str) = string.split('-')
    residue = name_to_residue[res_str]
    return residue.atom(atom_str)


def hydrogens_attached_to(atom, topology):
    """Find all hydrogens bonded to a given atom.
    """
    bonded_atoms = set(sum([[a1, a2] for a1, a2 in topology.bonds
                            if atom in [a1, a2]], [])) - set({atom})
    bonded_hydrogens = [a for a in bonded_atoms if a.element.symbol == 'H']
    return bonded_hydrogens


def list_possible_triples(residue_1, residue_2, motif, topology):
    """Get H-bond triples for a pair of residues, given H-bonding motif.
    """
    name_to_res = {res.name : res for res in [residue_1, residue_2]}
    hbonds = motif[frozenset(name_to_res.keys())]
    triples = {}
    for (bond_num, (donor_str, acceptor_str)) in enumerate(hbonds):
        donor_atom, acceptor_atom = [
            get_atom_from_string(s, name_to_res)
            for s in [donor_str, acceptor_str]
        ]
        hydrogen_atoms = hydrogens_attached_to(donor_atom, topology)
        local_triples = list(itertools.product([donor_atom], hydrogen_atoms,
                                               [acceptor_atom]))
        triples[bond_num] = local_triples
    return triples


def possible_motif_hydrogen_bonds(topology, motif, chainid, other_chains):
    """
    Returns
    -------
    dict :
        {residue_pair: {motif_bond_number: list_of_tripels}}
    """
    try:
        _ = len(other_chains)
    except TypeError:
        other_chains = [other_chains]

    possible_hbonds = {}
    for residue in topology.chain(chainid).residues:
        eligible_partners = [
            res for idx in other_chains
            for res in topology.chain(idx).residues
            if res.name == base_pairing_partner[residue.name]
        ]
        residue_possible_hbonds = {
            frozenset({residue, partner}): \
                list_possible_triples(residue, partner, motif, topology)
            for partner in eligible_partners
        }
        possible_hbonds.update(residue_possible_hbonds)
    return possible_hbonds


def possible_hbonds_to_indices(possible_hbonds):
    possible_bonds_idx = {}
    for residue_pair, motif_bonds in possible_hbonds.items():
        res_pair_idx = frozenset([r.index for r in residue_pair])
        motif_bonds_idx = {num: [tuple([a.index for a in triple])
                                 for triple in triples]
                           for num, triples in motif_bonds.items()}
        possible_bonds_idx[res_pair_idx] = motif_bonds_idx
    return possible_bonds_idx

def unnest_triples(nested_triples):
    return sum([triples for motif_bonds in nested_triples.values()
                for triples in motif_bonds.values()], [])

def baker_hubbard_from_triples(trajectory, triples, d_cut=0.25,
                               theta_cut=120.0):
    dist_pairs = [(t[1], t[2]) for t in triples]
    distances = md.compute_distances(trajectory, dist_pairs)
    angles = md.compute_angles(trajectory, triples)
    theta_cut *= np.pi / 180.0
    result = np.logical_and(distances < d_cut, angles > theta_cut)
    return result

def check_motif(possible_hbonds_idx, results):
    has_hbond = []
    for res_pair_idx, motif_hbonds in possible_hbonds_idx.items():
        intermediate = {hbond_num: any([results[t] for t in triples])
                        for hbond_num, triples in motif_hbonds.items()}
        if all(intermediate.values()):
            has_hbond += [res_pair_idx]
    return has_hbond
