"""
maclandrol: 22/07/19
This is an attempt to reverse engineer the BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures)
approach for molecule fragmentation and use it as an heuristic for assembling molecules. The original paper on BRICS can be found here:
http://dacemirror.sci-hub.tw/journal-article/93060992e8d889318b77b562c0e5b75f/degen2008.pdf.
This makes senses from a methodological point of view, but I can't either guarantee that its is working as expected of if it's the best
way to tackle this problem. The goal here is to reconstruct a set of original molecules, which if they were to be fragmented using BRICS, should yield
the same fragment set in input. Thus, in theory fragments obtained using BRICS CAN be assembled into the original molecules with this method.
This differs from rdkit BRICSBuild implementation that requires the presence of dummy indicator atoms added by a prior BRICS fragmentation.
"""

from typing import Optional

import copy
import json
import itertools
import random
import re
import pkg_resources

from functools import lru_cache

import numpy as np

from rdkit import Chem
from rdkit.Chem import rdChemReactions

import datamol as dm
from ..types import Mol

CCQ = "[$([#6][!#6;!#1]):1]!@[#6;!a:2]>>[*:1].[*:2]"
CCQ_RETRO = "[$([#6;!H0][!#6;!#1]):1].[#6;!a;!H0:2]>>[*:1][*:2]"

BRICS_ENVIRONS = {
    "L1": "[C;D3]([#0,#6,#7,#8])(=O)",
    "L3": "[O;D2]-;!@[#0,#6,#1]",
    "L4": "[C;!D1;!$(C=*)]-;!@[#6]",
    # 'L5':'[N;!D1;!$(N*!-*);!$(N=*);!$(N-[!C;!#0])]-[#0,C]',
    "L5": "[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]",
    "L6": "[C;D3;!R](=O)-;!@[#0,#6,#7,#8]",
    "L7": "[C;D2,D3]-[#6]",
    "#L8": "[C;!R;!D1]-;!@[#6]",
    "L8": "[C;!R;!D1;!$(C!-*)]",
    "L9": "[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]",
    "L10": "[N;R;$(N(@C(=O))@[C,N,O,S])]",
    "L11": "[S;D2](-;!@[#0,#6])",
    "L12": "[S;D4]([#6,#0])(=O)(=O)",
    "L13": "[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]",
    "L14": "[c;$(c(:[c,n,o,s]):[n,o,s])]",
    "L15": "[C;$(C(-;@C)-;@C)]",
    "L16": "[c;$(c(:c):c)]",
}

BRICS_RETRO_ENVIRONS = {
    "L1": "[C;!D4;!D3](=O)",
    "L3": "[O;D1]-;!@[#0,#6,#1]",
    "L4": "[C;!D4;!$(C=*)]-;!@[#6]",
    # 'L5':'[N;!D1;!$(N*!-*);!$(N=*);!$(N-[!C;!#0])]-[#0,C]',
    "L5": "[N;!D3;!X4;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]",
    "L6": "[C;!D3;!D4;!R](=O)-;!@[#0,#6,#7,#8]",
    "L7": "[C;D1,D2]-[#6]",
    "#L8": "[C;!R;!D4]-;!@[#6]",
    "L8": "[C;!R;!D4;!$(C!-*)]",
    # I updated this one to require an implicit hydrogen
    "L9": "[n;+0;h;$(n(:[c,n,o,s]):[c,n,o,s])]",
    "L10": "[N;D2;R;$(N(@C(=O))@[C,N,O,S])]",
    "L11": "[S;D1](-;!@[#0,#6])",
    "L12": "[S;D3]([#6,#0])(=O)(=O)",
    "L13": "[C;h;$(C(-;@[C,N,O,S,F])-;@[N,O,S])]",
    "L14": "[c;h;$(c(:[c,n,o,s,F]):[n,o,s])]",
    "L15": "[C;h;$(C(-;@C)-;@C)]",
    "L16": "[c;h;$(c(:c):c)]",
}

REACTIONSDEFS = (
    # L1
    [
        ("1", "3", "-"),
        ("1", "5", "-"),
        ("1", "10", "-"),
    ],
    # L3
    [
        ("3", "4", "-"),
        ("3", "13", "-"),
        ("3", "14", "-"),
        ("3", "15", "-"),
        ("3", "16", "-"),
    ],
    # L4
    [
        ("4", "5", "-"),
        ("4", "11", "-"),
    ],
    # L5
    [
        ("5", "12", "-"),
        ("5", "14", "-"),
        ("5", "16", "-"),
        ("5", "13", "-"),
        ("5", "15", "-"),
    ],
    # L6
    [
        ("6", "13", "-"),
        ("6", "14", "-"),
        ("6", "15", "-"),
        ("6", "16", "-"),
    ],
    # L7
    [
        ("7", "7", "="),
    ],
    # L8
    [
        ("8", "9", "-"),
        ("8", "10", "-"),
        ("8", "13", "-"),
        ("8", "14", "-"),
        ("8", "15", "-"),
        ("8", "16", "-"),
    ],
    # L9
    [
        ("9", "13", "-"),
        ("9", "14", "-"),
        ("9", "15", "-"),
        ("9", "16", "-"),
    ],
    # L10
    [
        ("10", "13", "-"),
        ("10", "14", "-"),
        ("10", "15", "-"),
        ("10", "16", "-"),
    ],
    # L11
    [
        ("11", "13", "-"),
        ("11", "14", "-"),
        ("11", "15", "-"),
        ("11", "16", "-"),
    ],
    # L13
    [
        ("13", "14", "-"),
        ("13", "15", "-"),
        ("13", "16", "-"),
    ],
    # L14
    [
        ("14", "14", "-"),
        ("14", "15", "-"),
        ("14", "16", "-"),
    ],
    # L15
    [
        ("15", "16", "-"),
    ],
    # L16
    [
        ("16", "16", "-"),
    ],
)


@lru_cache(maxsize=None)
def getbrics_list():
    all_brics = []
    all_brics_retro = []
    all_brics_type = []
    smartsGps = copy.deepcopy(REACTIONSDEFS)
    for gp in smartsGps:
        for defn in gp:
            g1, g2, bnd = defn
            r1 = BRICS_ENVIRONS["L" + g1]
            r2 = BRICS_ENVIRONS["L" + g2]
            rr1 = BRICS_RETRO_ENVIRONS["L" + g1]
            rr2 = BRICS_RETRO_ENVIRONS["L" + g2]
            g1 = re.sub("[a-z,A-Z]", "", g1)
            g2 = re.sub("[a-z,A-Z]", "", g2)
            # sma = '[$(%s):1].[$(%s):2]>>[$(%s):1]%s[$(%s):2]' % (r1, r2, r1, bnd, r2)
            # should be equivalent to above ...
            retrosma = "[$(%s):1].[$(%s):2]>>[*:1]%s[*:2]" % (rr1, rr2, bnd)
            sma = "[$(%s):1]%s;!@[$(%s):2]>>[%s*]-[*:1].[%s*]-[*:2]" % (
                r1,
                bnd,
                r2,
                g1,
                g2,
            )
            try:
                rx = rdChemReactions.ReactionFromSmarts(sma)
                retro_rx = rdChemReactions.ReactionFromSmarts(retrosma)
                rx.Initialize()
                retro_rx.Initialize()
                all_brics_type.append("L{}-L{}".format(g1, g2))
                all_brics.append(rx)
                all_brics_retro.append(retro_rx)
            except Exception:
                pass
    return all_brics, all_brics_retro, all_brics_type


@lru_cache(maxsize=None)
def get_reactions_list():
    all_rxns = []
    all_rxns_retro = []
    all_rxns_type = []
    with pkg_resources.resource_stream("datamol", "data/reactions.json") as IN:
        rxns = json.load(IN)
        for k, data in rxns.items():
            try:
                rxn = data.get("syn_smarts")
                retro_rxn = data.get("retro_smarts")
                retro_rxn = rdChemReactions.ReactionFromSmarts(retro_rxn)
                rxn = rdChemReactions.ReactionFromSmarts(rxn)
                if (rxn.GetNumReactantTemplates() == 2 and rxn.GetNumProductTemplates() == 1) and (
                    retro_rxn.GetNumReactantTemplates() == 1
                    and retro_rxn.GetNumProductTemplates() == 2
                ):
                    all_rxns.append(retro_rxn)
                    all_rxns_retro.append(rxn)
                    all_rxns_type.append(k)

            except Exception:
                pass
    return all_rxns, all_rxns_retro, all_rxns_type


ALL_RXNS, ALL_RXNS_RETRO, ALL_RXNS_TYPE = get_reactions_list()
ALL_BRICS, ALL_BRICS_RETRO, ALL_BRICS_TYPE = getbrics_list()


def _can_continue_with(mol, rxns):
    for rxn in rxns:
        if rxn.IsMoleculeReactant(mol):
            return True
    return False


def _run_at_all_rct(rxn, mol1, mol2):
    library = []
    rxn = rdChemReactions.ReactionFromSmarts(rdChemReactions.ReactionToSmarts(rxn))
    # display(rxn)
    m1 = rxn.GetReactantTemplate(0)
    m2 = rxn.GetReactantTemplate(1)
    mol1_valid = mol1 is not None
    mol2_valid = mol2 is not None
    isR1 = mol1_valid and mol1.HasSubstructMatch(m1)
    isR2 = mol1_valid and mol1.HasSubstructMatch(m2)
    if isR1 and mol2_valid and mol2.HasSubstructMatch(m2):
        library.extend(rxn.RunReactants((mol1, mol2)))
    if isR2 and mol2_valid and mol2.HasSubstructMatch(m1):
        library.extend(rxn.RunReactants((mol2, mol1)))
    if library:
        library = list(itertools.chain(*library))
    for m in library:
        mol = None
        mSmi = ""
        try:
            mSmi = dm.to_smiles(m)
            mol = dm.to_mol(mSmi)
        except:
            pass
        if mol is None:
            try:
                mol.UpdatePropertyCache()
                mol = dm.sanitize_mol(mol)
                mSmi = dm.to_smiles(m)
                mol = dm.to_mol(mSmi)
            except:
                pass
        if mSmi:
            yield mol, mSmi


def break_mol(
    mol: Chem.rdchem.Mol,
    minFragmentSize: int = 1,
    silent: bool = True,
    onlyUseReactions: list = [],
    randomize: bool = False,
    mode: str = "brics",
    returnTree: bool = False,
):
    """Breaks a molecules into a list of fragment."""

    if mode.lower() == "brics":
        all_reactions = ALL_BRICS
        all_reactions_type = ALL_BRICS_TYPE
    elif mode.lower() == "rxn":
        all_reactions = ALL_RXNS
        all_reactions_type = ALL_RXNS_TYPE
    else:
        all_reactions = ALL_BRICS + ALL_RXNS
        all_reactions_type = ALL_BRICS_TYPE + ALL_RXNS_TYPE
    if randomize:
        p = np.random.permutation(len(all_reactions))
        all_reactions = [all_reactions[ind] for ind in p]
        all_reactions_type = [all_reactions_type[ind] for ind in p]

    nx = dm.graph._get_networkx()
    mSmi = dm.to_smiles(mol, isomeric=True)
    G = nx.DiGraph()
    node_num = 0
    G.add_node(node_num, smiles=mSmi, mol=mol)
    allNodes = set()
    activePool = {mSmi: node_num}
    allNodes.add(mSmi)
    while activePool:
        nSmi = list(activePool.keys())[0]
        parent = activePool.pop(nSmi)
        node = G.nodes[parent]
        mol = node["mol"]
        for rxnIdx, reaction in zip(all_reactions_type, all_reactions):
            if onlyUseReactions and rxnIdx not in onlyUseReactions:
                continue
            ps = reaction.RunReactants((mol,))
            if ps:

                all_pass = [
                    all([prod.GetNumAtoms(onlyExplicit=True) > minFragmentSize for prod in p_])
                    for p_ in ps
                ]
                nz_i = 0
                while nz_i < len(all_pass) and not all_pass[nz_i]:
                    nz_i += 1
                if not silent:
                    print(nSmi, "->", len(ps), "products and selected ", nz_i)
                    # display(MolsToGridImage(list(itertools.chain(*list(ps))), molsPerRow=2))
                prodSeq = ps[nz_i % len(all_pass)]
                seqOk = True
                # we want to disqualify small fragments, so sort the product sequence by size
                prodSeq = [(prod.GetNumAtoms(onlyExplicit=True), prod) for prod in prodSeq]
                prodSeq.sort(key=lambda x: x[0])
                for _, prod in prodSeq:
                    prod.sanitized = True
                    try:
                        Chem.SanitizeMol(prod)
                    except:
                        if dm.sanitize_mol(prod) is None:
                            seqOk = False
                            break
                        continue
                    pSmi = dm.to_smiles(prod, isomeric=True)
                    seqOk = seqOk and (dm.to_mol(pSmi) is not None)

                    notDummies = sum([atm.GetSymbol() != "*" for atm in prod.GetAtoms()])
                    # nDummies = pSmi.count('*')
                    # if minFragmentSize > 0 and (nats - nDummies < minFragmentSize):
                    if minFragmentSize > 0 and notDummies < minFragmentSize:
                        seqOk = False
                        break
                    prod.pSmi = pSmi

                if seqOk:
                    for _, prod in prodSeq:
                        if not prod.sanitized:
                            continue
                        pSmi = prod.pSmi
                        node_num += 1
                        usmi = dm.to_smiles(dm.fix_mol(prod), isomeric=True)
                        G.add_node(node_num, smiles=usmi, mol=prod)
                        G.add_edge(parent, node_num)
                        if usmi not in allNodes:
                            activePool[pSmi] = node_num
                            allNodes.add(usmi)
                    G.nodes[parent]["rxn"] = rxnIdx
                    break  # at least one reaction matches

    leaves_smiles = [
        G.nodes[n]["smiles"] for n in G.nodes() if G.in_degree(n) != 0 and G.out_degree(n) == 0
    ]
    if returnTree:
        return leaves_smiles, allNodes, G
    return leaves_smiles, allNodes


def build(ll_mols, max_n_mols=float("inf"), mode="brics", frag_rxn=None, ADD_RNXS=[]):
    """Build a super molecule from a list of fragments"""

    seen = set()
    stop = False
    CUR_RXNS = []
    CUR_RXNS_TYPE = []

    if mode == "brics":
        CUR_RXNS = ALL_BRICS_RETRO
        CUR_RXNS_TYPE = ALL_BRICS_TYPE
    elif mode == "rxn":
        CUR_RXNS = ALL_RXNS_RETRO
        CUR_RXNS_TYPE = ALL_RXNS_TYPE
    elif mode is not None:
        CUR_RXNS = ALL_BRICS_RETRO + ALL_RXNS_RETRO
        CUR_RXNS_TYPE = ALL_BRICS_TYPE + ALL_RXNS_TYPE

    if ADD_RNXS is not None:
        ADD_RNXS_TYPE = [f"RXN-{i}" for i in range(len(ADD_RNXS))]
        if isinstance(ADD_RNXS, dict):
            ADD_RNXS_TYPE = ADD_RNXS.keys()
            ADD_RNXS = ADD_RNXS.values()
        CUR_RXNS += list(ADD_RNXS)
        CUR_RXNS_TYPE += list(ADD_RNXS_TYPE)

    for i, rxn_type in enumerate(CUR_RXNS_TYPE):
        if (frag_rxn is not None) and (frag_rxn.strip('"') == rxn_type):
            CUR_RXNS = [CUR_RXNS[i]]
            break

    for fraglist in itertools.product(*ll_mols):
        if stop:
            break

        fraglist = list(fraglist)
        for rxn in CUR_RXNS:  # should be size==1 if frag_rxn is provided
            ps = []
            try:
                ps = _run_at_all_rct(rxn, fraglist[0], fraglist[1])
            except Exception:
                pass
            for m, mSmi in ps:
                if len(seen) >= max_n_mols:
                    stop = True
                    break
                if mSmi not in seen:
                    seen.add(mSmi)
                    yield m


def assemble_fragment_order(
    fragmentlist: list,
    seen: Optional[Mol] = None,
    allow_incomplete: bool = False,
    max_n_mols: float = float("inf"),
    RXNS=None,
):
    """Assemble a list of fragment into a set of possible molecules under rules defined by the brics algorithm

    We are of course assuming:

        1. that the order in the fragmentlist matter :D !
        2. that none of the fragment has explicitly defined hydrogen atoms.
        3. only a list of unique molecule is internally maintained

    Args:
        fragmentlist: list of original fragments to grow
        seen: original molecules used as base. If none, the first element of fragment list will be poped out
        allow_incomplete: Whether to accept assembled molecules with missing fragment
    """

    if RXNS is None:
        RXNS = ALL_BRICS_RETRO

    fragmentlist = list(fragmentlist)
    yield_counter = 0
    if seen is None:
        seen = fragmentlist.pop(0)
    seen = [dm.to_smiles(seen)]  # only one molecule to assemble
    while yield_counter < max_n_mols and len(fragmentlist) > 0:
        # find all the way to add this fragment to seen
        frag = fragmentlist.pop(0)
        level_set = [dm.to_mol(x) for x in seen]
        seen = set()
        for sm in level_set:
            try:
                # there is no point in even trying something on molecules that cannot be kekulized
                for rxn in RXNS:
                    for m, mSmi in _run_at_all_rct(rxn, frag, sm):
                        if allow_incomplete and mSmi not in seen:
                            yield m
                            yield_counter += 1
                        seen.add(mSmi)
            except Exception as e:
                print(e)
                pass

    for m in seen:
        if yield_counter < max_n_mols:
            yield dm.to_mol(m)
            yield_counter += 1


def assemble_fragment_iter(
    fragmentlist,
    seens=None,
    scrambleReagents=False,
    max_n_mols=float("inf"),
    maxdepth=3,
    as_smiles=True,
    RXNS=None,
):
    """Perform an assembly from fragment given all potential RXNS transformation."""

    if RXNS is None:
        RXNS = ALL_BRICS_RETRO

    seen = set()
    if max_n_mols <= 0:
        return
    if not seens:
        seens = list(fragmentlist)
    if scrambleReagents:
        seens = list(seens)
        random.shuffle(seens, random=random.random)

    for seen in seens:
        nextSteps = []
        for rxn in RXNS:
            for fg in fragmentlist:
                for m, pSmi in _run_at_all_rct(rxn, fg, seen):
                    if pSmi not in seen:
                        seen.add(pSmi)
                        yield m if not as_smiles else pSmi
                    if _can_continue_with(m, rxn):
                        nextSteps.append(m)

        if nextSteps and len(seen) <= max_n_mols and maxdepth > 0:
            for p in assemble_fragment_iter(
                fragmentlist,
                seens=nextSteps,
                scrambleReagents=scrambleReagents,
                max_n_mols=(max_n_mols - len(seen)),
                maxdepth=maxdepth - 1,
            ):
                pSmi = dm.to_smiles(p, True)
                if pSmi not in seen:
                    seen.add(pSmi)
                    yield p if not as_smiles else pSmi
                    if len(seen) >= max_n_mols:
                        return
