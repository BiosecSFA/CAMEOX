import numpy as np
import torch
import torch.nn as nn
import importlib
from pathlib import Path

potts_is_installed = importlib.util.find_spec("potts") is not None
if not (potts_is_installed):
    error_msg = '''Potts library is not installed. 
             Run: `git clone git@github.com:hnisonoff/potts.git`
                  `cd potts`
                  `pip install .`
    '''
    raise ImportError(error_msg)
else:
    from potts import Potts
    from potts.mcmc import GWGCategoricalSampler, OneHotCategorical

ALPHABET = 'ARNDCQEGHILKMFPSTWYV-'
AA_TO_I = {aa: i for i, aa in enumerate(ALPHABET)}
I_TO_AA = {i: aa for aa, i in AA_TO_I.items()}

BASE_TO_I = {b: i for i, b in enumerate("ACGT")}
I_TO_BASE = {i: b for b, i in BASE_TO_I.items()}


def _get_codon_table():
    bases = [l.upper() for l in 'tcag']
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    return codon_table


CODON_TABLE = _get_codon_table()


def load_codon_to_aa():
    hidden = 100
    codon_to_aa_model = nn.Sequential(nn.Linear(12, hidden), nn.ReLU(),
                                      nn.Linear(hidden, hidden), nn.ReLU(),
                                      nn.Linear(hidden, 21))
    model_weights_fn = Path(
        __file__).parent.parent.parent / 'models' / 'codon_to_aa_model.pt'

    codon_to_aa_model.load_state_dict(torch.load(model_weights_fn))
    return codon_to_aa_model


def load_ccmpred(ccmpred_fn):
    '''
    Create a Potts model from ccmpred outptut file
    '''
    lines = [l.strip() for l in open(ccmpred_fn)]
    idx_first_comment = next(
        (i for i, l in enumerate(lines) if l.startswith("#")))
    L = idx_first_comment

    rows = []
    for l in lines[:L]:
        # add 0 for gap?
        row = np.asarray(list(map(float, l.split())) + [0])
        rows.append(row)
    h = np.asarray(rows)

    idxs_comments = [i for i, l in enumerate(lines) if l.startswith("#")]

    # get residue-residue interaction matrix
    W = np.zeros((L, L, 21, 21))
    num_pairs = (L * (L - 1) / 2)
    counter = 0
    for idx in range(len(idxs_comments)):
        i, j = map(int, lines[idxs_comments[idx]].split()[1:])
        start = idxs_comments[idx] + 1
        end = idxs_comments[idx +
                            1] if idx + 1 < len(idxs_comments) else len(lines)
        block = lines[start:end]
        sub_W = np.asarray(
            [np.asarray(list(map(float, l.split()))) for l in block])
        W[i, j] = sub_W
        W[j, i] = sub_W.T
        counter += 1
        if counter >= num_pairs:
            break

    h = torch.tensor(-h, dtype=torch.float64)
    W = torch.tensor(-W * 2, dtype=torch.float64)
    model = Potts(h=h, W=W.transpose(1, 2).reshape((L * 21, L * 21)))
    return model


class NucGapHelper():

    def __init__(self, aa_seq, gap_idx=20, device='cpu'):
        super().__init__()
        self.aa_seq = aa_seq
        self.L = len(self.aa_seq)  # length with gaps
        self.gap_idx = gap_idx
        self.non_gap_indices = [i for i, aa in enumerate(aa_seq) if aa != '-']
        self.device = device

    def add_gaps(self, design):
        bs = design.shape[0]
        X = torch.zeros(bs, self.L, 21, device=self.device)
        X[..., self.gap_idx] = 1
        X[:, self.non_gap_indices, :] = design
        return X


class NucGWGSampler(GWGCategoricalSampler):

    def __init__(self, bs, device, codon_to_aa_model, nuc_gap_helper):
        super().__init__(bs, 0, 0, device)
        self.codon_to_aa_model = codon_to_aa_model
        self.nuc_gap_helper = nuc_gap_helper
        self.temp = 1.0

    def compute_neg_energy_and_proposal(self, nts, model):
        nts.requires_grad_()
        bs = nts.shape[0]
        codons = nts.reshape((bs, -1, 12))
        aas = self.codon_to_aa_model(codons)
        aas_with_gaps = self.nuc_gap_helper.add_gaps(aas)
        energies = model(aas_with_gaps)
        f_x = -energies / self.temp
        grad_f_x = torch.autograd.grad(f_x.sum(), nts, retain_graph=True)[0]
        with torch.no_grad():
            d_tilde = (grad_f_x -
                       (nts * grad_f_x).sum(dim=-1).unsqueeze(dim=-1))
        probs = torch.softmax(d_tilde.reshape(d_tilde.shape[0], -1) / 2,
                              dim=-1)
        #dist = OneHotCategorical(logits=d_tilde.reshape(d_tilde.shape[0], -1) / 2)
        dist = OneHotCategorical(probs=probs)
        return f_x, dist


class EntNucGWGSampler(GWGCategoricalSampler):

    def __init__(self,
                 bs,
                 device,
                 codon_to_aa_model,
                 ent_nt_seq,
                 nuc_gap_helper_large,
                 nuc_gap_helper_small,
                 codon_table,
                 offset=1,
                 nt_idx_fix=None,
                 large_wt=0.5,
                 small_wt=0.5):
        super().__init__(bs, 0, 0, device)
        uppercase_idxs = [
            i for i, nt in enumerate(ent_nt_seq) if nt == nt.upper()
        ]
        self.smaller_gene_start = uppercase_idxs[offset]
        self.smaller_gene_end = uppercase_idxs[-1] + 1
        nts_smaller = ent_nt_seq[self.smaller_gene_start:self.smaller_gene_end]
        num_nts_smaller = len(nts_smaller)
        assert num_nts_smaller % 3 == 0

        self.codon_to_aa_model = codon_to_aa_model
        self.nuc_gap_helper_large = nuc_gap_helper_large
        self.nuc_gap_helper_small = nuc_gap_helper_small
        self.temp = 1.0
        self.large_wt = large_wt
        self.small_wt = small_wt

        stop_codons = [c for c, aa in codon_table.items() if aa == "*"]
        stop_codons_enc = torch.tensor([[BASE_TO_I[b] for b in codon]
                                        for codon in stop_codons])
        stop_codons_oh = torch.nn.functional.one_hot(
            stop_codons_enc,
            num_classes=4).to(torch.float).reshape(-1, 12).to(device)
        self.stop_codons_oh = stop_codons_oh

        constrain_mask = torch.zeros(len(ent_nt_seq), dtype=torch.bool)
        if not (nt_idx_fix is None):
            constrain_mask[nt_idx_fix] = True
        self.constrain_mask = constrain_mask.unsqueeze(1).unsqueeze(0).repeat(
            bs, 1, 4).to(device)

    def nt_to_aas(self, nts, nuc_gap_helper):
        bs = nts.shape[0]
        codons = nts.reshape((bs, -1, 12))
        aas = self.codon_to_aa_model(codons)
        aas_with_gaps = nuc_gap_helper.add_gaps(aas)
        return aas_with_gaps

    def compute_neg_energy_and_proposal(self, nts, model):
        model_large, model_small = model

        nts.requires_grad_()
        aas_with_gaps_large = self.nt_to_aas(nts, self.nuc_gap_helper_large)
        if False:
            energies_large = model_large(aas_with_gaps_large)
        plls_large = model_large.pseudolikelihood(aas_with_gaps_large)
        nts_small = nts[:, self.smaller_gene_start:self.smaller_gene_end, :]
        stop_codon = nts_small[:, -3:, :]  # last codon is stop
        penalty = -torch.prod(torch.abs(
            self.stop_codons_oh @ stop_codon.reshape(-1, 12).T - 3.0),
                              dim=0) * 1000

        aas_with_gaps_small = self.nt_to_aas(nts_small[:, :-3, :],
                                             self.nuc_gap_helper_small)
        if False:
            energies_small = model_small(aas_with_gaps_small)
        plls_small = model_small.pseudolikelihood(aas_with_gaps_small)
        if False:
            f_x = (((-energies_large * self.large_wt) +
                    (-energies_small * self.small_wt)) / self.temp) + penalty
        else:
            f_x = (((plls_large * self.large_wt) +
                    (plls_small * self.small_wt)) / self.temp) + penalty
        grad_f_x = torch.autograd.grad(f_x.sum(), nts, retain_graph=True)[0]
        with torch.no_grad():
            d_tilde = (grad_f_x -
                       (nts * grad_f_x).sum(dim=-1).unsqueeze(dim=-1))
            d_tilde -= (self.constrain_mask.to(torch.float) * 1000)

        probs = torch.softmax(d_tilde.reshape(d_tilde.shape[0], -1) / 2,
                              dim=-1)

        #dist = OneHotCategorical(logits=d_tilde.reshape(d_tilde.shape[0], -1) / 2)
        dist = OneHotCategorical(probs=probs)
        return f_x, dist


def add_methionine_constraints(model):
    h = model.h.reshape(model.L, model.A).clone()
    # for each amino acid index
    for i in range(model.A):
        # if the index is not methionine, add large penalty
        if I_TO_AA[i] != 'M':
            h[0, i] += 2000
    model_with_constraints = Potts(h=h, W=model.W.weight)
    return model_with_constraints


def get_pll_scores(nts_oh, sampler, model_large, model_small):
    with torch.no_grad():
        aas_large = sampler.nt_to_aas(nts_oh, sampler.nuc_gap_helper_large)
        # energies_large = model_large(aas_large)
        plls_large = model_large.pseudolikelihood(aas_large)

        nts_oh_small = nts_oh[:, sampler.smaller_gene_start:sampler.
                              smaller_gene_end, :]
        aas_small = sampler.nt_to_aas(nts_oh_small[:, :-3, :],
                                      sampler.nuc_gap_helper_small)
        # energies_small = model_small(aas_small)
        plls_small = model_small.pseudolikelihood(aas_small)
    return plls_large, plls_small


def optimize_ent_nuc(nts_oh,
                     sampler,
                     model_large,
                     model_small,
                     large_wt,
                     small_wt,
                     to_keep=10,
                     bs=1000):
    if nts_oh.ndim == 2:
        nts_oh = nts_oh.unsqueeze(0)
        nts_oh = nts_oh.repeat((bs, 1, 1)).to(device)

    model = (model_large, model_small)
    # anneal
    sampler.temp = 0.5
    for _ in range(100 * 2):
        sampler.temp *= 0.99
        nts_oh, _ = sampler.sample(nts_oh, model)

    # get_top `to_keep` sequences
    with torch.no_grad():
        aas_large = sampler.nt_to_aas(nts_oh, sampler.nuc_gap_helper_large)
        #energies_large = model_large(aas_large)
        plls_large = model_large.pseudolikelihood(aas_large)

        nts_oh_small = nts_oh[:, sampler.smaller_gene_start:sampler.
                              smaller_gene_end, :]
        aas_small = sampler.nt_to_aas(nts_oh_small[:, :-3, :],
                                      sampler.nuc_gap_helper_small)
        #energies_small = model_small(aas_small)
        plls_small = model_small.pseudolikelihood(aas_small)

    idxs = torch.argsort(-((large_wt * plls_large) +
                           (small_wt + plls_small)))[:to_keep]
    nts_oh = nts_oh[idxs].repeat(bs // len(idxs), 1, 1)

    # anneal again
    sampler.temp = 0.5
    for _ in range(100 * 2):
        sampler.temp *= 0.99
        nts_oh, _ = sampler.sample(nts_oh, model)
    return nts_oh


def check_constraints(nts_oh, sampler, codon_table):
    nts_oh_small = nts_oh[:, sampler.smaller_gene_start:sampler.
                          smaller_gene_end, :]
    first_codon_large = set(
        map(tuple, nts_oh[:, :3, :].argmax(dim=-1).cpu().numpy()))
    first_codon_small = set(
        map(tuple, nts_oh_small[:, :3, :].argmax(dim=-1).cpu().numpy()))
    for c in first_codon_large:
        c = ''.join([I_TO_BASE[i] for i in c])
        assert codon_table[c] == "M"
    for c in first_codon_small:
        c = ''.join([I_TO_BASE[i] for i in c])
        assert codon_table[c] == "M"

    last_codons_oh = set(
        map(tuple, nts_oh_small[:, -3:, :].argmax(dim=-1).cpu().numpy()))
    last_codons = [''.join([I_TO_BASE[i] for i in c]) for c in last_codons_oh]
    for c in last_codons:
        assert codon_table[c] == "*"
