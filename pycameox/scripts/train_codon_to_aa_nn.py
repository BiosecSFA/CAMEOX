import torch
import torch.nn as nn
import torch.nn.functional as F
from tqdm import tqdm
from torch.utils.data import TensorDataset
from superpolyak import (
    AlternatingProjections,
    BundleLinearSystemSolver,
    NewtonCG,
    ProxGradient,
    SuperPolyak,
)
from superpolyak.util import superpolyak_coupled_with_fallback

bases = [l.upper() for l in 'tcag']
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

base_to_i = {base: i for i, base in enumerate("ACGT")}
i_to_base = {i: base for i, base in enumerate("ACGT")}
oh_codon_aa_pairs = []
for codon in codons:
    aa = codon_table[codon]
    if aa == "*":
        continue
    aa_enc = torch.tensor([aa_to_i[aa]])
    aa_oh = torch.nn.functional.one_hot(aa_enc, num_classes=len(ALPHABET)).to(
        torch.float)
    codon_enc = torch.tensor(
        [base_to_i[codon[0]], base_to_i[codon[1]], base_to_i[codon[2]]])
    codon_oh = torch.nn.functional.one_hot(codon_enc,
                                           num_classes=4).to(torch.float)
    oh_codon_aa_pairs.append((codon_oh.reshape(-1), aa_oh))
codons_oh, aas_oh = zip(*oh_codon_aa_pairs)
train_ds = TensorDataset(torch.vstack(codons_oh), torch.vstack(aas_oh))

hidden = 100
codon_to_aa_model = nn.Sequential(nn.Linear(12, hidden), nn.ReLU(),
                                  nn.Linear(hidden, hidden), nn.ReLU(),
                                  nn.Linear(hidden, 21))

codon_to_aa_model.train()
'''
First we train with SGD
'''

lr = 1e-3
optim = torch.optim.Adam(codon_to_aa_model.parameters(), lr=lr)
n_epoch = 1000
X, y = train_ds[:]
X, y = X.cuda(), y.cuda()
codon_to_aa_model.cuda()
loss_fn = nn.L1Loss()
epoch_iter = tqdm(range(n_epoch))
for _ in epoch_iter:
    optim.zero_grad()
    y_hat = codon_to_aa_model(X)
    loss = F.l1_loss(y_hat, y)
    loss += F.mse_loss(y_hat, y)
    loss.backward()
    optim.step()
    epoch_iter.set_postfix(loss=loss.item())
'''
Then we finish the optimization with superpolyak
'''
codon_to_aa_model.cpu()
max_elts = 200
eta_est = 0.5
max_oracle_calls = 10000

linsys_solver = BundleLinearSystemSolver.QR

# optimizer = torch.optim.LBFGS(codon_to_aa_model.parameters(),
#                               history_size=100,
#                               max_iter=40,
#                               line_search_fn="strong_wolfe")
optimizer = SuperPolyak([p for p in codon_to_aa_model.parameters()],
                        max_elts=max_elts,
                        eta_est=eta_est,
                        linsys_solver=linsys_solver)


# Closure function to allow us to call backward.
def closure():
    optimizer.zero_grad()
    y_hat = codon_to_aa_model(X)
    loss = F.l1_loss(y_hat, y)
    loss += F.mse_loss(y_hat, y)
    loss.backward()
    return loss


X, y = X.cpu(), y.cpu()
codon_to_aa_model.cpu()
epoch_iter = tqdm(range(40))
for _ in epoch_iter:
    loss = optimizer.step(closure)
    epoch_iter.set_postfix(loss=loss)
torch.save(codon_to_aa_model.state_dict(), "codon_to_aa_model.pt")
