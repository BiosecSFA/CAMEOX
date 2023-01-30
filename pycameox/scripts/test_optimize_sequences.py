import argparse
from pathlib import Path
import statistics as stats
from potts import Potts
import torch
import pandas as pd
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed

import copy
from pycameox.optimization import (load_ccmpred, load_codon_to_aa,
                                   NucGapHelper, add_methionine_constraints,
                                   EntNucGWGSampler, get_pll_scores,
                                   optimize_ent_nuc, check_constraints)

from pycameox.optimization import (AA_TO_I as aa_to_i, I_TO_AA as i_to_aa,
                                   BASE_TO_I as base_to_i, I_TO_BASE as
                                   i_to_base, CODON_TABLE as codon_table)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fn',
                        type=Path,
                        help="Input CSV File from outparser.jl")
    parser.add_argument('--out_fn', type=Path, help="Output CSV File")
    return parser.parse_args()


if __name__ == "__main__":
    columns = [
        "orig_fullseq", "weighted_score_before", "pll_large_before",
        "pll_small_before", "opt_fullseq", "weighted_score_after",
        "pll_large_after", "pll_small_after", "opt_fullseq_unc",
        "weighted_score_after_unc", "pll_large_after_unc",
        "pll_small_after_unc"
    ]
    args = parse_args()
    input_fn = args.input_fn
    out_fn = args.out_fn
    cameox_results = pd.read_csv(input_fn)

    run_dir = (input_fn.parent.parent.parent)
    PROTa = 'aroB_pf5_uref100'
    PROTb = 'infA_pf5_uref100'

    ccmpred_fn_infA = (run_dir / PROTb / Path(PROTb).with_suffix('.raw'))
    model_small = load_ccmpred(ccmpred_fn_infA)

    ccmpred_fn_aroB = (run_dir / PROTa / Path(PROTa).with_suffix('.raw'))
    model_large = load_ccmpred(ccmpred_fn_aroB)

    codon_to_aa_model = load_codon_to_aa()

    large_wt = 0.5
    small_wt = 0.5
    device = 'cuda'
    bs = 1000
    model_large.to(device)
    model_large_with_constraints = add_methionine_constraints(model_large)
    model_large_with_constraints.to(device)

    model_small.to(device)
    model_small_with_constraints = add_methionine_constraints(model_small)
    model_small_with_constraints.to(device)

    model = (model_large_with_constraints, model_small_with_constraints)

    num_sequences = len(cameox_results)
    pbar = tqdm(range(num_sequences))
    data = []
    for seq_idx in pbar:

        cameox_result = cameox_results.iloc[seq_idx]

        ent_nt_seq = cameox_result.full_seq
        if len(ent_nt_seq) != 1056:
            print()
            print(f"Skipping idx {seq_idx} due to weird issue with large")
            continue
        aas_large = cameox_result.loc[f"{PROTa}_seq"]
        aas_small = cameox_result.loc[f"{PROTb}_seq"]
        nuc_gap_helper_large = NucGapHelper(aas_large,
                                            gap_idx=aa_to_i['-'],
                                            device=device)
        nuc_gap_helper_small = NucGapHelper(aas_small,
                                            gap_idx=aa_to_i['-'],
                                            device=device)

        lowercase_idxs = [
            i for i, aa in enumerate(ent_nt_seq) if aa != aa.upper()
        ]
        nt_idx_fix = lowercase_idxs
        sampler = EntNucGWGSampler(bs,
                                   device,
                                   codon_to_aa_model.to(device),
                                   ent_nt_seq,
                                   nuc_gap_helper_large,
                                   nuc_gap_helper_small,
                                   codon_table,
                                   offset=1,
                                   nt_idx_fix=nt_idx_fix,
                                   large_wt=large_wt,
                                   small_wt=small_wt)

        nts_enc = torch.tensor([base_to_i[b] for b in ent_nt_seq.upper()])
        nts_oh = torch.nn.functional.one_hot(nts_enc, num_classes=4).to(
            torch.float).unsqueeze(0).to(device)

        nts_oh_small = nts_oh[:, sampler.smaller_gene_start:sampler.
                              smaller_gene_end, :]
        num_codons_small = nts_oh_small.shape[1] // 3
        num_non_gap_aas_small = len(nuc_gap_helper_small.non_gap_indices)
        if num_codons_small - 1 != num_non_gap_aas_small:
            print()
            print(f"Skipping idx {seq_idx} due to weird issue with small")
            continue

        plls_large_before, plls_small_before = get_pll_scores(
            nts_oh, sampler, model_large, model_small)
        plls_large_before, plls_small_before = plls_large_before.item(
        ), plls_small_before.item()
        score_before = ((large_wt * plls_large_before) +
                        (small_wt + plls_small_before))
        #print(f"Score: {score_before}", )

        nts_oh = optimize_ent_nuc(nts_oh.repeat((bs, 1, 1)),
                                  sampler,
                                  model_large_with_constraints,
                                  model_small_with_constraints,
                                  large_wt,
                                  small_wt,
                                  to_keep=10,
                                  bs=1000)

        plls_large, plls_small = get_pll_scores(nts_oh, sampler, model_large,
                                                model_small)
        idx = torch.argmax((large_wt * plls_large) + (small_wt + plls_small))

        plls_large_after = plls_large[idx].item()
        plls_small_after = plls_small[idx].item()
        score_after = ((large_wt * plls_large_after) +
                       (small_wt + plls_small_after))
        check_constraints(nts_oh, sampler, codon_table)

        optimized_nt = ''.join([
            i_to_base[i].lower()
            for i in nts_oh[idx].cpu().numpy().argmax(axis=-1)
        ])
        optimized_nt = ''.join([
            nt if i in lowercase_idxs else nt.upper()
            for i, nt in enumerate(optimized_nt)
        ])

        ############################################################
        ######## without constraining mutation region ##############
        ############################################################
        sampler = EntNucGWGSampler(bs,
                                   device,
                                   codon_to_aa_model.to(device),
                                   ent_nt_seq,
                                   nuc_gap_helper_large,
                                   nuc_gap_helper_small,
                                   codon_table,
                                   offset=1,
                                   nt_idx_fix=None,
                                   large_wt=large_wt,
                                   small_wt=small_wt)

        nts_enc = torch.tensor([base_to_i[b] for b in ent_nt_seq.upper()])
        nts_oh = torch.nn.functional.one_hot(nts_enc, num_classes=4).to(
            torch.float).unsqueeze(0).to(device)
        nts_oh = optimize_ent_nuc(nts_oh.repeat((bs, 1, 1)),
                                  sampler,
                                  model_large_with_constraints,
                                  model_small_with_constraints,
                                  large_wt,
                                  small_wt,
                                  to_keep=10,
                                  bs=1000)

        plls_large, plls_small = get_pll_scores(nts_oh, sampler, model_large,
                                                model_small)
        idx = torch.argmax((large_wt * plls_large) + (small_wt + plls_small))

        plls_large_after_unc = plls_large[idx].item()
        plls_small_after_unc = plls_small[idx].item()
        score_after_unc = ((large_wt * plls_large_after_unc) +
                           (small_wt + plls_small_after_unc))
        check_constraints(nts_oh, sampler, codon_table)

        optimized_nt_unc = ''.join([
            i_to_base[i].lower()
            for i in nts_oh[idx].cpu().numpy().argmax(axis=-1)
        ])
        optimized_nt_unc = ''.join([
            nt if i in lowercase_idxs else nt.upper()
            for i, nt in enumerate(optimized_nt)
        ])

        #print(f"Score: {score_after}")
        pbar.set_description(
            f"Score Before: {score_before:.3f}    Score After: {score_after:.3f}    Score After Unc: {score_after_unc:.3f}",
            refresh=True)

        result = (ent_nt_seq, score_before, plls_large_before,
                  plls_small_before, optimized_nt, score_after,
                  plls_large_after, plls_small_after, optimized_nt_unc,
                  score_after_unc, plls_large_after_unc, plls_small_after_unc)

        data.append(result)
        df = pd.DataFrame(data, columns=columns)
        df.to_csv(out_fn, index=False)
    df = pd.DataFrame(data, columns=columns)
    df.to_csv(out_fn, index=False)
