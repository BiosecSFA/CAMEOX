"""Python module for CAMEOX downstream pipeline"""

# ## Initialization
# ### Dependencies
from collections import Counter
import difflib
import re
import typing
from typing import Set, List, Any, Dict, Iterator, Optional, NamedTuple

from pycameox.config import RunsSet, SampleSet

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.offline as py
import scipy


# ### General methods
# #### Get ~ratio
def get_ratio(seqA, seqB):
    return difflib.SequenceMatcher(None, seqA, seqB, autojunk=False).ratio()


# #### Get redundancy dataframe for a set of CAMEOS runs
def get_redundancy(datasets: RunsSet,
                   protA: str,
                   protB: str,
                   proteins: Dict[str, str],
                   source: str = 'supp',
                   verbose: bool = False) -> Optional[pd.DataFrame]:
    """Get redundancy dataframe for a set of CAMEOX runs.

    In the obtained dataframe:
     - Column 'max_multy' is the maximum multiplicity of the variant accounting
     for all the samples.
     - Column 'abundance' is the count of the variant accounting for all the
     samples where it appears twice or more times.
     - Column 'spread' is the fraction of datasets (CAMEOX runs) where the
     variant appears twice or more times.
    """

    def vprint(*arguments, **kargs) -> None:
        """Print only if verbose mode is enabled"""
        if verbose:
            print(*arguments, **kargs)

    seqA: str = protA + '_' + ('trans' if source == 'supp' else 'seq')
    seqB: str = protB + '_' + ('trans' if source == 'supp' else 'seq')

    # Get set of dataframes with redundant variants for each run
    redundant: RunsSet = RunsSet({})
    for name, df in datasets.items():
        # vardup contains only duplicated rows
        vardup = df.loc[df.duplicated(subset=[seqA, seqB], keep=False)].copy()
        if vardup.empty:
            vprint(f' No redundant variants for {name}')
            continue
        dupcnt = vardup.value_counts(subset=[seqA, seqB])
        vardup['multiplicity'] = vardup.apply(
            lambda row: dupcnt[row[seqA],
                               row[seqB]],
            axis='columns', raw=False)
        redundant[name] = vardup.drop_duplicates(
            subset=[seqA, seqB],
            ignore_index=True).sort_values(
            by='multiplicity', ascending=False, ignore_index=True)
        vprint(
            f' For {name}, {len(redundant[name])} redundant variants with '
            f'maximum multiplicity of {redundant[name]["multiplicity"][0]}')
    if not redundant:
        vprint('NOTE: No redundancy found!')
        return(None)

    # Get single redundancy dataframe including total abundance and relative
    # spread of redundant variants
    redundancy: pd.DataFrame = pd.concat(
        [df for df in redundant.values()],
        ignore_index=True)
    redundancy['abundance'] = redundancy.groupby(
        by=[seqA, seqB])['multiplicity'].transform(sum)
    redundancy['spread'] = redundancy.groupby(
        by=[seqA, seqB])['multiplicity'].transform(len) / len(datasets)
    redundancy.sort_values(
        by='multiplicity',
        ascending=False,
        ignore_index=True,
        inplace=True)
    redundancy.drop_duplicates(
        subset=[
            seqA,
            seqB],
        keep='first',
        ignore_index=True,
        inplace=True)
    redundancy.sort_values(
        by='abundance',
        ascending=True,
        ignore_index=True,
        inplace=True)
    redundancy.rename(
        columns={
            'multiplicity': 'max_multy'},
        inplace=True,
        errors='raise')
    vprint('TOP maximum multiplicity:', max(redundancy['max_multy']))
    vprint('TOP abundance:', max(redundancy['abundance']))
    vprint(f'TOP pread: {max(redundancy["spread"]):.2%}')

    # Calculate ~ratio
    vprint('Now calculating similarity ratios:')
    vprint(f'\tCalculating ratios of {protA}... ', end='', flush=True)
    ratios_A: pd.Series = redundancy[seqA].apply(
        get_ratio, args=(proteins[protA],))
    vprint(f'OK! \n\tCalculating ratios of {protB}... ', end='', flush=True)
    ratios_B: pd.Series = redundancy[seqB].apply(
        get_ratio, args=(proteins[protB],))
    vprint('OK!')
    redundancy[f'ratio_{protA}'] = ratios_A
    redundancy[f'ratio_{protB}'] = ratios_B

    # Calculate ERP
    renucl = re.compile('A|T|C|G')
    vprint('Finally, calculating entanglement relative position...')
    erp: pd.Series
    try:
        erp = redundancy['full_seq'].apply(
            lambda fullseq: round(
                (renucl.search(fullseq).start() / len(fullseq)), 3))
        redundancy['ERP'] = erp
    except AttributeError:
        raise

    return(redundancy)


# #### Interactive plot for CAMEOS results extracted with outparser.jl
def general_iplot(datasets: RunsSet,
                  protA: str,
                  protB: str,
                  source: str = 'supp',
                  # Dict with keys 'name', 'x', and 'y'
                  density: Dict[str, pd.Series] = None,
                  title: str = None,
                  fout: str = None,
                  auto_open: bool = True) -> None:

    pslsA: str = protA + '_psls' + ('_S6' if source == 'supp' else '')
    pslsB: str = protB + '_psls' + ('_S6' if source == 'supp' else '')
    simrA: str = 'ratio_' + protA
    simrB: str = 'ratio_' + protB

    # Prepare data for density plot if not explicitly entered
    if density is None:
        density = {}
        density['name'] = 'Aggregate'
        density['x'] = [
            pd.concat(
                [dataset[pslsA] for dataset in datasets.values()],
                ignore_index=True)]
        density['y'] = [
            pd.concat(
                [dataset[pslsB] for dataset in datasets.values()],
                ignore_index=True)]
    elif 'name' not in density or 'x' not in density or 'y' not in density:
        raise KeyError(
            'ERROR: If set, density dict must have keys "name", "x", and "y"!')

    data: List[Any] = []
    for name, dataset in datasets.items():
        trace = go.Scatter(
            x=dataset[pslsA],
            y=dataset[pslsB],
            text=(
                #                f'<i>Dataset</i>: <b>{name}</b>' + '<br>'
                '<i>Id</i>: ' + dataset.index.astype(str) + '<br>'
                + f'<i>apll {protA}</i>: ' + dataset[pslsA].apply(
                    lambda v: f'{v:.2f}<br>')
                + f'<i>apll {protB}</i>: ' + dataset[pslsB].apply(
                    lambda v: f'{v:.2f}<br>')
                + f'<i>~ratio {protA}</i>: ' + \
                dataset[simrA].apply(lambda v: f'{v:.0%}<br>')
                + f'<i>~ratio {protB}</i>: ' + \
                dataset[simrB].apply(lambda v: f'{v:.0%}')
            ),
            name=name,
            hovertemplate=('%{text}'),
            opacity=0.7,
            mode='markers',
            marker={
                'symbol': 'circle',
                'size': 8,
                'coloraxis': 'coloraxis',
                'showscale': True,
                'line': {
                    'color': 'rgba(0, 0, 0, 0.5)',
                    'width': 0},
            },
        )
        if source == 'supp':  # Special format for supp figure
            if 'success' in name:
                trace['marker']['line']['color'] = 'rgba(0, 1.0, 0, 1.0)'
                trace['marker']['line']['width'] = 2.0
            elif '7500' in name:
                trace['marker']['symbol'] = 'square'
                trace['marker']['line']['color'] = 'rgba(0, 0, 0, 0.4)'
                trace['marker']['line']['width'] = 1.0
            elif '2500' in name:
                trace['marker']['symbol'] = 'diamond'
                trace['marker']['line']['color'] = 'rgba(0, 0, 0, 0.4)'
                trace['marker']['line']['width'] = 1.0
        if name.startswith('intersect'):
            trace['marker']['symbol'] = 'x'
            trace['marker']['line']['color'] = 'rgba(0, 0, 0, 0.4)'
            trace['marker']['line']['width'] = 1.0
        data.append(trace)

    layout: Dict[str,
                 Any] = {'legend': {'x': +0.1,
                                    'y': -0.2},
                         'legend_orientation': 'h',
                         'title': title,
                         'xaxis': {'title': f'anti-pseudo-log-likelihood {protA} (rel. to WT)',
                                   'type': 'linear',
                                   },
                         'yaxis': {'title': f'anti-pseudo-log-likelihood {protB} (rel. to WT)',
                                   'type': 'linear',
                                   },
                         'coloraxis': {'colorscale': ('plasma' if source == 'supp' else 'Hot'),
                                       'colorbar': {'showticklabels': False,
                                                    'title': '',
                                                    },
                                       },
                         'hovermode': 'closest',
                         'updatemenus': [{'buttons': [{"args": [{'marker.color': None},
                                                                {'coloraxis.colorscale': None,
                                                                 'coloraxis.colorbar.title': '',
                                                                 'coloraxis.colorbar.showticklabels': False,
                                                                 }],
                                                       "label": "None",
                                                       "method": "update"},
                                                      {"args": [{'marker.color': [dataset[[simrA,
                                                                                           simrB]].mean(axis='columns') for dataset in datasets.values()]},
                                                                {'coloraxis.colorscale': 'Hot',
                                                                 'coloraxis.colorbar.title': '~ratio avg',
                                                                 'coloraxis.colorbar.showticklabels': True,
                                                                 }],
                                                       "label": "By ~ratio avg",
                                                       "method": "update"},
                                                      ],
                                          'direction': "right",
                                          'pad': {"r": 10,
                                                  "t": 10},
                                          'active': 0,
                                          'showactive': True,
                                          'type': "buttons",
                                          'x': 0.75,
                                          'xanchor': "left",
                                          'y': -0.17,
                                          'yanchor': "bottom"},
                                         {"buttons": [{"args": [{'type': 'scatter',
                                                                 'name': [name for name in datasets.keys()],
                                                                 'x': [dataset[pslsA] for dataset in datasets.values()],
                                                                 'y': [dataset[pslsB] for dataset in datasets.values()],
                                                                 'hovertemplate':'%{text}',
                                                                 },
                                                                {'xaxis.title': f'Anti-pseudo-log-likelihood for {protA} (rel. to WT)',
                                                                 'xaxis.type': 'lin',
                                                                 'yaxis.title': f'Anti-pseudo-log-likelihood for {protB} (rel. to WT)',
                                                                 'yaxis.type': 'lin',
                                                                 'updatemenus[0].buttons[0].visible': True,
                                                                 'updatemenus[0].buttons[1].visible': True,
                                                                 'updatemenus[0].buttons[2].visible': True,
                                                                 'annotations[1].visible': True,
                                                                 }],
                                                       "label": "APLL",
                                                       "method": "update"},
                                                      {"args": [{'type': 'scattergl',
                                                                 'name': [name for name in datasets.keys()],
                                                                 'x': [dataset[pslsA].apply(lambda s: np.exp(-s)) for dataset in datasets.values()],
                                                                 'y': [dataset[pslsB].apply(lambda s: np.exp(-s)) for dataset in datasets.values()],
                                                                 'hovertemplate':'%{text}',
                                                                 },
                                                       {'xaxis.title': f'pseudo-likelihood for {protA} (rel. to WT)',
                                                        'xaxis.type': 'log',
                                                        'yaxis.title': f'pseudo-likelihood for {protB} (rel. to WT)',
                                                        'yaxis.type': 'log',
                                                        'updatemenus[0].buttons[0].visible': True,
                                                        'updatemenus[0].buttons[1].visible': True,
                                                        'updatemenus[0].buttons[2].visible': True,
                                                        'annotations[1].visible': True,
                                                        }],
                                                       "label": "PL",
                                                       "method": "update"},
                                                      {"args": [{'type': 'scatter',
                                                                 'name': [name for name in datasets.keys()],
                                                                 'x': [dataset[simrA] for dataset in datasets.values()],
                                                                 'y': [dataset[simrB] for dataset in datasets.values()],
                                                                 'hovertemplate':'%{text}',
                                                                 },
                                                       {'xaxis.title': f'Similarity ratio for {protA}',
                                                        'xaxis.type': 'lin',
                                                        'yaxis.title': f'Similarity ratio for {protB}',
                                                        'yaxis.type': 'lin',
                                                        'updatemenus[0].buttons[0].visible': True,
                                                        'updatemenus[0].buttons[1].visible': True,
                                                        'updatemenus[0].buttons[2].visible': True,
                                                        'annotations[1].visible': True,
                                                        }],
                                                       "label": "~ratio",
                                                       "method": "update"},
                                                      {"args": [{'type': 'histogram2dcontour',
                                                                 'name': density['name'],
                                                                 'x': density['x'],
                                                                 'y': density['y'],
                                                                 'histnorm':'probability density',
                                                                 'hovertemplate':'',
                                                                 },
                                                       {'xaxis.title': f'Normalized density of variants in the anti-pseudo-log-likelihood space for {protA} (rel. to WT)',
                                                        'xaxis.type': 'lin',
                                                        'yaxis.title': f'Normalized density of variants in the anti-pseudo-log-likelihood space for {protB} (rel. to WT)',
                                                        'yaxis.type': 'lin',
                                                        'updatemenus[0].buttons[0].visible': False,
                                                        'updatemenus[0].buttons[1].visible': False,
                                                        'updatemenus[0].buttons[2].visible': False,
                                                        'annotations[1].visible': False,
                                                        }],
                                                       "label": "Density",
                                                       "method": "update"},
                                                      ],
                                          'direction': "right",
                                          'pad': {"r": 10,
                                                  "t": 10},
                                          'active': 0,
                                          'showactive': True,
                                          'type': "buttons",
                                          'x': 0.05,
                                          'xanchor': "auto",
                                          'y': -0.17,
                                          'yanchor': "bottom"},
                                         ],
                         'annotations': [{'text': 'Plot:',
                                          'showarrow': False,
                                          'x': 0.015,
                                          'xref': 'paper',
                                          'xanchor': "auto",
                                          'y': -0.16,
                                          'yref': 'paper',
                                          'yanchor': "bottom"},
                                         {'text': 'Color-scale:',
                                          'showarrow': False,
                                          'x': 0.75,
                                          'xref': 'paper',
                                          'xanchor': "auto",
                                          'y': -0.16,
                                          'yref': 'paper',
                                          'yanchor': "bottom"},
                                         ],
                         }

    # Add ERP colorscale button if there is ERP data
    if all('ERP' in dataset for dataset in datasets.values()):
        layout['updatemenus'][0]['buttons'].append({
            "args": [{'marker.color': [dataset['ERP'] for dataset in datasets.values()]},
                     {'coloraxis.colorscale': 'Bluered', 'coloraxis.colorbar.title': 'ERP',
                     'coloraxis.colorbar.showticklabels': True, }],
            "label": "By ERP",
            "method": "update"
        },)

    # Deal with fout choices
    if fout is None:
        layout['autosize'] = False
        layout['width'] = 960
        layout['height'] = 800

    pyfig: Dict[str, Any] = {'data': data, 'layout': layout}

    if fout is None:
        py.iplot(pyfig)
    else:
        py.plot(pyfig, filename=fout, auto_open=auto_open)


# #### Interactive plot for redundancy/multiplicity analysis
def redundancy_iplot(redundancy: pd.DataFrame,
                     protA: str,
                     protB: str,
                     name: str = 'redundancy',
                     title: str = None,
                     fout: str = None,
                     auto_open: bool = True) -> None:

    pslsA: str = protA + '_psls'
    pslsB: str = protB + '_psls'
    simrA: str = 'ratio_' + protA
    simrB: str = 'ratio_' + protB

    data: List[Any] = []
    trace = go.Scatter(
        x=redundancy[pslsA],
        y=redundancy[pslsB],
        text=('<b><i>Abundance</i></b>: '
              + redundancy['abundance'].apply(lambda v: f'{v:d}<br>')
              + '<b><i>Spread</i></b>: '
              + redundancy['spread'].apply(lambda v: f'{v:.0%}<br>')
              + f'<i>apll {protA}</i>: ' + redundancy[pslsA].apply(
                  lambda v: f'{v:.2f}<br>')
              + f'<i>apll {protB}</i>: ' + redundancy[pslsB].apply(
                  lambda v: f'{v:.2f}<br>')
                + f'<i>~ratio {protA}</i>: ' + redundancy[simrA].apply(
                    lambda v: f'{v:.0%}<br>')
                + f'<i>~ratio {protB}</i>: '
              + redundancy[simrB].apply(lambda v: f'{v:.0%}')
              ),
        name=name,
        hovertemplate=('%{text}'),
        #        opacity=0.7,
        mode='markers',
        marker={
            'symbol': 'circle',
            'color': redundancy['abundance'],
            'size': redundancy['spread'],
            'sizeref': 0.02,
            'sizemode': 'diameter',
            'coloraxis': 'coloraxis',
            'showscale': True,
        },
    )
    data.append(trace)

    layout: Dict[str, Any] = {
        'legend': {'x': +0.1, 'y': -0.2},
        'legend_orientation': 'h',
        'title': title,
        'xaxis': {
            'title': f'anti-pseudo-log-likelihood {protA} (rel. to WT)',
            'type': 'linear',
        },
        'yaxis': {
            'title': f'anti-pseudo-log-likelihood {protB} (rel. to WT)',
            'type': 'linear',
        },
        'coloraxis': {
            'colorscale': ('Reds'),
            'colorbar': {
                'showticklabels': True,
                'title': 'Abundance',
            },
        },
        'hovermode': 'closest',
        'updatemenus': [
            {
                'buttons': [
                    {
                        "args": [{'marker.color': [redundancy['abundance']]},
                                 {'coloraxis.colorscale': 'Reds', 'coloraxis.colorbar.title': 'Abundance',
                                  'coloraxis.colorbar.showticklabels': True, }],
                        "label": "Multiplicity",
                        "method": "update"
                    },
                    {
                        "args": [{'marker.color': [redundancy[[simrA, simrB]].mean(axis='columns')]},
                                 {'coloraxis.colorscale': 'Hot', 'coloraxis.colorbar.title': '~ratio avg',
                                  'coloraxis.colorbar.showticklabels': True, }],
                        "label": "By ~ratio avg",
                        "method": "update"
                    },
                ],
                'direction': "right",
                'pad': {"r": 10, "t": 10},
                'active': 0,
                'showactive': True,
                'type': "buttons",
                'x': 0.70,
                'xanchor': "left",
                'y': -0.17,
                'yanchor': "bottom"
            },
        ],
        'annotations': [
            {
                'text': 'Color-scale:',
                'showarrow': False,
                'x': 0.70,
                'xref': 'paper',
                'xanchor': "auto",
                'y': -0.16,
                'yref': 'paper',
                'yanchor': "bottom"
            },
        ],
    }

    # Add ERP colorscale button if there is ERP data
    if 'ERP' in redundancy:
        layout['updatemenus'][0]['buttons'].append({
            "args": [{'marker.color': [redundancy['ERP']]},
                     {'coloraxis.colorscale': 'Bluered', 'coloraxis.colorbar.title': 'ERP',
                     'coloraxis.colorbar.showticklabels': True, }],
            "label": "By ERP",
            "method": "update"

},)

    # Deal with fout choices
    if fout is None:
        layout['autosize'] = False
        layout['width'] = 960
        layout['height'] = 800

    pyfig: Dict[str, Any] = {'data': data, 'layout': layout}

    if fout is None:
        py.iplot(pyfig)
    else:
        py.plot(pyfig, filename=fout, auto_open=auto_open)


# ## Sample variants for experimental validation
# ### Sampling criterium: Pareto's front


def sample_pareto_front(dset: pd.DataFrame,
                        protA: str,
                        protB: str,
                        size: int = 500,
                        erp: float = None,
                        verbose: bool = False) -> pd.DataFrame:
    """Sample variants on Pareto front for a set of CAMEOX runs and for a given ERP.
    """

    def vprint(*arguments, **kargs) -> None:
        """Print only if verbose mode is enabled"""
        if verbose:
            print(*arguments, **kargs)

    pslsA: str = protA + '_psls'
    pslsB: str = protB + '_psls'

    if erp is not None:  # Restrict by erp if given
        dset = dset.loc[dset['ERP'] == erp]
    print("Number of different CAMEOX solutions for this "
          f"{'entanglement' if erp is None else 'ERP'}: {len(dset)}")

    # Get "top pair" (defined as the one with minimum pslsA + pslsB)
    pd.options.mode.chained_assignment = None
    dset.loc[:,'psls_sum'] = dset[pslsA] + dset[pslsB]
    dset.sort_values(
        by='psls_sum',
        ascending=True,
        ignore_index=True,
        inplace=True)
    pd.options.mode.chained_assignment = 'warn'
    vprint(
        f"Top pair: {pslsA}={dset.loc[0][pslsA]:.2f}, {pslsB}={dset.loc[0][pslsB]:.2f}")

    # Define aux variables and dataframes sorted by pslsA and pslsB
    bypslsA = dset.sort_values(by=pslsA, ascending=True)
    top_idx_bypslsA = bypslsA.index.get_loc(0)
    vprint(
        f'There are {top_idx_bypslsA} variants with better {pslsA} than the top pair')
    bypslsB = dset.sort_values(by=pslsB, ascending=True)
    top_idx_bypslsB = bypslsB.index.get_loc(0)
    vprint(
        f'There are {top_idx_bypslsB} variants with better {pslsB} than the top pair')

    # **Algorithm core**
    # Explore Pareto's front from center (top pair) towards extremes ("above" and "below") in psls (aka APLL/NPLL)
    #
    # Initialization
    pareto_points: Dict[int, pd.Series] = {}
    cntr_above: int = 0  # Counter of variants in Pareto's front "above" the top
    cntr_below: int = 0  # Counter of variants in Pareto's front "below" the top
    idxA: int
    idxB: int
    idx_bypslsA: int = top_idx_bypslsA
    idx_bypslsB: int = top_idx_bypslsB
    vprint(
        f' 0\t{idx_bypslsA}!\t{idx_bypslsB}^:\t{bypslsA.iloc[idx_bypslsA][pslsA]:.2f}\t{bypslsB.iloc[idx_bypslsB][pslsB]:.2f}')
    # Add "top pair" as 1st (zero) result
    pareto_points[0] = bypslsA.iloc[idx_bypslsA]
    # Main loop: go until too many variants or reach limits for the APLL values
    while cntr_above - cntr_below < size and (
            idx_bypslsA + 1 < len(bypslsA) or idx_bypslsB + 1 < len(bypslsB)):
        # Get pair in front "below"
        for idxA in range(idx_bypslsA + 1, len(bypslsA)):
            if bypslsA.iloc[idxA][pslsB] < bypslsA.iloc[idx_bypslsA][pslsB]:
                cntr_below -= 1
                vprint(
                    f'{cntr_below}\t{idx_bypslsA}\t{idxA}:'
                    f'\t{bypslsA.iloc[idxA][pslsA]:.2f}\t{bypslsA.iloc[idxA][pslsB]:.2f}'
                    f'\t\t{idx_bypslsA+1} of {len(bypslsA)}')
                pareto_points[cntr_below] = bypslsA.iloc[idxA]
                break
        idx_bypslsA = idxA

        # Get pair in front "above"
        for idxB in range(idx_bypslsB + 1, len(bypslsB)):
            if bypslsB.iloc[idxB][pslsA] < bypslsB.iloc[idx_bypslsB][pslsA]:
                cntr_above += 1
                vprint(
                    f'{cntr_above}\t{idx_bypslsB}\t{idxB}:'
                    f'\t{bypslsB.iloc[idxB][pslsA]:.2f}\t{bypslsB.iloc[idxB][pslsB]:.2f}'
                    f'\t\t{idx_bypslsB+1} of {len(bypslsB)}')
                pareto_points[cntr_above] = bypslsB.iloc[idxB]
                break
        idx_bypslsB = idxB
    print(f'Found {cntr_above - cntr_below + 1} variants in Pareto front: top, {cntr_above} above, {-cntr_below} below')

    return pd.DataFrame.from_dict(pareto_points, orient='index')


# ### Sampling criterium: multiplicity
def sample_multiplicity(dset: pd.DataFrame,
                        protA: str,
                        protB: str,
                        size: int = 500,
                        erp: float = None,
                        sampled: SampleSet = None,
                        verbose: bool = False) -> pd.DataFrame:
    """Sample variants with multiplicity.
    """

    def vprint(*arguments, **kargs) -> None:
        """Print only if verbose mode is enabled"""
        if verbose:
            print(*arguments, **kargs)

    if erp is not None:  # Restrict by erp if given
        dset = dset[dset['ERP'] == erp]
    print("Number of different redundant CAMEOX solutions for this "
          f"{'entanglement' if erp is None else 'ERP'}: {len(dset)}")

    if sampled is None:
        return dset.iloc[:size]
    else:
        all_sampled: pd.DataFrame = pd.concat(
            [df for df in sampled.values()],
            ignore_index=True)
        all_sampled_full_seqs: Set[str] = set(all_sampled['full_seq'])
        print(
            f'INFO: Checking against {len(all_sampled_full_seqs)} already sampled variants')
        num_rows: List[int] = []
        for num, row in enumerate(dset.itertuples(index=False)):
            if row.full_seq not in all_sampled_full_seqs:
                num_rows.append(num)
            else:
                vprint(
                    f'  Found already sampled var for redundant var {num+1} with abundance={row.abundance} '
                    f'and max_multiplity={row.max_multy}. Skipping...')
            if len(num_rows) >= size:
                break
        print(
            f'Succesfully sampled {len(num_rows)} unique redundant variants.')
        return dset.take(num_rows)


# ### Sampling criterium: overdensities in APLL space
def inv_avg_pdist_com(num_stds,
                      H: np.ndarray,
                      H_mean: float,
                      H_std: float,
                      desired_features: int = 3,
                      penalization: float = 0.1,
                      debug: bool = False,
                      results: List[Dict[str, Any]] = None,
                      ) -> float:
    """Min target: inverse averaged pairwise distance for centers of mass"""
    H_zeroclipped = np.clip(
        H - (H_mean + num_stds * H_std),
        0, 1)  # Negative values to zero
    (L, num_feat) = scipy.ndimage.label(H_zeroclipped)
    com_avdist: float = 0.0
    if num_feat < 1:
        if results is not None:
            results.append(
                {'num_feat': num_feat, 'weight': 0.0,
                 'com_avdist': com_avdist, 'L': L})
        if debug:
            print(f' inv_avg_pdist_com: #feat={num_feat} -> return max value!')
        # Deliver maximum value for undesired num_features
        return (float(np.finfo(float).max) / 10)
    com = scipy.ndimage.center_of_mass(
        H_zeroclipped, L, range(
            1, num_feat + 1))
    com_pdist = scipy.spatial.distance.pdist(com, 'seuclidean')
    com_avdist = com_pdist.mean()
    # Add penalization for each diverging feature from the desired number
    weight: float = 1 + penalization * abs(desired_features - num_feat)
    if debug:
        print(
            f' inv_avg_pdist_com: #feat={num_feat}\t weight={weight}\t com_avdist={com_avdist}')
    if results is not None:
        results.append(
            {'num_feat': num_feat, 'weight': weight,
             'com_avdist': com_avdist, 'L': L})
    return weight / com_avdist


def search_overdensities(H: np.ndarray,
                         overdensities: int = 3,
                         penalization: float = 0.1,
                         std_max: float = 5.0,  # Max std from mean to explore
                         verbose: bool = False,
                         debug: bool = False) -> pd.DataFrame:
    """Search for different number of overdensities"""

    # Mask bins with no variants (0.0) in the APLL space so they are not
    # included in statistics
    H_ma = np.ma.masked_values(H, 0.0)
    # Prepare args for the objective function
    H_mean: float = H_ma.mean()
    H_std: float = H_ma.std()
    desired_features: int = overdensities
    results: List[Dict[str, Any]] = []
    # Call to the minimizator
    res = scipy.optimize.minimize_scalar(
        inv_avg_pdist_com, bounds=(0, 5), method='bounded',
        options={'disp': (3 if verbose else 0)},
        args=(H, H_mean, H_std, desired_features, penalization, debug, results))

    # Get the best results from the optimum
    best_results: List[Dict[str, Any]] = []
    inv_avg_pdist_com(
        res.x,
        H,
        H_mean,
        H_std,
        desired_features,
        penalization,
        debug,
        results=best_results)
    # num_features = best_results[0]['num_features']
    L = best_results[0]['L']
    locs: List[slice] = scipy.ndimage.find_objects(L)

    # Return the OptimizeResult object and a list of dicts with results for
    # every call and locs for the best
    return res, results, locs


def plot_overdensities(X: np.ndarray,
                       Y: np.ndarray,
                       H: np.ndarray,
                       locs: List[List[slice]],
                       protA: str,
                       protB: str,
                       axlim: List[float] = None,
                       ) -> plt.Axes:
    """Plot selection of overdensities."""

    num_features: int = len(locs)
    qmesh = plt.pcolormesh(X, Y, H)
    ax = plt.gca()
    ax.set_aspect('equal')
    cbar = plt.colorbar(qmesh, ax=ax)
    cbar.formatter.set_powerlimits((0, 0))
    for f in range(num_features):
        sliceX = locs[f][1]
        sliceY = locs[f][0]
        pslsA_min = X[0][sliceX.start]
        pslsA_max = X[0][sliceX.stop]
        pslsB_min = Y[:, 0][sliceY.start]
        pslsB_max = Y[:, 0][sliceY.stop]
        plt.plot([pslsA_min, pslsA_min, pslsA_max, pslsA_max, pslsA_min],
                 [pslsB_min, pslsB_max, pslsB_max, pslsB_min, pslsB_min],
                 color="red")
    plt.title("2D histogram of the APLL with detected\noverdensity regions")
    plt.xlabel(f'APLL {protA} (rel. to WT)')
    plt.ylabel(f'APLL {protB} (rel. to WT)')
    if axlim is not None:
        plt.axis(axlim)
    plt.tight_layout()
    plt.show()
    return ax


def sample_overdensities(dset: pd.DataFrame,
                         protA: str,
                         protB: str,
                         size: int = 500,
                         bins: int = 25,
                         overdensities: int = 3,
                         penalization: float = 0.1,
                         erp: float = None,
                         sampled: SampleSet = None,
                         verbose: bool = False,
                         debug: bool = False) -> pd.DataFrame:
    """Sample variants from APLL overdensities.
    """

    def vprint(*arguments, **kargs) -> None:
        """Print only if verbose mode is enabled"""
        if verbose:
            print(*arguments, **kargs)

    pslsA: str = protA + '_psls'
    pslsB: str = protB + '_psls'

    if erp is not None:  # Restrict by erp if given
        dset = dset[dset['ERP'] == erp]
    print("Number of different CAMEOX solutions for this "
          f"{'entanglement' if erp is None else 'ERP'}: {len(dset)}")

    # Limit with iloc was needed to avoid not enough memory before rolling
    # update TOSS-3.7-17 on mammoth
    nrpsls = dset[[pslsA, pslsB]]  # .iloc[:50000]

    H, xedges, yedges = np.histogram2d(
        nrpsls[pslsA],
        nrpsls[pslsB],
        density=True, bins=bins)
    X, Y = np.meshgrid(xedges, yedges)
    H = H.T  # Get H prepared for plotting

    # Search for the right number of overdensities
    (res, results, locs) = search_overdensities(
        H, overdensities, penalization, verbose=verbose, debug=debug)

    vprint('>> Optimization results <<')
    vprint(res)
    iters_df = pd.DataFrame.from_records(
        results, index=range(1, len(results) + 1))
    if verbose:
        print(iters_df)

    # Plot selection of overdensities
    axlim: List[float] = [
        min(dset[pslsA]),
        max(dset[pslsA]),
        min(dset[pslsB]),
        max(dset[pslsB])]
    ax = plot_overdensities(X, Y, H, locs, protA, protB, axlim)

    # Retrieve variants from detected overdensities
    vprint('>> Results for overdensities <<')
    num_features: int = len(locs)
    ods: SampleSet = SampleSet({})
    for f in range(num_features):
        # Filter variants in the detected overdensity region (rectangle)
        sliceX = locs[f][1]
        sliceY = locs[f][0]
        pslsA_min = xedges[sliceX.start]
        pslsA_max = xedges[sliceX.stop]
        pslsB_min = yedges[sliceY.start]
        pslsB_max = yedges[sliceY.stop]
        od = dset.query(
            f'{pslsA_min} < {pslsA} < {pslsA_max} and {pslsB_min} < {pslsB} < {pslsB_max}')
        od_psls = od[[pslsA, pslsB]]
        # Estimate (local) max density point within overdensity region as
        # center of top bin in 2D histogram
        od_H, od_xedges, od_yedges = np.histogram2d(
            od_psls[pslsA], od_psls[pslsB], bins=bins)
        inds = np.unravel_index(np.argmax(od_H, axis=None), od_H.shape)
        # x-coord of estimated local max density
        od_pslsA = (od_xedges[inds[0]] + od_xedges[inds[0] + 1]) / 2
        # y-coord of estimated local max density
        od_pslsB = (od_yedges[inds[1]] + od_yedges[inds[1] + 1]) / 2
        # Sort variants by distance to local max density point
        od_cdist = scipy.spatial.distance.cdist(
            od_psls, np.array([[od_pslsA, od_pslsB]]))
        pd.options.mode.chained_assignment = None
        od.loc[:,'cdist'] = od_cdist
        print(
            f"Number of different CAMEOX solutions for overdensity {f+1} circa ({int(np.rint(od_pslsA))},{int(np.rint(od_pslsB))}): {len(od)}")
        vprint(
            f'Overdensity {f+1}: Max 2Dhist is {int(od_H.max())}@{inds} with od_pslsA={od_pslsA:.2f} and od_pslsB={od_pslsB:.2f}')
        od.sort_values(
            'cdist',
            ascending=True,
            inplace=True,
            ignore_index=True)
        pd.options.mode.chained_assignment = 'warn'
        ods[f] = od
        # Plot top sampling points within overdensity regions
        od_psls = od[[pslsA, pslsB]].iloc[:int(len(od) / 100)]
        plt.scatter(od_psls[pslsA], od_psls[pslsB], marker='.')
        plt.plot([pslsA_min, pslsA_min, pslsA_max, pslsA_max, pslsA_min],
                 [pslsB_min, pslsB_max, pslsB_max, pslsB_min, pslsB_min],
                 color="red")

    # Show final plot
    ax = plt.gca()
    ax.set_aspect('equal')
    # plt.axis(axlim)
    plt.title(
        "Overdensities in the APLL space\nand 1% of top variants within each of them")
    plt.xlabel(f'APLL {protA} (rel. to WT)')
    plt.ylabel(f'APLL {protB} (rel. to WT)')
    plt.show()

    # Retrieve variants from detected overdensities
    vprint('>> Sampling from detected overdensities sets...')
    all_sampled_full_seqs: Set[str] = set()
    if sampled is not None:
        all_sampled: pd.DataFrame = pd.concat(
            [df for df in sampled.values()],
            ignore_index=True)
        all_sampled_full_seqs = set(all_sampled['full_seq'])
        print(
            f'Sampling while checking against {len(all_sampled_full_seqs)} already sampled variants...')
    ods_iterrows: Dict[int, Iterator] = {  # Dict with initialized itertuples for the overdensities sets
        f: ods[f].itertuples() for f in range(num_features)}
    od_num_rows: Dict[int, List[int]] = {  # Dict with overdensity map to list of selected columns
        f: [] for f in range(num_features)}
    cnt: typing.Counter[int] = Counter(  # Keeps track of the number of checked variants per overdensity region
        {f: 0 for f in range(num_features)})
    last_cnt_total: int = -1  # Sum of the counts in the previous iteration

    def get_num_sampled(od_num_rows: Dict[int, List[int]] = od_num_rows,
                        num_features: int = num_features,
                        ) -> int:
        """Get number of sampled variants"""
        return np.array([len(od_num_rows[f])
                        for f in range(num_features)], dtype=int).sum()

    # Sampling main loop: get variants from overdensities in a round-robin
    # fashion
    while get_num_sampled() < size and last_cnt_total < sum(
            cnt.values()):  # 2nd condition checks for all itertuples exhausted
        last_cnt_total = sum(cnt.values())
        for f in range(num_features):
            row: NamedTuple
            try:
                row = next(ods_iterrows[f])
            except StopIteration:  # Overdensity set exhausted, jump to next overdensity set
                vprint(
                    f'  Overdensity {f+1} itertuple exhausted! Continuing loop...')
                continue
            if row.full_seq not in all_sampled_full_seqs:  # The set will be empty if sampled is None
                od_num_rows[f].append(cnt[f])
            else:
                vprint(
                    f'  Found already sampled var for overdensity {f+1} var {cnt[f]+1}. Skipping...')
            cnt[f] += 1
            if get_num_sampled() >= size:
                break

    print(
        f'Sampled {get_num_sampled()} unique redundant variants from a total of {sum(cnt.values())}:')
    print(
        ' '.join(
            [f'{len(od_num_rows[f])} from OD{f+1};'
             for f in range(num_features)]))

    # Build and return the joint dataframe
    sampling: pd.DataFrame
    sampling = pd.concat([ods[f].take(od_num_rows[f])
                          for f in range(num_features)],
                         ignore_index=True)
    return sampling


# ### Sampling criterium: random (control)
def sample_random(dset: pd.DataFrame,
                  protA: str,
                  protB: str,
                  size: int = 500,
                  erp: float = None,
                  sampled: SampleSet = None,
                  seed: int = None,
                  verbose: bool = False) -> pd.DataFrame:
    """Randomly sample variants
    """

    def vprint(*arguments, **kargs) -> None:
        """Print only if verbose mode is enabled"""
        if verbose:
            print(*arguments, **kargs)

    if erp is not None:  # Restrict by erp if given
        dset = dset[dset['ERP'] == erp]
    print("Number of different redundant CAMEOX solutions for this "
          f"{'entanglement' if erp is None else 'ERP'}: {len(dset)}")

    rng = np.random.default_rng(seed=seed)
    rnd_idx: List[int] = list(
        rng.choice(
            len(dset),
            size=len(dset),
            replace=False,
            shuffle=True))
    num_rows: List[int] = []

    # Retrieve variants from detected overdensities
    vprint('>> Randomly sampling...')
    all_sampled_full_seqs: Set[str] = set()
    if sampled is not None:
        all_sampled: pd.DataFrame = pd.concat(
            [df for df in sampled.values()],
            ignore_index=True)
        all_sampled_full_seqs = set(all_sampled['full_seq'])
        print(
            f'Sampling while checking against {len(all_sampled_full_seqs)} already sampled variants...')

    # Sampling main loop
    while len(num_rows) < size:
        try:
            raw_idx = rnd_idx.pop()
        # Random index list exhausted (may happen if size large and many
        # already sampled)
        except IndexError:
            print('\n WARNING: Random index list exhausted! Finishing...')
            break
        full_seq: str = dset.iloc[raw_idx]['full_seq']
        if full_seq not in all_sampled_full_seqs:  # The set would be initially empty if sampled is None
            num_rows.append(raw_idx)
            vprint(raw_idx, end=', ')
            # dset is normally nonredundant, but just in case...
            all_sampled_full_seqs.add(full_seq)
        else:
            print(
                f'  Found already sampled variant in row number {raw_idx}. Skipping...')

    print(f'\nSampled {len(num_rows)} unique random variants.')
    return dset.take(num_rows)


# ### Save sampled variants
def save_sampled(sampled: SampleSet,
                 protA: str,
                 protB: str) -> None:
    """Save sampled variants"""

    total: int = 0
    for name, dset in sampled.items():
        print(
            f"Saving sampled dataset '{name}' with {len(dset)} variants... ",
            end='')
        dset[['full_seq', f'{protA}_seq', f'{protB}_seq',
         f'{protA}_psls', f'{protB}_psls','ERP']].to_csv(
            f'{name}_CAMEOX_{protA}_{protB}_pairs.csv', index=False, header=True)
        total += len(dset)
        print('OK!')
    print(
        f"Saved {total} total different variants distributed in {len(sampled)} datasets")
