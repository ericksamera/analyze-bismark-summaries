#!/usr/bin/env python3
__description__ =\
"""
Purpose: Script generates a box plot of methylation from a list of files that are produced from `process-bismark-coverage.py`
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "Stable enough."
# -----------------------------------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
from pathlib import Path
# -----------------------------------------------------------------------------
from collections import Counter
import pandas as pd
import numpy as np
from skbio import DistanceMatrix, stats
from skbio.tree import nj
import plotly.graph_objects as go
import plotly.express as px
# -----------------------------------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=f"{__description__}",
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter,
        allow_abbrev=False)

    parser.add_argument('input_path',
        metavar="PATH",
        type=Path,
        nargs='+',
        help=f"path of input file")
    parser.add_argument('-o', '--out-dir', dest='output_dir',
        metavar="PATH",
        type=Path,
        default=Path.cwd(),
        help=f"directory for output")

    group_csv_filtering = parser.add_argument_group('csv filtering options')
    group_csv_filtering.add_argument('-m', '--min-reads',
        metavar="INT",
        type=int,
        default=1000,
        help=f"minumum coverage or read depth to consider CpG position. (DEFAULT=1000)")
    group_csv_filtering.add_argument('-r', '--min-pos-rep',
        metavar="FLOAT",
        type=float,
        default=1,
        help=f"minimum representation at CpG position relative to samples in primer (DEFAULT=0.25)")

    group_primer_mut_ex = parser.add_argument_group('mutually exclusive primer options')
    mut_ex_primers = group_primer_mut_ex.add_mutually_exclusive_group()
    mut_ex_primers.add_argument('--include-primer', dest='primer_inclusions',
        metavar="STR",
        type=str,
        help=f"include primer(s) (ex: 'RTL1-BN1X;SNRPN-BP1X')")
    mut_ex_primers.add_argument('--exclude-primer', dest='primer_exclusions',
        metavar="STR",
        type=str,
        help=f"exclude primer(s) (ex: 'RTL1-BN1X;SNRPN-BP1X')")

    group_sample_mut_ex = parser.add_argument_group('mutually exclusive sample options')
    mut_ex_samples = group_sample_mut_ex.add_mutually_exclusive_group()
    mut_ex_samples.add_argument('--include-sample', dest='sample_inclusions',
        metavar="STR",
        type=str,
        help=f"include primer(s) (ex: 'BSX-AMP-BS-05-Monarch;BSX-GC-21222')")
    mut_ex_samples.add_argument('--exclude-sample', dest='sample_exclusions',
        metavar="STR",
        type=str,
        help=f"exclude primer(s) (ex: 'BSX-AMP-BS-05-Monarch;BSX-GC-21222')")

    group_export = parser.add_argument_group('export options')
    group_export.add_argument('--show',
        action='store_true',
        help=f"show HTML visualization")
    group_export.add_argument('--export-by-primer',
        action='store_true',
        help=f"export HTML visualization by primer (currently DEFUNCT)")

    args = parser.parse_args()
    # parser errors and processing
    # -------------------------------------------------------------------------
    args.input_path = [file for file in args.input_path if file.is_file()]
    if not args.input_path: parser.error('Invalid file specified!')

    args.output_dir.mkdir(exist_ok=True)

    return args
# -----------------------------------------------------------------------------
def _filename_filter(_paths_list: list, _filter_index: int, _criteria: str, _inclusion: bool) -> list:
    """
    Modular function to filter based on some filename critera.
    """
    if _inclusion: return [path for path in _paths_list if path.stem.split('_')[_filter_index] in _criteria.split(';')]
    else: return [path for path in _paths_list if path.stem.split('_')[_filter_index] not in _criteria.split(';')]
def _apply_args(args: Namespace) -> Namespace:
    """
    Function to lazily apply this filtering so that this process is less messy.
    """
    if args.primer_exclusions or args.primer_inclusions:
        inclusion_var = True if args.primer_inclusions else False
        criteria_var = args.primer_inclusions if args.primer_inclusions else args.primer_exclusions
        filter_index=-2
        args.input_path: list = _filename_filter(
            _paths_list=args.input_path,
            _filter_index=filter_index,
            _criteria=criteria_var,
            _inclusion=inclusion_var)
    if args.sample_exclusions or args.sample_inclusions:
        inclusion_var = True if args.sample_inclusions else False
        criteria_var = args.sample_inclusions if args.sample_inclusions else args.sample_exclusions
        filter_index=0
        args.input_path: list = _filename_filter(
            _paths_list=args.input_path,
            _filter_index=filter_index,
            _criteria=criteria_var,
            _inclusion=inclusion_var)
    
    return args
# -----------------------------------------------------------------------------
def _filter_positions_by_rep(_input_dict: dict, args: Namespace) -> dict:
    """
    """
    output_dict: dict = {}
    for primer, samples_dict in _input_dict.items():
        output_dict[primer] = {}

        position_represenation = []

        for sample in samples_dict.keys():
            position_represenation += list(samples_dict[sample].keys())
        
        position_counts = Counter(position_represenation)
        max_rep = max(position_counts.values())
        allowed_positions = [key for key, value in position_counts.items() if value/max_rep >= args.min_pos_rep]

        for sample, position_dict in samples_dict.items():
            output_dict[primer][sample] = {}
            for position in position_dict:
                if position not in allowed_positions: continue
                output_dict[primer][sample][position] = _input_dict[primer][sample][position]

    return output_dict
def _get_primers(_input_csv_paths: list) -> list:
    """
    Function to get the number of primers from the list of input CSV files.
    """
    return list(set([primer.stem.split('_')[1] for primer in _input_csv_paths]))
def _generate_dicts(_input_csv_paths: list, args: Namespace) -> tuple:
    
    _methylation_dict: dict = {}
    _coverage_dict: dict = {}

    for csv_path in _input_csv_paths:
        sample_id, primer, _ = csv_path.stem.split('_')

        if primer not in _methylation_dict:
            _methylation_dict[primer] = {}
            _coverage_dict[primer] = {}
        if sample_id not in _methylation_dict[primer]: 
            _methylation_dict[primer][sample_id] = {}
            _coverage_dict[primer][sample_id] = {}

        with open(csv_path, encoding='UTF-8') as input_file:
            headers: list = [header_val for header_val in input_file.readline().strip().split(',')]
            for line in input_file.readlines():
                line_dict: dict = {key: value for key, value in zip(headers, line.strip().split(','))}
                
                coverage = int(line_dict['count_methylated']) + int(line_dict['count_unmethylated'])
                if not coverage > args.min_reads: continue
                position_key: str = f"{line_dict['chromosome']}:{line_dict['start']}"
                _methylation_dict[primer][sample_id][position_key]: float = float(line_dict['methylation_percentage'])
                _coverage_dict[primer][sample_id][position_key]: int = int(coverage)
    return _filter_positions_by_rep(_methylation_dict, args), _filter_positions_by_rep(_coverage_dict, args)
def _transform_dict(_input_primer_dict: dict, add_reference_vals: bool=True) -> dict:
    """
    """
    _reference_vals = {
        'methylated': 100,
        'hemi-methylated': 50,
        'un-methylated': 0
    }

    transformed_dict: dict = {}
    for sample, positions_dict in _input_primer_dict.items():
        for position, value in positions_dict.items():
            if position not in transformed_dict: transformed_dict[position] = {}
            transformed_dict[position][sample] = value
    
    if add_reference_vals:
        for sample in _reference_vals:
            for position, position_dict in transformed_dict.items():
                transformed_dict[position][sample] = _reference_vals[sample]

    return transformed_dict
def _generate_distance_matrix(_input_dict: dict) -> list:
    """
    """
    matrix = []
    for key1, value1 in _input_dict.items():
        distances = []
        for key2, value2 in _input_dict.items():
            cohens_h = abs(2 * ( np.arcsin((value1/100) ** (1/2)) - np.arcsin((value2/100) ** (1/2))))
            distances.append(cohens_h)
        matrix.append(distances)
    return matrix
def _generate_matrices_dict(_input_dict: dict, add_reference_vals: bool=True) -> tuple:
    """
    """
    matrices = {}
    samples = {}
    for primer in _input_dict:
        matrices[primer], samples[primer] = {}, {}
        positions_dict: dict = _transform_dict(_input_dict[primer])
        samples_list = [sample_name for sample_name in _input_dict[primer] if _input_dict[primer][sample_name]]
        if len(samples_list) < 2: continue
        samples[primer] = samples_list
        if add_reference_vals:
            samples[primer] += ['methylated', 'hemi-methylated', 'un-methylated']
        for position in positions_dict:
            matrices[primer][position] =_generate_distance_matrix(positions_dict[position])
    return matrices, samples
def _sum_matrices_dict(_input_matrices: list, _sample_num: int) -> np.array:
    """
    Function sums the input and returns a summed aray.
    """
    initialized_matrix: list = np.array([[0*_sample_num]*_sample_num])
    for matrix in _input_matrices:
        initialized_matrix = np.add(initialized_matrix, np.array(matrix))
    return initialized_matrix
def _generate_pcoa_dict(_matrices_dict: dict, _samples_dict: dict) -> tuple:
    """
    """
    per_primer_pcoa_dict = {}
    for primer, matrix_dicts in _matrices_dict.items():
        if not list(matrix_dicts.values()): continue
        summated_distance_matrix = _sum_matrices_dict(matrix_dicts.values(), len(_samples_dict[primer]))
        distance_matrix_object = DistanceMatrix(summated_distance_matrix, _samples_dict[primer])
        pcoa_result = stats.ordination.pcoa(distance_matrix_object)
        per_primer_pcoa_dict[primer] = pcoa_result
    return per_primer_pcoa_dict
def _generate_truth_list(_offset_dict: dict, max_value: int) -> list:
    """
    """
    if not _offset_dict['start'] > 1: truth_list = []
    else: truth_list = [False]*_offset_dict['start']
    for i in range(_offset_dict['end']-_offset_dict['start']):
        truth_list += [True]
    while len(truth_list) < max_value:
        truth_list += [False]    
    return(truth_list)
def _plot_PCoA(_input_dict: dict, _metadata_dict: dict, _output_dir: Path, _output_prefix: str='', _export: bool=True, _show: bool=False) -> None:
    """
    """

    _reference_vals = ['methylated', 'hemi-methylated', 'un-methylated']

    colors = []
    for i in range(5): colors += list(px.colors.qualitative.Bold)

    _plot_offset = {}
    _list_of_plots = []

    overall_counter = 0
    for i, primer in enumerate(_input_dict):
        if i==0: _plot_offset[primer] = {'start': 0}
        elif primer not in _plot_offset: _plot_offset[primer] = {'start': overall_counter}
        _pca_by_sample: dict = {}
        colors_dict = {key: value for key, value in zip(_metadata_dict[primer], colors)}
        for sample, x, y, z in zip(_input_dict[primer].samples.index, list(_input_dict[primer].samples['PC1']), list(_input_dict[primer].samples['PC2']), list(_input_dict[primer].samples['PC3'])):
            overall_counter += 1
            if sample not in _pca_by_sample: 
                _pca_by_sample[sample] = {}
                _pca_by_sample[sample]['x'], _pca_by_sample[sample]['y'], _pca_by_sample[sample]['z'] = [], [], []
            _pca_by_sample[sample]['x'].append(x)
            _pca_by_sample[sample]['y'].append(y)
            _pca_by_sample[sample]['z'].append(z)
        _plot_offset[primer]['end'] = overall_counter
        for sample in _pca_by_sample:
            markerstyle = 'circle' if sample not in _reference_vals else 'x'
            markersize = 10 if sample not in _reference_vals else 5
            _list_of_plots.append(go.Scatter3d(
                x=_pca_by_sample[sample]['x'],
                y=_pca_by_sample[sample]['y'],
                z=_pca_by_sample[sample]['z'],
                mode='markers',
                name=sample,
                visible=False,
                marker=dict(
                    size=markersize,
                    symbol=markerstyle,
                    color=colors_dict[sample])))

    fig = go.Figure(_list_of_plots)

    buttons_list = []
    for primer in _plot_offset:
        buttons_list.append(
            dict(
                label = primer,
                method = 'update',
                args = [
                    {'visible': _generate_truth_list(_plot_offset[primer], overall_counter)},
                    {'title':
                        f'({primer}) PCoA across all CpG positions<br>'+
                        f'<sup>v{__version__} : {__author__} | {__comments__}</sup>',
                    'scene.xaxis': {'title': f'PCA1 ({_input_dict[primer].proportion_explained[0]*100:.2f} %)'},
                    'scene.yaxis': {'title': f'PCA2 ({_input_dict[primer].proportion_explained[1]*100:.2f} %)'},
                    'scene.zaxis': {'title': f'PCA3 ({_input_dict[primer].proportion_explained[3]*100:.2f} %)'},
                    'showlegend':True}]))
    

    fig.update_layout(
        updatemenus=[go.layout.Updatemenu(
            active=0,
            buttons=buttons_list)])

    camera = dict(eye=dict(x=0., y=0., z=1))
    fig.update_layout(
        title_text=f'PCoA of samples across all CpG positions<br><sup>v{__version__} : {__author__} | {__comments__}   </sup>',
        autosize=True, scene_camera=camera)
    if _show: fig.show()
    if _export: fig.write_html(_output_dir.joinpath(f'summary-methylation-3D-PCoA.html'))
# -----------------------------------------------------------------------------
def pcoa_wrapper(args: Namespace) -> None:
    """
    """
    methylation_dict, metadata_dict = _generate_dicts(args.input_path, args)
    matrices_dict, samples_dict = _generate_matrices_dict(methylation_dict)
    pcoa_dict = _generate_pcoa_dict(matrices_dict, samples_dict)
    _plot_PCoA(
        _input_dict=pcoa_dict, 
        _metadata_dict=samples_dict,
        _output_dir=args.output_dir)
# -----------------------------------------------------------------------------
def main() -> None:
    """ Insert docstring here """
    args = _apply_args(get_args())
    pcoa_wrapper(args)
    return None
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()