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
import plotly.graph_objects as go
from plotly.subplots import make_subplots
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
        default=0.25,
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
    group_export.add_argument('--no-summary-export',
        action='store_true',
        help=f"specifically don't export the summary HTML visualization")
    group_export.add_argument('--export-by-primer',
        action='store_true',
        help=f"export HTML visualization by primer")

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
def _plot_dict(_input_dict: dict, _metadata_dict: dict, _output_dir: Path, _output_prefix: str, _export: bool=True, _show: bool=False) -> None:
    """
    """
    colors = []
    for i in range(5):
        colors += list(px.colors.qualitative.Bold)

    fig = make_subplots(
        rows=len(_input_dict.keys()), 
        cols=1,
        subplot_titles=list(_input_dict.keys()))

    for i, primer in enumerate(_input_dict):
        _input_df = pd.DataFrame.from_dict(_input_dict[primer], orient='index')
        _input_df = _input_df.reindex(columns=sorted(_input_df.columns))

        _metadata_df = pd.DataFrame.from_dict(_metadata_dict[primer], orient='index')
        _metadata_df = _metadata_df.reindex(columns=sorted(_metadata_df.columns))

        headers = _input_df.columns.values
        samples = list(_input_df.index)
        colors_dict = {key: value for key, value in zip(samples, colors)}

        for header in headers:
            methylation_per_header = []
            for sample in _input_df.index:
                methylation_per_header.append(_input_df[header].loc[sample])
            fig.add_trace(go.Box(
                boxpoints='suspectedoutliers',
                y=methylation_per_header,
                name=header,
                pointpos=0,
                width=0.9,
                hoverinfo='skip',
                showlegend=False,
                marker=dict(
                    color='grey',
                    symbol="diamond",
                    size=12)
                ),
                row=i+1, col=1)

        show_in_legend = False if i > 0 else True

        for sample in _input_df.index:
            methylation_percent = list(_input_df.loc[sample])
            fig.add_trace(go.Scatter(
                x=headers,
                y=methylation_percent,
                name=sample,
                mode='markers',
                showlegend=show_in_legend,
                marker=dict(color=colors_dict[sample]),
                hovertemplate = \
                    '<b>%{y}</b>'+
                    '<br>%{text}',
                text=['Depth {:.0f}'.format(i) for i in _metadata_df.loc[sample]],
                ),
                row=i+1, col=1)

    fig.update_yaxes(range=[-10, 110])
    fig.update_layout(
        title_text=f'Boxplot of methylation per primer per CpG<br><sup>v{__version__} : {__author__} | {__comments__}   </sup>', 
        boxmode='group', autosize=True, height=400*len(_input_dict.keys()))
    if _show: fig.show()
    if _export: fig.write_html(_output_dir.joinpath(f'methylation-boxplots-{_output_prefix}.html'))
# -----------------------------------------------------------------------------
def boxplot_wrapper(args: Namespace) -> None:
    """
    """
    methylation_dict, metadata_dict = _generate_dicts(args.input_path, args)
    summary_export_arg = True if not args.no_summary_export else False
    show_arg = True if args.no_summary_export else False
    _plot_dict(
        _input_dict=methylation_dict,
        _metadata_dict=metadata_dict,
        _output_dir=args.output_dir,
        _output_prefix='summary',
        _export=summary_export_arg,
        _show=show_arg)
    if args.export_by_primer:
        for primer in methylation_dict:
            empty_methylation_dict, empty_metadata_dict = {}, {}
            empty_methylation_dict[primer] = methylation_dict[primer]
            empty_metadata_dict[primer] = metadata_dict[primer]
            _plot_dict(
                _input_dict=empty_methylation_dict,
                _metadata_dict=empty_metadata_dict,
                _output_dir=args.output_dir,
                _output_prefix=primer,
                _export=True,
                _show=False)
# -----------------------------------------------------------------------------
def main() -> None:
    """ Insert docstring here """
    args = _apply_args(get_args())
    boxplot_wrapper(args)
    return None
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()