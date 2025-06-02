from typing import Union, List, Tuple,Optional
import pyranges as pr
from pathlib import Path
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pyBigWig

__all__ = ['coveragePlot', 'genebodyPlot','linkagePlot','read_and_filter_gtf','regionPlot']

def create_colormap(values: List[float], cmap: str = 'viridis') -> np.ndarray:
    """
    Create a colormap based on a list of continuous values.

    Parameters
    ----------
    values : List[float]
        A list of continuous values to be mapped to colors.
    cmap : str, optional
        The name of the colormap to use (default is 'viridis').

    Returns
    -------
    np.ndarray
        An array of RGBA color values corresponding to the input continuous values.

    Notes
    -----
    The function normalizes the input values between 0 and 1, then uses the specified colormap to map the values
    to corresponding colors.

    Examples
    --------
    >>> values = [0.1, 0.5, 1.0, 1.5, 2.0]
    >>> colors_mapped = create_colormap(values, cmap='plasma')
    >>> print(colors_mapped)
    array([[0.050383, 0.029803, 0.527975, 1. ],
           [0.983868, 0.904867, 0.136897, 1. ],
           [0.940015, 0.975158, 0.131326, 1. ],
           [0.744136, 0.691368, 0.102124, 1. ],
           [0.323998, 0.229806, 0.506349, 1. ]])
    """
    # Normalize the continuous values between 0 and 1
    norm = matplotlib.colors.Normalize(vmin=min(values), vmax=max(values))

    # Get the colormap (e.g., 'viridis', 'plasma', etc.)
    cmap = plt.get_cmap(cmap)

    # Map normalized values to colors
    color_values = cmap(norm(values))

    return color_values


def getcoverage(
    bw_file: Path, chrom: str, start: int, end: int
) -> Tuple[np.array, np.array]:
    """
    Retrieves coverage data from a BigWig file for a specified genomic region.

    Parameters
    ----------
    bw_file : Path
        The path to the BigWig file.
    chrom : str
        The chromosome of interest.
    start : int
        The start position of the genomic region.
    end : int
        The end position of the genomic region.

    Returns
    -------
    Tuple[np.array, np.array]
        A tuple containing two numpy arrays: positions and their corresponding coverage values.
    """
    with pyBigWig.open(bw_file) as bw:
        y = bw.values(chrom, start, end)
        y = np.nan_to_num(y)
    x = np.array(range(start, end, 1))

    return x, y


def coveragePlot(
    bw_file: Path,
    chrom: str,
    start: int,
    end: int,
    ylim_max:float=3,
    highlight_regions: Union[None, List[Tuple[int, int]]] = None,
    label: str = 'celltype',
    color: str = 'red',
    ax: Union[None,  matplotlib.axes.Axes] = None,
) -> matplotlib.axes.Axes:
    """
    Plots the coverage from a BigWig file over a specified genomic region.

    Parameters
    ----------
    bw_file : Path
        The path to the BigWig file.
    chrom : str
        The chromosome of interest.
    start : int
        The start position of the genomic region.
    end : int
        The end position of the genomic region.
    ylim_max:float
        the maximum value of peaks
    highlight_regions : Union[None, List[Tuple[int, int]]], optional
        A list of regions to highlight, specified as tuples of start and end positions. Default is None.
    label : str, optional
        The label for the plot. Default is 'celltype'.
    color : str, optional
        The color for the coverage plot. Default is 'red'.
    ax : Union[None, matplotlib.axes], optional
        A matplotlib axes object to plot on. If None, a new axes is created.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the plotted coverage.
    """
    assert start < end, 'start > end'
    x, y = getcoverage(bw_file, chrom, start, end)
    if ax == None:
        ax = plt.gca()
    ax.fill_between(x, y1=y, y2=0, step="mid", linewidth=0, color=color)
    ax.patch.set_alpha(0)
    ax.get_yaxis().set_ticks([])
    ax.set_xticks([])
    ax.set_ylabel(label, rotation=0, ha='right')
    ax.spines[['left','right', 'top']].set_visible(False)
    if highlight_regions:
        for x1, x2 in highlight_regions:
            assert (
                x1 >= start and x2 <= end
            ), f'the highlight region ({x1},{x2})is out of plotting region'
            ax.axvspan(x1, x2, alpha=0.05, color='gray')
    ax.set_xlim(start, end)
    ax.set_ylim(0,ylim_max)
    return ax


def read_and_filter_gtf(gtf_file: Path) -> pr.PyRanges:
    """
    Reads a GTF file and filters out rows without gene name to create a PyRanges object.

    Parameters
    ----------
    gtf_file : Path
        The path to the GTF file to be read.

    Returns
    -------
    pr.PyRanges
        A PyRanges object containing the filtered gene annotations.
    """
    gtf = pr.read_gtf(gtf_file)

    # Filter out rows where 'gene_name' is None or empty
    filtered_gtf = gtf[gtf.gene_name.notnull() & (gtf.gene_name != '')]

    return filtered_gtf


def genebodyPlot(
    pr_gtf: pr.PyRanges,
    chrom: str,
    start: int,
    end: int,
    gene_height: int = 1,
    exon_height: int = 2,
    gene_gap_ratio:float=0.5,
    gene_fontsize: int = 15,
    ax: Union[None, matplotlib.axes.Axes] = None,
) -> matplotlib.axes:
    """
    Plots the gene body structure within a specified genomic region.

    Parameters
    ----------
    pr_gtf : pr.PyRanges
        A PyRanges object containing gene annotations.
    chrom : str
        The chromosome of interest.
    start : int
        The start position of the genomic region.
    end : int
        The end position of the genomic region.
    gene_height : int, optional
        The height of the gene representation. Default is 1.
    exon_height : int, optional
        The height of the exon representation. Default is 2.
    gene_gap_ratio: float,optional
        The gap relative to exon height between two genes. Default is 0.5.
    gene_fontsize : int, optional
        The font size for gene labels. Default is 15.
    ax : Union[None, matplotlib.axes], optional
        A matplotlib axes object to plot on. If None, a new axes is created.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the plotted gene body structure.
    """
    pr_region = pr.PyRanges(chromosomes=[chrom], starts=[start], ends=[end])
    # draw the genes of interest, from our gtf
    # intersect genes gtf with the region of interest
    if ax == None:
        ax = plt.gca()
    ax.set_xlim(start, end)
    # intersect gtf with region and get first 9 columns
    gtf_region_intersect = pr_gtf.intersect(pr_region).df
    # only keep exon and gene info
    gtf_region_intersect = gtf_region_intersect.loc[
            gtf_region_intersect.Feature.isin(
                ['gene', 'exon', 'start_codon', 'stop_codon']
            )
            & (gtf_region_intersect.gene_type == 'protein_coding'),
            :,
        ]

    genes_in_window = set(gtf_region_intersect.gene_name)
    n_genes_in_window = len(genes_in_window)

    # set ylim otherwise the plot wouldn't show up
    ax.set_ylim(0, exon_height * n_genes_in_window * (1+gene_gap_ratio))

    for idx, _gene in enumerate(genes_in_window):
        _exon_bottom = idx * (exon_height * (1+gene_gap_ratio))
        _gene_bottom = _exon_bottom + gene_height / 2

        start_point = gtf_region_intersect.loc[
            gtf_region_intersect['gene_name'] == _gene, 'Start'
        ].min()
        end_point = gtf_region_intersect.loc[
            gtf_region_intersect['gene_name'] == _gene, 'Start'
        ].max()
        for _, part in gtf_region_intersect.loc[
            gtf_region_intersect['gene_name'] == _gene
        ].iterrows():
            # make exons thick
            if part['Feature'] == 'exon':
                exon_start = part['Start']
                exon_end = part['End']
                # draw rectangle for exon
                rect = mpatches.Rectangle(
                    (exon_start, _exon_bottom),
                    exon_end - exon_start,
                    exon_height,
                    color='navy',
                )
                ax.add_patch(rect)
                prev_exon_end = exon_end
            # using arrow to note direction
            elif part['Feature'] == 'gene':
                gene_start = part['Start']
                gene_end = part['End']
                x =np.linspace(gene_start,gene_end, 20, endpoint=True)
                y = np.ones(x.shape[0])* (_gene_bottom + gene_height/2)
                ax.plot(x,y, marker='4' if part['Strand'] =='+' else '3',ms=15,color='navy')

            # note the start site
            elif part['Feature'] == 'start_codon':
                gene_start = max(part['Start'] - 1500, start_point)
                gene_end = part['End']
                rect = mpatches.Rectangle(
                    (gene_start, _gene_bottom),
                    gene_end - gene_start,
                    gene_height,
                    fill=True,
                    color="darkgreen",
                )
                ax.add_patch(rect)
            # note the end site
            elif part['Feature'] == 'stop_codon':
                gene_start = part['Start']
                gene_end = min(part['End'] + 1500, end_point)
                rect = mpatches.Rectangle(
                    (gene_start, _gene_bottom),
                    gene_end - gene_start,
                    gene_height,
                    fill=True,
                    color="k",
                )
                ax.add_patch(rect)
                

        ax.text(
            start_point, _exon_bottom + exon_height * 1.1, _gene, fontsize=gene_fontsize
        )
        
    ax.spines[['left', 'right', 'top', 'bottom']].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    return ax


def linkagePlot(
    pr_interact: pr.PyRanges,
    chrom: str,
    start: int,
    end: int,
    label: str = 'Coaccessibility',
    cbar_label:str='score',
    rad:float=.3,
    ax: Union[None, matplotlib.axes.Axes] = None,
    cax:Union[None, matplotlib.axes.Axes] = None,
    cmap: str = 'plasma',
) -> matplotlib.axes.Axes:
    """
    Generate a linkage plot for a specified chromosomal region using PyRanges data.

    Parameters
    ----------
    pr_interact : pr.PyRanges
        A PyRanges object containing interaction data or genomic ranges to plot.
    chrom : str
        The chromosome on which the region is located (e.g., 'chr1').
    start : int
        The start position of the region of interest on the chromosome.
    end : int
        The end position of the region of interest on the chromosome.
    rad : float, optional
        The rad of two point connection (default is .3).
    label : str, optional
        The label for the plot (default is 'Coaccessibility').
    cbar_label : str, optional
        The label for the color bar (default is 'score').
    ax : Union[None, matplotlib.axes.Axes], optional
        A Matplotlib Axes object on which to plot. If None, a new Axes object is created (default is None).
    cax : Union[None, matplotlib.axes.Axes], optional
        A Matplotlib Axes object on which to plot color bar. If None, the color bar will not be drawn (default is None).
    cmap : str, optional
        The colormap to use for the plot (default is 'plasma').

    Returns
    -------
    matplotlib.axes.Axes
        The Matplotlib Axes object containing the plot.

    Notes
    -----
    This function visualizes the interactions or genomic regions provided in `pr_interact` for the specified region
    on the chromosome. If `ax` is None, a new plot is created. The y-axis is scaled according to `ylim_min`, and the
    plot is labeled using the `label` parameter.

    Examples
    --------
    >>> import pyranges as pr
    >>> # Assuming pr_interact is a PyRanges object with interaction data
    >>> linkagePlot(pr_interact, 'chr1', 100000, 200000, ylim_min=-0.5, label='Gene Interaction')
    <matplotlib.axes._subplots.AxesSubplot object at 0x...>
    """
    if ax is None:
        ax = plt.gca()
    pr_region = pr.PyRanges(chromosomes=[chrom], starts=[start], ends=[end])
    region_interact_intersect = pr_interact.intersect(pr_region).df.sort_values('score')

    continuous_values = region_interact_intersect['score'].tolist()
    colors_mapped = create_colormap(continuous_values, cmap=cmap)
    for i, row in region_interact_intersect.iterrows():
        if int(row['targetStart']) > end:
            continue
        posA = (int(row['sourceStart']), 0)
        posB = (int(row['targetStart']), 0)
        # this to ensure arcs are always down
        sign = '-' if posA[0] > posB[0] else ''
        color = colors_mapped[i]
        arrow = mpatches.FancyArrowPatch(
            posA=posA, posB=posB, connectionstyle=f"arc3,rad={sign}{rad}", color=color,
        )
        ax.add_patch(arrow)
    ax.set_xlim(start, end)
    ax.set_ylim(-1, 0)
    ax.spines[['left', 'right', 'bottom']].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.set_ylabel(label, rotation=0, ha='right')
    if cax is not None:
        plt.colorbar(
            plt.cm.ScalarMappable(
                norm=matplotlib.colors.Normalize(
                    vmin=min(continuous_values), vmax=max(continuous_values)
                ),
                cmap=cmap,
            ),
            cax=cax,
            shrink=0.7,
            label=cbar_label,
        )
    else:
        plt.colorbar(
            plt.cm.ScalarMappable(
                norm=matplotlib.colors.Normalize(
                    vmin=min(continuous_values), vmax=max(continuous_values)
                ),
                cmap=cmap,
            ),
            ax=ax,
            pad=0,
            label=cbar_label,
            orientation='horizontal', fraction=.2
            )
    return ax


def regionPlot(
    regions: pr.PyRanges,
    chrom: str,
    start: int,
    end: int,
    highlight_regions:List[Tuple[int]],
    label: str = 'Peaks',
    color: str = 'k',
    region_bed_height: int = 2,
    ax: Optional[matplotlib.axes.Axes] = None
) -> matplotlib.axes.Axes:
    """
    Plots genomic regions on a specified chromosome and position range.

    Parameters:
    -----------
    regions : pr.PyRanges
        A PyRanges object containing genomic regions to plot.
    chrom : str
        Chromosome identifier (e.g., 'chr1', 'chr2').
    start : int
        Start position of the genomic range to plot.
    end : int
        End position of the genomic range to plot.
    label : str, optional
        Label for the plotted regions. Default is 'Peaks'.
    color : str, optional
        Color of the plotted regions. Default is 'k' (black).
    region_bed_height : int, optional
        Height of the region bed in the plot. Default is 2.
    ax : Optional[matplotlib.axes.Axes], optional
        Matplotlib axis on which to plot. If None, a new axis is created.

    Returns:
    --------
    matplotlib.axes.Axes
        The axis with the plotted genomic regions.

    """

    if ax == None:
        ax = plt.gca()
    pr_region = pr.PyRanges(chromosomes=[chrom], starts=[start], ends=[end])
    regions_intersect = regions.intersect(pr_region).df
    ax.set_ylim(0, region_bed_height*2)
    ax.set_xlim(start,end)
    for _, r in regions_intersect.iterrows():
        bed_start = r['Start']
        bed_end = r['End']
        rect = mpatches.Rectangle(
            (bed_start,0.5*region_bed_height), bed_end-bed_start, region_bed_height, fill=True, color=color, linewidth=1)
        ax.add_patch(rect)
    if highlight_regions:
        for x1, x2 in highlight_regions:
            assert (
                x1 >= start and x2 <= end
            ), f'the highlight region ({x1},{x2})is out of plotting region'
            ax.axvspan(x1, x2, alpha=0.05, color='grey')
            
    ax.spines[['left', 'right', 'top', 'bottom']].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel(label, rotation=0, ha='right')
    return ax