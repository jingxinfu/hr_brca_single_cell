from plotnine import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Union,Tuple,List,Dict
import matplotlib
import seaborn as sns
def plot_proportion_barplot(
    adata,
    yaxis,
    fill,
    fill_breakdown=None,
    yaxis_label=None,
    y_label_func=None,
    xaxis_label='Proportions',
    fill_label=None,
    percent_limit=5.,
    show_percent=True,
    height_scale=1.,
    width_scale=1.,
    legend_position=None,
    return_df=False,
    normalize_by=None,
    format_x_as_percent=True,
    remove_x_axis_ticks=False,
    swap_axes=False,
    external_percent_label=None,
    bar_position='fill',
    fill_groups=None,
    ):

    import mizani
    import matplotlib.patheffects as pe

    if yaxis_label is None:
        yaxis_label = yaxis
    if fill_label is None:
        fill_label = fill

    adata._sanitize()

    fill_dict = {k: v for k, v in zip(
        adata.obs[fill].cat.categories, adata.uns[f'{fill}_colors'])}

    groupby = [yaxis, fill] + ([external_percent_label]
                               if external_percent_label is not None else [])
    if fill_breakdown is not None and (fill_breakdown not in groupby):
        groupby.append(fill_breakdown)
    if normalize_by is not None and (normalize_by not in groupby):
        print('normalizing by a factor that is not plotted, please make sure you want to do that.')
        groupby.append(normalize_by)

    df_level1 = pd.DataFrame(adata.obs.groupby(
        groupby, observed=True).size(), columns=['counts'])

    if normalize_by:
        scales = df_level1.reset_index().groupby(
            normalize_by)[['counts']].sum()
        scales = scales.sum().div(scales)

        df_level1 = df_level1.multiply(scales)

    df_level0 = df_level1.reset_index().groupby(yaxis)[['counts']].sum()
    df = df_level1.div(df_level0, level=yaxis).reset_index()

    df[fill] = pd.Categorical(
        df[fill], categories=reversed(adata.obs[fill].cat.categories))
    if swap_axes:  # do not reverse order if swap_axes is given
        def rev_func(x): return x
    else:
        rev_func = reversed
    df[yaxis] = pd.Categorical(
        df[yaxis], categories=rev_func(adata.obs[yaxis].cat.categories))

    df['counts_coarse'] = df.groupby([yaxis, fill], observed=True)[
        'counts'].transform('sum')
    df['counts_coarse_round_percent'] = (
        df.counts_coarse*100).round().astype(int)

    df['_show_text'] = df.counts_coarse_round_percent >= percent_limit
    df['_show_breakdown'] = (
        df.counts_coarse_round_percent >= percent_limit) if fill_breakdown else False

    # collapse breakdown of small groups
    if fill_breakdown:
        df = df[(~df.duplicated([yaxis, fill])) | (df._show_breakdown)].copy()
        df.loc[~df._show_breakdown,
               'counts'] = df.loc[~df._show_breakdown, 'counts_coarse']
        df['_show_breakdown'] = True

    cs = df.sort_values([yaxis, fill], ascending=False).drop_duplicates(
        [yaxis, fill]).groupby(yaxis, observed=True)['counts_coarse'].transform(pd.Series.cumsum)
    df['cumsum_mean'] = cs - df['counts_coarse'] + (df['counts_coarse']/2)

    figure_width = 8*width_scale
    figure_height = 0.4*df[yaxis].nunique()*height_scale

    if swap_axes:
        figure_width, figure_height = figure_height, figure_width
        
    data=df[df[fill].isin(fill_groups)] if fill_groups is not None else df
    g = (
        ggplot(aes(x=yaxis, y='counts', fill=fill, group=fill), data) +
        geom_bar(position=bar_position, stat='identity', mapping=aes(color='_show_breakdown'), size=0.08) +
        (scale_y_continuous(labels=mizani.formatters.percent) if format_x_as_percent else geom_blank()) +
        (scale_x_discrete(labels=y_label_func) if y_label_func is not None else geom_blank()) +
        (coord_flip() if not swap_axes else geom_blank()) +
        theme_minimal() +
        theme(
            text=element_text(),
            figure_size=(figure_width, figure_height),
            legend_position=legend_position,
            axis_text_x=element_blank() if remove_x_axis_ticks else element_text(
                angle=90 if swap_axes else 0),
            axis_ticks_major_x=element_blank() if remove_x_axis_ticks else None,
            axis_ticks_minor_x=element_blank() if remove_x_axis_ticks else None,
            panel_grid_major=element_blank(), panel_grid_minor=element_blank(),
        ) +
        scale_color_manual(values={True: 'black', False: 'none'}) +
        scale_fill_manual(values=fill_dict) +
        labs(x=yaxis_label, y=xaxis_label, fill=fill_label) +
        guides(fill=guide_legend(reverse=True), color=None)
    )

    if show_percent:
        if external_percent_label is not None:
            label = external_percent_label
        else:
            label = 'counts_coarse_round_percent'
        g += geom_text(aes(label=label, y='cumsum_mean'), data=df[df._show_text],
                       color='white', size=8, fontweight='bold',
                       path_effects=(pe.Stroke(linewidth=1, foreground='black'), pe.Normal()))

    if return_df:
        return g, df
    else:
        return g

    
################## heatmap
import anndata
import scanpy as sc
from typing import List,Tuple,Union
from PyComplexHeatmap import \
ClusterMapPlotter,HeatmapAnnotation,anno_simple,anno_label,DotClustermapPlotter,anno_scatterplot

def aggregate_heatmap(adata:anndata.AnnData,
                      features:Union[list,pd.DataFrame],
                      palette:dict,
                      top_anno_columns:List[str],
                      output_path:Union[str,None]=None,
                      col_split:Union[str,None]=None,
                      row_split:Union[str,None]=None,
                      standard_normalization:bool=True,
                      figsize: Union[None,Tuple[int]]=(10,10),
                      vmin:float=-3,
                      vmax:float=3,
                      cmap:str='RdBu_r',
                      cbar_label:str='Expression',
                      col_cluster:Union[None,bool]=True,
                      row_cluster:Union[None,bool]=True,
                      delta_col: Union[None,Tuple[str]]=None,
                      delta_order: Union[None,List[str]]=None,
                      **kwargs):
    ## plot
    plt.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size'] = 7
    plt.rcParams['axes.grid'] = False 
    plt.figure(figsize=figsize,dpi=150)


    markers = features if isinstance(features,list) else features.index.tolist()
    top_anno_columns = [delta_col]+top_anno_columns if delta_col is not None else top_anno_columns
    data = sc.get.obs_df(adata,markers+top_anno_columns)
    
    ## calculate the average value in each group
    data = data.groupby(top_anno_columns).mean()
    data.dropna(axis=0,how='all',inplace=True)

    if delta_col is not None:
        ## calculate the mean difference between two condition
        assert len(delta_order)==2, "Please provide a list with two elements to indicate the comparison order"
        delta0_data = data.loc[delta_order[0],:]
        delta1_data = data.loc[delta_order[1],:]
        overlap_rows = delta0_data.index.intersection(delta1_data.index)
        data = delta1_data.loc[overlap_rows,:] - delta0_data.loc[overlap_rows,:]
        del top_anno_columns[top_anno_columns.index(delta_col)]

    ## standard normalization
    if standard_normalization:
        data = data - data.mean()
        data = data / data.std()
    data.reset_index(inplace=True)

    ## top annotation
    col_dict={}
    for col in top_anno_columns:
        colors = { k:v for k,v in palette[col].items() if k in data[col].unique()}
        if col == col_split:
            col_dict['label']=anno_label(data[col],
                                         colors='k',
                                 # colors=colors,
                                 merge=True,rotation=20,size=15)
        col_dict[col] = anno_simple(data[col],
                                    legend_kws=dict(frameon=False),
                                    colors=colors,
                                    legend=False if col == col_split else True)

   
    col_ha = HeatmapAnnotation(**col_dict,verbose=0,axis=1)
    ## left annotation
    row_ha = None
    row_dict={}
    row_orders = markers
    if isinstance(features,pd.DataFrame):
        row_orders = features.sort_values(features.columns.tolist()).index
        for row in features:
            colors = { k:v for k,v in palette[row].items() if k in features[row].unique()}
            if row == row_split:
                row_dict['label']=anno_label(features[row],
                                             colors='k',
                                 # colors=colors,
                                 merge=True,rotation=20,size=15)
            row_dict[row] = anno_simple(features[row],
                                        legend_kws=dict(frameon=False),
                                        legend=False if row == row_split else True,
                                        colors=colors)
        row_ha = HeatmapAnnotation(**row_dict,verbose=0,axis=0,text_kws={'fontsize':28})     
    col_orders = data.sort_values(top_anno_columns).index
    
    cm = ClusterMapPlotter(data=data[markers].T.loc[row_orders,col_orders],
                           top_annotation=col_ha,
                           left_annotation=row_ha,
                           # top annotation
                           col_split=None if col_split is None else data[col_split],
                           col_cluster=col_cluster,
                           # bottom annotation
                           row_names_side='right',
                           show_rownames=True,
                           row_split=row_split if not isinstance(row_split,str)else features[row_split],
                           row_cluster=row_cluster,
                           # row_split_order=None if row_split is None else features[row_split].unique().tolist(), 
                           vmin=vmin,
                           vmax=vmax,
                           label=cbar_label, 
                           cmap=cmap,
                           **kwargs)
    if output_path is not None:
        plt.savefig(output_path, bbox_inches='tight',dpi=200)

def small_sample_box_visual(dt:pd.DataFrame,x:str,y:str,palette:Dict[str,str],ax:Union[None,matplotlib.axes.Axes]=None,**kwargs):
    if ax is None:
        _,ax=plt.subplots(1,1,figsize=(3,6),dpi=200)
    PROPS = {
        'boxprops':{'facecolor':'none', 'edgecolor':'black'},
        'medianprops':{'color':'black'},
        'whiskerprops':{'color':'black'},
        'capprops':{'color':'black'}
    }
    PROPS.update(kwargs)
    sns.stripplot(data=dt,x=x,y=y,ax=ax,palette=palette,**kwargs)
    sns.boxplot(data=dt,x=x,y=y,ax=ax,palette=palette,showfliers=False,linewidth=1,**PROPS)
    xticklabels=[]
    for handle in ax.get_xticklabels():
        text = handle.get_text()
        xticklabels.append(f"{text}\n(N={(dt[x]==text).sum()})")
    ax.set_xticklabels(xticklabels)
    return ax




def check_if_matplotlib(return_mpl=False):
    if not return_mpl:
        try:
            import matplotlib.pyplot as plt
        except Exception:
            raise ImportError('matplotlib is not installed. Please install it with: pip install matplotlib')
        return plt
    else:
        try:
            import matplotlib as mpl
        except Exception:
            raise ImportError('matplotlib is not installed. Please install it with: pip install matplotlib')
        return mpl


def check_if_adjustText():
    try:
        import adjustText as at
    except Exception:
        raise ImportError('adjustText is not installed. Please install it with: pip install adjustText')
    return at

def filter_limits(df, sign_limit=None, lFCs_limit=None):

    # Define limits if not defined
    if sign_limit is None:
        sign_limit = np.inf
    if lFCs_limit is None:
        lFCs_limit = np.inf

    # Filter by absolute value limits
    msk_sign = df['pvals'] < np.abs(sign_limit)
    msk_lFCs = np.abs(df['logFCs']) < np.abs(lFCs_limit)
    df = df.loc[msk_sign & msk_lFCs]

    return df
def plot_volcano_df(data, x, y, top=5, sign_thr=0.05, lFCs_thr=0.5, sign_limit=None, lFCs_limit=None, 
                    show_list = None, color_pos='#D62728',
                    color_neg='#1F77B4', color_null='gray', figsize=(7, 5), dpi=100, ax=None):
    """
    Plot logFC and p-values from a long formated data-frame.

    Parameters
    ----------
    data : pd.DataFrame
        Results of DEA in long format.
    x : str
        Column name of data storing the logFCs.
    y : str
        Columns name of data storing the p-values.
    top : int
        Number of top differentially expressed features to show.
    sign_thr : float
        Significance threshold for p-values.
    lFCs_thr : float
        Significance threshold for logFCs.
    sign_limit : float
        Limit of p-values to plot in -log10.
    lFCs_limit : float
        Limit of logFCs to plot in absolute value.
    color_pos: str
        Color to plot significant positive genes.
    color_neg: str
        Color to plot significant negative genes.
    color_null: str
        Color to plot rest of the genes.
    figsize : tuple
        Figure size.
    dpi : int
        DPI resolution of figure.
    ax : Axes, None
        A matplotlib axes object. If None returns new figure.
    return_fig : bool
        Whether to return a Figure object or not.
    save : str, None
        Path to where to save the plot. Infer the filetype if ending on {``.pdf``, ``.png``, ``.svg``}.

    Returns
    -------
    fig : Figure, None
        If return_fig, returns Figure object.
    """

    # Load plotting packages
    plt = check_if_matplotlib()
    at = check_if_adjustText()

    # Transform sign_thr
    sign_thr = -np.log10(sign_thr)

    # Extract df
    df = data.copy()
    df['logFCs'] = df[x]
    df['pvals'] = -np.log10(df[y])

    # Filter by limits
    df = filter_limits(df, sign_limit=sign_limit, lFCs_limit=lFCs_limit)

    # Define color by up or down regulation and significance
    df['weight'] = color_null
    up_msk = (df['logFCs'] >= lFCs_thr) & (df['pvals'] >= sign_thr)
    dw_msk = (df['logFCs'] <= -lFCs_thr) & (df['pvals'] >= sign_thr)
    df.loc[up_msk, 'weight'] = color_pos
    df.loc[dw_msk, 'weight'] = color_neg

    # Plot
    fig = None
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    df.plot.scatter(x='logFCs', y='pvals', c='weight', sharex=False, ax=ax)
    ax.set_axisbelow(True)

    # Draw sign lines
    ax.axhline(y=sign_thr, linestyle='--', color="black")
    ax.axvline(x=lFCs_thr, linestyle='--', color="black")
    ax.axvline(x=-lFCs_thr, linestyle='--', color="black")

    # Plot top sign features
    if show_list is None:
        signs = df[up_msk | dw_msk].sort_values('pvals', ascending=False)
        signs = signs.iloc[:top]
    else:
        signs = df.loc[df.index.intersection(show_list),:].sort_values('pvals', ascending=False)

    # Add labels
    ax.set_ylabel('-log10(pvals)')
    texts = []
    for x, y, s in zip(signs['logFCs'], signs['pvals'], signs.index):
        texts.append(ax.text(x, y, s))
    if len(texts) > 0:
        at.adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'), ax=ax)


    return fig,ax