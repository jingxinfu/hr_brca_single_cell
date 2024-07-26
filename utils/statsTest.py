import pandas as pd
import numpy as np
from typing import List
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

def lmTest(dt:pd.DataFrame,Y:List[str],x:str,random_var:str,Z:List[str],explored_z:str=None)->pd.DataFrame:
    result =[]
    l0,l1 = dt[x].cat.categories
    if explored_z is not None and explored_z not in Z:
        Z.append(explored_z)
    data = dt[Y+[x]+Z+[random_var]].copy()
    for y in Y:
        reg = smf.mixedlm(f"{y}~0+{x}+{'+'.join(Z)}".strip('+'),data=data,groups=data[random_var]).fit()
        coef = reg.params[f'{x}[{l1}]'] - reg.params[f'{x}[{l0}]']
        logfc =  np.log2( data.loc[data[x]==l1,:].groupby(random_var)[y].mean().mean() / data.loc[data[x]==l0,:].groupby(random_var)[y].mean().mean())
        ## set params other than compared two coefs to 0
        ## since there is an additional params 'random_var', the we should subtract 3 other than 2
        # pvalue = reg.t_test(np.array([[-1,1]+[0]*(reg.params.size-3)])).pvalue
        pvalue = reg.wald_test(f"{x}[{l1}]={x}[{l0}]").pvalue
        if explored_z is None:
            result.append(pd.Series({
                 "Ctrl":l0,
            "Experiment":l1,
            "Y":y,
            "Coef":coef,
            "log2FC":logfc,
            "Pvalue":pvalue
            }).to_frame().T)
        else:
            result.append(pd.Series({
                "Ctrl":l0,
            "Experiment":l1,
            "Y":y,
            "Coef":coef,
            "log2FC":logfc,
            "Pvalue":pvalue,
                explored_z:'+'.join(data[explored_z].unique().tolist())
            }).to_frame().T)
            Z_minus = [z for z in Z if z !=explored_z]
            for gv,subdata in data.groupby(explored_z):
                reg = smf.mixedlm(f"{y}~0+{x}+{'+'.join(Z_minus)}".strip('+'),data=subdata,groups=subdata[random_var]).fit()
                coef = reg.params[f'{x}[{l1}]'] - reg.params[f'{x}[{l0}]']
                logfc =  np.log2( data.loc[data[x]==l1,y].mean() / data.loc[data[x]==l0,y].mean())
                ## set params other than compared two coefs to 0
                ## since there are additional params 'random_var', the we should subtract 3 other than 2
                pvalue = reg.wald_test(f"{x}[{l1}]={x}[{l0}]").pvalue
                result.append(pd.Series({
                    "Ctrl":l0,
            "Experiment":l1,
                    "Y":y,
                    "Coef":coef,
                    "log2FC":logfc,
                    "Pvalue":pvalue,
                    explored_z:gv
                }).to_frame().T)
    
    result = pd.concat(result,axis=0).sort_values('log2FC')
    ## FDR Pvalue adjust
    if explored_z is None:
        result['FDR'] =  multipletests(result['Pvalue'].values,method='fdr_bh')[1]
    else:
        for gv,_ in result.groupby(explored_z):
            result.loc[result[explored_z]==gv,'FDR'] =  multipletests(result.loc[result[explored_z]==gv,'Pvalue'].values,
                                                                      method='fdr_bh')[1]
    return result