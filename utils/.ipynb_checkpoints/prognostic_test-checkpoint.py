import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.fitters.coxph_fitter import CoxPHFitter

def prognostic_test(df:pd.DataFrame,test_var:str,duration_col='xpfs_mos_s1',event_col='xpfs_ind_s1',confounders:list=[])->pd.DataFrame:
    """
    Test the prognostic association for single variable with or without ajusting confounders

    Parameters
    ----------
    df : pd.DataFrame
        patient clinical information data table
    test_var : str
        column name in ``df`` to test
    duration_col : str, optional
        duration column, by default 'xpfs_mos_s1'
    event_col : str, optional
        event column, by default 'xpfs_ind_s1'
    covariates: list, optional
        using multiva
    Returns
    -------
    pd.DataFrame
       testing result
    """
    dt = df[[duration_col,event_col,test_var]+confounders].copy()
    for x in confounders:
        if dt[x].dtype =='O' or dt[x].dtype.name =='category':
            dt[x],_ = pd.factorize(dt[x])

    if dt[test_var].dtype =='O':
        dt[test_var] = dt[test_var].str.replace('\n',' ')
        categories = dt[test_var].unique()
        dt = pd.get_dummies(dt,columns=[test_var])
        result = []
        for category in categories:
            sub_test_var = test_var+'_'+category
            # if category in ['Unknown','Not\nAssessed']:
            #     continue
            tmp_dt = dt[[duration_col,event_col,sub_test_var]+confounders]
            cph = CoxPHFitter(alpha=.1)
            res = cph.fit(tmp_dt,duration_col=duration_col,event_col = event_col).summary.loc[sub_test_var]
            res['label'] = category
            res['group'] = test_var
            res['N'] = (tmp_dt[sub_test_var]==1).sum()
            if len(confounders)==0: # record the estimated median survival time 
                kmf = kmf = KaplanMeierFitter()
                kmf.fit(tmp_dt.loc[tmp_dt[sub_test_var]==1,duration_col],
                        tmp_dt.loc[tmp_dt[sub_test_var]==1,event_col],
                        alpha=.1)
                res['Median PFS'] = f"{kmf.median_survival_time_:.2}"
            result.append(res)
        result = pd.concat(result,axis=1).T

    else:
        cph = CoxPHFitter(alpha=.1)
        cph.fit(dt,duration_col=duration_col,event_col = event_col)
        result = cph.summary.loc[test_var].to_frame().T
        result['label'] = 'continuous'
        result['group'] = test_var
        result['N'] = dt.shape[0]
        if len(confounders)==0: # record the estimated median survival time 
            result['Median PFS'] = ''
    result['est_and_ci'] = result[['exp(coef) lower 90%','exp(coef) upper 90%','exp(coef)']].apply(lambda r:f"{r[2]:.2f}[{r[0]:.2f},{r[1]:.2f}]",axis=1)
    
    return result