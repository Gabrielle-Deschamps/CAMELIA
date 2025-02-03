import numpy as np
import pandas as pd
from sklearn.metrics import *

## Define score function
def generate_scorer(nn_model, data_in, metric_func=None, **kwargs):
    X, X1 = data_in
    def score(X, y):
        n = X.shape[0]
        y_pred = nn_model.predict((X,np.zeros((n,1))))['species']
        
        if metric_func is None:
            return roc_auc_score(y, y_pred, average=None)
        
        else:
            return metric_func(y, y_pred, **kwargs)
    
    return score


def compute_pimp(nn_model, data_in, data_out, splist, tlist, nrepet=10, log_period=10, metric_func=None, **kwargs):
    scorer = generate_scorer(nn_model, data_in, metric_func=metric_func, **kwargs)
    
    nspecies = data_out['species'].shape[1]
    nfeats = data_in[0].shape[1]
    
    ref_scores = scorer(X=data_in[0],y=data_out['species'])
    
    all_repet = []
    for repet in range(nrepet):
        if repet%log_period==0:
            print('Done %d repetitions'%repet)
            
        pimp = []
        for c in range(nspecies,nfeats):
            X_perm = data_in[0].copy()
            np.random.shuffle(X_perm[:,c])

            scores = scorer(X_perm,data_out['species'])-ref_scores
            pimp.append(scores)

        all_repet.append(pimp)
        
        
    pimp_df = pd.DataFrame(data=np.stack([np.stack(arr,axis=1) for arr in all_repet], axis=0).mean(axis=0),
                           columns=['richness']+['cm_%s'%tr for tr in tlist]+['cstd_%s'%tr for tr in tlist],
                           index=splist)
    
    return pimp_df