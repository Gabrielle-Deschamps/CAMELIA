import numpy as np
import pandas as pd

from sklearn.metrics import recall_score, roc_auc_score, average_precision_score, precision_score, f1_score, roc_curve, jaccard_score
from sklearn.metrics import r2_score, mean_squared_error, max_error
from scipy.stats.stats import pearsonr, spearmanr
from utilities.data_utilities import *

import matplotlib.pyplot as plt
import seaborn as sns

def optimize_thres(Y, Y_hat):
    opt_ths = []
    splist = Y.columns.tolist()
    for sp in splist:
        # calculate roc curves
        fpr, tpr, thresholds = roc_curve(Y[sp], Y_hat[sp])
        # get the best threshold
        J = tpr - fpr

        opt_th = thresholds[np.argmax(J)]
        
        opt_ths.append(opt_th)
    
    return opt_ths


def evaluate_regression_model(Y_true,Y_pred, targets, algo='model', do_plot=False):
    perfs = {}
    
    r2_vals = {}
    for sp, score in zip(targets,r2_score(y_true=Y_true,y_pred=Y_pred,multioutput='raw_values')):
        r2_vals[sp]=score
            
    r2_vals['average']=r2_score(y_true=Y_true,y_pred=Y_pred,multioutput='uniform_average')
    r2_vals['weighted']=r2_score(y_true=Y_true,y_pred=Y_pred,multioutput='variance_weighted')
    
    perfs['r2']=r2_vals
    
    mse_vals = {}
    for sp, score in zip(targets,mean_squared_error(y_true=Y_true,y_pred=Y_pred,multioutput='raw_values')):
        mse_vals[sp]=score
    
    mse_vals['average']=mean_squared_error(y_true=Y_true,y_pred=Y_pred,multioutput='uniform_average')
    
    perfs['mse']=mse_vals   
    
    ### Correlation and residuals
    Y_true = Y_true.apply(pd.to_numeric, errors='coerce')
    Y_pred = Y_pred.apply(pd.to_numeric, errors='coerce')
    error = Y_true-Y_pred
    pears_vals = {}
    spear_vals = {}
    pears_pvals = {}
    spear_pvals = {}
    
    for rv in targets:
        pears = pearsonr(Y_true[rv],Y_pred[rv])
        spear = spearmanr(Y_true[rv],Y_pred[rv])
        
        pears_vals[rv]= pears[0]
        spear_vals[rv]= spear[0]
        
        pears_pvals[rv]= pears[1]
        spear_pvals[rv]= spear[1]        
        
        if do_plot:
            fig,ax=plt.subplots(1,2,figsize=(8,4))
            ### Prediction
            sns.scatterplot(x=Y_true[rv],y=Y_pred[rv],ax=ax[0])
            ax[0].axline((0, 0), slope=1)
            ax[0].set_xlabel('True')
            ax[0].set_ylabel('Predicted')
            ax[0].set_title('%s, R2:%.2f'%(rv,r2_vals[rv]))

            ### Residuals
            error[rv].hist(ax=ax[1])
            ax[1].set_title('%s residuals'%rv)

            fig.suptitle('%s'%algo)
        
    
    perfs['pearson']=pears_vals
    perfs['pearson_pv']=pears_pvals
    perfs['spearman']=spear_vals
    perfs['spearman_pv']=spear_pvals
    
    return perfs

def eval_species_perfs(Ys_hat,Ys_true, th_list=0.5):
    aucs = []
    praucs = []
    
    splist = Ys_true.columns.tolist()
    
    for sp in splist:
        try:
            aucs.append(roc_auc_score(y_score=Ys_hat[sp],y_true=Ys_true[sp],average=None))
        except:
            aucs.append(np.nan)
            
        try:
            praucs.append(average_precision_score(y_score=Ys_hat[sp],y_true=Ys_true[sp],average=None))
        except:
            praucs.append(np.nan)            

    sens = recall_score(y_pred=(Ys_hat>th_list),y_true=Ys_true,average=None)
    specs = recall_score(y_pred=1-(Ys_hat>th_list),y_true=1-Ys_true,average=None)
    precs = precision_score(y_pred=1-(Ys_hat>th_list),y_true=1-Ys_true,average=None)
    f1s = f1_score(y_pred=1-(Ys_hat>th_list),y_true=1-Ys_true,average=None)

    support = Ys_true.sum(axis=0)

    sp_perfs = pd.DataFrame(data=np.stack([aucs,praucs,sens,specs,precs,f1s,support],axis=1),
                            columns=['auc','prauc','sensitivity','specificity','precision','f1','support'])
    
    sp_perfs['species']=splist
    sp_perfs['tss']=sp_perfs['sensitivity']+sp_perfs['specificity']-1
    
    return sp_perfs

comm_indices = ['richness','cm','cstd','fskew','fkurt']

def plot_community_pred(Y_test,pred, tlist, test_idx=None, algo='mtnn', plot_residual=False, plot_perfs=False):
    sel_indices = [att for att in comm_indices if (att in Y_test.keys()) & (att in pred.keys())]
    trait_indices = [c for c in sel_indices if c!='richness']
    Y_comm = np.concatenate([pred[att] for att in sel_indices]+[Y_test[att] for att in sel_indices],axis=1)
    
    cols = []
    for att in sel_indices:
        if att=='richness':
            cols.append(att)
        else:
            cols += ['%s_%s'%(tr,att) for tr in tlist]

    obs_cols = ['obs_%s'%c for c in cols]
    pred_cols = ['pred_%s'%c for c in cols]
    Y_comm = pd.DataFrame(data=Y_comm, columns=pred_cols+obs_cols,index=test_idx)
    
    if plot_perfs:
        if 'richness' in sel_indices:
            fig, ax = plt.subplots(1,1)
            Y_comm.plot.scatter(x='obs_richness',y='pred_richness',ax=ax)
            r, p = pearsonr(Y_comm['obs_richness'], Y_comm['pred_richness'])
            sns.regplot(data=Y_comm,x='obs_richness',y='pred_richness',ax=ax)
            ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),transform=ax.transAxes)
            fig.suptitle('Richness') 

        fig, ax = plt.subplots(len(trait_indices),len(tlist),figsize=(10,10))
        for i, ind in enumerate(trait_indices):
            for j, tr in enumerate(tlist):
                r, p = pearsonr(Y_comm['obs_%s_%s'%(tr,ind)], Y_comm['pred_%s_%s'%(tr,ind)])
                sns.regplot(data=Y_comm,x='obs_%s_%s'%(tr,ind),y='pred_%s_%s'%(tr,ind),ax=ax[i,j])
                ax[i,j].text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),transform=ax[i,j].transAxes)

                ax[i,j].set_xlabel('observed')
                ax[i,j].set_ylabel('predicted')
                ax[i,j].set_title('%s - %s'%(ind,tr))

        fig.suptitle('Functional indices prediction')    
        fig.tight_layout()
    
    Y_comm_true = Y_comm[obs_cols]
    Y_comm_hat = Y_comm[pred_cols]
    
    Y_comm_true.columns = Y_comm_hat.columns = cols
    perfs = evaluate_regression_model(Y_comm_true,Y_comm_hat, cols, algo=algo, do_plot=plot_residual)
    perf_df = pd.DataFrame.from_dict(perfs)
    
    return perf_df


def evaluate_prediction(plant_model, data_train, data_test, T, splist, tlist, method='raw', do_residual=False, eps=1e-6):
    data_train_in, data_train_out = data_train
    data_test_in, data_test_out = data_test
    
    ### Train prediction
    obs_train = data_train_out
    obs_train['richness'] = data_train_out['species'].sum(axis=1).reshape(-1,1)
    Ys_true_train = pd.DataFrame(data=obs_train['species'], columns=splist)
    
    pred_train = plant_model.predict(data_train_in)
    Ys_hat_train = pd.DataFrame(data=pred_train['species'], columns=splist)
    
    ### Test prediction
    obs_test = data_test_out
    obs_test['richness'] = data_test_out['species'].sum(axis=1).reshape(-1,1)
    Ys_true_test = pd.DataFrame(data=obs_test['species'], columns=splist) 
    
    prediction = plant_model.predict(data_test_in)
    pred_test = prediction
    Ys_hat_test = pd.DataFrame(data=pred_test['species'], columns=splist)   
    
    if method=='raw':
        pred_test['richness']= Ys_hat_test.sum(axis=1).values.reshape(-1,1)
        prediction = Ys_hat_test
        try:
            comm_perfs = plot_community_pred(obs_test,pred_test,tlist,algo='mtnn_raw',plot_residual=do_residual)
            species_perfs = None
        except:
            comm_perfs = None
        
    if method=='opth':
        th_list = optimize_thres(Ys_true_train,Ys_hat_train)
        species_perfs = eval_species_perfs(Ys_hat_test, Ys_true_test, th_list=th_list)
        rich, cm, cstd, _, _ = compute_functional_indices((Ys_hat_test>th_list)*1,T)
        prediction = (Ys_hat_test>th_list)*1
        pred_test = {
            'species': (Ys_hat_test>th_list)*1,
            'cm': cm,
            'cstd': cstd,
            'richness': rich
        }
        
        try:
            comm_perfs = plot_community_pred(obs_test,pred_test,tlist,algo='mtnn_opth',plot_residual=do_residual)
            comm_perfs['jaccard']=jaccard_score(y_true=obs_test['species'],y_pred=pred_test['species'], average='micro')
        except:
            comm_perfs = None
            
    return prediction, species_perfs, comm_perfs


def evaluate_pfal(P_train, P_test, Y_train, Y_test, splist, tlist, T, method='raw', do_residual=False, eps=1e-6):
   
    ### Train prediction
    obs_train = {}
    obs_train['species'] = Y_train
    obs_train['richness'] = Y_train.sum(axis=1).values.reshape(-1,1)
    Ys_true_train = Y_train
    
    pred_train ={}
    pred_train['species'] = P_train
    
    Ys_hat_train = P_train
    
    ### Test prediction
    rich, cm, cstd, _, _ = compute_functional_indices(Y_test,T)
    obs_test = {
            'species': Y_test,
            'cm': cm,
            'cstd': cstd,
            'richness': rich
    }
    
    Ys_true_test = Y_test
    pred_test ={}
    pred_test['species'] = P_test
    
    Ys_hat_test=P_test
    
    if method=='raw':
        pred_test['richness']= Ys_hat_test.sum(axis=1).values.reshape(-1,1)
        rich, cm, cstd, _, _ = compute_functional_indices(Ys_hat_test,T)
        pred_test = {
            'species': Ys_hat_test,
            'cm': cm,
            'cstd': cstd,
            'richness': rich
        }
        try:
            comm_perfs = plot_community_pred(obs_test,pred_test,tlist,algo='sdm_raw',plot_residual=do_residual)
            species_perfs = None
        except:
            comm_perfs = None
    
    if method=='th05':
        th_list = 0.5
        species_perfs = eval_species_perfs(Ys_hat_test, Ys_true_test, th_list=th_list)
        rich, cm, cstd, _, _ = compute_functional_indices((Ys_hat_test>th_list)*1,T)
        pred_test = {
            'species': (Ys_hat_test>th_list)*1,
            'cm': cm,
            'cstd': cstd,
            'richness': rich
        }
        
        try:
            comm_perfs = plot_community_pred(obs_test,pred_test,tlist,algo='sdm_th05',plot_residual=do_residual)
            comm_perfs['jaccard']=jaccard_score(y_true=obs_test['species'],y_pred=pred_test['species'], average='micro')
        except:
            comm_perfs = None
        
    if method=='opth':
        th_list = optimize_thres(Ys_true_train,Ys_hat_train)
        species_perfs = eval_species_perfs(Ys_hat_test, Ys_true_test, th_list=th_list)
        rich, cm, cstd, _, _ = compute_functional_indices((Ys_hat_test>th_list)*1,T)
        pred_test = {
            'species': (Ys_hat_test>th_list)*1,
            'cm': cm,
            'cstd': cstd,
            'richness': rich
        }
        
        try:
            comm_perfs = plot_community_pred(obs_test,pred_test,tlist,algo='sdm_opth',plot_residual=do_residual)
            comm_perfs['jaccard']=jaccard_score(y_true=obs_test['species'],y_pred=pred_test['species'], average='micro')
        except:
            comm_perfs = None
            
    if method=='prr':
        filter_prr = (Ys_hat_test.rank(axis=1,method='min',ascending=False)<=rich_test)*1
        Y_hat_prr = Ys_hat_test * filter_prr
        species_perfs = eval_species_perfs(Y_hat_prr, Ys_true_test, th_list=eps)
        
        rich, cwm, fdis, _, _ = compute_functional_indices((Y_hat_prr>eps)*1,T_scaled)
        pred_test = {
            'species': (Y_hat_prr>eps)*1,
            'cwm': cwm,
            'fdis': fdis,
            'richness': rich
        } 
        
        try:
            comm_perfs = plot_community_pred(obs_test,pred_test,tlist,algo='sdm_prr',plot_residual=do_residual)
            comm_perfs['jaccard']=jaccard_score(y_true=obs_test['species'],y_pred=pred_test['species'], average='micro')
        except:
            comm_perfs = None
        y_true = Ys_true_test.values
        y_pred = (Y_hat_prr>eps).values

    return comm_perfs