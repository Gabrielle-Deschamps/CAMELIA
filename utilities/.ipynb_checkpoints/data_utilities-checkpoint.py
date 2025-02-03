import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler

def compute_richness(Y):
    A = Y.values.reshape((nsite,nspecies,1)) 
    return A.sum(axis=1)

def compute_functional_indices(Y,T):
    nsite, nspecies = Y.shape
    _, ntrait = T.shape
    
    # First, we'll build the site x species x trait matrix
    A = Y.values.reshape((nsite,nspecies,1)) 
    B = T.values.reshape((1,nspecies,ntrait))
    Y_trait = A * B
    Y_rich = Y.sum(axis=1).values.reshape(-1,1)
    
    # Community weighted mean trait
    T_mean = Y_trait.sum(axis=1)/Y_rich
    T_mean_ext = np.expand_dims(T_mean,axis=1)
    T_cent = Y_trait - T_mean_ext
    T_var = (A*(T_cent ** 2)).sum(axis=1) / Y_rich
    
    # Community weighted standard deviation of trait
    T_std = np.sqrt(T_var)
    
    Z_stat = A*T_cent / T_std.reshape((nsite,1,ntrait))
    
    # Community weighted skewness of trait
    T_skew = np.power(Z_stat, 3).sum(axis=1)/Y_rich
    
    # Community weighted kurtosis of trait
    T_kurt = np.power(Z_stat, 4).sum(axis=1)/Y_rich
    
    return Y_rich, T_mean, T_std, T_skew, T_kurt

def get_species_data(metadata,Y,fold=2,eval_mode='extrapolation'):
    nsite, nspecies = Y.shape
    
    if eval_mode=="interpolation":
        P = pd.read_csv('data/species_models_pa_inter.csv')
    else:
        P = pd.read_csv('data/species_models_pa_extra.csv',index_col=0)
    
    Y.index = P.index = metadata.index = np.arange(nsite)
    
    print('Using %s fold %d'%(eval_mode,fold))
    test_idx=metadata.query('group_%s==%d'%(eval_mode,fold)).index
    train_idx=metadata.query('group_%s!=%d'%(eval_mode,fold)).index
    
    ### Species output
    Y_train = Y.loc[train_idx]
    Y_test = Y.loc[test_idx]
    
    ### Species input
    P_train = P.loc[train_idx]
    P_test = P.loc[test_idx] 
    
    P_train.columns = Y_train.columns
    P_test.columns = Y_test.columns
    
    return P_train, Y_train, P_test, Y_test


def generate_dataset(metadata, Y, T, fold=2,eval_mode='extrapolation',in_comm='obs',out_comm='obs', trait_scaler=None):
    print('Using %s fold %d'%(eval_mode,fold))
    nsite, nspecies = Y.shape
    
    if eval_mode=="interpolation":
        P = pd.read_csv('data/species_models_pa_inter.csv')
    else:
        P = pd.read_csv('data/species_models_pa_extra.csv',index_col=0)
    
    Y.index = P.index = metadata.index = np.arange(nsite)
    test_idx=metadata.query('group_%s==%d'%(eval_mode,fold)).index
    train_idx=metadata.query('group_%s!=%d'%(eval_mode,fold)).index
    
    ### Species output
    Y_train = Y.loc[train_idx]
    Y_test = Y.loc[test_idx] 
        
    ### Species input
    P_train = P.loc[train_idx]
    P_test = P.loc[test_idx]
    
    ### Scaling traits (eventually)
    if trait_scaler is not None:
        T = pd.DataFrame(data=trait_scaler.transform(T),columns=T.columns,index=T.index)
    
    ### Community indices
    #### Observed
    Y_rich_obs_train, T_mean_obs_train, T_std_obs_train, _, _ = compute_functional_indices(Y_train,T)
    Y_rich_obs_test, T_mean_obs_test, T_std_obs_test, _, _ = compute_functional_indices(Y_test,T)
    
    #### MEM predictions
    if eval_mode=='interpolation':
        T_MEM_cwm = pd.read_csv('data/mem/CM_obs_degrad.csv',sep=';',decimal=',')
        T_MEM_fdis = pd.read_csv('data/mem/uFDis_obs_degrad.csv',sep=';',decimal=',')
        Y_MEM_rich = pd.read_csv('data/mem/SR_obs_degrad.csv',sep=';',decimal=',')  
    else:
        T_MEM_cwm = pd.read_csv('data/mem/CM_trait_based_extra.csv',sep=';',decimal=',')
        T_MEM_fdis = pd.read_csv('data/mem/uFDis_trait_based_extra.csv',sep=';',decimal=',')
        Y_MEM_rich = pd.read_csv('data/mem/SR_extra.csv',sep=';',decimal=',')
    
    if trait_scaler is not None:
        T_MEM_cwm = pd.DataFrame(data=trait_scaler.transform(T_MEM_cwm),columns=T_MEM_cwm.columns,index=T_MEM_cwm.index)
        T_MEM_fdis = pd.DataFrame(data=trait_scaler.transform(T_MEM_fdis),columns=T_MEM_fdis.columns,index=T_MEM_fdis.index)
                   
    Y_rich_mem_train = Y_MEM_rich.loc[train_idx,['SR']].values.astype(float)
    T_mean_mem_train = T_MEM_cwm.loc[train_idx].values
    T_std_mem_train = T_MEM_fdis.loc[train_idx].values

    Y_rich_mem_test = Y_MEM_rich.loc[test_idx,['SR']].values.astype(float)
    T_mean_mem_test = T_MEM_cwm.loc[test_idx].values
    T_std_mem_test = T_MEM_fdis.loc[test_idx].values      
    
    ### Format outputs    
    if out_comm=='obs':
        print('Output: observed community indices')
        train_data_out = {
            'species': Y_train.values,
            'richness': Y_rich_obs_train,
            'cwm': T_mean_obs_train,
            'fdis': T_std_obs_train
        }
        
        test_data_out = {
            'species': Y_test.values,
            'richness': Y_rich_obs_test,
            'cwm': T_mean_obs_test,
            'fdis': T_std_obs_test
        }
        
    else:
        print('Output: MEM predicted community indices')
        train_data_out = {
            'species': Y_train.values,
            'richness': Y_rich_mem_train,
            'cwm': T_mean_mem_train,
            'fdis': T_std_mem_train
        }
        
        test_data_out = {
            'species': Y_test.values,
            'richness': Y_rich_mem_test,
            'cwm': T_mean_mem_test,
            'fdis': T_std_mem_test
        }
    
    ### Format inputs
    rich_std = StandardScaler()
    cwm_std = StandardScaler() if trait_scaler is None else None
    fdis_std = StandardScaler() if trait_scaler is None else None
    
    if in_comm=='obs': 
        print('Input: observed community indices')
        rich_std.fit(Y_rich_obs_train)
        
        if trait_scaler is None:
            cwm_std.fit(T_mean_obs_train)
            fdis_std.fit(T_std_obs_train)
    
        train_data_in = (np.concatenate([P_train, 
                                        rich_std.transform(Y_rich_obs_train),
                                        cwm_std.transform(T_mean_obs_train) if trait_scaler is None else T_mean_obs_train, 
                                        fdis_std.transform(T_std_obs_train) if trait_scaler is None else T_std_obs_train], axis=1),
                         
                         Y_rich_obs_train.astype(np.int32))

        test_data_in = (np.concatenate([P_test, 
                                       rich_std.transform(Y_rich_obs_test),
                                       cwm_std.transform(T_mean_obs_test) if trait_scaler is None else T_mean_obs_test, 
                                       fdis_std.transform(T_std_obs_test) if trait_scaler is None else T_std_obs_test], axis=1),
                        
                        Y_rich_obs_test.astype(np.int32))
        
    else:
        print('Input: MEM predicted community indices')
        rich_std.fit(Y_rich_mem_train)
        
        if trait_scaler is None:
            cwm_std.fit(T_mean_mem_train)
            fdis_std.fit(T_std_mem_train)
        
        train_data_in = (np.concatenate([P_train,
                                        rich_std.transform(Y_rich_mem_train),
                                        cwm_std.transform(T_mean_mem_train) if trait_scaler is None else T_mean_mem_train, 
                                        fdis_std.transform(T_std_mem_train) if trait_scaler is None else T_std_mem_train], axis=1),
                         
                         Y_rich_mem_train.astype(np.int32))

        test_data_in = (np.concatenate([P_test, 
                                       rich_std.transform(Y_rich_mem_test),
                                       cwm_std.transform(T_mean_mem_test) if trait_scaler is None else T_mean_mem_test, 
                                       fdis_std.transform(T_std_mem_test) if trait_scaler is None else T_std_mem_test], axis=1),
                        
                        Y_rich_mem_test.astype(np.int32))
        
        
    return (train_data_in, train_data_out), (test_data_in, test_data_out), {'richness':rich_std, 'cwm':cwm_std, 'fdis': fdis_std}





def generate_dataset_2(metadata, Y, T, fold=2,eval_mode='extrapolation',in_comm='obs',out_comm='obs', trait_scaler=None):
    print('Using %s fold %d'%(eval_mode,fold))
    nsite, nspecies = Y.shape
    
    if eval_mode=="interpolation":
        P = pd.read_csv('data/species_models_pa_inter.csv')
    else:
        P = pd.read_csv('data/species_models_pa_extra.csv',index_col=0)
    
    Y.index = P.index = metadata.index = np.arange(nsite)
    test_idx=metadata.query('group_%s==%d'%(eval_mode,fold)).index
    train_idx=metadata.query('group_%s!=%d'%(eval_mode,fold)).index
    
    ### Species output
    Y_train = Y.loc[train_idx]
    Y_test = Y.loc[test_idx] 
        
    ### Species input
    P_train = P.loc[train_idx]
    P_test = P.loc[test_idx]
    
    ### Scaling traits (eventually)
    if trait_scaler is not None:
        T = pd.DataFrame(data=trait_scaler.transform(T),columns=T.columns,index=T.index)
    
    ### Community indices
    #### Observed
    Y_rich_obs_train, T_mean_obs_train, T_std_obs_train, _, _ = compute_functional_indices(Y_train,T)
    Y_rich_obs_test, T_mean_obs_test, T_std_obs_test, _, _ = compute_functional_indices(Y_test,T)
    
    #### MEM predictions
    if eval_mode=='interpolation':
        T_MEM_cwm = pd.read_csv('data/mem/CM_trait_based_inter.csv',sep=';',decimal=',')
        T_MEM_fdis = pd.read_csv('data/mem/uFDis_trait_based_inter.csv',sep=';',decimal=',')
        Y_MEM_rich = pd.read_csv('data/mem/SR_inter.csv',sep=';',decimal=',')  
    else:
        T_MEM_cwm = pd.read_csv('data/mem/CM_trait_based_extra.csv',sep=';',decimal=',')
        T_MEM_fdis = pd.read_csv('data/mem/uFDis_trait_based_extra.csv',sep=';',decimal=',')
        Y_MEM_rich = pd.read_csv('data/mem/SR_extra.csv',sep=';',decimal=',')
    
    if trait_scaler is not None:
        T_MEM_cwm = pd.DataFrame(data=trait_scaler.transform(T_MEM_cwm),columns=T_MEM_cwm.columns,index=T_MEM_cwm.index)
        T_MEM_fdis = pd.DataFrame(data=trait_scaler.transform(T_MEM_fdis),columns=T_MEM_fdis.columns,index=T_MEM_fdis.index)
                   
    Y_rich_mem_train = Y_MEM_rich.loc[train_idx,['SR']].values.astype(float)
    T_mean_mem_train = T_MEM_cwm.loc[train_idx].values
    T_std_mem_train = T_MEM_fdis.loc[train_idx].values

    Y_rich_mem_test = Y_MEM_rich.loc[test_idx,['SR']].values.astype(float)
    T_mean_mem_test = T_MEM_cwm.loc[test_idx].values
    T_std_mem_test = T_MEM_fdis.loc[test_idx].values      
    
    ### Format outputs    
    if out_comm=='obs':
        print('Output: observed community indices')
        train_data_out = {
            'species': Y_train.values,
            'richness': Y_rich_obs_train,
            'cwm': T_mean_obs_train,
            'fdis': T_std_obs_train
        }
        
        test_data_out = {
            'species': Y_test.values,
            'richness': Y_rich_obs_test,
            'cwm': T_mean_obs_test,
            'fdis': T_std_obs_test
        }
        
    else:
        print('Output: MEM predicted community indices')
        train_data_out = {
            'species': Y_train.values,
            'richness': Y_rich_mem_train,
            'cwm': T_mean_mem_train,
            'fdis': T_std_mem_train
        }
        
        test_data_out = {
            'species': Y_test.values,
            'richness': Y_rich_mem_test,
            'cwm': T_mean_mem_test,
            'fdis': T_std_mem_test
        }
    
    ### Format inputs
    if in_comm=='obs':
        print('Input: observed community indices')
        train_data_in = {
            'species': P_train.values,
            'richness': Y_rich_obs_train,
            'cwm': T_mean_obs_train,
            'fdis': T_std_obs_train
        }
        
        test_data_in = {
            'species': P_test.values,
            'richness': Y_rich_obs_test,
            'cwm': T_mean_obs_test,
            'fdis': T_std_obs_test
        }
        
    else:
        print('Input: MEM predicted community indices')
        train_data_in = {
            'species': P_train.values,
            'richness': Y_rich_mem_train,
            'cwm': T_mean_mem_train,
            'fdis': T_std_mem_train
        }
        
        test_data_in = {
            'species': P_test.values,
            'richness': Y_rich_mem_test,
            'cwm': T_mean_mem_test,
            'fdis': T_std_mem_test
        }
   
    return (train_data_in, train_data_out), (test_data_in, test_data_out)