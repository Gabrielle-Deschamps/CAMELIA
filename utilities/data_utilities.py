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

def get_species_data(metadata, Y, P, fold=2):
    Y = Y.reindex(metadata.index)
    P = P.reindex(metadata.index)
    
    print('Using interpolation fold %d'%(fold))
    test_idx = metadata[metadata['fold'] == fold].index
    train_idx = metadata[metadata['fold'] != fold].index
    
    ### Species output
    Y_train = Y.loc[train_idx]
    Y_test = Y.loc[test_idx]
    
    ### Species input
    P_train = P.loc[train_idx]
    P_test = P.loc[test_idx] 
    
    P_train.columns = Y_train.columns
    P_test.columns = Y_test.columns
    
    return P_train, Y_train, P_test, Y_test


def generate_dataset(metadata, Y, P, T, fold=2, in_comm='obs',out_comm='obs'):
    print('Using interpolation fold %d'%(fold))    
    Y = Y.reindex(metadata.index)
    P = P.reindex(metadata.index)    
    test_idx = metadata[metadata['fold'] == fold].index
    train_idx = metadata[metadata['fold'] != fold].index
    
    ### Species output
    Y_train = Y.loc[train_idx]
    Y_test = Y.loc[test_idx] 
        
    ### Species input
    P_train = P.loc[train_idx]
    P_test = P.loc[test_idx]
    
    ### Community indices
    #### Observed
    Y_rich_obs_train, T_mean_obs_train, T_std_obs_train, _, _ = compute_functional_indices(Y_train,T)
    Y_rich_obs_test, T_mean_obs_test, T_std_obs_test, _, _ = compute_functional_indices(Y_test,T)
    
    #### community_indices
    if in_comm!='obs':
        T_input_cm = pd.read_csv(f'data/community_indices/CM_{in_comm}.csv', index_col=0)
        T_input_cstd = pd.read_csv(f'data/community_indices/CSTD_{in_comm}.csv', index_col=0)
        Y_input_rich = pd.read_csv(f'data/community_indices/SR_{in_comm}.csv', index_col=0)
                   
        Y_rich_input_train = Y_input_rich.loc[train_idx,['SR']].values.astype(float)
        T_mean_input_train = T_input_cm.loc[train_idx].values
        T_std_input_train = T_input_cstd.loc[train_idx].values

        Y_rich_input_test = Y_input_rich.loc[test_idx,['SR']].values.astype(float)
        T_mean_input_test = T_input_cm.loc[test_idx].values
        T_std_input_test = T_input_cstd.loc[test_idx].values      
    
    ### Format outputs    
    if out_comm=='obs':
        print('Output: observed community indices')
        train_data_out = {
            'species': Y_train.values,
            'richness': Y_rich_obs_train,
            'cm': T_mean_obs_train,
            'cstd': T_std_obs_train
        }
        
        test_data_out = {
            'species': Y_test.values,
            'richness': Y_rich_obs_test,
            'cm': T_mean_obs_test,
            'cstd': T_std_obs_test
        }
        
    else:
        print(f'Output: {in_comm} community indices')
        train_data_out = {
            'species': Y_train.values,
            'richness': Y_rich_input_train,
            'cm': T_mean_input_train,
            'cstd': T_std_input_train
        }
        
        test_data_out = {
            'species': Y_test.values,
            'richness': Y_rich_input_test,
            'cm': T_mean_input_test,
            'cstd': T_std_input_test
        }
    
    ### Format inputs
    rich_std = StandardScaler()
    cm_std = StandardScaler()
    cstd_std = StandardScaler()
    
    if in_comm=='obs': 
        print('Input: observed community indices')
        rich_std.fit(Y_rich_obs_train)
        
        cm_std.fit(T_mean_obs_train)
        cstd_std.fit(T_std_obs_train)
    
        train_data_in = (np.concatenate([P_train, 
                                        rich_std.transform(Y_rich_obs_train),
                                        cm_std.transform(T_mean_obs_train),
                                        cstd_std.transform(T_std_obs_train)], axis=1),
                         
                         Y_rich_obs_train.astype(np.int32))

        test_data_in = (np.concatenate([P_test, 
                                       rich_std.transform(Y_rich_obs_test),
                                       cm_std.transform(T_mean_obs_test),
                                       cstd_std.transform(T_std_obs_test)], axis=1),
                        
                        Y_rich_obs_test.astype(np.int32))
        
    else:
        print(f'Input: {in_comm} community indices')
        rich_std.fit(Y_rich_input_train)
        cm_std.fit(T_mean_input_train)
        cstd_std.fit(T_std_input_train)
        
        train_data_in = (np.concatenate([P_train,
                                        rich_std.transform(Y_rich_input_train),
                                        cm_std.transform(T_mean_input_train), 
                                        cstd_std.transform(T_std_input_train)], axis=1),
                         
                         Y_rich_input_train.astype(np.int32))

        test_data_in = (np.concatenate([P_test, 
                                       rich_std.transform(Y_rich_input_test),
                                       cm_std.transform(T_mean_input_test), 
                                       cstd_std.transform(T_std_input_test)], axis=1),
                        
                        Y_rich_input_test.astype(np.int32))
        
        
    return (train_data_in, train_data_out), (test_data_in, test_data_out), {'richness':rich_std, 'cm':cm_std, 'cstd': cstd_std}