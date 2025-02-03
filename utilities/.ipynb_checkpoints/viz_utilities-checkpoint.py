import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_history(history, use_skew=True, use_kurt=True, use_rich=True):
    learn_history = pd.DataFrame(history)
    learn_history['epoch']=np.arange(learn_history.shape[0])
    
    fig, ax = plt.subplots(3,2,figsize=(8,8))
    learn_history.plot(x='epoch',y='species_loss', ax=ax[0,0])
    learn_history.plot(x='epoch',y='val_species_loss', ax=ax[0,0])
    
    if use_rich:
        learn_history.plot(x='epoch',y='richness_loss', ax=ax[0,1])
        learn_history.plot(x='epoch',y='val_richness_loss', ax=ax[0,1])
        
    learn_history.plot(x='epoch',y='cwm_loss', ax=ax[1,0])
    learn_history.plot(x='epoch',y='val_cwm_loss', ax=ax[1,0])
    learn_history.plot(x='epoch',y='fdis_loss', ax=ax[1,1])
    learn_history.plot(x='epoch',y='val_fdis_loss', ax=ax[1,1])
    
    if use_skew:
        learn_history.plot(x='epoch',y='fskew_loss', ax=ax[2,0])
        learn_history.plot(x='epoch',y='val_fskew_loss', ax=ax[2,0])
    
    if use_kurt:
        learn_history.plot(x='epoch',y='fkurt_loss', ax=ax[2,1])
        learn_history.plot(x='epoch',y='val_fkurt_loss', ax=ax[2,1])
    
    fig.suptitle('Learning history - loss')
    fig.tight_layout()
    
    fig, ax = plt.subplots(4,2,figsize=(8,10))
    learn_history.plot(x='epoch',y='species_auc', ax=ax[0,0])
    learn_history.plot(x='epoch',y='val_species_auc', ax=ax[0,0])
    
    learn_history.plot(x='epoch',y='species_aupr', ax=ax[0,1])
    learn_history.plot(x='epoch',y='val_species_aupr', ax=ax[0,1]) 
    
    learn_history.plot(x='epoch',y='cwm_r_square', ax=ax[1,0])
    learn_history.plot(x='epoch',y='val_cwm_r_square', ax=ax[1,0])   
    
    learn_history.plot(x='epoch',y='fdis_r_square', ax=ax[1,1])
    learn_history.plot(x='epoch',y='val_fdis_r_square', ax=ax[1,1]) 
    
    learn_history.plot(x='epoch',y='species_recall', ax=ax[2,0])
    learn_history.plot(x='epoch',y='val_species_recall', ax=ax[2,0]) 
    
    learn_history.plot(x='epoch',y='species_precision', ax=ax[2,1])
    learn_history.plot(x='epoch',y='val_species_precision', ax=ax[2,1])
    
    learn_history.plot(x='epoch',y='richness_r_square', ax=ax[3,0])
    learn_history.plot(x='epoch',y='val_richness_r_square', ax=ax[3,0]) 
    
    
    return fig
    