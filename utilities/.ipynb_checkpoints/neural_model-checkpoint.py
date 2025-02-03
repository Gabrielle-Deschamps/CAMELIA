import os
import numpy as np
import random
import tensorflow as tf
#import tensorflow_addons as tfa
from tensorflow.keras import layers, models, regularizers
tfk = tf.keras

def set_seed(seed: int = 1234) -> None:
    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)
    # When running on the CuDNN backend, two further options must be set
#     os.environ['TF_CUDNN_DETERMINISTIC'] = '1'
#     os.environ['TF_DETERMINISTIC_OPS'] = '1'
    # Set a fixed value for the hash seed
    os.environ["PYTHONHASHSEED"] = str(seed)
    print(f"Random seed set as {seed}")
    
    
class FCBlock(layers.Layer):
    def __init__(self, units, activation='relu', dropout_rate=0.0, l1=0.0, l2=0.0, name="fclayer", **kwargs):
        super().__init__(name=name, **kwargs)
        self.units = units
        self.activation = activation
        self.dropout_rate = dropout_rate
        self.l1 = l2
        self.l2 = l2
        
    def build(self, input_shape):
        self.dense = layers.Dense(self.units,name='dense',dtype=tf.float32)
        self.activ = layers.Activation(activation=self.activation,name='activation')
        
        self.activ_reg = layers.ActivityRegularization(l1=self.l1, l2=self.l2)
        if self.dropout_rate > 0.0:
            self.dropout = layers.Dropout(self.dropout_rate, name='dropout',dtype=tf.float32)
        else:
            self.dropout = None
            
    def call(self, inputs):
        x = self.dense(inputs)
        x = self.activ(x)
        x = self.activ_reg(x)
        
        if self.dropout is not None:
            x = self.dropout(x)
            
        return x
    
    def get_config(self):
        return {"units": self.units, "activation":self.activation, "dropout_rate":self.dropout_rate, "l1": self.l1, "l2":self.l2}
    
    
class MLPBlock(layers.Layer):
    def __init__(self, hidden_units, activation='relu', dropout_rate=0.0, l1=0.0, l2=0.0, name="mlpblock", **kwargs):
        super().__init__(name=name, **kwargs)
        
        self.hidden_units = hidden_units
        self.activation = activation
        self.dropout_rate = dropout_rate
        self.l1 = l1
        self.l2 = l2
        
    def build(self, input_shape):
        self.hidden_layers = []
        for nl, nn in enumerate(self.hidden_units):
            hidden = FCBlock(units=nn,activation=self.activation, dropout_rate=self.dropout_rate,l1=self.l1,l2=self.l2,name='fclayer%d'%nl)
            self.hidden_layers.append(hidden)
            
    def call(self, inputs):
        x = inputs
        
        for hidden_layer in self.hidden_layers:
            x = hidden_layer(x)
            
        return x
    
    def get_config(self):
        return {"hidden_units": self.hidden_units, "activation":self.activation, "dropout_rate":self.dropout_rate, "l1": self.l1, "l2":self.l2}    
    
    
class ProbaRanking(layers.Layer):
    def __init__(self, name="prr", **kwargs):
        super().__init__(name=name, **kwargs)
        
    def build(self, input_shape):
        self.topr_layer = layers.Lambda(lambda x: tf.cast(tf.less(tf.argsort(x[0],axis=1,direction='DESCENDING'),x[1]),tf.float32), name='target_prr')
        self.mask_layer = layers.Lambda(lambda x: x[0] * x[1],name='mask_prr')
            
    def call(self, inputs):
        in_proba, in_total = inputs
        topr_mask = self.topr_layer([in_proba, in_total])
        masked_prr = self.mask_layer([in_proba, topr_mask])
   
        return masked_prr

    def get_config(self):
        return {}


class WeightedMoments(layers.Layer):
    def __init__(self, trait_matrix, name="weightedmom", **kwargs):
        super().__init__(name=name, **kwargs)
        self.trait_matrix = trait_matrix
        
    def build(self, input_shape):
        self.species_trait = tf.constant(self.trait_matrix,dtype=tf.float32, name='trait_constant')
        self.funcdiv_layer = layers.Lambda(lambda x: tf.nn.weighted_moments(x[0], axes=[1], keepdims=False, frequency_weights=x[1], name='weighted_moments'))
          
    def call(self, inputs):
        weights = tf.expand_dims(inputs,axis=2, name='species_weights')
        #trait_output = tf.multiply(weights, self.species_trait, name='trait_weights') # to check here if not weighted 2 times du to 
        # the use of tf.nn.weighted_moments 

        cwm, cwvar = self.funcdiv_layer([self.species_trait, weights])
        cwstd = tf.sqrt(cwvar)
        
        return cwm, cwstd
    
    def get_config(self):
        return {"trait_matrix": self.trait_matrix}
    
    
class MTNN(tfk.Model):
    def __init__(self, num_species, num_traits, trait_matrix = None , hidden_units=[], apply_prr=False, optim_rich = False, hidden_activation='relu', out_activation='sigmoid',dropout_rate=0.0,l1=0.0,l2=0.0,name='mtnn',**kwargs):
        super().__init__(name=name, **kwargs)
        
        ### Problem setting
        self.num_species = num_species
        self.num_traits = num_traits
        
        ### Architecture
        self.hidden_units = hidden_units
        self.hidden_activation = hidden_activation
        self.out_activation=out_activation
        self.apply_prr = apply_prr
        self.optim_rich = optim_rich
        
        ## Regularization
        self.dropout_rate = dropout_rate
        self.l1 = l1
        self.l2 = l2
        
        #### Components
        self.trait_matrix = trait_matrix if trait_matrix is not None else None #np.expand_dims(trait_matrix,axis=0) 
        self.mlp_block = MLPBlock(hidden_units=self.hidden_units, activation=self.hidden_activation, dropout_rate=self.dropout_rate, l1=0.0, l2=0.0, name='hidden')
        self.classifier = FCBlock(units=self.num_species, activation=self.out_activation, dropout_rate=0.0, l1=self.l1, l2=self.l2, name='species_output')

        self.prr_layer = ProbaRanking(name='species_prr') if self.apply_prr else None
        self.richness_layer = layers.Lambda(lambda x: tf.reduce_sum(x,keepdims=True,axis=1),
                                            name='richness',dtype=tf.float32) if optim_rich else None
        self.funcdiv_layer = WeightedMoments(trait_matrix=self.trait_matrix, name='funcdiv')
        
    def call(self, inputs):
        x_in, rich_in = inputs
        
        features = self.mlp_block(x_in)
        species_output = self.classifier(features)
        out_layers = {'species': species_output}
        
        if self.optim_rich:
            richness_output = self.richness_layer(species_output)
            out_layers.update({'richness': richness_output})
        
        if self.trait_matrix is not None:
            if self.apply_prr:
                species_prr = self.prr_layer([species_output, rich_in])
                traitcwm, traitstd = self.funcdiv_layer(species_prr)
            else:
                traitcwm, traitstd = self.funcdiv_layer(species_output)
                
            out_layers.update({'cwm':traitcwm, 'fdis': traitstd})
        
        return out_layers    
    
    def get_config(self):
        return {"num_species":self.num_species, 
                "num_traits":self.num_traits, 
                "trait_matrix":self.trait_matrix, 
                "hidden_units":self.hidden_units, 
                "apply_prr":self.apply_prr, 
                "optim_rich":self.optim_rich, 
                "hidden_activation":self.hidden_activation, 
                "out_activation":self.out_activation,
                "dropout_rate":self.dropout_rate,
                "l1": self.l1,
                "l2": self.l2
               }    