import pandas as pd
import numpy as np

def read_dataset(dataset):
    '''
    Summary:
    Reads a given dataset and outputs the distance and redshift
    to that galaxy as well as the uncertanties in each measurement.
    
    Parameters
    ----------
    dataset : a dataset in the working directory. 
    
    NOTE: New datasets can easily be added simply by adding an 
    if statement below corresponding to your new dataset. You
    will need to work out how to extract the distance and redshift
    from the dataset. Datasets tend to be incredibly different from
    each other so there is no one-size-fits-all approach here.
    '''
    
    
    if dataset == 'type_1a.txt':
        
        ## Read the dataset and convert it to a numpy array. 
        ## I think the pandas read_csv function works much
        ## better than numpy's loadtxt. So I use it and convert.
        type_1a = pd.read_csv('type_1a.txt',delimiter='\t')
        type_1a = type_1a.to_numpy()
        
        ## Extract the distance and add some uncertainty to it
        ## Unforunately, some datasets do not provide uncertainties
        ## alongside their data so I have to provide some uniform
        ## uncertainties.
        distance = type_1a[:,1]
        distance_uncert = np.random.uniform(0,10,len(distance))
        
        ## Extract the redshift and redshift uncertainty.
        redshift = type_1a[:,2]
        redshift_uncert = type_1a[:,3]
        
        return distance, distance_uncert, redshift, redshift_uncert
    
    if dataset == 'SDSS_lum_distance_cepheids.csv':
        
        ## Read the dataset and convert it to a numpy array. 
        ## I think the pandas read_csv function works much
        ## better than numpy's loadtxt. So I use it and convert.
        SDSS = pd.read_csv('SDSS_lum_distance_cepheids.csv',delimiter='\t')
        SDSS = SDSS.to_numpy()
        
        ## Extract the distance and add some uncertainty to it
        ## Unforunately, some datasets do not provide uncertainties
        ## alongside their data so I have to provide some uniform
        ## uncertainties.
        distance = SDSS[:,0]
        distance_uncert = np.random.uniform(0,10,len(distance))
        
        ## Extract the redshift and redshift uncertainty.
        redshift = SDSS[:,1]
        redshift_uncert = np.random.uniform(0,0.6,len(redshift))
        
        return distance, distance_uncert, redshift, redshift_uncert
        
    
    if dataset == 'leda_distance.csv':
        
        ## Read the dataset and convert it to a numpy array. 
        ## I think the pandas read_csv function works much
        ## better than numpy's loadtxt. So I use it and convert.
        leda = pd.read_csv('leda_distance.csv',delimiter=',')
        leda = leda.to_numpy()
        
        ## Extract the distance and add some uncertainty to it
        ## Unforunately, some datasets do not provide uncertainties
        ## alongside their data so I have to provide some uniform
        ## uncertainties.
        
        ## In this case, the dataset contains the distance
        ## modulus. Therefore, it has to be converted to
        ## the actual distance in megaparsecs via the conversion
        ## d in pc = 10^(distance_mod/5 + 1)
        distance_mod = leda[:,1]
        distance = 10**(distance_mod/5 + 1)/(10**6)
        distance_uncert = np.random.uniform(0,1,len(distance))
        
        ## In this dataset, the redshift isn't explicitly given but
        ## the recessional velocity is. We calculate the redshift as 
        ## z = v/c
        vsqr = leda[:,2]
        redshift = vsqr/299792
        redshift = np.random.uniform(0,5e-4,len(distance))
        
        return distance, distance_uncert, redshift, redshift_uncert
    
    
        
        
        
        
        
        
        
        
        
        