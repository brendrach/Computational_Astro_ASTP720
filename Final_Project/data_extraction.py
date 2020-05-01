import pandas as pd
import numpy as np

c = 299792 ## km/s

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
        type_1a = type_1a.dropna()
        type_1a = type_1a.to_numpy()
        
        ## Sort the data by the first column so the 
        ## galaxies are ordered by their distance.
        type_1a = type_1a[type_1a[:,0].argsort()]
        
        ## Extract the distance and add some uncertainty to it
        ## Unforunately, some datasets do not provide uncertainties
        ## alongside their data so I have to provide some uniform
        ## uncertainties.
        distance = type_1a[:,1]
        distance_uncert = np.random.uniform(0,10,len(distance))
        
        ## Extract the redshift and compute the 
        ## velocity as v = z*c.
        redshift = type_1a[:,2]
        velocity = redshift * c
        velocity_uncert = np.random.uniform(0,10,len(velocity))
        
        distance = distance.astype('float64')
        distance_uncert = distance_uncert.astype('float64')
        velocity = velocity.astype('float64')
        velocity_uncert = velocity_uncert.astype('float64')
        
        
        return distance, distance_uncert, velocity, velocity_uncert
    
    if dataset == 'SDSS_lum_distance_cepheids.csv':
        
        ## Read the dataset and convert it to a numpy array. 
        ## I think the pandas read_csv function works much
        ## better than numpy's loadtxt. So I use it and convert.
        SDSS = pd.read_csv('SDSS_lum_distance_cepheids.csv',delimiter='\t')
        SDSS = SDSS.dropna()
        SDSS = SDSS.to_numpy()
        
        ## Sort the data by the first column so the 
        ## galaxies are ordered by their distance.
        SDSS = SDSS[SDSS[:,0].argsort()]
        
        ## Extract the distance and add some uncertainty to it
        ## Unforunately, some datasets do not provide uncertainties
        ## alongside their data so I have to provide some uniform
        ## uncertainties.
        distance = SDSS[:,0]
        distance_uncert = np.random.uniform(0,10,len(distance))
        
        ## Extract the redshift and compute the 
        ## velocity as v = z*c.
        redshift = SDSS[:,1]
        velocity = redshift * c
        velocity_uncert = np.random.uniform(0,10,len(velocity))
        
        distance = distance.astype('float64')
        distance_uncert = distance_uncert.astype('float64')
        velocity = velocity.astype('float64')
        velocity_uncert = velocity_uncert.astype('float64')
        
        return distance, distance_uncert, velocity, velocity_uncert
        
    
    if dataset == 'leda_distance.csv':
        
        ## Read the dataset and convert it to a numpy array. 
        ## I think the pandas read_csv function works much
        ## better than numpy's loadtxt. So I use it and convert.
        leda = pd.read_csv('leda_distance.csv',delimiter=',')
        leda = leda.dropna()
        leda = leda.to_numpy()
        
        ## Sort the data by the first column so the 
        ## galaxies are ordered by their distance.
        leda = leda[leda[:,0].argsort()]
        
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
        velocity = leda[:,2]
        velocity_uncert = np.random.uniform(0,10,len(velocity))
        
        distance = distance.astype('float64')
        distance_uncert = distance_uncert.astype('float64')
        velocity = velocity.astype('float64')
        velocity_uncert = velocity_uncert.astype('float64')
        
        return distance, distance_uncert, velocity, velocity_uncert
    
    
        
        
        
        
        
        
        
        
        
        