import numpy as np
import astropy.constants as const
import astropy.units as u
from cosmolopy import constants

## The value of G in proper units. I've used
## the package cosmolopy which has this particular
## unit system built in. This is more convenient than
## converting astropy's G to the correct unit system.
G = constants.G_const_Mpc_Msun_s

## Compute step size. In order to be consistent in our units, 
## we will have to use a step_size in seconds. However, it is built
## to accept a simple input of step
step_size_in_years = 10000000
step_size = (step_size_in_years * u.yr).to(u.s).value

## Galaxy mass - we are assuming all galaxies are 1e12 Solar Masses
M = 1e12

## The limit of L/D where L is the dimensions of the node and 
## D is the distance from a given galaxy to that Node.
## If D * theta_lim < L, we will consider individual members of the group when we 
## perform our acceleration calculations. If not, we assume we 
## are far enough away to calculate the node as a single mass point.
theta_lim = 0.7

## The threshold for force softening. If the distance to a neighboring node or
## galaxies COM is less than this value, we will compute the acceleration 
## with a softening term, eps, to avoid singularities or runaway accelerations.
eps = 0.001


class Tree:
    def __init__(self, origin, size, masses, points, ids, leaves=[]):
        '''
        Summary:
        Initializes a tree structure that will serve as a global node
        within which all galaxies will reside. We will add galaxies to the 
        structure and reproduce nodes until each node only holds one
        galaxy. Said another way, we will reproduce until we reach only leaves.
        
        Parameters
        ----------
        origin : the coordinates of the center of the node.
        masses : the array of masses. In our case it is easy, we are assuming all galaxies
                 have a mass of 1e12. However, it is still good to include this parameter for
                 future generalization.
        points : An array of coordinate points.
        ids : the identifiers of each node in the Tree. Helps to keep track of where
              each parent and child belongs.
        leaves : keeping track of when a node reach a leaf. This is our goal after all.
    
        '''
        
        ## The coordinates of the center of the node. 
        self.origin = origin
        
        ## The side length of the node.
        self.size = size
        
        ## How many children each node has. Starting from 0.
        self.children = []
 
        ## Check if we only have 1 point in the node. If we do, we can easily calculate
        ## the center of mass, total mass, attribute an id, and calculate the gravitational
        ## acceleration at that point.
        if len(points) == 1:
            
            ## Store the fact that we've reached a leaf
            leaves.append(self)
            
            ## Once we've reached a leaf, the center of mass, mass, 
            ## and acceleration do not even require computations. 
            self.COM = points[0]
            self.mass = masses[0]
            self.id = ids[0]
            self.g = np.zeros(3)       
        else:
            ## If there are more than 1 galaxy in the node, we have to 
            ## make the node reproduce. 
            self.Reproduce(points, masses, ids, leaves)
 
            ## We can now sum the mass and position vectors to compute the 
            ## total mass and center of mass for each node. 
            com_total = np.zeros(3)
            m_total = 0.  
            for child in self.children:
                ## Find the mass of each child.
                m = child.mass
                
                ## Find the center of mass of each child.
                com = child.COM
                
                ## Sum all the masses in a node.
                m_total += m
                
                ## Find the center of mass vector for a node.
                com_total += com * m  
            
            self.mass = m_total
            self.COM = com_total / self.mass  
 
    def Reproduce(self, points, masses, ids, leaves):
        '''
        Summary:
        Takes a Node and causes it to reproduce. i.e. add children. 
        
        Parameters
        ----------
        masses : the array of masses. In our case it is easy, we are assuming all galaxies
                 have a mass of 1e12. However, it is still good to include this parameter for
                 future generalization.
        points : An array of coordinate points.
        ids : the identifiers of each node in the Tree. Helps to keep track of where
              each parent and child belongs.
        leaves : keeping track of when a node reach a leaf. This is our goal after all.
    
        '''
        
        ## This is a cool vectorization that seems to be a popular implementation
        ## if Tree-codes. 
        
        ## We identify all points that are above the midpoint or below the midpoint.
        ## For example, if the origin is 5 and the first point is 4, octant_index
        ## will store a False in place of the actual coordinate. This way, we 
        ## can easily and quickly narrow down where each galaxy belongs in the 
        ## reproduction structure.
        child_index = (points > self.origin)
        
        ## Loop over the 8 children.
        for i in range(2): 
            for j in range(2):
                for k in range(2):
                    
                    ## Check to see if the child_index is True or False. 
                    ## If True, the points at that index belong in the child node.
                    ## If False, the points belong in a different child node.
                    in_child_node = np.all(child_index == np.bool_([i,j,k]), axis=1)
                    
                    ## Some child nodes may not have any galaxies.
                    if not np.any(in_child_node): 
                        continue
                    
                    ## We now need to set up the children. 
                    ## In order to initialize them, we need to know where their 
                    ## origin will be. 
                    child_offset = 0.5*self.size*(np.array([i,j,k])-0.5) 
                    
                    ## Initialize the child nodes.
                    self.children.append(Tree(self.origin+child_offset,
                                                 self.size/2,
                                                 masses[in_child_node],
                                                 points[in_child_node],
                                                 ids[in_child_node],
                                                 leaves))
  
def calc_contribution(Node, galaxy):
    '''
    Summary:
    Calculates the gravitational acceleration on a galaxy from
    a given node. The tree class will walk this structure 
    to decide if we are far enough away to calculate the COM of a 
    whole node or if considering individual galaxies is necessary.
    
    Parameters
    ----------
    Node : A Node object that may contain children or be a leaf.
    galaxy : the galaxy object that is being accelerated due to the mass
             of the node object.
    '''
    
    ## The vector that defines the distance between the node object
    ## and the galaxy.
    r_vector = Node.COM - galaxy.COM
    
    ## The magnitude of the r vector.
    r_mag = np.sqrt(np.sum(r_vector**2)) 
    
    ## If r_mag == 0, we would be calculating the galaxies influence on itself
    ## That is useless computation. Therefore, we are only concerned with cases 
    ## where r_mag != 0.
    if r_mag>0:
        
        # if the node only has one particle or theta is small enough,
        #  add the field contribution to value stored in node.g
        
        ## If the node only has one galaxy, we can add the contribution
        ## without further reproduction of the Node.
        ## We also incorporate force softening. 
        if (len(Node.children)==0): 
            ## If r_mag is less than eps, we apply force softening. If not, we
            ## calculating the acceleration like normal. 
            if r_mag <= eps:
                galaxy.g += G * Node.mass * r_vector/(r_mag * (r_mag**2 + eps**2))
            elif r_mag > eps:
                galaxy.g += G * Node.mass * r_vector/r_mag**3
        
        ## If the ratio of the length of the node to the distance between the
        ## the galaxy and the node is below a threshold set above, we can treat
        ## the entire node like a single mass point concentrated at its COM.
        ## We also incorporate force softening.
        elif (Node.size/r_mag < theta_lim):
            if r_mag <= eps:
                galaxy.g += G * Node.mass * r_vector/(r_mag * (r_mag**2 + eps**2))
            elif r_mag > eps:
                galaxy.g += G * Node.mass * r_vector/r_mag**3
        
        ## If neither of the above criteria are met, the node needs to 
        ## split further into children.
        else:
            for child in Node.children: calc_contribution(child, galaxy)
                              

def calc_accel(points, masses):
    '''
    Summary:
    Calculates the acceleration on each point from all other points. 
    
    Parameters
    ----------
    points : a list of galaxy coordinates.
    masses : the mass of each node. 
    '''
    
    ## Calculate the origin in each dimension. 
    origin = (np.max(points,axis=0)+np.min(points,axis=0))/2
    
    ## Calculate the dimensions of the bounding box.
    dimensions = np.max(np.max(points,axis=0)-np.min(points,axis=0))
    
    ## Keeping track of leaves.
    leaves = []
    
    ## Establish the root node and build the tree.
    root = Tree(origin, dimensions, masses, points, np.arange(len(masses)), leaves)
 
    ## Store the accelerations
    a = np.empty_like(points)
    
    ## Loop over all leaves in the Tree.
    for i,leaf in enumerate(leaves):
        
        ## Sum all the accelerations on each galaxy in the array.
        calc_contribution(root, leaf)  
        a[leaf.id] = leaf.g
 
    return a


def calc_position(galaxies, history_galaxies, accel):
    '''
    Summary:
    Calculate the position of each galaxy in the array after each iteration.

    Parameters
    ----------
    galaxies : the coordinate of every galaxy in the current timestep.
    history_galaxies : the coordinate of every galaxy in the last timestep.
    accel : the acceleration of each galaxy. Calculated in barnes_hut.py
    '''
    
    new_positions = step_size**2 * accel + 2 * galaxies - history_galaxies
    
    return new_positions
