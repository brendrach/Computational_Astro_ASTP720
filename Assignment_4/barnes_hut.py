import numpy as np
import astropy.constants as const
import astropy.units as u

galaxies_t0 = np.load('galaxies0.npy')

max_x = np.max(galaxies_t0[:,0]) + 1
min_x = np.min(galaxies_t0[:,0]) - 1

max_y = np.max(galaxies_t0[:,1]) + 1
min_y = np.min(galaxies_t0[:,1]) - 1

max_z = np.max(galaxies_t0[:,2]) + 1
min_z = np.min(galaxies_t0[:,2]) - 1

init_side_len = (max_x, max_y, max_z)



class galaxy:
    
    
    def __init__(self, x_coord, y_coord, z_coord, mass = 1e12):
        '''
        Summary:
        Initializes a galaxy with given coordinates and mass.
        
        Parameters
        ----------
        xyz_coord : the coordinates for the galaxy. 
        mass : the mass of the galaxy.
        '''
        self.x = x_coord
        self.y = y_coord
        self.z = z_coord
        self.mass = mass
        self.coords = [self.x, self.y, self.z]
        self.accel = [0,0,0]
        
        
class Node:
    
    def __init__(self, x_min, y_min, z_min, x_max, y_max, z_max, parent=None, children=None):
        '''
        Summary:
        Initializes a node (cube) with corner (x,y,z) that extends in
        each direction a side length (side_len).
        
        Parameters
        ----------
        xyz_min : the minimum range of the global tree
        xyz_max : the maximum range of the global tree
        parent : the parent node.
        children : the number of children the node has.
        '''
        
        
        self.x_min = x_min
        self.y_min = y_min
        self.z_min = z_min
        self.x_max = x_max
        self.y_max = y_max
        self.z_max = z_max
        
        ## The mass of the node.
        self.mass = 0
        
        ## The center of mass of the node. 
        self.com = [0,0,0]
        
        ## The list of children in the node, i.e. subnodes.
        self.children = []
        
        ## The list of galaxies in the node.
        self.galaxies = []
        
        if parent is not None:
            self.parent = parent
        
        
    def insert_galaxy(self, galaxy):
        '''
        Summary:
        Inserts a galaxy into the node object.
        
        Parameters
        ----------
        galaxy: the galaxy we are adding to the node object.
        '''
        
        ## Append the individual galaxy to the list of galaxies in the node.
        self.galaxies.append(galaxy)
        
    
    def insert_node(self, sub_node):
        '''
        Summary:
        Inserts a sub node into the parent node.
        
        Parameters
        ----------
        sub_node : the node object that will be added to the parent node.
        
        '''
              
        ## Append the child node to the list of nodes under the parent.
        self.children.append(sub_node)
        
        
    def center_mass_calc(self):
        '''
        Summary:
        
        Calculates the center of mass of the node object.
        '''
        
        ## If our node has children, we have to calculate the 
        ## center of mass of all the child nodes.
        if len(self.children) > 0:
            
            ## Loop through the sub nodes in the parent node.
            for sub_nodes in self.children:
                
                ## In each node, recursively calculate the center of mass
                ## of the node.
                sub_node_com = sub_nodes.center_mass_calc()
                
                ## Store the center of mass coordinates and mass
                ## of each sub_node. 
                x_sub_com = sub_node_com[0]
                y_sub_com = sub_node_com[1]
                z_sub_com = sub_node_com[2]
                M_sub = sub_node_com[3]
                
                ## Use the center of mass of each sub node to calculate
                ## the center of mass of the parent node.
                self.com[0] += x_sub_com * M_sub
                self.com[1] += y_sub_com * M_sub
                self.com[2] += z_sub_com * M_sub
                self.mass += M_sub
            
            ## Right now, self.com isn't truly the center of mass. We
            ## still need to divide by the total mass. 
            self.com[0] = self.com[0] / self.mass
            self.com[1] = self.com[1] / self.mass
            self.com[2] = self.com[2] / self.mass
            
            return self.com[0], self.com[1], self.com[2], self.mass
        
        ## If our node has no children, then 
        ## the center of mass, com, and coords, x, y, z, all will be 
        ## the properties of the single galaxy in the node or 0 if
        ## no galaxy is present.
        elif len(self.children) == 0:
            self.com[0] = self.galaxies[0].x
            self.com[1] = self.galaxies[0].y
            self.com[2] = self.galaxies[0].z
            
            return self.galaxies[0].x, self.galaxies[0].y, self.galaxies[0].z, self.galaxies[0].mass
        
class Tree:

    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max, galaxy_population):
        '''
        Summary:
        Initializes a tree structure that will serve as a global node
        within which all galaxies will reside. It will then be parsed
        with the node class.
        
        Parameters
        ----------
        xyz_min : the minimum range of the global tree
        xyz_max : the maximum range of the global tree
        galaxy_population : a list of galaxies and their coordinates.
    
        '''
        
        ## Described above
        self.x_min = x_min
        self.y_min = y_min
        self.z_min = z_min
        self.x_max = x_max
        self.y_max = y_max
        self.z_max = z_max
        
        # The length of the Tree. Must be a cube!
        self.l = x_max - x_min
        
        ## The coordinates of all galaxies in the population of galaxies.
        self.galaxies_x = galaxy_population[:,0]
        self.galaxies_y = galaxy_population[:,1]
        self.galaxies_z = galaxy_population[:,2]
        
        ## Initialize the global or root node. The side length parameter
        ## is defined here as x_max. Before we can use the assumption
        ## that the side length should be x_max, we should first make
        ## sure that the node will be square as it should be.
        if x_max != y_max or x_max != z_max or y_max != z_max:
            exit('You are trying to initialize a non-square node. This \
                 code will not support non-square nodes')
            
        self.root = Node(x_min, y_min, z_min, x_max)
    
    def make_subnodes(self, parent):
        '''
        Summary:
        Makes the tree structure reproduce into subnodes based on whether
        a galaxy is present in the node. 
        
        Parameters
        ----------
        parent: a Tree structure that will serve as the parent node, that we
                will make subnodes from.
    
        '''
        
        ## The minimum value of x in the tree structure.
        x_min = parent.x_min
        y_min = parent.y_min
        z_min = parent.z_min
        
        ## The length of the node.
        half_length = parent.l / 2
        
        ## The midpoint of each node.
        midpoint_x = x_min + half_length
        midpoint_y = y_min + half_length
        midpoint_z = z_min + half_length
        
        ## in the cube, I will refer to the x axis as going from left to right
        ## y axis in to out
        ## z axis up to down.
        
        ## so we will have 8 quadrants, denoted 
        ## l_i_d
        ## l_i_u
        ## r_i_d
        ## r_i_u
        ## l_o_d
        ## l_o_u
        ## r_o_d
        ## r_o_u
        
        l_i_d = 0 #Node(x, y, z, half_length)
        
        for i in range(len()) 
        
        l_i_u = 0 #Node(x, y, z+half_length, half_length)
        r_i_d = 0 #Node(x+half_length, y, z, half_length)
        r_i_u = 0 #Node(x+half_length, y, z+half_length, half_length)
        l_o_d = 0 #Node(x, y+half_length, z, half_length)
        l_o_u = 0 #Node(x, y+half_length, z+half_length, half_length)
        r_o_d = 0 #Node(x+half_length, y+half_length, z, half_length)
        r_o_u = 0 #Node(x+half_length, y+half_length, z+half_length, half_length)
        
        
        
    
        
        