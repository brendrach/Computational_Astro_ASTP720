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
    
    def __init__(self, x, y, z, side_len):
        '''
        Summary:
        Initializes a node (cube) with corner (x,y,z) that extends in
        each direction a side length (side_len).
        
        Parameters
        ----------
        x, y, z : the coordinates of the corner of a node. 
        side_len : the length of the cube that extends in each direction
                   a distance, side_len.
        '''
        
        ## Can be thought of as x_min, y_min, z_min
        self.x = x
        self.y = y
        self.z = z
        
        ## The mass of the node.
        self.mass = 0
        
        ## The center of mass of the node. 
        self.com = [0,0,0]
        
        ## The list of children in the node, i.e. subnodes.
        self.children = []
        
        ## The length of the node. This extends from x, y, z.
        self.l = side_len
        
        ## The list of galaxies in the node.
        self.galaxies = []
        
        
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
        
        ## All galaxies in the tree.
        self.galaxies = galaxy_population
        
        ### EVERYTHING BELOW THIS LINE NEEDS UPDATED!###
        ################################################
        ################################################
        ################################################
        ################################################
        ################################################
        ################################################
        
        
    
    def make_subnodes(self):
        
        x = self.x
        y = self.y
        z = self.z
        
        half_length = self.l / 2
        
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
        
        l_i_d = Node(x, y, z, half_length)
        l_i_u = Node(x, y, z+half_length, half_length)
        r_i_d = Node(x+half_length, y, z, half_length)
        r_i_u = Node(x+half_length, y, z+half_length, half_length)
        l_o_d = Node(x, y+half_length, z, half_length)
        l_o_u = Node(x, y+half_length, z+half_length, half_length)
        r_o_d = Node(x+half_length, y+half_length, z, half_length)
        r_o_u = Node(x+half_length, y+half_length, z+half_length, half_length)
        
        self.children = [l_i_d, l_i_u, r_i_d, r_i_u, l_o_d, l_o_u, r_o_d, r_o_u]
        
        for i in range(0, len(self.children)):
            if self.children[i].galaxy_presence(self.galaxies):
                self.children[i].add_galaxy(self.galaxies)
        
    
        
        