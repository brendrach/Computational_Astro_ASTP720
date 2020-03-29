import numpy as np
import astropy.constants as const
import astropy.units as u
from cosmolopy import constants

## The value of G in proper units. I've used
## the package cosmolopy which has this particular
## unit system built in. This is more convenient than
## converting astropy's G to the correct unit system.
G = constants.G_const_Mpc_Msun_s

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
    
    def __init__(self, x_min, y_min, z_min, x_max, y_max, z_max, parent, gals):
        '''
        Summary:
        Initializes a node (cube) with corner (x,y,z) that extends in
        each direction a side length (side_len).
        
        Parameters
        ----------
        xyz_min : the minimum range of the global tree
        xyz_max : the maximum range of the global tree
        parent : the parent node.
        gals : the number of galaxies the node has.
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
        
        self.parent = parent
            
        self.gals = gals
        
        
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
        
        l = parent.l
        
        ## Load the galaxies in the parent. 
        gal_x = self.galaxies_x 
        gal_y = self.galaxies_y
        gal_z = self.galaxies_z 
        
        ## We will be looping over the list of galaxies a lot 
        ## so it is useful to define this.
        len_gal = len(galaxy_x_coords)
        
        ## In the cube, I will denote the quantities in the following way: 
        ## x axis as going from left to right.
        ## y axis as going from in to out.
        ## z axis as going from down to up.
        ## If this is confusing, trying reading my description again
        ## while looking at a 3D cube. Good example here:
        ## https://www.researchgate.net/figure/Cube-ABCD-EFGH-on-3D-coordinate-axis-with-a-length-of-the-edge-in-4-units_fig1_323231999
        ## Therefore, a description like l_i_d means the quadrant of the 
        ## cube that is on the (x = left, y = inner, z = down) boundary. Or more
        ## aptly written, the sub-node on the lower left side of the parent.
        
        ## so we will have 8 quadrants, denoted 
        ## l_i_d
        ## l_i_u
        ## r_i_d
        ## r_i_u
        ## l_o_d
        ## l_o_u
        ## r_o_d
        ## r_o_u
        
        ## Initially, I spent some time trying to vectorize this. 
        ## I'm sure it could be done, but it seems the brute force
        ## method is the best way for me to attack this. 
        
        ## This function is quite verbose. I'll comment the first instance
        ## and leave the rest uncommented. The functionality is identical so
        ## this shouldn't cause too much confusion.
        
        ## Initially, this node has 0 galaxies.
        l_i_d = 0 
        for i in range(len_gal):
            ## Loop over the galaxies and place them into the 
            ## proper nodes based on the galaxies coordinates.
            if x_min <= gal_x[i] and gal_x[i] < x_min + l/2 and \
               y_min <= gal_y[i] and gal_y[i] < y_min + l/2 and \
               z_min <= gal_z[i] and gal_z[i] < z_min + l/2 :
                   l_i_d = l_i_d + 1   
        ## If there are galaxies in this node, we need to define it as a 
        ## node object so we can calculate important quantites, com, accel, etc...
        if l_i_d != 0:
            parent.insert_node(Node(x_min, y_min, z_min, x_min+l/2, y_min+l/2, z_min+l/2, parent=parent, children=l_i_d))
            parent.gals += l_i_d
            
            
        r_i_d = 0
        for i in range(len_gal):
            if x_min + l/2 <= gal_x[i] and gal_x[i] < x_min + l and \
               y_min <= gal_y[i] and gal_y[i] < y_min + l/2 and \
               z_min <= gal_z[i] and gal_z[i] < z_min + l/2 :
                   r_i_d = r_i_d + 1
        if r_i_d != 0:
            parent.insert_node(Node(x_min+l/2, y_min, z_min, x_min+l, y_min+l/2, z_min+l/2, parent=parent, children=r_i_d))
            parent.gals += r_i_d
            
           
        r_i_u = 0
        for i in range(len_gal):
            if x_min + l/2 <= gal_x[i] and gal_x[i] < x_min + l and \
               y_min <= gal_y[i] and gal_y[i] < y_min + l/2 and \
               z_min + l/2 <= gal_z[i] and gal_z[i] < z_min + l :
                   r_i_u = r_i_u + 1
        if r_i_u != 0:
            parent.insert_node(Node(x_min+l/2, y_min, z_min+l/2, x_min+l, y_min+l/2, z_min+l, parent=parent, children=r_i_u))
            parent.gals += r_i_u
        

        l_o_d = 0
        for i in range(len_gal):
            if x_min <= gal_x[i] and gal_x[i] < x_min + l/2 and \
               y_min + l/2 <= gal_y[i] and gal_y[i] < y_min + l and \
               z_min <= gal_z[i] and gal_z[i] < z_min + l/2 :
                   l_o_d = l_o_d + 1
        if l_o_d != 0:
            parent.insert_node(Node(x_min, y_min+l/2, z_min, x_min+l/2, y_min+l, z_min+l/2, parent=parent, children=l_o_d))
            parent.gals += l_o_d
            
            
        l_o_u = 0
        for i in range(len_gal):
            if x_min <= gal_x[i] and gal_x[i] < x_min + l/2 and \
               y_min + l/2 <= gal_y[i] and gal_y[i] < y_min + l and \
               z_min + l/2 <= gal_z[i] and gal_z[i] < z_min + l :
                   l_o_u = l_o_u + 1
        if l_o_u != 0:
            parent.insert_node(Node(x_min, y_min+l/2, z_min+l/2, x_min+l/2, y_min+l, z_min+l, parent=parent, children=l_o_u))
            parent.gals += l_o_u
            
        
        r_o_d = 0
        for i in range(len_gal):
            if x_min + l/2 <= gal_x[i] and gal_x[i] < x_min + l and \
               y_min + l/2 <= gal_y[i] and gal_y[i] < y_min + l and \
               z_min <= gal_z[i] and gal_z[i] < z_min + l/2 :
                   r_o_d = r_o_d + 1
        if r_o_d != 0:
            parent.insert_node(Node(x_min+l/2, y_min+l/2, z_min, x_min+l, y_min+l, z_min+l/2, parent=parent, children=r_o_d))
            parent.gals += r_o_d
        
        
        r_o_u = 0
        for i in range(len_gal):
            if x_min + l/2 <= gal_x[i] and gal_x[i] < x_min + l and \
               y_min + l/2 <= gal_y[i] and gal_y[i] < y_min + l and \
               z_min + l/2 <= gal_z[i] and gal_z[i] < z_min + l :
                   r_o_u = r_o_u + 1
        if r_o_u != 0:
            parent.insert_node(Node(x_min+l/2, y_min+l/2, z_min+l/2, x_min+l, y_min+l, z_min+l, parent=parent, children=r_o_u))
            parent.gals += r_o_u
            
            
            
            
        def update_acceleration(self, Node, galaxy):
            '''
            Summary:
            Takes in a node and a particular galaxy and calculates the gravitational
            acceleration of the galaxy due to the mass of the node.
            
            Parameters
            ----------
            node : A node object that may contain children or be a leaf.
            galaxy : the galaxy object that is being accelerated due to the mass
                     of the node object.
        
            '''
            
            ## Getting the previous acceleration of the galaxy.
            gal_accel = galaxy.accel
            
            ## If the node has no children, it must be a leaf with 1 galaxy.
            if len(Node.children) == 0:
                
                ## Calculate r_vector in each direction.
                r_x = Node.com[0] - galaxy.x
                r_y = Node.com[1] - galaxy.y
                r_z = Node.com[2] - galaxy.z
                
                ## Calculate the magnitude of r.
                r = np.sqrt(r_x**2.0 + r_y**2.0 + r_z**2.0)
                
                ## The mass of the node.
                M = Node.mass
                
                ## The threshold for force softening. I chose this to be 1 pc.
                eps = 0.001
                
                ## Implementing force softening as described in Eq. 13 of Michael Lam's
                ## N-body interactions notes.
                if r == 0:
                    gal_accel[0] += 0
                    gal_accel[1] += 0
                    gal_accel[2] += 0
                
                elif r < eps and r > 0:
                    gal_accel[0] += G*M*r_x/(r * (r**2 + eps**2))
                    gal_accel[1] += G*M*r_y/(r * (r**2 + eps**2))
                    gal_accel[2] += G*M*r_z/(r * (r**2 + eps**2))
                   
                elif r > eps:
                    eps = 0
                    gal_accel[0] += G*M*r_x/(r * (r**2 + eps**2))
                    gal_accel[1] += G*M*r_y/(r * (r**2 + eps**2))
                    gal_accel[2] += G*M*r_z/(r * (r**2 + eps**2))
                    
            
            ## If the node has children, we must either treat the node and all of its 
            ## children as a single mass point by calculating its center of mass. If the 
            ## node is close enough, we will calculate the individual masses of the galaxies
            ## in it or its sub nodes (children).
            elif len(Node.children) > 0:
                for sub_node in Node.children:
                    
                    ## Calculate r_vector in each direction.
                    r_x = sub_node.com[0] - galaxy.x
                    r_y = sub_node.com[1] - galaxy.y
                    r_z = sub_node.com[2] - galaxy.z
                    
                    ## Calculate the magnitude of r.
                    r = np.sqrt(r_x**2.0 + r_y**2.0 + r_z**2.0)
                    
                    ## The mass of the node.
                    M = sub_node.mass
                    
                    ## The threshold for force softening. I chose this to be 1 pc.
                    eps = 0.001
                    
                    ## The distance from the galaxy to the center of mass of the node.
                    ## Simply a redefinition of r to match the descriptions in 
                    ## the Barnes-Hut notes.
                    D = r
                    
                    ## The length of the node. 
                    L = sub_node.x_max - sub_node.x_min
                    
                    ## The limit of L/D. In other words, if D * theta_lim < L, 
                    ## we will consider individual members of the group when we 
                    ## perform our acceleration calculations. If not, we assume we 
                    ## are far enough away to calculate the node as a single mass point.
                    theta_lim = 0.5
                    if L/D > theta_lim:
                        self.update_acceleration(sub_node, galaxy)
                        
                    
                    else:
                        if r == 0:
                            gal_accel[0] += 0
                            gal_accel[1] += 0
                            gal_accel[2] += 0
                        
                        elif r < eps and r > 0:
                            gal_accel[0] += G*M*r_x/(r * (r**2 + eps**2))
                            gal_accel[1] += G*M*r_y/(r * (r**2 + eps**2))
                            gal_accel[2] += G*M*r_z/(r * (r**2 + eps**2))
                           
                        elif r > eps:
                            eps = 0
                            gal_accel[0] += G*M*r_x/(r * (r**2 + eps**2))
                            gal_accel[1] += G*M*r_y/(r * (r**2 + eps**2))
                            gal_accel[2] += G*M*r_z/(r * (r**2 + eps**2))
                        
                        