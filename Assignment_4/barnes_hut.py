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
        self.x = x_coord
        self.y = y_coord
        self.z = z_coord
        self.mass = mass
        self.coords = np.array([self.x, self.y, self.z])
        
        
        
        
class Node:
    
    def __init__(self, x, y, z, side_len):
        
        self.x = x
        self.y = y
        self.z = y
        self.com = None
        self.children = []
        self.l = side_len
        self.galaxies = []
        
    def insert_galaxy(self, galaxy):
        if len(self.galaxies) == 0:
            self.galaxies.append(galaxy)
            
        elif len(self.galaxies) > 0:
            self.galaxies.append(galaxy)
            
            if len(self.children) == 0:
                self.make_subnodes()
                
            for i in range(len(self.children)):
                
                if 
                
            
        
    def galaxy_presence(self, galaxy_coords):
        
        x_pos = galaxy_coords[0]
        y_pos = galaxy_coords[1]
        z_pos = galaxy_coords[2]
        
        if x_pos < self.x + self.side_len and x_pos > self.x:
            if y_pos < self.y + self.side_len and y_pos > self.y:
                if z_pos < self.z + self.side_len and z_pos > self.z:
                    return bool(1)
                
        else:
            return bool(0)
            
        
    def com_calc(self):
        if len(self.galaxies) == 0:
            self.com = [0,0,0]
            self.node_mass = 0
            return 0 ###Just need to exit the function here
            
        elif len(self.galaxies) == 1:
            self.com = self.galaxies[0].coords
            self.node_mass = self.galaxies[0].mass
            return 0
        
        com = [0,0,0]
        mass = 0
        
        for i in self.galaxies:
            if i.com == None:
                i.com_calc()
            com += i.com * i.node_mass
            mass += i.node_mass
        
        self.com = com * (1 / mass)
        self.node_mass = mass
        
    
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
        
    
        
        