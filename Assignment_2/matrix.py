import numpy as np

class Matrix:
    
    def __init__(self , shape , elems = []):
        '''
        Summary:
        Initializes a matrix with a given shape. At the init
        phase, the user can pass elements if they would like. Else, 
        they will be zeros. 
        
        Parameters
        ----------
        shape : a given 2-dimensional shape of rows and columns. Given
                in the form shape = [int, int]
        elems : float, fill the matrix at initializion.
        '''
        
        ## If we don't fill the matrix, initialize it.
        if elems == []:
            self.elems = []
            
            ## Loop over each row and along each column. 
            
            ## Initialize the rows
            row = 0
            while row < shape[0]:
                
                ## Initialize the columns
                col = 0
                self.elems.append([])
                while col < shape[1]:
                    self.elems[-1].append(0)
                    
                    ## Count along the columns
                    col = col + 1
                ## Count along the rows
                row = row + 1
        
        ## If we pass it elements, it just keeps those elements.
        else:
            self.elems = elems
            
        self.shape = shape
       
        
        
    def add(self , other):
        '''
        Summary:
        Adds an array other to our initialized array, self. 
        
        Parameters
        ----------
        other : a Matrix of identical dimenions to self.
        '''
        
        output = Matrix(self.shape)
        
        ## Make sure the shapes are the same.
        if other.shape[0] != self.shape[0] and other.shape[1] != self.shape[1]:
             exit('The matrices you are trying to add are not the same size')
        
        ## Count along the rows and columns and add each element.
        for rows in range(len(self.elems)):
            for cols in range(len(self.elems[rows])):
                output.elems[rows][cols] = self.elems[rows][cols] + other.elems[rows][cols]
                
        return output
    
    
    
    def transpose(self):
        '''
        Summary:
        Computes the transpose of a matrix.
        '''
        
        output = Matrix(self.shape)
        
        ## Set each matrix[i,j] = output[j,i]
        ## Modeled after Eq. 6 in Michael Lam's ASTP720 Matrix Notes. 
        for rows in range(len(self.elems)):
            for cols in range(len(self.elems[rows])):
                output.elems[cols][rows] = self.elems[rows][cols]
        
        return output
    
    
    def residual(self, ith_row, jth_col):
        '''
        Summary: 
        Returns self with the ith_row and jth_column removed.
        Parameters
        ----------
        ith_row : int, the row you want removed.
        jth_col : int, the column you want removed.
        '''
        
        ## Initialize the new Matrix to be smaller than self.
        output = Matrix([self.shape[0]-1, self.shape[1]-1])
        
        ## Empty the elements so we can append to them.
        output.elems = []
        
        ## Remove the row and column
        for rows in range(len(self.elems)):
            
            ## Check if we're at the ith_row
                if rows == ith_row:  
                    ## Do nothing, jump to beginning
                    continue
                
                ## Add a new row
                output.elems.append([])
                
                for cols in range(len(self.elems)):
                    ## Check if we're at the ith_row
                    if cols == jth_col:
                        ## Do nothing
                        continue
                    
                    ## Output the value at self[row][cols] to the 
                    ## newly initalized row. 
                    output.elems[-1].append(self.elems[rows][cols])
                    
        return output
        

    def multiply_matrix(self, other):
        '''
        Summary:
        Computes the product of two matrices. Returns a matrix.
        
        Parameters
        ----------
        other : the matrix we will multiply self by.
        '''
        
        
        if len(self.elems[0]) != len(other.elems[1]):
            exit("This matrix multiplication is not valid, the number of rows \
                 in Matrix A must equal the number of columns in Matrix B. ")
        
        output = Matrix([self.shape[0], other.shape[1]], elems=[])
        
        ## Loop over the rows and columns. This loop is modeled after 
        ## Eq. 4 in Michael Lam's ASTP720 Matrix notes.
        for rows in range(len(output.elems)):
            for cols in range(len(output.elems[rows])):
                
                index_value = 0
                for k in range(0, len(self.elems[0])):
                    index_value = index_value + self.elems[rows][k] * other.elems[k][cols]                    
                
                output.elems[rows][cols] = index_value
        
        return output
    
    
    
    def multiply_constant(self, other):
        '''
        Summary: 
        Computes the product of our number and self. Returns a matrix.
        Parameters
        ----------
        other : float/int, a constant we will multiply our matrix with.
        '''
        
        output = Matrix(self.shape)
        
        ## Multiply each element by other.
        for rows in range(len(self.elems)):
            for cols in range(len(self.elems[rows])):
                output.elems[rows][cols] = self.elems[rows][cols] * other
                
        return output


    def determinant(self):
        '''
        Summary: 
        Computes the determinant of a matrix. Returns a value.
        Parameters
        ----------
        None.
        '''
        
        ## Initialize the output and the row we will use to calculate
        ## the determinant.
        output = 0
        starting_row = 0
        
        ## If the matrix is a 2x2 - we will calculate
        ## the determinant analytically.
        ## This was suggested by a source found here:
        ## https://integratedmlai.com/find-the-determinant-of-a-matrix-with-pure-python-without-numpy-or-scipy/
        if self.shape[0] == 2 and self.shape[1] == 2:
            output = self.elems[0][0] * self.elems[1][1] - self.elems[1][0] * self.elems[0][1]
            return output
        
        
        ## Compute the determinant 
        for cols in range(len(self.elems)):
            ## Calculate the residual
            residual = self.residual(starting_row, cols)
            
            ## Calculate the determinant of the residual.
            residual_determinant = residual.determinant()
            
            ## Calculate the cofactor.
            cofactor = (-1)**(int(starting_row+1+cols+1)) * residual_determinant
            
            ## Compute the determinant.
            output = output + cofactor * self.elems[starting_row][cols]
                
        return output
    
    
    def inverse(self):
        '''
        Summary: 
        Computes the inverse of a matrix. Returns a matrix.
        Parameters
        ----------
        None.
        '''
        
        ## Initialize the output matrix.
        output = Matrix(self.shape)
        
        ## Calculate the determinant
        determinant = self.determinant()
        
        for rows in range(len(self.elems)):
            for cols in range(len(self.elems[rows])):
                ## Calculate the transpose of the residual
                residual = self.transpose().residual(rows, cols)
                
                ## Calculate the determinant of the transposed residual.
                residual_determinant = residual.determinant()
            
                ## Calculate the tranposed cofactor.
                cofactor = (-1)**(int(rows+1+cols+1)) * residual_determinant
            
                ## Calculate the inverse of the matrix. 
                output.elems[rows][cols] = cofactor / determinant
                
        return output
                
    
    def trace(self):
        '''
        Summary: 
        Computes the trace of a matrix. Returns a value.
        Parameters
        ----------
        None.
        '''
    
        ## Initialize the output.
        output = 0
        
        for rows in range(len(self.elems)):
            for cols in range(len(self.elems[rows])):
                
                ## Sum the diagonal
                if rows == cols:
                    output = output + self.elems[rows][cols]
         
        return output
        
    
    def LU_decomposition(self):
        '''
        Summary: 
        Computes the LU decomposition. Returns two matrices - L and U.
        Parameters
        ----------
        None.
        '''
        
        ## Initialize the matrices.
        L = Matrix(self.shape)
        U = Matrix(self.shape)
        
        for rows in range(len(self.elems)):
            
            for cols in range(rows, len(self.elems[rows])):
                    
                ## Compute u_ij from Eq. 29 of Michael Lam's Matrix notes.
                temp_sum_U = 0
                k = 0
                while k < (rows):
                    temp_sum_U = temp_sum_U + L.elems[rows][k] * U.elems[k][cols]
                    k = k + 1
                U.elems[rows][cols] = (self.elems[rows][cols] - temp_sum_U)
                
            for cols in range(rows, len(self.elems[rows])):
                
                ## If we are on the diagonal, l_ii = 1.
                if rows == cols:
                    L.elems[rows][cols] = 1
                    continue
                
                ## Compute l_ij from Eq. 28 of Michael Lam's Matrix notes.
                temp_sum_L = 0
                k = 0
                while k < (rows):
                    temp_sum_L = temp_sum_L + L.elems[cols][k] * U.elems[k][rows]
                    k = k+1
                L.elems[cols][rows] = 1/U.elems[rows][rows] * (self.elems[cols][rows] - temp_sum_L)

                
        return L, U
    

    
    
    
    
    
    
    
    
    
    
    
    
    