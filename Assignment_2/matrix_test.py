import unittest
from matrix import Matrix
import numpy as np
import scipy

class TestMatrix(unittest.TestCase):

    def test_add(self):
        a = [[1,2,3],[3,5,6],[6,1,1]]
        b = [[4,2,5],[1,7,5],[3,9,2]]
        
        A = Matrix([3,3], elems = a)
        B = Matrix([3,3], elems = b)
        correct = np.add(a,b)
        self.assertEqual(A.add(B).elems[1][1] , correct[1][1])
        
        
    def test_multiply(self):
        a = [[1,2,3],[3,5,6],[6,1,1]]
        b = [[4,2,5],[1,7,5],[3,9,2]]
        
        A = Matrix([3,3], elems = a)
        B = Matrix([3,3], elems = b)
        correct = np.matmul(a,b)
        self.assertEqual(A.multiply_matrix(B).elems[1][1] , correct[1][1])
        
        
    def test_transpose(self):
        a = [[1,2,3],[3,5,6],[6,1,1]]
        A = Matrix([3,3], elems = a)
        
        correct = np.transpose(a)
        self.assertEqual(A.transpose().elems[1][1] , correct[1][1])
        
        
    def test_inverse(self):
        a = [[1,2,3],[3,5,6],[6,1,1]]
        A = Matrix([3,3], elems = a)
        
        correct = np.linalg.inv(a)
        self.assertAlmostEqual(A.inverse().elems[1][1], correct[1][1])
        
        
    def test_trace(self):
        a = [[1,2,3],[3,5,6],[6,1,1]]
        A = Matrix([3,3], elems = a)
        
        self.assertAlmostEqual(A.trace(), np.trace(a))
        
        
    def determinant_test(self):
        a = [[1,2,3],[3,5,6],[6,1,1]]
        A = Matrix([3,3], elems = a)
        
        self.assertAlmostEqual(A.determinant().elems[1][1], np.linalg.det(a)[1][1])
        
        
    def LU_decomp_test(self):
        a = [[1,2,3],[3,5,6],[6,1,1]]
        A = Matrix([3,3], elems = a)
        
        correct = scipy.linalg.lu(a)
        self.assertAlmostEqual(A.LU_decomposition()[0].elems[1][1], correct[1][1])
        
        
		
if __name__ == "__main__":
	unittest.main()