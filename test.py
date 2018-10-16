#!/usr/bin/env python

from laplace import LaplaceSolver
import numpy as np
import unittest
import time


class TestSolution(unittest.TestCase):

    def test_top_bc(self):

        solver = LaplaceSolver(nx=4, ny=3)
        solver.set_boundary_condtion('top', lambda x,y: 10)
        solver.numba_solve(quiet=True)
        np.testing.assert_allclose(solver.get_solution(), np.array([[ 0.        ,  0.        ,  0.        ,  0.        ],
                                                                 [ 0.        ,  4.09090909,  4.09090909,  0.        ],
                                                                 [10.        , 10.        , 10.        , 10.        ]], dtype=np.double), atol=0.01)


    def test_left_bc(self):

        solver = LaplaceSolver(nx=4,ny=4)
        solver.set_boundary_condtion('left', lambda x,y: 7)
        solver.numba_solve(quiet=True)
        np.testing.assert_allclose(solver.get_solution(),  np.array([[7.   , 0.   , 0.   , 0.   ],
                                                                     [7.   , 2.625, 0.875, 0.   ],
                                                                     [7.   , 2.625, 0.875, 0.   ],
                                                                     [7.   , 0.   , 0.   , 0.   ]], dtype=np.double) , atol=0.01)


    def test_right_bc(self):

        solver = LaplaceSolver(nx=4,ny=3)
        solver.set_boundary_condtion('right', lambda x,y: 5)
        solver.numba_solve(quiet=True)
        np.testing.assert_allclose(solver.get_solution(), np.array([[0.        , 0.        , 0.        , 5.        ],
                                                                    [0.        , 0.12121212, 0.78787879, 5.        ],
                                                                    [0.        , 0.        , 0.        , 5.        ]], dtype=np.double), atol=0.01)


    def test_bottom_bc(self):
        solver = LaplaceSolver(nx=3,ny=3)
        solver.set_boundary_condtion('bottom', lambda x,y: 14)
        solver.numba_solve(quiet=True)
        np.testing.assert_allclose(solver.get_solution(), np.array([[14. , 14. , 14. ],
                                                                    [ 0. ,  3.5,  0. ],
                                                                    [ 0. ,  0. ,  0. ]], dtype=np.double), atol=0.01)
        
        
    def test_timing(self):
        
        solver = LaplaceSolver()
        solver.set_boundary_condtion('top', lambda x,y: 10)
        solver.set_boundary_condtion('bottom', lambda x,y: 10)
        start = time.time()
        solver.solve(quiet=True)
        end = time.time()
        t1 = end - start
        solver.reset()
        solver.set_boundary_condtion('top', lambda x,y: 10)
        solver.set_boundary_condtion('bottom', lambda x,y: 10)
        start = time.time()
        solver.numba_solve(quiet=True)
        end = time.time()
        t2 = end - start
        
        assert t1 / t2 > 10.0
        
        
        
if __name__ == '__main__':
        unittest.main()
