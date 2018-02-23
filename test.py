#!/usr/bin/env python

from laplace import LaplaceSolver
import numpy as np


def test_top_bc():
    
    solver = LaplaceSolver(nx=4, ny=3)
    solver.set_boundary_condtion('top', lambda x,y: 10)
    solver.swig_solve(quiet=True)
    np.testing.assert_allclose(solver.get_solution(), np.array([[0.,0.,0.],[0.,0.,2.35294118],[2.35294118,0.,0.],[10.,10., 10.]]), atol=0.01) 

    
def test_left_bc():
    
    solver = LaplaceSolver(nx=4,ny=4)
    solver.set_boundary_condtion('left', lambda x,y: 7)
    solver.swig_solve(quiet=True)
    np.testing.assert_allclose(solver.get_solution(),  np.array([[7., 0., 0., 0.],[7., 2.625, 0.875, 0.], [7., 2.625, 0.875, 0.   ],[7., 0., 0., 0.]]) , atol=0.01)
    
    
def test_right_bc():
    
    solver = LaplaceSolver(nx=4,ny=3)
    solver.set_boundary_condtion('right', lambda x,y: 5)
    solver.swig_solve(quiet=True)
    np.testing.assert_allclose(solver.get_solution(), np.array([[0., 0., 5.],[0., 0., 0.30252101], [0.87394958, 0., 5.], [0., 0., 5.]]), atol=0.01)
    
    
def test_bottom_bc():
    solver = LaplaceSolver(nx=3,ny=3)
    solver.set_boundary_condtion('bottom', lambda x,y: 14)
    solver.swig_solve(quiet=True)
    np.testing.assert_allclose(solver.get_solution(), np.array([[14., 14., 14.], [0.,3.5,0.],[0.,0.,0.]]), atol=0.01)
