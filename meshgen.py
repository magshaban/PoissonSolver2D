# The part of <possionSolver2D.py> 
# This file generates the mesh with defining a data structure for storing the mesh.
# This file under development, 
# It will be generated automatically in the next code's version based on Delaunay Triangulation.
# you can also use a third party file to generate the mesh.
# t is the connectivity matrix 
# e is the edge matrix 
# p is the nodes coordinates matrix. 
# points is an array that generates the coordinates of the nodes.  

import numpy as np 

def meshgen():
    #p = np.array(([0,1,2,2,0],[0,0,0,1,1]))
    #t =np.array(([1,1,2],[4,2,3],[5,4,4]))


    p = np.array([
        [-1,-0.5,0,0.5,1,-1,-0.5,0,0.5,1,-1,-0.5,0,
         0.5,1,-1,-0.5,0,0.5,1,-1,-0.5,0,0.5,1],
        [-1,-1,-1,-1,-1,-0.5,-0.5,-0.5,-0.5,-0.5,0,0,
         0,0,0,0.5,0.5,0.5,0.5,0.5,1,1,1,1,1]
        ])
    t = np.array([
        [1,2,2,2,3,4,4,4,
         6,6,12,12,13,14,14,14,
         11,17,17,13,13,19,19,15,
         16,16,17,18,18,18,19,20
         ],
        
        [2,7,8,3,4,9,9,5,
         11,12,7,7,8,8,9,15,
         12,12,12,12,14,14,14,14,
         21,17,18,22,23,19,20,24
         ],
        
        [6,6,7,8,8,8,10,10,
         12,7,8,13,14,9,10,10,
         16,16,18,18,18,18,20,20,
         22,22,22,23,24,24,24,25
         ]
        ])
    E = np.array([
                  [1,2],[2,3],[3,4],[4,5],
                  [5,10],[10,15],[15,20],[20,25],
                  [25,24],[24,23],[23,22],[22,21],
                  [21,16],[16,11],[11,6],[6,1]
                  ])
    e = E.T

    return p,t,e

def delpoints():
    points = np.array([[0,0],[1,0],[-1,0],[0,1],
                   [0,-1],[-1, -1], [1, -1], 
                   [1, 1], [-1, 1],[-0.5,1],
                   [0.5,1],[-1,0.5],[0,0.5],
                   [1,0.5],[-0.5,0],[0.5,0],
                   [-1,-0.5],[0,-0.5],[1,-0.5],
                   [-0.5,-1],[0.5,-1],[0.5,0.5],
                   [-0.5,-0.5],[0.5,-0.5],[-0.5,0.5] ])
    return points
