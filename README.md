# PoissonSolver2D
This is a Poisson 2D equation Python solver using FEM.<br />

The main file is <PoissonSolver2D.py><br />

The Poisson eguation (with a Robin general boundry condition) given as following: 
 
    -\nabla.(a(xy)\nabla{u}) = f(x,y)         in \Omega
    -n.(a(xy)\nabla{u}) = k(u - g_D) - g_N    in \partial{\Omega}

Where a(x,y) > 0, f(x,y), g_D and g_N are given functions.


     
# Results: 
With

    The domain is [-1,1]^2
    f(x,y) = 0 
    a(x,y) = 1
    k = 1
    g_D sin(x+y)^2
    g_N =1 on x = -1  and 0 elsewhere
    
we have the approximated solution u_h: <br />

![result](/solution1.png)
<br />  
![result](/solution2.png)
<br />
The stiffness  matrix: <br />
![result55](/stiffmatrix.png)
<br />

which evaluated on the following mesh: <br /> 
![result55](/mesh.png)
<br />    

 

# Author:
 Maged Shaban <br />
 magshaban[at]gmail.com <br />


