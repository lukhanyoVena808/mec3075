import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import os
import time
import logging

def create_output_directory(m: int,
                            n_steps: int,
                            direct_solve: bool) -> str:
    """
    Creates a directory for storing the output within the './results' directory. The directory in the results directory
    must use the following naming convention: 'm-{m}-n_steps-{n_steps}-direct-{direct_solve}'
    :param m: The number of nodes in either the x or y direction
    :param n_steps: The number of time steps
    :param direct_solve: Whether a direct or iterative linear solver is to be used
    :return: The output path
    """
    # TODO: Create name for simulation using prescribed naming convention
    # TODO: Create variable for output path
    # TODO: Create directory for storing results using prescribed naming convention (see os.makedirs() documentation)
    # TODO: return string of output path
    directory = f"m-{m}-n_steps-{n_steps}-direct-{direct_solve}"
    parent_dir="./results"
    path = os.path.join(parent_dir,directory)
    dir = os.mkdir(path)
    print(dir)
    return dir
#pathb=create_output_directory(100,100, False)
#print(pathb)

def create_and_save_grid(m: int,
                         out_path: str) -> (np.ndarray, np.ndarray, float):
    """
    Creates X and Y grids on the domain [0, 1]^2 with m nodes in a given row. The X and Y arrays are saved to the output
    directory as "X.npy" and "Y.npy", respectively.
    :param m: The number of nodes in a given row
    :param out_path: The path to save the X and Y arrays to
    :return: The X array, Y array, and h, that is the distance between neighbouring nodes.
    """
    # TODO create X and Y mesh grid arrays (look up the numpy.meshgrid documentation)
    # TODO save X and Y arrays
    X = np.linspace(1/(m+1),m*(m+1),m)
    Y = np.linspace(1/(m+1),m*(m+1),m)
    
    X_out_path = out_path + "\\X" 
    Y_out_path = out_path + "\\Y"
    X, Y = np.meshgrid(X,Y)
    
    np.save(out_path+"/x.npy", Y)
    np.save(out_path+"/Y.npy", Y)

    h = 1/(m+1)
    return X_out_path,Y_out_path,np.linalg.norm((int(X[0][0])-int(X[0][1])))
#    return X_out_path,Y_out_path,

# raise NotlmplementedError


    
    

 #   raise NotImplementedError


def get_T0(X: np.ndarray,
           Y: np.ndarray,
           a: float,
           b: float,
           c: float,
           d: float) -> np.ndarray:
    """
    Defines the initial temperature distribution in the plate.
    :param X: A meshgrid array (shape=[m, m]) of X positions
    :param Y: A meshgrid array (shape=[m, m]) of Y positions
    :param a: Student specific parameter
    :param b: Student specific parameter
    :param c: Student specific parameter
    :param d: Student specific parameter
    :return: The initial temperature distribution in the plate.
    """
    # TODO Complete as described above.
    T_Int = np.multiply(a,np.sin(np.multiply(b,X))) + np.multiply(c,np.cos(np.multiply(d,Y)))
    return T_Int



def create_logger(out_path: str) -> logging.Logger:
    """
    Create a Logger object for logging values of interest at each time step.
    The log file should be named info.log
    :param out_path: The path where the log file is to be stored
    :return: A Logger object
    """
    # TODO instantiate logger (see logging.getLogger() documentation)
    # TODO set logging level to logging.INFO
    # TODO create a logging.FileHandler object
    # TODO set logging.FileHandler object logging level to logging.INFO
    # TODO set format for file handler (see logging.FileHandler.setFormatter() documentation)
    # TODO add file handler to logger
    # TODO return logger
    
    log = logging.getLogger("mylogger")
    log.setLevel(level=logging.INFO)
    f =logging.FileHandler("handler.log")
    f.setLevel(level=logging.INFO)
    f.setFormatter(logging.Formatter)
    log.addHandler(f)
    return log

   


def create_K(m: int,
             hdt: float) -> sp.csc_matrix:
             
    k = sp.diags(4+hdt,-1,m,format='csc')
    for j in range(m):
        if ((k[j][0]==0)) and (k[j][m-1]==0):
            k[j][j+1]=0
            k[j+1][j]=0
    return k        
    
    """
    Creates the sparse matrix in the system of linear equations to be solved.
    :param m: The number of nodes in a given row of the gird.
    :param hdt: h^2/\\Delta t
    :return: A scipy sparse matrix
    """
    
    # TODO Create K using scipy.sparse.diags (see online documentation). Use a 'csc' format.
    # TODO use a for-loop to set appropriate entries on the first diagonal to zero
    raise NotImplementedError


def get_solver(K: sp.csc_matrix,
               direct_solve: bool) -> callable:
    """
    Returns a function for solving the linear system of equations represented by K. That is, one should be able to use
    this function as follows to find d that satisfies Kd=f:
        solver = get_solver(K, true)
        d = solver(f)
    The argument "direct_solve" determines if the resultant solver uses a direct method or an iterative method. In the
    case of the direct method scipy.sparse.linalg.factorized should be used and in the case of the iterative method
    scipy.sparse.linalg.cg should be used (see documentation).
    :param K: Sparse matrix for linear system of equations.
    :param direct_solve: If true, uses a direct solver. If false, uses an iterative solver.
    :return: A function that returns d that satisfies Kd=f
    """
    # TODO Use an if statement to implement the functionality described above
    if direct_solve == True:
        solver = sp.linalg.factorizated(K)
    elif direct_solve == False:
        solver = lambda f: sp.linalg.cg(K,f, atol=0)[0]

    return solver

    


def apply_time_step(t: float,
                    dt: float,
                    solver: callable,
                    step: int,
                    hdt: float,
                    logger: logging.Logger,
                    start_time: float,
                    out_path: str,
                    Tn: np.ndarray) -> (np.ndarray, float, int):
    """
    Determine the solution at the next time step (among other things).
    :param t: The time at the start of the time step
    :param dt: The time increment
    :param solver: The function used to solve the linear system of equations
    :param step: The number of step n
    :param hdt: h^2/\\Delta t
    :param logger: Self explanatory
    :param start_time: The time at which the solving procedure started
    :param out_path: The path to where results are saved
    :param Tn: The solution at step n
    :return: Tn1, t at time step n+1, n+1
    """
    # TODO determine Tn1, that is, the solution at the next time step
    k = 4 + hdt

    # kTn1 = hdt*Tn[n]
    A = sp.csc_matrix([k],dtype=int)
    B = sp.csc_matrix([hdt*Tn[n]], dtype=int)
    Tn1 = solver(A, B)
   
        
    # TODO determine the time at Tn1
    # TODO log the information for the current time step.
    #  Use the format "step: {step}\t t:{t} s\t T_max: {maximum temperature in Tn1}\t T_min: {minimum}\t compute time: {amount of time taken to compute Tn1}
    # TODO Save the solution Tn1 using the format "T-{step}.npy". The shape should be that of the meshgrids X, Y ([m,m])
    #  and the positions in the arrays should correspond.
    # TODO increment the step number
    
    # T_Int[m+1] = T_Int[m]-(T_Int[m+1]-T_Int[m])
    
    
    raise NotImplementedError


def save_summary(m: int,
                 n_steps_completed: int,
                 out_path: str,
                 T_max: float,
                 T_min: float,
                 a: float,
                 b: float,
                 c: float,
                 d: float) -> None:
    """
    Save a summary of the simulation a json file format (use Google to look this up if it is unfamiliar)
    :param m: The number of nodes in a given row
    :param n_steps_completed: The number of time steps completed
    :param out_path: The path to where the simulation information is stored
    :param T_max: The maximum temperature throughout the simulation
    :param T_min: The minimum temperature throught the sumulation
    :param a: Student specific parameter
    :param b: Student specific parameter
    :param c: Student specific parameter
    :param d: Student specific parameter
    :return:
    """
    
    # TODO Implement the function as described above. Use the provided example output as a guide.
    raise NotImplementedError


def solve(m: int,
          n_steps: int,
          direct_solve: bool,
          a: float = 1,
          b: float = 1,
          c: float = 1,
          d: float = 1,
          t_end: float = 0.4) -> float:
    """
    Solve the problem described in the project hand out
    :param m: The number of nodes in a given row of the mesh
    :param n_steps: The number of time steps to be used
    :param direct_solve: Whether a direct or iterative solver is used
    :param a: Student specific parameter
    :param b: Student specific parameter
    :param c: Student specific parameter
    :param d: Student specific parameter
    :param t_end: The end time of the simulation
    :return: The total compute time taken for this function to run
    """

    # TODO: Calculate time step size
    # TODO: create output directory
    # TODO: Create logger
    # TODO: Create mesh grid for node positions and save arrays in output directory
    # TODO: Create K (see scipy.sparse.diags documentation)
    # TODO: get solver
    # TODO: Assign initial T values
    # TODO: save initial T values to output directory as "T-0.npy"
    # TODO: Initialize t, step, T_max, T_min, and start_time variables
    # TODO use a loop that iteratively calls apply_time_step to calculate Tn1, t and step. Remember to update Tn in
    #  each iteration.
    # TODO: Save summary.json for simulation
    # TODO: Return total time taken to solve problem (start_time - time.time())
    raise NotImplementedError
