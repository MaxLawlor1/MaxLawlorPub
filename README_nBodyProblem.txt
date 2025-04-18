# N-BODY PROBLEM SIMULATION USING RUNGE-KUTTA

The project approximates the N-body problem using user inputs and the 4th order Runge-Kutta method, along with STL and object-oriented programming. This problem is the prediction of the movement of N planets, or bodies, in space that are affected by each other's gravitational pull. The user should write an input file but has no interaction with the program during execution, and receives an output file with positional and velocity data at each timestep for the N bodies. 

# Files, Functions and Classes

Coursework.project: This is the corresponding project file for "main.ccp"
main.ccp: This file contains the C++ code including library inclusion; the class "Bodies" to store planetary position, acceleration and mass; and the following functions:
	initcons: Reads the input file "parameters.txt" to the "parameters" vector and throws an error if it is formatted incorrectly. Requires no inputs. Returns the "parameters" vector.
	sumr: This function checks whether the bodies will try to occupy the same position in space on the next timestep, and throws an error if so. Requires the number of bodies N, the current body's index i, Runge-Kutta constants kx and ky, the vector bodies in the class Bodies (containing positional and velocity data for each body), and the timestep deltat. There are no outputs returned.
	rk4: Integrates the 4th order Runge-Kutta equation and writes a results file "output.txt". Takes the vector bodies from the class Bodies, the gravitational constant  G, the total simulation runtime T, the timestep deltat, and the total number of bodies N. There are no outputs returned as the function writes the necessary information to the file "output.txt".
	main: The main function that calls the above functions as needed to run the N body problem simulation. Is responsible for looping through the N bodies present in the problem.

# Inputs

"parameters.txt" file: The first line is the gravitational constant G, the total simulation time T, and the timestep deltat separated by spaces (not commas). The next lines are the initial condititions for each body (x-position, y-position, x-velocity, y-velocity, mass). These values can be dimensionalised in any units, as long as they are consistent.

# Outputs

"outputs.txt" file: contains the body index i, time t, positions of the bodies xi and yi, and their velocity ˙xi and ˙yi of each body at each time step in the same format as the input file.

# Usage

The specified "parameters.txt" file in the same directory as "main.cpp" and "Coursework.project" should be appropriately formatted.
eg.
1.0 5.0 0.001
3.0 3.0 0.0 0.0 1.0

The project can then be compiled and run, and the user will receive the "output.txt" file which can then be used to plot velocity or position of the different bodies in a program of their choice.

# Max Lawlor, https://github.com/MaxLawlor1/MaxLawlorPub, 17/04/2025