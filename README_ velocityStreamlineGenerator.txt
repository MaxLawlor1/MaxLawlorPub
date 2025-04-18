# STREAMLINE AND VECTOR GENERATOR FOR 2D NACA AEROFOIL

This set of functions uses a 2D panel code to discretise an aerofoil into a selectable number of panels, from which the velocity of the airflow in an area around the aerofoil can be determined. The panel code assumes the airflow to be inviscid, irrotational, incompressible and steady. 

If a NACA2412 aerofoil is input by the user, the program will run a test case instead, and output lift curve slopes (lift coefficient vs angle of attack) for the aerofoil using 50, 100 and 200 panels at 10 degrees angle of attack to demonstrate the convergence of this code to data computationally determined using the XFOIL program developed by MIT. 

# Functions

cdoublet.m
	Function to compute the u and v velocity components produced at an input point p with coordinates p=[xp,zp] by a constant strength doublet panel with input end points at p1=[x1,z1] and p2=[x2,z2] and unit strenth. Outputs are x-velocity u and y-velocity v.

liftco.m
	Function to calculate the lift coefficient cl and effective dynamic viscosity mu at each aerofoil panel. This is done by solving a linear matrix system of the panel positions and relative distances and enforcing the Kutta condition at the aerofoil trailing edge. Inputs are panel coordinates x and z, number of panels N, angle of attack AoA (degrees), free-stream velocity uinf (m/s). Outputs are lift coefficient cl and effective dynamic viscosity mu.

panelgen.m
	Function to generate and discretise the aerofoil into panels. Takes the aerofoil code NACA, number of panels N and angle of attack AoA (deg), and outputs 2 arrays listing the x and y coordinates of the panels xf and z. 

velocities.m
	Function to calculate the velocity at each flow gridpoint and apply the necessary angle of attack. Takes the panel midpoints x and z, dynamic viscosity mu, free-stream velocity uinf, number of panels N and angle of attack AoA, and outputs 4 arrays - the flow gridpoints x and y position (x2, z2) and the x and y velocity at each gridpoint (U, V)

# Inputs
	Free-stream velocity (in m/s, eg. 10)
	Aerofoil Panel Code (the 4-digit NACA aerofoil code, eg. 2412)
	Angle of attack (in degrees, eg. 5)
	Number of panels (integer, eg. 100)

# Outputs

A plot of the aerofoil and airflow showing 2D velocity streamlines. The plot title indicates the aerofoil used, the number of panels, the angle of attack in degrees and the free-stream velocity in m/s.

A plot of the aerofoil showing 2D velocity vectors at each flow gridpoint. A larger size of arrow indicates a larger velocity magnitude at that point. The plot title indicates the aerofoil used, the number of panels, the angle of attack in degrees and the free stream velocity in m/s.

Additional plots are generated if the user inputs the test case of NACA2412 aerofoil:
	Lift coefficient (C_L) vs angle of attack (alpha, deg) for 50 panels
	Lift coefficient (C_L) vs angle of attack (alpha, deg) for 100 panels
	Lift coefficient (C_L) vs angle of attack (alpha, deg) for 200 panels


# Usage

The user should run the file "RUNME streamline and vector generator.m" in MATLAB and input the requested information in the command window when prompted.
The (non-test case) prompts should appear as follows: 
	Free stream velocity = 
	Aerofoil Panel Code = 
	Angle of attack = 
	number of panels = 

Example inputs are:
	15
	2410
	5
	100

To activate the test case, the user should input "2412" when prompted to input "Aerofoil Panel Code = ".


# Max Lawlor, https://github.com/MaxLawlor1/MaxLawlorPub, 17/04/2025
