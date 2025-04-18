# ORBITAL MECHANICS COURSEWORK

This single MATLAB file OrbitalMechanicsCoursework.m calculates and displays the answers to the 2024 coursework set at Imperial College London to identify the orbital mechanics in certain scenarios for the Magnetospheric Composition Explorer (MACEX) ESA mission, which consists of a spacecraft orbiting the Earth on a nominal science orbit. This is done using orbital relations as derived in lectures.

# Inputs

The user is not required to input any data, but the code uses the orbital elements of the MACEX as follows:
	a_1 = 15272 km
	e_1 = 0.3187
	i_1 = 88.39 degrees
	Omega_1 = 65.93 degrees
	omega_1 = 337.62 degrees
	theta_1 = 248.68 degrees

The code also uses initial conditions in the Earth-Centred Inertial (ECI) frame as follows for "Part 1":
	position, r_0 = [2031; 5599; 2775] km
	velocity, v_0 = [-0.8419; -3.204; 7.081] km/s 

And as a transer orbit for "Part 2":
	position, r_1 = [-4085.1; -9918; -11215.3] km
	velocity, v_1 = [1.911; 4.131; -2.135] km/s 
	transfer time Delta_t = 6.753 hours

# Outputs

The code outputs the following in the command window:

Part 1(i)
Position vector in perifocal frame, r_1,p
Velocity vector in perifocal frame, v_1,p


Part 1(ii)
Perifocal reference frame unit vectors:
i_e
i_w
i_p


Part 2(i)
Minimum-energy semi-major axis, a_m


Part 2(ii)
Transfer time for parabolic orbit, Delta_t_p
Transfer time provided, Delta_t
Identification of an elliptical transfer orbit


Part 2(iii)
Transfer time for minimum energy orbit, Delta_t_m


Part 2(iv)
Semi-major axis of transfer orbit, a_01


Part 2(v)
Departure velocity impulse vector, Delta_v_0
Arrival velocity impulse vector, Delta_v_1

# Usage

The user must only run the file OrbitalMechanicsCoursework.m in MATLAB


# Max Lawlor, https://github.com/MaxLawlor1/MaxLawlorPub, 17/04/2025