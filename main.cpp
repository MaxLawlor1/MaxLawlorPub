#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream> // include appropriate libraries
#include <vector>
using namespace std;

vector<float> initcons() {          // function to read file "parameters.txt"
  ifstream vfile("parameters.txt"); // start reading the file "parameters.txt"
  vector<float> parameters; // initialise a vector to store each number from the file
  float num; // initialise a temporary file to store each number in turn from the file

  if (vfile.good()) {
    while (true) { // read the file
      if (vfile.eof()) {
        break;
      }      // break the loop when the end of the file is reached
      else { // run this if the end of the file has not been reached
        vfile >> num;
        if (!vfile) {
          throw logic_error("Error in parameters.txt input -\nFile must only contain numbers"); // throw error in case file contains string characters
		} else {parameters.push_back(num);} // put each number in the variable "num" then add to the end of the vector "parameters"
      }
    }
    vfile.close(); // close the file
  } else {cout << "Error opening file" << endl;} // error message if the file cannot be opened

  return parameters; // return the vector of each number in the file
}

class Bodies { // create "Bodies" class to house data for each body
public:
  float x; // define data to be housed: location in space (x,y), acceleration (xd,yd) and mass (m)
  float y;
  float xd;
  float yd;
  float m;
  Bodies(float a, float b, float c, float d, float e);
  
  void sumr(int N, int i, float kx, float ky, vector<Bodies> bodies, float deltat, float &sumx, float &sumy) {
  float d;  // initialise variable for distance between bodies
  sumx = 0; // initialise sums as zero
  sumy = 0;

  for (int j = 0; j < N; j++) {     // loop through each body for comparison with body i
    if (i != j) { // check if we are trying to compare body i to itself
      d = sqrt(((bodies[i].x - bodies[j].x) * (bodies[i].x - bodies[j].x)) + (bodies[i].y - bodies[j].y) * (bodies[i].y - bodies[j].y));
	  // ^ calculate distance between body i and j
      if (d == 0) {
        throw logic_error("Error in simulation: bodies attempting to occupy the same spatial coordinate");} 
		// ^ throw an error if the bodies are on top of each other (to avoid dividing by zero)
      else {
        sumx += (bodies[j].m / (d * d * d)) * (bodies[j].x - bodies[i].x - (kx * deltat)); // sum for calculation of x acceleration
        sumy += (bodies[j].m / (d * d * d)) * (bodies[j].y - bodies[i].y - (ky * deltat)); // sum for calculation of y acceleration
      }
    } else {continue;} // continue with the loop if trying to find the effect the body has on itself
  }
}

void rk4(vector<Bodies> bodies, float G, float T, float deltat, int N) { // void function to integrate using Runge-Kutta 4th and write a file of results 
  ofstream outp("output.txt", ios::out | ios::trunc); // start writing the file "output.txt"
  if (outp.good()) { // check if file creation has succeeded
    float kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4, kxd1, kxd2, kxd3, kxd4, kyd1, kyd2, kyd3, kyd4, sumx, sumy;
	// ^ initialise float numbers to use
	float atemp[5*N]; // initialise float array for temporary values
    for (float t = 0; t <= T; t += deltat) { // loop for increasing time with step delta t
      for (int i = 0; i < N; i++) { // loop through each body
        outp << i + 1 << " " << t << " " << bodies[i].x << " " << bodies[i].y << " " << bodies[i].xd << " " << bodies[i].yd << endl;
		// ^ write data to file

        kx1 = deltat * bodies[i].xd; // calculate Runge-Kutta constants for x position
        kx2 = deltat * (bodies[i].xd + (0.5 * kx1));
        kx3 = deltat * (bodies[i].xd + (0.5 * kx2));
        kx4 = deltat * (bodies[i].xd + kx3);

        ky1 = deltat * bodies[i].yd; // calculate Runge-Kutta constants for y position
        ky2 = deltat * (bodies[i].yd + (0.5 * ky1));
        ky3 = deltat * (bodies[i].yd + (0.5 * ky2));
        ky4 = deltat * (bodies[i].yd + ky3);

        sumr(N, i, 0.0, 0.0, bodies, deltat, sumx, sumy);
		// ^ call function to calculate multiplier for first RKconstants for x'' and y''
        kxd1 = deltat * G * sumx; // calculate first Runge-Kutta constants for x'' and y''
        kyd1 = deltat * G * sumy;

        sumr(N, i, (0.5 * kxd1), (0.5 * kyd1), bodies, deltat, sumx, sumy);
		// ^ call function to calculate multiplier for second RK constants for x'' and y''
        kxd2 = deltat * G * sumx; // calculate second Runge-Kutta constants for x'' and y''
        kyd2 = deltat * G * sumy;

        sumr(N, i, (0.5 * kxd2), (0.5 * kyd2), bodies, deltat, sumx, sumy); 
		// ^ call function to calculate multiplier for third RK constants for x'' and y''
        kxd3 = deltat * G * sumx; // calculate third Runge-Kutta constants for x'' and y''
        kyd3 = deltat * G * sumy;

        sumr(N, i, kxd3, kyd3, bodies, deltat, sumx, sumy);
		// ^ call function to calculate multiplier for fourth RK constants for x'' and y''
        kxd4 = deltat * G * sumx; // calculate fourth Runge-Kutta constants for x'' and y''
        kyd4 = deltat * G * sumy;

        // calculate next terms in each data point for this body and store in a temporary array
		atemp[5*i] = (bodies[i].x + ((kx1 / 6.0) + (kx2 / 3.0) + (kx3 / 3.0) + (kx4 / 6.0)));
        atemp[(5*i)+1] = (bodies[i].y + ((ky1 / 6.0) + (ky2 / 3.0) + (ky3 / 3.0) + (ky4 / 6.0)));
		atemp[(5*i)+2] = (bodies[i].xd + ((kxd1 / 6.0) + (kxd2 / 3.0) + (kxd3 / 3.0) + (kxd4 / 6.0)));
        atemp[(5*i)+3] = (bodies[i].yd + ((kyd1 / 6.0) + (kyd2 / 3.0) + (kyd3 / 3.0) + (kyd4 / 6.0)));
		atemp[(5*i)+4] = (bodies[i].m);
		// these values will replace the old ones at the end of the timestep, once calculations are completed
      }
	bodies.erase(bodies.begin(),bodies.end()); // clear vector for new values
	for (int k = 0; k < 5*N; k+=5) { // replace old values in class vector at the end of the timestep
		bodies.push_back(Bodies(atemp[k], atemp[k+1], atemp[k+2], atemp[k+3], atemp[k+4]));
	}
    }

    outp.close(); // close file once writing is completed
  }
  else {cout << "Error writing file" << endl;} // return an error message if unable to write the file
}
};

Bodies::Bodies(float a, float b, float c, float d, float e) { // define class constructor to accept certain inputs to defined data
  x = a, y = b, xd = c, yd = d, m = e;
}



int main() {
  try {
    vector<float> parameters = initcons(); // store the values from "parameters.txt" in the vector parameters
    if ((parameters.size() - 3) % 5 != 0 ||
        (parameters.size() - 3) == 0) { // throw an error if the wrong number of parameters has been input
      throw logic_error("Error in parameters.txt input -\nWrong number of parameters input - file must contain G, T, dt and 5 data points for each body");
    } 
	else {
      float G = parameters[0], T = parameters[1], deltat = parameters[2]; // store constants in separate variables

      if (G <= 0 || T < 0 || deltat <= 0) { // throw an error if parameters are impossible values
        throw logic_error("Error in parameters.txt input -\nMust have:\nG > 0\nT >= 0\ndt > 0");
      } else if (deltat > T) { // throw an error if delta t > T
        throw logic_error("Error in parameters.txt input -\ndt must be less than T");
      } else {
        parameters.erase(parameters.begin(), parameters.begin() + 3);
		// ^ remove first 3 values so "parameters" is a vector of body parameters only
        int N = parameters.size() / 5; // find number of bodies to be considered

        vector<Bodies> body; // create a class type vector to store variables for each body
        body.reserve(N); // preallocate space for vector

        for (int i = 1; i <= N; i++) { // loop through vector "body" and insert each variable from our original vector
          body.push_back(Bodies(parameters[(5 * i) - 5], parameters[(5 * i) - 4], parameters[(5 * i) - 3], parameters[(5 * i) - 2], parameters[(5 * i) - 1]));
		  if (parameters[(5 * i) - 1] < 0) { // throw an error if mass is less than zero
			throw logic_error("Error in parameters.txt input -\nBodies must have positive mass");
		  }
		  else{continue;}
        }

        body[0].rk4(body, G, T, deltat, N); // call "rk4" function to integrate with Runge-Kutta 4th and write file
      }
    }
  } catch (const exception &e) { // catching errors
    cout << "Error:" << endl << e.what() << endl << "\n(output.txt file has been removed)" << endl; // error message output
    remove("output.txt"); // deletes ouput file that was run using the error
    return 1;
  }

  return 0;
}