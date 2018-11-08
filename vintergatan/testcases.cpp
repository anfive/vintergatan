/** testcases.cpp
Author: Andrea Ferrario

Implementation of test cases for the simulation algorithm.
Individual test cases can be run with -testcase <n>.

*/
#include "testcases.h"
#define _USE_MATH_DEFINES
#include <math.h> 
#include <random>

bool setupTestCase(const unsigned int testCase, std::unique_ptr<Particle[]>& particleData, unsigned int& particleCount, double& deltaT) {

  switch (testCase) {
  case 0: { // Two galaxies

    constexpr int n1 = 600, n2 = 600;
    deltaT = 1e-6;
    particleCount = n1 + n2;
    particleData.reset(new Particle[particleCount]);
    std::default_random_engine e(0);
    std::uniform_real_distribution<double> d(0.0, 1.0);

    int i;
    {
      constexpr double a = 0.4;
      constexpr double b = 0.1;
      constexpr double v = 50.0;

      for (i = 0; i < n1; ++i) {

        const double r = 0.1 + a * d(e) * 0.99;
        const double rn = d(e);
        const double alpha = 2.0 * M_PI * d(e);

        particleData[i].pos.x = r * sin(alpha);
        particleData[i].pos.y = r * cos(alpha);
        particleData[i].pos.z = -0.4;
        
        particleData[i].vel.x += v * r * cos(alpha);
        particleData[i].vel.y += -v * r * sin(alpha);
        particleData[i].vel.z += 20.0;
        
      }
    }

    {
      constexpr double a = 0.4;
      constexpr double b = 0.2;
      constexpr double v = 50.0;
      for (; i < n1 + n2; ++i) {

        const double r = 0.1 + a * d(e) * 0.99;

        const double alpha = 2.0 * M_PI * d(e);

        particleData[i].pos.x = r * sin(alpha);
        particleData[i].pos.y = r * cos(alpha);
        particleData[i].pos.z = 0.4;
        
        particleData[i].vel.x = -v * r * cos(alpha);
        particleData[i].vel.y = v * r * sin(alpha);
        particleData[i].vel.z = -20.0;
        

      }
    }

  } break;
  case 1: { // Cube

    int side = 25;

    particleCount = side * side * side;
    deltaT = 1e-5;
    particleData.reset(new Particle[particleCount]);
    constexpr double l = 0.8;
    constexpr double v = 200.0;
    int c = 0;
    for (auto i = 0; i < side; ++i) {
      for (auto j = 0; j < side; ++j) {
        for (auto k = 0; k < side; ++k) {

          particleData[c].pos.x = -l / 2.0 + (l * i) / side;
          particleData[c].pos.y = -l / 2.0 + (l * j) / side;
          particleData[c].pos.z = -l / 2.0 + (l * k) / side;

          const double r = sqrt(particleData[c].pos.x * particleData[c].pos.x + particleData[c].pos.y * particleData[c].pos.y);
          particleData[c].vel.x = v * particleData[c].pos.y * r;
          particleData[c].vel.y = -v * particleData[c].pos.x * r;
          particleData[c].vel.z = 0.0;

          ++c;

        }
      }
    }

  } break;
  case 2: { 

    int side = 25;

    particleCount = 3;
    deltaT = 1e-5;
    particleData.reset(new Particle[particleCount]);
    
    particleData[0].pos = Vector3(0.0, 0.0, -0.4);
    particleData[0].vel = Vector3(-5.0, 0.0, 0.0);

    particleData[1].pos = Vector3(0.0, 0.0, -0.42);
    particleData[1].vel = Vector3(5.0, 0.0, 0.0);

    particleData[2].pos = Vector3(0.0, 0.0, 0.4);
    particleData[2].vel = Vector3(0.0, 0.0, 0.0);

  } break;
  default: {
    return false;
  }
  }
  return true;

}

