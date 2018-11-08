/** testcases.h
Author: Andrea Ferrario

Header file for the function to setup test cases for the simulation algorithm.

*/
#pragma once

#include <memory>

#include "simulation.h"

bool setupTestCase(const unsigned int testCase, std::unique_ptr<Particle[]>& particleData, unsigned int& particleCount, double& deltaT);
