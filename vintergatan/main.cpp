/** main.cpp
Author: Andrea Ferrario

Entry point for the program (main function).
Initializes and runs the graphics and the simulation.

*/

#include <iostream>
#include <memory>
#include <thread>
#include <string>
#include <cstring>
#include "omp.h"

#include "graphics.h"
#include "simulation.h"
#include "testcases.h"

/** Prototype for octree unit tests defined in octree_test.cpp */
void octreeTests();

/** Prints usage information*/
void printUsage() {
  std::cout <<
    "Usage: vintergatan [options]" << std::endl <<
    "Options: " << std::endl <<
    "  -testcase <n>      Run test case n (default: 0)." << std::endl <<
    "  -nthreads <n>      Use n threads for the simulation (default: cores - 1)." << std::endl <<
    "  -tspslimit <n>     Limit time integration to n time steps per second (default: off)." << std::endl;
}

/** Checks if the command line argument in args at position i is equal to the specified string,
in which case it reads the value from the next command line argument.

Returns false if the command line arguments are erroneously specified. */
bool readArgument(const char* option, char** const args, int& i, const int nargs, unsigned int& value) {
  if (strcmp(args[i], option) == 0) {
    if (i == nargs - 1) {
      printUsage();
      return false;
    }
    else {
      ++i;
      try {
        value = std::stoi(args[i]);
        return true;
      }
      catch (...) {
        printUsage();
        return false;
      }
    }
  }
  return true;
}

/** Entry point for the program. */
int main(int nargs, char** args) {

  std::cout <<
    "vintergatan v.0.1 by Andrea Ferrario" << std::endl <<
    "Simple n-body gravitational simulator." << std::endl;

  // Default test case
  unsigned int testcase = 0;
  // Default number of threads
  unsigned int numThreads = std::thread::hardware_concurrency() - 1;
  // Default timestep limit (unlimited)
  unsigned int timeStepLimit = 0;

  // Parse command line arguments.

  if (nargs == 2) {
    if (strcmp(args[1], "-octreetests") == 0) {
      // Run octree unit tests.
      octreeTests();
      return 0;
    }
  }

  for (int i = 1; i < nargs; ++i) {
    if (!readArgument("-testcase", args, i, nargs, testcase))
      return -1;
    if (!readArgument("-nthreads", args, i, nargs, numThreads))
      return -1;
    if (!readArgument("-tspslimit", args, i, nargs, timeStepLimit))
      return -1;
  
  }

  // Validate number of threads.
  if (numThreads < 1) {
    printUsage();
    return -1;
  }

  // Set up simulation data for the specified test case.

  unsigned int particleCount;
  std::unique_ptr<Particle[]> particleData;
  double deltaT;
  if (!setupTestCase(testcase, particleData, particleCount, deltaT)) {
    std::cout << "Test case number " << testcase << " invalid." << std::endl;
    return -1;
  }

  // Initialize OpenMP. Use one thread for rendering and numThreads threads for the simulation.
  omp_set_num_threads(numThreads + 1);
  omp_set_nested(1);

  // Allocate rendering data, that is shared between the simulation and the graphics.
  std::unique_ptr<float[]> renderingData(new float[3 * particleCount]);
  std::fill(renderingData.get(), renderingData.get() + particleCount * 3, 0.0f);

  // Initialize the simulation.
  Simulation simulation(numThreads, renderingData.get(), particleData.get(), particleCount, deltaT, timeStepLimit);

  // Start the rendering thread and the simulation master thread.
#pragma omp parallel num_threads(2)
  {

    if (omp_get_thread_num() == 0) {
      // This is the rendering thread.
      // Initialize the graphics and open the window before starting the simulations
      initializeGraphics(renderingData.get(), particleCount);
    }

#pragma omp barrier

    if (omp_get_thread_num() == 0) {
      // This is the rendering thread.
      // Run the graphics. The following function call blocks until the window is closed.
      runGraphics();

      // Stop the simulation.
      simulation.stop();
    }
    else {

      // This is the simulation master thread.
      // Run the simulation. The following function call blocks until the simulation is
      // stopped when the rendering thread calls Simulation::stop().
      simulation.runParallel();

    }

  }

  return 0;
}