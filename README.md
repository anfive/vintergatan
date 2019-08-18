# Vintergatan
Multi-threaded octree-based multi-body gravitational dynamics simulator.
Author: Andrea "AnFive" Ferrario.

The program was written as a work sample in about ~40 hours. No validation or guarantee of any accuracy is provided - on the contrary, it most definitely does not handle well the dynamics of close particles due to the constant time step.

The program requires at least two threads, one for the simulation and one for the rendering.

## Build

Dependencies are included as submodules.

* On Windows, use CMake to create the Visual Studio solution.
The Windows build was tested with VS 2019 Community on Windows 10.
* On Linux (Ubuntu), install the dependencies:
   `apt install build-essential cmake libomp-dev xorg-dev libglu1-mesa-dev`
   Then use CMake with your favourite generator. The Ubuntu build was tested on Ubuntu 18.04 LTS.
