/** graphics.h
Author: Andrea Ferrario

Header file containing the prototypes for the functions
that control the graphics.

*/
#pragma once

/** Initialize the graphics, opening the window.
particleData is a pointer to an array of length 3*particleCount, containing the three-dimensional
coordinates of all particles.
*/
void initializeGraphics(const float* const particleData, const int particleCount);

/** Runs the graphics.
This function blocks until the window is closed.
*/
void runGraphics();
