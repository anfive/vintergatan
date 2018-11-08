/** settings.h
Author: Andrea Ferrario

Global settings for the simulation.

*/

#pragma once

#include <chrono>

/** Size of the simulation area, a cube from -BOX_SIZE/2 to BOX_SIZE/2. */
constexpr double BOX_SIZE = 2.0;

/** "Gravitational constant", a global constant that scales the computed
forces between particles. */
constexpr double G = 1.0;

/** Maximum number of octree levels allowed.
This also determines the particle collision distance, since two particles will collide
when they occupy the same finest-division octant.
*/
constexpr unsigned int MAX_OCTREE_DEPTH = 21; // Must be < 64 / 3, since it determines the size of the Morton index, which is a 64-bit integer.

/** Value of the "softening distance" used to mitigate the singularity in the gravitational force
between close particles.
*/
constexpr double SOFTENING_DISTANCE = 1.0 / 1024.0;

/** Interval at which performance information is written to the console.
*/
constexpr std::chrono::milliseconds PERFORMANCE_LOGGING_INTERVAL(1000);