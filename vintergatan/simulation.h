/** simulation.h
Author: Andrea Ferrario

Header file containing declarations of data structures used for the simulation.

*/
#pragma once

#include <memory>
#include <cassert>

#include "settings.h"

/** A three-dimensional vector, stored as Cartesian coordinates.
Works pretty much as you would expect.
*/
struct Vector3 {
  /** Cartesian components. */
  double x, y, z;

  /** Constructor. */
  Vector3(const double x = 0.0,
    const double y = 0.0,
    const double z = 0.0)
    : x(x), y(y), z(z) {
  }

  /** Vector sum. */
  Vector3 operator+(const Vector3& other) const {
    return Vector3(x + other.x, y + other.y, z + other.z);
  }

  /** Vector subtraction. */
  Vector3 operator-(const Vector3& other) const {
    return Vector3(x - other.x, y - other.y, z - other.z);
  }

  /** Multiplication by a scalar. */
  Vector3 operator*(const double& scalar) const {
    return Vector3(x * scalar, y * scalar, z * scalar);
  }

  /** Division by a scalar. */
  Vector3 operator/(const double& scalar) const {
    return Vector3(x / scalar, y / scalar, z / scalar);
  }

  /** Returns the dot product of this vector with another vector. */
  double dot(const Vector3& other) const {
    return x * other.x + y * other.y + z * other.z;
  }

  /** Assignment operator. */
  Vector3& operator=(const double& scalar) {
    x = scalar;
    y = scalar;
    z = scalar;
    return *this;
  }

  /** Returns the square of the length of the vector. */
  double absSq() const {
    return dot(*this);
  }
};

/** Class for a Morton index.
A Morton index identifies an octant in the simulation space.
An octant is a cubic region of space obtained by recursively dividing the simulation cell in eight equal subcells.
The level of the Morton index is the number of recursive subdivisions, e.g.:
- a Morton index of level 0 (the highest level) identifies the undivided cell (so it has a single value);
- a Morton index of level 1 identifies one of the 8 octants (so it requires 3 bits),
- a Morton index of level 2 identifies one of the 64 octants obtained by further dividing each of the level 1 octants (6 bits),
etc.

The Morton index has the following properties:
1. Given the index of an octant at a certain level, the index of the octant at the higher (coarser) level containing the original
octant can be obtained by discarding the least significant 3 bits; in practice, each group of 3 bits in the index identifies
the octant at the corresponding level.
2. The three bits in each group (corresponding to a level) can be interpreted as the Cartesian coordinates of the identified octant;
for example:
- the octant at level 1 with index 101 is the octant with z > 0, y < 0, x > 0.
- the octant at level 2 with index 101001 is obtained by dividing the previous octant (101) and taking the one with smaller z and y and larger x.

The index is used in the simulation to efficiently construct the octree. From the coordinates of a point it is possible to efficiently
compute its Morton index at a desired level. For a given maximum length of the index, it is then possible to sort 3-dimensional points
according to their Morton index, as long as there are no two points in the same lowest level (finest) octant.

This class can store a Morton index at any level (from 0 to MAX_OCTREE_DEPTH), and provides utilities for manipulating indices.
The actual numerical value of the index is not accessible and kept private.
Indices are created from coordinate vectors and used only in comparison and sorting.

*/
class MortonIndex {

  /** The numerical value of the Morton index.
  The index is stored in the 3*level least significant bits.
  */
  size_t index;

  /** The level (of subdivision) this index corresponds to.
  level = 0 means no subdivisions (the simulation cell).
  The level value also indicates the number of 3-bits group in the index.
  */
  size_t level;

  /** Private constructor from the actual values. */
  explicit MortonIndex(const size_t index, const size_t level)
    : index(index), level(level) {
    // Verify at compile time that the size of the ::index member is enough to store the coordinates at the specified
    // maximum octree depth.
    static_assert(sizeof(decltype(index)) * 8 >= MAX_OCTREE_DEPTH * 3, "Size of MortonIndex::index is not sufficient for the speci");
  }

public:

  /** Default constructor; creates a Morton index at level 0. */
  MortonIndex() : index(0), level(0) {}

  /** Constructor from a spatial position vector.
  Creates a Morton index at the lowest (finest) level (i.e. MAX_OCTREE_DEPTH) for the specified coordinates.*/
  explicit MortonIndex(const Vector3& position) {

    // Construct the index.
    index = 0;

    // Loop through the levels starting from the coarsets.
    // For each level, determine the position of the point in the current subdivision and compute the 3-bit
    // group that identifies the subdivision.
    double octantSize = BOX_SIZE;
    Vector3 octantPosition(-BOX_SIZE / 2.0, -BOX_SIZE / 2.0, -BOX_SIZE / 2.0);
    for (auto level = 0; level < MAX_OCTREE_DEPTH; ++level) {

      // Shift the index to add a new 3-bits group.
      index <<= 3;
      // Divide the current octant and determine which one of the 8 subdivisions
      // contains the point.
      octantSize /= 2.0;
      if (position.x - octantPosition.x >= octantSize) {
        index |= 0x1; // First bit (LSB)
        octantPosition.x += octantSize;
      }
      if (position.y - octantPosition.y >= octantSize) {
        index |= 0x2; // Second bit
        octantPosition.y += octantSize;
      }
      if (position.z - octantPosition.z >= octantSize) {
        index |= 0x4; // Third bit
        octantPosition.z += octantSize;
      }
    }

    // indices created from coordinates are at the lowest (finest) level.
    level = MAX_OCTREE_DEPTH;
  }

  /** Use the default copy constructor */
  MortonIndex(const MortonIndex&) = default;

  /** Returns the level of this index. */
  size_t getLevel() const {
    return level;
  }

  /** Returns the index of the octant at the specified level that contains this octant.
  The specified otherLevel must be equal or higher (coarser = with a lower numeric value) than the
  level of this index. */
  MortonIndex getIndexAtLevel(const size_t otherLevel) const {
    assert(level >= otherLevel);
    // The index is computed easily by shifting of a number of groups
    // equal to the difference in levels.
    const size_t newIndex = index >> (3 * (level - otherLevel));
    return MortonIndex(newIndex, otherLevel);
  }

  /** Equality operator. It only makes sense to compare indices at the same level. */
  bool operator==(const MortonIndex& other) const {
    assert(level == other.level);
    return index == other.index;
  }

  /** Unequality operator. It only makes sense to compare indices at the same level. */
  bool operator!=(const MortonIndex& other) const {
    assert(level == other.level);
    return index != other.index;
  }

  /** Comparison operator. It only makes sense to compare indices at the same level. */
  bool operator<(const MortonIndex& other) const {
    assert(level == other.level);
    return index < other.index;
  }

};

/** Particle data.
The struct contains (almost) all the data for a particle used in the simulation.
Certain transient data used in time integration (acceleration and jerk) is stored separately for performance reasons.

The Particle data is stored in a vector shared among all threads, so the Particle struct is cache-aligned (at 64-bytes boundaries)
to avoid false sharing.
*/
struct alignas(64) Particle {

  /** The position of the particle at the last computed time step. */
  Vector3 pos;
  /** The velocity of the particle at the last computed time step. */
  Vector3 vel;

  /** The position of the particle at the next time step predicted by the time integrator. */
  Vector3 predictedPos;
  /** The velocity of the particle at the next time step predicted by the time integrator. */
  Vector3 predictedVel;

  /** The Morton index associated to the particle's predicted position. */
  MortonIndex index;

  /** The mass of the particle.
  Inactive particles are assigned a mass of NaN. */
  double mass;

  /** Flag controlling if the particle is active in the simulation.
  Particles are set to inactive when they collide with other particles or leave the simulation area.*/
  bool active;

  /** Computational weight of the particle in the last computed time step.
  The computational weight is the number of force calculations performed for the particle.
  This value is used to perform load balancing across threads in the force calculation part of the next step. */
  unsigned int computationalWeight;

  /** Default constructor.*/
  Particle() : mass(1.0), computationalWeight(1), active(true) {}

};

struct SimulationImpl;

/** Simulation object.
Contains the state of the simulation and runs it.
*/
class Simulation {

  /** Pointer to implementation. */
  std::unique_ptr<SimulationImpl> impl;

public:

  /** Constructor. Initializes the simulation.
  Arguments:

  threadCount     The number of simulation threads to use (at least one)
  renderingData   Pointer to the rendering data array shared with the graphics. Must be of size 3 * particleCount.
  At each time step the array is filled with the coordinates (x, y, z) of all particles. Coordinates of inactive particles
  are set to NAN.
  particleData    Pointer to the particle data. Must be of size particleCount.
  The caller must initialize the particle's positions and velocities (Particle::pos and Particle::vel). The other members of Particle
  are initialized in this constructor.
  deltaT          Time step to use for time integration.
  timeStepLimit   Number of time steps to compute per second (on average). A value of 0 means that the simulation will run at full speed.
  */
  Simulation(const unsigned int threadCount,
    float* const renderingData,
    Particle* const particleData,
    const unsigned int particleCount,
    const double deltaT,
    const unsigned int timeStepLimit);

  /** Destructor. */
  ~Simulation();

  /** Runs the simulation using the specified number of threads.
  This method blocks and continues running until the stop() method is called
  (from another thread). */
  void runParallel();

  /** Stops a running simulation. */
  void stop();

};

