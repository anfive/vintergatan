/** octree.h
Author: Andrea Ferrario

Header file containing declarations of the octree data structure and
classes used for its manipulation.

*/
#pragma once

#include <vector>
#include <cassert>
#include "simulation.h"

/** Octree node.
An octree node represents an octant, i.e. a cubic recursive subdivision of the cubic simulation area.
The level of the node is the number of recursive subdivisions.
The octree node is identified by its Morton index (see documentation for MortonIndex in simulation.h).

Octree nodes can be "particle" nodes (in which case they have no children), or "parent" nodes, in which
case they have up to eight children nodes.
Each node also has a "sibling" pointer, that points to:
- the next child of this node's parent (in order of Morton indices), or
- the parent's sibling pointer, if this node's parent has no other children.
The sibling pointer is used to efficiently iterate through the tree when performing force calculations.

Since octree nodes are stored linearly in the octree (as a vector), it is not necessary to store the pointers
to the children, since this node's first child is the next node in the vector.
Similarly, the sibling pointer is just the index of the node in the vector.


*/
class OctreeNode {

  friend class Octree;
  friend class OctreeIterator;

  /** Morton index of this octree node. */
  MortonIndex index;
  
  /** Sibling pointer (linear position in the octree's array)*/
  size_t sibling;

  /** Constructor for a node.
  If a particle pointer is specified, this node is constructed as a particle node, and its center of
  mass position and velocity and the mass are taken from the particle. Otherwise it is considered
  a "parent" node and its values are initialized to zero (to be computed in Octree::computeCOMs()).
  */
  OctreeNode(const MortonIndex index, const double openingThreshold, const Particle* particle = nullptr);

  /** Returns the level of this node. */
  size_t getLevel() const {
    return index.getLevel();
  }

public :

  /** Pointer to a particle, if this is a particle node, or nullptr if this is a parent node. */
  const Particle* particle;

  /** Position of this node's center of mass.
  If this node is a particle node, this position corresponds to the position of the particle.
  Otherwise it is computed from the centers of mass of its children.
  */
  Vector3 centerOfMass;

  /** Velocity of this node's center of mass.
  If this node is a particle node, this velocity corresponds to the velocity of the particle.
  Otherwise it is computed from the velocities of the centers of mass of its children.
  */
  Vector3 centerOfMassVelocity;

  /** Mass of this node.
  If this node is a particle node, this is the particle's mass.
  Otherwise it is the total mass of all its children.
  */
  double mass;

  /** Node opening threshold. If the distance of a particle to the center of mass squared is larger
  than this value, the force calculation algorithm will use the node's center of mass and total mass to
  compute the force (reducing the number of force calculations).
  Otherwise, the algorithm will proceed with investigating the children. */
  double openingThreshold;

};

/** Barnes-Hut octree data structure.
The octree is used to organize all the particles in the simulation according to their position (i.e. using
their Morton index) and computing the center of mass of groups of particles in order to reduce the number of
force calculations required at each time step.

An octree is obtained by dividing recursively the simulation cell in octants and creating a tree structure
from this subdivision. The root node is the simulation cell and each node has eight children corresponding to
the eight octants in which the node can be subdivided. The recursion depth is called the "level" of the node.
This means that e.g. a node at level 3 is a cube of size 1/(2^3) = 1/8 of the simulation cell.

The idea behing the Barnes-Hut tree is to only store some of these nodes according to the following criteria:
- If a node contains a single particle, it is considered a "particle" node with no children.
- If a node contains multiple particle, it is considered a "parent" node and it will have one or more child nodes,
  as necessary, to contain all the particles.

When constructing the tree, the center of mass for each parent node is computed from its children nodes and stored.

Using a Barnes-Hut tree allows significantly speeding up the calculation of the forces between particles by only adding
the contribution of the center of mass of distant nodes (instead of the contributions of all particles under it).
The algorithm has a complexity of O(n*log(n)) in the number of particles (instead of O(n^2) of a direct sum).

* Octree construction and node Storage

The octree is constructed by adding one particle at a time. The particles must be added in order of their Morton index.
The nodes are stored linearly in a vector with the following criteria:
- The root node is at the position 0.
- The first child node immediately follows the parent node.
- Other children are added after the first child and all its descendants.

With these criteria:
- the operation of adding a particle to the node is an operation linear in the depth of the tree, since it requires just checking
  the particle's index against the last added node's index and possibly adding children to it;
- tree iteration takes advantage of the memory locality of the nodes stored in the vector.

When running parallelized, each thread constructs a different octree and the contributions of all octrees are added for each particle.
The Octree class is therefore cache-aligned at 64 bytes boundaries to avoid false sharing.

*/
class alignas(64) Octree {

  friend class OctreeIterator;

  /** Linear storage of all the octree nodes. */
  std::vector<OctreeNode> nodes;

  /** Indices + 1 of the last added node at each octree level,
  used to update the sibling pointers. */
  std::vector<size_t> lastNodesInLevels;
  
  /** Debug method that checks the consistency of the octree.
  Will throw assert failures if the structure is inconsistent. */
  bool verify() const;

  /** Copy constructor and assignment operators are deleted. */

  Octree(const Octree&) = delete;
  Octree& operator=(const Octree&) = delete;

public:

  /** Default constructor. */
  Octree();

  /** Move constructor. */
  Octree(Octree&& other);

  /** Clears the octree, removing all the nodes.
  Internally, this method does not deallocate the vectors, so the reserved memory is kept.
  Costructing a new octree with the same size will then require no allocations. */
  void clear();

  /** Adds a particle to the octree.
  Particles must be added in order of their Morton index */
  void addParticle(const Particle* p);

  /** Computes the centers of mass position and velocities and the total mass for all nodes.
  Must be called once the tree has been filled with particles. */
  void computeCOMs();

  /** Returns true if the octree is empty (has no nodes). */
  bool empty() const;

};

/** Class used to iterate through an octree during force calculation.

The getNextNode() method returns the next node in depth-first order.
The getSiblingNode() returns the sibling of the last returned node, i.e. the next child of the
node's parent, or the parent's sibling.

Both methods return nullptr if there are no more nodes.

The typical usage of this class is:
1. Get a node with getNextNode().
2. Determine if the node should be "opened", i.e. its children investigated.
3. If yes, go to 1.
4. If no, process the current node, call getSiblingNode() and go to 2.

*/
class OctreeIterator {

  /** Reference to the octree being iterated on. */
  const Octree& octree;
  /** Linear index of the current node in the node vector. */
  size_t currentNodeIndex;
  /** Linear index of the sibling of the current node in the node vector. */
  size_t currentSiblingIndex;

public :

  /** Constructor for an object that iterates on the given octree. */
  explicit OctreeIterator(const Octree& octree);

  /** Returns the next node in depth-first order.
  Returns nullptr if there are no more nodes. */
  const OctreeNode* getNextNode();

  /** Returns the sibling of the last returned node in depth-first order.
  Returns nullptr if there are no more nodes. */
  const OctreeNode* getSiblingNode();

};