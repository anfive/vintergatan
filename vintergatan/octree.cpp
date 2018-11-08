/** octree.cpp
Author: Andrea Ferrario

Implementation of the Barnes-Hut octree data structure.

*/

#include <math.h>

#include "octree.h"


/** Numerical value for the sibling pointer that indicates that the current node has no siblings. */
constexpr size_t NO_SIBLING = std::numeric_limits<size_t>::max();

OctreeNode::OctreeNode(const MortonIndex index, const double openingThreshold, const Particle * particle)
  : index(index), openingThreshold(openingThreshold), sibling(NO_SIBLING),
  particle(particle),
  mass(particle ? particle->mass : 0.0),
  centerOfMass(particle ? particle->predictedPos * particle->mass : Vector3()),
  centerOfMassVelocity(particle ? particle->predictedVel * particle->mass : Vector3()) {

}

bool Octree::verify() const {

  // Check structure consistency.

  // Store parent of the current node at each level
  const OctreeNode* parents[MAX_OCTREE_DEPTH + 1];
  std::fill(parents, parents + MAX_OCTREE_DEPTH + 1, nullptr);
  
  for (auto i = 0; i < nodes.size(); ++i) {

    const auto& node = nodes[i];
    auto level = node.getLevel();
    
    if (level == 0) {
      // Only the root node can be at level 0.
      assert(i == 0);
    }
    else {
      
      auto parent = parents[level - 1];
      // All nodes except the root node must have a parent
      assert(parent);
      // The parent must be at one level higher (coarser)
      assert(level == parent->getLevel() + 1);
      // The Morton index of the parent must be consistent with the index of this node.
      assert(node.index.getIndexAtLevel(parent->getLevel()) == parent->index);
    }

    if (!node.particle) {
      // This is a parent node.
      // Store it for this level.
      parents[level] = &node;
      // There must be child nodes.
      assert(i < nodes.size() - 1);
      // The next node in the vector must be its first child.
      assert(nodes[i + 1].getLevel() == level + 1);
    }
    // Check the sibling pointer
    if (node.sibling != NO_SIBLING) {
      // Root node has no siblings
      assert(level > 0);

      auto parent = parents[level - 1];
      assert(parent);
      // The sibling must be another child of the parent or must be the parent's sibling.
      assert(nodes[node.sibling].getLevel() == level || (nodes[node.sibling].getLevel() < level && node.sibling == parent->sibling));
    }
  }
  
  return true; // Always return true; any errors are caught by the asserts.
}

Octree::Octree() {
  // Reserve some memory for the nodes.
  nodes.reserve(1024);
  lastNodesInLevels.resize(MAX_OCTREE_DEPTH + 1, 0);
}

Octree::Octree(Octree && other) 
  : nodes(std::move(other.nodes)), lastNodesInLevels(std::move(other.lastNodesInLevels)) {
}

void Octree::clear() {
  // Clear the nodes, retaining capacity.
  nodes.clear();
  std::fill(lastNodesInLevels.begin(), lastNodesInLevels.end(), 0);
}

void Octree::addParticle(const Particle * p) {

  if (nodes.empty()) {
    // This is the first particle being added.
    // Create the root node at the top level.

    const auto index = p->index.getIndexAtLevel(0);
    nodes.push_back(OctreeNode(index, BOX_SIZE * BOX_SIZE, p));
    return;
  }

  // Check the last added node to determine if
  // - the particle is in the same octant, or
  // - the particle is in another octant.

  auto& lastNode = nodes.back();
  assert(lastNode.particle); // The last node is the deepest child, so it must be a particle node.
  const auto lastNodeLevel = lastNode.getLevel();
  const auto lastNodeIndex = lastNode.index;

  const auto newParticleIndex = p->index;

  // Compare Morton indices at the level of the last added node.
  if (lastNodeIndex != newParticleIndex.getIndexAtLevel(lastNodeLevel)) {
    // The particle is in another octant.

    // Move to higher level until we find a common parent.
    auto currentLevel = lastNodeLevel;
    // Index of the ancestor of the last node at the current level.
    MortonIndex currentIndex;
    do {
      --currentLevel;
      assert(currentLevel + 1 != 0);
    } while (newParticleIndex.getIndexAtLevel(currentLevel) != lastNodeIndex.getIndexAtLevel(currentLevel));

    // We found a common parent.
    // Add the new node as a child of the common parent.

    const double openingThreshold = pow(BOX_SIZE / pow(2.0, currentLevel + 1), 2.0);
    const auto indexOfAddedNode = nodes.size();
    nodes.push_back(OctreeNode(newParticleIndex.getIndexAtLevel(currentLevel + 1), openingThreshold, p));

    // Update the sibling pointers of all the last nodes at each lower level to point
    // to the newly added node. NOTE: no need to update the values of the
    // lastNodesInLevels: they are overwritten when nodes of the appropriate level are added.

    for (auto level = lastNodeLevel; level >= currentLevel + 1; --level) {
      assert(nodes[lastNodesInLevels[level]].sibling == NO_SIBLING);
      nodes[lastNodesInLevels[level]].sibling = indexOfAddedNode;
    }

    // Set the index of the new node as the current parent at this level.
    lastNodesInLevels[currentLevel + 1] = indexOfAddedNode;

    // If debug asserts are enabled, perform a consistency check on the tree.
    // This significantly slows down this method, so asserts must be disabled
    // in production.
    assert(verify());
   
  }
  else {
    // The particle is in the same octant.

    // Transform the current last node into a parent node (set the particle pointer to null).

    auto otherParticle = lastNode.particle;
    lastNode.particle = nullptr;
    lastNode.mass = 0.0;

    // Particles must be added in increasing Morton index, and must have distict indices.
    assert(otherParticle->index < p->index);

    // Now move to lower levels until we find the finest common octant.

    auto currentLevel = lastNodeLevel + 1;
    const double openingThreshold = pow(BOX_SIZE / pow(2.0, currentLevel), 2.0);
    // Index of the new node for the last particle in the tree
    auto currentIndex1 = otherParticle->index.getIndexAtLevel(currentLevel);
    // Index of the new node for the new particle (being added).
    auto currentIndex2 = newParticleIndex.getIndexAtLevel(currentLevel);
    while (currentIndex1 == currentIndex2) {
      // The two particles are in the same octant for this level as well.
      assert(currentLevel < MAX_OCTREE_DEPTH + 1);

      // Create a parent node for the two particles
      lastNodesInLevels[currentLevel] = nodes.size();
      nodes.push_back(OctreeNode(currentIndex1, openingThreshold));

      // Proceed to the lower level
      ++currentLevel;
      currentIndex1 = otherParticle->index.getIndexAtLevel(currentLevel);
      currentIndex2 = newParticleIndex.getIndexAtLevel(currentLevel);
    }

    // Add two new sibling nodes for the two particles.
    nodes.push_back(OctreeNode(currentIndex1, openingThreshold, otherParticle));

    const auto indexOfAddedNode = nodes.size();
    lastNodesInLevels[currentLevel] = indexOfAddedNode;
    nodes.push_back(OctreeNode(currentIndex2, openingThreshold, p));
    nodes[indexOfAddedNode - 1].sibling = indexOfAddedNode;
  }

  // If debug asserts are enabled, perform a consistency check on the tree.
  // This significantly slows down this method, so asserts must be disabled
  // in production.
  assert(verify());

}

void Octree::computeCOMs() {

  if (nodes.empty()) {
    return;
  }

  // Create work arrays for the calculation.
  // These arrays collect the contributions to the center of mass positions and velocities
  // for all children at each level. Contributions are stored as a sum of vector * mass, so
  // the final value of position and velocity can be computed by just dividing by the mass.
  double totalMass[MAX_OCTREE_DEPTH + 1];
  std::fill(totalMass, totalMass + MAX_OCTREE_DEPTH + 1, 0.0);
  Vector3 centerOfMass[MAX_OCTREE_DEPTH + 1];
  std::fill(centerOfMass, centerOfMass + MAX_OCTREE_DEPTH + 1, Vector3());
  Vector3 centerOfMassVelocity[MAX_OCTREE_DEPTH + 1];
  std::fill(centerOfMassVelocity, centerOfMassVelocity + MAX_OCTREE_DEPTH + 1, Vector3());

  // Iterate on the nodes in reverse order (bottom-up).
  // For each node, compute its center of mass from its children's and add its
  // contribution to the parent's center of mass.
  // The node order ensures that all children of a node are processed together
  // at their own level before moving on to the children of another node, and before
  // processing the parent.

  // Note that centers of mass of particle nodes are automatically set when the nodes
  // are created. Here we just need to update the center of mass of the parent nodes.

  for (auto i = nodes.size(); i > 0; --i) {

    auto& node = nodes[i - 1];
    auto currentLevel = node.getLevel();

    if (!node.particle) {
      // This is a parent node.
      assert(currentLevel <= MAX_OCTREE_DEPTH);

      // Compute the center of mass position and velocity by dividing all contributions from the children
      // by the total mass.
      node.centerOfMass = centerOfMass[currentLevel + 1] / totalMass[currentLevel + 1];
      node.centerOfMassVelocity = centerOfMassVelocity[currentLevel + 1] / totalMass[currentLevel + 1];
      node.mass = totalMass[currentLevel + 1];

      // Reset the collected contributions (we will process children of another node in the next loop).
      centerOfMass[currentLevel + 1] = Vector3();
      centerOfMassVelocity[currentLevel + 1] = Vector3();
      totalMass[currentLevel + 1] = 0.0;

    }

    // Add this node's contributions to the parent level's data.
    totalMass[currentLevel] += node.mass;
    centerOfMass[currentLevel] = centerOfMass[currentLevel] + node.centerOfMass * node.mass;
    centerOfMassVelocity[currentLevel] = centerOfMassVelocity[currentLevel] + node.centerOfMassVelocity * node.mass;

  }

}

bool Octree::empty() const {
  return nodes.empty();
}

OctreeIterator::OctreeIterator(const Octree& octree)
  : octree(octree),
    currentNodeIndex(static_cast<size_t>(-1)) {
  // currentNodeIndex is initialized to a value that gives zero
  // once incremented by one.

}

const OctreeNode* OctreeIterator::getNextNode() {
  if (currentNodeIndex == octree.nodes.size() - 1) {
    // No more nodes.
    return nullptr;
  }
  // Move to the next node.
  ++currentNodeIndex;
  return &octree.nodes[currentNodeIndex];
}

const OctreeNode* OctreeIterator::getSiblingNode() {
  if (currentNodeIndex == octree.nodes.size() - 1) {
    // No more nodes.
    return nullptr;
  }
  // Get sibling pointer of the currnet node.
  const auto& currentNode = octree.nodes[currentNodeIndex];
  auto sibling = currentNode.sibling;
  if (sibling != NO_SIBLING) {
    // Move the current node pointer to the sibling and proceed.
    currentNodeIndex = sibling;
    return &octree.nodes[currentNodeIndex];
  }
  else {
    // The current node (nor any of its ancestor) has any sibling.
    // The iteration is concluded.
    currentNodeIndex = octree.nodes.size() - 1;
    return nullptr;
  }

}
