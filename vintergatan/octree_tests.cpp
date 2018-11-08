/* octree_tests.cpp
Author: Andrea Ferrario

Some quick and simple tests for the octree data structure and its methods.
They mostly rely on Octree's internal verify() method to check the results
(or I ran them with the debugger and checked manually - no time to write better checks...)

*/
#include <iostream>
#include <algorithm>

#include "octree.h"

#ifdef _MSC_VER

void test1() {

  Octree o;

  std::vector<Particle> p(4);
  p[0].pos = Vector3(-.5, -.5, -1.0);
  p[1].pos = Vector3(.5, -.5, -1.0);
  p[2].pos = Vector3(-.5, .5, -1.0);
  p[3].pos = Vector3(.5, .5, -1.0);

  for (int i = 0; i < p.size(); ++i) {
    p[i].index = MortonIndex(p[i].pos);
  }

  std::sort(p.begin(), p.end(), [](const Particle& p1, const Particle& p2) { return p1.index < p2.index; });

  for (int i = 0; i < p.size(); ++i) {
    o.addParticle(&p[i]);
  }

}

void test2() {

  auto mi0 = MortonIndex(Vector3(-1.0, -1.0, -1.0));
  auto mi1 = MortonIndex(Vector3(1.0, 1.0, 1.0));
  auto mi2 = MortonIndex(Vector3(0.0, 0.0, 0.0));
  auto mi3 = MortonIndex(Vector3(0.001, 0.0, 0.0));
  auto mi4 = MortonIndex(Vector3(0.0, 0.001, 0.0));
  auto mi5 = MortonIndex(Vector3(0.0, 0.0, 0.001));
  auto mi6 = MortonIndex(Vector3(-1.0+0.001, -1.0, -1.0));
  auto mi7 = MortonIndex(Vector3(-1.0, -1.0+0.001, -1.0));
  auto mi8 = MortonIndex(Vector3(-1.0, - 1.0, -1.0+0.001));

  auto mi2_1 = mi2.getIndexAtLevel(1);
  auto mi2_2 = mi2.getIndexAtLevel(2);
  assert(mi2_1 == mi2_2.getIndexAtLevel(1));

  auto eq1 = MortonIndex(Vector3(-1.0, -1.0, -1.0));
  auto eq2 = MortonIndex(Vector3(-1.0 + 1.999 / pow(2.0, MAX_OCTREE_DEPTH), -1.0, -1.0));
  assert(eq1 == eq2);

  eq1 = MortonIndex(Vector3(-1.0, -1.0, -1.0));
  eq2 = MortonIndex(Vector3(-1.0 + 2.0 / pow(2.0, MAX_OCTREE_DEPTH), -1.0, -1.0));
  assert(eq1 != eq2);

}

void test3() {

  Octree o;
  int ix = 0;
  
  std::vector<Particle> p(6);
  p[ix++].pos = Vector3(-.5, -.5, -1.0);
  p[ix++].pos = Vector3(.5, -.5, -1.0);
  p[ix++].pos = Vector3(-.5, .5, -1.0);
  p[ix++].pos = Vector3(-.5 + 0.0001, .5, -1.0);
  p[ix++].pos = Vector3(-.5 + 0.0002, .5, -1.0);
  p[ix++].pos = Vector3(.5, .5, -1.0);

  for (int i = 0; i < p.size(); ++i) {
    p[i].index = MortonIndex(p[i].pos);
  }

  std::sort(p.begin(), p.end(), [](const Particle& p1, const Particle& p2) { return p1.index < p2.index; });

  for (int i = 0; i < p.size(); ++i) {
    o.addParticle(&p[i]);
  }

  o.computeCOMs();

}

void test4() {

  Octree o;
  int ix = 0;

  std::vector<Particle> p(4);
  
  p[ix++].pos = Vector3(-1.0, -1.0, -1.0);
  p[ix++].pos = Vector3(-1.0 + 2.0 / pow(2.0, MAX_OCTREE_DEPTH), -1.0, -1.0);
  p[ix++].pos = Vector3(1.0, 1.0, 1.0);
  p[ix++].pos = Vector3(1.0 - 2.0001 / pow(2.0, MAX_OCTREE_DEPTH), 1.0, 1.0);

  for (int i = 0; i < p.size(); ++i) {
    p[i].index = MortonIndex(p[i].pos);
  }

  std::sort(p.begin(), p.end(), [](const Particle& p1, const Particle& p2) { return p1.index < p2.index; });

  for (int i = 0; i < p.size(); ++i) {
    o.addParticle(&p[i]);
  }

  o.computeCOMs();

}

void octreeTests() {

  std::cout << "Running octree unit tests" << std::endl;
  
  test1();
  test2();
  test3();
  test4();

  std::cout << "Done." << std::endl;

}

#else

void octreeTests() {
  std::cout << "Octree unit tests not available." << std::endl;
}

#endif