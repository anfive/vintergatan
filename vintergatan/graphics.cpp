/** graphics.cpp
Author: Andrea Ferrario

Implementation of the functions for graphics initialization and
for the rendering loop.

*/

#include <exception>
#include <iostream>
#include "GL/glew.h"
#include "GLFW/glfw3.h"

#include "graphics.h"
#include "settings.h"

/** 3-dimensional coordinates of the particles.
Updated by the simulation algorithm. */
const float* particleData;

/** Number of particles. */
unsigned int particleCount;

/** The graphics window. */
GLFWwindow* window;

/** VAO for the simulation box. */
GLuint box;

/** Position of the cursor. */
double cursorX, cursorY;

/** Angles used for view rotation. */
GLfloat a = 25.0, b = 70.0;

/** Vertices of the simulation region. */
constexpr float simCube[] = {
  -BOX_SIZE / 2.0, -BOX_SIZE / 2.0, -BOX_SIZE / 2.0,
  -BOX_SIZE / 2.0,  BOX_SIZE / 2.0, -BOX_SIZE / 2.0,
  BOX_SIZE / 2.0, -BOX_SIZE / 2.0, -BOX_SIZE / 2.0,
  BOX_SIZE / 2.0,  BOX_SIZE / 2.0, -BOX_SIZE / 2.0,
  -BOX_SIZE / 2.0, -BOX_SIZE / 2.0, BOX_SIZE / 2.0,
  -BOX_SIZE / 2.0,  BOX_SIZE / 2.0, BOX_SIZE / 2.0,
  BOX_SIZE / 2.0, -BOX_SIZE / 2.0, BOX_SIZE / 2.0,
  BOX_SIZE / 2.0,  BOX_SIZE / 2.0, BOX_SIZE / 2.0,
};

/** Indices into simCube to draw the wireframe simulation cube as strip quads. */
constexpr GLuint simCubeIndices[] = {
  0, 1, 2, 3, 6, 7, 4, 5, 0, 1
};

void onCursorPositionChange(GLFWwindow* window, double x, double y) {
  if (glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED) {

    a += static_cast<GLfloat>(x - cursorX) / 10.0f;
    b += static_cast<GLfloat>(y - cursorY) / 10.0f;

    cursorX = x;
    cursorY = y;
  }
}

void onMouseButton(GLFWwindow* window, int button, int action, int mods) {
  if (button != GLFW_MOUSE_BUTTON_LEFT) {
    return;
  }

  if (action == GLFW_PRESS)
  {
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwGetCursorPos(window, &cursorX, &cursorY);
  }
  else {
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
  }
}

void initializeGraphics(const float* const particleData, const int particleCount) {
  ::particleCount = particleCount;
  ::particleData = particleData;

  // Initialize openGL

  if (!glfwInit()) {
    std::cerr << "GLFW initialization failed" << std::endl;
    throw std::exception();
  }

  window = glfwCreateWindow(1200, 1200, "Vintergatan", NULL, NULL);
  if (!window) {
    std::cerr << "Windows creation failed" << std::endl;
    glfwTerminate();
    throw std::exception();
  }

  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  if (glewInit() != GLEW_OK) {
    std::cerr << "GLEW initialization failed" << std::endl;
    glfwTerminate();
    throw std::exception();
  }

  const GLubyte* renderer = glGetString(GL_RENDERER);
  const GLubyte* version = glGetString(GL_VERSION);
  std::cout << "Starting rendering thread." << std::endl;
  std::cout << "Renderer: " << renderer << std::endl;
  std::cout << "OpenGL version supported: " << version << std::endl;

  // Register mouse event handlers
  glfwSetMouseButtonCallback(window, onMouseButton);
  glfwSetCursorPosCallback(window, onCursorPositionChange);

  // Create the box

  box = 0;
  glGenVertexArrays(1, &box);
  glBindVertexArray(box);
  
  //glEnableVertexAttribArray(box);

  GLuint boxBuffer = 0;
  glGenBuffers(1, &boxBuffer);
  glBindBuffer(GL_ARRAY_BUFFER, boxBuffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(simCube), simCube, GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  glEnableClientState(GL_VERTEX_ARRAY);

  glEnable(GL_DEPTH_CLAMP);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

}

void runGraphics() {

  while (!glfwWindowShouldClose(window)) {

    glClear(GL_COLOR_BUFFER_BIT);

    // Rotate the view

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glRotatef(b, 1.0, 0.0, 0.0);
    glRotatef(a, 0.0, 0.0, 1.0);

    // Draw the box

    glBindVertexArray(box);
    
    glDrawElements(GL_QUAD_STRIP, 10, GL_UNSIGNED_INT, simCubeIndices);
    glBindVertexArray(0);

    // Draw the particles

    glVertexPointer(3, GL_FLOAT, 3 * sizeof(float), particleData);
    glDrawArrays(GL_POINTS, 0, particleCount);

    glfwSwapBuffers(window);

    // Handle events
    glfwPollEvents();
  }

  glfwTerminate();
}
