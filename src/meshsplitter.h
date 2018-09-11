#pragma once
#include "voro++.hh"
#include "cinder/gl/gl.h"
#include <vector>
#include "vorohelpers.h"
#include "glm/gtc/epsilon.hpp"
#include "geometry.h"

using std::vector;
using ci::TriMesh;

void setTitleMessage(const std::string& str);


vector<TriMesh> splitMesh(const TriMesh& sourceMesh, voro::container& container);

/**
 * Broad test for mesh - cell intersection
 */
bool cellIntersectsMesh(const TriMesh& mesh, voro::voronoicell& cell, const glm::vec3& particlePos);

/**
 * Calculate intersection of a source mesh and a voronoi cell as a new mesh
 */
bool intersectMesh(const TriMesh& source, voro::voronoicell& cell, const glm::vec3& particlePos, TriMesh& outMesh);


/**
 * Split triangle by a face. Return part of the triangle that is in front of the face
 */
std::vector<Triangle> cutTriangleByFace(const Triangle& triangle, const Face& face);


/**
 * Split triangle by a line segment inside the triangle. Return new triangles
 */
std::vector<Triangle> splitTriangleBySegment(const Triangle&triangle, const LineSegment& segment,const Plane& halfspace);


std::vector<Triangle> getTriangles(const TriMesh& mesh);

std::vector<cinder::TriMesh> testSplit(const class TriMesh& mesh);