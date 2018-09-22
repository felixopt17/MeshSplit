#pragma once
#include "voro++.hh"
#include "cinder/gl/gl.h"
#include <vector>
#include "vorohelpers.h"
#include "glm/gtc/epsilon.hpp"
#include "geometry.h"
#include <unordered_set>

using std::vector;
using ci::TriMesh;

typedef std::vector<Face> Cell;

void setTitleMessage(const std::string& str);


vector<TriMesh> splitMesh(const TriMesh& sourceMesh, const std::vector<Cell>& cells);

/**
 * Broad test for mesh - cell intersection
 */
bool cellIntersectsMesh(const TriMesh& mesh, const std::vector<Face>& cellFaces);

/**
 * Calculate intersection of a source mesh and a voronoi cell as a new mesh
 */
bool intersectMesh(const TriMesh& source, const std::vector<Face>& cellFaces, TriMesh& outMesh);


/**
 * Split triangle by a line segment inside the triangle. Return new triangles
 */
std::vector<Triangle> splitTriangleBySegment(const Triangle&triangle, const LineSegment& segment,const Plane& halfspace);


std::vector<Triangle> getTriangles(const TriMesh& mesh);

std::vector<cinder::TriMesh> testSplit(const class TriMesh& mesh);

/// Attempt to create a cap for the holecut by segments
std::vector<Triangle> createCap(const std::vector<struct OrientedLineSegment>& splitSegments, const Face& face, std::unordered_set<glm::vec3>& usedVertices);


std::vector<std::string> geogebraExport(const std::vector<OrientedLineSegment>& segments);