#pragma once
#include "voro++.hh"
#include "cinder/gl/gl.h"
#include <vector>

/**
 * Helper functions for voro++
 */

struct Face
{
	Face() = default;
	Face(const struct Triangle& triangle);
	/// Get normal vector of the plane
	glm::vec3 getNormal() const;

	/// Find direction vectors that are inside the plane and are not parallel
	std::array<glm::vec3, 2> getDirectionVectors() const;

	/// Remove all duplicate vertices in the face
	void removeDuplicates();

	/// Order the vertices so that they are listed in CCW order when looking against normal direction
	void orient(const glm::vec3& normal);

	std::vector<glm::vec3> vertices;
};

/// Get all faces of a voronoi cell
std::vector<Face> getFaces(voro::voronoicell& cell);

/// Get all faces of a voronoi cell using the cells'edges (this avoids broken face data of voro++)
std::vector<Face> getFacesFromEdges(voro::voronoicell& cell);

ci::TriMesh meshFromFaces(const std::vector<Face>& faces);
