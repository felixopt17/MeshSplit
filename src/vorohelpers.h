#pragma once
#include "voro++.hh"
#include "cinder/gl/gl.h"
#include <vector>
#include <numeric>
#include <glm/gtx/hash.hpp>

/**
 * Helper functions for voro++
 */

struct Face
{
	Face() = default;
	Face(const struct Triangle& triangle);
	/// Get normal vector of the plane
	glm::vec3 getNormal() const;

	glm::vec3 calcCenterPoint() const
	{
		const glm::vec3 centerPoint = std::accumulate(vertices.begin(), vertices.end(), glm::vec3(0.f, 0.f, 0.f), std::plus<glm::vec3>())
			/ static_cast<float>(vertices.size());

		return centerPoint;
	}

	/// Get triangles making up the face. Make sure to orient the vertices in CCW order.
	std::vector<Triangle> triangulate() const;

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
std::vector<Face> getFacesFromEdges(voro::voronoicell& cell, const glm::vec3& particlePos);

ci::TriMesh meshFromFaces(const std::vector<Face>& faces);

