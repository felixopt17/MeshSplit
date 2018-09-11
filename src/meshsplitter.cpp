#include "vorohelpers.h"
#include "meshsplitter.h"
#include "glm/gtc/epsilon.hpp"

#ifndef _TEST_
#include <cinder/app/AppBase.h>
#endif


using std::vector;
using ci::TriMesh;
using glm::vec3;

void setTitleMessage(const std::string& str)
{
#ifndef _TEST_
	cinder::app::getWindow()->setTitle(str);
#endif
}


vector<TriMesh> splitMesh(const TriMesh& sourceMesh, voro::container& con)
{
	vector<TriMesh> result;

	// Loop over all cells and test for collision with the mesh
	voro::c_loop_all vLoop(con);
	vLoop.start();
	do
	{
		setTitleMessage("Splitting with cell " + std::to_string(vLoop.pid()) + "/" + std::to_string(con.total_particles()));

		voro::voronoicell vcell;
		double px, py, pz;
		vLoop.pos(px, py, pz);
		glm::vec3 particlePos(px, py, pz);

		if (con.compute_cell(vcell, vLoop))
		{
			if (cellIntersectsMesh(sourceMesh, vcell, particlePos))
			{
				TriMesh intersectionMesh;
				if (intersectMesh(sourceMesh, vcell, particlePos, intersectionMesh))
				{
					intersectionMesh.recalculateNormals();
					result.emplace_back(std::move(intersectionMesh));
				}
			}
		}

	} while (vLoop.inc());

	setTitleMessage("Splitting complete");
	return result;
}

cinder::AxisAlignedBox calculateCellBounds(voro::voronoicell& cell, const glm::vec3& particlePos)
{
	using glm::vec3;

	vec3 min(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	vec3 max(std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest());

	std::vector<double> vertPositions;
	cell.vertices(particlePos.x, particlePos.y, particlePos.z, vertPositions);

	for (int i = 0; i < cell.p; i++)
	{
		const glm::vec3 pos(vertPositions[3 * i], vertPositions[3 * i + 1], vertPositions[3 * i + 2]);

		min.x = glm::min(pos.x, min.x);
		min.y = glm::min(pos.y, min.y);
		min.z = glm::min(pos.z, min.z);

		max.x = glm::max(pos.x, max.x);
		max.y = glm::max(pos.y, max.y);
		max.z = glm::max(pos.z, max.z);

	}

	return cinder::AxisAlignedBox(min, max);
}



bool cellIntersectsMesh(const TriMesh& mesh, voro::voronoicell& cell, const glm::vec3& particlePos)
{
	auto bounds = mesh.calcBoundingBox();
	const auto cellBounds = calculateCellBounds(cell, particlePos);

	return bounds.intersects(cellBounds);
}

std::vector<Triangle> getTriangles(const TriMesh& mesh)
{
	std::vector<Triangle> triangles;
	const auto numTriangles = mesh.getNumTriangles();
	for (size_t tNum = 0; tNum < numTriangles; tNum++)
	{
		glm::vec3 a, b, c;
		mesh.getTriangleVertices(tNum, &a, &b, &c);
		triangles.emplace_back(Triangle(a, b, c));
	}

	return triangles;
}

std::vector<cinder::TriMesh> testSplit(const TriMesh& mesh)
{
	TriMesh outMesh;
	Plane plane(vec3(0.0f, 0.3f, 0.f), vec3(1, 1, 0.2));
	auto tris = getTriangles(mesh);
	std::vector<Triangle> clippedTris;
	for (auto& t : tris)
	{
		auto r = cutTriangleByPlane(t, plane);

		for (const Triangle& tri : r)
		{
			const auto vCount = static_cast<uint32_t>(outMesh.getNumVertices());

			outMesh.appendPosition(tri.a);
			outMesh.appendPosition(tri.b);
			outMesh.appendPosition(tri.c);

			outMesh.appendTriangle(vCount, vCount + 1, vCount + 2);
		}
	}

	outMesh.recalculateNormals();
	return { outMesh };
}

bool intersectMesh(const TriMesh& source, voro::voronoicell& cell, const glm::vec3& particlePos, TriMesh& outMesh)
{
	/*
	 * Intersect all triangles in source mesh with the cell
	 * Triangles that intersect are split, others that are fully inside the cell are kept whole
	 */

	auto triangles = getTriangles(source);

	auto faces = getFacesFromEdges(cell, particlePos);
	for (const auto& face : faces)
	{
		std::vector<Triangle> newTriangles;
		const Plane plane(face);
		for (const Triangle& tri : triangles)
		{
			auto splitResult = cutTriangleByPlane(tri, plane);
			newTriangles.insert(std::end(newTriangles), std::begin(splitResult), std::end(splitResult));
		}

		std::swap(triangles, newTriangles);
	}

	// Add resulting triangles to mesh
	for (const Triangle& tri : triangles)
	{
		const auto vCount = static_cast<uint32_t>(outMesh.getNumVertices());

		outMesh.appendPosition(tri.a);
		outMesh.appendPosition(tri.b);
		outMesh.appendPosition(tri.c);

		outMesh.appendTriangle(vCount, vCount + 1, vCount + 2);
	}

	return !triangles.empty();
}

/**
 * Cut triangle by face. Only keep the part of the triangle that is in front of the face
 */
std::vector<Triangle> cutTriangleByFace(const Triangle& triangle, const Face& face)
{
	const Plane triPlane(triangle);
	const Plane facePlane(face);

	// Triangles that are fully in front of the plane are kept unchanged
	if (facePlane.isInFrontStrict(triangle.a) && facePlane.isInFrontStrict(triangle.b) && facePlane.isInFrontStrict(triangle.c))
	{
		return { triangle };
	}

	// Find new triangles from intersection
	const auto intersection = triPlane.intersect(facePlane);
	if (intersection.isLine)
	{
		// Limit the intersection to the faces
		const auto faceSegment = cutSegmentByFaceEdges(LineSegment(intersection.line), face);
		const auto commonSegment = cutSegmentByFaceEdges(faceSegment, Face(triangle));

		if (commonSegment.getLength() > 0)
		{
			return splitTriangleBySegment(triangle, commonSegment, facePlane);
		}
		else
		{
			return {};
		}
	}

	//TODO is not a line, is plane
	return { triangle }; // bad temp "fix"
}


/// Get closest of two points
glm::vec3 getClosest(const glm::vec3& source, const glm::vec3& a, const::glm::vec3& b)
{
	return glm::distance(source, a) > glm::distance(source, b) ? b : a;
}

/**
 * Split triangle by a segment. Keep original vertices that are in the correct halfspace
 */
std::vector<Triangle> splitTriangleBySegment(const Triangle& triangle, const LineSegment& segment, const Plane& halfspace)
{
	std::vector<glm::vec3> verticesToKeep;
	if (halfspace.isInFrontStrict(triangle.a)) verticesToKeep.push_back(triangle.a);
	if (halfspace.isInFrontStrict(triangle.b)) verticesToKeep.push_back(triangle.b);
	if (halfspace.isInFrontStrict(triangle.c)) verticesToKeep.push_back(triangle.c);

	assert(!verticesToKeep.empty());

	if (verticesToKeep.size() == 1)
	{
		Triangle result(verticesToKeep[0], segment.getStart(), segment.getEnd());

		// Keep triangle pointing to the same direction
		if (glm::dot(result.getNormal(), triangle.getNormal()) < 0)
			result.changeWinding();

		return { result };
	}

	if (verticesToKeep.size() == 2)
	{
		Triangle first(verticesToKeep[0], segment.getStart(), segment.getEnd());
		Triangle second(verticesToKeep[0], verticesToKeep[1], getClosest(verticesToKeep[1], segment.getStart(), segment.getEnd()));

		// Keep triangles pointing to the same direction
		if (glm::dot(first.getNormal(), triangle.getNormal()) < 0)
			first.changeWinding();

		if (glm::dot(second.getNormal(), triangle.getNormal()) < 0)
			second.changeWinding();

		return { first, second };
	}

	return {}; //this should not happen
}
