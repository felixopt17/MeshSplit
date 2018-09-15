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
	TriMesh capMesh;
	Face face;
	face.vertices = { vec3(-3.0f, 0.3f, 3.f), vec3(3.0f, 0.3f, 3.f),
		 vec3(3.0f, 0.3f, -3.f), vec3(-3.0f, 0.3f, -3.f)
	};
	face.orient(vec3(1,- 1, 0.2));

	//Plane plane(vec3(0.0f, 0.3f, 0.f), vec3(1, 1, 0.2));
	Plane plane(face);

	auto tris = getTriangles(mesh);
	std::vector<Triangle> clippedTris;
	std::vector<OrientedLineSegment> segments;
	for (auto& t : tris)
	{
		auto r = cutTriangleByPlane(t, plane, segments);

		for (const Triangle& tri : r)
		{
			const auto vCount = static_cast<uint32_t>(outMesh.getNumVertices());

			outMesh.appendPosition(tri.a);
			outMesh.appendPosition(tri.b);
			outMesh.appendPosition(tri.c);

			outMesh.appendTriangle(vCount, vCount + 1, vCount + 2);
		}
	}

	auto capTriangles =  createCap(segments, face);
	for (const Triangle& tri : capTriangles)
	{
		const auto vCount = static_cast<uint32_t>(capMesh.getNumVertices());

		capMesh.appendPosition(tri.a);
		capMesh.appendPosition(tri.b);
		capMesh.appendPosition(tri.c);

		capMesh.appendTriangle(vCount, vCount + 1, vCount + 2);
	}

	outMesh.recalculateNormals();
	capMesh.recalculateNormals();
	return { outMesh, capMesh };
}

/// Try to merge a new segment with the old ones, possibly creating a triangle to fill
void tryMergeSegments(const OrientedLineSegment& segment, std::vector<OrientedLineSegment>& segments, std::vector<Triangle>& triangles, const Face& face)
{
	const auto faceNormal = face.getNormal();
	// Calculate a plane that contains the segment with a normal pointing inside the fill shape
	const auto calcSegmentPlane = [&](const OrientedLineSegment& segment)
	{
		auto normal = glm::cross(faceNormal, segment.segment.line.direction);
		if (glm::dot(normal, segment.insideDir) < 0.f)
			normal = -normal;

		return Plane(segment.segment.getStart(), normal);
	};

	// Look for old segments with one common point and attempt to merge with them into triangle and a segment
	for (auto otherSegmentIt = segments.begin(); otherSegmentIt != segments.end(); ++otherSegmentIt)
	{
		// Note:
		// Succesful call to mergeSegments() invalidates the iterator!
		const auto mergeSegments = [&](const vec3& commonVertex, const vec3& firstOther, const vec3& secondOther)
		{
			if (fVecEquals(firstOther, secondOther))
			{
				// 0-size triangle, delete both segments
				segments.erase(otherSegmentIt);
				return true;
			}

			auto plane1 = calcSegmentPlane(segment);
			auto plane2 = calcSegmentPlane(*otherSegmentIt);

			const auto newSideMidpoint = (firstOther + secondOther) / 2.f;
			if (plane1.isInFrontStrict(newSideMidpoint) && plane2.isInFrontStrict(newSideMidpoint))
			{
				OrientedLineSegment newSegment(firstOther, secondOther, glm::normalize(newSideMidpoint - commonVertex));
				Triangle tri(commonVertex, firstOther, secondOther);

				//keep correct winding order
				if (glm::dot(tri.getNormal(), faceNormal) < 0)
					tri.changeWinding();

				//save triangle
				triangles.push_back(tri);

				//try to merge newly created segment
				segments.erase(otherSegmentIt);
				tryMergeSegments(newSegment, segments, triangles, face);

				return true;
			}

			return false;
		};

		if (fVecEquals(segment.segment.getStart(), otherSegmentIt->segment.getStart()))
		{
			if (mergeSegments(segment.segment.getStart(), segment.segment.getEnd(), otherSegmentIt->segment.getEnd()))
				return;
		}
		else
			if (fVecEquals(segment.segment.getStart(), otherSegmentIt->segment.getEnd()))
			{
				if (mergeSegments(segment.segment.getStart(), segment.segment.getEnd(), otherSegmentIt->segment.getStart()))
					return;
			}
			else
				if (fVecEquals(segment.segment.getEnd(), otherSegmentIt->segment.getStart()))
				{
					if (mergeSegments(segment.segment.getEnd(), segment.segment.getStart(), otherSegmentIt->segment.getStart()))
						return;
				}
				else
					if (fVecEquals(segment.segment.getEnd(), otherSegmentIt->segment.getEnd()))
					{
						if (mergeSegments(segment.segment.getEnd(), segment.segment.getStart(), otherSegmentIt->segment.getEnd()))
							return;
					}
	}

	// Could not be merged, save it for later
	segments.push_back(segment);
}

std::vector<Triangle> createCap(const std::vector<struct OrientedLineSegment>& splitSegments, const Face& face)
{
	const vec3 faceNormal = face.getNormal();
	std::vector<Triangle> triangles;
	std::vector<OrientedLineSegment> filteredSegments;
	std::vector<glm::vec3> newFaceVertices;

	// Filter away segments outside the face and cut crossing segments
	// Attempt to merge segments into triangles
	for (const auto s : splitSegments)
	{
		auto f = OrientedLineSegment(cutSegmentByFaceEdges(s.segment, face), s.insideDir);
		if (f.segment.getLength() > 0)
		{
			tryMergeSegments(f, filteredSegments, triangles, face);

			// Segments that get cut create new face vertices
			if (!fVecEquals(s.segment.getStart(), f.segment.getStart()))
				newFaceVertices.push_back(f.segment.getStart());

			if (!fVecEquals(s.segment.getEnd(), f.segment.getEnd()))
				newFaceVertices.push_back(f.segment.getEnd());
		}
	}

	Face allEdgeVerts(face);
	allEdgeVerts.vertices.insert(allEdgeVerts.vertices.end(), newFaceVertices.begin(), newFaceVertices.end());
	allEdgeVerts.orient(face.getNormal());
	const vec3 centerPoint = allEdgeVerts.calcCenterPoint();

	// Merge remaining segments with next face edge to create more triangles
	for (auto intersectingSegment = filteredSegments.begin(); intersectingSegment != filteredSegments.end();)
	{
		//find it's intersecting vertex
		const auto intersectionVert = std::find_if(allEdgeVerts.vertices.begin(), allEdgeVerts.vertices.end(), [&](const vec3& v)
		{
			return fVecEquals(intersectingSegment->segment.getStart(), v) || fVecEquals(intersectingSegment->segment.getEnd(), v);
		});


		if (intersectionVert != allEdgeVerts.vertices.end())
		{
			// Find next vertex that is inside the triangle
			const auto itToNext = intersectionVert + 1 != allEdgeVerts.vertices.end() ? intersectionVert + 1 : allEdgeVerts.vertices.begin();
			const vec3 directionForward = *itToNext - *intersectionVert;

			if (glm::dot(directionForward, intersectingSegment->insideDir))
			{
				OrientedLineSegment newSegment(*intersectionVert, *itToNext, centerPoint - *itToNext);
				newSegment.fixInsideDir(faceNormal);
				tryMergeSegments(newSegment, filteredSegments, triangles, face);
			}
			else
			{
				const auto itToPrev = intersectionVert == allEdgeVerts.vertices.begin() ? --allEdgeVerts.vertices.end() : intersectionVert - 1;
				OrientedLineSegment newSegment(*intersectionVert, *itToPrev, centerPoint - *itToPrev);
				newSegment.fixInsideDir(faceNormal);
				tryMergeSegments(newSegment, filteredSegments, triangles, face);
			}

			// Iterator was invalidated, need to start loop from a good point
			intersectingSegment = filteredSegments.begin();
		}
		else
		{
			// A segment that is not intersecting with the edge remained inside
			// this should not happen for properly sealed meshes
			//assert(false);
			++intersectingSegment;
			continue;
		}
	}


	//TODO conthere

	return triangles;
}

bool intersectMesh(const TriMesh& source, voro::voronoicell& cell, const glm::vec3& particlePos, TriMesh& outMesh)
{
	/*
	 * Intersect all triangles in source mesh with the cell
	 * Triangles that intersect are split, others that are fully inside the cell are kept whole
	 */

	auto triangles = getTriangles(source);
	std::vector<Triangle> capTriangles;

	auto faces = getFacesFromEdges(cell, particlePos);
	for (const auto& face : faces)
	{
		std::vector<Triangle> newTriangles;
		std::vector<OrientedLineSegment> splitSegments; //segments created by splitting a triangle in half
		const Plane plane(face);

		for (const Triangle& tri : triangles)
		{
			auto splitResult = cutTriangleByPlane(tri, plane, splitSegments);
			newTriangles.insert(std::end(newTriangles), std::begin(splitResult), std::end(splitResult));
		}

		auto capResult = createCap(splitSegments, face);
		capTriangles.insert(std::end(capTriangles), std::begin(capResult), std::end(capResult));



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
	for (const Triangle& tri : capTriangles)
	{
		const auto vCount = static_cast<uint32_t>(outMesh.getNumVertices());

		outMesh.appendPosition(tri.a);
		outMesh.appendPosition(tri.b);
		outMesh.appendPosition(tri.c);

		outMesh.appendTriangle(vCount, vCount + 1, vCount + 2);
	}

	return !triangles.empty();
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
