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


vector<TriMesh> splitMesh(const TriMesh& sourceMesh, const std::vector<Cell>& cells)
{
	vector<TriMesh> result;

	// Loop over all cells and test for collision with the mesh
	for (size_t i = 0; i < cells.size(); i++)
	{
		const Cell& cell = cells[i];
		setTitleMessage("Splitting with cell " + std::to_string(i) + "/" + std::to_string(cells.size()));

		TriMesh intersectionMesh;
		if (cellIntersectsMesh(sourceMesh, cell))
		{
			if(intersectMesh(sourceMesh, cell, intersectionMesh))
			{
				intersectionMesh.recalculateNormals();
				
			}
		}
		result.emplace_back(std::move(intersectionMesh));

	}

	setTitleMessage("Splitting complete");
	return result;
}

cinder::AxisAlignedBox calculateCellBounds(const std::vector<Face>& cellFaces)
{
	using glm::vec3;

	vec3 min(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	vec3 max(std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest());

	for (const Face& face : cellFaces)
		for (const vec3& pos : face.vertices)
		{
			min.x = glm::min(pos.x, min.x);
			min.y = glm::min(pos.y, min.y);
			min.z = glm::min(pos.z, min.z);

			max.x = glm::max(pos.x, max.x);
			max.y = glm::max(pos.y, max.y);
			max.z = glm::max(pos.z, max.z);
		}

	return cinder::AxisAlignedBox(min, max);
}



bool cellIntersectsMesh(const TriMesh& mesh, const std::vector<Face>& cellFaces)
{
	auto bounds = mesh.calcBoundingBox();
	const auto cellBounds = calculateCellBounds(cellFaces);

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
	face.vertices = { vec3(-3.0f, 0.3f, 3.f), vec3(0.3f, 0.3f, 3.f),
		 vec3(0.3f, 0.3f, -3.f), vec3(-3.0f, 0.3f, -3.f)
	};
	face.orient(vec3(1, -1, 0.2));

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

	std::unordered_set<vec3> usedVerts;
	auto capTriangles = createCap(segments, face, usedVerts);
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
				// 0-size triangle, delete both segments (this one occurs naturally from the last triangle edge)
				segments.erase(otherSegmentIt);
				return true;
			}

			auto plane1 = calcSegmentPlane(segment);
			auto plane2 = calcSegmentPlane(*otherSegmentIt);

			const auto newSideMidpoint = (firstOther + secondOther) / 2.f;
			if (plane1.isInFrontOrInside(newSideMidpoint) && plane2.isInFrontOrInside(newSideMidpoint))
			{
				OrientedLineSegment newSegment(firstOther, secondOther, glm::normalize(newSideMidpoint - commonVertex));
				newSegment.fixInsideDir(faceNormal);
				Triangle tri(commonVertex, firstOther, secondOther);

				//make sure that triangle's normal points outside the cell
				if (glm::dot(tri.getNormal(), faceNormal) > 0)
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
					if (mergeSegments(segment.segment.getEnd(), segment.segment.getStart(), otherSegmentIt->segment.getEnd()))
						return;
				}
				else
					if (fVecEquals(segment.segment.getEnd(), otherSegmentIt->segment.getEnd()))
					{
						if (mergeSegments(segment.segment.getEnd(), segment.segment.getStart(), otherSegmentIt->segment.getStart()))
							return;
					}
	}

	// Could not be merged, save it for later
	segments.push_back(segment);
}

/**
 * Add new segment into the array, merging with other segments on the same line
 */
void addAndMergeSameLine(const OrientedLineSegment& segment, std::vector<struct OrientedLineSegment>& otherSegments, const Face& face)
{
	const auto faceNormal = face.getNormal();

	for (auto otherSegmentIt = otherSegments.begin(); otherSegmentIt != otherSegments.end(); ++otherSegmentIt)
	{
		if (segment.segment.canMerge(otherSegmentIt->segment) && glm::dot(segment.insideDir, otherSegmentIt->insideDir) > 0)
		{
			// Segments are along the same line, combine them into one
			OrientedLineSegment newSegment(segment.segment.merge(otherSegmentIt->segment), segment.insideDir);
			newSegment.fixInsideDir(faceNormal);
			otherSegments.erase(otherSegmentIt);
			addAndMergeSameLine(newSegment, otherSegments, face);
			return;
		}
	}

	otherSegments.push_back(segment);
}

void removeShortSegments(vector<OrientedLineSegment>& segments)
{
	for (size_t i = 0; i < segments.size();)
	{
		if (segments[i].segment.getLength() < 0.001f)
		{
			segments.erase(segments.begin() + i);
		}
		else
		{
			++i;
		}
	}
}

/**
 * Find all vertices between newFaceVertices
 */
std::vector<vec3> findEdgesBetweenNewVerts(const size_t startVertIdx, const std::vector<vec3>& faceVerts, const std::unordered_set<vec3>& newFaceVerts, std::unordered_set<vec3>& usedVerts, bool goForward)
{
	assert(startVertIdx < faceVerts.size());
	auto nextIt = [&](const auto& it) {return it + 1 != faceVerts.end() ? it + 1 : faceVerts.begin(); };
	auto prevIt = [&](const auto& it) {return it == faceVerts.begin() ? faceVerts.end() - 1 : it - 1; };

	std::vector<vec3> result;
	auto currentIt = faceVerts.begin() + startVertIdx;
	usedVerts.insert(*currentIt);
	result.emplace_back(*currentIt);
	while (true)
	{
		const auto secondVertIt = goForward ? nextIt(currentIt) : prevIt(currentIt);
		result.emplace_back(*secondVertIt);
		currentIt = secondVertIt;

		if (newFaceVerts.find(*currentIt) != newFaceVerts.end())
		{
			usedVerts.insert(*currentIt);
			break; //found another intersection, lets stop here
		}

	}

	return result;
}

/**
 * Remove segments that do not start and end at the face edges.
 * Those segments can come from bad input data, floating precission errors, or infinitely thin slices along the face boundary
 */
void removeBrokenSegments(vector<OrientedLineSegment>& segments, const Face& face)
{
	Plane facePlane(face);
	const vec3 faceNormal = face.getNormal();

	const vec3 faceCenter = std::accumulate(face.vertices.begin(), face.vertices.end(), vec3(0)) / static_cast<float>(face.vertices.size());

	for (size_t i = 0; i < segments.size();)
	{
		const auto& segment = segments[i];

		const bool isSegmentOnPlane = std::max(facePlane.pointDistance(segment.segment.getStart()),
			facePlane.pointDistance(segment.segment.getEnd())) < 0.001f;
		// Test segment for each face edge
		float closestStartDist = std::numeric_limits<float>::max();
		for (size_t j = 0; j < face.vertices.size(); j++)
		{
			const size_t nextIdx = (j + 1) % face.vertices.size();
			const vec3 edgeDir = face.vertices[nextIdx] - face.vertices[j];
			Plane edgePlane(face.vertices[j], face.vertices[nextIdx], face.vertices[j] + faceNormal);
			closestStartDist = std::min(edgePlane.pointDistance(segment.segment.getStart()), closestStartDist);
		}

		float closestEndDist = std::numeric_limits<float>::max();
		for (size_t j = 0; j < face.vertices.size(); j++)
		{
			const size_t nextIdx = (j + 1) % face.vertices.size();
			const vec3 edgeDir = face.vertices[nextIdx] - face.vertices[j];
			Plane edgePlane(face.vertices[j], face.vertices[nextIdx], face.vertices[j] + faceNormal);
			closestEndDist = std::min(edgePlane.pointDistance(segment.segment.getEnd()), closestEndDist);
		}

		const bool bothVerticesOnEdge = closestStartDist <= EPSILON_SCALE * glm::epsilon<float>() &&
			closestEndDist <= EPSILON_SCALE * glm::epsilon<float>();

		if (!bothVerticesOnEdge || !isSegmentOnPlane)
		{
			segments.erase(segments.begin() + i);
		}
		else
		{
			++i;
		}
	}
}

std::vector<Triangle> createCap(const std::vector<struct OrientedLineSegment>& splitSegments, const Face& face, std::unordered_set<glm::vec3>& globalUsedVerts)
{
	const vec3 faceNormal = face.getNormal();
	std::vector<Triangle> triangles;
	std::vector<OrientedLineSegment> cutSegments;
	std::vector<OrientedLineSegment> filteredSegments;

	// Filter away segments outside the face and cut crossing segments
	for (const auto s : splitSegments)
	{
		auto f = OrientedLineSegment(cutSegmentByFaceEdges(s.segment, face), s.insideDir);
		if (f.segment.getLength() > 0)
		{
			addAndMergeSameLine(f, cutSegments, face);
		}
	}

	std::vector<std::string> gbrSegments = geogebraExport(cutSegments);

	// Attempt to merge segments into triangles
	for (const auto& s : cutSegments)
	{
		tryMergeSegments(s, filteredSegments, triangles, face);
	}

	gbrSegments = geogebraExport(filteredSegments);

	removeShortSegments(filteredSegments);

	//removeBrokenSegments(filteredSegments, face);

	// Segments that did not get merged or triangulized must intersect the face edge at both points
	std::unordered_set<vec3> newFaceVerts;
	Face allEdgeVerts(face);
	for (const auto& s : filteredSegments)
	{
		allEdgeVerts.vertices.emplace_back(s.segment.getStart());
		allEdgeVerts.vertices.emplace_back(s.segment.getEnd());
		newFaceVerts.insert(s.segment.getStart());
		newFaceVerts.insert(s.segment.getEnd());
	}

	allEdgeVerts.orient(face.getNormal()); //Oriented to CW order from the outside
	const vec3 centerPoint = allEdgeVerts.calcCenterPoint();


	std::unordered_set<vec3> bannedVerts; //vertices used to create this face's cap
	std::vector<OrientedLineSegment> segmentsToAdd;
	// Merge with some face segments to complete the triangles
	for (const auto& segment : filteredSegments)
	{
		const auto intersectionStartVert = std::find_if(allEdgeVerts.vertices.begin(), allEdgeVerts.vertices.end(), [&](const vec3& v)
		{
			return fVecEquals(segment.segment.getStart(), v);
		});
		const auto intersectionEndVert = std::find_if(allEdgeVerts.vertices.begin(), allEdgeVerts.vertices.end(), [&](const vec3& v)
		{
			return fVecEquals(segment.segment.getEnd(), v);
		});


		const auto lineDirecton = segment.segment.a < segment.segment.b ? segment.segment.line.direction : -segment.segment.line.direction;
		const vec3 clockwiseDirection = glm::cross(-faceNormal, lineDirecton); //CW direction from the outside

		const bool goClockwise = glm::dot(clockwiseDirection, segment.insideDir) > 0;


		if (bannedVerts.find(*intersectionStartVert) == bannedVerts.end())
		{
			const size_t startIdx = intersectionStartVert - allEdgeVerts.vertices.begin();
			std::vector<vec3> pointsToAdd = findEdgesBetweenNewVerts(startIdx, allEdgeVerts.vertices, newFaceVerts, bannedVerts, goClockwise);
			globalUsedVerts.insert(pointsToAdd[0]);
			for (size_t i = 1; i < pointsToAdd.size(); i++)
			{
				OrientedLineSegment newSegment(pointsToAdd[i - 1], pointsToAdd[i], centerPoint - pointsToAdd[i]);
				newSegment.fixInsideDir(faceNormal);
				segmentsToAdd.emplace_back(newSegment);
				globalUsedVerts.insert(pointsToAdd[i]);
			}
		}

		if (bannedVerts.find(*intersectionEndVert) == bannedVerts.end())
		{
			const size_t startIdx = intersectionEndVert - allEdgeVerts.vertices.begin();
			std::vector<vec3> pointsToAdd = findEdgesBetweenNewVerts(startIdx, allEdgeVerts.vertices, newFaceVerts, bannedVerts, !goClockwise);
			globalUsedVerts.insert(pointsToAdd[0]);
			for (size_t i = 1; i < pointsToAdd.size(); i++)
			{
				OrientedLineSegment newSegment(pointsToAdd[i - 1], pointsToAdd[i], centerPoint - pointsToAdd[i]);
				newSegment.fixInsideDir(faceNormal);
				segmentsToAdd.emplace_back(newSegment);
				globalUsedVerts.insert(pointsToAdd[i]);
			}
		}

	}

	for (const auto& s : segmentsToAdd)
	{
		tryMergeSegments(s, filteredSegments, triangles, face);
	}

	if (!filteredSegments.empty())
	{
		//assert(filteredSegments.empty());
		1 + 1;
	}


	return triangles;
}

std::vector<std::string> geogebraExport(const std::vector<OrientedLineSegment>& segments)
{
	std::vector<std::string> gbrSegments;
	auto v3ToStr = [](const glm::vec3& v) {return std::string("(") + std::to_string(v.x) + "," + std::to_string(v.y) + "," + std::to_string(v.z) + ")"; };
	for (const auto& s : segments)
	{
		gbrSegments.emplace_back(std::string("Segment(") + v3ToStr(s.segment.getStart()) + "," + v3ToStr(s.segment.getEnd()) + ")");
	}

	return gbrSegments;
}

bool intersectMesh(const TriMesh& source, const std::vector<Face>& faces, TriMesh& outMesh)
{
	/*
	 * Intersect all triangles in source mesh with the cell
	 * Triangles that intersect are split, others that are fully inside the cell are kept whole
	 */

	auto triangles = getTriangles(source);
	std::vector<Triangle> capTriangles;
	std::unordered_set<glm::vec3> usedVertices; //face vertices that are part of the mesh

	std::vector<Face> emptyFaces;
	for (int i = 0; i < faces.size(); i++)
	{
		const auto& face = faces[i];

		std::vector<Triangle> newTriangles;
		std::vector<OrientedLineSegment> splitSegments; //segments created by splitting a triangle in half
		const Plane plane(face);

		for (const Triangle& tri : triangles)
		{
			auto splitResult = cutTriangleByPlane(tri, plane, splitSegments);
			newTriangles.insert(std::end(newTriangles), std::begin(splitResult), std::end(splitResult));
		}


		std::vector<std::string> gbrSegments = geogebraExport(splitSegments);


		auto capResult = createCap(splitSegments, face, usedVertices);

		capTriangles.insert(std::end(capTriangles), std::begin(capResult), std::end(capResult));

		if (triangles.empty())
		{
			emptyFaces.push_back(face);
		}

		std::swap(triangles, newTriangles);
	}

	// If face has no intersections and its vertices are used by other faces fill the whole face
	for (auto& face : emptyFaces)
	{
		if (usedVertices.find(face.vertices[0]) != usedVertices.end())
		{
			//assume the rest are used aswell
			face.orient(-face.getNormal()); //Orient in CCW order from the outside
			auto faceTris = face.triangulate();
			capTriangles.insert(capTriangles.end(), faceTris.begin(), faceTris.end());
		}
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
	if (halfspace.isInFrontOrInside(triangle.a)) verticesToKeep.push_back(triangle.a);
	if (halfspace.isInFrontOrInside(triangle.b)) verticesToKeep.push_back(triangle.b);
	if (halfspace.isInFrontOrInside(triangle.c)) verticesToKeep.push_back(triangle.c);

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
