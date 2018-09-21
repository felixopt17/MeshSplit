#pragma once
#include "cinder/gl/gl.h"
#include "vorohelpers.h"
#include "glm/gtc/epsilon.hpp"
#include <algorithm>

constexpr float EPSILON_SCALE = static_cast<float>(2<<9); //1024

template<typename vecType>
bool allFalse(const vecType& bVec)
{
	return bVec.x == false && bVec.y == false && bVec.z == false;
}

template<typename vecType>
bool allTrue(const vecType& bVec)
{
	return bVec.x && bVec.y && bVec.z;
}

/// Comnpare two float vectors for equality with eps. tolerance
template<typename fVecType>
bool fVecEquals(const fVecType& a, const fVecType& b)
{
	const auto epsilonVec3 = glm::vec3(glm::epsilon<float>());

	return allTrue<fVecType>(glm::epsilonEqual(a, b, epsilonVec3*static_cast<float>(EPSILON_SCALE)));
}



struct Line
{
	Line(const glm::vec3& origin, const glm::vec3& direction) :origin(origin), direction(glm::normalize(direction)) {}
	Line& operator=(const Line& other) = default;


	glm::vec3 origin;
	glm::vec3 direction;
};

struct LineSegment
{
	LineSegment(const Line& line, float a = -KINDA_BIG_NUMBER, float b = KINDA_BIG_NUMBER) :
	line(line), a(a), b(b) {}
	LineSegment(const glm::vec3& from, const glm::vec3& to):line(from, (to-from)),
	a(0), b(glm::length(to-from)){}

	glm::vec3 getStart() const { return line.origin + a * line.direction; }
	glm::vec3 getEnd() const { return line.origin + b * line.direction; }

	/// Can merge with this line segment to form a continuous segment?
	bool canMerge(const LineSegment& other) const;

	/// Merge with a line segment to create a new one. Both segments must be connected
	LineSegment merge(const LineSegment& other) const;

	LineSegment& operator=(const LineSegment& other) = default;

	float getLength() const { return std::max(b - a, 0.f); }

	LineSegment intersect(const LineSegment& other) const;

	/// Cut line segment by a plane, keeping the part of the segment in front of the plane
	LineSegment cut(const struct Plane& plane) const;

	Line line;

	float a;
	float b;

private:
	/// Starting value for segment limits. Starting at float limits would be too imprecise
	static const int KINDA_BIG_NUMBER = 10000; 
};

struct OrientedLineSegment
{
	OrientedLineSegment(const glm::vec3& from, const glm::vec3& to, const glm::vec3& insideDir) :segment(from, to), insideDir(insideDir){}
	OrientedLineSegment(const LineSegment& segment, const glm::vec3& insideDir):segment(segment), insideDir(insideDir){}

	/// Fix inside direction so that it is perpendicular to plane normal
	void fixInsideDir(const glm::vec3& planeNormal)
	{
		auto fixedDir = glm::cross(segment.line.direction, planeNormal);
		if (glm::dot(insideDir, fixedDir) > 0)
			insideDir = fixedDir;
		else
			insideDir = -fixedDir;
	}
	LineSegment segment;
	glm::vec3 insideDir;
};

/**
 * Plane in 3D space
 */
struct Plane
{
	Plane(const glm::vec3& point, const glm::vec3& normal) : normal(glm::normalize(normal)), scale(glm::dot(point, normal)) {}
	Plane(const glm::vec3& pointA, const glm::vec3& pointB, const glm::vec3& pointC)
		: Plane(pointA, glm::normalize(glm::cross((pointB - pointA), (pointC - pointA)))) {}

	Plane(const struct Triangle& tri);
	Plane(const struct Face& face);

	/// Calculate intersection with other plane
	struct PlaneIntersectionResult intersect(const Plane& other) const;

	/// Is point in front of the plane (does not use thickness)
	bool isInFrontStrict(const glm::vec3& point) const
	{
		return glm::dot(normal, point) >= scale;
	}

	bool isInside(const glm::vec3& point) const
	{
		return pointDistance(point) < THICKNESS;
	}

	float pointDistance(const glm::vec3& point) const
	{
		return glm::abs(glm::dot(normal, point) - scale);
	}

	glm::vec3 normal;
	float scale;
	static const float THICKNESS;
};

struct PlaneIntersectionResult
{
	PlaneIntersectionResult(const Line& line) :line(line), isLine(true), isPlane(false) {}
	PlaneIntersectionResult(const Plane& plane) :plane(plane), isLine(false), isPlane(true) {}
	PlaneIntersectionResult() :isLine(false), isPlane(false) {}

	union
	{
		Line line;
		Plane plane;
	};

	const bool isLine;
	const bool isPlane;
};


struct Triangle
{
	Triangle(const glm::vec3& a, const glm::vec3& b, const  glm::vec3& c) :a(a), b(b), c(c) {}

	glm::vec3 a;
	glm::vec3 b;
	glm::vec3 c;

	Plane getPlane()const { return Plane(a, b, c); }
	glm::vec3 getNormal() const { return glm::normalize(glm::cross(b - a, c - a)); }

	/// Change winding order of vertices
	void changeWinding() { std::swap(a, c); }
};

/// Cut line segment going in the plane of the face with face's edges
LineSegment cutSegmentByFaceEdges(const LineSegment& segment, const Face& face);

/// Is point in front of the plane (or inside)
bool isPointInFront(const Plane& plane, const glm::vec3& point);


/// Split triangle by a plane. Return part of the triangle in front of the plane. Output splitting segments into a segment array
std::vector<Triangle> cutTriangleByPlane(const Triangle& triangle, const Plane& plane, std::vector<struct OrientedLineSegment>& segments);