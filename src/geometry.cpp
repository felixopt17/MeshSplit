#include "geometry.h"

using glm::vec3;


LineSegment LineSegment::intersect(const LineSegment& other) const
{
	assert(fVecEquals(line.origin, other.line.origin));
	assert(fVecEquals(line.direction, other.line.direction));

	return LineSegment(line, std::max(a, other.a), std::min(b, other.b));
}

LineSegment LineSegment::cut(const Plane& plane) const
{
	// Line Segment - Plane intersection based on Real-Time Collision Detection by Christer Ericson
	const auto ab = getEnd() - getStart();
	auto t = (plane.scale - glm::dot(plane.normal, getStart())) / glm::dot(plane.normal, ab);

	if(t>=0.f && t<=1.0f)
	{
		LineSegment result(*this);

		// Keep the part in front
		if (plane.isInFrontStrict(getStart()))
		{
			const auto distFromA = t * result.getLength();
			result.b = result.a + distFromA;
		}
		else
		{
			const auto distFromB = (1 - t)*result.getLength();
			result.a = result.b - distFromB;
		}

		return result;
	}

	if (plane.isInside(line.origin))
		return LineSegment(*this);
	else
		return LineSegment(line, 0,0);

}

const float Plane::THICKNESS(std::numeric_limits<float>::epsilon()*EPSILON_SCALE);

Plane::Plane(const struct Triangle& tri) :Plane(tri.a, tri.b, tri.c) {}
Plane::Plane(const struct Face& face)
{
	assert(face.vertices.size() >= 3);
	if (face.vertices.size() < 3)
	{
		std::cerr << "Too few vertices in a face to make a plane! (" << face.vertices.size() << ")""\n";
		exit(1);
	}

	normal = face.getNormal();
	scale = glm::dot(face.vertices[0], normal);
}

PlaneIntersectionResult Plane::intersect(const Plane& other) const
{
	if (fVecEquals(normal, other.normal))
	{
		// The two planes parallel

		if (glm::epsilonEqual(scale, other.scale, EPSILON_SCALE*glm::epsilon<float>()))
			return PlaneIntersectionResult(*this);
		else
			return PlaneIntersectionResult();
	}

	// Intersection is a line

	// Line calculation based on method from Real-Time Collision Detection by Christer Ericson
	const auto lineDirection = glm::cross(normal, other.normal);
	const auto point = glm::cross(scale*other.normal - other.scale*normal, lineDirection);

	return PlaneIntersectionResult(Line(point, lineDirection));
}

LineSegment cutSegmentByFaceEdges(const LineSegment& originalSegment, const Face& face)
{
	if (glm::epsilonEqual(originalSegment.getLength(), 0.f, glm::epsilon<float>()*EPSILON_SCALE))
		return originalSegment;

	const Plane facePlane(face);


	LineSegment segment(originalSegment);
	const auto vertCount = face.vertices.size();
	for(size_t i=0; i<vertCount;i++)
	{
		const auto edgeStart = face.vertices[i];
		const auto edgeEnd = face.vertices[(i + 1) % vertCount];
		const auto edgeDirection = glm::normalize(edgeEnd - edgeStart);

		if (glm::epsilonEqual(abs(glm::dot(segment.line.direction, edgeDirection)), 1.f, EPSILON_SCALE*glm::epsilon<float>()))
			continue; // Line does not intersect the edge

		// Construct a temporary plane to test against
		const glm::vec3 outDirection = glm::normalize(glm::cross(edgeDirection, segment.line.direction));
		const Plane edgePlane(edgeStart, glm::cross(outDirection, edgeDirection));

		segment = segment.cut(edgePlane);
	}

	return segment;
}


bool isPointInFront(const Plane& plane, const glm::vec3& point)
{
	return plane.isInFrontStrict(point);
}


std::vector<Triangle> cutTriangleByPlane(const Triangle& triangle, const Plane& plane, std::vector<OrientedLineSegment>& segments)
{
	std::vector<vec3> pointsInFront;
	std::vector<vec3> pointsInBack;
	if (plane.isInFrontStrict(triangle.a)) pointsInFront.push_back(triangle.a); else pointsInBack.push_back(triangle.a);
	if (plane.isInFrontStrict(triangle.b)) pointsInFront.push_back(triangle.b); else pointsInBack.push_back(triangle.b);
	if (plane.isInFrontStrict(triangle.c)) pointsInFront.push_back(triangle.c); else pointsInBack.push_back(triangle.c);

	if (pointsInFront.empty())
		return {};

	if (pointsInFront.size() == 3)
		return { triangle };

	if (pointsInFront.size() == 1)
	{
		// Create one triangle

		// Cut sides of the triangle with the plane
		auto aToB = LineSegment(pointsInFront[0], pointsInBack[0]).cut(plane);
		auto aToC = LineSegment(pointsInFront[0], pointsInBack[1]).cut(plane);

		Triangle result(pointsInFront[0], aToB.getEnd(), aToC.getEnd());

		// Keep triangle pointing to the same direction
		if (glm::dot(result.getNormal(), triangle.getNormal()) < 0)
			result.changeWinding();

		segments.emplace_back(aToB.getEnd(), aToC.getEnd(), triangle.getNormal());
		return { result };
	}

	assert(pointsInFront.size() == 2);
	// Create two triangles

	auto aToC = LineSegment(pointsInFront[0], pointsInBack[0]).cut(plane);
	auto bToC = LineSegment(pointsInFront[1], pointsInBack[0]).cut(plane);

	vec3 intersections[2] = { aToC.getEnd(), bToC.getEnd() };

	Triangle result1(pointsInFront[0], pointsInFront[1], intersections[0]);
	Triangle result2(pointsInFront[1], intersections[0], intersections[1]);

	// Keep triangle pointing in the same direction
	if (glm::dot(result1.getNormal(), triangle.getNormal()) < 0)
		result1.changeWinding();
	if (glm::dot(result2.getNormal(), triangle.getNormal()) < 0)
		result2.changeWinding();

	segments.emplace_back(intersections[0], intersections[1], triangle.getNormal());
	return { result1, result2 };
}