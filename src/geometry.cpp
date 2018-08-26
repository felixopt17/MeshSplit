#include "geometry.h"


LineSegment LineSegment::intersect(const LineSegment& other) const
{
	assert(fVecEquals(line.origin, other.line.origin));
	assert(fVecEquals(line.direction, other.line.direction));

	return LineSegment(line, std::max(a, other.a), std::min(b, other.b));
}

Plane::Plane(const struct Triangle& tri) :Plane(tri.a, tri.b, tri.c) {}
Plane::Plane(const struct Face& face)
{
	assert(face.vertices.size() >= 3);
	if (face.vertices.size() <= 3)
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
	if (originalSegment.getLength() == 0)
		return originalSegment;

	const Plane facePlane(face);
	assert(face.vertices.size() >= 3);
	assert(glm::epsilonEqual(glm::dot(originalSegment.line.direction, facePlane.normal), 0.f, EPSILON_SCALE*glm::epsilon<float>())); //are perpendicular
	assert(glm::epsilonEqual(glm::dot(originalSegment.line.origin, facePlane.normal), facePlane.scale, EPSILON_SCALE*glm::epsilon<float>())); //is point on plane

	/**
	 * 1) Project line to the face
	 * 2) Cut 2d line segment by each of the edges
	 */


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
		glm::vec3 outDirection = glm::normalize(glm::cross(edgeDirection, segment.line.direction));
		Plane edgePlane(edgeStart, glm::cross(outDirection, edgeDirection));

		const float t = (edgePlane.scale - glm::dot(edgePlane.normal, segment.getStart())) / 
			glm::dot(edgePlane.normal, segment.line.direction*segment.getLength());

		if(t>=0.f && t<=1.0f)
		{
			const auto splitPointValue = segment.a + t * segment.getLength();
			//TODO test that we keep the correct side of the segment
			if(isPointInFront(edgePlane, segment.getStart()))
			{
				// Segment starts in front of the plane
				// Keep start -> split 
				segment.b = splitPointValue;
			}
			else
			{
				// Segment starts behind the plane
				// Keep split-> end
				segment.a = splitPointValue;
			}
		}
	}

	return segment;
}


bool isPointInFront(const Plane& plane, const glm::vec3& point)
{
	return plane.isInFront(point);
}
