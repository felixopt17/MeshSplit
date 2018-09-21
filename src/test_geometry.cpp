#include <cinder/app/App.h>
#include <cinder/app/RendererGl.h>
#ifdef _TEST_
#include "geometry.h"
#include <gtest/gtest.h>

// -----------
// Line
// -----------

using glm::vec3;

void assertVec(const glm::vec3& a, const glm::vec3& b)
{
	ASSERT_FLOAT_EQ(a.x, b.x);
	ASSERT_FLOAT_EQ(a.y, b.y);
	ASSERT_FLOAT_EQ(a.z, b.z);
}

bool lineEquals(const Line& a, const Line& b)
{
	const auto originDifference = b.origin - a.origin;
	auto directionRatio = originDifference / a.direction;

	bool equals = true;
	// Make sure direction has at least one valid component
	equals &= !(glm::epsilonEqual(glm::length2(a.direction), 0.f, glm::epsilon<float>()*EPSILON_SCALE));

	// Vector between origins must be a multiple of a direction vector
	// If direction vector has 0 avoid nan comparison



	if(glm::epsilonNotEqual(a.direction.x, 0.f, glm::epsilon<float>()*EPSILON_SCALE)
		&& glm::epsilonNotEqual(a.direction.y, 0.f, glm::epsilon<float>()*EPSILON_SCALE))
	{
		equals &= glm::epsilonEqual(directionRatio.x, directionRatio.y, glm::epsilon<float>()*EPSILON_SCALE);
	}

	if (glm::epsilonNotEqual(a.direction.y, 0.f, glm::epsilon<float>()*EPSILON_SCALE)
		&& glm::epsilonNotEqual(a.direction.z, 0.f, glm::epsilon<float>()*EPSILON_SCALE))
	{
		equals &= glm::epsilonEqual(directionRatio.y, directionRatio.z, glm::epsilon<float>()*EPSILON_SCALE);
	}

	if (glm::epsilonNotEqual(a.direction.x, 0.f, glm::epsilon<float>()*EPSILON_SCALE)
		&& glm::epsilonNotEqual(a.direction.z, 0.f, glm::epsilon<float>()*EPSILON_SCALE))
	{
		equals &= glm::epsilonEqual(directionRatio.x, directionRatio.z, glm::epsilon<float>()*EPSILON_SCALE);
	}

	return equals;
}

TEST(Asserts, LineEqualsTest)
{
	vec3 a(5.f, 2.f, 3.f);
	vec3 b(6.f, 3.f, 3.f);
	vec3 d(1.f, 1.f, 0.f);

	const Line line1(a, d);
	const Line line2(b, d);
	const Line line3(b, -d); //reversed direction vector
	ASSERT_TRUE(lineEquals(line1, line2)); //same direction, different origin
	ASSERT_TRUE(lineEquals(line1, line3)); //reverse direction, different origin

	vec3 c(2.f, 0.f, 1.f);
	const Line otherLine(c, d);
	ASSERT_FALSE(lineEquals(line1, otherLine));
}

TEST(LineTest, DirectionNormalized)
{
	const Line line(glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f, 0.f, 2.f));
	assertVec(line.direction, glm::vec3(0.f, 0.f, 1.f));
}

TEST(LineSegmentTest, LineStartEnd)
{
	const float LINE_SIZE = 3.f;
	Line line(glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f, 0.f, 1.f));
	LineSegment segment(line, 0, LINE_SIZE);

	auto startPoint = segment.getStart();
	assertVec(startPoint, glm::vec3(0.f, 0.f, 0.f));
	assertVec(segment.getEnd(), glm::vec3(0.f, 0.f, LINE_SIZE));
	ASSERT_FLOAT_EQ(segment.getLength(), LINE_SIZE);
}

TEST(LineSegmentTest, IntersectReducesInterval)
{

	const Line line(glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f, 0.f, 1.f));
	const LineSegment segment(line, -10.f, 8.f);
	const LineSegment segment2(line, -2.f, 6.f);

	LineSegment intersection = segment.intersect(segment2);
	ASSERT_FLOAT_EQ(intersection.getLength(), 8.f);
	ASSERT_FLOAT_EQ(intersection.a, -2.f);
	ASSERT_FLOAT_EQ(intersection.b, 6.f);
}

TEST(LineSegmentTest, IntersectPreservesInterval)
{
	const Line line(glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f, 0.f, 1.f));
	const LineSegment segment(line, -2.f, 6.f);
	const LineSegment segment2(line, -10.f, 8.f);

	LineSegment intersection = segment.intersect(segment2);
	ASSERT_FLOAT_EQ(intersection.getLength(), 8.f);
	ASSERT_FLOAT_EQ(intersection.a, -2.f);
	ASSERT_FLOAT_EQ(intersection.b, 6.f);
}

TEST(LineSegmentTest, IntersectPlaneStraight)
{
	const Plane horizontalPlane(glm::vec3(0.f, 5.f, 0.f), glm::vec3(0.f, 1.f, 0.f));
	const Line line(vec3(0.f, 10.f, 0.f), vec3(0.f, 1.f, 0.f));
	const Line reverseLine(vec3(0.f, 10.f, 0.f), vec3(0.f, -1.f, 0.f));
	const LineSegment segment(line, -10.f, 10.f);
	const LineSegment reverseSegment(reverseLine, -10.f, 10.f);

	// Test segment starting behind the plane
	auto cutSegment = segment.cut(horizontalPlane);
	assertVec(cutSegment.getStart(), vec3(0.f, 5.0f, 0.f));
	assertVec(cutSegment.getEnd(), vec3(0.f, 20.0f, 0.f));

	// Test segment starting in front of the plane
	auto cutSegmentReverse = reverseSegment.cut(horizontalPlane);
	assertVec(cutSegmentReverse.getStart(), vec3(0.f, 20.0f, 0.f));
	assertVec(cutSegmentReverse.getEnd(), vec3(0.f, 5.0f, 0.f));
}

TEST(LineSegmentTest, IntersectPlane)
{
	const Plane horizontalPlane(glm::vec3(0.f, 5.f, 0.f), glm::vec3(0.f, 1.f, 0.f));
	const Line line(vec3(0.f, 0.f, 0.f), vec3(1.f, 1.f, 1.f));
	const LineSegment segment(line, -10.f, 10.f);

	auto cutSegment = segment.cut(horizontalPlane);
	assertVec(cutSegment.getStart(), vec3(5.f, 5.0f, 5.f));
	assertVec(cutSegment.getEnd(),segment.getEnd());
}

TEST(LineSegmentTest, IntersectPlaneRegression1)
{
	const Plane plane(vec3(-3.f, 0.f, 0.f), vec3(1.f, 0.f, 0.f));
	const LineSegment segment(vec3(0.12f, 0.8f, 0.f), vec3(0.0f, 0.8f, 0.f));

	auto cutSegment = segment.cut(plane);
	assertVec(cutSegment.getStart(), segment.getStart());
	assertVec(cutSegment.getEnd(), segment.getEnd());
}

TEST(LineSegmentTest, MergeSameLine)
{
	const Line line1(vec3(0.f, 0.f, 0.f), vec3(1.f, 1.f, 1.f));

	const LineSegment a(line1, 0.f, 1.f);
	const LineSegment b(line1, 0.5f, 2.f);
	ASSERT_TRUE(a.canMerge(b));
	ASSERT_TRUE(b.canMerge(a));

	const auto result = a.merge(b);
	EXPECT_FLOAT_EQ(result.a , 0.0f);
	EXPECT_FLOAT_EQ(result.b,  2.0f);

	assertVec(result.line.direction, line1.direction);
	assertVec(result.line.origin, line1.origin);

	// Make sure result is correct from either side
	const auto result2 = b.merge(a);
	assertVec(result.getStart(), result2.getStart());
	assertVec(result.getEnd(), result2.getEnd());
}

TEST(LineSegmentTest, MergeFailTest)
{
	const Line line1(vec3(0.f, 0.f, 0.f), vec3(1.f, 1.f, 1.f));
	const Line line2(vec3(0.f, 0.f, 0.f), vec3(1.f, 0.f, -1.f));

	const LineSegment a(line1, 0.f, 1.f);
	const LineSegment b(line2, 0.5f, 2.f);
	EXPECT_FALSE(a.canMerge(b));
	EXPECT_FALSE(b.canMerge(a));
}

TEST(LineSegmentTest, MergeOppositeLine)
{
	const Line line1(vec3(0.f, 0.f, 0.f), vec3(1.f, 1.f, 1.f));
	const Line line2(vec3(0.f, 0.f, 0.f), vec3(1.f, 0.f, -1.f));

	const LineSegment a(line1, 0.f, 1.f);

	/// Origin on the correct line, reverse direction vector
	const LineSegment d(Line(line1.origin, -line1.direction), 0.0f, 2.f);
	ASSERT_TRUE(a.canMerge(d));
	ASSERT_TRUE(d.canMerge(a));

	const auto result = a.merge(d);
	EXPECT_FLOAT_EQ(result.a, -2.f);
	EXPECT_FLOAT_EQ(result.b, 1.f);
}

TEST(PlaneTest, PlaneConstructorEquality)
{
	const glm::vec3 a(0.f, 0.f, 0.f);
	const glm::vec3 b(0.f, 0.f, 1.f);
	const glm::vec3 c(1.f, 0.f, 1.f);
	const glm::vec3 n(0.f, 1.f, 0.f);

	// Make sure a plane constructed in different way is the same
	Plane plane1(a, b, c);
	Plane plane2(a, n);
	assertVec(plane1.normal, plane2.normal);
	ASSERT_FLOAT_EQ(plane1.scale, plane2.scale);

	Plane plane3(Triangle(a, b, c));
	assertVec(plane1.normal, plane3.normal);
	ASSERT_FLOAT_EQ(plane1.scale, plane3.scale);

	Face f;
	f.vertices = { a,b,c };

	Plane plane4(f);
	assertVec(plane1.normal, plane4.normal);
	ASSERT_FLOAT_EQ(plane1.scale, plane4.scale);
}

TEST(PlaneTest, PlaneIntersectionLineBasic)
{
	const Plane horizontalPlane(glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f, 1.f, 0.f));
	const Plane verticalPlane(glm::vec3(0.f, 0.f, 0.f), glm::vec3(1.f, 0.f, 0.f));

	PlaneIntersectionResult intersection = horizontalPlane.intersect(verticalPlane);

	ASSERT_EQ(intersection.isLine, true);
	ASSERT_EQ(intersection.isPlane,false);

	const Line correctIntersection(vec3(0.f, 0.f, 0.f), vec3(0.f, 0.f, 1.f));
	ASSERT_TRUE(lineEquals(correctIntersection, intersection.line));
}

TEST(PlaneTest, PlaneIntersectionLineBasic2)
{
	const vec3 planeOrigin(3.f, 5.f, -4.f);
	const Plane horizontalPlane(planeOrigin, glm::vec3(0.f, 1.f, 0.f));
	const Plane verticalPlane(glm::vec3(0.f, -2.f, 3.f), glm::vec3(1.f, 0.f, 0.f));

	PlaneIntersectionResult intersection = horizontalPlane.intersect(verticalPlane);

	ASSERT_EQ(intersection.isLine, true);
	ASSERT_EQ(intersection.isPlane, false);

	const Line correctIntersection(planeOrigin, vec3(0.f, 0.f, 1.f));
	ASSERT_TRUE(lineEquals(correctIntersection, intersection.line));
}


TEST(PlaneTest, PointInFrontBasic)
{
	const vec3 planeOrigin(3.f, 5.f, -4.f);
	const Plane horizontalPlane(planeOrigin, glm::vec3(0.f, 1.f, 0.f));

	const glm::vec3 a(0.f, 0.f, 0.f);
	ASSERT_FALSE(horizontalPlane.isInFrontStrict(a)); //Point below

	const glm::vec3 b(1.f, 5.5f, -4.f);
	ASSERT_TRUE(horizontalPlane.isInFrontStrict(b)); //Point above

	ASSERT_TRUE(horizontalPlane.isInFrontStrict(planeOrigin)); //Point in plane
}

TEST(TriangleTest, TrianglePlaneCutNoChange)
{
	const Triangle tri(vec3(0.f, 0.f, 0.f), vec3(1.f, 0.f, 1.f), vec3(1.f, 0.f, 0.f));
	const Plane plane(vec3(2.f, 0.f, 0.f), vec3(-1.f, 0.f, 0.f));
	std::vector<OrientedLineSegment> segments;
	auto cutResult = cutTriangleByPlane(tri, plane, segments);
	// Triangle should remain unchanged
	ASSERT_TRUE(cutResult.size() == 1);
	const auto& resultTri = cutResult[0];
	assertVec(tri.a, resultTri.a);
	assertVec(tri.b, resultTri.b);
	assertVec(tri.c, resultTri.c);
}

TEST(TriangleTest, TrianglePlaneCutRemove)
{
	const Triangle tri(vec3(0.f, 0.f, 0.f), vec3(1.f, 0.f, 1.f), vec3(1.f, 0.f, 0.f));
	const Plane plane(vec3(-2.f, 0.f, 0.f), vec3(-1.f, 0.f, 0.f));
	std::vector<OrientedLineSegment> segments;
	auto cutResult = cutTriangleByPlane(tri, plane, segments);
	// Triangle should be completely removed
	ASSERT_TRUE(cutResult.empty());
}

TEST(TriangleTest, TrianglePlaneCut1Vec)
{
	const Triangle tri(vec3(0.f, 0.f, 0.f), vec3(1.f, 0.f, 1.f), vec3(1.f, 0.f, 0.f));
	const Plane plane(vec3(0.5f, 0.f, 0.f), vec3(-1.f, 0.f, 0.f));
	std::vector<OrientedLineSegment> segments;
	auto cutResult = cutTriangleByPlane(tri, plane, segments);
	// Triangle should be partially clipped
	ASSERT_TRUE(cutResult.size() == 1);
	const auto& resultTri = cutResult[0];
	assertVec(resultTri.a, vec3(0.f, 0.f, 0.f));
	assertVec(resultTri.b, vec3(0.5f, 0.f, 0.5f));
	assertVec(resultTri.c, vec3(0.5f, 0.f, 0.f));
}

TEST(TriangleTest, TrianglePlaneCut2Vec)
{
	const Triangle tri(vec3(0.f, 0.f, 0.f), vec3(1.f, 0.f, 1.f), vec3(1.f, 0.f, 0.f));
	const Plane plane(vec3(0.5f, 0.f, 0.f), vec3(1.f, 0.f, 0.f));
	std::vector<OrientedLineSegment> segments;

	auto cutResult = cutTriangleByPlane(tri, plane, segments);
	// Triangle should be partially clipped, resulting in two triangles
	ASSERT_TRUE(cutResult.size() == 2);
	assertVec(cutResult[0].a, vec3(1.f, 0.f, 1.f));
	assertVec(cutResult[0].b, vec3(1.f, 0.f, 0.f));
	assertVec(cutResult[0].c, vec3(0.5f, 0.f, 0.5f));

	assertVec(cutResult[1].a, vec3(0.5f, 0.f, 0.f));
	assertVec(cutResult[1].b, vec3(0.5f, 0.f, 0.5f));
	assertVec(cutResult[1].c, vec3(1.f, 0.f, 0.f));
}




#include <windows.h>
#include <WinBase.h>
#include <shellapi.h>
#include <vector>

int __stdcall WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
	int argc = 0;
	LPWSTR* str = CommandLineToArgvW(GetCommandLineW(), &argc);
	std::vector < std::unique_ptr<char[]>> argBuffer;
	for(int i=0;i<argc;i++)
	{
		const size_t stringSize = wcslen(str[i]);
		argBuffer.emplace_back(std::make_unique<char[]>(stringSize));
		size_t bytesWritten = 0;
		wcstombs(argBuffer.at(i).get(), str[i], stringSize);
	}

	std::vector<char*> argv;
	for (auto& a : argBuffer)
		argv.push_back(a.get());

	::testing::InitGoogleTest(&argc,argv.data());
	return RUN_ALL_TESTS();
}


#endif