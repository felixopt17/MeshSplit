#include "vorohelpers.h"
#include "geometry.h"
#include <numeric>
#include "glm/gtx/vector_angle.hpp"

Face::Face(const Triangle& triangle)
{
	vertices.push_back(triangle.a);
	vertices.push_back(triangle.b);
	vertices.push_back(triangle.c);
}

glm::vec3 Face::getNormal() const
{
	auto directions = getDirectionVectors();
	return glm::normalize(glm::cross(directions[0], directions[1]));
}

std::array<glm::vec3, 2> Face::getDirectionVectors() const
{
	assert(vertices.size() >= 3);

	std::array<glm::vec3, 2> directions;
	directions[0] = glm::normalize(vertices[1] - vertices[0]);

	// Need to make sure that the second direction is not on the same line as the first one
	for (size_t i = 2; i < vertices.size(); i++)
	{
		directions[1] = glm::normalize(vertices[i] - vertices[0]);
		const auto dot = glm::dot(directions[0], directions[1]);

		const auto tolerance = EPSILON_SCALE * glm::epsilon<float>();
		if (dot > -1.f + tolerance && dot < 1.f - tolerance)
			break;
	}

	return directions;
}

void Face::removeDuplicates()
{
	std::vector<glm::vec3> newVerts;
	for (const glm::vec3& oldVert : vertices)
	{
		bool match = false;
		for (const glm::vec3& newVert : newVerts)
		{
			if (fVecEquals(oldVert, newVert))
			{
				match = true;
				break;
			}
		}

		if (!match)
			newVerts.push_back(oldVert);
	}

	std::swap(newVerts, vertices);
}

void Face::orient(const glm::vec3& normal)
{
	if (vertices.empty())
		return;

	const glm::vec3 centerPoint = std::accumulate(vertices.begin(), vertices.end(), glm::vec3(0.f, 0.f, 0.f), std::plus<glm::vec3>())
		/ static_cast<float>(vertices.size());


	// Calculate angle between vertices and center point
	const glm::vec3& a = glm::normalize(vertices[0] - centerPoint);
	std::vector<std::pair<glm::vec3, float>> vertAngles;
	for (int i = 1; i < vertices.size(); i++)
	{
		const glm::vec3& b = glm::normalize(vertices[i] - centerPoint);
		float angle = glm::orientedAngle(a, b, normal);
		if (angle < 0)
			angle += 2 * glm::pi<float>(); //keep angle positive

		vertAngles.emplace_back(std::make_pair(vertices[i], angle));
	}

	// Sort vertices by angle to the first vert
	std::sort(vertAngles.begin(), vertAngles.end(), [](const auto& a, const auto& b) {return a.second < b.second; });

	// Get vertices from the sorted array
	std::vector<glm::vec3> orderedVertices;
	orderedVertices.push_back(vertices[0]);
	for (const auto& vertAnglePair : vertAngles)
	{
		orderedVertices.push_back(vertAnglePair.first);
	}

	std::swap(orderedVertices, vertices);
}

std::vector<Face> getFacesFromEdges(voro::voronoicell& cell, const glm::vec3& particlePos)
{
	std::vector<int> vertexOrders;
	cell.vertex_orders(vertexOrders);


	const auto vertToPoint = [&](int vertex) {
		return glm::vec3(
			static_cast<float>(cell.pts[3 * vertex]),
			static_cast<float>(cell.pts[3 * vertex + 1]),
			static_cast<float>(cell.pts[3 * vertex + 2]));
	};

	std::set<std::set<int>> faceSets;
	const auto numVertices = cell.p;
	for (int vertex = 0; vertex < numVertices; vertex++) //For each vertex in the cell
	{
		const auto numEdges = cell.nu[vertex];
		for (int vertEdge = 0; vertEdge < numEdges; vertEdge++)
		{
			const auto neighborVert = cell.ed[vertex][vertEdge]; // For their neighbour

			const auto neighborEdges = cell.nu[neighborVert];
			for (int neighborEdge = 0; neighborEdge < neighborEdges; neighborEdge++) //For all their edges
			{
				const auto distantNeighbor = cell.ed[neighborVert][neighborEdge];
				if (distantNeighbor == vertex)
					continue; //ignore the edge we came from



				std::set<int> visitedVertices{ vertex, neighborVert };
				const Plane facePlane(vertToPoint(vertex), vertToPoint(neighborVert), vertToPoint(distantNeighbor));
				// Find next vertices in the cell that belong to this plane

				int currentVert = distantNeighbor;
				while (visitedVertices.find(currentVert) == visitedVertices.end())
				{
					visitedVertices.insert(currentVert);

					// Find a new neighbor on the same plane to visit
					for (int i = 0; i < cell.nu[currentVert]; i++)
					{
						const auto neighbor = cell.ed[currentVert][i];
						if (visitedVertices.find(neighbor) != visitedVertices.end())
							continue; //test new vertices only

						// Test if neighbor is close enough to the common plane
						if (facePlane.pointDistance(vertToPoint(neighbor)) < glm::epsilon<float>()*EPSILON_SCALE)
						{
							currentVert = neighbor;
							break;
						}
					}

				}

				faceSets.insert(visitedVertices);
			}
		}
	}




	glm::vec3 centerPoint(0.f, 0.f, 0.f);
	for (int i = 0; i < cell.p; i++)
		centerPoint += vertToPoint(i);
	centerPoint /= static_cast<float>(cell.p);

	// Create faces from the faceset
	std::vector<Face> result;
	for (const auto& faceSet : faceSets)
	{
		Face face;
		for (const int vert : faceSet)
		{
			face.vertices.push_back(vertToPoint(vert)+particlePos);
		}

		auto normal = face.getNormal();
		if(glm::dot(normal, centerPoint - (face.vertices[0]-particlePos)) < 0)
		{
			normal = -normal; //make normal face inward
		}

		face.orient(normal);
		result.push_back(face);
	}

	return result;
}

ci::TriMesh meshFromFaces(const std::vector<Face>& faces)
{
	cinder::TriMesh mesh;

	for(const Face& face :faces)
	{
		if (face.vertices.size() < 3)
			continue;


		const auto originIdx = mesh.getNumVertices();
		mesh.appendPosition(face.vertices[0]);

		for (uint32_t i = 2; i < static_cast<uint32_t>(face.vertices.size()); i++)
		{
			const uint32_t cornerIdx = static_cast<uint32_t>(mesh.getNumVertices());
			mesh.appendPosition(face.vertices[i - 1]);
			mesh.appendPosition(face.vertices[i]);
			mesh.appendTriangle(originIdx, cornerIdx + 1, cornerIdx + 0);
		}
	}

	/*for(size_t i=0;i<mesh.getNumVertices();i++)
		mesh.appendColorRgba(ci::ColorA(0.8f, 0.f, 0.f, 0.2f));*/

	mesh.recalculateNormals(true);
	return mesh;
}
