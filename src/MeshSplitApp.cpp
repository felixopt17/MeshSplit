#include "cinder/app/App.h"
#include "cinder/app/RendererGl.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/Shader.h"


#include "voro++.hh"
#include "vorohelpers.h"
#include <random>
#include "meshsplitter.h"

#ifndef _TEST_

/// Scale mesh to fit into a uniform cube
void scaleMesh(TriMesh& original)
{
	glm::vec3 maxPos(0);


	glm::vec3* positions = original.getPositions<3>();

	const auto vertCount = original.getNumVertices();
	for (size_t i = 0; i < vertCount; i++)
	{
		maxPos.x = std::max(maxPos.x, positions[i].x);
		maxPos.y = std::max(maxPos.y, positions[i].y);
		maxPos.z = std::max(maxPos.z, positions[i].z);
	}

	const float scale = 1 / std::max(maxPos.x, std::max(maxPos.y, maxPos.z));

	for (size_t i = 0; i < vertCount; i++)
	{
		positions[i] *= scale;
	}
}

using namespace ci;
using namespace ci::app;
using namespace std;

std::default_random_engine re;

class MeshSplitApp : public App {
public:
	MeshSplitApp();
	void generateVoronoiCells();
	void setup() override;
	void mouseDown(MouseEvent event) override;
	void update() override;
	void draw() override;
	void keyDown(KeyEvent event) override;

	void drawFaceLines(Camera& cam);
private:
	TriMesh mesh;
	TriMesh faceMesh;
	voro::container con;

	bool drawFace = false;

	std::vector<std::vector<Face>> cells;
	std::vector<TriMesh> meshParts;
	std::vector<vec3> particlePositions;
	static const int MAX_SIZE = 2;

	int currentPart = 0;
	int currentFace = 0;

private:

};

MeshSplitApp::MeshSplitApp() :con(-MAX_SIZE, MAX_SIZE, -MAX_SIZE, MAX_SIZE, -MAX_SIZE, MAX_SIZE, 5, 5, 5, false, false, false, 8)
{

}


void MeshSplitApp::generateVoronoiCells()
{
	getWindow()->setTitle("Generating voronoi cells");

	const auto particleCount = 60;
	uniform_real_distribution<double> posDist(-MAX_SIZE, MAX_SIZE);
	for (int particleID = 0; particleID < particleCount; particleID++)
	{
		con.put(particleID++, posDist(re), posDist(re), posDist(re));
	}


	getWindow()->setTitle("Voronoi cells generated");

	// Loop over all particles and save computed voronoi cells
	voro::c_loop_all vLoop(con);
	vLoop.start();
	do
	{
		getWindow()->setTitle("Generating cell faces " + to_string(vLoop.pid()) + "/" + to_string(con.total_particles()));
		voro::voronoicell vcell;

		if (con.compute_cell(vcell, vLoop))
		{
			double px, py, pz;
			vLoop.pos(px, py, pz);
			const glm::vec3 particlePos(px, py, pz);

			cells.emplace_back(getFacesFromEdges(vcell, particlePos));
			particlePositions.emplace_back(particlePos);
		}

	} while (vLoop.inc());

	getWindow()->setTitle("Cell faces generated");
}

void MeshSplitApp::setup()
{
	cells.clear();
	meshParts.clear();
	particlePositions.clear();

	mesh = geom::Cube();
	scaleMesh(mesh);
	generateVoronoiCells();

	meshParts = splitMesh(mesh, cells);
}

void MeshSplitApp::mouseDown(MouseEvent event)
{
}



void MeshSplitApp::drawFaceLines(Camera& cam)
{
	gl::draw(faceMesh);
}

void MeshSplitApp::update()
{

}




void MeshSplitApp::draw()
{
	gl::enableDepth(true);
	gl::clear(Color(0.3f, 1.0f, 0.3f));
	gl::enableFaceCulling(true);

	gl::lineWidth(1);

	CameraPersp cam;
	cam.lookAt(vec3(4, 4, 4), vec3(0));
	gl::setMatrices(cam);

	const auto lambert = gl::ShaderDef().lambert();
	const auto shader = gl::getStockShader(lambert);
	shader->bind();

	float offsetScale = static_cast<float>(fmod(getElapsedSeconds(), 3)) - 1.0f;
	offsetScale = std::max(offsetScale, 0.f);

	gl::pushModelMatrix();
	//gl::scale(2, 2, 2);
	gl::rotate(static_cast<float>(3.1415*2.0*getElapsedSeconds() / 8), vec3(0, 1, 0));
	//gl::draw(mesh);

	for (size_t i = 0; i < meshParts.size(); i++)
	{
		gl::pushModelMatrix();

		gl::translate(particlePositions[i] * offsetScale);
		const auto& part = meshParts[i];
		if (i == currentPart % meshParts.size())
		{
			const bool oldState = gl::isWireframeEnabled();
			if (!oldState)
				gl::enableWireframe();
			gl::draw(part);
			if (!oldState)
				gl::disableWireframe();
		}
		else
		{
			gl::draw(part);
		}
		gl::popModelMatrix();
	}


	if (drawFace)
		drawFaceLines(cam);


	gl::popModelMatrix();

}

void MeshSplitApp::keyDown(KeyEvent event)
{
	if (event.getChar() == '+')
	{
		currentPart++;
		currentFace = 0;
		if (!meshParts.empty())
			mesh = meshParts[currentPart%meshParts.size()];

		auto& faces = cells[currentPart%cells.size()];
		faceMesh = meshFromFaces({ faces });
	}

	if (event.getChar() == 'p')
	{
		currentFace++;
		auto& faces = cells[currentPart%cells.size()];

		faceMesh = meshFromFaces({ faces[currentFace%faces.size()] });
	}

	if (event.getChar() == 'f')
	{
		drawFace = !drawFace;
		auto& faces = cells[currentPart%cells.size()];
		faceMesh = meshFromFaces({ faces[currentFace%faces.size()] });
	}

	if (event.getChar() == 'w')
	{
		static bool wireframe;
		if (wireframe)
			gl::enableWireframe();
		else
			gl::disableWireframe();

		wireframe = !wireframe;
	}

	if (event.getChar() == 'r')
	{
		re.seed(getElapsedFrames());
		setup();
	}
}



CINDER_APP(MeshSplitApp, RendererGl)

#endif