#include "cinder/app/App.h"
#include "cinder/app/RendererGl.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/Shader.h"


#include "voro++.hh"
#include "vorohelpers.h"
#include <random>
#include "meshsplitter.h"

#ifndef _TEST_


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

	void draw3DLine(Camera& cam, vec3 a, vec3 b);
	void drawFaceLines(Camera& cam);
private:
	TriMesh mesh;
	TriMesh faceMesh;
	voro::container con;

	bool drawFace = false;

	std::vector<std::vector<Face>> cells;
	std::vector<TriMesh> meshParts;
	static const int MAX_SIZE = 5;

	int currentPart = 0;
	int currentFace = 0;

private:

};

MeshSplitApp::MeshSplitApp() :con(-MAX_SIZE, MAX_SIZE, -MAX_SIZE, MAX_SIZE, -MAX_SIZE, MAX_SIZE, 1, 1, 1, false, false, false, 10)
{

}


void MeshSplitApp::generateVoronoiCells()
{
	getWindow()->setTitle("Generating voronoi cells");

	const auto scale = 2;
	uniform_real_distribution<double> offsetDist(-0.5, 0.5);
	int particleID = 0;
	for (int x = -scale; x <= scale; x++)
		for (int y = -scale; y <= scale; y++)
			for (int z = -scale; z <= scale; z++)
			{
				con.put(particleID++, (x + offsetDist(re)) / scale, (y + offsetDist(re)) / scale, (z + offsetDist(re)) / scale);
			}

	con.compute_all_cells();

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
			glm::vec3 particlePos(px, py, pz);

			cells.emplace_back(getFacesFromEdges(vcell, particlePos));
		}

	} while (vLoop.inc());

	getWindow()->setTitle("Cell faces generated");
}

void MeshSplitApp::setup()
{
	mesh = geom::Cube();
	generateVoronoiCells();

	meshParts = splitMesh(mesh, cells);
	//meshParts = testSplit(mesh);
}

void MeshSplitApp::mouseDown(MouseEvent event)
{
}


void MeshSplitApp::draw3DLine(Camera& cam, vec3 a, vec3 b)
{
	const auto screenA = cam.worldToScreen(a, getWindowWidth(), getWindowHeight());
	const auto screenB = cam.worldToScreen(b, getWindowWidth(), getWindowHeight());
	ci::gl::drawLine(screenA, screenB);
}

void MeshSplitApp::drawFaceLines(Camera& cam)
{
	/*auto& faces = cells[currentPart%meshParts.size()];

	auto& face = faces[currentFace%faces.size()];
	for(size_t i=1; i<=face.vertices.size();i++)
	{
		draw3DLine(cam, face.vertices[i - 1], face.vertices[i%face.vertices.size()]);
	}*/
	gl::draw(faceMesh);
}

void MeshSplitApp::update()
{

}




void MeshSplitApp::draw()
{
	gl::enableDepth(true);
	gl::clear(Color(0.5, 0.5, 0.5));
	gl::enableFaceCulling(false);

	gl::lineWidth(3);

	CameraPersp cam;
	cam.lookAt(vec3(2, 2, 2), vec3(0));
	gl::setMatrices(cam);

	const auto lambert = gl::ShaderDef().lambert();
	const auto shader = gl::getStockShader(lambert);
	shader->bind();

	gl::lineWidth(1);

	//gl::drawSphere(vec3(), 1.0f, 40);

	gl::pushModelMatrix();
	//gl::scale(2, 2, 2);
	gl::rotate(static_cast<float>(3.1415*2.0*getElapsedSeconds() / 8), vec3(0, 1, 0));
	gl::draw(mesh);

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
}



CINDER_APP(MeshSplitApp, RendererGl)

#endif