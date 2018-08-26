#include "cinder/app/App.h"
#include "cinder/app/RendererGl.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/Shader.h"


#include "voro++.hh"
#include "vorohelpers.h"
#include <random>


using namespace ci;
using namespace ci::app;
using namespace std;

std::default_random_engine re;

class MeshSplitApp : public App {
  public:
	  MeshSplitApp();
	void setup() override;
	void mouseDown( MouseEvent event ) override;
	void update() override;
	void draw() override;
	void keyDown(KeyEvent event) override;

private:
	TriMesh mesh;
	voro::container con;

	std::vector<std::vector<Face>> cells;
	static const int MAX_SIZE = 5;

	int currentCell = 0;

private:

};

MeshSplitApp::MeshSplitApp() :con(-MAX_SIZE, MAX_SIZE, -MAX_SIZE, MAX_SIZE, -MAX_SIZE, MAX_SIZE, 1, 1, 1, false, false, false, 10)
{
	
}


void MeshSplitApp::setup()
{
	mesh = geom::Teapot();

	uniform_real_distribution<double> offsetDist(-0.5, 0.5);
	int particleID = 0;
	for(int x=-1; x<=1;x++)
		for(int y=-1;y<=1;y++)
			for (int z = -1; z <= 1; z++)
		{
				con.put(particleID++, x+offsetDist(re), y + offsetDist(re), z + offsetDist(re));
		}

	con.compute_all_cells();

	// Loop over all particles and save computed voronoi cells
	voro::c_loop_all vLoop(con);
	vLoop.start();
	do
	{
		voro::voronoicell vcell;

		if(con.compute_cell(vcell, vLoop))
		{
			cells.emplace_back(getFacesFromEdges(vcell));
		}

	} while (vLoop.inc());


}

void MeshSplitApp::mouseDown( MouseEvent event )
{
}

void MeshSplitApp::update()
{

}




void MeshSplitApp::draw()
{
	gl::enableDepth(true);
	gl::clear(Color(0.5, 0.5, 0.5));
	gl::enableFaceCulling(false);

	CameraPersp cam;
	cam.lookAt(vec3(3, 3, 3), vec3(0));
	gl::setMatrices(cam);

	const auto lambert = gl::ShaderDef().lambert();
	const auto shader = gl::getStockShader(lambert);
	shader->bind();

	gl::lineWidth(3);
	gl::color(ColorA(1.0f, 0.0f, 0.0f, 0.5f));

	//gl::drawSphere(vec3(), 1.0f, 40);

	gl::pushModelMatrix();
	gl::rotate(3.1415f*2.f*getElapsedSeconds() / 8, vec3(0,1,0));
	gl::draw(mesh);

	gl::popModelMatrix();

}

void MeshSplitApp::keyDown(KeyEvent event)
{
	if(event.getChar()=='+')
	{
		currentCell++;
		mesh = meshFromFaces(cells[currentCell%cells.size()]);
	}
}



CINDER_APP( MeshSplitApp, RendererGl )
