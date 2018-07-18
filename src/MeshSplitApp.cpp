#include "cinder/app/App.h"
#include "cinder/app/RendererGl.h"
#include "cinder/gl/gl.h"

using namespace ci;
using namespace ci::app;
using namespace std;

class MeshSplitApp : public App {
  public:
	void setup() override;
	void mouseDown( MouseEvent event ) override;
	void update() override;
	void draw() override;
};

void MeshSplitApp::setup()
{
}

void MeshSplitApp::mouseDown( MouseEvent event )
{
}

void MeshSplitApp::update()
{
}

void MeshSplitApp::draw()
{
	gl::clear( Color( 0, 0, 0 ) ); 
}

CINDER_APP( MeshSplitApp, RendererGl )
