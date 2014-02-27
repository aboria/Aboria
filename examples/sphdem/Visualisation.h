/*
 * visualisation.h
 *
 *  Created on: 1 Feb 2014
 *      Author: mrobins
 */

#ifndef VISUALISATION_H_
#define VISUALISATION_H_

#include <vtkVersion.h>
#include <vtkProperty.h>
#include <vtkPlaneSource.h>
#include <vtkCubeSource.h>
#include <vtkVertexGlyphFilter.h>
#include "vtkActor2D.h"
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRegularPolygonSource.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkAxesActor.h>
#include <vtkCamera.h>
#include <vtkCallbackCommand.h>

#include <vtkAxis.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <thread>

void TimerCallbackFunction ( vtkObject* caller,
							long unsigned int vtkNotUsed(eventId),
							void* clientData,
							void* vtkNotUsed(callData) ) {
  vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =
      static_cast<vtkVertexGlyphFilter*>(clientData);

  vtkRenderWindowInteractor *iren =
    static_cast<vtkRenderWindowInteractor*>(caller);
  vertexFilter->Modified();
  iren->Render();

}

void interactor_loop(vtkSmartPointer<vtkRenderWindowInteractor> arg) {
	std::cout << "starting interactor_loop"<<std::endl;
	arg->Start();
}

class Visualisation {
public:
	Visualisation(const Vect3d& min, const Vect3d& max) {

		renderer = vtkSmartPointer<vtkRenderer>::New();
		renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		vertexGlyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();

		renderWindow->AddRenderer(renderer);
		renderWindowInteractor->SetRenderWindow(renderWindow);

		renderWindow->SetSize(800,600);
		renderer->SetBackground(.1,.2,.3); // Background color dark blue

		vtkSmartPointer<vtkAxesActor> axes =
				vtkSmartPointer<vtkAxesActor>::New();

		orientationWidget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
		//orientationWidget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
		orientationWidget->SetOrientationMarker( axes );
		orientationWidget->SetInteractor( renderWindowInteractor );
		//orientationWidget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
		orientationWidget->SetEnabled( 1 );

		// Create a cube.
		vtkSmartPointer<vtkCubeSource> cubeSource =
				vtkSmartPointer<vtkCubeSource>::New();
		cubeSource->SetBounds(min[0],max[0],min[1],max[1],min[2],max[2]);

		// Create a mapper and actor.
		vtkSmartPointer<vtkPolyDataMapper> mapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(cubeSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor =
				vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetRepresentationToWireframe();
		renderer->AddActor(actor);

		const Vect3d mid = 0.5*(max+min);
		renderer->GetActiveCamera()->SetFocalPoint(mid[0],mid[1],mid[2]);
		renderer->GetActiveCamera()->SetPosition(mid[0],mid[1] + max[1]-min[1],mid[2]);
		renderer->GetActiveCamera()->SetViewUp(0,0,1);
		renderer->ResetCamera();

		renderWindow->Render();
		renderWindowInteractor->Initialize();

		//renderWindowInteractor->SetStillUpdateRate(100);
		renderWindowInteractor->CreateRepeatingTimer(200);
		vtkSmartPointer<vtkCallbackCommand> timerCallback =
				vtkSmartPointer<vtkCallbackCommand>::New();
		timerCallback->SetCallback ( TimerCallbackFunction );
		timerCallback->SetClientData(vertexGlyphFilter);
		renderWindowInteractor->AddObserver ( vtkCommand::TimerEvent, timerCallback );
	}
	void glyph_points(vtkSmartPointer<vtkUnstructuredGrid> grid) {

#if VTK_MAJOR_VERSION <= 5
		vertexGlyphFilter->AddInput(grid);
#else
		vertexGlyphFilter->AddInputData(grid);
#endif
		vertexGlyphFilter->Update();

		vtkSmartPointer<vtkPolyDataMapper> mapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(vertexGlyphFilter->GetOutputPort());
		mapper->Update();

		vtkSmartPointer<vtkActor> actor =
				vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		actor->GetProperty()->SetPointSize(10);
		renderer->AddActor(actor);
		renderWindow->Render();
	}

	void start_render_loop() {
		//std::thread render_thread(interactor_loop,renderWindowInteractor);
		//render_thread.detach();
	}
	void restart_render_loop() {
		vertexGlyphFilter->Modified();
		renderWindowInteractor->EnableRenderOn();
	}
	void stop_render_loop() {
		renderWindowInteractor->EnableRenderOff();
	}

	static void vtkWriteGrid(const char *name, int timestep, vtkSmartPointer<vtkUnstructuredGrid> grid) {

		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
				vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

#if VTK_MAJOR_VERSION > 5
		writer->SetInputData(grid);
#else
		writer->SetInput(grid);
#endif

		writer->SetDataModeToBinary();

		char buffer[100];
		sprintf(buffer,"%s%05d.vtu",name,timestep);

		writer->SetFileName(buffer);
		writer->Write();
	}


private:
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindow> renderWindow;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
	vtkSmartPointer<vtkOrientationMarkerWidget> orientationWidget;
	vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter;
};


#endif /* VISUALISATION_H_ */
