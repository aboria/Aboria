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

#include <vtkAxis.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>


class Visualisation {
public:
	Visualisation() {
		renderer = vtkSmartPointer<vtkRenderer>::New();
		renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

		renderWindow->SetSize(800,600);

		renderWindow->AddRenderer(renderer);
		renderWindowInteractor->SetRenderWindow(renderWindow);

		renderer->SetBackground(.1,.2,.3); // Background color dark blue

		vtkSmartPointer<vtkAxesActor> axes =
				vtkSmartPointer<vtkAxesActor>::New();

		orientationWidget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
		//orientationWidget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
		orientationWidget->SetOrientationMarker( axes );
		orientationWidget->SetInteractor( renderWindowInteractor );
		//orientationWidget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
		orientationWidget->SetEnabled( 1 );

		renderWindowInteractor->Initialize();

		//orientationWidget->InteractiveOn();

		//renderWindowInteractor->Initialize();
	}
	void glyph_points(vtkSmartPointer<vtkUnstructuredGrid> grid) {
		vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
				vtkSmartPointer<vtkVertexGlyphFilter>::New();
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
		actor->GetProperty()->SetPointSize(3);
		renderer->AddActor(actor);
	}

	void start_render_loop() {
		renderWindowInteractor->Start();
		renderWindowInteractor->EnableRenderOn();
	}
	void stop_render_loop() {
		renderWindowInteractor->EnableRenderOff();

	}

private:
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindow> renderWindow;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
	vtkSmartPointer<vtkOrientationMarkerWidget> orientationWidget;
};


#endif /* VISUALISATION_H_ */
