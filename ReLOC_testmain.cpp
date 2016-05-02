#include "ReLOC.h"

#include "ReLOC.h"

#include <stdexcept>


using namespace std;



int main(int argc, char** argv)
{
	string absDatasetPath = argv[1];	//folder containing the partial views to register
	float Rx = atof(argv[2]);			//radius (in unit of mesh resolution) of the support used for computing the X axis

	//load meshes
	vtkPolyDataCollection* cMeshes = ReLOC::ReadMeshes(absDatasetPath); 
	

	ReLOC::ReLOC reLOC;

	//Set parameters of ReLOC
	reLOC.Verbose(true);
	
	////use random detector
	//reLOC.UseRandomDetector(true);

	//set FLARE parameters
	reLOC.GetFLAREParams().radiusInmeshRes_tangent = Rx;
	//reLOC.GetFLAREParams().selector = "FLARE_Sampled_T";	//Perform a sub-sampling of the support used for the estimation of FLARE x axis (It reduces CPUtime without affecting meaningfully the number of correctly aligned view pairs)
	reLOC.GetFLAREParams().selector = "FLARE";

	////For High resolution datasets
	//reLOC.GetFlatPoints3DDetectorParams().radius = 8;
	//reLOC.GetFlatPoints3DDetectorParams().excludedPointsRadius = 8;
	//reLOC.GetFlatPoints3DDetectorParams().partitionRadius = 8;
	//reLOC.GetFlatPoints3DDetectorParams().coverageAreaPercentage_Hierarchical = 0.5;
	//reLOC.SetHoughBinSizeInMeshRes(6); 
	//reLOC.GetRansac().threshInUnitOfMeshRes = 7;

	//Apply ReLOC
	vector<vtkTransform*> vPairwiseRigidMotions = reLOC.Align(cMeshes);

	cout << "Registration completed" << endl;
}

