#include "ReLOC.h"

#include "Pairwise3DRegistrationEvaluation.h"
#include "ReLOC.h"

#include <stdexcept>


using namespace std;



/** \brief Wrapper of ReLOC algorithm. Implement PairwiseRegistration() method.
*/
class ReLOCTest : public Pairwise3DRegistrationEvaluation::IPairwise3DRegistrationAlgorithm
{
	public:
	ReLOC::ReLOC reLOC;

	vtkTransform* PairwiseRegistration(vtkPolyData* cloud_trg, vtkPolyData* cloud_ref, const int idx_trg, const int idx_ref, double &secCPUtime)
	{
		Pairwise3DRegistrationEvaluation::CpuTimeProfiler cpuTime;
		
		vtkTransform* rigidMotion = reLOC.Align(idx_trg, idx_ref);

		secCPUtime = cpuTime.GetElapsedSecs();

		return rigidMotion;	
	};

};



int main(int argc, char** argv)
{
	string absDatasetPath = argv[1];	//folder containing the partial views to register
	float Rx = atof(argv[2]);			//radius (in unit of mesh resolution) of the support used for computing the X axis
	string absResultFilename = argv[3];	//filename of the output csv file reporting the result of the registration

	//Instantiate the wrapper of the algorithm to test
	ReLOCTest* ptrReLOCAlgorithm = new ReLOCTest();

	//Instantiate the benchmark
	Pairwise3DRegistrationEvaluation::Pairwise3DRegistrationBenchmark benchmark; 

	//Set the parameters of the benchmark
	benchmark.SetDatasetPath(absDatasetPath);
	benchmark.SetAlgorithm(ptrReLOCAlgorithm);
	//optional params
	//vector<pair<int, int>> vViewPairs;
	//vViewPairs.push_back(pair<int,int>(1,2));
	//vViewPairs.push_back(pair<int,int>(1,3));
	//vViewPairs.push_back(pair<int,int>(2,3));
	//benchmark.SetViewPairs(vViewPairs);
	//benchmark.SetMeshFileExtension("ply");
	//benchmark.m_lowMemory = false;
	//benchmark.m_RMSEthresh = 5.0;
	//benchmark.m_skipICP = false;
	//benchmark.m_ICPnMaxIters = 10;
	//benchmark.m_ICPeps = 0.1;
	//benchmark.m_ICPmaxRadius = 8.0;
	//benchmark.m_ICPsamplingPerc = 1.0;

	//Load partial views and ground truth
	benchmark.PrepareEvaluation();

	
	//Set parameters of ReLOC
	ptrReLOCAlgorithm->reLOC.Verbose(true);
	

	////use random detector
	//ptrReLOCAlgorithm->reLOC.UseRandomDetector(true);

	//set FLARE parameters
	ptrReLOCAlgorithm->reLOC.GetFLAREParams().radiusInmeshRes_tangent = Rx;
	//ptrReLOCAlgorithm->reLOC.GetFLAREParams().selector = "FLARE_Sampled_T";	//Perform a sub-sampling of the support used for the estimation of FLARE x axis (It reduces CPUtime without affecting meaningfully the number of correctly aligned view pairs)
	ptrReLOCAlgorithm->reLOC.GetFLAREParams().selector = "FLARE";

	////For High resolution datasets
	//ptrReLOCAlgorithm->reLOC.GetFlatPoints3DDetectorParams().radius = 8;
	//ptrReLOCAlgorithm->reLOC.GetFlatPoints3DDetectorParams().excludedPointsRadius = 8;
	//ptrReLOCAlgorithm->reLOC.GetFlatPoints3DDetectorParams().partitionRadius = 8;
	//ptrReLOCAlgorithm->reLOC.GetFlatPoints3DDetectorParams().coverageAreaPercentage_Hierarchical = 0.5;
	//ptrReLOCAlgorithm->reLOC.SetHoughBinSizeInMeshRes(6); 
	//ptrReLOCAlgorithm->reLOC.GetRansac().threshInUnitOfMeshRes = 7;
	//benchmark.m_ICPsamplingPerc = 0.05;

	//Start ReLOC registration. Compute detection and description 
	Pairwise3DRegistrationEvaluation::CpuTimeProfiler cpuTimeProfiler_Prepare;
	ptrReLOCAlgorithm->reLOC.PrepareRegistration(benchmark.GetMeshes(), benchmark.GetMeanMeshRes() );
	double cpuPrepare = cpuTimeProfiler_Prepare.GetElapsedSecs();
	
	//Run the benchmark on all the view pairs
	benchmark.Evaluate();

	benchmark.PrintResults(absResultFilename, benchmark.GetNRegistrations(), benchmark.GetRMSE(), benchmark.GetTotalCPUtime() + cpuPrepare );

	delete ptrReLOCAlgorithm;
}

