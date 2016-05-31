
#include "Pairwise3DRegistrationEvaluation.h"

#include <stdexcept>

#include <pcl/registration/ia_ReLOC.h>

using namespace std;

using namespace pcl;
using namespace pcl::io;
using namespace pcl::registration;


/** \brief Wrapper of Pcl version fo ReLOC algorithm.
*/
class PclReLOCTest : public Pairwise3DRegistrationEvaluation::IPairwise3DRegistrationAlgorithm
{
	public:
  ReLOCInitialAlignment <PointXYZ, PointXYZ> reLOC;
  int idx_trg_prev;
  int idx_ref_prev;

  PclReLOCTest():
    idx_trg_prev(-1),
    idx_ref_prev(-1)
  {};

	vtkTransform* PairwiseRegistration(vtkPolyData* cloud_trg, vtkPolyData* cloud_ref, const int idx_trg, const int idx_ref, double &secCPUtime)
	{
    if(idx_trg != idx_trg_prev)
    {
		  //access cloud_trg so to fill your data
		  float* ptrPoints_trg = Pairwise3DRegistrationEvaluation::GetPolyDataPointsPointer(cloud_trg);	//x0,y0,z0,x1,y1,z1,x2,y2,z2. etc...
		  vtkIdType nPoints_trg = Pairwise3DRegistrationEvaluation::GetPolyDataNumberOfPoints(cloud_trg);

		  float* ptrNormals_trg = Pairwise3DRegistrationEvaluation::GetPolyDataPointNormalsPointer(cloud_trg);	//nx0,ny0,nz0,nx1,ny1,nz1,nx2,ny2,nz2. etc...

      PointCloud <PointXYZ> target;
      target.resize(nPoints_trg);
      PointCloud <Normal> target_normals;
      target_normals.resize(nPoints_trg);
      for(int po=0; po<nPoints_trg; po++)
      {
        target.at(po).x = *ptrPoints_trg++;
        target.at(po).y = *ptrPoints_trg++;
        target.at(po).z = *ptrPoints_trg++;

        target_normals.at(po).normal_x = *ptrNormals_trg++;
        target_normals.at(po).normal_y = *ptrNormals_trg++;
        target_normals.at(po).normal_z = *ptrNormals_trg++;
      }

      reLOC.setInputTarget (target.makeShared());
      reLOC.setTargetNormals (target_normals.makeShared()); 
    }

    if(idx_ref != idx_ref_prev)
    {
      //access cloud_ref so to fill your data
		  float* ptrPoints_ref = Pairwise3DRegistrationEvaluation::GetPolyDataPointsPointer(cloud_ref);	//x0,y0,z0,x1,y1,z1,x2,y2,z2. etc...
		  vtkIdType nPoints_ref = Pairwise3DRegistrationEvaluation::GetPolyDataNumberOfPoints(cloud_ref);

		  float* ptrNormals_ref = Pairwise3DRegistrationEvaluation::GetPolyDataPointNormalsPointer(cloud_ref);	//nx0,ny0,nz0,nx1,ny1,nz1,nx2,ny2,nz2. etc...

      PointCloud <PointXYZ> source;
      source.resize(nPoints_ref);
      PointCloud <Normal> source_normals;
      source_normals.resize(nPoints_ref);
      for(int po=0; po<nPoints_ref; po++)
      {
        source.at(po).x = *ptrPoints_ref++;
        source.at(po).y = *ptrPoints_ref++;
        source.at(po).z = *ptrPoints_ref++;

        source_normals.at(po).normal_x = *ptrNormals_ref++;
        source_normals.at(po).normal_y = *ptrNormals_ref++;
        source_normals.at(po).normal_z = *ptrNormals_ref++;
      }

      reLOC.setInputSource (source.makeShared());
      reLOC.setSourceNormals (source_normals.makeShared()); 
    }

       
 

    //perform alignment
		Pairwise3DRegistrationEvaluation::CpuTimeProfiler cpuTime;
		
    PointCloud <PointXYZ> source_aligned;
		reLOC.align(source_aligned);

		secCPUtime = cpuTime.GetElapsedSecs();


    //get rigid motion
    Eigen::Matrix4f rigidMotion_4f = reLOC.getFinalTransformation();

		vtkTransform* rigidMotion = vtkTransform::New();
		rigidMotion->GetMatrix()->Element[0][0] = rigidMotion_4f(0,0);
		rigidMotion->GetMatrix()->Element[0][1] = rigidMotion_4f(0,1);
		rigidMotion->GetMatrix()->Element[0][2] = rigidMotion_4f(0,2);
		rigidMotion->GetMatrix()->Element[0][3] = rigidMotion_4f(0,3);
		rigidMotion->GetMatrix()->Element[1][0] = rigidMotion_4f(1,0);
		rigidMotion->GetMatrix()->Element[1][1] = rigidMotion_4f(1,1);
		rigidMotion->GetMatrix()->Element[1][2] = rigidMotion_4f(1,2);
		rigidMotion->GetMatrix()->Element[1][3] = rigidMotion_4f(1,3);
		rigidMotion->GetMatrix()->Element[2][0] = rigidMotion_4f(2,0);
		rigidMotion->GetMatrix()->Element[2][1] = rigidMotion_4f(2,1);
		rigidMotion->GetMatrix()->Element[2][2] = rigidMotion_4f(2,2);
		rigidMotion->GetMatrix()->Element[2][3] = rigidMotion_4f(2,3);
		rigidMotion->GetMatrix()->Element[3][0] = rigidMotion_4f(3,0);
		rigidMotion->GetMatrix()->Element[3][1] = rigidMotion_4f(3,1);
		rigidMotion->GetMatrix()->Element[3][2] = rigidMotion_4f(3,2);
		rigidMotion->GetMatrix()->Element[3][3] = rigidMotion_4f(3,3);


    idx_trg_prev = idx_trg;
    idx_ref_prev = idx_ref;

		return rigidMotion;	
	};

};





int main(int argc, char** argv)
{
	string absDatasetPath = argv[1];	//folder containing the partial views to register
	float Rx = atof(argv[2]);			//radius (in unit of mesh resolution) of the support used for computing the X axis
	string absResultFilename = argv[3];	//filename of the output csv file reporting the result of the registration

	//Instantiate the wrapper of the algorithm to test
	PclReLOCTest* ptrReLOCAlgorithm = new PclReLOCTest();

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

  double meshRes = benchmark.GetMeanMeshRes();
	
	//Set parameters of ReLOC
  //ptrReLOCAlgorithm->reLOC.useRandomDetector(true);  // use random detector
  ptrReLOCAlgorithm->reLOC.setSeed (0xFFFFFFFF);
  ptrReLOCAlgorithm->reLOC.setFlatKeypointRf (3.0*meshRes); //5 in the paper
  ptrReLOCAlgorithm->reLOC.setFlatKeypointRdiscard (2.0*meshRes);
  ptrReLOCAlgorithm->reLOC.setFlatKeypointR1search (2.0*meshRes);
  ptrReLOCAlgorithm->reLOC.setFlatKeypointR2search (20.0*meshRes);
  ptrReLOCAlgorithm->reLOC.setFlareNormalRadius (5.0*meshRes);
  ptrReLOCAlgorithm->reLOC.setFlareTangentRadius (Rx*meshRes);
  //ptrReLOCAlgorithm->reLOC.setFlareXsupportSamplingPerc(0.2);  //Perform a sub-sampling of the support used for the estimation of FLARE x axis (It reduces CPUtime without affecting meaningfully the number of correctly aligned view pairs)
  ptrReLOCAlgorithm->reLOC.setHoughSbin (2.0*meshRes);
  ptrReLOCAlgorithm->reLOC.setRansacT (8.0*meshRes);

   
  
 // //For High resolution datasets
 // ptrReLOCAlgorithm->reLOC.setFlatKeypointRf (8.0*meshRes);
 // ptrReLOCAlgorithm->reLOC.setFlatKeypointRdiscard (8.0*meshRes);
 // ptrReLOCAlgorithm->reLOC.setFlatKeypointR1search (8.0*meshRes);
 // ptrReLOCAlgorithm->reLOC.setFlatKeypointT2search (0.5);
 // ptrReLOCAlgorithm->reLOC.setHoughSbin (6.0*meshRes);
 // ptrReLOCAlgorithm->reLOC.setRansacT (7.0*meshRes);
	//benchmark.m_ICPsamplingPerc = 0.05;

	//Run the benchmark on all the view pairs
	benchmark.Evaluate();

	benchmark.PrintResults(absResultFilename, benchmark.GetNRegistrations(), benchmark.GetRMSE(), benchmark.GetTotalCPUtime() );

	delete ptrReLOCAlgorithm;
}

