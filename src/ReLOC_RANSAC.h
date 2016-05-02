#ifndef RELOC_RANSAC_H
#define RELOC_RANSAC_H


#include "ReLOC_LRF.h"

#include <vector>
#include <limits>

#include "vtkPolyData.h"
#include "vtkPolyDataCollection.h"
#include "vtkTransform.h"


namespace ReLOC
{

	class AbsoluteOrientation
	{
	protected:
		vtkTransform*		m_transform;
		double				m_matrTransform[16];
		double				m_rmse;
		bool				m_calcRMSE;
		bool				m_calcScale;

	public:
		AbsoluteOrientation(bool calcScale = true, bool calcRMSE = false);
		virtual				~AbsoluteOrientation();

		vtkTransform*		GetTransform(){return m_transform;}
		double*				GetMatrTransform(){return m_matrTransform;}

		//every polydata must be allocated, vtkTransform mustn't be deallocated
		virtual vtkTransform*		Orientate(vtkPolyData* polyData1In, vtkPolyData* polyData2In);

		static int ExtractPoints(vtkPolyData* poly, Feature3D* features3D, int numFeatures3D, float* &outPoints );
		static vtkPolyData* ExtractPoints(vtkPolyData* poly, Feature3D* features3D, int numFeatures3D );
		static bool ExtractPoints(FeatureMatch* featureMatches, int numFeatureMatches, Feature3D* trgFeatures3D, int numTrgFeatures3D, Feature3D* refFeatures3D, int numRefFeatures3D, float* &trgPoints, float* &refPoints);

		int	Orientate(FeatureMatch* featureMatches, int numFeatureMatches, Feature3D* trgFeatures3D, int numTrgFeatures3D, Feature3D* refFeatures3D, int numRefFeatures3D, double* &matrTransform);

		vtkTransform* Orientate(FeatureMatch* featureMatches, int numFeatureMatches, Feature3D* trgFeatures3D, int numTrgFeatures3D, Feature3D* refFeatures3D, int numRefFeatures3D);

		virtual int Orientate(float *points1, float* points2, int nPoints, double* &matrTransform) = 0;

		static double CalcRMSE(float* points1, float* points2, int nPoints);
		static double CalcRMSE(vtkPolyData* polyData1, vtkPolyData* polyData2);

		double GetRMSE(){return m_rmse;}
		void CalcScale(bool calcScale){m_calcScale = calcScale;}
		void CalcRMSE(bool calcRMSE){m_calcRMSE = calcRMSE;}
		bool GetCalcScale(){return m_calcScale;}
		bool GetCalcRMSE(){return m_calcRMSE;}
	};


	class HornOrientation : public AbsoluteOrientation
	{
	public:
		HornOrientation(bool calcScale = true, bool calcRMSE = false);

		virtual int Orientate(float *points1, float* points2, int nPoints, double* &matrTransform);

	};


	class PairwiseRegistration
	{
	protected:
		vtkTransform*		m_transform;
		double				m_rmse;
		bool				m_calcRMSE;
		bool				m_calcScale;
		bool				m_verbose;

	public:

		PairwiseRegistration(bool calcScale = true, bool calcRMSE = false, bool verbose = false);
		virtual ~PairwiseRegistration();

		vtkTransform*		GetTransform(){return m_transform;}

		vtkTransform*	Register(vtkPolyData* polyData1In, vtkPolyData* polyData2In){return m_transform;};

		static double CalcRMSE(float* points1, float* points2, int nPoints);
		static double CalcRMSE(vtkPolyData* polyData1, vtkPolyData* polyData2);

		double GetRMSE(){return m_rmse;}
		void CalcScale(bool calcScale){m_calcScale = calcScale;}
		void CalcRMSE(bool calcRMSE){m_calcRMSE = calcRMSE;}
		bool GetCalcScale(){return m_calcScale;}
		bool GetCalcRMSE(){return m_calcRMSE;}
	};


	class FeaturePointsPairwiseRegistration : public PairwiseRegistration
	{
	public:

		FeaturePointsPairwiseRegistration(bool m_calcScale = true, bool calcRMSE = false, bool verbose = false);

		virtual vtkTransform*	Register(vtkPolyData* polyData1In, vtkPolyData* polyData2In) = 0;

		virtual vtkTransform*	Register(vtkPolyData* polyData1In, vtkPolyData* polyData2In, std::vector<vtkIdType> &polyPata1Indices, std::vector<vtkIdType> &polyPata2Indices);

		static bool CreatePolyData(FeatureMatch* featureMatches, int numFeatureMatches, Feature3D* trgFeatures3D, int numTrgFeatures3D, Feature3D* refFeatures3D, int numRefFeatures3D, vtkPolyData* &trgPolyData, vtkPolyData* &refPolyData);

		virtual vtkTransform*	Register(FeatureMatch* featureMatches, int numFeatureMatches, Feature3D* trgFeatures3D, int numTrgFeatures3D, Feature3D* refFeatures3D, int numRefFeatures3D);

		static bool CreatePolyData(vtkPolyDataCollection* polyDataCollection, const int cloudToRegisterIdx, int* nearCloudMatchIdx, int numNearCloudMatches, const std::vector<Feature3D*> &vFeatures3D, const std::vector<int> &vNFeatures3D, const std::vector<FeatureMatch*> &vFeatureMatches, const std::vector<int> &vNFeatureMatches, const std::vector<std::pair<int, int> > &vMatchPairs, vtkPolyData* &nearCloudsFeaturePoints, vtkPolyData* &cloudToRegisterFeaturePoints, const float threshDistance);

		virtual vtkTransform*	Register(vtkPolyDataCollection* polyDataCollection, const int cloudToRegisterIdx, int* nearCloudMatchIdx, int numNearCloudMatches, const std::vector<Feature3D*> &vFeatures3D, const std::vector<int> &vNFeatures3D, const std::vector<FeatureMatch*> &vFeatureMatches, const std::vector<int> &vNFeatureMatches, const std::vector<std::pair<int, int> > &vMatchPairs, const float threshDistance );
	};


	class RANSACPairwiseRegistration : public FeaturePointsPairwiseRegistration
	{
	protected:
		AbsoluteOrientation						&m_absOrientation;
		float									m_outlierDistance;
		int										m_maxConsensusSetSize;
		double									m_probGoodResult;
		int										m_nMaxIter;
		int*									m_maxConsensusSet;

	public:
		RANSACPairwiseRegistration(AbsoluteOrientation &absOrientation, float outlierDistance = 0.5, float probGoodResult = 0.99, int nMaxIter = std::numeric_limits<int>::max(), bool calcRMSE = false, bool verbose = false);
		virtual  ~RANSACPairwiseRegistration();
		virtual vtkTransform*	Register(vtkPolyData* polyData1In, vtkPolyData* polyData2In);
		virtual float GetOutlierDistance(){return m_outlierDistance;}
		virtual void SetOutlierDistance(float outlierDistance){m_outlierDistance = outlierDistance;}
		virtual int GetMaxConsensusSetSize(){return m_maxConsensusSetSize;}
		virtual int* GetMaxConsensusSet(){return m_maxConsensusSet;}
		virtual void CalcScale(bool calcScale){m_calcScale = calcScale; m_absOrientation.CalcScale(calcScale);}
	};


	class Ransac
	{
	public:
		FeaturePointsPairwiseRegistration*	ptr;
		AbsoluteOrientation*			aoptr;

		float							threshInUnitOfMeshRes;
		int								maxIters;
		double							prob;

		std::vector<FeatureMatch>		outMatches;
		int								consensusSetSize;
		int*							ptrMaxConsensusSet;
		vtkTransform*					outTransform;

		bool							verbose;

		Ransac(const float ransacThreshInUnitOfMeshRes = 8.0f, const int ransacMaxIters = 1000, const double ransacProb = 0.99);

		void Allocate(const double meanMeshRes);

		vtkTransform* Apply(const std::vector<FeatureMatch> &matches, const Feature3D* feat3D_trg, const int nFeat3D_trg, const Feature3D* feat3D_ref, const int nFeat3D_ref);

		~Ransac();
	};

}


#endif



