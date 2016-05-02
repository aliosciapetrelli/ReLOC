#include "ReLOC.h"

#include "ReLOC_utils.h"


#include <stdexcept>

using namespace::std;


ReLOC::ReLOC::ReLOC():
	m_cMeshes(NULL),
	m_idx_trg_prev(-1),
	m_meanMeshRes(-1),
	m_verbose(false),
	m_useRandomDetector(false),
	m_ptrDetector(NULL),
	m_ptrDescriptor(NULL),
	m_ptrMatcher(NULL),
	m_houghSpace(NULL),
	m_houghSigmaFactor(2.0f),
	m_houghBinSizeInMeshRes(2.0f),
	m_houghBBoxSizeFactor(1.4f)
{
}


ReLOC::ReLOC::~ReLOC()
{
	deallocate();
}

void ReLOC::ReLOC::deallocate()
{
	DeleteTransforms(m_vPairwiseTransforms);	

	if(m_ptrDetector)
	{
		delete m_ptrDetector;
		m_ptrDetector = NULL;
	}

	if(m_ptrDescriptor)
	{
		delete m_ptrDescriptor;
		m_ptrDescriptor = NULL;
	}

	if(m_ptrMatcher)
	{
		delete m_ptrMatcher;
		m_ptrMatcher = NULL;
	}

	DestroyHoughSpace(m_houghSpace);

	m_vViewPairs.clear();

	m_meanMeshRes = -1.0;
}


void ReLOC::ReLOC::PrepareRegistration(vtkPolyDataCollection* cMeshes, const double meanMeshRes, std::vector< std::pair<int,int> > &vViewPairs)
{
	deallocate();

	m_cMeshes = cMeshes;
	int nMeshes = m_cMeshes->GetNumberOfItems();

	//Compute mesh resolution
	m_meanMeshRes = meanMeshRes;
	if(m_meanMeshRes == -1)
	{
		m_meanMeshRes = ComputeMeshResolution(m_cMeshes);
	}




	//Detection
	if(m_useRandomDetector)
	{
		m_randomDetectorParams.radius = m_randomDetectorParams.radius * m_meanMeshRes;
		m_ptrDetector = new Random3DDetector(m_randomDetectorParams);
	}
	else
	{
		m_flatPoints3DDetectorParams.radius = m_flatPoints3DDetectorParams.radius * m_meanMeshRes;
		m_flatPoints3DDetectorParams.partitionRadius = m_flatPoints3DDetectorParams.partitionRadius * m_meanMeshRes;
		m_flatPoints3DDetectorParams.excludedPointsRadius = m_flatPoints3DDetectorParams.excludedPointsRadius * m_meanMeshRes;
		m_flatPoints3DDetectorParams.partitionRadius_Hierarchical = m_flatPoints3DDetectorParams.partitionRadius_Hierarchical * m_meanMeshRes;
		m_ptrDetector = new FlatPoints3DDetector(m_flatPoints3DDetectorParams);
	}
	
	//extract features
	for(int fi=0; fi<nMeshes; fi++)
	{
		if(m_verbose)
		{
			cout << "Detection - Mesh: " << fi << "/" << nMeshes-1 << endl;
		}
		vtkPolyData* mesh = (vtkPolyData*)(cMeshes->GetItemAsObject(fi));
		

		Feature3D* feat;
		int numFeats = m_ptrDetector->extract(mesh, feat);

		//push features of current mesh in FeaturesVect
		m_ptrDetector->getNumFeaturesVect().push_back(numFeats);
		if(numFeats>0)
		{
			Feature3D* feat3D = new Feature3D[numFeats];
			memcpy(feat3D, m_ptrDetector->getFeatures(), sizeof(Feature3D)*numFeats);
			m_ptrDetector->getFeaturesVect().push_back(feat3D);
		}
		else
		{
			m_ptrDetector->getFeaturesVect().push_back(NULL);
		}

	}


	//Description
	m_LRFDescriptor_1_Params.meshRes = m_meanMeshRes;
	m_ptrDescriptor = new LRFDescriptor_1(m_LRFDescriptor_1_Params);

	//perform description
	m_vDescs.clear();
	for(int fi=0; fi<nMeshes; fi++)
	{
		if(m_verbose)
		{
			cout << "Description - Mesh: " << fi << "/" << nMeshes-1 << endl;
		}
		vtkPolyData* mesh = (vtkPolyData*)(cMeshes->GetItemAsObject(fi));

		Desc desc;

		Feature3D* feat = m_ptrDetector->getFeaturesVect()[fi];
		int numFeats = m_ptrDetector->getNumFeaturesVect()[fi];

		m_ptrDescriptor->describe(mesh, feat, desc, numFeats );

		m_vDescs.push_back(desc);
	}



	//Matching
	m_ptrMatcher = new LRFMatcher_1(m_LRFMatcher_1_Params);

	m_vViewPairs = vViewPairs;
	if(m_vViewPairs.size() == 0)
	{
		CalcAllViewPairs(m_vViewPairs, m_cMeshes->GetNumberOfItems());
	}

	m_idx_trg_prev = -1;

	m_vPairwiseTransforms.resize(m_vViewPairs.size(), NULL);


	//Hough voting
	GetCentroids_cMeshes(m_cMeshes, m_vMeshCentroids);				//compute centroids of each partial view
	GetBBox_StdDev_Clouds(m_cMeshes, m_vMeshBBoxes, m_vMeshCentroids, m_houghSigmaFactor); //estimate hough space extension for each partial view
	GetCentroidsWrtLRFs(m_vMeshCentroids, m_ptrDetector->getFeaturesVect(), m_ptrDetector->getNumFeaturesVect(), m_vvCentroidsWrtLRF);	//express each centroid w.r.t. each LRF for each partial view 

	//ransac
	m_RANSAC.Allocate(m_meanMeshRes);

}

std::vector<vtkTransform*> &ReLOC::ReLOC::Align(vtkPolyDataCollection* cMeshes, const double meanMeshRes, std::vector< std::pair<int,int> > &vViewPairs)
{
	PrepareRegistration(cMeshes, meanMeshRes, vViewPairs);

	for(size_t pa=0; pa<m_vViewPairs.size(); pa++)
	{
		int idx_trg = m_vViewPairs[pa].first;
		int idx_ref = m_vViewPairs[pa].second;

		if(m_verbose)
			cout << "Align pair: " << pa << " / " <<  m_vViewPairs.size() << "\ttrg: " << idx_trg << "\tref: " << idx_ref << " ... ";
		
		m_vPairwiseTransforms[pa] = Align(idx_trg, idx_ref);
	}

	DestroyHoughSpace(m_houghSpace);

	return m_vPairwiseTransforms;
}


vtkTransform* ReLOC::ReLOC::Align(const int idx_trg, const int idx_ref)
{

	if(m_vDescs[idx_ref].getNumDesc() == 0)
	{
		if(m_verbose)
			cout << "No ref Desc" << endl;
		return NULL;
	}

	if(m_vDescs[idx_trg].getNumDesc() == 0)
	{
		if(m_verbose)
			cout << "No trg Desc" << endl;
		return NULL;
	}

	FeatureMatch* matches = NULL;
	int nMatches = 0;
	m_ptrMatcher->setReference(m_vDescs[idx_ref]);
	nMatches = m_ptrMatcher->match(m_vDescs[idx_trg], matches);			
	if(nMatches == 0)
	{
		if(m_verbose)
			cout << "No matches" << endl;
		return NULL;
	}

	const Feature3D* feat3D_trg = m_ptrDetector->getFeaturesVect()[idx_trg];
	const int nFeat3D_trg = m_ptrDetector->getNumFeaturesVect()[idx_trg];
	const Feature3D* feat3D_ref = m_ptrDetector->getFeaturesVect()[idx_ref];
	const int nFeat3D_ref = m_ptrDetector->getNumFeaturesVect()[idx_ref];


	vector< MaximumContainer<unsigned short int> > vMaxima_afterRansac;

	//hough voting

	//only if mesh_trg has changed allocate a new hough space
	if(m_idx_trg_prev != idx_trg)
	{
		DestroyHoughSpace(m_houghSpace);
		CreateHoughSpace(m_houghSpace, &m_vMeshBBoxes[idx_trg*6], m_houghBinSizeInMeshRes * m_meanMeshRes, m_houghBBoxSizeFactor);
	}
	m_idx_trg_prev = idx_trg;

	//given the coordinates of the centroid of mesh_ref represented w.r.t. LRFs of mesh_ref, represents the coordinates w.r.t the LRFs of mesh_trg
	vector<float> transformedCentroidSpace;
	vector<unsigned short int> vWeights;
	CreateTransformedCentroidSpace(matches, nMatches, feat3D_trg, m_vvCentroidsWrtLRF[idx_ref], transformedCentroidSpace, vWeights, false );

	vector< MaximumContainer<unsigned short int> > vMaxima;
	PopulateHoughSpaceAndFindMaxima(m_houghSpace, transformedCentroidSpace, vWeights, vMaxima);

	vector<FeatureMatch> houghMatches;
	ObtainHoughMatches(vMaxima, matches, houghMatches);

	if(!houghMatches.size())
	{
		if(m_verbose)
			cout << "No matches after Hough voting" << endl;
		return NULL;
	}


	//RANSAC
	vtkTransform* ransacTransform = m_RANSAC.Apply(houghMatches, feat3D_trg, nFeat3D_trg, feat3D_ref, nFeat3D_ref);

	if( (!m_RANSAC.outTransform) || (m_RANSAC.consensusSetSize<3) )
	{
		if(m_verbose)
			cout << "RANSAC failed: No consensus set" << endl;
		return NULL;
	}

	if(m_verbose)
		cout << "Consensus set size: " << m_RANSAC.consensusSetSize << endl;

	return Duplicate(m_RANSAC.outTransform);

}

vtkTransform* ReLOC::ReLOC::Align(vtkPolyData* mesh_trg, vtkPolyData* mesh_ref)
{
	vtkPolyDataCollection* cMeshes = vtkPolyDataCollection::New();
	cMeshes->AddItem(mesh_trg);
	cMeshes->AddItem(mesh_ref);

	vector<vtkTransform*> vTransforms = Align(cMeshes);

	return vTransforms[0];
}