#include "ReLOC_RANSAC.h"

#include "ReLOC_utils.h"

#include <algorithm>


using namespace std;

#include "vtkMath.h"



//******************** PairwiseRegistration ****************


ReLOC::PairwiseRegistration::PairwiseRegistration(bool calcScale, bool calcRMSE, bool verbose)
	:m_calcRMSE(calcRMSE)
	,m_calcScale(calcScale)
	,m_verbose(verbose)
{
	m_transform = vtkTransform::New();
	m_rmse = -1.0;
}

double ReLOC::PairwiseRegistration::CalcRMSE(float* points1, float* points2, int nPoints)
{
	double rmse = calcDistance2BetweenPoints(points1, points2, nPoints);
	rmse /=(double)nPoints;
	rmse = sqrt(rmse);

	return rmse;
}

double ReLOC::PairwiseRegistration::CalcRMSE(vtkPolyData* polyData1, vtkPolyData* polyData2)
{
	if(polyData1->GetNumberOfPoints() != polyData2->GetNumberOfPoints())
	{
		return -1.0;
	}

	return CalcRMSE( (float*)polyData1->GetPoints()->GetVoidPointer(0), (float*)polyData2->GetPoints()->GetVoidPointer(0), polyData1->GetNumberOfPoints() );
}

ReLOC::PairwiseRegistration::~PairwiseRegistration()
{
	if(m_transform)
	{
		m_transform->Delete();
		m_transform = NULL;
	}
}




//******************** FeaturePointsPairwiseRegistration ****************


ReLOC::FeaturePointsPairwiseRegistration::FeaturePointsPairwiseRegistration(bool calcScale, bool calcRMSE, bool verbose)
	:PairwiseRegistration(calcScale, calcRMSE, verbose)
{
}


vtkTransform* ReLOC::FeaturePointsPairwiseRegistration::Register(vtkPolyData* polyData1In, vtkPolyData* polyData2In, std::vector<vtkIdType> &polyPata1Indices, std::vector<vtkIdType> &polyPata2Indices)
{
	if(polyPata1Indices.size() != polyPata2Indices.size())
	{
		return NULL;
	}

	vtkPolyData* poly1 = ReLOC::CreatePolyData(polyData1In, &polyPata1Indices[0], (vtkIdType)polyPata1Indices.size());
	vtkPolyData* poly2 = ReLOC::CreatePolyData(polyData2In, &polyPata2Indices[0], (vtkIdType)polyPata2Indices.size());

	m_transform = Register(poly1, poly2);

	poly1->Delete();
	poly2->Delete();

	return NULL;
}

vtkTransform* ReLOC::FeaturePointsPairwiseRegistration::Register(FeatureMatch* featureMatches, int numFeatureMatches, Feature3D* trgFeatures3D, int numTrgFeatures3D, Feature3D* refFeatures3D, int numRefFeatures3D)
{
	vtkPolyData* refPolyData;
	vtkPolyData* trgPolyData;
	if( !CreatePolyData(featureMatches, numFeatureMatches, trgFeatures3D, numTrgFeatures3D, refFeatures3D, numRefFeatures3D, trgPolyData, refPolyData) )
	{
		return NULL;
	}

	m_transform = Register(trgPolyData, refPolyData );

	refPolyData->Delete();
	trgPolyData->Delete();

	return m_transform;
}


bool ReLOC::FeaturePointsPairwiseRegistration::CreatePolyData(FeatureMatch* featureMatches, int numFeatureMatches, Feature3D* trgFeatures3D, int numTrgFeatures3D, Feature3D* refFeatures3D, int numRefFeatures3D, vtkPolyData* &trgPolyData, vtkPolyData* &refPolyData)
{

	refPolyData = NULL;
	trgPolyData = NULL;
	if(numFeatureMatches==0)
	{
		return false;
	}

	refPolyData = vtkPolyData::New();
	trgPolyData = vtkPolyData::New();

	float* ptrPolyRef = AllocatePoints(refPolyData, numFeatureMatches);
	float* ptrPolyTrg = AllocatePoints(trgPolyData, numFeatureMatches);

	FeatureMatch* ptrFeatureMatch = &featureMatches[0];
	for(int po=0; po<numFeatureMatches; po++)
	{
		if(ptrFeatureMatch->refIndex >= numRefFeatures3D)
		{
			refPolyData->Delete();
			trgPolyData->Delete();
			refPolyData = NULL;
			trgPolyData = NULL;
			return false;
		}
		if(ptrFeatureMatch->trgIndex >= numTrgFeatures3D)
		{
			refPolyData->Delete();
			trgPolyData->Delete();
			refPolyData = NULL;
			trgPolyData = NULL;
			return false;
		}
		*ptrPolyRef++ = refFeatures3D[ptrFeatureMatch->refIndex].x;
		*ptrPolyRef++ = refFeatures3D[ptrFeatureMatch->refIndex].y;
		*ptrPolyRef++ = refFeatures3D[ptrFeatureMatch->refIndex].z;
		*ptrPolyTrg++ = trgFeatures3D[ptrFeatureMatch->trgIndex].x;
		*ptrPolyTrg++ = trgFeatures3D[ptrFeatureMatch->trgIndex].y;
		*ptrPolyTrg++ = trgFeatures3D[ptrFeatureMatch->trgIndex].z;
		ptrFeatureMatch++;
	}
	return true;
}


vtkTransform*	ReLOC::FeaturePointsPairwiseRegistration::Register(vtkPolyDataCollection* polyDataCollection, const int cloudToRegisterIdx, int* nearCloudMatchIdx, int numNearCloudMatches, const std::vector<Feature3D*> &vFeatures3D, const std::vector<int> &vNFeatures3D, const std::vector<FeatureMatch*> &vFeatureMatches, const std::vector<int> &vNFeatureMatches, const std::vector<std::pair<int, int> > &vMatchPairs, const float threshDistance )
{
	vtkPolyData* nearCloudsFeaturePoints;
	vtkPolyData* cloudToRegisterFeaturePoints;
	if( !CreatePolyData(polyDataCollection, cloudToRegisterIdx, nearCloudMatchIdx, numNearCloudMatches, vFeatures3D, vNFeatures3D, vFeatureMatches, vNFeatureMatches, vMatchPairs, nearCloudsFeaturePoints, cloudToRegisterFeaturePoints, threshDistance) )
	{
		return NULL;
	}

	m_transform = Register(nearCloudsFeaturePoints, cloudToRegisterFeaturePoints);

	nearCloudsFeaturePoints->Delete();
	cloudToRegisterFeaturePoints->Delete();

	return m_transform;
}

bool ReLOC::FeaturePointsPairwiseRegistration::CreatePolyData(vtkPolyDataCollection* polyDataCollection, const int cloudToRegisterIdx, int* nearCloudMatchIdx, int numNearCloudMatches, const std::vector<Feature3D*> &vFeatures3D, const std::vector<int> &vNFeatures3D, const std::vector<FeatureMatch*> &vFeatureMatches, const std::vector<int> &vNFeatureMatches, const std::vector<std::pair<int, int> > &vMatchPairs, vtkPolyData* &nearCloudsFeaturePoints, vtkPolyData* &cloudToRegisterFeaturePoints, const float threshDistance)
{
	nearCloudsFeaturePoints = NULL;
	cloudToRegisterFeaturePoints = NULL;
	if(numNearCloudMatches<=0)
	{
		return false;
	}

	//calc number of clouds points
	int numCloudPoints=0;
	for(int ma=0; ma<numNearCloudMatches; ma++)
	{
		if( nearCloudMatchIdx[ma]>=0 )
		{
			numCloudPoints += vNFeatureMatches[ nearCloudMatchIdx[ma] ];
		}
	}

	float* newCloudToRegisterPoints = new float[numCloudPoints*3];
	float* newNearCloudsPoints = new float[numCloudPoints*3];

	int nearCloudIdx;
	bool isCloudToRegisterFirst;
	FeatureMatch* ptrFeatureMatch;

	float* ptrCloudToRegisterPoints = (float*)((vtkPolyData*)(polyDataCollection->GetItemAsObject(cloudToRegisterIdx)))->GetPoints()->GetVoidPointer(0);
	float* ptrNearCloudPoints;

	Feature3D* ptrCloudToRegisterFeature3D = vFeatures3D[cloudToRegisterIdx];
	int numCloudToRegisterFeatures3D = vNFeatures3D[cloudToRegisterIdx];
	int numNearCloudFeatures3D;
	Feature3D* ptrNearCloudFeature3D;

	int cloudToRegisterMatchId, nearCloudMatchId;
	int cloudToRegisterPointId, nearCloudPointId;

	//double point[3];
	float* ptrNewCloudToRegisterPoints = newCloudToRegisterPoints;
	float* ptrNewNearCloudsPoints = newNearCloudsPoints;
	int numNewPoints = 0;
	for(int ma=0; ma<numNearCloudMatches; ma++)
	{
		if( nearCloudMatchIdx[ma]>=0 )
		{
			//find id of near cloud
			if( vMatchPairs[ nearCloudMatchIdx[ma] ].first == cloudToRegisterIdx)
			{
				nearCloudIdx = vMatchPairs[ nearCloudMatchIdx[ma] ].second;
				isCloudToRegisterFirst = true;
			}
			else
			{
				nearCloudIdx = vMatchPairs[ nearCloudMatchIdx[ma] ].first;
				isCloudToRegisterFirst = false;
			}
			ptrNearCloudPoints = (float*)((vtkPolyData*)(polyDataCollection->GetItemAsObject(nearCloudIdx)))->GetPoints()->GetVoidPointer(0);

			ptrNearCloudFeature3D = vFeatures3D[nearCloudIdx];
			numNearCloudFeatures3D = vNFeatures3D[nearCloudIdx];

			ptrFeatureMatch = vFeatureMatches[ nearCloudMatchIdx[ma] ];

			char escKey = 'e';
			for(int po=0; po<vNFeatureMatches[ nearCloudMatchIdx[ma] ]; po++)
			{
				if(isCloudToRegisterFirst)
				{
					cloudToRegisterMatchId = ptrFeatureMatch->refIndex;
					nearCloudMatchId = ptrFeatureMatch->trgIndex;
				}
				else
				{
					cloudToRegisterMatchId = ptrFeatureMatch->trgIndex;
					nearCloudMatchId = ptrFeatureMatch->refIndex;
				}

				if(cloudToRegisterMatchId >= numCloudToRegisterFeatures3D)
				{
					printf("ERROR: cloudToRegisterMatchId >= numCloudToRegisterFeatures3D\n");
					return false;
				}
				if(nearCloudMatchId >= numNearCloudFeatures3D)
				{
					printf("ERROR: nearCloudMatchId >= numNearCloudFeatures3D\n");
					return false;
				}
				cloudToRegisterPointId = ptrCloudToRegisterFeature3D[cloudToRegisterMatchId].index;
				nearCloudPointId = ptrNearCloudFeature3D[nearCloudMatchId].index;

				if(EuclideanDistance_3D(ptrCloudToRegisterPoints + 3*cloudToRegisterPointId, ptrNearCloudPoints + 3*nearCloudPointId)<threshDistance)
				{

					float prova1[3];
					float prova2[3];
					memcpy(prova1, ptrCloudToRegisterPoints + 3*cloudToRegisterPointId, 3*sizeof(float));
					memcpy(prova2, ptrNearCloudPoints + 3*nearCloudPointId, 3*sizeof(float));

					memcpy(ptrNewCloudToRegisterPoints, ptrCloudToRegisterPoints + 3*cloudToRegisterPointId, 3*sizeof(float));
					memcpy(ptrNewNearCloudsPoints, ptrNearCloudPoints + 3*nearCloudPointId, 3*sizeof(float));

					//if(escKey == 'e')
					//{
					//	vtkMan.AddSphere(*ptrCloudToRegister, *(ptrCloudToRegister+1), *(ptrCloudToRegister+2), 1, 8, 8, 0.0, 1.0, 0.0, 0.5);
					//	vtkMan.AddSphere(*ptrNearClouds, *(ptrNearClouds+1), *(ptrNearClouds+2), 1, 8, 8, 1.0, 0.0, 1.0, 0.5);
					//	vtkMan.PrintNumber(*ptrCloudToRegister, 0, 0);
					//	vtkMan.PrintNumber(*(ptrCloudToRegister+1), 0, 1);
					//	vtkMan.PrintNumber(*(ptrCloudToRegister+2), 0, 2);
					//	vtkMan.PrintNumber(*ptrNearClouds, 1, 0);
					//	vtkMan.PrintNumber(*(ptrNearClouds+1), 1, 1);
					//	vtkMan.PrintNumber(*(ptrNearClouds+2), 1, 2);
					//	escKey = vtkMan.DrawAndWait();
					//}


					ptrNewCloudToRegisterPoints+=3;
					ptrNewNearCloudsPoints+=3;
					numNewPoints++;
				}
				ptrFeatureMatch++;
			}
			//vtkMan.RemoveAllPolyData();
		}

	}
	if(!numNewPoints)
	{
		return false;
	}

	nearCloudsFeaturePoints = ReLOC::CreatePolyData( newNearCloudsPoints, numNewPoints);
	cloudToRegisterFeaturePoints = ReLOC::CreatePolyData( newCloudToRegisterPoints, numNewPoints);
	delete[] newCloudToRegisterPoints;
	delete[] newNearCloudsPoints;

	return true;
}





ReLOC::Ransac::Ransac(const float ransacThreshInUnitOfMeshRes, const int ransacMaxIters, const double ransacProb)
{
	threshInUnitOfMeshRes = ransacThreshInUnitOfMeshRes;
	maxIters = ransacMaxIters;
	prob = ransacProb;

	ptr = NULL;
	aoptr = NULL;

	consensusSetSize = 0;

	verbose = false;

	vtkObject::GlobalWarningDisplayOff();	//to disable warnings produced by the JacobiN method 
}

void ReLOC::Ransac::Allocate(const double meanMeshRes)
{
	if(aoptr)
	{
		delete aoptr;
	}
	aoptr = new HornOrientation(false, false);

	if(ptr)
	{
		delete ptr;
	}
	ptr = new RANSACPairwiseRegistration(*aoptr, threshInUnitOfMeshRes * meanMeshRes, prob, maxIters, false, true);
}

vtkTransform* ReLOC::Ransac::Apply(const vector<FeatureMatch> &matches, const Feature3D* feat3D_trg, const int nFeat3D_trg, const Feature3D* feat3D_ref, const int nFeat3D_ref)
{
	consensusSetSize = 0;
	ptrMaxConsensusSet = NULL;
	outMatches.clear();

	outTransform = ptr->Register((FeatureMatch*)&matches[0], (int)matches.size(), (Feature3D*)feat3D_trg, nFeat3D_trg, (Feature3D*)feat3D_ref, nFeat3D_ref);
	if(outTransform)
	{
		consensusSetSize = ((RANSACPairwiseRegistration*)ptr)->GetMaxConsensusSetSize();
		outMatches.resize(consensusSetSize);

		ptrMaxConsensusSet = ((RANSACPairwiseRegistration*)ptr)->GetMaxConsensusSet();
		int* ptrPtrConsensusSet = ptrMaxConsensusSet;
		for(int co=0; co<consensusSetSize; co++)
		{
			outMatches[co] = matches[*ptrPtrConsensusSet];
			ptrPtrConsensusSet++;
		}
	}
	else
	{
		if(verbose)
		{
			cout << "Impossible to find rotoTranslation" << endl;
		}
	}

	return outTransform;
}

ReLOC::Ransac::~Ransac()
{
	if(ptr)
	{
		delete ptr;
	}
	if(aoptr)
	{
		delete aoptr;
	}
}




ReLOC::AbsoluteOrientation::AbsoluteOrientation(bool calcScale, bool calcRMSE)
	:m_calcRMSE(calcRMSE)
	,m_calcScale(calcScale)
{
	m_transform = vtkTransform::New();
	m_rmse = -1.0;
}

vtkTransform* ReLOC::AbsoluteOrientation::Orientate(vtkPolyData* polyData1In, vtkPolyData* polyData2In )
{
	if(polyData1In->GetNumberOfPoints() != polyData2In->GetNumberOfPoints())
	{
		if(m_transform)
		{
			m_transform->Delete();
			m_transform = NULL;
		}
		return NULL;
	}

	double* ptrMatr = NULL;
	int ret = Orientate( (float*)(polyData1In->GetPoints()->GetVoidPointer(0)), (float*)(polyData2In->GetPoints()->GetVoidPointer(0)), polyData1In->GetNumberOfPoints(), ptrMatr);

	if(m_transform)
	{
		m_transform->Delete();
		m_transform = NULL;
	}
	if(ret>=0)
	{
		m_transform = vtkTransform::New();
		m_transform->SetMatrix(ptrMatr);
	}

	return m_transform;
}

int ReLOC::AbsoluteOrientation::ExtractPoints(vtkPolyData* poly, Feature3D* features3D, int numFeatures3D, float* &outPoints )
{
	outPoints = new float[numFeatures3D*3];

	float* points = (float*)poly->GetPoints()->GetVoidPointer(0);
	for(int fe = 0; fe < numFeatures3D; fe++)
	{
		memcpy(outPoints + 3*fe, points + 3*features3D[fe].index, 3*sizeof(float));
	}

	return numFeatures3D;
}

vtkPolyData* ReLOC::AbsoluteOrientation::ExtractPoints(vtkPolyData* poly, Feature3D* features3D, int numFeatures3D )
{
	float* outPoints;
	int numPoints = ExtractPoints(poly, features3D, numFeatures3D, outPoints );
	vtkPolyData* outPoly = CreatePolyData(outPoints, numPoints);
	delete[] outPoints;
	return outPoly;
}

bool ReLOC::AbsoluteOrientation::ExtractPoints(FeatureMatch* featureMatches, int numFeatureMatches, Feature3D* trgFeatures3D, int numTrgFeatures3D, Feature3D* refFeatures3D, int numRefFeatures3D, float* &trgPoints, float* &refPoints)
{
	refPoints = new float[numFeatureMatches*3];
	trgPoints = new float[numFeatureMatches*3];

	float* ptrRefPoints = refPoints;
	float* ptrTrgPoints = trgPoints;
	FeatureMatch* ptrFeatureMatch = &featureMatches[0];
	for(int po=0; po<numFeatureMatches; po++)
	{
		if(ptrFeatureMatch->refIndex >= numRefFeatures3D)
		{
			return false;
		}
		if(ptrFeatureMatch->trgIndex >= numTrgFeatures3D)
		{
			return false;
		}
		*ptrRefPoints++ = refFeatures3D[ptrFeatureMatch->refIndex].x;
		*ptrRefPoints++ = refFeatures3D[ptrFeatureMatch->refIndex].y;
		*ptrRefPoints++ = refFeatures3D[ptrFeatureMatch->refIndex].z;
		*ptrTrgPoints++ = trgFeatures3D[ptrFeatureMatch->trgIndex].x;
		*ptrTrgPoints++ = trgFeatures3D[ptrFeatureMatch->trgIndex].y;
		*ptrTrgPoints++ = trgFeatures3D[ptrFeatureMatch->trgIndex].z;
		ptrFeatureMatch++;
	}
	return true;
}


int ReLOC::AbsoluteOrientation::Orientate(FeatureMatch* featureMatches, int numFeatureMatches, Feature3D* trgFeatures3D, int numTrgFeatures3D, Feature3D* refFeatures3D, int numRefFeatures3D, double* &matrTransform)
{
	float* refPoints;
	float* trgPoints;
	if( !ExtractPoints(featureMatches, numFeatureMatches, trgFeatures3D, numTrgFeatures3D, refFeatures3D, numRefFeatures3D, trgPoints, refPoints) )
	{
		delete[] trgPoints;
		delete[] refPoints;
		return -1;
	}

	int result = Orientate(trgPoints, refPoints, numFeatureMatches, matrTransform);

	delete[] trgPoints;
	delete[] refPoints;

	return result;
}

vtkTransform* ReLOC::AbsoluteOrientation::Orientate(FeatureMatch* featureMatches, int numFeatureMatches, Feature3D* trgFeatures3D, int numTrgFeatures3D, Feature3D* refFeatures3D, int numRefFeatures3D)
{
	double* matrTransform;
	if(Orientate(featureMatches, numFeatureMatches, trgFeatures3D, numTrgFeatures3D, refFeatures3D, numRefFeatures3D, matrTransform) == -1)
	{
		return NULL;
	}

	if(m_transform)
	{
		m_transform->Delete();
		m_transform = NULL;
	}
	m_transform = vtkTransform::New();
	m_transform->SetMatrix(matrTransform);

	return m_transform;
}

double ReLOC::AbsoluteOrientation::CalcRMSE(float* points1, float* points2, int nPoints)
{
	double rmse = 0;
	float* ptrPoints1 = points1;
	float* ptrPoints2 = points2;
	for(int po=0; po<nPoints; po++)
	{
		rmse += vtkMath::Distance2BetweenPoints(ptrPoints1, ptrPoints2);
		ptrPoints1+=3;
		ptrPoints2+=3;
	}
	rmse /=(double)nPoints;
	rmse = sqrt(rmse);

	return rmse;
}

double ReLOC::AbsoluteOrientation::CalcRMSE(vtkPolyData* polyData1, vtkPolyData* polyData2)
{
	if(polyData1->GetNumberOfPoints() != polyData2->GetNumberOfPoints())
	{
		return -1.0;
	}

	return CalcRMSE( (float*)polyData1->GetPoints()->GetVoidPointer(0), (float*)polyData2->GetPoints()->GetVoidPointer(0), polyData1->GetNumberOfPoints() );
}


ReLOC::AbsoluteOrientation::~AbsoluteOrientation()
{
	if(m_transform)
	{
		m_transform->Delete();
		m_transform = NULL;
	}
}





ReLOC::HornOrientation::HornOrientation(bool calcScale, bool calcRMSE)
	:AbsoluteOrientation(calcScale, calcRMSE)
{
}

int ReLOC::HornOrientation::Orientate(float *points1, float* points2, int nPoints, double* &matrTransform)
{
	if(nPoints < 3)
	{
		return -1;
	}

	float* points1Out = new float[nPoints*3];
	float* points2Out = new float[nPoints*3];
	memcpy(points1Out, points1, sizeof(float)*nPoints*3);
	memcpy(points2Out, points2, sizeof(float)*nPoints*3);

	float barycenter1[3];
	float barycenter2[3];
	GetCentroid(points1Out, nPoints, barycenter1);
	GetCentroid(points2Out, nPoints, barycenter2);

	//translate points in barycenters
	float invbarycenter1[3];
	float invbarycenter2[3];
	invbarycenter1[0] = -barycenter1[0];
	invbarycenter1[1] = -barycenter1[1];
	invbarycenter1[2] = -barycenter1[2];
	invbarycenter2[0] = -barycenter2[0];
	invbarycenter2[1] = -barycenter2[1];
	invbarycenter2[2] = -barycenter2[2];
	Translate(points1Out, nPoints, invbarycenter1);
	Translate(points2Out, nPoints, invbarycenter2);

	float* ptr1;
	float* ptr2;

	//calc scale
	float scale=1.0;
	if(m_calcScale)
	{
		double sum1 = 0.0;
		double sum2 = 0.0;
		ptr1 = points1Out;
		ptr2 = points2Out;
		for(int po=0; po<nPoints; po++)
		{
			sum1 += sqrt( (ptr1[0] * ptr1[0]) + (ptr1[1] * ptr1[1]) + (ptr1[2] * ptr1[2]) );
			sum2 += sqrt( (ptr2[0] * ptr2[0]) + (ptr2[1] * ptr2[1]) + (ptr2[2] * ptr2[2]) );
			ptr1+=3;
			ptr2+=3;
		}
		if(!IsZero( sum2 ))
		{
			scale = sum1 / sum2;
		}
	}

	//calc rotation
	float Sxx=0, Sxy=0, Sxz=0;
	float Syx=0, Syy=0, Syz=0;
	float Szx=0, Szy=0, Szz=0;
	ptr1 = points1Out;
	ptr2 = points2Out;
	for (int po=0; po<nPoints ;po++)
	{
		Sxx += ptr1[0] * ptr2[0];
		Sxy += ptr1[0] * ptr2[1];
		Sxz += ptr1[0] * ptr2[2];
		Syx += ptr1[1] * ptr2[0];
		Syy += ptr1[1] * ptr2[1];
		Syz += ptr1[1] * ptr2[2];
		Szx += ptr1[2] * ptr2[0];
		Szy += ptr1[2] * ptr2[1];
		Szz += ptr1[2] * ptr2[2];
		ptr1 +=3;
		ptr2 +=3;
	}


	//calc N
	double *n[4];
	n[0] = new double[4];
	n[1] = new double[4];
	n[2] = new double[4];
	n[3] = new double[4];
	n[0][0] = Sxx + Syy + Szz;
	n[0][1] = Syz - Szy;
	n[0][2] = Szx - Sxz;
	n[0][3] = Sxy - Syx;

	n[1][0] = Syz - Szy;
	n[1][1] = Sxx - Syy - Szz;
	n[1][2] = Sxy + Syx;
	n[1][3] = Szx + Sxz;

	n[2][0] = Szx - Sxz;
	n[2][1] = Sxy + Syx;
	n[2][2] = - Sxx + Syy - Szz;
	n[2][3] = Syz + Szy;

	n[3][0] = Sxy - Syx;
	n[3][1] = Szx + Sxz;
	n[3][2] = Syz + Szy;
	n[3][3] = - Sxx - Syy + Szz;

	double *evect[4];
	evect[0] = new double[4];
	evect[1] = new double[4];
	evect[2] = new double[4];
	evect[3] = new double[4];
	memset(evect[0], 0.0, sizeof(double)*4);
	memset(evect[1], 0.0, sizeof(double)*4);
	memset(evect[2], 0.0, sizeof(double)*4);
	memset(evect[3], 0.0, sizeof(double)*4);

	double evalues[] = {0, 0, 0, 0};

	int ret = vtkMath::JacobiN(n, 4, evalues, evect);
	if(ret == 0)
	{
		delete[] points1Out;
		delete[] points2Out;

		for (int i = 0; i < 4; i++)
		{
			delete [] evect[i];
			delete [] n[i];
		}

		return -1;
	}



	//calc final transform
	//initialize last row
	for(int co=0; co<3; co++)
	{
		m_matrTransform[3*4 + co] = 0;
	}
	m_matrTransform[3*4 + 3] = 1;

	//translate to barycenter1
	for(int ro=0; ro<3; ro++)
	{
		m_matrTransform[ro*4+3] = barycenter1[ro];
	}

	//first eigenvector is the rotation quaternion
	//scale and rotate
	double rotMatr[3][3];


	//as  JacobiN returns evect in column wise:
	double maxEvect[] = { - evect[0][0], evect[1][0], evect[2][0], evect[3][0]};
	vtkMath::QuaternionToMatrix3x3(maxEvect, rotMatr);
	for(int ro=0; ro<3; ro++)
	{
		for(int co=0; co<3; co++)
		{
			m_matrTransform[ro*4+co] = rotMatr[ro][co] * scale;
		}
	}

	//translate to invbarycenter2 (It's like m_transform->Translate(invbarycenter2) )
	double matrixTrasl[16];
	vtkMatrix4x4::Identity(matrixTrasl);
	matrixTrasl[4*0+3] = invbarycenter2[0];
	matrixTrasl[4*1+3] = invbarycenter2[1];
	matrixTrasl[4*2+3] = invbarycenter2[2];
	vtkMatrix4x4::Multiply4x4(m_matrTransform, matrixTrasl, m_matrTransform);


	if(m_calcRMSE)
	{
		Translate(points2Out, nPoints, barycenter2);
		Translate(points1Out, nPoints, barycenter1);
		ApplyTransform(points2Out, nPoints, m_matrTransform);
		m_rmse = CalcRMSE(points1Out, points2Out, nPoints);
	}

	delete[] points1Out;
	delete[] points2Out;

	for (int i = 0; i < 4; i++)
	{
		delete [] evect[i];
		delete [] n[i];
	}

	matrTransform = m_matrTransform;
	return 0;

}



//*************  RANSAC *********

ReLOC::RANSACPairwiseRegistration::RANSACPairwiseRegistration(AbsoluteOrientation &absOrientation, float outlierDistance, float probGoodResult, int nMaxIter, bool calcRMSE, bool verbose)
:FeaturePointsPairwiseRegistration(absOrientation.GetCalcScale(), calcRMSE, verbose),
m_absOrientation(absOrientation),
m_outlierDistance(outlierDistance),
m_maxConsensusSetSize(0),
m_maxConsensusSet(NULL),
m_probGoodResult(probGoodResult),
m_nMaxIter(nMaxIter)
{
}

ReLOC::RANSACPairwiseRegistration::~RANSACPairwiseRegistration()
{
	if(m_maxConsensusSet)
	{
		delete[] m_maxConsensusSet;
		m_maxConsensusSet = NULL;
	}
}

vtkTransform*	ReLOC::RANSACPairwiseRegistration::Register(vtkPolyData* polyData1In, vtkPolyData* polyData2In)
{

	
	m_maxConsensusSetSize = 0;
	m_rmse = -1;

	if(polyData1In->GetNumberOfPoints() != polyData2In->GetNumberOfPoints())
	{
		if(m_transform)
		{
			m_transform->Delete();
			m_transform = NULL;
			m_rmse = -1;
		}
		return NULL;
	}

	int nPoints = polyData1In->GetNumberOfPoints();
	if(nPoints < 3)
	{
		if(m_transform)
		{
			m_transform->Delete();
			m_transform = NULL;
			m_rmse = -1;
		}
		return NULL;
	}

	if(m_transform)
	{
		m_transform->Delete();
		m_transform = NULL;
	}
	m_transform = vtkTransform::New();
	m_transform->Identity();


	float* ptrPolyData1In = (float*)(polyData1In->GetPoints()->GetVoidPointer(0));
	float* ptrPolyData2In = (float*)(polyData2In->GetPoints()->GetVoidPointer(0));

	float* points1Temp = new float[nPoints*3];
	float* points2Temp = new float[nPoints*3];
	memcpy(points1Temp, ptrPolyData1In, sizeof(float)*nPoints*3);
	memcpy(points2Temp, ptrPolyData2In, sizeof(float)*nPoints*3);

	int consensusSetSize;
	int maxConsensusSetSize = 0;
	int* consensusSet = new int[nPoints];

	if(m_maxConsensusSet)
	{
		delete[] m_maxConsensusSet;
		m_maxConsensusSet = NULL;
	}
	m_maxConsensusSet = new int[nPoints];


	double inliersPercent = 0.001;

	// Calc attempts number
	long long N = m_nMaxIter; //CVPipes::Core::Impl::Min(m_nMaxIter, ceil(  log(1.0 - m_probGoodResult) / log(1.0 - pow(inliersPercent, 3))  ) );

	int rnd[3];
	srand(5);
	float randPoints1[9];	//3 polyData1In random points
	float randPoints2[9];	//3 polyData2In random points

	double* ptrTempMatr;
	int nTestAtMaxConsensus = 0;
	int test = 0;
	for (test = 0; test < N; test++)
	{
		//insert points in temp array
		memcpy(points1Temp, ptrPolyData1In, sizeof(float)*3*nPoints);
		memcpy(points2Temp, ptrPolyData2In, sizeof(float)*3*nPoints);

		//select 3 random points
		rnd[0] = (int)( ((float)rand()/(float)RAND_MAX) * nPoints);
		do
		{
			rnd[1] = (int)(  ((float)rand()/(float)RAND_MAX) * nPoints);
		}while (rnd[0] == rnd[1]);

		do
		{
			rnd[2] = (int)(  ((float)rand()/(float)RAND_MAX) * nPoints);
		} while ((rnd[2] == rnd[0])||(rnd[2] == rnd[1]));

		for (int i = 0; i < 3; i++)
		{
			memcpy(randPoints1+i*3, ptrPolyData1In + 3*rnd[i], sizeof(float)*3) ;
			memcpy(randPoints2+i*3, ptrPolyData2In + 3*rnd[i], sizeof(float)*3) ;
		}

		int ret = m_absOrientation.Orientate(randPoints1, randPoints2, 3, ptrTempMatr);
		if(ret < 0)
		{
			continue;
		}

		ApplyTransform(points2Temp, nPoints, ptrTempMatr);

		consensusSetSize = 0;
		float euclideanDistance;
		for (int po=0; po<nPoints; po++)
		{
			euclideanDistance = EuclideanDistance_3D(points1Temp + po*3, points2Temp + po*3);
			if (euclideanDistance<m_outlierDistance )
			{
				consensusSet[consensusSetSize] = po;
				consensusSetSize++;
			}
		}

		if (consensusSetSize > maxConsensusSetSize)
		{
			//new best consensus set
			maxConsensusSetSize = consensusSetSize;
			memcpy(m_maxConsensusSet, consensusSet, maxConsensusSetSize*sizeof(int));

			// adaptive update of iteractions number
			if(maxConsensusSetSize == nPoints)
			{
				N = 0;
			}
			else
			{
				N = Min(m_nMaxIter, (int)ceil(  log(1.0 - m_probGoodResult) / log(1.0 - pow(maxConsensusSetSize / (double)nPoints, 3))  ) );
			}

			nTestAtMaxConsensus = test + 1;
		}
	}

	m_maxConsensusSetSize = maxConsensusSetSize;
	if(m_maxConsensusSetSize < 3)
	{
		m_transform->Delete();
		m_transform = NULL;
		m_rmse = -1;
	}
	else
	{
		for (int i=0; i<maxConsensusSetSize; i++)
		{
			memcpy(points1Temp+i*3, ptrPolyData1In + 3*m_maxConsensusSet[i], sizeof(float)*3) ;
			memcpy(points2Temp+i*3, ptrPolyData2In + 3*m_maxConsensusSet[i], sizeof(float)*3) ;
		}

		//get final transform
		int ret = m_absOrientation.Orientate(points1Temp, points2Temp, maxConsensusSetSize, ptrTempMatr);
		if(ret < 0)
		{
			m_transform->Delete();
			m_transform = NULL;
			m_rmse = -1;
		}
		
		m_transform->SetMatrix(ptrTempMatr);

		if(m_calcRMSE)
		{
			ApplyTransform(points2Temp, maxConsensusSetSize, ptrTempMatr);
			m_rmse = CalcRMSE(points1Temp, points2Temp, maxConsensusSetSize);
		}
	}


	delete[] points1Temp;
	delete[] points2Temp;
	delete[] consensusSet;
	return m_transform;
}

