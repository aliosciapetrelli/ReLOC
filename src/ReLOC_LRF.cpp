#include "ReLOC_LRF.h"

#include <algorithm>


using namespace std;

#include "vtkMath.h"
#include "vtkPointData.h"
#include <set>
#include "vtkIdList.h"


ReLOC::Feature3D::Feature3D()
: x(0)
, y(0)
, z(0)
, scale(0)
, index(-1)
, score(0)
{
	memset(rf, 0, 9*sizeof(float));
}

ReLOC::Feature3D::Feature3D(const Feature3D & other)
: x(other.x)
, y(other.y)
, z(other.z)
, scale(other.scale)
, index(other.index)
, score(other.score)
{
	memcpy(rf, other.rf, 9*sizeof(float));
}

ReLOC::Feature3D& ReLOC::Feature3D::operator =(const Feature3D & other)
{
  x = other.x;
  y = other.y;
  z = other.z;
  scale = other.scale;
  index = other.index;
  score = other.score;
  memcpy(rf, other.rf, 9*sizeof(float));
	return *this;
}

bool ReLOC::Feature3D::operator !=(const Feature3D & other) const
{
	return ! (*this == other);
}

bool ReLOC::Feature3D::operator ==(const Feature3D & other) const
{
	if (x != other.x)
		return false;

	if (y != other.y)
		return false;

	if (z != other.z)
		return false;

	if (scale != other.scale)
		return false;

	if (index != other.index)
		return false;

	if (score != other.score)
		return false;

	for (int i=0; i<9; i++)
		if (rf[i] != other.rf[i])
			return false;

	return true;
}




ReLOC::I3DDetector::I3DDetector()
{
	m_feat = NULL;
	m_numFeat = 0;
	m_featCapacity = 0;
	m_borderDistance = -1;
}

ReLOC::I3DDetector::~I3DDetector()
{
	if (m_feat != NULL)
	{
		delete [] m_feat;
		m_feat = NULL;
	}

	if(m_vFeat.size() > 0)
	{
		for(unsigned int i = 0; i < m_vFeat.size(); i++)
		{
			if(m_vFeat[i] != NULL)
			{
				delete [] m_vFeat[i];
				m_vFeat[i] = NULL;
			}
		}
	}
}

void ReLOC::I3DDetector::updateFeatCapacity(Feature3D * & feat, int & capacity, const int nFeat)
{
	if (capacity < nFeat)
	{
		if (feat != NULL)
		{
			delete [] feat;
			feat = NULL;
		}
	}

	if(feat == NULL && nFeat > 0) {
		feat = new Feature3D[nFeat];
		capacity = nFeat;
	}
}

void ReLOC::I3DDetector::setFeatures(Feature3D* feat, const int numFeatures, const bool copyFeats)
{
	if(copyFeats)
	{
		updateFeatCapacity(m_feat, m_featCapacity, numFeatures);
		memcpy(m_feat, feat, sizeof(Feature3D)*numFeatures);
	}
	else
	{
		if (m_feat != NULL)
		{
			delete [] m_feat;
			m_feat = NULL;
		}
		m_feat = feat;
		m_featCapacity = numFeatures;
	}
	m_numFeat = numFeatures;
};

void ReLOC::I3DDetector::resizeVect(const unsigned int newLength)
{
	if(m_vFeat.size() < newLength)
	{
		for(unsigned int i = 0; i < m_vFeat.size(); i++)
		{
			if(m_vFeat[i] != NULL)
			{
				delete [] m_vFeat[i];
				m_vFeat[i] = NULL;
			}
		}
	}

	if(newLength < m_vFeat.size())
	{
		for(unsigned int i = newLength; i < m_vFeat.size(); i++)
		{
			if(m_vFeat[i] != NULL)
			{
				delete [] m_vFeat[i];
				m_vFeat[i] = NULL;
			}
		}
	}

	m_vFeatCapacity.resize(newLength);
	m_vFeat.resize(newLength);
	m_vNumFeat.resize(newLength);

	for(unsigned int i = 0; i < m_vFeat.size(); i++)
	{
		m_vFeat[i] = NULL;
	}
}


int ReLOC::I3DDetector::extract(vtkPolyData* cloud, Feature3D* & feat)
{
	int nFeat;
	extractImpl(cloud, m_feat, nFeat, m_featCapacity);

	feat = m_feat;
	m_numFeat = nFeat;

	return nFeat;
}

void ReLOC::I3DDetector::extract(vtkPolyDataCollection* polyDataCollection )
{
 	resizeVect(polyDataCollection->GetNumberOfItems());
	polyDataCollection->InitTraversal();
	std::vector<int>::iterator itN = m_vNumFeat.begin();
	std::vector<int>::iterator itFeatCapacity = m_vFeatCapacity.begin();
	int count = 1;
	printf("\n");
	for(std::vector<Feature3D*>::iterator itFeat = m_vFeat.begin(); itFeat != m_vFeat.end(); ++itFeat )
	{
		vtkPolyData* poly = polyDataCollection->GetNextItem();
		extractImpl(poly, *itFeat, *itN, *itFeatCapacity);
		printf("\tmodel %d / %d : %d\r", count, polyDataCollection->GetNumberOfItems(), *itN);
		itN++;
		itFeatCapacity++;
		count++;
	}
}

bool ReLOC::I3DDetector::save(const char *path, const int numFeat, const Feature3D* feat, const bool binaryFile)
{
	CreateFullPath(GetPathDir(path));

	if(binaryFile)
	{
		//binary version
		fstream binary_file(path,ios::out|ios::binary);
		if (!binary_file.is_open() )
			return false;

		binary_file.write((char*)(&numFeat), sizeof(numFeat));

		for(int i=0; i<numFeat; i++)
		{
			binary_file.write((char*)(&feat[i].x), sizeof(feat[i].x));
			binary_file.write((char*)(&feat[i].y), sizeof(feat[i].y));
			binary_file.write((char*)(&feat[i].z), sizeof(feat[i].z));
			binary_file.write((char*)(&feat[i].scale), sizeof(feat[i].scale));
			binary_file.write((char*)(&feat[i].index), sizeof(feat[i].index));
			binary_file.write((char*)(&feat[i].score), sizeof(feat[i].score));
			for(int j=0; j<9; j++)
				binary_file.write((char*)(&feat[i].rf[j]), sizeof(feat[i].rf[j]));
		}
		binary_file.close();
	}
	else
	{
		//ascii version
		std::ofstream file(path);
		file.precision(15);

		if (!file.is_open() )
			return false;

		file << numFeat << "\n";

		for(int i=0; i<numFeat; i++)
		{
			file << feat[i].x << " " << feat[i].y << " " << feat[i].z << " " << feat[i].scale << " " << feat[i].index << " " << feat[i].score << " ";
			for(int j=0; j<9; j++)
				file << feat[i].rf[j] << " ";
			file << "\n";
		}

		file.close();
	}

	return true;
}

bool ReLOC::I3DDetector::save(const char *path, const bool binaryFile)
{
	return save(path, m_numFeat, m_feat, binaryFile);
}

bool ReLOC::I3DDetector::save(const char *path, const char *fileTemplate, const bool binaryFile)
{
	char filename[300];
	for(size_t de=0; de<m_vNumFeat.size(); de++)
	{
		sprintf(filename, "%s/%s%03d.feat", path, fileTemplate, de);
		if( !save(filename, m_vNumFeat[de], m_vFeat[de], binaryFile))
		{
			return false;
		}
	}
	return true;
}

int ReLOC::I3DDetector::load(const char *path, Feature3D* &feat, const bool binaryFile){

	load(path, m_feat, m_numFeat, m_featCapacity, binaryFile );

	feat = m_feat;

	return m_numFeat;
}

int ReLOC::I3DDetector::load(const char *path, Feature3D* &feat, int &numFeat, int &featCapacity, const bool binaryFile )
{
	if(binaryFile)
	{
		//binary version
		fstream binary_file(path, ios::binary|ios::in);

		if (!binary_file.is_open() )
		{
			numFeat = 0;
			updateFeatCapacity(feat, featCapacity, numFeat);
			return -1;
		}

		binary_file.read((char*)(&numFeat), sizeof(numFeat));

		updateFeatCapacity(feat, featCapacity, numFeat);

		for(int i=0; i<numFeat; i++)
		{
			binary_file.read((char*)(&feat[i].x), sizeof(feat[i].x));
			binary_file.read((char*)(&feat[i].y), sizeof(feat[i].y));
			binary_file.read((char*)(&feat[i].z), sizeof(feat[i].z));
			binary_file.read((char*)(&feat[i].scale), sizeof(feat[i].scale));
			binary_file.read((char*)(&feat[i].index), sizeof(feat[i].index));
			binary_file.read((char*)(&feat[i].score), sizeof(feat[i].score));
			for(int j=0; j<9; j++)
				binary_file.read((char*)(&feat[i].rf[j]), sizeof(feat[i].rf[j]));
		}

		binary_file.close();
	}
	else
	{
		//ascii version
		std::ifstream file(path);

		if (!file.is_open() )
		{
			numFeat = 0;
			updateFeatCapacity(feat, featCapacity, numFeat);
			return -1;
		}

		file >> numFeat;

		updateFeatCapacity(feat, featCapacity, numFeat);

		for(int i=0; i<numFeat; i++){
			file >> feat[i].x >> feat[i].y >> feat[i].z >> feat[i].scale >> feat[i].index >> feat[i].score;
			for(int j=0; j<9; j++)
				file >> feat[i].rf[j];
		}

		file.close();
	}

	return numFeat;
}

void ReLOC::I3DDetector::load(const char *path, const char *fileTemplate, const bool binaryFile)
{
	std::vector<std::string> foundFiles;
	FindFilesStartWith(path, fileTemplate, foundFiles);
	resizeVect((const unsigned int)foundFiles.size());
	for(size_t de=0; de<foundFiles.size(); de++)
	{
		load(foundFiles[de].c_str(), m_vFeat[de], m_vNumFeat[de], m_vFeatCapacity[de], binaryFile);
	}
}


int ReLOC::I3DDetector::compute3DIndex(vtkPolyData* poly, KdTree & kdt, Feature3D & feat, float meshResolution)
{
	int res;
	float point[3];
	point[0]=feat.x;
	point[1]=feat.y;
	point[2]=feat.z;
	double distance = kdt.FindNearestPoint(point,res);
	if(meshResolution>0)
	{
		if(distance>2*meshResolution)
			cout << "warning: get3DIndex -> distance > 2 * meshResolution";
	}
	return res;
}



bool ReLOC::I3DDetector::save(vector<string> &vAbsFileNames, const bool binaryFile)
{
	for(size_t de=0; de<m_vNumFeat.size(); de++)
	{
		if( !save(ChangeExtension(vAbsFileNames[de].c_str(), "feat").c_str(), m_vNumFeat[de], m_vFeat[de], binaryFile))
		{
			return false;
		}
	}
	return true;
}

bool ReLOC::I3DDetector::save(const char *path, vector<string> &vFileNames, const bool binaryFile)
{
	char filename[300];
	for(size_t de=0; de<m_vNumFeat.size(); de++)
	{
		sprintf(filename, "%s/%s", path, ChangeExtension(vFileNames[de].c_str(), "feat").c_str());
		if( !save(filename, m_vNumFeat[de], m_vFeat[de], binaryFile))
		{
			return false;
		}
	}
	return true;
}

void ReLOC::I3DDetector::load(vector<string> &vAbsFileNames, const bool binaryFile)
{
	resizeVect((const unsigned int)vAbsFileNames.size());
	for(size_t de=0; de<vAbsFileNames.size(); de++)
	{
		load(ChangeExtension(vAbsFileNames[de].c_str(), "feat").c_str(), m_vFeat[de], m_vNumFeat[de], m_vFeatCapacity[de], binaryFile);
	}
}

void ReLOC::I3DDetector::load(const char *path, vector<string> &vFileNames, const bool binaryFile)
{
	char filename[300];
	resizeVect((const unsigned int)vFileNames.size());
	for(size_t de=0; de<vFileNames.size(); de++)
	{
		sprintf(filename, "%s/%s", path, ChangeExtension(vFileNames[de].c_str(), "feat").c_str());
		load(filename, m_vFeat[de], m_vNumFeat[de], m_vFeatCapacity[de], binaryFile);
	}
}



/// **************************************** FlatPoints3DDetector ****************************************

ReLOC::FlatPoints3DDetector::FlatPoints3DDetector(bool random, double radius, double partitionRadius, double excludedPointsRadius, double partitionRadius_Hierarchical, float coverageAreaPercentage, float coverageAreaPercentage_Hierarchical, int minNeigh, int64_t seed) : I3DDetector()
{
	m_minNeigh = minNeigh;
	m_radius = radius;
	m_partitionRadius = partitionRadius;
	m_excludedPointsRadius = excludedPointsRadius;
	m_partitionRadius_Hierarchical = partitionRadius_Hierarchical;
	m_coverageAreaPercentage = coverageAreaPercentage;
	m_coverageAreaPercentage_Hierarchical = coverageAreaPercentage_Hierarchical;

	m_seed = seed;
	if(random)
		srand(time(NULL));
	else
		srand(seed);
}

ReLOC::FlatPoints3DDetector::FlatPoints3DDetector(FlatPoints3DDetectorParams &params)
{
	m_minNeigh = params.numMinNeighbors;
	m_radius = params.radius;
	m_partitionRadius = params.partitionRadius;
	m_excludedPointsRadius = params.excludedPointsRadius;
	m_partitionRadius_Hierarchical = params.partitionRadius_Hierarchical;
	m_coverageAreaPercentage = params.coverageAreaPercentage;
	m_coverageAreaPercentage_Hierarchical = params.coverageAreaPercentage_Hierarchical;
		
	m_seed = params.seed;
	if(params.random)
		srand(time(NULL));
	else
		srand(params.seed);
}

ReLOC::FlatPoints3DDetector::~FlatPoints3DDetector()
{
}

void ReLOC::FlatPoints3DDetector::extractImpl(vtkPolyData* cloud, Feature3D* & feat, int &numFeatures, int & featCapacity)
{
	if(m_radius <= 0.0 || m_minNeigh <= 0)
	{
		cout << "ERROR (FlatPoints3DDetector::extractImpl(): (m_radius <= 0.0 || m_minNeigh <= 0)";
		exit(-1);
	}


	int n_cloud_points = cloud->GetNumberOfPoints();

	std::vector<float> vFeaturePointScores(n_cloud_points, -2.0f);

	//first scan of the 2-step process
	std::vector<vtkIdType> vFeaturePointIndexes;
	ApplyFlatPointsDetector(cloud, vFeaturePointScores, m_excludedPointsRadius, m_partitionRadius, m_coverageAreaPercentage, vFeaturePointIndexes);

	//create polydata comprising the flat points extrated at step 1.
	vtkPolyData* mesh_step2 = CreatePolyData(cloud, &vFeaturePointIndexes[0], (vtkIdType)vFeaturePointIndexes.size(), false, false, false);

	//Get the scores of points of reduced mesh mesh_step2 from the scores of the points of the original cloud
	std::vector<float> vFeaturePointScores_H(vFeaturePointIndexes.size(), -2.0f);
	for(size_t i = 0; i < vFeaturePointIndexes.size(); i++)
	{
		vFeaturePointScores_H[i] = vFeaturePointScores[vFeaturePointIndexes[i]];
	}

	//second scan of the 2-step process
	std::vector<vtkIdType> vFeaturePointIndexes_H;
	ApplyFlatPointsDetector(mesh_step2, vFeaturePointScores_H, m_excludedPointsRadius, m_partitionRadius_Hierarchical, m_coverageAreaPercentage_Hierarchical, vFeaturePointIndexes_H);

	//allocate array of feature points
	numFeatures = (int)vFeaturePointIndexes_H.size();
	updateFeatCapacity(feat, featCapacity, numFeatures);

	//Get the indices of the points of the original mesh that have been selected as feature points from the the reduced mesh mesh_step2
	for(size_t i = 0; i < vFeaturePointIndexes_H.size(); i++)
	{
		vFeaturePointIndexes_H[i] = vFeaturePointIndexes[vFeaturePointIndexes_H[i]];
	}

	//populate the array of feature points
	for(size_t fe = 0; fe < vFeaturePointIndexes_H.size(); fe++)
	{			
		float* points = GetPolyDataPointsPointer(cloud, vFeaturePointIndexes_H[fe]);
		feat[fe].x = points[0];
		feat[fe].y = points[1];
		feat[fe].z = points[2];
		feat[fe].scale = m_radius;
		feat[fe].index = vFeaturePointIndexes_H[fe];
		feat[fe].score = vFeaturePointScores[ vFeaturePointIndexes_H[fe] ];			
	}


	//delete mesh
	mesh_step2->Delete();
}



/// *****************************************************************************************

void ReLOC::FlatPoints3DDetector::ApplyFlatPointsDetector(vtkPolyData* cloud, std::vector<float> & vFeaturePointScores, float exPRadius, float partRadius, float covAreaP, std::vector<vtkIdType> & vFeaturePointIndexes)
{
	float* neighNormal;
	float cosine= 0.0f;
	bool rightLeft = true;

	int n_cloud_points = cloud->GetNumberOfPoints();

	std::vector<bool> selectablePointSeed;
	std::vector<bool> selectablePointBest;
	selectablePointSeed.resize(n_cloud_points, true);
	selectablePointBest.resize(n_cloud_points, true);

	vFeaturePointIndexes.resize(n_cloud_points);
	int numFeatures = 0;  


	//compute Kdtree
	KdTree kdtree(cloud);

	vtkIdList* NNpoints_partRadius = vtkIdList::New();

	//get polydata point Normals
	float* ptrNormals = GetPolyDataPointNormalsPointer(cloud);


	float coveragePercentage_best = 0.0f;
	int numUnselectablePoints_best = 0;
	float coveragePercentage_seed = 0.0f;
	int numUnselectablePoints_seed = 0;
	//search for flat points until the algorithm covers covAreaP percentage of the mesh
	while((coveragePercentage_best <= covAreaP) && (coveragePercentage_seed <= covAreaP))
	{
		//select a random point
		int randomNum = rand()%n_cloud_points;
		int pointIndex = randomNum;

		//if the random point has already been extracted, search for the points at left or right of the random point until find a not-already-selected point
		while((!selectablePointSeed[pointIndex]))
		{
			if(rightLeft)
			{
				pointIndex = (++pointIndex) % n_cloud_points;
			}
			else
			{
				pointIndex--;
				if(pointIndex == -1)
				{
					pointIndex = n_cloud_points - 1;
				}
			}
		}
		rightLeft = !rightLeft;


		//remove the selected point (and its neighbours) from the list of seed points (dark green points in Fig.6 of paper)
		vtkIdList* NNpoints_remCenter = kdtree.FindPointsWithinRadius(pointIndex, exPRadius);
		int nNeighbors_remCenter = NNpoints_remCenter->GetNumberOfIds();
		vtkIdType* ptrIds_rem = NNpoints_remCenter->GetPointer(0);
		for(size_t iRC = 0; iRC < nNeighbors_remCenter; iRC++)
		{
			int indexRC = *ptrIds_rem;
			if(selectablePointSeed[indexRC])
			{
				numUnselectablePoints_seed++;
			}
			selectablePointSeed[indexRC] = false;

			ptrIds_rem++;
		}

		//find neighbour points (yellow points in Fig.6 of paper)
		if(partRadius == exPRadius)
		{
			NNpoints_partRadius->DeepCopy( NNpoints_remCenter );
		}
		else
		{
			NNpoints_partRadius->DeepCopy( kdtree.FindPointsWithinRadius(pointIndex, partRadius) );
		}
		vtkIdType* ptrIds_part = NNpoints_partRadius->GetPointer(0);
		int nNeighbors_partRadius = NNpoints_partRadius->GetNumberOfIds();

		if(nNeighbors_partRadius < m_minNeigh)
		{
			continue;					
		}

		//iterate on yellow points 
		float bestFlatness = -1.0f;
		int bestFlatnessIndex = -1;		
		float *candidateNormal = NULL;
		for(size_t index = 0; index < nNeighbors_partRadius; index++)
		{		
			int iPartPoint = *ptrIds_part;

			//check if the point is selectable as the best
			if(!selectablePointBest[iPartPoint])
			{
				++ptrIds_part;	
				continue;
			}

			//compute flatness if it hasn't already been computed
			float flatness = 0.0f;
			if(vFeaturePointScores[iPartPoint] != -2.0f)
			{					
				flatness = vFeaturePointScores[iPartPoint];
			}
			else
			{
				vtkIdList* NNpoints_detRadius = kdtree.FindPointsWithinRadius(iPartPoint, m_radius); 
				int nNeighbors_detRadius = NNpoints_detRadius->GetNumberOfIds();

				//check the minimum number of neighbors
				if( nNeighbors_detRadius < m_minNeigh)
				{
					vFeaturePointScores[iPartPoint] = -1.0f;
					++ptrIds_part;						
					continue;					
				}

				//compute flatness index
				vtkIdType* ptrIds_det = NNpoints_detRadius->GetPointer(0);

				candidateNormal = ptrNormals + 3 * iPartPoint;

				ptrIds_det = NNpoints_detRadius->GetPointer(0);	
				float meanCosine = 0.0f;

				//compute cosine index
				for(int ne = 0; ne < nNeighbors_detRadius; ++ne)
				{
					neighNormal = ptrNormals + 3* *ptrIds_det;

					cosine = candidateNormal[0] * neighNormal[0] + candidateNormal[1] * neighNormal[1] + candidateNormal[2] * neighNormal[2];
					meanCosine += cosine;

					++ptrIds_det;
				}				
				flatness = meanCosine / (float)nNeighbors_detRadius;


				vFeaturePointScores[iPartPoint] = flatness;	
			}


			if(bestFlatness < flatness)
			{
				bestFlatness = flatness;
				bestFlatnessIndex = iPartPoint;
			}

			++ptrIds_part;

		}

		//add flattest feature point index into the vector
		if(bestFlatnessIndex != -1)
		{
			vFeaturePointIndexes[numFeatures] = bestFlatnessIndex;
			numFeatures++;
		}


		//remove the flattest point (and its neighbours) from the list of best points (dark fuchsia points in Fig.6 of paper)
		vtkIdList* NNpoints_remBest = kdtree.FindPointsWithinRadius(bestFlatnessIndex, exPRadius);
		int nNeighbors_remBest = NNpoints_remBest->GetNumberOfIds();
		vtkIdType* ptrIdBest = NNpoints_remBest->GetPointer(0);
		for(size_t iRB = 0; iRB < nNeighbors_remBest; iRB++)
		{
			int indexRB = *ptrIdBest;					
			if(selectablePointBest[indexRB])
			{
				numUnselectablePoints_best++;
			}
			selectablePointBest[indexRB] = false;
			ptrIdBest++;

		}

		coveragePercentage_best = numUnselectablePoints_best/(float)n_cloud_points;
		coveragePercentage_seed = numUnselectablePoints_seed/(float)n_cloud_points;
		//cout << "coveragePercentage: " << coveragePercentage_best << "\t" << coveragePercentage_seed << " \r";

	}// while(coveragePercentage <= covAreaP)


	vFeaturePointIndexes.resize(numFeatures);


	NNpoints_partRadius->Delete();
}



ReLOC::Random3DDetector::Random3DDetector(int nPoints, bool random, double radius,  int minNeigh, double borderDistance, int64_t seed) : I3DDetector()
{
	m_requestedNumFeat = nPoints;
	m_minNeigh = minNeigh;
	m_radius = radius;
	m_borderDistance = borderDistance;


	m_seed = seed;
	if(random)
		srand(time(NULL));
	else
		srand(seed);

}

ReLOC::Random3DDetector::Random3DDetector(RandomDetectorParams params)
{
	m_requestedNumFeat = params.nFeatures;
	m_minNeigh = params.numMinNeighbors;
	m_radius = params.radius;
	m_borderDistance = params.borderDistance;

	m_seed = params.seed;
	if(params.random)
		srand(time(NULL));
	else
		srand(params.seed);
}

ReLOC::Random3DDetector::~Random3DDetector()
{
}



void ReLOC::Random3DDetector::extractImpl(vtkPolyData* cloud, Feature3D* & feat, int &numFeatures, int & featCapacity)
{

	double point[3];
	int randomNum;

	int n_cloud_points = cloud->GetNumberOfPoints();
	int detectablePointsLeft = n_cloud_points;

	bool *edgePointList;
	if(m_borderDistance > 0.0)
		GetBoundaryPoints(cloud, edgePointList);

	numFeatures = Min(m_requestedNumFeat, n_cloud_points);
	updateFeatCapacity(feat, featCapacity, numFeatures);

	bool checkCurv = false;
	double* ptrPrincCurv = NULL;
	if(cloud->GetPointData()->GetArray("Principal Curvatures"))
	{
		checkCurv = true;
		ptrPrincCurv = (double *)(cloud->GetPointData()->GetArray("Principal Curvatures")->GetVoidPointer(0));

	}

	KdTree kdtree;
	if(m_radius > 0.0 && m_radius != m_borderDistance && m_minNeigh > 0)
	{
		kdtree.SetPolyData(cloud);
	}
	vtkIdList* NNpoints;

	int extractedPoints = 0;

	std::set<int> ids;

	while ( extractedPoints < numFeatures){


		numFeatures = Min(m_requestedNumFeat, detectablePointsLeft);

		do
		{
			randomNum = rand()%n_cloud_points;
		}
		while( (checkCurv) && ( *(ptrPrincCurv + 2*randomNum) == std::numeric_limits<double>::infinity() ) );

		cloud->GetPoint(randomNum, point);

		if (ids.find(randomNum) != ids.end()) // found
			continue;
		ids.insert(randomNum);


		//Check the minimum number of neighbors (if radius != borderDistance)
		if(m_radius > 0.0 && m_radius != m_borderDistance && m_minNeigh > 0){

			NNpoints = kdtree.FindPointsWithinRadius(point, m_radius);
			int nNeighbors = NNpoints->GetNumberOfIds();

			if( nNeighbors-1 < m_minNeigh){
				detectablePointsLeft--;
				continue;
				//goto label_newpoint;
			}
		}

		//Check if the point is not too close to the border
		if(m_borderDistance > 0.0){

			NNpoints = kdtree.FindPointsWithinRadius(point, m_borderDistance);
			int nNeighbors = NNpoints->GetNumberOfIds();

			for(int j=0; j<nNeighbors; j++){
				if( edgePointList[ NNpoints->GetId(j) ]){
					detectablePointsLeft--;
					continue;
				}
			}

			//Check the minimum number of neighbors (if radius = borderDistance)
			if( nNeighbors < m_minNeigh && m_radius == m_borderDistance){
				detectablePointsLeft--;
				continue;
			}
		}


		feat[extractedPoints].x = point[0];
		feat[extractedPoints].y = point[1];
		feat[extractedPoints].z = point[2];


		feat[extractedPoints].scale = m_radius;
		feat[extractedPoints].index = randomNum;
		feat[extractedPoints].score = 1.0;

		extractedPoints++;

	}

	numFeatures = extractedPoints;

	if(m_borderDistance > 0.0)
		delete [] edgePointList;

}




ReLOC::I3DDescriptor::I3DDescriptor()
{
	m_descLength = 0;
}

ReLOC::I3DDescriptor::~I3DDescriptor()
{

}

void ReLOC::I3DDescriptor::describe(vtkPolyData* cloud, Feature3D* & feat, Desc & desc, int nFeat)
{
	describeImpl(cloud, feat, desc, nFeat);
}


void ReLOC::I3DDescriptor::describe(vtkPolyDataCollection* polyDataCollection, std::vector<Feature3D*> & vFeat, std::vector<Desc> & vDesc, std::vector<int> vNumFeat)
{
	vDesc.resize(polyDataCollection->GetNumberOfItems());
	polyDataCollection->InitTraversal();
	std::vector<Desc>::iterator itDesc = vDesc.begin();
	std::vector<int>::iterator itN = vNumFeat.begin();
	int currentLength = m_descLength;
	for(std::vector<Feature3D*>::iterator itFeat = vFeat.begin(); itFeat != vFeat.end(); ++itFeat )
	{
		vtkPolyData* poly = polyDataCollection->GetNextItem();
		describeImpl(poly, *itFeat, *itDesc, *itN);
		m_descLength = currentLength;
		printf("\nNum descriptors... %d\n", *itN);
		itDesc++;
		itN++;
	}
}



ReLOC::Desc::Desc()
{
	init();
}

ReLOC::Desc::Desc(const char *path)
{
	init();
	load(path);
}

ReLOC::Desc::Desc(const Desc & other)
{
	init();
	*this = other;
}

ReLOC::Desc::Desc(const Desc & other, bool forceDeepCopy)
{
	init();

	if(other.getNumDesc() == 0)
	{
		return;
	}

	if(!other.m_ptrOwner && !forceDeepCopy)
	{
		m_size = other.m_size;
		m_length = other.m_length;
		m_capacity = 0;
		m_ptrOwner = false;
		m_doubleInit = other.m_doubleInit;
		m_bitsetInit = other.m_bitsetInit;
		m_floatInit = other.m_floatInit;
		m_doubleDesc = other.m_doubleDesc;
		m_floatDesc = other.m_floatDesc;
		return;
	}

	if(other.isDouble())
	{
		double **desc_d = updateDoubleInternal(other.getNumDesc(), other.getLength());
		double const* const* other_d = other.getDoubleReadOnly();
		memcpy(desc_d[0], other_d[0], sizeof(double)*m_length*m_size);
	}
	if(other.isFloat())
	{
		float **desc_f = updateFloatInternal(other.getNumDesc(), other.getLength());
		float const* const* other_f = other.getFloatReadOnly();
		memcpy(desc_f[0], other_f[0], sizeof(float)*m_length*m_size);
	}
}

ReLOC::Desc::~Desc()
{
	release();
}

void ReLOC::Desc::init()
{
	m_size = 0;
	m_length = 0;
	m_capacity = 0;
	m_ptrOwner = true;
	m_doubleInit = false;
	m_bitsetInit = false;
	m_floatInit = false;
	m_doubleDesc = NULL;
	m_floatDesc = NULL;
}

void ReLOC::Desc::release()
{
  if(m_ptrOwner)
  {
    if(m_doubleInit)
    {
      if (m_doubleDesc != NULL)
      {
        delete [] m_doubleDesc[0];
        delete [] m_doubleDesc;
      }
    }

    if(m_floatInit)
    {
      if (m_floatDesc != NULL)
      {
        delete [] m_floatDesc[0];
        delete [] m_floatDesc;
      }
    }

  }

  init();
}

ReLOC::Desc & ReLOC::Desc::operator =(const Desc & other)
{
	if(this != &other)
	{
		release();

		if(other.getNumDesc() == 0)
		{
			return *this;
		}

		if(!other.m_ptrOwner)
		{
			m_size = other.m_size;
			m_length = other.m_length;
			m_capacity = 0;
			m_ptrOwner = false;
			m_doubleInit = other.m_doubleInit;
			m_bitsetInit = other.m_bitsetInit;
			m_floatInit = other.m_floatInit;
			m_doubleDesc = other.m_doubleDesc;
			m_floatDesc = other.m_floatDesc;
			return *this;
		}

		if(other.isDouble())
		{
			double **desc_d = updateDoubleInternal(other.getNumDesc(), other.getLength());
			double const* const * other_d = other.getDoubleReadOnly();
			memcpy(desc_d[0], other_d[0], sizeof(double)*m_length*m_size);
		}
		if(other.isFloat())
		{
			float **desc_f = updateFloatInternal(other.getNumDesc(), other.getLength());
			float const* const* other_f = other.getFloatReadOnly();
			memcpy(desc_f[0], other_f[0], sizeof(float)*m_length*m_size);
		}
	}

	return *this;
}

ReLOC::Desc & ReLOC::Desc::operator +=(const Desc & other)
{
	if(this != &other)
	{
		if(!m_doubleInit && !m_bitsetInit && !m_floatInit)
		{
			deepCopyFrom(other);
		}
		else if(m_doubleInit && other.isDouble())
		{
			// copy original descriptors
			Desc desc_cpy(*this, true);

			// reallocate desc for the new size
			double **desc_d = updateDouble(desc_cpy.getNumDesc() + other.getNumDesc(), desc_cpy.getLength());

			// merge
			double const* const * desc_cpy_d = desc_cpy.getDoubleReadOnly();
			double const* const * other_d = other.getDoubleReadOnly();
			memcpy(desc_d[0], desc_cpy_d[0], sizeof(double) * desc_cpy.getNumDesc() * desc_cpy.getLength());
			memcpy(desc_d[0] + desc_cpy.getNumDesc() * desc_cpy.getLength(), other_d[0], sizeof(double) * other.getNumDesc() * other.getLength());

			// release old desc
			desc_cpy.release();
		}
		else if(m_floatInit && other.isFloat())
		{
			// copy original descriptors
			Desc desc_cpy(*this, true);

			// reallocate desc for the new size
			float **desc_f = updateFloat(desc_cpy.getNumDesc() + other.getNumDesc(), desc_cpy.getLength());

			// merge
			float const* const * desc_cpy_f = desc_cpy.getFloatReadOnly();
			float const* const * other_f = other.getFloatReadOnly();
			memcpy(desc_f[0], desc_cpy_f[0], sizeof(float) * desc_cpy.getNumDesc() * desc_cpy.getLength());
			memcpy(desc_f[0] + desc_cpy.getNumDesc() * desc_cpy.getLength(), other_f[0], sizeof(float) * other.getNumDesc() * other.getLength());

			// release old desc
			desc_cpy.release();
		}
		else
		{
			std::cout << std::endl << "[ERROR] Trying to merge descriptors of different type." << std::endl;
			throw std::logic_error("Trying to merge descriptors of different type.");
		}
	}

	return *this;
}

double** ReLOC::Desc::updateDoubleInternal(int size, int length)
{
	if(!m_doubleInit)
		m_doubleInit = true;

	if (m_capacity < size || length > m_length)
	{
		if (m_doubleDesc != NULL && m_ptrOwner)
		{
			delete [] m_doubleDesc[0];
			delete [] m_doubleDesc;
		}
		m_doubleDesc = NULL;
	}

	if (m_doubleDesc == NULL && size > 0)
	{
		m_doubleDesc = new double*[size];
		m_doubleDesc[0] = new double[size*length];
		m_ptrOwner = true;

		for (int i = 1; i < size; ++i)
			m_doubleDesc[i] = m_doubleDesc[i-1] + length;

		m_capacity = size;
	}

	for (int i = 0; i < size; i++)
		for (int j = 0; j < length; j++)
			m_doubleDesc[i][j] = 0.0;

	m_size = size;
	m_length = length;
	return m_doubleDesc;
}

double** ReLOC::Desc::updateDouble(int size, int length)
{
	if(m_floatInit)
	{
		if (m_floatDesc != NULL && m_ptrOwner)
		{
			delete [] m_floatDesc[0];
			delete [] m_floatDesc;
		}
		m_floatInit = false;
		m_floatDesc = NULL;
	}

	return updateDoubleInternal(size, length);
}

float** ReLOC::Desc::updateFloatInternal(int size, int length)
{
	if(!m_floatInit)
		m_floatInit = true;

	if (m_capacity < size || length > m_length)
	{
		if (m_floatDesc != NULL && m_ptrOwner)
		{
			delete [] m_floatDesc[0];
			delete [] m_floatDesc;
		}
		m_floatDesc = NULL;
	}

	if (m_floatDesc == NULL && size > 0)
	{
		m_floatDesc = new float*[size];
		m_floatDesc[0] = new float[size*length];
		m_ptrOwner = true;

		for (int i = 1; i < size; ++i)
			m_floatDesc[i] = m_floatDesc[i-1] + length;

		m_capacity = size;
	}

	for (int i = 0; i < size; i++)
		for (int j = 0; j < length; j++)
			m_floatDesc[i][j] = 0.0;

	m_size = size;
	m_length = length;
	return m_floatDesc;
}

float** ReLOC::Desc::updateFloat(int size, int length)
{
	if(m_doubleInit)
	{
		if (m_doubleDesc != NULL && m_ptrOwner)
		{
			delete [] m_doubleDesc[0];
			delete [] m_doubleDesc;
		}
		m_doubleInit = false;
		m_doubleDesc = NULL;
	}

	return updateFloatInternal(size, length);
}

double** ReLOC::Desc::getDouble()
{
	if (m_floatInit && !m_doubleInit)
	{
		bool wasOwner = m_ptrOwner;
		updateDoubleInternal(m_size, m_length);

		for(int i = 0; i < m_size; i++)
			for(int j = 0; j < m_length; j++)
				m_doubleDesc[i][j] = m_floatDesc[i][j];

		if (!wasOwner)
		{
			m_floatInit = false;
			m_floatDesc = NULL;
		}
	}

	return m_doubleDesc;
}

float** ReLOC::Desc::getFloat()
{
	if (!m_floatInit && m_doubleInit)
	{
    bool wasOwner = m_ptrOwner;
		updateFloatInternal(m_size, m_length);

		for(int i = 0; i < m_size; i++)
			for(int j = 0; j < m_length; j++)
				m_floatDesc[i][j] = (float)m_doubleDesc[i][j];

    if (!wasOwner)
		{
      m_doubleInit = false;
      m_doubleDesc = NULL;
    }
	}

	return m_floatDesc;
}

double const * const * ReLOC::Desc::getDoubleReadOnly() const
{
	if (!m_doubleInit)
		return NULL;

	return m_doubleDesc;
}

float const * const * ReLOC::Desc::getFloatReadOnly() const
{
	if (!m_floatInit)
		return NULL;

	return m_floatDesc;
}


bool ReLOC::Desc::save(const char *path, const bool binaryFile, const bool append)
{
	if(m_doubleInit && m_bitsetInit || m_floatInit && m_bitsetInit || m_doubleInit && m_floatInit)
	{
		std::cout << std::endl << "[FATAL ERROR] ReLOC::Desc::save(): more than one type of descriptor defined." << std::endl;
		throw std::logic_error("ReLOC::Desc::save(): more than one type of descriptor defined.");
	}

	CreateFullPath(GetPathDir(path));

	if(binaryFile)
	{
		//binary version
		fstream binary_file;

		if(append)
		{
			binary_file.open(path, ios::app|ios::binary);
		}
		else
		{
			binary_file.open(path, ios::out|ios::binary);
		}

		if (!binary_file.is_open() )
		{
			return false;
		}

		if(m_doubleInit)
		{
			binary_file << "double ";
		}
		else if(m_floatInit)
		{
			binary_file << "float ";
		}
		else if(m_bitsetInit)
		{
			binary_file << "bitset ";
		}

		binary_file << m_size << " " << m_length << "\n";
		//binary_file.write((char*)(&m_size), sizeof(m_size));
		//binary_file.write((char*)(&m_length), sizeof(m_length));

		if(m_doubleInit)
		{
			for(int i=0; i<m_size; i++)
			{
				for(int j=0; j<m_length; j++)
				{
					binary_file.write((char*)(&(m_doubleDesc[i][j])), sizeof(m_doubleDesc[i][j]));
				}
			}
		}
		else if(m_floatInit)
		{
			for(int i=0; i<m_size; i++)
			{
				for(int j=0; j<m_length; j++)
				{
					binary_file.write((char*)(&(m_floatDesc[i][j])), sizeof(m_floatDesc[i][j]));
				}
			}
		}

		binary_file.close();
	}
	else
	{
		//ascii version
		std::ofstream file;

		if(append)
		{
			file.open(path, ios::app);
		}
		else
		{
			file.open(path);
		}

		file.precision(30);

		if (!file.is_open())
		{
			return false;
		}

		if(m_doubleInit)
		{
			file << "double ";
		}
		else if(m_floatInit)
		{
			file << "float ";
		}
		else if(m_bitsetInit)
		{
			file << "bitset ";
		}

		file << m_size << " " << m_length << "\n";

		for(int i = 0; i < m_size; i++)
		{
			if(m_doubleInit)
			{
				for(int j = 0; j < m_length; j++)
				{
					file << m_doubleDesc[i][j] << " ";
				}
			}
			else if(m_floatInit)
			{
				for(int j = 0; j < m_length; j++)
				{
					file << m_floatDesc[i][j] << " ";
				}
			}
		}

		file.close();
	}

	return true;
}

bool ReLOC::Desc::save(const char *path, const char *fileTemplate, std::vector<Desc> & vDesc, const bool binaryFile)
{
	std::stringstream filename;
	for(size_t de=0; de<vDesc.size(); de++)
	{
		filename.str("");
		filename << path << "/" << fileTemplate << std::setw(3) << std::setfill('0') << de << ".desc";
		if( !vDesc[de].save(filename.str().c_str(), binaryFile))
		{
			return false;
		}
	}
	return true;
}

bool ReLOC::Desc::save(vector<string> &vAbsFileNames, std::vector<Desc> & vDesc, const bool binaryFile)
{
	for(size_t de=0; de<vDesc.size(); de++)
	{
		if( !vDesc[de].save(ChangeExtension(vAbsFileNames[de].c_str(), "desc").c_str(), binaryFile))
		{
			return false;
		}
	}
	return true;
}

bool ReLOC::Desc::save(const char *path, vector<string> &vFileNames, std::vector<Desc> & vDesc, const bool binaryFile)
{
	std::stringstream filename;
	for(size_t de=0; de<vDesc.size(); de++)
	{
		filename.str("");
		filename << path << "/" << ChangeExtension(vFileNames[de].c_str(), "desc").c_str();
		if( !vDesc[de].save(filename.str().c_str(), binaryFile))
		{
			return false;
		}
	}
	return true;
}

bool ReLOC::Desc::load(const char *path, const bool binaryFile)
{
	release();

	if(binaryFile)
	{
		//binary version
		fstream binary_file(path, ios::binary|ios::in);

		if (!binary_file.is_open() )
		{
			return false;
		}

		std::string type;
		binary_file >> type;
		binary_file >> m_size;
		binary_file >> m_length;
		char name[10];
		binary_file.getline(name,10);

		//binary_file.read((char*)(&m_size), sizeof(m_size));
		//binary_file.read((char*)(&m_length), sizeof(m_length));

		if(type.compare("double") == 0)
		{
			updateDoubleInternal(m_size, m_length);
		}
		else if(type.compare("float") == 0)
		{
			updateFloatInternal(m_size, m_length);
		}

		if(m_doubleInit)
		{
			for(int i=0; i<m_size; i++)
			{
				for(int j=0; j<m_length; j++)
				{
					binary_file.read((char*)(&(m_doubleDesc[i][j])), sizeof(m_doubleDesc[i][j]));
				}
			}
		}
		else if(m_floatInit)
		{
			for(int i=0; i<m_size; i++)
			{
				for(int j=0; j<m_length; j++)
				{
					binary_file.read((char*)(&(m_floatDesc[i][j])), sizeof(m_floatDesc[i][j]));
				}
			}
		}

		binary_file.close();
	}
	else
	{
		//ascii version
		std::ifstream file(path);

		if (!file.is_open() )
		{
			return false;
		}

		std::string type;
		file >> type;
		file >> m_size;
		file >> m_length;
		std::string strLine;
		std::getline(file, strLine);

		if(type.compare("double") == 0)
		{
			updateDoubleInternal(m_size, m_length);
		}
		else if(type.compare("float") == 0)
		{
			updateFloatInternal(m_size, m_length);
		}

		std::string bitset_str;

		for(int i = 0; i < m_size; i++)
		{
			if(m_doubleInit)
			{
				for(int j = 0; j < m_length; j++)
				{
					file >> m_doubleDesc[i][j];
				}
			}
			else if(m_floatInit)
			{
				for(int j = 0; j < m_length; j++)
				{
					file >> m_floatDesc[i][j];
				}
			}
		}

		file.close();
	}

	return true;
}

bool ReLOC::Desc::load(const char *path, const char *fileTemplate, std::vector<Desc> & vDesc, const bool binaryFile)
{
	std::vector<std::string> foundFiles;
	FindFilesStartWith(path, fileTemplate, foundFiles);
	vDesc.resize((const unsigned int)foundFiles.size());
	for(size_t de=0; de<foundFiles.size(); de++)
	{
		if(!vDesc[de].load(foundFiles[de].c_str(), binaryFile))
			return false;
	}
	return true;
}

bool ReLOC::Desc::load(vector<string> &vAbsFileNames, std::vector<Desc> & vDesc, const bool binaryFile)
{
	vDesc.resize((const unsigned int)vAbsFileNames.size());
	vDesc.resize( vAbsFileNames.size() );
	for(size_t de=0; de<vAbsFileNames.size(); de++)
	{
		if(!vDesc[de].load(ChangeExtension(vAbsFileNames[de].c_str(), "desc").c_str(), binaryFile))
		//if(!load(ChangeExtension(vAbsFileNames[de].c_str(), "desc").c_str(), m_vDesc[de], m_vNumDesc[de], m_descLength, m_vDescCapacity[de], binaryFile))
		{
			return false;
		}
	}
	return true;
}

bool ReLOC::Desc::load(const char *path, vector<string> &vFileNames, std::vector<Desc> & vDesc, const bool binaryFile)
{
	std::stringstream filename;
	vDesc.resize( vFileNames.size() );
	for(size_t de=0; de<vFileNames.size(); de++)
	{
		filename.str("");
		filename << path << "/" << ChangeExtension(vFileNames[de].c_str(), "desc").c_str();

		if(!vDesc[de].load(filename.str().c_str(), binaryFile))
		{
			return false;
		}
	}
	return true;
}


ReLOC::Desc ReLOC::Desc::segment(unsigned start, unsigned n) const
{
  Desc d_copy;
  d_copy.m_size = n;
  d_copy.m_length = m_length;
  d_copy.m_capacity = 0;
  d_copy.m_ptrOwner = false;
  d_copy.m_doubleInit = m_doubleInit;
  d_copy.m_bitsetInit = m_bitsetInit;
  d_copy.m_floatInit = m_floatInit;
  d_copy.m_doubleDesc = NULL;
  d_copy.m_floatDesc = NULL;
  if(m_doubleDesc != NULL)
    d_copy.m_doubleDesc = &m_doubleDesc[start];
  if(m_floatDesc != NULL)
    d_copy.m_floatDesc = &m_floatDesc[start];
  return d_copy;
}

ReLOC::Desc& ReLOC::Desc::deepCopyFrom(const Desc& other)
{
  if(this != &other)
	{
		release();

		if(other.getNumDesc() == 0)
		{
			return *this;
		}

		if(other.isDouble())
		{
			double **desc_d = updateDoubleInternal(other.getNumDesc(), other.getLength());
			double const* const * other_d = other.getDoubleReadOnly();
			memcpy(desc_d[0], other_d[0], sizeof(double)*m_length*m_size);
		}
		if(other.isFloat())
		{
			float **desc_f = updateFloatInternal(other.getNumDesc(), other.getLength());
			float const* const* other_f = other.getFloatReadOnly();
			memcpy(desc_f[0], other_f[0], sizeof(float)*m_length*m_size);
		}
	}

	return *this;
}




ReLOC::LRFDescriptor_1::LRFDescriptor_1(LRFDescriptor_1_Params params) : I3DDescriptor()
{
	setParams(params);
}

ReLOC::LRFDescriptor_1::~LRFDescriptor_1()
{

}


void ReLOC::LRFDescriptor_1::describeImpl(vtkPolyData *cloud, Feature3D *&feat, Desc & desc, int nPoints)
{
	KdTree kdtree(cloud);

	float *ptrDesc = *(desc.updateFloat(nPoints, m_descLength));

	Feature3D* ptrFeat = feat;
	if(m_params.flareParams.selector == "FLARE_Sampled_T")
	{
		KdTree kdTree_sampled;
		PopulateKdTree_NaiveSampling(cloud, m_params.flareParams.samplingPerc_tangent, kdTree_sampled);

		for(int po=0; po<nPoints; po++)
		{
			*ptrDesc = getLocalRF_FLARE_Sampled_T(cloud, kdtree, kdTree_sampled, ptrFeat->index, ptrFeat->rf, m_params.flareParams, m_params.meshRes);
			ptrDesc++;
			ptrFeat++;
		}
	}
	else
	{
		for(int po=0; po<nPoints; po++)
		{
			*ptrDesc = getLocalRF_FLARE(cloud, kdtree, ptrFeat->index, ptrFeat->rf, m_params.flareParams, m_params.meshRes);
			ptrDesc++;
			ptrFeat++;
		}
	}
}


float ReLOC::getLocalRF_FLARE_Sampled_T(vtkPolyData *cloud, KdTree &kdTree, KdTree &kdTree_sampled, const int index, float *rfc, FLAREparams &params, const float meshRes)
{
	vtkIdType* ptrIds;
	float* ptrPoints = GetPolyDataPointsPointer(cloud);
	float* ptrNormals = GetPolyDataPointNormalsPointer(cloud );
	float* neighPoint;
	float* neighNormal;


	//Z axis
	vtkIdList* NNpoints = kdTree.FindPointsWithinRadius(index, (float)(params.radiusInmeshRes_normal * meshRes));
	int nNeighbours = NNpoints->GetNumberOfIds();

	if (nNeighbours < params.minNumNeighbours_normal)
	{
		memcpy(rfc+6, ptrNormals + 3*index, sizeof(float)*3);
	}
	else
	{
		//find centroid for plane fitting
		float centroidP[] = {0.0f, 0.0f, 0.0f};
		float meanNormal[] = {0.0f, 0.0f, 0.0f};
		ptrIds = NNpoints->GetPointer(0);
		for(int ne = 0; ne < nNeighbours; ++ne)
		{
			neighPoint = ptrPoints + 3* *ptrIds;
			neighNormal = ptrNormals + 3* *ptrIds;
			
			centroidP[0] += neighPoint[0];
			centroidP[1] += neighPoint[1];
			centroidP[2] += neighPoint[2];

			meanNormal[0] += neighNormal[0];
			meanNormal[1] += neighNormal[1];
			meanNormal[2] += neighNormal[2];

			++ptrIds;
		}
		centroidP[0] /= (float)nNeighbours;
		centroidP[1] /= (float)nNeighbours;
		centroidP[2] /= (float)nNeighbours;

		vtkMath::Normalize(meanNormal);

		//plane fitting
		float noCentroid[3];
		memset(params.covM[0], 0, sizeof(float)*3);
		memset(params.covM[1], 0, sizeof(float)*3);
		memset(params.covM[2], 0, sizeof(float)*3);

		ptrIds = NNpoints->GetPointer(0);
		float temp;
		for(int ne = 0; ne < nNeighbours; ++ne)
		{
			neighPoint = ptrPoints + 3* *ptrIds;

			// Difference between current point and origin
			noCentroid[0] = neighPoint[0] - centroidP[0];
			noCentroid[1] = neighPoint[1] - centroidP[1];
			noCentroid[2] = neighPoint[2] - centroidP[2];

			//populate covariance matrix
			params.covM[0][0] += noCentroid[0] * noCentroid[0];
			params.covM[1][1] += noCentroid[1] * noCentroid[1];
			params.covM[2][2] += noCentroid[2] * noCentroid[2];

			temp =  noCentroid[0] * noCentroid[1];
			params.covM[0][1] += temp;
			params.covM[1][0] += temp;
			temp =  noCentroid[0] * noCentroid[2];
			params.covM[0][2] += temp;
			params.covM[2][0] += temp;
			temp =  noCentroid[1] * noCentroid[2];
			params.covM[1][2] += temp;
			params.covM[2][1] += temp;

			++ptrIds;
		}


		float eval[3];
		int resJ = vtkMath::Jacobi(params.covM, eval, params.evect);
		if( resJ != 1)
		{
			cout << "ERROR (getLocalRF_Fast): resJ != 1 " << resJ;
		}

		rfc[6] = params.evect[0][2];
		rfc[7] = params.evect[1][2];
		rfc[8] = params.evect[2][2];

		//disambiguate Z axis with mean normal
		if (rfc[6]*meanNormal[0] + rfc[7]*meanNormal[1] + rfc[8]*meanNormal[2] < 0)
		{
			rfc[6] = -rfc[6];
			rfc[7] = -rfc[7];
			rfc[8] = -rfc[8];
		}

	}



	//X axis
	float* featurePoint = ptrPoints + 3 * index;
	float* ptrPoints_sampled = GetPolyDataPointsPointer( kdTree_sampled.GetPolyData() );

	NNpoints = kdTree_sampled.FindPointsWithinRadius(featurePoint, params.radiusInmeshRes_tangent * meshRes);
	nNeighbours = NNpoints->GetNumberOfIds();

	if(nNeighbours < params.minNumNeighbours_tangent)
	{
		RandomOrthogonalAxis(rfc+6, rfc);

		vtkMath::Cross(&rfc[6], &rfc[0], &rfc[3]);
		return std::numeric_limits<float>::max();
	}

	float shapeScore;
	float bestShapeScore = -std::numeric_limits<float>::max();
	int bestShapeIndex = -1;

	float neighDistance;
	
	float radius2 = params.radiusInmeshRes_tangent * meshRes * params.radiusInmeshRes_tangent * meshRes;
	float crownDistance2 = params.crownThresh * params.crownThresh * radius2;

	ptrIds = NNpoints->GetPointer(0);
	float euclSupp_x;
	float euclSupp_y;
	float euclSupp_z;
	for(int ne = 0; ne < nNeighbours; ++ne)
	{
		neighPoint = ptrPoints_sampled + 3* *ptrIds;
		
		//EuclideanDistance2_3D
		euclSupp_x = neighPoint[0]-featurePoint[0];
		euclSupp_y = neighPoint[1]-featurePoint[1];
		euclSupp_z = neighPoint[2]-featurePoint[2];
		neighDistance = euclSupp_x*euclSupp_x + euclSupp_y*euclSupp_y + euclSupp_z*euclSupp_z;

		if(neighDistance > crownDistance2)
		{
			shapeScore = neighPoint[0] * rfc[6] + neighPoint[1] * rfc[7] + neighPoint[2] * rfc[8];

			if(shapeScore > bestShapeScore)
			{
				bestShapeIndex = *ptrIds;
				bestShapeScore = shapeScore;
			}
		
		}

		++ptrIds;
	}

	if(bestShapeIndex == -1)
	{
		RandomOrthogonalAxis(rfc+6, rfc);
		vtkMath::Cross(&rfc[6], &rfc[0], &rfc[3]);
		return std::numeric_limits<float>::max();
	}

	DirectedOrthogonalAxis(rfc+6, featurePoint, ptrPoints_sampled + 3 * bestShapeIndex, rfc);
	vtkMath::Cross(rfc+6, rfc+0, rfc+3);

	bestShapeScore -= (rfc[6]*featurePoint[0] + rfc[7]*featurePoint[1] + rfc[8]*featurePoint[2]);
	return bestShapeScore/meshRes;

}



float ReLOC::getLocalRF_FLARE(vtkPolyData *cloud, KdTree &kdTree, const int index, float *rfc, FLAREparams &params, const float meshRes)
{
	vtkIdType* ptrIds;
	float* ptrPoints = GetPolyDataPointsPointer(cloud);
	float* ptrNormals = GetPolyDataPointNormalsPointer(cloud );
	float* neighPoint;
	float* neighNormal;



	//Z axis
	vtkIdList* NNpoints = kdTree.FindPointsWithinRadius(index, (float)(params.radiusInmeshRes_normal * meshRes));
	int nNeighbours = NNpoints->GetNumberOfIds();

	if (nNeighbours < params.minNumNeighbours_normal)
	{
		memcpy(rfc+6, ptrNormals + 3*index, sizeof(float)*3);
	}
	else
	{
		//find centroid for plane fitting
		float centroidP[] = {0.0f, 0.0f, 0.0f};
		float meanNormal[] = {0.0f, 0.0f, 0.0f};
		ptrIds = NNpoints->GetPointer(0);
		for(int ne = 0; ne < nNeighbours; ++ne)
		{
			neighPoint = ptrPoints + 3* *ptrIds;
			neighNormal = ptrNormals + 3* *ptrIds;
			
			centroidP[0] += neighPoint[0];
			centroidP[1] += neighPoint[1];
			centroidP[2] += neighPoint[2];

			meanNormal[0] += neighNormal[0];
			meanNormal[1] += neighNormal[1];
			meanNormal[2] += neighNormal[2];

			++ptrIds;
		}
		centroidP[0] /= (float)nNeighbours;
		centroidP[1] /= (float)nNeighbours;
		centroidP[2] /= (float)nNeighbours;

		vtkMath::Normalize(meanNormal);

		//plane fitting
		float noCentroid[3];
		memset(params.covM[0], 0, sizeof(float)*3);
		memset(params.covM[1], 0, sizeof(float)*3);
		memset(params.covM[2], 0, sizeof(float)*3);

		ptrIds = NNpoints->GetPointer(0);
		float temp;
		for(int ne = 0; ne < nNeighbours; ++ne)
		{
			neighPoint = ptrPoints + 3* *ptrIds;

			// Difference between current point and origin
			noCentroid[0] = neighPoint[0] - centroidP[0];
			noCentroid[1] = neighPoint[1] - centroidP[1];
			noCentroid[2] = neighPoint[2] - centroidP[2];

			//populate covariance matrix
			params.covM[0][0] += noCentroid[0] * noCentroid[0];
			params.covM[1][1] += noCentroid[1] * noCentroid[1];
			params.covM[2][2] += noCentroid[2] * noCentroid[2];

			temp =  noCentroid[0] * noCentroid[1];
			params.covM[0][1] += temp;
			params.covM[1][0] += temp;
			temp =  noCentroid[0] * noCentroid[2];
			params.covM[0][2] += temp;
			params.covM[2][0] += temp;
			temp =  noCentroid[1] * noCentroid[2];
			params.covM[1][2] += temp;
			params.covM[2][1] += temp;

			++ptrIds;
		}


		float eval[3];
		int resJ = vtkMath::Jacobi(params.covM, eval, params.evect);
		if( resJ != 1)
		{
			cout << "ERROR (getLocalRF_Fast): resJ != 1 " << resJ;
		}

		rfc[6] = params.evect[0][2];
		rfc[7] = params.evect[1][2];
		rfc[8] = params.evect[2][2];

		//disambiguate Z axis with mean normal
		if (rfc[6]*meanNormal[0] + rfc[7]*meanNormal[1] + rfc[8]*meanNormal[2] < 0)
		{
			rfc[6] = -rfc[6];
			rfc[7] = -rfc[7];
			rfc[8] = -rfc[8];
		}

	}




	//X axis
	NNpoints = kdTree.FindPointsWithinRadius(index, params.radiusInmeshRes_tangent * meshRes);
	nNeighbours = NNpoints->GetNumberOfIds();

	if(nNeighbours < params.minNumNeighbours_tangent)
	{
		//RandomOrthogonalAxis(rfc+6, rfc);
		//vtkMath::Cross(&rfc[6], &rfc[0], &rfc[3]);
		rfc[0] = std::numeric_limits<float>::quiet_NaN();
		return std::numeric_limits<float>::max();
	}

	float shapeScore;
	float bestShapeScore = -std::numeric_limits<float>::max();
	int bestShapeIndex = -1;

	float neighDistance;
	
	float radius2 = params.radiusInmeshRes_tangent * meshRes * params.radiusInmeshRes_tangent * meshRes;
	float crownDistance2 = params.crownThresh * params.crownThresh * radius2;

	ptrIds = NNpoints->GetPointer(0);
	float* featurePoint = ptrPoints + 3 * index;
	float euclSupp_x;
	float euclSupp_y;
	float euclSupp_z;
	for(int ne = 0; ne < nNeighbours; ++ne)
	{
		neighPoint = ptrPoints + 3* *ptrIds;
		
		//EuclideanDistance2_3D
		euclSupp_x = neighPoint[0]-featurePoint[0];
		euclSupp_y = neighPoint[1]-featurePoint[1];
		euclSupp_z = neighPoint[2]-featurePoint[2];
		neighDistance = euclSupp_x*euclSupp_x + euclSupp_y*euclSupp_y + euclSupp_z*euclSupp_z;

		if(neighDistance > crownDistance2)
		{
			shapeScore = neighPoint[0] * rfc[6] + neighPoint[1] * rfc[7] + neighPoint[2] * rfc[8];

			if(shapeScore > bestShapeScore)
			{
				bestShapeIndex = *ptrIds;
				bestShapeScore = shapeScore;
			}
		
		}

		++ptrIds;
	}

	if(bestShapeIndex == -1)
	{
		//RandomOrthogonalAxis(rfc+6, rfc);
		//vtkMath::Cross(&rfc[6], &rfc[0], &rfc[3]);
		rfc[0] = std::numeric_limits<float>::quiet_NaN();
		return std::numeric_limits<float>::max();
	}

	DirectedOrthogonalAxis(rfc+6, featurePoint, GetPolyDataPointsPointer(cloud, bestShapeIndex), rfc);
	vtkMath::Cross(rfc+6, rfc+0, rfc+3);

	bestShapeScore -= (rfc[6]*featurePoint[0] + rfc[7]*featurePoint[1] + rfc[8]*featurePoint[2]);
	return bestShapeScore/meshRes;

}




//********************* LRFMatcher_1 *********************

ReLOC::LRFMatcher_1::LRFMatcher_1(float signedDistanceThresh): IMatcher()
{
	m_signedDistanceThresh = signedDistanceThresh;
}

ReLOC::LRFMatcher_1::LRFMatcher_1(LRFMatcher_1_Params params):   IMatcher(),
	m_signedDistanceThresh(params.signedDistanceThresh)
{
}

ReLOC::LRFMatcher_1::~LRFMatcher_1()
{

}

void ReLOC::LRFMatcher_1::Sort_wrtFirstFeat(Desc & desc, vector<int> &vId_sorted)
{
	int descLength = desc.getLength();
	int nDescs = desc.getNumDesc();

	vector<float> firstFeats(nDescs);

	float* ptrDesc_1 = &firstFeats[0];
	float* ptrDescs = *(desc.getFloat());
	for(int de=0; de<nDescs; de++)
	{
		*ptrDesc_1 = *ptrDescs;
		ptrDescs += descLength;
		ptrDesc_1++;
	}

	//allocate ids vector for sorting
	vId_sorted.resize(nDescs);
	for(int id=0; id<nDescs; id++)
	{
		vId_sorted[id] = id; 
	}

	//sort features w.r.t. first feature
	std::sort(vId_sorted.begin(), vId_sorted.end(), SortIds_Wrt_FloatDecr(firstFeats) );
}

void ReLOC::LRFMatcher_1::ReplaceIds(FeatureMatch* &featureMatch, const int nMatches, const vector<int> &vNewIds_trg, const vector<int> &vNewIds_ref)
{

	FeatureMatch* ptrFeatureMatch = featureMatch;
	int const * ptrvNewIds_trg = &vNewIds_trg[0];
	int const * ptrvNewIds_ref = &vNewIds_ref[0];
	for(int ma=0; ma<nMatches; ma++)
	{
		ptrFeatureMatch->refIndex = vNewIds_ref[ ptrFeatureMatch->refIndex ];
		ptrFeatureMatch->trgIndex = vNewIds_trg[ ptrFeatureMatch->trgIndex ];
		ptrFeatureMatch++;
		ptrvNewIds_trg++;
		ptrvNewIds_ref++;
	}

}

bool ReLOC::LRFMatcher_1::setReferenceImpl(Desc & refDescriptors)
{
	m_desc_ref = refDescriptors;

	//find sorting ids w.r.t. descriptor first feature
	Sort_wrtFirstFeat(m_desc_ref, m_vId_sorted_descRef);

	//sort w.r.t. indeces
	m_desc_ref.sort(m_vId_sorted_descRef);

	return true;
}

int ReLOC::LRFMatcher_1::matchImpl(Desc & trgDescriptors, FeatureMatch* & featureMatch)
{
	if(featureMatch != NULL)
	{
		delete [] featureMatch;
		featureMatch = NULL;
	}

	//sort trgDescriptors
	Desc trgDescriptors_sorted = trgDescriptors;
	
	//find sorting ids w.r.t. descriptor first feature
	vector<int> vId_sorted_descTrg;
	Sort_wrtFirstFeat(trgDescriptors_sorted, vId_sorted_descTrg);

	//sort w.r.t. indeces
	trgDescriptors_sorted.sort(vId_sorted_descTrg);


	float* desc_trg = *(trgDescriptors_sorted.getFloat());
	float* desc_ref = *(m_desc_ref.getFloat());
	int numDesc_trg = trgDescriptors_sorted.getNumDesc();
	int numDesc_ref = m_desc_ref.getNumDesc();

	//find maxDiff_ref_trg for viewpair
	float maxScore_trg, maxScore_ref, minScore_trg, minScore_ref;

	float* ptrDesc_trg = desc_trg;
	int iDesc = 0;
	while( (maxScore_trg = *(ptrDesc_trg++) ) == std::numeric_limits<float>::max() )
	{
		iDesc++;
	}
	if(iDesc == numDesc_trg)
	{
		//no valid descs in trg
		return 0;
	}

	float* ptrDesc_ref = desc_ref;
	iDesc = 0;
	while( (maxScore_ref = *(ptrDesc_ref++) ) == std::numeric_limits<float>::max() )
	{
		iDesc++;
	}
	if(iDesc == numDesc_ref)
	{
		//no valid descs in ref
		return 0;
	}

	ptrDesc_trg = desc_trg + numDesc_trg - 1;
	while( (minScore_trg = *(ptrDesc_trg--) ) == std::numeric_limits<float>::max() ) {}

	ptrDesc_ref = desc_ref + numDesc_ref - 1;
	while( (minScore_ref = *(ptrDesc_ref--) ) == std::numeric_limits<float>::max() ) {}

	float maxDiff_ref_trg = Max(fabs(maxScore_trg - minScore_ref), fabs(maxScore_ref - minScore_trg));


	float minMin = Min(minScore_ref, minScore_trg);
	float maxMax = Max(maxScore_ref, maxScore_trg);

	float delta = maxDiff_ref_trg * m_signedDistanceThresh;
	//float flatThresh = maxMax * flatThreshPerc;

	float minDelta, maxDelta;


	
	//find num correct matches
	int goodMatches = 0;
	int fe1 = 0;
	bool endScan = false;

	//skip max values
	ptrDesc_trg = desc_trg;
	int fe0_start = 0;
	while( *(ptrDesc_trg++)  == std::numeric_limits<float>::max() ) {fe0_start++;}
	ptrDesc_trg--;
	

	//skip max values
	int iMinDelta = 0;
	while( desc_ref[iMinDelta++] == std::numeric_limits<float>::max() ) {}
	iMinDelta--;

	for(int fe0=fe0_start; fe0 < numDesc_trg; fe0++)
	{
		minDelta = *ptrDesc_trg + delta;
		maxDelta = *ptrDesc_trg - delta;

		while(desc_ref[iMinDelta] > minDelta)
		{
			iMinDelta++;
			if(iMinDelta == numDesc_ref)
			{
				//desc_ref[numDesc_ref-1] does not enter in trg window and all scores in trg will be lower than desc_ref[iMinDelta]
				endScan = true;
				break;
			}
		}
		if(endScan)
		{
			//desc_ref[numDesc_ref-1] does not enter in trg window and all scores in trg will be lower than desc_ref[iMinDelta]
			break;
		}
		fe1 = iMinDelta;
		while( (fe1<numDesc_ref) && (desc_ref[fe1] > maxDelta) )
		{
			//float score = (*ptrDesc_trg + desc_ref[fe1])/2.0;

			//if(fabs(score) < flatThresh)
			//{
			//	fe1++;
			//	continue;
			//}

			goodMatches++;
			fe1++;
		}
		ptrDesc_trg++;
	}

	if(goodMatches == 0)
	{
		return 0;
	}

	//resize matches
	featureMatch = new FeatureMatch[goodMatches];


	//populate matches
	FeatureMatch* ptrMatches = featureMatch;
	fe1 = 0;
	endScan = false;

	//skip max values
	ptrDesc_trg = desc_trg;
	while( *(ptrDesc_trg++)  == std::numeric_limits<float>::max() ) {}
	ptrDesc_trg--;

	//skip max values
	iMinDelta = 0;
	while( desc_ref[iMinDelta++] == std::numeric_limits<float>::max() ) {}
	iMinDelta--;

	for(int fe0=fe0_start; fe0 < numDesc_trg; fe0++)
	{
		minDelta = *ptrDesc_trg + delta;
		maxDelta = *ptrDesc_trg - delta;

		while(desc_ref[iMinDelta] > minDelta)
		{
			iMinDelta++;
			if(iMinDelta == numDesc_ref)
			{
				//desc_ref[numDesc_ref-1] does not enter in trg window and all scores in trg will be lower than desc_ref[iMinDelta]
				endScan = true;
				break;
			}
		}
		if(endScan)
		{
			//desc_ref[numDesc_ref-1] does not enter in trg window and all scores in trg will be lower than desc_ref[iMinDelta]
			break;
		}
		fe1 = iMinDelta;
		while( (fe1<numDesc_ref) && (desc_ref[fe1] > maxDelta) )
		{

			ptrMatches->matchScore = 1.0f - fabs(*ptrDesc_trg - desc_ref[fe1] )/maxDiff_ref_trg;
			ptrMatches->refIndex = fe1;
			ptrMatches->trgIndex = fe0;


			ptrMatches++;
			fe1++;
		}
		ptrDesc_trg++;
	}

	ReplaceIds(featureMatch, goodMatches, vId_sorted_descTrg, m_vId_sorted_descRef);


	return goodMatches;
}






ReLOC::IMatcher::IMatcher()
{
	m_featureMatch = NULL;
	m_numMatches = 0;

	m_nTotalDesc = 0;

	m_matcherStrategy = E_MATCHER_STRATEGY_RATIO;
	m_knn = -1;
	m_featureDistance = -1.0;

	m_referenceInitialized = false;
	m_numRefDescriptors = 0;
}

ReLOC::IMatcher::~IMatcher()
{
	if(m_featureMatch != NULL)
	{
		delete [] m_featureMatch;
		m_featureMatch = NULL;
	}

	if(m_vFeatureMatch.size() > 0)
	{
		for(unsigned int i = 0; i < m_vFeatureMatch.size(); i++)
		{
			if(m_vFeatureMatch[i] != NULL)
			{
				delete [] m_vFeatureMatch[i];
				m_vFeatureMatch[i] = NULL;
			}
		}
	}

	// reset single index data
	reset();
}

void ReLOC::IMatcher::SetSingleIndexRatioMatching(const std::vector<int> & refModelIndices)
{
	if(refModelIndices.empty())
	{
		std::cerr << "[IMatcher::SetSingleIndexRatioMatching] ERROR: THE REFERENCE INDEX VECTOR IS EMPTY!!" << std::endl;
		throw std::logic_error("[IMatcher::SetSingleIndexRatioMatching] ERROR: THE REFERENCE INDEX VECTOR IS EMPTY!!");
	}

	m_matcherStrategy = E_MATCHER_STRATEGY_SINGLE_INDEX_RATIO;
	m_refModelIndices = refModelIndices;
}

void ReLOC::IMatcher::SetSingleIndexRatioMatching()
{
	if(m_refModelIndices.empty())
	{
		std::cerr << "[IMatcher::SetSingleIndexRatioMatching] ERROR: Model indices vector empty! This mode can be used only if a SINGLE INDEX is built using the 'add' and 'train' methods!" << std::endl;
		throw std::logic_error("[IMatcher::SetSingleIndexRatioMatching] ERROR: Model indices vector empty! This mode can be used only if a SINGLE INDEX is built using the 'add' and 'train' methods!");
	}

	m_matcherStrategy = E_MATCHER_STRATEGY_SINGLE_INDEX_RATIO;
}

bool ReLOC::IMatcher::setReference(Desc & refDescriptors)
{
	m_numRefDescriptors = refDescriptors.getNumDesc();

	m_referenceInitialized = (m_numRefDescriptors != 0);
	if(m_matcherStrategy == E_MATCHER_STRATEGY_RATIO)
	{
		m_referenceInitialized = (m_numRefDescriptors > 1);
	}

	if (m_referenceInitialized)
	{
		m_referenceInitialized = setReferenceImpl(refDescriptors);
	}

	return m_referenceInitialized;
}

void ReLOC::IMatcher::reset()
{

	for(int md = 0; md < (int) m_vDescForSingleIndex.size(); md++)
	{
		m_vDescForSingleIndex[md].release();
	}

	m_vDescForSingleIndex.clear();
	m_nTotalDesc = 0;
}


int ReLOC::IMatcher::getModelMatches(const FeatureMatch * matches, const int nMatches, const int id, FeatureMatch * & modelFm)
{
	if(modelFm != NULL)
	{
		delete [] modelFm;
		modelFm = NULL;
	}

	int mm = 0;

	for(int ma = 0; ma < nMatches; ma++)
	{
		if( m_refModelIndices[matches[ma].refIndex] == id)
		{
			mm++;
		}
	}

	if(mm == 0)
	{
		return 0;
	}

	modelFm = new FeatureMatch[mm];

	mm = 0;
	
	for(int ma = 0; ma < nMatches; ma++)
	{
		if( m_refModelIndices[matches[ma].refIndex] == id)
		{
			modelFm[mm] = (matches[ma]);
			mm++;
		}
	}

	return mm;
}

int ReLOC::IMatcher::match(Desc & trgDescriptors, FeatureMatch *&featureMatch)
{
	if (!m_referenceInitialized || trgDescriptors.getNumDesc() == 0)
		return 0;

	m_numMatches = matchImpl(trgDescriptors, m_featureMatch);
	featureMatch = m_featureMatch;
	return m_numMatches;
}

void ReLOC::IMatcher::match(std::vector<Desc> &vDesc, std::vector< std::pair<int,int> > & pairs)
{
	resizeVect((const unsigned int)(pairs.size()));

	PairSort sortPredicate;
	std::sort(pairs.begin(), pairs.end(), sortPredicate);

	std::vector<FeatureMatch*>::reverse_iterator itFeatMatch = m_vFeatureMatch.rbegin();
	std::vector<int>::reverse_iterator itNumMatches = m_vNumMatches.rbegin();
	std::vector< std::pair<int,int> >::reverse_iterator itPairs = pairs.rbegin();

	int currentRef = itPairs->second;
	setReference(vDesc[currentRef]);

//	FeatureMatch* ptrFeatureMatch = NULL;
	for(; itPairs != pairs.rend(); ++itPairs )
	{
		if(itPairs->second != currentRef) {
			currentRef = itPairs->second;
			setReference(vDesc[currentRef]);
		}
		printf("Matches %d %d ... ", itPairs->first, itPairs->second);

		if (!m_referenceInitialized || vDesc[itPairs->first].getNumDesc() == 0)
			*itNumMatches = 0;
		else
			*itNumMatches = matchImpl(vDesc[itPairs->first], *itFeatMatch);

		printf("OK  Num Matches: %d\n", *itNumMatches);

		itFeatMatch++;
		itNumMatches++;
	}
}

void ReLOC::IMatcher::resizeVect(const unsigned int newLength)
{
	for(unsigned int i = 0; i < m_vFeatureMatch.size(); i++)
	{
		if(m_vFeatureMatch[i] != NULL)
		{
			delete [] m_vFeatureMatch[i];
			m_vFeatureMatch[i] = NULL;
		}
	}

	m_vNumMatches.resize(newLength);
	m_vFeatureMatch.resize(newLength);

	for(size_t i = m_vFeatureMatch.size(); i < newLength; i++)
	{
		m_vFeatureMatch[i] = NULL;
	}
}

// Removes matches with same refIndex leaving the one with the best score (since the matches
// are in descending order, we take the first for each refIndex as it is the best one)
// The number of reference features is required in order to speed up the computation
int ReLOC::IMatcher::cleanFeatureMatch()
{
	int reducedNumGoodMatches = 0;

	bool *checked = new bool[m_numRefDescriptors];
	memset(checked, 0, sizeof(bool) * m_numRefDescriptors);

	for(int ma = 0; ma < m_numMatches; ma++)
	{
		if(!checked[m_featureMatch[ma].refIndex])
		{
			checked[m_featureMatch[ma].refIndex] = true;
			m_featureMatch[reducedNumGoodMatches] = m_featureMatch[ma];
			reducedNumGoodMatches++;
		}
	}

	delete [] checked;

	m_numMatches = reducedNumGoodMatches;
	return m_numMatches;
}

void ReLOC::Desc::sort(const vector<int> &vId_sorted)
{
	float *ptrDesc = getFloat()[0];
	
	int descLength = getLength();
	int nDescs = getNumDesc();

	if((size_t)nDescs != vId_sorted.size())
	{
		cout << "ERROR: (Desc::sort) nDescs != vId_sorted.size()   " << nDescs << " != " << vId_sorted.size();
		getchar();
		exit(-1);
	}

	if(isFloat())
	{
		vector<float>		vDesc_beforeSorting(nDescs*descLength);

		//copy not sorted descriptors 
		memcpy(&vDesc_beforeSorting[0], ptrDesc, sizeof(float)*nDescs*descLength);
	
		//populate descriptors w.r.t. ids
		float* ptrDesc_beforeSorting = &vDesc_beforeSorting[0];
		const int *ptrId = &vId_sorted[0];

		for(int po=0; po<nDescs; po++)
		{
			memcpy(ptrDesc, ptrDesc_beforeSorting + (*ptrId) * descLength, sizeof(float)*descLength);
			ptrId++;
			ptrDesc += descLength;
		}
	}

	if(isDouble())
	{
		vector<double>		vDesc_beforeSorting(nDescs*descLength);

		//copy not sorted descriptors 
		memcpy(&vDesc_beforeSorting[0], ptrDesc, sizeof(double)*nDescs*descLength);
	
		//populate descriptors w.r.t. ids
		double* ptrDesc_beforeSorting = &vDesc_beforeSorting[0];
		const int *ptrId = &vId_sorted[0];

		for(int po=0; po<nDescs; po++)
		{
			memcpy(ptrDesc, ptrDesc_beforeSorting + (*ptrId) * descLength, sizeof(double)*descLength);
			ptrId++;
			ptrDesc += descLength;
		}
	}

	if(isBitset())
	{
		cout << "Desc::sort for bitset: Do it yourself";
		getchar();
	}
}

