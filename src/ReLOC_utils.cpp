#include "ReLOC_utils.h"


#include <algorithm>


using namespace std;




#include <io.h>   // For access().
#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().
//
//#include <cassert>
//
//
//#include "vtkTransformFilter.h"
#include "vtkIdTypeArray.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
//#include "vtkPolyDataWriter.h"
#include "vtkErrorCode.h"
#include "vtkPLYReader.h"
//#include "vtkOBJReader.h"
//#include "vtkSTLReader.h"
//#include "vtkPolyDataReader.h"
//#include "vtkSimplePointsReader.h"
#include "vtkPolyDataNormals.h"
#include "vtkCleanPolyData.h"
#include "vtkKdTreePointLocator.h"
#include "vtkPointLocator.h"
#include "vtkMath.h"
#include "vtkPointData.h"
//#include "vtkTriangle.h"
#include <vtkVersion.h>


void ReLOC::GetCentroid(vtkPolyData *polyData, float centroid[3])
{
	double centroid_double[] = {0.0, 0.0, 0.0};
	float* ptrPoints = (float*)(polyData->GetPoints()->GetVoidPointer(0));
	int nPoints = polyData->GetNumberOfPoints();
	for(int po=0; po<nPoints; po++)
	{
		centroid_double[0] += *ptrPoints;
		ptrPoints++;
		centroid_double[1] += *ptrPoints;
		ptrPoints++;
		centroid_double[2] += *ptrPoints;
		ptrPoints++;
	}
	centroid[0] = (float)(centroid_double[0]/nPoints);
	centroid[1] = (float)(centroid_double[1]/nPoints);
	centroid[2] = (float)(centroid_double[2]/nPoints);
}


void ReLOC::GetCentroids_cMeshes(vtkPolyDataCollection* coll, std::vector<float> &vCentroids)
{
	vCentroids.resize(coll->GetNumberOfItems()*3);
	coll->InitTraversal();
	vtkPolyData* poly;
	for(int po=0; po<coll->GetNumberOfItems(); po++)
	{
		poly = coll->GetNextItem();
		GetCentroid(poly, &(vCentroids[3*po]));
	}
}


void ReLOC::GetBBox_StdDev_Clouds(vtkPolyDataCollection* cMeshes, std::vector<double> & vBBoxes, const std::vector<float> &vCentroids, const float sigmaFactor)
{
	vBBoxes.resize(cMeshes->GetNumberOfItems()*6);
	int me=0;
	cMeshes->InitTraversal();
	vtkPolyData* poly;
	while( (poly = cMeshes->GetNextItem()) )
	{
		GetBBox_StdDev_Cloud(poly, &vBBoxes[me*6], &vCentroids[me*3], sigmaFactor);
		me++;
	}
}

void ReLOC::GetBBox_StdDev_Cloud(vtkPolyData* cloud, double* bbox, const float* centroid, const float sigmaFactor)
{
	float* ptrPoints = GetPolyDataPointsPointer( cloud );
	vtkIdType nPoints = cloud->GetNumberOfPoints();
	vector<double> vXs(nPoints);
	vector<double> vYs(nPoints);
	vector<double> vZs(nPoints);
	double* ptrXs = &vXs[0];
	double* ptrYs = &vYs[0];
	double* ptrZs = &vZs[0];
	for(vtkIdType po=0; po<nPoints; po++)
	{
		*ptrXs++ = *ptrPoints++;
		*ptrYs++ = *ptrPoints++;
		*ptrZs++ = *ptrPoints++;
	}

	//find sigma for each dimension
	double sigma_X = StdDev(vXs, (double)centroid[0], 0.0);
	double sigma_Y = StdDev(vYs, (double)centroid[1], 0.0);
	double sigma_Z = StdDev(vZs, (double)centroid[2], 0.0);

	bbox[0] = centroid[0] - sigma_X * sigmaFactor;
	bbox[1] = centroid[0] + sigma_X * sigmaFactor;
	bbox[2] = centroid[1] - sigma_Y * sigmaFactor;
	bbox[3] = centroid[1] + sigma_Y * sigmaFactor;
	bbox[4] = centroid[2] - sigma_Z * sigmaFactor;
	bbox[5] = centroid[2] + sigma_Z * sigmaFactor;
}

float* ReLOC::GetPolyDataPointsPointer(vtkPolyData* polyData, const unsigned int iPoint)
{
	vtkPoints* Array = polyData->GetPoints();

	if(Array->GetNumberOfPoints() == 0 )
		return NULL;
	else
		return (float*)(Array->GetVoidPointer(iPoint*3));
}



double ReLOC::calcDistance2BetweenPoints(float* points1, float* points2, int nPoints)
{
	double rmse = 0;
	float* ptrPoints1 = points1;
	float* ptrPoints2 = points2;
	for(int po=0; po<nPoints; po++)
	{
		rmse += EuclideanDistance2_3D(ptrPoints1, ptrPoints2);
		ptrPoints1 += 3;
		ptrPoints2 += 3;
	}

	return rmse;
}


vtkPolyData* ReLOC::CreatePolyData(vtkPolyData* polyData, vtkIdType* pointIds, vtkIdType nPointIds, const bool insertNormals, const bool insertColors, const bool insertTriangles)
{
	vtkPolyData* poly = vtkPolyData::New();

	vtkIdType* newPointsPos = new vtkIdType[polyData->GetNumberOfPoints()];
	for(int po=0; po<polyData->GetNumberOfPoints(); po++)
	{
		newPointsPos[po] = -1;
	}

	//insert points
	AllocatePoints(poly, nPointIds);
	float* ptrNewPoints = (float*)(poly->GetPoints()->GetVoidPointer(0));
	float* ptrOldPoints = (float*)(polyData->GetPoints()->GetVoidPointer(0));
	for(int po=0; po< nPointIds; po++)
	{
		memcpy(ptrNewPoints, ptrOldPoints+ 3 * pointIds[po], sizeof(float)*3);
		ptrNewPoints+=3;
		newPointsPos[(pointIds[po])] = po;
	}

	//insert point normals
	if(insertNormals)
	{
		float* polyDataNormals = GetPolyDataPointNormalsPointer(polyData);
		if(polyDataNormals != NULL)
		{
			float* polyNormals = AllocateNormals(poly);
			for(int i=0; i<nPointIds; i++)
			{
				polyNormals[3*i] = polyDataNormals[3*pointIds[i]];
				polyNormals[3*i + 1] = polyDataNormals[3*pointIds[i] + 1];
				polyNormals[3*i + 2] = polyDataNormals[3*pointIds[i] + 2];
			}
		}
	}

	//insert point colors
	if(insertColors)
	{
		unsigned char* polyDataColors = GetPolyDataColorsPointer(polyData);
		if(polyDataColors != NULL)
		{
			unsigned char* polyColors = AllocateColors(poly);
			for(int i=0; i<nPointIds; i++)
			{
				polyColors[3*i] = polyDataColors[3*pointIds[i]];
				polyColors[3*i + 1] = polyDataColors[3*pointIds[i] + 1];
				polyColors[3*i + 2] = polyDataColors[3*pointIds[i] + 2];
			}
		}
	}


	//insert triangles
	if(insertTriangles)
	{
		if(polyData->GetNumberOfPolys() > 0)
		{
			vtkIdTypeArray* newCells = vtkIdTypeArray::New();
			newCells->SetNumberOfComponents(4);
			newCells->SetNumberOfTuples(polyData->GetNumberOfPolys());
			vtkIdType* ptrOldCells = (vtkIdType *)(polyData->GetPolys()->GetPointer());
			vtkIdType* ptrNewCells = (vtkIdType *)(newCells->GetPointer(0));
			int nNewCells = 0;
			for(int ce=0; ce<polyData->GetNumberOfPolys(); ce++)
			{
				if( (newPointsPos[ptrOldCells[1]]!=-1) && (newPointsPos[ptrOldCells[2]]!=-1) && (newPointsPos[ptrOldCells[3]]!=-1) )
				{
					*ptrNewCells = 3;
					ptrNewCells++;
					*ptrNewCells = newPointsPos[ptrOldCells[1]];
					ptrNewCells++;
					*ptrNewCells = newPointsPos[ptrOldCells[2]];
					ptrNewCells++;
					*ptrNewCells = newPointsPos[ptrOldCells[3]];
					ptrNewCells++;
					nNewCells++;
				}
				ptrOldCells += 4;
			}
			newCells->Resize(nNewCells);

			vtkCellArray* cellArray = vtkCellArray::New();
			cellArray->SetCells(nNewCells, newCells);
			AllocateTriangles(poly, nNewCells);
			poly->SetPolys(cellArray);
			cellArray->Delete();
			newCells->Delete();
		}
	}


	delete[] newPointsPos;
	return poly;
}




float* ReLOC::AllocatePoints(vtkPolyData *polyData, vtkIdType nPoints)
{
	vtkFloatArray* floatArray = vtkFloatArray::New();
	floatArray->SetNumberOfComponents(3);
	floatArray->SetNumberOfTuples(nPoints);

	vtkPoints* points3D = vtkPoints::New();
	points3D->SetData(floatArray);

	polyData->SetPoints(points3D);

	floatArray->Delete();
	points3D->Delete();

	return (float*)(polyData->GetPoints()->GetVoidPointer(0));

}

float* ReLOC::GetPolyDataPointNormalsPointer(vtkPolyData* polyData, const unsigned int iNormal)
{
	vtkDataArray* Array = polyData->GetPointData()->GetNormals();

	if(Array == NULL)
		return NULL;
	else
		return (float*)(Array->GetVoidPointer(iNormal*3));
}


float* ReLOC::AllocateNormals(vtkPolyData *polyData)
{

	vtkFloatArray* normals = vtkFloatArray::New();
	normals->SetNumberOfComponents( 3 );
	normals->SetNumberOfTuples( polyData->GetNumberOfPoints() );
	normals->SetName("Normals");
	polyData->GetPointData()->SetNormals(normals);
	normals->Delete();
	return (float *)(polyData->GetPointData()->GetNormals()->GetVoidPointer(0));
}

unsigned char* ReLOC::GetPolyDataColorsPointer(vtkPolyData* polyData, const unsigned int iColor)
{
	vtkDataArray* Array = polyData->GetPointData()->GetArray("RGB");

	if(Array == NULL)
		return NULL;
	else
		return (unsigned char*)(Array->GetVoidPointer(iColor*3));
}

unsigned char* ReLOC::AllocateColors(vtkPolyData *polyData, unsigned char* color, bool* mask)
{
	vtkUnsignedCharArray* colorArray = vtkUnsignedCharArray::New();
	colorArray->SetNumberOfComponents( 3 );
	colorArray->SetName("RGB");
	colorArray->SetNumberOfTuples( polyData->GetNumberOfPoints() );
	if(color)
	{
		unsigned char* ptrColors = (unsigned char *)(colorArray->GetVoidPointer(0));
		unsigned char noColor[3] = {200, 200, 200};
		for(int co=0; co<polyData->GetNumberOfPoints(); co++)
		{
			if(mask && !mask[co])
			{
				memcpy(ptrColors, noColor, sizeof(unsigned char)*3);
			}
			else
			{
				memcpy(ptrColors, color, sizeof(unsigned char)*3);
			}
			ptrColors+=3;
		}
	}
	polyData->GetPointData()->AddArray(colorArray);
	polyData->GetPointData()->SetActiveScalars("RGB");
	colorArray->Delete();
	return (unsigned char *)(polyData->GetPointData()->GetScalars("RGB")->GetVoidPointer(0));
}


vtkIdType* ReLOC::AllocateTriangles(vtkPolyData *polyData, const vtkIdType nTriangles)
{
	vtkIdTypeArray *idTypeTrianglesArray = vtkIdTypeArray::New();
	idTypeTrianglesArray->SetNumberOfComponents(4);
	idTypeTrianglesArray->SetNumberOfTuples(nTriangles);

	vtkCellArray *trianglesCellArray = vtkCellArray::New();
	trianglesCellArray->SetCells( nTriangles, idTypeTrianglesArray);

	polyData->SetPolys(trianglesCellArray);

	idTypeTrianglesArray->Delete();
	trianglesCellArray->Delete();

	return polyData->GetPolys()->GetPointer();

}


vtkPolyData* ReLOC::CreatePolyData(const float* points, const vtkIdType nPoints, const unsigned char* colors, const vtkIdType* triangles, const vtkIdType nTriangles)
{
	vtkPolyData* polyData = vtkPolyData::New();

	//insert points
	AllocatePoints(polyData, nPoints);
	memcpy(polyData->GetPoints()->GetVoidPointer(0), points, nPoints*3 * sizeof(float) );

	if(colors)
	{
		//insert colors
		AllocateColors(polyData);
		memcpy(polyData->GetPointData()->GetArray("RGB")->GetVoidPointer(0), colors, nPoints* 3 * sizeof(unsigned char) );
	}

	if(triangles)
	{
		//insert triangles
		AllocateTriangles(polyData, nTriangles);

		vtkIdType* ptrTrianglesArray = (vtkIdType *)polyData->GetPolys()->GetPointer();
		const vtkIdType* ptrTriangles = triangles;
		for(vtkIdType tr=0; tr<nTriangles; tr++)
		{
			*(ptrTrianglesArray++)= 3;
			memcpy(ptrTrianglesArray, ptrTriangles, sizeof(vtkIdType)*3);
			ptrTriangles+=3;
			ptrTrianglesArray+=3;
		}
	}
	return polyData;
}

void ReLOC::Translate(float* points, int nPoints, const float point[3])
{
	float* ptrPoints = points;
	for(int po=0; po<nPoints; po++)
	{
		*ptrPoints += point[0];
		ptrPoints++;
		*ptrPoints += point[1];
		ptrPoints++;
		*ptrPoints += point[2];
		ptrPoints++;
	}
}



void ReLOC::ApplyTransform(float* points, int nPoints, const double* matrTransf)
{
	float* ptrPoints = points;
	float pointIn[4];
	float pointOut[4];
	pointIn[3] = 1;
	for(int po=0; po<nPoints; po++)
	{
		memcpy(pointIn, ptrPoints, sizeof(float)*3);
		vtkMatrix4x4::MultiplyPoint(matrTransf, pointIn, pointOut);
		memcpy(ptrPoints, pointOut, sizeof(float)*3);
		ptrPoints+=3;
	}
}


void ReLOC::GetCentroid(float* points, int nPoints, float centroid[3])
{
	double centroid_double[] = {0.0, 0.0, 0.0};
	float* ptrPoints = points;
	for(int po=0; po<nPoints; po++)
	{
		centroid_double[0] += *ptrPoints;
		ptrPoints++;
		centroid_double[1] += *ptrPoints;
		ptrPoints++;
		centroid_double[2] += *ptrPoints;
		ptrPoints++;
	}
	centroid[0] = (float)(centroid_double[0]/nPoints);
	centroid[1] = (float)(centroid_double[1]/nPoints);
	centroid[2] = (float)(centroid_double[2]/nPoints);
}




ReLOC::KdTree::KdTree(const bool useTrueKdTree)
{
	if(useTrueKdTree)
	{
		m_kdTree = vtkKdTreePointLocator::New();
	}
	else
	{
		m_kdTree = vtkPointLocator::New();
	}
	m_pointsList = vtkIdList::New();
}


ReLOC::KdTree::KdTree(vtkPolyData* polyData, const bool useTrueKdTree)
{
	if(useTrueKdTree)
	{
		m_kdTree = vtkKdTreePointLocator::New();
	}
	else
	{
		m_kdTree = vtkPointLocator::New();
	}
	m_pointsList = vtkIdList::New();
	SetPolyData(polyData);
}

void ReLOC::KdTree::SetPolyData(vtkPolyData* polyData)
{
	m_polyData = polyData;
	m_kdTree->SetDataSet(polyData);
	m_kdTree->Update();
}

ReLOC::KdTree::~KdTree()
{
	m_kdTree->Delete();
	m_pointsList->Delete();

}


vtkIdList* ReLOC::KdTree::FindPointsWithinRadius(const float* const point, float radius )
{
	const double p[3] = {point[0], point[1], point[2]};
	m_kdTree->FindPointsWithinRadius(radius, p, m_pointsList);
	return m_pointsList;
}

vtkIdList* ReLOC::KdTree::FindPointsWithinRadius(const int pointIndex, float radius )
{
	float* point = (float*)(m_polyData->GetPoints()->GetVoidPointer(pointIndex*3));
	const double p[3] = {point[0], point[1], point[2]};
	m_kdTree->FindPointsWithinRadius(radius, p, m_pointsList);
	return m_pointsList;
}



vtkIdList* ReLOC::KdTree::FindPointsWithinRadius(const double* point, double radius )
{
	const double p[3] = {point[0], point[1], point[2]};

	m_kdTree->FindPointsWithinRadius(radius, p, m_pointsList);
	return m_pointsList;
}

int ReLOC::KdTree::FindNearestPoint(double* point)
{
	return m_kdTree->FindClosestPoint(point);
}

int ReLOC::KdTree::FindNearestPoint(float* point)
{
	const double p[3] = {point[0], point[1], point[2]};
	return m_kdTree->FindClosestPoint(p);
}

int ReLOC::KdTree::FindNearestPoint(const int pointIndex)
{
	float* point = (float*)(m_polyData->GetPoints()->GetVoidPointer(pointIndex*3));
	const double p[3] = {point[0], point[1], point[2]};
	return m_kdTree->FindClosestPoint(p);
}

int ReLOC::KdTree::FindNearestPointWithinRadius(double* point, double radius, double & dist)
{
	return m_kdTree->FindClosestPointWithinRadius(radius, point, dist);
}

int ReLOC::KdTree::FindNearestPointWithinRadius(float* point, float radius, double & dist)
{
	double doublePoint[] = {point[0],point[1],point[2]};
	return m_kdTree->FindClosestPointWithinRadius(radius, doublePoint, dist);
}

int ReLOC::KdTree::FindNearestPointWithinRadius(const int pointIndex, float radius, double & dist)
{
	float* point = (float*)(m_polyData->GetPoints()->GetVoidPointer(pointIndex*3));
	double doublePoint[] = {point[0],point[1],point[2]};
	return m_kdTree->FindClosestPointWithinRadius(radius, doublePoint, dist);
}


double ReLOC::KdTree::FindNearestPoint(double* point, int &nearestPointId)
{
	nearestPointId = m_kdTree->FindClosestPoint(point);
	double nearestPoint[3];
	m_polyData->GetPoint(nearestPointId, nearestPoint);
	return sqrt(vtkMath::Distance2BetweenPoints(nearestPoint, point));
}

double ReLOC::KdTree::FindNearestPoint(float* point, int &nearestPointId)
{
	const double p[3] = {point[0], point[1], point[2]};
	nearestPointId = m_kdTree->FindClosestPoint(p);
	float* nearestPoint = GetPolyDataPointsPointer(m_polyData, nearestPointId);
	return sqrt(vtkMath::Distance2BetweenPoints(nearestPoint, point));
}

double ReLOC::KdTree::FindNearestPoint(const int pointIndex, int &nearestPointId)
{
	float* point = (float*)(m_polyData->GetPoints()->GetVoidPointer(pointIndex*3));
	const double p[3] = {point[0], point[1], point[2]};
	nearestPointId = m_kdTree->FindClosestPoint(p);
	float* nearestPoint = GetPolyDataPointsPointer(m_polyData, nearestPointId);
	return sqrt(vtkMath::Distance2BetweenPoints(nearestPoint, point));
}

vtkIdList* ReLOC::KdTree::FindNearestNPoints(int n, const double* point)
{
	const double p[3] = {point[0], point[1], point[2]};
	m_kdTree->FindClosestNPoints(n, p, m_pointsList);
	return m_pointsList;
}

vtkIdList* ReLOC::KdTree::FindNearestNPoints(int n, const float* point)
{
	const double p[3] = {point[0], point[1], point[2]};
	m_kdTree->FindClosestNPoints(n, p, m_pointsList);
	return m_pointsList;
}

vtkIdList* ReLOC::KdTree::FindNearestNPoints(int n, const int pointIndex)
{
	float* point = (float*)(m_polyData->GetPoints()->GetVoidPointer(pointIndex*3));
	const double p[3] = {point[0], point[1], point[2]};

	m_kdTree->FindClosestNPoints(n, p, m_pointsList);
	return m_pointsList;
}


std::string ReLOC::ChangeExtension(const std::string &filename, const std::string &newExt)
{
	string strTemp = filename;
	size_t dotPos = strTemp.find_last_of(".");
	if(dotPos == std::string::npos)
	{
		strTemp = strTemp + "." + newExt;
	}
	else
	{
		strTemp.replace(dotPos+1, strTemp.size(), newExt);
	}
	return strTemp;
}


bool ReLOC::ExistsDir(const string &directory)
{
	if ( _access( directory.c_str(), 0 ) == 0 )
    {
        struct stat status;
        stat( directory.c_str(), &status );

        if ( status.st_mode & S_IFDIR )
        {
            return true;
        }

		//The path you entered is a file.
		return false;
    }

	return false;
}

void ReLOC::CreateDir(const string &directory)
{
	#ifdef WIN32
		CreateDirectory(directory.c_str(), NULL);
	#elif defined (__GNUC__)
	/*  read/write/exe for owner
        read/exe for group owner
        read for other              */
        mkdir(directory.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH);
	#elif 
		#ifndef (BOOST_IS_NOT_INCLUDED)
        create_directory(path(directory.c_str()));
		#endif
	#endif
}

void ReLOC::CreateFullPath(const string &directory)
{
	string strDirectory = directory;
	if( (strDirectory[strDirectory.size()-1] != ':') && (!ExistsDir(directory)) )
	{
		string subDir = GetPathDir(directory);
		CreateFullPath(subDir.c_str());
		CreateDir(directory);
	}
}

string ReLOC::GetPathDir(const string &filename)
{
	string strFinal = filename;
	size_t posSlash = strFinal.rfind('/');
	size_t posBackSlash = strFinal.rfind('\\');
	if( (posSlash == string::npos) && (posBackSlash == string::npos) )
	{
		return ".";
	}

	if (posSlash == string::npos)
	{
		return strFinal.substr(0, posBackSlash);
	}

	if (posBackSlash == string::npos)
	{
		return strFinal.substr(0, posSlash);
	}

	return strFinal.substr(0, (posSlash > posBackSlash)?posSlash:posBackSlash);

}


bool ReLOC::FindFilesStartWith(const std::string & path, const std::string & startingString, std::vector<std::string> & foundFiles, bool getOnlyFileName, const int nMaxItems )
{
	foundFiles.clear();

#ifdef _MSC_VER
	// Required structs for searching for files and directories
	WIN32_FIND_DATA FindFileData;
	HANDLE hFind = INVALID_HANDLE_VALUE;

	// Build the file search string...
	char searchDir[2048] = {0};
	char fullpath[2048] = {0};

	// ...if it already is a path that ends with \ or /, add '*'...
	if(path.at(path.length() - 1) == '\\' || path.at(path.length() - 1) == '/')
	{
		_snprintf(searchDir, 2047, "%s*", path.c_str());
		_snprintf(fullpath, 2047, "%s", path.c_str()); // just copy path
	}
	// ...otherwise, add '\*' to the end of the path.
	else
	{
		_snprintf(searchDir, 2047, "%s/*", path.c_str());
		_snprintf(fullpath, 2047, "%s/", path.c_str()); // copy path and add slash (required when building filenames)
	}

	// Find the first file in the directory.
	hFind = FindFirstFile(searchDir, &FindFileData);

	// If there is no file, return
	if (hFind == INVALID_HANDLE_VALUE)
	{
		return false;
	}

	// loop
	do
	{
		// Skip ".", ".." and all directories
		if(FindFileData.dwFileAttributes == FILE_ATTRIBUTE_DIRECTORY || strcmp(FindFileData.cFileName, ".") == 0 || strcmp(FindFileData.cFileName, "..") == 0)
		{
			continue;
		}

		// Store filename into the vector if it starts with the given string
		if(startingString.size() > 0)
		{
			if(std::string(FindFileData.cFileName).find(startingString) == 0)
			{
				if(getOnlyFileName)
				{
					foundFiles.push_back(FindFileData.cFileName);
				}
				else
				{
					// File found: create a path to the file
					char filePath[2048] = {0};
					_snprintf(filePath, 2047, "%s%s", fullpath, FindFileData.cFileName);
					// Add it to vector of files found
					foundFiles.push_back(filePath);
				}
			}
		}
		else // Always store filename if a starting string has not been provided
		{
			if(getOnlyFileName)
			{
				foundFiles.push_back(FindFileData.cFileName);
			}
			else
			{
				// File found: create a path to the file
				char filePath[2048] = {0};
				_snprintf(filePath, 2047, "%s%s", fullpath, FindFileData.cFileName);
				// Add it to vector of files found
				foundFiles.push_back(filePath);
			}
		}
	}
	// Loop while we find more files
	while(FindNextFile(hFind, &FindFileData) != 0);

	// Release
	FindClose(hFind);

	sort(foundFiles.begin(), foundFiles.end(), StringCompare_Smart_Incr());

	foundFiles.resize( Min(nMaxItems, (int)foundFiles.size() ) );

	return true;
#else // not _MSC_VER
    DIR* directory = opendir(path.c_str());
    if(directory)
    {
        string parent(path);
        if(parent[parent.length()-1] != '/')
            parent.append("/");

        struct dirent dirEntry;
        struct dirent* res = &dirEntry;
        while((readdir_r(directory, &dirEntry, &res) == 0) && (res)) // thread-safe
            if((dirEntry.d_type == DT_REG) &&
                    (strncmp(dirEntry.d_name, startingString.c_str(), startingString.length()) == 0) &&
                    (strcmp(dirEntry.d_name, ".") != 0) &&
                    (strcmp(dirEntry.d_name, "..") != 0))
			if(getOnlyFileName)
			{
				foundFiles.push_back(dirEntry.d_name);
			}
			else
			{
                foundFiles.push_back(parent + dirEntry.d_name);
			}
        closedir(directory);

		sort(foundFiles.begin(), foundFiles.end(), StringCompare_Smart_Incr());

		foundFiles.resize( Min(nMaxItems, (int)foundFiles.size() ) );

        return true;
    }

    return false;
#endif // _MSC_VER
}


int ReLOC::GetBoundaryPoints(vtkPolyData *polydata, bool* &boundaryPointsIds)
{
	boundaryPointsIds = new bool[polydata->GetNumberOfPoints()];
	for(int po=0; po<polydata->GetNumberOfPoints(); po++)
	{
		boundaryPointsIds[po] = false;
	}
	vtkIdType* ptrCells = polydata->GetPolys()->GetPointer();
	int nVertex, ve1, ve2;
	vtkIdList* idList = vtkIdList::New();
	polydata->BuildLinks();
	//for every cell in polydata
	for(int ce=0; ce<polydata->GetNumberOfPolys(); ce++)
	{
		nVertex = *ptrCells;
		ptrCells++;
		//for every edge in polydata
		for(int ve=0; ve<nVertex; ve++)
		{
			ve1 = *(ptrCells + ve);
			ve2 = *(ptrCells + ((ve+1)%nVertex) );
			polydata->GetCellEdgeNeighbors(ce, ve1, ve2,idList);
			if(idList->GetNumberOfIds()<1)
			{
				assert(ve1 >= 0 && ve1 < polydata->GetNumberOfPoints());
				assert(ve2 >= 0 && ve2 < polydata->GetNumberOfPoints());
				boundaryPointsIds[ve1] = true;
				boundaryPointsIds[ve2] = true;
			}
		}
		ptrCells += nVertex;
	}
	idList->Delete();
	return polydata->GetNumberOfPoints();
}


void ReLOC::PopulateKdTree_NaiveSampling(vtkPolyData* cloud, const float samplingPerc, KdTree &kdTree)
{
	vtkPolyData* sampledCloud = CreatePolyData_NaiveSampling(cloud, samplingPerc);

	kdTree.SetPolyData(sampledCloud);

	sampledCloud->Delete();
}

vtkPolyData* ReLOC::CreatePolyData_NaiveSampling(vtkPolyData* cloud, const float samplingPerc)
{
	float samplingIncr = 1.0f/samplingPerc;

	int nPoints = cloud->GetNumberOfPoints();

	vtkIdType* pointIds = new vtkIdType[nPoints];
	vtkIdType* ptrpointIds = pointIds;
	vtkIdType nPointIds = 0;
	for(float ne = 0.0f; ne < (float)nPoints; ne+= samplingIncr)
	{
		*ptrpointIds = (vtkIdType)ne;
		++ptrpointIds;

		nPointIds++;
	}

	vtkPolyData* sampledCloud = CreatePolyData(cloud, pointIds, nPointIds, false, false, false);
	
	delete[] pointIds;

	return sampledCloud;
}



void ReLOC::RandomOrthogonalAxis(float* axis, float* randOrthoAxis)
{
	if(!AreEqual(axis[2], 0.0f))
	{
		randOrthoAxis[0] = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
		randOrthoAxis[1] = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
		randOrthoAxis[2] = - (axis[0]*randOrthoAxis[0] + axis[1]*randOrthoAxis[1]) / axis[2];
	}
	else if(!AreEqual(axis[1], 0.0f))
	{
		randOrthoAxis[0] = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
		randOrthoAxis[2] = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
		randOrthoAxis[1] = - (axis[0]*randOrthoAxis[0] + axis[2]*randOrthoAxis[2]) / axis[1];
	}
	else if(!AreEqual(axis[0], 0.0f))
	{
		randOrthoAxis[1] = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
		randOrthoAxis[2] = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
		randOrthoAxis[0] = - (axis[1]*randOrthoAxis[1] + axis[2]*randOrthoAxis[2]) / axis[0];
	}

	vtkMath::Normalize(randOrthoAxis);

	// check if the computed x axis is orthogonal to the normal
	assert(AreEqual(randOrthoAxis[0]*axis[0] + randOrthoAxis[1]*axis[1] + randOrthoAxis[2]*axis[2], 0.0f, 1E-6f));

}


void ReLOC::DirectedOrthogonalAxis(float* axis, float* axisOrigin, float* point, float* directedOrthoAxis)
{
	float projection[3];
	ProjectPointOnPlane(point, axisOrigin, axis, projection);
	directedOrthoAxis[0] = projection[0] - axisOrigin[0];
	directedOrthoAxis[1] = projection[1] - axisOrigin[1];
	directedOrthoAxis[2] = projection[2] - axisOrigin[2];

	vtkMath::Normalize(directedOrthoAxis);

	// check if the computed x axis is orthogonal to the normal
	assert(AreEqual((float)(directedOrthoAxis[0]*axis[0] + directedOrthoAxis[1]*axis[1] + directedOrthoAxis[2]*axis[2]), 0.0f, 2E-3f));

}

void ReLOC::ProjectPointOnPlane(float* x, float* origin, float* normal, float* xproj)
{
	float t, xo[3];

	xo[0] = x[0] - origin[0];
	xo[1] = x[1] - origin[1];
	xo[2] = x[2] - origin[2];

	t = Dot3D(normal,xo);

	xproj[0] = x[0] - t * normal[0];
	xproj[1] = x[1] - t * normal[1];
	xproj[2] = x[2] - t * normal[2];
}

float ReLOC::Dot3D(const float pt1[3], const float pt2[3])
{
    return (pt1[0]*pt2[0] + pt1[1]*pt2[1] + pt1[2]*pt2[2]);
}


void ReLOC::CalcAllViewPairs(std::vector<std::pair<int, int> > &vMatchPairs, int nElements)
{
	int nPairs = (int)( (nElements) * ((nElements-1)/2.0) );

	vMatchPairs.resize(nPairs);

	std::vector<std::pair<int, int> >::iterator itPa = vMatchPairs.begin();
	for(int i=0; i<nElements-1; i++)
	{
		for(int j=i+1; j<nElements; j++)
		{
			itPa->first = i;
			itPa->second = j;
			itPa++;
		}
	}
}


vtkTransform* ReLOC::Duplicate(vtkTransform* source)
{
	if(!source)
		return NULL;

	vtkTransform* transf = vtkTransform::New();
	transf->DeepCopy(source);
	return transf;
}

void ReLOC::DeleteTransforms(std::vector<vtkTransform*> &transforms)
{
	for(size_t tr=0; tr<transforms.size(); tr++)
	{
		if(transforms[tr])
		{
			transforms[tr]->Delete();
			transforms[tr] = NULL;
		}
	}
}


double ReLOC::ComputeMeshResolution(vtkPolyDataCollection* polyDataCollection)
{
	double meshRes = 0;

	vtkPolyData* poly;
	polyDataCollection->InitTraversal();
	int n=0;
	while( (poly = polyDataCollection->GetNextItem()) )
	{
		meshRes += ComputeMeshResolution(poly);
		n++;
	}
	meshRes /= (double)polyDataCollection->GetNumberOfItems();
	
	return meshRes;
}

double ReLOC::ComputeMeshResolution(vtkPolyData* cloud)
{
	double meshResolution = 0;
	int numdistances = 0;

	if (cloud == NULL)
		return 0.0;

	vtkCellArray* polys = cloud->GetPolys();

	if (polys == NULL)
		return 0.0;

	vtkIdType* ptrCellIds = polys->GetPointer();

	if (ptrCellIds == NULL)
		return 0.0;

	int nEdges;
	int firstVer;
	int secondVer;
	double firstPoint[3], secondPoint[3];
	for(int ce=0; ce<cloud->GetNumberOfPolys(); ce++)
	{
		nEdges = *ptrCellIds;
		ptrCellIds++;
		for(int ed=0; ed<nEdges; ed++)
		{
			firstVer = *(ptrCellIds+ed);
			secondVer = *(ptrCellIds+((ed+1)%nEdges));
			cloud->GetPoint(firstVer, firstPoint);
			cloud->GetPoint(secondVer, secondPoint);
			meshResolution += EuclideanDistance_3D(firstPoint, secondPoint );
			numdistances++;
		}
		ptrCellIds+=nEdges;
	}
	if(numdistances!=0)
	{
		meshResolution/=numdistances;
	}

	return meshResolution;
}



vtkPolyDataCollection* ReLOC::ReadMeshes(const string &meshesPath, const string meshExt)
{
	vector<string> vAbsMeshFileNames;
	FindFilesEndWith(meshesPath, meshExt, vAbsMeshFileNames, false);
	return ReadMeshes(vAbsMeshFileNames);
}

vtkPolyDataCollection* ReLOC::ReadMeshes(const vector<string> &vAbsMeshFileNames)
{
	vtkPolyDataCollection* meshes;

	//load meshes in ply format
	cout << "Loading clouds in ply format...";
	meshes = LoadPolyData(vAbsMeshFileNames);
	if(!meshes)
	{
		throw runtime_error("ERROR: Impossible to load ply files"); 
	}

	Preprocessing(meshes);


	return meshes;
}


bool ReLOC::FindFilesEndWith(const std::string & path, const std::string & endingString, std::vector<std::string> & foundFiles, bool getOnlyFileName, const int nMaxItems )
{
	foundFiles.clear();

#ifdef _MSC_VER
	// Required structs for searching for files and directories
	WIN32_FIND_DATA FindFileData;
	HANDLE hFind = INVALID_HANDLE_VALUE;

	// Build the file search string...
	char searchDir[2048] = {0};
	char fullpath[2048] = {0};

	// ...if it already is a path that ends with \ or /, add '*'...
	if(path.at(path.length() - 1) == '\\' || path.at(path.length() - 1) == '/')
	{
		_snprintf(searchDir, 2047, "%s*", path.c_str());
		_snprintf(fullpath, 2047, "%s", path.c_str()); // just copy path
	}
	// ...otherwise, add '\*' to the end of the path.
	else
	{
		_snprintf(searchDir, 2047, "%s/*", path.c_str());
		_snprintf(fullpath, 2047, "%s/", path.c_str()); // copy path and add slash (required when building filenames)
	}

	// Find the first file in the directory.
	hFind = FindFirstFile(searchDir, &FindFileData);

	// If there is no file, return
	if (hFind == INVALID_HANDLE_VALUE)
	{
		return false;
	}

	// loop
	do
	{
		// Skip ".", ".." and all directories
		if( ((FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0) || strcmp(FindFileData.cFileName, ".") == 0 || strcmp(FindFileData.cFileName, "..") == 0)
		{
			continue;
		}

		if(endingString.size() > std::string(FindFileData.cFileName).size())
		{
			continue;
		}

		// Store filename into the vector if it starts with the given string
		if(endingString.size() > 0)
		{
			if(std::string(FindFileData.cFileName).rfind(endingString) == (std::string(FindFileData.cFileName).size() - endingString.size()) )
			{
				if(getOnlyFileName)
				{
					foundFiles.push_back(FindFileData.cFileName);
				}
				else
				{
					// File found: create a path to the file
					char filePath[2048] = {0};
					_snprintf(filePath, 2047, "%s%s", fullpath, FindFileData.cFileName);
					// Add it to vector of files found
					foundFiles.push_back(filePath);
				}
			}
		}
		else // Always store filename if a starting string has not been provided
		{
			if(getOnlyFileName)
			{
				foundFiles.push_back(FindFileData.cFileName);
			}
			else
			{
				// File found: create a path to the file
				char filePath[2048] = {0};
				_snprintf(filePath, 2047, "%s%s", fullpath, FindFileData.cFileName);
				// Add it to vector of files found
				foundFiles.push_back(filePath);
			}
		}
	}
	// Loop while we find more files
	while(FindNextFile(hFind, &FindFileData) != 0);

	// Release
	FindClose(hFind);

	sort(foundFiles.begin(), foundFiles.end(), StringCompare_Smart_Incr());

	foundFiles.resize( Min(nMaxItems, (int)foundFiles.size() ) );

	return true;
#else // not _MSC_VER
    DIR* directory = opendir(path.c_str());
    if(directory)
    {
        string parent(path);
        if(parent[parent.length()-1] != '/')
            parent.append("/");

        struct dirent dirEntry;
        struct dirent* res = &dirEntry;
        while((readdir_r(directory, &dirEntry, &res) == 0) && (res)) // thread-safe
            if((dirEntry.d_type == DT_REG) &&
                    (strncmp(dirEntry.d_name+(d_namlen-endingString.size()), endingString.c_str(), endingString.length()) == 0) &&
                    (strcmp(dirEntry.d_name, ".") != 0) &&
                    (strcmp(dirEntry.d_name, "..") != 0))
			if(getOnlyFileName)
			{
				foundFiles.push_back(dirEntry.d_name);
			}
			else
			{
                foundFiles.push_back(parent + dirEntry.d_name);
			}
        closedir(directory);

		sort(foundFiles.begin(), foundFiles.end(), StringCompare_Smart_Incr());

		foundFiles.resize( Min(nMaxItems, (int)foundFiles.size() ) );

        return true;
    }

    return false;
#endif // _MSC_VER
}



vtkPolyDataCollection* ReLOC::LoadPolyData(const std::vector<std::string> &filenames)
{
	vtkPolyDataCollection* polyDataCollection = vtkPolyDataCollection::New();
	vtkPolyData* poly;
	for(std::vector<std::string>::const_iterator iter = filenames.begin(); iter!= filenames.end(); ++iter)
	{
		poly = LoadPolyData(iter->data());
		if(poly == NULL)
		{
			polyDataCollection->Delete();
			return NULL;
		}
		polyDataCollection->AddItem(poly);
		poly->Delete();
	}

	return polyDataCollection;
}



vtkPolyData* ReLOC::LoadPolyData(const std::string &absFilename)
{
	std::string strFileName = absFilename;
	#ifndef BOOST_IS_NOT_INCLUDED
		boost::filesystem::path filepath(filename);
    #if defined (BOOST_FILESYSTEM_VERSION) && BOOST_FILESYSTEM_VERSION >= 3
      std::string ext = filepath.extension().string();
    #else
      std::string ext = filepath.extension();
    #endif
	#else
		std::string ext = strFileName.substr(strFileName.find_last_of("."));
	#endif
	std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

	if (ext == ".ply")
		return LoadPly(strFileName.c_str());
	//else if (ext == ".off")
	//	return LoadOff(strFileName.c_str());
	//else if (ext == ".obj")
	//	return LoadOBJ(strFileName.c_str());
	//else if (ext == ".stl")
	//	return LoadSTL(strFileName.c_str());
	//else if (ext == ".vtk")
	//	return LoadVtk(strFileName.c_str());
	//else if (ext == ".vert")
	//	return LoadVertTri(strFileName.substr(0, strFileName.length() - ext.length()).c_str());
	//else if (ext == ".tri")
	//	return LoadVertTri(strFileName.substr(0, strFileName.length() - ext.length()).c_str());
	//else if (ext == ".pcd")
	//	return LoadPcd(strFileName.c_str());
	//else if (ext == ".txt")
	//	return LoadSimpleTxt(strFileName.c_str());
	//else if (ext == ".xyz")
	//	return LoadXYZ(strFileName.c_str());
	else
		throw runtime_error("ERROR Type not supported");

	return NULL;
}


vtkPolyData* ReLOC::LoadPly(const char* filename)
{
	vtkPLYReader  *plyReader = vtkPLYReader::New();
	plyReader->SetFileName(filename);
	plyReader->Update();
	int errCod = plyReader->GetErrorCode();
	if(errCod != vtkErrorCode::NoError)
	{
		const char* errString = vtkErrorCode::GetStringFromErrorCode(errCod);
		printf("\nERROR LoadPly: %s\n", errString);
		plyReader->Delete();
		return NULL;
	}

	vtkPolyData* polyData = vtkPolyData::New();
	polyData->ShallowCopy(plyReader->GetOutput());
	plyReader->Delete();
	return polyData;
}


void ReLOC::Preprocessing(vtkPolyDataCollection* &polyDataCollection)
{
	cout << "Computing normals...";
	CalcMeshNormals(polyDataCollection);
	CleanPolyData(polyDataCollection, true, 0.0);
	cout << "OK" << endl;
}


void ReLOC::CalcMeshNormals(vtkPolyDataCollection *polyDataCollection, const double dFeatureAngle)
{
	vtkPolyData* iter=NULL;
	polyDataCollection->InitTraversal();
	while((iter = polyDataCollection->GetNextItem()))
	{
		CalcMeshNormals(iter, dFeatureAngle);
	}
}

void ReLOC::CalcMeshNormals(vtkPolyData *polyData, const double dFeatureAngle)
{
	if(polyData == NULL)
	{
		return;
	}
	vtkPolyDataNormals* polydataNormals = vtkPolyDataNormals::New();
	if(dFeatureAngle < 0)
	{
		polydataNormals->SplittingOff();
	}
	else
	{
		polydataNormals->SplittingOn();
		polydataNormals->SetFeatureAngle(dFeatureAngle);
	}
#if VTK_MAJOR_VERSION <= 5
	polydataNormals->SetInput(polyData);
#else
	polydataNormals->SetInputData(polyData);
#endif
	polydataNormals->SetInput(polyData);
	polydataNormals->ComputePointNormalsOn();
	polydataNormals->AutoOrientNormalsOff();
	polydataNormals->ConsistencyOff();
	polydataNormals->ComputeCellNormalsOff();

	polydataNormals->Update();

	polyData->ShallowCopy( polydataNormals->GetOutput() );
	polydataNormals->Delete();
}


void ReLOC::CleanPolyData(vtkPolyDataCollection *polyDataCollection, const bool mergePoints, const double tolerance, bool removeNotPolysCells)
{
	vtkPolyData* iter=NULL;
	polyDataCollection->InitTraversal();
	while((iter = polyDataCollection->GetNextItem()))
	{
		CleanPolyData(iter, mergePoints, tolerance, removeNotPolysCells);
	}
}

void ReLOC::CleanPolyData(vtkPolyData *polyData, const bool mergePoints, const double tolerance, bool removeNotPolysCells)
{
	if(polyData == NULL)
	{
		return;
	}

	if(removeNotPolysCells)
	{
		RemoveNotPolysCells(polyData);
	}

	vtkCleanPolyData* cleanPolyData = vtkCleanPolyData::New();
#if VTK_MAJOR_VERSION <= 5
	cleanPolyData->SetInputConnection(polyData->GetProducerPort());
#else
	cleanPolyData->SetInputData(polyData);
#endif
	cleanPolyData->SetAbsoluteTolerance(tolerance);
	cleanPolyData->SetPointMerging(mergePoints);
	cleanPolyData->ToleranceIsAbsoluteOn();
	//cleanPolyData->ConvertLinesToPointsOff();
	//cleanPolyData->ConvertPolysToLinesOff();
	//cleanPolyData->ConvertStripsToPolysOff();
	//cleanPolyData->PieceInvariantOff();

	cleanPolyData->Update();
	polyData->ShallowCopy( cleanPolyData->GetOutput() );
	cleanPolyData->Delete();

	if(removeNotPolysCells)
	{
		RemoveNotPolysCells(polyData);
	}
}

void ReLOC::RemoveNotPolysCells(vtkPolyData *polyData)
{
	polyData->GetLines()->Reset();
	polyData->GetStrips()->Reset();
	polyData->GetVerts()->Reset();
}

