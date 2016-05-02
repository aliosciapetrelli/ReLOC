#ifndef RELOC_UTILS_H
#define RELOC_UTILS_H

#define BOOST_IS_NOT_INCLUDED

#include <vector>
#include <sstream>
#include "vtkPolyData.h"
#include "vtkPolyDataCollection.h"
#include "vtkTransform.h"
#include "vtkAbstractPointLocator.h"

#define BOOST_IS_NOT_INCLUDED

#ifdef _MSC_VER
	#ifndef _CRT_SECURE_NO_WARNINGS
		#define _CRT_SECURE_NO_WARNINGS
	#endif
	#define NOMINMAX
	#include "windows.h"
#endif


namespace ReLOC
{

	/*!
	* \brief Compute minimum.
	*
	* \param a,b instances to compare.
	*/
	template <typename Type> inline Type Min(Type a, Type b) { return (a <= b)? a : b; };

	/*!
	* \brief Compute maximum.
	*
	* \param a,b instances to compare.
	*/
	template <typename Type> inline Type Max(Type a, Type b) { return (a >= b)? a : b; };


	float* GetPolyDataPointsPointer(vtkPolyData* polyData, const unsigned int iPoint=0);
	float* GetPolyDataPointNormalsPointer(vtkPolyData* polyData, const unsigned int iNormal=0);
	vtkIdType* GetPolyDataTrianglesPointer(vtkPolyData* polyData, const unsigned int iTriangle=0);
	unsigned char* GetPolyDataColorsPointer(vtkPolyData* polyData, const unsigned int iColor=0);
	vtkIdType GetPolyDataNumberOfPoints(vtkPolyData* polyData);

	float* AllocateNormals(vtkPolyData *polyData);
	float* AllocatePoints(vtkPolyData *polyData, const vtkIdType nPoints);
	unsigned char* AllocateColors(vtkPolyData *polyData, unsigned char* color = NULL, bool* mask = NULL);
	vtkIdType* AllocateTriangles(vtkPolyData *polyData, const vtkIdType nTriangles);


	vtkPolyData* CreatePolyData(const float* points, const vtkIdType nPoints, const unsigned char* colors = NULL, const vtkIdType* triangles = NULL, const vtkIdType nTriangles = 0);

	vtkPolyData* CreatePolyData(vtkPolyData* polyData, vtkIdType* pointIds, vtkIdType nPointIds, const bool insertNormals = true, const bool insertColors = true, const bool insertTriangles = true);

	

	/*!
	* \brief Compute squared euclidean distance between two 3D points.
	*
	* \param point1 3D point
	* \param point2 3D point
	* \return squared euclidean distance
	*/
	inline float EuclideanDistance2_3D(float* point1, float* point2){return (point1[0]-point2[0])*(point1[0]-point2[0]) + (point1[1]-point2[1])*(point1[1]-point2[1]) + (point1[2]-point2[2])*(point1[2]-point2[2]);  };

	/*!
	* \brief Compute squared euclidean distance between two 3D points.
	*
	* \param point1 3D point
	* \param point2 3D point
	* \return squared euclidean distance
	*/
	inline double EuclideanDistance2_3D(double* point1, double* point2){return (point1[0]-point2[0])*(point1[0]-point2[0]) + (point1[1]-point2[1])*(point1[1]-point2[1]) + (point1[2]-point2[2])*(point1[2]-point2[2]);  };


	/*!
	* \brief Compute euclidean distance between two 3D points.
	*
	* \param point1 3D point
	* \param point2 3D point
	* \return euclidean distance
	*/
	inline float EuclideanDistance_3D(float* point1, float* point2){ return sqrt(EuclideanDistance2_3D(point1, point2)); };

	/*!
	* \brief Compute euclidean distance between two 3D points.
	*
	* \param point1 3D point
	* \param point2 3D point
	* \return euclidean distance
	*/
	inline double EuclideanDistance_3D(double* point1, double* point2){ return sqrt(EuclideanDistance2_3D(point1, point2)); };

	double calcDistance2BetweenPoints(float* points1, float* points2, int nPoints);

	template<typename Value> inline Value StdDev(Value* vValues, const Value average, int nValues, Value init)
	{
		Value strDev = init;
		for(int i=0; i<nValues; i++)
		{
			strDev+=(vValues[i]*vValues[i]);
		}
		if(nValues)
		{
			strDev/=nValues;
		}
		strDev -= (average*average);
		strDev = sqrt(strDev);
		return strDev;
	}

	template<typename Value> inline Value StdDev(const std::vector<Value> &vValues, const Value average, const Value &init)
	{
		return StdDev((Value*)(&vValues[0]), average, (int)vValues.size(), init);
	}

	void GetCentroid(float* points, int nPoints, float centroid[3]);
	void GetCentroid(vtkPolyData *points, float centroid[3]);
	void GetCentroids_cMeshes(vtkPolyDataCollection* coll, std::vector<float> &vCentroids);

	void Translate(float* points, int nPoints, const float point[3]);

	void GetBBox_StdDev_Cloud(vtkPolyData* cloud, double* bbox, const float* centroid, const float sigmaFactor = 3.0f);
	void GetBBox_StdDev_Clouds(vtkPolyDataCollection* cMeshes, std::vector<double> & vBBoxes, const std::vector<float> &vCentroids, const float sigmaFactor = 3.0f);

	const float zeroFloatEps8 = 1E-8f;
	inline bool IsZero(float val, float zeroFloatEps = zeroFloatEps8){return (abs(val)<zeroFloatEps);};

	/*!
	* \brief Check if val1 and val2 are equals.
	*
	* \param val1 first number to check.
	* \param val2 second number to check.
	* \return true if val1 is equal to val2, false otherwise.
	*/
	inline bool AreEqual(float val1, float val2, float zeroFloatEps = zeroFloatEps8){return (abs(val1 - val2)<zeroFloatEps);};

	void ApplyTransform(float* points, int nPoints, const double* matrTransf);


	class KdTree
	{
	protected:
		vtkAbstractPointLocator*	m_kdTree;
		vtkIdList*				m_pointsList;
		vtkPolyData*			m_polyData;

	public:
		KdTree(const bool useTrueKdTree = false);
		KdTree(vtkPolyData* polyData, const bool useTrueKdTree = false);
		~KdTree();

		void SetPolyData(vtkPolyData* polyData);
		vtkPolyData* GetPolyData(){return m_polyData;};

		vtkAbstractPointLocator* GetKdTree(){return m_kdTree;}
		vtkIdList* GetFoundPoints(){return m_pointsList;}

		//FindPointsWithinRadius
		vtkIdList* FindPointsWithinRadius(const float* const point, float radius);
		vtkIdList* FindPointsWithinRadius(const double* point, double radius);
		/*!
		* \param pointIndex index of m_polydata point.
		*/
		vtkIdList* FindPointsWithinRadius(const int pointIndex, float radius);

		//FindNearestPoint
		int FindNearestPoint(double* point);
		int FindNearestPoint(float* point);
		/*!
		* \param pointIndex index of m_polydata point.
		*/
		int FindNearestPoint(const int pointIndex);

		//FindNearestPoint, return distance
		double FindNearestPoint(double* point, int &nearestPointId);
		double FindNearestPoint(float* point, int &nearestPointId);
		/*!
		* \param pointIndex index of m_polydata point.
		*/
		double FindNearestPoint(const int pointIndex, int &nearestPointId);

		//FindNearestPointWithinRadius
		int FindNearestPointWithinRadius(double* point, double radius, double & dist);
		/*!
		* \param point
		* \param radius
		* \param dist output parameter
		* \return the index of the nearest point in the mesh of the polydata, -1 if not found
		*/
		int FindNearestPointWithinRadius(float* point, float radius, double & dist);
		/*!
		* \param pointIndex index of m_polydata point.
		*/
		int FindNearestPointWithinRadius(const int pointIndex, float radius, double & dist);

		//FindNearestNPoints
		vtkIdList* FindNearestNPoints(int n, const double* point);
		vtkIdList* FindNearestNPoints(int n, const float* point);
		vtkIdList* FindNearestNPoints(int n, const int pointIndex);
	};


	std::string ChangeExtension(const std::string &filename, const std::string &newExt);

	bool ExistsDir(const std::string &directory);

	void CreateDir(const std::string &directory);
	void CreateFullPath(const std::string &directory);
	std::string GetPathDir(const std::string &filename);

	bool FindFilesStartWith(const std::string & path, const std::string & startingString, std::vector<std::string> & foundFiles, bool getOnlyFileName = false, const int nMaxItems = std::numeric_limits<int>::max());


	template<typename T> T LexicalCast(const std::string& s)
	{
		std::stringstream ss(s);

		T result;
		if ((ss >> result).fail() || !(ss >> std::ws).eof())
		{
			//throw std::bad_cast();
			cout << "ERROR:Impossible to cast " << s;
			getchar();
			exit(-1);
		}

		return result;
	}

	//logical sorting (e.g. Windows explorer)
	class StringCompare_Smart_Incr
	{
	public:
		inline bool operator() (const std::string& a, const std::string& b) const
		{
			unsigned posStr = 0;
			while( (posStr < a.size() ) && (posStr < b.size()) )
			{
				unsigned tkn_idx_a = (unsigned int) a.find_first_of("0123456789", posStr);
				unsigned tkn_idx_b = (unsigned int) b.find_first_of("0123456789", posStr);
				std::string suba = a.substr(posStr, tkn_idx_a - posStr );
				std::string subb = b.substr(posStr, tkn_idx_b - posStr );
				if(suba == subb)
				{
					//same substring

					if(tkn_idx_a == a.size())
					{
						//end of a and at least of b
						return true;
					}

					if(tkn_idx_a == a.size())
					{
						//end of b but not of a
						return false;
					}

					unsigned numberEnd_a = (unsigned int) a.find_first_not_of("0123456789", tkn_idx_a+1);
					unsigned numberEnd_b = (unsigned int) b.find_first_not_of("0123456789", tkn_idx_b+1);
					//check number
					long long number_a = LexicalCast<long long>(a.substr(tkn_idx_a, numberEnd_a - tkn_idx_a));
					long long number_b = LexicalCast<long long>(b.substr(tkn_idx_b, numberEnd_b - tkn_idx_b));
					//long number_a = std::atol(a.substr(tkn_idx_a).c_str());
					//long number_b = std::atol(b.substr(tkn_idx_b).c_str());
					if(number_a != number_b)
					{
						return (number_a < number_b);
					}
				}
				else
				{
					//different substring
					return (suba < subb);
				}
				posStr = (unsigned int) a.find_first_not_of("0123456789", tkn_idx_a + 1);
			}

			return ( a.size() < b.size() );
		}
	};


	int GetBoundaryPoints(vtkPolyData *polydata, bool* &boundaryPointsIds);



	void PopulateKdTree_NaiveSampling(vtkPolyData* cloud, const float samplingPerc, KdTree &kdTree);
	vtkPolyData* CreatePolyData_NaiveSampling(vtkPolyData* cloud, const float samplingPerc);


	/*!
	* \brief Given an axis (with origin axisOrigin), return the orthogonal axis directed to point.
	*
	* axis must be normalized.
	*/
	void DirectedOrthogonalAxis(float* axis, float* axisOrigin, float* point, float* directedOrthoAxis);

	void RandomOrthogonalAxis(float* axis, float* randOrthoAxis);

	/*!
	* \brief Given a plane (origin and normal) and a point x, return the projection xproj of x on plane
	*
	* Equivalent to vtkPlane::ProjectPoint()
	*/
	void ProjectPointOnPlane(float* x, float* origin, float* normal, float* xproj);


	/*!
	* \brief Compute dot product between two 3D vectors
	*
	* \param pt1 3D vector
	* \param pt2 3D vector
	* \return dot product between pt1 and pt2
	*/
	float Dot3D(const float pt1[3], const float pt2[3]);


	//sort a vector of indeces w.r.t. values in vector f
	class SortIds_Wrt_FloatDecr
	{
	public:
		SortIds_Wrt_FloatDecr(std::vector<float> & f) : _f(f) {}
		bool operator()(int i, int j) {
			return _f[i] > _f[j];
		}
	private:
		std::vector<float> & _f;
	};


	void CalcAllViewPairs(std::vector<std::pair<int, int> > &vMatchPairs, int nElements);


	template<typename T> void ExtendBBox(T* bounds, const T factor_X, const T factor_Y, const T factor_Z)
	{
		T size_X = bounds[1] - bounds[0];
		T size_Y = bounds[3] - bounds[2];
		T size_Z = bounds[5] - bounds[4];

		T ext2_X = (size_X * factor_X - size_X)/2.0;
		T ext2_Y = (size_Y * factor_Y - size_Y)/2.0;
		T ext2_Z = (size_Z * factor_Z - size_Z)/2.0;

		bounds[0] = bounds[0] - ext2_X;
		bounds[1] = bounds[1] + ext2_X;
		bounds[2] = bounds[2] - ext2_Y;
		bounds[3] = bounds[3] + ext2_Y;
		bounds[4] = bounds[4] - ext2_Z;
		bounds[5] = bounds[5] + ext2_Z;

	}

	vtkTransform* Duplicate(vtkTransform* source);

	void DeleteTransforms(std::vector<vtkTransform*> &transforms);

	double ComputeMeshResolution(vtkPolyData* cloud);
	double ComputeMeshResolution(vtkPolyDataCollection* polyDataCollection);



	vtkPolyDataCollection* ReadMeshes(const std::vector<std::string> &vAbsMeshFileNames );
	vtkPolyDataCollection* ReadMeshes(const std::string &meshesPath, const std::string meshExt = "ply" );

	bool FindFilesEndWith(const std::string & path, const std::string & endingString, std::vector<std::string> & foundFiles, bool getOnlyFileName = false, const int nMaxItems = std::numeric_limits<int>::max());

	vtkPolyData* LoadPolyData(const std::string &absFilename);
	vtkPolyDataCollection* LoadPolyData(const std::vector<std::string> &filenames);
	vtkPolyData* LoadPly(const char* filename);

	void Preprocessing(vtkPolyDataCollection* &polyDataCollection);

	void CalcMeshNormals(vtkPolyData *polyData, const double dFeatureAngle = -1.0);
	void CalcMeshNormals(vtkPolyDataCollection *polyDataCollection, const double dFeatureAngle = -1.0);

	void RemoveNotPolysCells(vtkPolyData *polyData);

	void CleanPolyData(vtkPolyData *polyData, const bool mergePoints = true, const double tolerance = 0.0, bool removeNotPolysCells = true);
	void CleanPolyData(vtkPolyDataCollection *polyDataCollection, const bool mergePoints = true, const double tolerance = 0.0, bool removeNotPolysCells = true);
}

#endif



