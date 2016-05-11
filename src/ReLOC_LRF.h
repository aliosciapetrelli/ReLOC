#ifndef RELOC_LRF_H
#define RELOC_LRF_H

#include "ReLOC_utils.h"

#include <stdint.h>

#include <string>
#include <vector>
#include "vtkPolyData.h"
#include "vtkPolyDataCollection.h"


namespace ReLOC
{

	//************************ DETECTION *******************************

	struct Feature3D
	{
		float x,y,z;
		int i,j,k;	//for structured data
		float scale;
		int index;
		float rf[9];
		float score;

		Feature3D();
		Feature3D(const Feature3D & other);

		Feature3D& operator=(const Feature3D & other);
		bool operator==(const Feature3D & other) const;
		bool operator!=(const Feature3D & other) const;

	};


	class I3DDetector
	{
	public:
		I3DDetector();
		virtual ~I3DDetector();
		virtual int extract(vtkPolyData* cloud, Feature3D* & feat);
		virtual void extract(vtkPolyDataCollection* polyDataCollection);
		int getNumFeatures() { return m_numFeat; } ;
		void setNumFeatures(int numFeatures) { m_numFeat = numFeatures; } ;
		Feature3D* getFeatures() { return m_feat; };
		void setFeatures(Feature3D* feat) { m_feat = feat; };
		void setFeatures(Feature3D* feat, const int numFeatures, const bool copyFeats = false);
		std::vector<int> & getNumFeaturesVect() { return m_vNumFeat; } ;
		std::vector<Feature3D*> & getFeaturesVect() { return m_vFeat; };
		void setNumFeaturesVect(std::vector<int> &vNumFeat){m_vNumFeat = vNumFeat;};
		void setFeaturesVect(std::vector<Feature3D*> &vFeat){m_vFeat = vFeat;};
		virtual bool save(const char *path, const bool binaryFile = false);
		virtual bool save(const char *path, const char *fileTemplate, const bool binaryFile = false);
		virtual bool save(std::vector<std::string> &vAbsFileNames, const bool binaryFile = false);
		virtual bool save(const char *path, std::vector<std::string> &vFileNames, const bool binaryFile = false);
		virtual int load(const char *path, Feature3D* &feat, const bool binaryFile = false);
		virtual void load(const char *path, const char *fileTemplate, const bool binaryFile = false);
		virtual void load(std::vector<std::string> &vAbsFileNames, const bool binaryFile = false);
		virtual void load(const char *path, std::vector<std::string> &vFileNames, const bool binaryFile = false);
		static bool save(const char *path, const int numFeat, const Feature3D* feat, const bool binaryFile = false);

		static int compute3DIndex(vtkPolyData* poly, KdTree & kdt, Feature3D & feat, float meshResolution=-1.0f);

		virtual void resizeVect(const unsigned int newLength);

	protected:
		double m_borderDistance;
		virtual void extractImpl(vtkPolyData* cloud, Feature3D* & feat, int & numFeatures, int & featCapacity) = 0;
		static void updateFeatCapacity(Feature3D * & feat, int & capacity, const int nFeat);

	private:
		Feature3D* m_feat;
		int m_numFeat;
		int m_featCapacity;
		std::vector<Feature3D*> m_vFeat;
		std::vector<int> m_vNumFeat;
		std::vector<int> m_vFeatCapacity;


		static int load(const char *path, Feature3D* &feat, int &numFeat, int &featCapacity, const bool binaryFile = false );

	};

	/// ***********************  Random3DDetector  ***********************

	struct RandomDetectorParams
	{
		int			nFeatures;
		int			numMinNeighbors;
		bool		random;
		double		radius;
		double		borderDistance;			///< min distance from the border
		int64_t		seed;

		RandomDetectorParams()	// default values
		{
			nFeatures = 2000;
			random = false;
			numMinNeighbors = 5;
			radius = 0.0;
			borderDistance = -1.0;
			seed = 0xFFFFFFFF;
		};
	};


	class Random3DDetector : public I3DDetector
	{
	public:
		//Random3DDetector(int nPoints, bool random = false);
		Random3DDetector(int nPoints, bool random = false, double radius = -1.0,  int minNeigh = 0, double borderDistance = -1.0, int64_t seed = 0xFFFFFFFF); //computes "nPoints" randomly that have at least "minNeigh" neighboring points in a sphere of radius "radius" and that are not closer to any border than "borderDistance"
		Random3DDetector(RandomDetectorParams params); //computes "nPoints" randomly that have at least "minNeigh" neighboring points in a sphere of radius "radius" and that are not closer to any border than "borderDistance"

		virtual ~Random3DDetector();

		virtual void extractImpl(vtkPolyData* cloud, Feature3D* & feat, int & numFeatures, int & descCapacity);

	private:
		int64_t m_seed;
		int m_minNeigh;
		double m_radius;
		int m_requestedNumFeat;
	};



	/// ***********************  FlatPoints3DDetector  ***********************



	struct FlatPoints3DDetectorParams
	{
		int				numMinNeighbors;
		bool			random;
		
		/** \brief Rf */
		double			radius;

		/** \brief R1search */
		double			partitionRadius;

		/** \brief R2search */
		double			partitionRadius_Hierarchical;

		/** \brief Rdiscard */
		double			excludedPointsRadius;

		/** \brief T1search */
		float			coverageAreaPercentage;

		/** \brief T2search */
		float			coverageAreaPercentage_Hierarchical;

		int64_t			seed;	

		FlatPoints3DDetectorParams()	// default values
		{
			random = false;
			numMinNeighbors = 5;
			radius = 3;   //5 in the paper
			partitionRadius = 2;
			excludedPointsRadius = 2;
			partitionRadius_Hierarchical = 20;
			coverageAreaPercentage = 0.9f;
			coverageAreaPercentage_Hierarchical = 0.9f;
			seed = 0xFFFFFFFF;	
		};
	};


	/** \brief Flat point detector proposed in:
	*  Petrelli A., Di Stefano L., "Pairwise registration by local orientation cues", Computer Graphics Forum, 2015
	*
	* \author Alioscia Petrelli
	*/
	class FlatPoints3DDetector : public I3DDetector
	{			
	public:
		FlatPoints3DDetector(bool random = false, double radius = 3, double partitionRadius = 2, double excludedPointsRadius = 2, double partitionRadius_Hierarchical = 20, float coverageAreaPercentage = 0.9f, float coverageAreaPercentage_Hierarchical = 0.9f, int minNeigh = 5, int64_t seed = 0xFFFFFFFF); 
		FlatPoints3DDetector(FlatPoints3DDetectorParams &params); 

		void ApplyFlatPointsDetector(vtkPolyData* cloud, std::vector<float> & vFeaturesScore, float exRadius, float partRadius, float covAreaP, std::vector<vtkIdType> & vFeaturePointIndexes);

		virtual ~FlatPoints3DDetector();

		virtual void extractImpl(vtkPolyData* cloud, Feature3D* & feat, int & numFeatures, int & descCapacity);

	private:
		int64_t			m_seed;
		int				m_minNeigh;

		/** \brief Rf */
		double			m_radius;

		/** \brief R1search */
		double			m_partitionRadius;
		
		/** \brief R2search */
		double			m_partitionRadius_Hierarchical;

		/** \brief Rdiscard */
		double			m_excludedPointsRadius;
		
		/** \brief T1search */
		float			m_coverageAreaPercentage;

		/** \brief T2search */
		float			m_coverageAreaPercentage_Hierarchical;
	};



	//*******************  DESCRIPTION *********************


	struct Desc
	{
	public:
		Desc();
		Desc(const char *path);
		// It does a shallow copy if and only if (!other.isOwner() && !forceDeepCopy)
		// Call shallowCopy() if you want a shallow copy when (other.isOwner()==TRUE)
		Desc(const Desc& other, bool forceDeepCopy);
		// It calls *this = other
		Desc(const Desc& other);
		~Desc();

		// It does a deep copy if and only if (other.isOwner()==TRUE)
		// Call shallowCopy() if you want a shallow copy when (other.isOwner()==TRUE)
		Desc & operator =(const Desc & other);
		Desc & operator +=(const Desc & other);

		double** updateDouble(int size, int length);
		float** updateFloat(int size, int length);
		double** getDouble();
		double const * const * getDoubleReadOnly() const;
		float** getFloat();
		float const * const * getFloatReadOnly() const;
		
		// Shallow copy of first/last n descriptors
		Desc head(unsigned n) const { return segment(0, n); }
		Desc tail(unsigned n) const { return segment(m_size-n, n); }
		// Shallow copy of n descriptors from start;
		Desc segment(unsigned start, unsigned n) const;

		Desc shallowCopy() const { return segment(0, m_size); }
		// Always force a deep copy of other
		// eg: d1 = d2.segment(a,b) --> d1 is a shallow copy
		//     d1.deepCopyFrom(d2.segment(a,b)) --> d1 is a deep copy
		Desc& deepCopyFrom(const Desc& other);

		bool save(const char *path, const bool binaryFile = false, const bool append = false);
		static bool save(const char *path, const char *fileTemplate, std::vector<Desc> &vDesc, const bool binaryFile = false);
		static bool save(std::vector<std::string> &vAbsFileNames, std::vector<Desc> &vDesc, const bool binaryFile = false);
		static bool save(const char *path, std::vector<std::string> &vFileNames, std::vector<Desc> &vDesc, const bool binaryFile = false);

		bool load(const char *path, const bool binaryFile = false);
		static bool load(const char *path, const char *fileTemplate, std::vector<Desc> &vDesc, const bool binaryFile = false);
		static bool load(std::vector<std::string> &vAbsFileNames, std::vector<Desc> &vDesc, const bool binaryFile = false);
		static bool load(const char *path, std::vector<std::string> &vFileNames, std::vector<Desc> &vDesc, const bool binaryFile = false);

		int getLength() const { return m_length; }
		int getNumDesc() const { return m_size; }
		void setNumDesc(int new_size) { m_size = new_size; }
		bool isDouble() const { return m_doubleInit; }
		bool isBitset() const { return m_bitsetInit; }
		bool isFloat() const { return m_floatInit; }
		bool isOwner() const { return m_ptrOwner; }

		void release();

		void sort(const std::vector<int> &vId_sorting);

	private:
		bool m_ptrOwner; // if TRUE the data pointer is deleted on destruction

		// supported descriptors types
		bool m_doubleInit, m_bitsetInit, m_floatInit;
		int m_size, m_length, m_capacity;

		double** m_doubleDesc;
		float** m_floatDesc;

		double** updateDoubleInternal(int size, int length);
		float** updateFloatInternal(int size, int length);

		void init();
	};



	class I3DDescriptor
	{
	public:
		I3DDescriptor();
		virtual ~I3DDescriptor();
		virtual void describe(vtkPolyData* cloud, Feature3D* & feat, Desc & desc, int nFeat);
		virtual void describe(vtkPolyDataCollection* polyDataCollection, std::vector<Feature3D*> & vFeat, std::vector<Desc> & vDesc, std::vector<int> vNFeat);

	protected:
		int m_descLength;

		virtual void describeImpl(vtkPolyData* cloud, Feature3D* & feat, Desc & desc, int nFeat) = 0;
	};



	///  ***********************  Flare ***********************

	/* \brief FLARE method proposed in Petrelli A., Di Stefano L., "A repeatable and efficient canonical reference for surface matching", 3DimPvt 2012.
	*/
	struct FLAREparams
	{
		/* \brief Select between the standard FLARE ("FLARE") and a faster version that performs a subsampling of the support used for the computation of X axis ("FLARE_Sampled_T"). */
		std::string	selector;

		/* \brief radius (Rz) in unit of mesh resolutions of the support used for the computation of Z axis. */
		float	radiusInmeshRes_normal;

		/* \brief radius (Rx) in unit of mesh resolutions of the support used for the computation of X axis. */
		float	radiusInmeshRes_tangent;

		int		minNumNeighbours_normal;
		int		minNumNeighbours_tangent;

		/* \brief Threshold that defines the thickness of the support crown used for the computation of X axis. */
		float	crownThresh;

		float samplingPerc_normal;

		/* \brief Percentage of points actually used for the computation of X axis. */
		float samplingPerc_tangent;

		//for plane fitting computation
		float *covM[3];
		float *evect[3];

		void Init()
		{
			selector = "FLARE";				//"FLARE_Sampled_T"
			radiusInmeshRes_normal = 5.0f;
			radiusInmeshRes_tangent = 50.0f;
			minNumNeighbours_normal = 6;
			minNumNeighbours_tangent = 6;

			crownThresh = 0.85f;

			samplingPerc_normal = 1.0f;
			samplingPerc_tangent = 0.2f;

			covM[0] = new float[3];
			covM[1] = new float[3];
			covM[2] = new float[3];
			evect[0] = new float[3];
			evect[1] = new float[3];
			evect[2] = new float[3];

		}

		FLAREparams& operator=(const FLAREparams &flareParams)
		{
			selector = flareParams.selector;
			radiusInmeshRes_normal = flareParams.radiusInmeshRes_normal;
			radiusInmeshRes_tangent = flareParams.radiusInmeshRes_tangent;
			minNumNeighbours_normal = flareParams.minNumNeighbours_normal;
			minNumNeighbours_tangent = flareParams.minNumNeighbours_tangent;

			crownThresh = flareParams.crownThresh;

			samplingPerc_normal = flareParams.samplingPerc_normal;
			samplingPerc_tangent = flareParams.samplingPerc_tangent;

			return *this;
		}

		FLAREparams()
		{
			Init();
		}

		FLAREparams(const FLAREparams &flareParams)
		{
			Init();

			*this = flareParams;
		}

		std::string GetNickname()
		{
			std::stringstream name;

			name << selector;
			if(radiusInmeshRes_normal != 5.0f)
			{
				name << "_Rn-" << radiusInmeshRes_normal;
			}
			name << "_Rt-" << radiusInmeshRes_tangent;
			if(minNumNeighbours_normal != 6)
			{
				name << "_nNeign-" << minNumNeighbours_normal;
			}
			if(minNumNeighbours_tangent != 6)
			{
				name << "_nNeigt-" << minNumNeighbours_tangent;
			}

			if(crownThresh != 0.85f)
			{
				name << "_cT-" << crownThresh;
			}

			if(selector != "FLARE")
			{
				if(samplingPerc_normal != 1.0f)
				{
					name << "_sPn-" << samplingPerc_normal;
				}
				if(samplingPerc_tangent != 0.2f)
				{
					name << "_sPt-" << samplingPerc_tangent;
				}
			}

			return name.str();
		}

		~FLAREparams()
		{
			for (int n = 0; n < 3; n++)
			{
				delete [] evect[n];
				evect[n] = 0;
			}

			for (int n = 0; n < 3; n++)
			{
				delete [] covM[n];
				covM[n] = 0;
			}
		}
	};


	/*!
	* \brief FLARE method proposed in Petrelli A., Di Stefano L., "A repeatable and efficient canonical reference for surface matching", 3DimPvt 2012.
	* 
	* \param cloud input mesh (Point normals must be computed before the method is applied)
	* \param kdTree kd-tree indexing cloud points
	* \param index idx oof the feture point for which FLARE is computed. 
	* \param rfc 9-element array representing X,Y and Z axis of the computed LRF.
	* \param params FLARE params
	* \param meshRes mesh resolution of cloud 
	*/
	float getLocalRF_FLARE(vtkPolyData *cloud, KdTree &kdTree, const int index, float *rfc, FLAREparams &params, const float meshRes);


	/*!
	* \brief Fast version of FLARE that considers sampled mesh for X axis estimation.
	* 
	* \param cloud input mesh (Point normals must be computed before the method is applied)
	* \param kdTree kd-tree indexing cloud points
	* \param kdTree_sampled kd-tree indexing a sampled version of the cloud 
	* \param index idx oof the feture point for which FLARE is computed. 
	* \param rfc 9-element array representing X,Y and Z axis of the computed LRF.
	* \param params FLARE params
	* \param meshRes mesh resolution of cloud 
	*/
	float getLocalRF_FLARE_Sampled_T(vtkPolyData *cloud, KdTree &kdTree, KdTree &kdTree_sampled, const int index, float *rfc, FLAREparams &params, const float meshRes);


	///  ***********************  LRFdescriptor_1  ***********************


	/*!
	* \brief Parameters of descriptor based on FLARE method as described in Petrelli A., Di Stefano L., "Pairwise registration by local orientation cues", Computer Graphics Forum, 2015
	*/
	struct LRFDescriptor_1_Params
	{
		FLAREparams flareParams;

		/* \brief mesh resolution of cloud on which the descriptor is computed*/ 
		double	meshRes;		

		LRFDescriptor_1_Params()
		{
			meshRes = 1.0;
		};
	};


	/*!
	* \brief Descriptor based on FLARE method as described in Petrelli A., Di Stefano L., "Pairwise registration by local orientation cues", Computer Graphics Forum, 2015
	*/
	class LRFDescriptor_1 : public I3DDescriptor			
	{			
	public:				
		LRFDescriptor_1(LRFDescriptor_1_Params params);
		virtual ~LRFDescriptor_1();

		virtual void describeImpl(vtkPolyData* cloud, Feature3D* & feat, Desc & desc, int nPoints);

		void setParams(LRFDescriptor_1_Params params){m_params = params; m_descLength = 1;};	

	private:	
		LRFDescriptor_1_Params m_params;
	};


	//************************ MATCHING *******************************

	struct FeatureMatch
	{
		int refIndex, trgIndex;
		double matchScore;
	};


	enum EmatcherStrategy
	{
		E_MATCHER_STRATEGY_RATIO,
		E_MATCHER_STRATEGY_KNN_DISTANCE,
		E_MATCHER_STRATEGY_SINGLE_INDEX_RATIO
	};


	// Predicate to sort pairs in increasing order
	struct PairSort
	{
		inline bool operator () (const std::pair<int,int>& m1, const std::pair<int,int>& m2) const {return (m1.first == m2.first)? (m1.second < m2.second):(m1.first < m2.first);/*return m1.first < m2.first;*/ };
	};

	class IMatcher
	{
	public:

		IMatcher();
		virtual ~IMatcher();
		bool setReference(Desc & refDescriptors);
		int match(Desc & trgDescriptors, FeatureMatch* & featureMatch);
		void match(std::vector<Desc> &vDesc, std::vector< std::pair<int,int> > & pairs);

		void reset();
		int getModelMatches(const FeatureMatch * matches, const int nMatches, const int id, FeatureMatch * & fm);

		FeatureMatch * getFeatureMatch() { return m_featureMatch; };
		void setFeatureMatch(FeatureMatch* featureMatch) { m_featureMatch = featureMatch; };
		int getNumMatch() { return m_numMatches; };
		void setNumMatch(int numMatch) { m_numMatches = numMatch; } ;
		std::vector<FeatureMatch*> & getFeatureMatchVect() { return m_vFeatureMatch; };
		std::vector<int> & getNumMatchVect() { return m_vNumMatches; };
		std::vector<int> & getModelIndices() { return m_refModelIndices; };
		void setNumMatchVect(std::vector<int> &vNumMatch){m_vNumMatches = vNumMatch;};
		void setFeatureMatchVect(std::vector<FeatureMatch*> &vFeatureMatch){m_vFeatureMatch = vFeatureMatch;};

		int cleanFeatureMatch();
		virtual std::string toString() = 0;

		virtual void SetMaxMatchesNum(int maxKNN, float featureDistance = std::numeric_limits<float>::max() ){m_featureDistance = featureDistance; m_knn = maxKNN; m_matcherStrategy = E_MATCHER_STRATEGY_KNN_DISTANCE;};
		virtual void SetRatioMatching(){m_matcherStrategy = E_MATCHER_STRATEGY_RATIO;};
		virtual void SetSingleIndexRatioMatching(const std::vector<int> & refModelIndices);
		virtual void SetSingleIndexRatioMatching();

	protected:
		EmatcherStrategy m_matcherStrategy;
		int m_knn;
		float m_featureDistance;
		bool m_referenceInitialized;
		int m_numRefDescriptors;

		// single index matching vars
		std::vector<Desc> m_vDescForSingleIndex;
		std::vector<int> m_refModelIndices;
		int m_nTotalDesc;

		// output vars
		FeatureMatch *m_featureMatch;
		int m_numMatches;
		std::vector<int> m_vNumMatches;
		std::vector<FeatureMatch*> m_vFeatureMatch;

		virtual void resizeVect(const unsigned int newLength);
		virtual int matchImpl(Desc & trgDescriptors, FeatureMatch* & featureMatch) = 0;
		virtual bool setReferenceImpl(Desc & refDescriptors) = 0;
	};


	//********************* LRFMatcher_1 *********************

	/*!
	* \brief Parameters of matcher based on FLARE method as described in Petrelli A., Di Stefano L., "Pairwise registration by local orientation cues", Computer Graphics Forum, 2015
	*/
	struct LRFMatcher_1_Params
	{
		/* \brief Threshold T_d used to prune false LRF correspondes */
		float signedDistanceThresh;

		LRFMatcher_1_Params()	// default values
		{
			signedDistanceThresh = 0.01f;
		}
	};

	/*!
	* \brief Matcher based on FLARE method as described in Petrelli A., Di Stefano L., "Pairwise registration by local orientation cues", Computer Graphics Forum, 2015
	*/
	class LRFMatcher_1 : public IMatcher
	{
	public:
		LRFMatcher_1(float signedDistanceThresh = 0.01f);
		LRFMatcher_1(LRFMatcher_1_Params params);
		virtual ~LRFMatcher_1();

		virtual std::string toString()
		{
			std::stringstream str;
			str << "LRFMatcher_1";
			str << "  signedDistanceThresh = " << m_signedDistanceThresh;
			return str.str();
		};

		virtual void SetSignedDistanceThresh(float signedDistanceThresh){m_signedDistanceThresh = signedDistanceThresh;};

		static void Sort_wrtFirstFeat(Desc & desc, std::vector<int> &vId_sorted);
		static void ReplaceIds(FeatureMatch* &featureMatch, const int nMatches, const std::vector<int> &vNewIds_trg, const std::vector<int> &vNewIds_ref);

	protected:
		float						m_signedDistanceThresh;

		Desc						m_desc_ref;
		std::vector<int>			m_vId_sorted_descRef;

		int matchImpl(Desc & trgDescriptors, FeatureMatch* & featureMatch);
		bool setReferenceImpl(Desc & refDescriptors);
	};

	


}


#endif



