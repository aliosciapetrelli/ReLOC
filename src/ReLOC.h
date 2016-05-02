#ifndef RELOC
#define RELOC

#include "ReLOC_RANSAC.h"
#include "ReLOC_Hough.h"

#include <vector>

#include "vtkTransform.h"
#include "vtkPolyDataCollection.h"



namespace ReLOC
{

	/** \brief Pairwise 3D registration algorithm "Registration by Local Orientation Cues" (ReLOC) proposed in:
	*  Petrelli A., Di Stefano L., "Pairwise registration by local orientation cues", Computer Graphics Forum, 2015
	*
	* \author Alioscia Petrelli
	*/
	class ReLOC
	{
	public:
		ReLOC();
		virtual ~ReLOC();

		/** \brief Return the rigid motion that moves mesh_ref against mesh_trg.
		*
		*/
		vtkTransform* Align(vtkPolyData* mesh_trg, vtkPolyData* mesh_ref);


		/** \brief perform the registration of each pair of views listed in vViewPairs.
		* \param cMeshes partial views to align.
		* \param meanMeshRes The mean mesh resolution of cMeshes. If not provided, the method computes the mean mesh resolution of cMeshes.
		* \param vViewPairs Pairs of views the algorithm attempts to align. If not provided the algorithm attempts to align each pair af views.
		* \return List of pairwise rigid motions that align the pair of views lited in vViewPairs. Each rigid motion represents the rigid motion that moves the second view against the first.
		*/
		std::vector<vtkTransform*> &Align(vtkPolyDataCollection* cMeshes, const double meanMeshRes = -1.0, std::vector< std::pair<int,int> > &vViewPairs = std::vector< std::pair<int,int> >());

		/** \brief Initialize the algorithm, extract feature points and computes local reference frames by applying FLARE algorithm.
		* \param cMeshes partial views to align.
		* \param meanMeshRes The mean mesh resolution of cMeshes. If not provided, the method computes the mean mesh resolution of cMeshes.
		* \param vViewPairs indices of the pairs of views the algorithm attempts to align. If not provided the algorithm attempts to align each pair af views.
		*/
		void PrepareRegistration(vtkPolyDataCollection* cMeshes, const double meanMeshRes = -1.0, std::vector< std::pair<int,int> > &vViewPairs = std::vector< std::pair<int,int> >());

		/** \brief Return the rigid motion that moves views idx_ref against idx_trg.
		*
		* PrepareRegistration() method has to be previously called.
		*/
		vtkTransform* Align(const int idx_trg, const int idx_ref);


		void Verbose(const bool verbose){m_verbose = verbose;}

		/** \brief Use random detector instead of the detector that extracts flat points. */
		void UseRandomDetector(const bool useRandomDetector){m_useRandomDetector = useRandomDetector;}

		RandomDetectorParams &GetRandomDetectorParams(){return m_randomDetectorParams;}

		/** \brief Get the flat points detector. Useful to set R_f, R_discard, R1_search, R2_search, T1_search, T2_search parameters. */
		FlatPoints3DDetectorParams &GetFlatPoints3DDetectorParams(){return m_flatPoints3DDetectorParams;}
		
		/** \brief Get the FLARE object. Useful to set R_z and R_x parameters. */
		FLAREparams &GetFLAREParams(){return m_LRFDescriptor_1_Params.flareParams;}

		/** \brief Get the RANSAC object. Useful to set T_ransac, N_ransac and P_ransac parameters. */
		Ransac &GetRansac(){return m_RANSAC;}

		/** \brief Get the matcher of LRF computed with FLARE algorithm. Useful to set T_D parameter. */
		LRFMatcher_1_Params &GetFLAREMatcherParams(){return m_LRFMatcher_1_Params;}

		/** \brief Set the extent of the Hough space in unit of the standard deviation of points coordinates.*/
		void SetHoughSigmaFactor(const float houghSigmaFactor){m_houghSigmaFactor = houghSigmaFactor;}

		/** \brief Set  the bin size of the Hough space in unit of mesh resolution (Sbin).*/
		void SetHoughBinSizeInMeshRes(const float houghBinSizeInMeshRes){m_houghBinSizeInMeshRes = houghBinSizeInMeshRes;}

		/** \brief Set the factor of magnification of the Hough space (f_hough).*/
		void SetHoughBBoxSizeFactor(const float houghBBoxSizeFactor){m_houghBBoxSizeFactor = houghBBoxSizeFactor;}


	protected:
		vtkPolyDataCollection*				m_cMeshes;
		int									m_idx_trg_prev;
		double								m_meanMeshRes;
		bool								m_verbose;
		

		std::vector<vtkTransform*>			m_vPairwiseTransforms;

		bool								m_useRandomDetector;
		I3DDetector*						m_ptrDetector;
		RandomDetectorParams				m_randomDetectorParams;
		FlatPoints3DDetectorParams			m_flatPoints3DDetectorParams;

		I3DDescriptor*						m_ptrDescriptor;
		LRFDescriptor_1_Params				m_LRFDescriptor_1_Params;
		std::vector<Desc>					m_vDescs;

		IMatcher*							m_ptrMatcher;
		LRFMatcher_1_Params					m_LRFMatcher_1_Params;
		std::vector< std::pair<int,int> >	m_vViewPairs;

		pstHough<unsigned short int, bool, pstHoughContinuousVoteStrategy, pstHoughNeighborhoodFindMaximaStrategy >* m_houghSpace;
		std::vector< std::vector<float> >	m_vvCentroidsWrtLRF;
		std::vector<double>					m_vMeshBBoxes;
		float								m_houghSigmaFactor;
		std::vector<float>					m_vMeshCentroids;
		float								m_houghBinSizeInMeshRes;
		float								m_houghBBoxSizeFactor;
		
		Ransac								m_RANSAC;
		


		void deallocate();

	};

}

#endif