#ifndef RELOC_HOUGH_H
#define RELOC_HOUGH_H

#include "ReLOC_LRF.h"


#include <vector>


namespace ReLOC
{

	template <typename CodomainType> struct ContinuousVote
	{
		ContinuousVote(int _id) : id(_id)
		{
		};

		//ContinuousVote(double* _coordinate, int nDim, CodomainType _vote, int _id)
		//{

		//	for (int i = 0; i < nDim; i++)
		//		coordinate.push_back(_coordinate[i]);
		//	vote = _vote;
		//	id = _id;
		//};

		//std::vector<double> coordinate;
		//CodomainType vote;

		int id;
	};


	template <typename CodomainType> struct MaximumContainer
	{
		MaximumContainer()
		{
			vote = 0;
			additionalInfoIndex = -1;
		};

		MaximumContainer(std::vector <double> &coordinates, CodomainType _vote, const std::vector <ContinuousVote<CodomainType> > &_continuousVote = std::vector <ContinuousVote<CodomainType> >(), unsigned int _additionalInfoIndex = -1)
		{
			coordinate = coordinates;
			vote = _vote;
			continuousVote = _continuousVote;
			additionalInfoIndex = _additionalInfoIndex;
		};

		std::vector<double> coordinate;
		CodomainType vote;
		std::vector < ContinuousVote<CodomainType> > continuousVote;
		unsigned int additionalInfoIndex;
	};




	template <typename CodomainType, typename AdditionalInfoType = char, class VoteStrategy = pstHoughSimpleVoteStrategy, class FindMaximaStrategy = pstHoughSimpleFindMaximaStrategy> class pstHough{

		friend class pstHoughSimpleVoteStrategy;
		friend class pstHoughContinuousVoteStrategy;
		friend class pstHoughContinuousInterpolatedVoteStrategy;

		friend class pstHoughNeighborhoodFindMaximaStrategy;
		friend class pstHoughFindContinuousMaximaStrategy;
		friend class pstHoughSimpleFindMaximaStrategy;
		friend class pstHoughSumVotesFindMaximaStrategy;

	public:

		pstHough(int nDim, const double *Mins, const double *Steps, const int *Bins, bool useAdditionalInfo=true){

			m_nDim = nDim;

			m_Mins = new double [nDim];
			memcpy(m_Mins, Mins, nDim*sizeof(double));

			m_Steps = new double [nDim];
			memcpy(m_Steps, Steps, nDim*sizeof(double));

			m_Bins = new int [nDim];
			memcpy(m_Bins, Bins, nDim*sizeof(int));

			m_partialBinProducts = new int [nDim+1];
			m_partialBinProducts[0] = 1;

			for(int i=1; i<=nDim; i++)
				m_partialBinProducts[i] = m_Bins[i-1]*m_partialBinProducts[i-1];

			m_HoughSpace = new CodomainType [m_partialBinProducts[nDim]];
			m_ContinuousVoteSpace = new std::vector < ContinuousVote<CodomainType> >[m_partialBinProducts[nDim]];
			memset(m_HoughSpace, 0, m_partialBinProducts[nDim] * sizeof(CodomainType));

			if(useAdditionalInfo)
				m_AdditionalInfoSpace = new std::vector <AdditionalInfoType> [m_partialBinProducts[nDim]];
			else
				m_AdditionalInfoSpace = NULL;

		};


		pstHough(int nDim, double *Mins, int *Bins, double *Maxs, bool useAdditionalInfo=true){

			m_nDim = nDim;

			m_Mins = new double [nDim];
			memcpy(m_Mins, Mins, nDim*sizeof(double));

			m_Steps = new double [nDim];

			m_Bins = new int [nDim];
			memcpy(m_Bins, Bins, nDim*sizeof(int));

			for(int i=0; i<nDim; i++){
				m_Steps[i] = (Maxs[i]-Mins[i]) / Bins[i];
			}

			m_partialBinProducts = new int [nDim+1];
			m_partialBinProducts[0] = 1;
			for(int i=1; i<=nDim; i++)
				m_partialBinProducts[i] = m_Bins[i-1]*m_partialBinProducts[i-1];

			m_HoughSpace = new CodomainType [m_partialBinProducts[nDim]];
			m_ContinuousVoteSpace = new std::vector < ContinuousVote<CodomainType> >[m_partialBinProducts[nDim]];

			memset(m_HoughSpace, 0, m_partialBinProducts[nDim] * sizeof(CodomainType));

			if(useAdditionalInfo)
				m_AdditionalInfoSpace = new std::vector <AdditionalInfoType> [m_partialBinProducts[nDim]];
			else
				m_AdditionalInfoSpace = NULL;

		};

		pstHough(int nDim, const double *Mins, const double *Steps, const double *Maxs, bool useAdditionalInfo = true){

			m_nDim = nDim;

			m_Mins = new double [nDim];
			memcpy(m_Mins, Mins, nDim*sizeof(double));

			m_Steps = new double [nDim];
			memcpy(m_Steps, Steps, nDim*sizeof(double));

			m_Bins = new int [nDim];

			for(int i=0; i<nDim; i++){
				m_Bins[i] = (int)floor( (Maxs[i]-Mins[i]) / m_Steps[i] );
			}

			m_partialBinProducts = new int [nDim+1];
			m_partialBinProducts[0] = 1;
			for(int i=1; i<=nDim; i++)
				m_partialBinProducts[i] = m_Bins[i-1]*m_partialBinProducts[i-1];

			m_HoughSpace = new CodomainType [m_partialBinProducts[nDim]];
			m_ContinuousVoteSpace = new std::vector < ContinuousVote<CodomainType> >[m_partialBinProducts[nDim]];
			memset(m_HoughSpace, 0, m_partialBinProducts[nDim] * sizeof(CodomainType));

			if(useAdditionalInfo)
				m_AdditionalInfoSpace = new std::vector <AdditionalInfoType> [m_partialBinProducts[nDim]];
			else
				m_AdditionalInfoSpace = NULL;

		};


		~pstHough(){

			delete [] m_HoughSpace;
			delete [] m_ContinuousVoteSpace;

			if(m_AdditionalInfoSpace)
				delete [] m_AdditionalInfoSpace;
			delete [] m_Mins;
			delete [] m_Bins;
			delete [] m_Steps;

		};

		void Reset() {
			memset(m_HoughSpace, 0, m_partialBinProducts[m_nDim] * sizeof(CodomainType));

			for(int k=0; k < m_partialBinProducts[m_nDim]; k++)
				m_ContinuousVoteSpace[k].clear();

			if (m_AdditionalInfoSpace)
				for(int k=0; k < m_partialBinProducts[m_nDim]; k++)
					m_AdditionalInfoSpace[k].clear();

		};

		void findMaxima(CodomainType Th, std::vector <MaximumContainer<CodomainType> > &maximumContainer){
			m_findMaximaStrategy.findMaxima(this, Th, maximumContainer);
		};

		void findMaxima(CodomainType Th, std::vector < MaximumContainer<CodomainType> > &maximumContainer, std::vector <std::vector <AdditionalInfoType> > &AdditionalInfo) {
			m_findMaximaStrategy.findMaxima(this, Th, maximumContainer, AdditionalInfo);
		};

		void findMaximum(MaximumContainer<CodomainType>  &maximum){
			m_findMaximaStrategy.findMaximum(this, maximum);
		};

		int Vote(double *singleVote, CodomainType weight, int id){

			return m_voteStrategy.Vote(this, singleVote, weight, id);

		};

		int Vote(double *singleVote, CodomainType weight, int id, AdditionalInfoType &AdditionalInfo){

			return m_voteStrategy.Vote(this, singleVote, weight, id, AdditionalInfo);

		};

		template <typename CodomainType>
		void FilterHoughSpace(int kSize = 3)
		{
			int hk = (kSize - 1) / 2;
			cv::Mat gKernel = cv::getGaussianKernel(kSize, 1.0, CV_32F);

			float *K = new float[kSize];
			for(int i = 0; i < kSize; i++)
				K[i] =  gKernel.at<float>(i);

			int currBin = 0;
			int moduledIndex;
			CodomainType temp;
			CodomainType * tempHough = NULL;
			CodomainType * outHough = new CodomainType [m_partialBinProducts[m_nDim]];
			double* coordCheck1 = new double[m_nDim];
			double* coordCheck2 = new double[m_nDim];
			bool outOfBounds;

			for(int d = 0; d < m_nDim; d++)
			{
				for(int i = 0; i < m_partialBinProducts[m_nDim]; i++)
				{
					temp = 0;

					// get d-th coordinate of i-th bin
					moduledIndex = i;

					for(int k = m_nDim - 1; k >= 0; k--)
					{
						moduledIndex = moduledIndex % m_partialBinProducts[k+1];
						coordCheck1[k] = moduledIndex / m_partialBinProducts[k];
					}

					for(int h = - hk; h <= hk; h++)
					{
						// get d-th coordinate of (i-h)th bin
						currBin = i - h * m_partialBinProducts[d];
						moduledIndex = currBin;

						for(int k = m_nDim - 1; k >= 0; k--)
						{
							moduledIndex = moduledIndex % m_partialBinProducts[k+1];
							coordCheck2[k] = moduledIndex / m_partialBinProducts[k];
						}

						// only the d-th coordinate can change, if the other coordinates change, it means that we have reached the maximum size for that dimension
						outOfBounds = false;
						for(int c = 0; c < m_nDim; c++)
							if(c != d && coordCheck1[c] != coordCheck2[c])
								outOfBounds = true;

						// accumulate
						if(!outOfBounds && currBin > 0 && currBin < m_partialBinProducts[m_nDim])
							temp += (CodomainType) (m_HoughSpace[currBin] * K[hk - h]);
					}

					outHough[i] = temp;
				}

				// swap
				tempHough = m_HoughSpace;
				m_HoughSpace = outHough;
				outHough = tempHough;
			}

			// releases
			delete [] K;
			delete [] outHough;
			delete [] coordCheck1;
			delete [] coordCheck2;
		}

		inline double getMin(int dimension)
		{
			return (dimension >= m_nDim) ? 0 : m_Mins[dimension];
		};

		inline int getBin(int dimension)
		{
			return (dimension >= m_nDim) ? 0 : m_Bins[dimension];
		};

		inline double getBinLength(int dimension)
		{
			return (dimension >= m_nDim) ? 0 : m_Steps[dimension];
		};

	private:

		int m_nDim;
		double *m_Mins;
		double *m_Steps;
		int *m_Bins;
		int *m_partialBinProducts;

		CodomainType *m_HoughSpace;
		std::vector <AdditionalInfoType> *m_AdditionalInfoSpace;
		std::vector < ContinuousVote<CodomainType> > *m_ContinuousVoteSpace;

		VoteStrategy m_voteStrategy;
		FindMaximaStrategy m_findMaximaStrategy;

		void findMaxima(CodomainType Th, std::vector< std::vector < int > > &MaximaCoord, std::vector<int> & MaximaBin)
		{
			MaximaCoord.clear();
			MaximaBin.clear();
			double* indexes = new double[m_nDim];

			for(int i=0; i<m_partialBinProducts[m_nDim]; i++){

				if (m_HoughSpace[i] < Th)
					continue;

				//Check with neighbors
				bool isMaximum = true;
				int moduledIndex = i;

				for(int k = m_nDim - 1; k >= 0; k--){

					moduledIndex = moduledIndex % m_partialBinProducts[k+1];
					indexes[k] = moduledIndex / m_partialBinProducts[k];

					if( indexes[k] > 0 && m_HoughSpace[i]<m_HoughSpace[i-m_partialBinProducts[k]]){
						isMaximum = false;
						break;
					}
					if(	indexes[k] < m_Bins[k]-1 && m_HoughSpace[i]<m_HoughSpace[i+m_partialBinProducts[k]] ){
						isMaximum = false;
						break;
					}
				}

				if(isMaximum)
				{
					std::vector<int> temp;
					for(int k=0; k < m_nDim; k++)
						temp.push_back((int)indexes[k]);

					MaximaCoord.push_back(temp);
					MaximaBin.push_back(i);
				}
			}
			delete [] indexes;
		};

	};


	//put in vvCentroidWrtLRF the centroids with respect to the basis of LRFs
	void GetCentroidsWrtLRFs(const std::vector<float> &vCentroids, const std::vector<Feature3D*> &vFeatures3D, const std::vector<int> &vNumFeatures3D, std::vector< std::vector<float> > &vvCentroidWrtLRF);


	template<typename CodomainType, class VoteStrategy, class MaximaStrategy> static void DestroyHoughSpace(pstHough<CodomainType, bool, VoteStrategy, MaximaStrategy >* &houghSpace)
	{
#ifndef _DEBUG
		if(houghSpace)
		{
			delete houghSpace;
			houghSpace = NULL;
		}
#endif
	}

	template<typename CodomainType, class VoteStrategy, class MaximaStrategy> void CreateHoughSpace(pstHough<CodomainType, bool, VoteStrategy, MaximaStrategy >* &houghSpace, double* bbox_trg, const float stepSize, const double bboxSizeFactor = 1.3)
	{
		double bounds_ext[6];
		memcpy(bounds_ext, bbox_trg, sizeof(double)*6);
		ExtendBBox(bounds_ext, bboxSizeFactor, bboxSizeFactor, bboxSizeFactor);

		int nDims = 3;
		double* vMin = new double[nDims];
		double* vMax = new double[nDims];
		double* vStep = new double[nDims];
		for(int bi = 0; bi < nDims; bi++)
		{
			vMin[bi] = bounds_ext[2*bi] - (bounds_ext[2*bi+1]-bounds_ext[2*bi])/2.0;
			vMax[bi] = bounds_ext[2*bi+1] + (bounds_ext[2*bi+1]-bounds_ext[2*bi])/2.0;
			vStep[bi] = stepSize;
		}

		//cout << "Hough Size: " << endl;
		//cout << vMin[0] << "\t" << vMax[0] << "\t" << (vMax[0] - vMin[0]) << "\t" << (int)((vMax[0] - vMin[0])/vStep[0]) << endl;
		//cout << vMin[1] << "\t" << vMax[1] << "\t" << (vMax[1] - vMin[1]) << "\t" << (int)((vMax[1] - vMin[1])/vStep[1]) << endl;
		//cout << vMin[2] << "\t" << vMax[2] << "\t" << (vMax[2] - vMin[2]) << "\t" << (int)((vMax[2] - vMin[2])/vStep[2]) << endl;
		//int size = (int)((vMax[0] - vMin[0])/vStep[0]);
		//size *= (int)((vMax[1] - vMin[1])/vStep[1]);
		//size *= (int)((vMax[2] - vMin[2])/vStep[2]);
		//cout << "Size: " << size << endl;

		houghSpace = new pstHough<CodomainType, bool, VoteStrategy, MaximaStrategy >(nDims, vMin, vStep, vMax, false);

		delete[] vMin;
		delete[] vMax;
		delete[] vStep;

	}




	class pstHoughContinuousVoteStrategy {

	public:

		template <typename CodomainType, typename AdditionalInfoType, typename FindMaximaStrategy> 
		int Vote(pstHough<CodomainType, AdditionalInfoType, pstHoughContinuousVoteStrategy, FindMaximaStrategy> *hough, double *singleVote, CodomainType weight, int id){

			int index = 0;

			for(int i=0; i<hough->m_nDim; i++){
				int currentBin = (int)floor((singleVote[i]-hough->m_Mins[i])/hough->m_Steps[i]);
				if(currentBin < 0 || currentBin >= hough->m_Bins[i]){
					//printf("Current Vote goes out of bounds in the Hough Table!\nDimension: %d, Value inserted: %f, Min value: %f, Max value: %f\n", i, singleVote[i], m_Mins[i], m_Mins[i] + m_Steps[i]*m_Bins[i]);
					return -1;
				}
				index += hough->m_partialBinProducts[i] * currentBin;
			}
			hough->m_HoughSpace[index] += weight;

			hough->m_ContinuousVoteSpace[index].push_back(ContinuousVote<CodomainType>(id));
			//hough->m_ContinuousVoteSpace[index].push_back(ContinuousVote<CodomainType>(singleVote, hough->m_nDim, weight, id));

			return index;
		};

		template <typename CodomainType, typename AdditionalInfoType, typename FindMaximaStrategy> 
		int Vote(pstHough<CodomainType, AdditionalInfoType, pstHoughContinuousVoteStrategy, FindMaximaStrategy> *hough, double *singleVote, CodomainType weight, int id, AdditionalInfoType &AdditionalInfo){

			int index = Vote(hough, singleVote, weight, id);
			if(index != -1)
				hough->m_AdditionalInfoSpace[index].push_back( AdditionalInfo);
			return index;
		};

		//private:

	};



	class pstHoughNeighborhoodFindMaximaStrategy {

	public:

		template <typename CodomainType, typename AdditionalInfoType, typename VoteStrategy>
		void findMaxima(pstHough<CodomainType, AdditionalInfoType, VoteStrategy, pstHoughNeighborhoodFindMaximaStrategy> *hough, CodomainType Th, std::vector <MaximumContainer<CodomainType> > &maximumContainer)
		{
			maximumContainer.clear();

			std::vector< std::vector<int> > maximaCoord;
			std::vector<int> maximaBin;

			hough->findMaxima(1, maximaCoord, maximaBin);
			//hough->findMaxima(Th, maximaCoord, maximaBin);

			for (unsigned int i = 0; i < maximaBin.size(); i++)
			{
				std::vector <double> temp;
				temp.resize(hough->m_nDim);
				double voteSum = 0.0;

				voteSum = hough->m_ContinuousVoteSpace[maximaBin[i]].size();

				//for (unsigned int j = 0; j < hough->m_ContinuousVoteSpace[maximaBin[i]].size(); j++)
				//{
				//	for(int k=0; k < hough->m_nDim; k++)

				//		temp[k] += hough->m_ContinuousVoteSpace[maximaBin[i]][j].vote * hough->m_ContinuousVoteSpace[maximaBin[i]][j].coordinate[k];

				//	voteSum += hough->m_ContinuousVoteSpace[maximaBin[i]][j].vote;
				//}

				for(int k=0; k < hough->m_nDim; k++)
					temp[k] /= voteSum;

				std::vector < ContinuousVote<CodomainType> > neighbourhoodContinuousVoteSpace;
				CodomainType nVotes = sumVotes(hough, maximaCoord[i], neighbourhoodContinuousVoteSpace);

				if (nVotes < Th)
					continue;

				maximumContainer.push_back( MaximumContainer<CodomainType>(temp, nVotes, neighbourhoodContinuousVoteSpace));
			}
		};

		template <typename CodomainType, typename AdditionalInfoType, typename VoteStrategy>
		void findMaxima(pstHough<CodomainType, AdditionalInfoType, VoteStrategy, pstHoughNeighborhoodFindMaximaStrategy> *hough, CodomainType Th, std::vector <MaximumContainer<CodomainType> > &maximumContainer, std::vector <std::vector <AdditionalInfoType> > &AdditionalInfo)
		{
			AdditionalInfo.clear();
			maximumContainer.clear();

			std::vector< std::vector<int> > maximaCoord;
			std::vector<int> maximaBin;

			hough->findMaxima(Th, maximaCoord, maximaBin);

			for (unsigned int i = 0; i < maximaBin.size(); i++)
			{
				std::vector <double> temp;
				temp.resize(hough->m_nDim);
				double voteSum = 0.0;

				for (unsigned int j = 0; j < hough->m_ContinuousVoteSpace[maximaBin[i]].size(); j++)
				{
					for(int k=0; k < hough->m_nDim; k++)

						temp[k] += hough->m_ContinuousVoteSpace[maximaBin[i]][j].vote * hough->m_ContinuousVoteSpace[maximaBin[i]][j].coordinate[k];

					voteSum += hough->m_ContinuousVoteSpace[maximaBin[i]][j].vote;
				}

				for(int k=0; k < hough->m_nDim; k++)
					temp[k] /= voteSum;

				std::vector < ContinuousVote<CodomainType> > neighbourhoodContinuousVoteSpace;
				std::vector < AdditionalInfoType > neighbourhoodAdditionalInfoSpace;
				CodomainType nVotes = sumVotes(hough, maximaCoord[i], neighbourhoodContinuousVoteSpace, neighbourhoodAdditionalInfoSpace);

				maximumContainer.push_back( MaximumContainer<CodomainType>(temp, nVotes, neighbourhoodContinuousVoteSpace, (unsigned int) AdditionalInfo.size()));
				AdditionalInfo.push_back(neighbourhoodAdditionalInfoSpace);
			}
		};

		template <typename CodomainType, typename AdditionalInfoType, typename VoteStrategy>
		void findMaximum(pstHough<CodomainType, AdditionalInfoType, VoteStrategy, pstHoughNeighborhoodFindMaximaStrategy> *hough, MaximumContainer<CodomainType> &maximum)
		{
			std::vector< std::vector<int> > maximaCoord;
			std::vector<int> maximaBin;

			//find all displacements
			int nNeighs = (int)pow((float)3, hough->m_nDim);
			vector<int> vDispl(nNeighs);
			vector<vector<int>> vvActiveDim(nNeighs, vector<int>(hough->m_nDim));
			vDispl[0] = -1;
			vvActiveDim[0][0] = -1;
			vvActiveDim[0][1] = 0;
			vvActiveDim[0][2] = 0;
			vDispl[1] = 0;
			vvActiveDim[1][0] = 0;
			vvActiveDim[1][1] = 0;
			vvActiveDim[1][2] = 0;
			vDispl[2] = +1;
			vvActiveDim[2][0] = 1;
			vvActiveDim[2][1] = 0;
			vvActiveDim[2][2] = 0;
			size_t vectorSize;
			int iNeigh = 3;
			for(int di = 1; di<hough->m_nDim; di++)
			{
				vectorSize = iNeigh;
				for(size_t si=0; si<vectorSize; si++)
				{
					vDispl[iNeigh] = -hough->m_partialBinProducts[di] + vDispl[si];
					vvActiveDim[iNeigh] = vvActiveDim[si];
					vvActiveDim[iNeigh][di] = -1;
					iNeigh++;

					vDispl[iNeigh] = hough->m_partialBinProducts[di] + vDispl[si];
					vvActiveDim[iNeigh] = vvActiveDim[si];
					vvActiveDim[iNeigh][di] = 1;
					iNeigh++;
				}
			}
			vectorSize = vDispl.size();

			CodomainType maxSum = 0;

			if(hough->m_nDim == 3)
			{
				CodomainType *ptrHoughSpace = hough->m_HoughSpace;
				int *ptrDispl;

				int width = hough->m_Bins[0] - 1;
				int height = hough->m_Bins[1] - 1;
				int depth = hough->m_Bins[2] - 1;
				int maxWidth = 0;
				int maxHeight = 0;
				int maxDepth = 0;

				ptrHoughSpace += hough->m_partialBinProducts[2];
				for(int de=1; de<depth; de++)
				{
					ptrHoughSpace += hough->m_partialBinProducts[1];
					for(int he=1; he<height; he++)
					{
						ptrHoughSpace++;
						for(int wi=1; wi<width; wi++)
						{
							if(*ptrHoughSpace != 0)
							{
								//sum all the 27 neighs
								ptrDispl = &vDispl[0];
								CodomainType voteSum = 0;
								for(size_t si=0; si<vectorSize; si++)
								{
									voteSum += *(ptrHoughSpace + *ptrDispl);
									ptrDispl++;
								}
								if(voteSum > maxSum)
								{
									maxSum = voteSum;
									maxWidth = wi;
									maxHeight = he;
									maxDepth = de;
								}
							}
							ptrHoughSpace++;
						}
						ptrHoughSpace++;
					}
					ptrHoughSpace += hough->m_partialBinProducts[1];
				}

				//consider borders
				int deOut=0;
				vector<int>* ptrActiveDim;
				ptrHoughSpace = hough->m_HoughSpace;
				for(int he=0; he<(height+1); he++)
				{
					for(int wi=0; wi<(width+1); wi++)
					{
						if(*ptrHoughSpace != 0)
						{
							//sum all the 27 neighs
							ptrDispl = &vDispl[0];
							ptrActiveDim = &vvActiveDim[0];
							CodomainType voteSum = 0;
							for(size_t si=0; si<vectorSize; si++)
							{
								if( (*ptrActiveDim)[2] != -1)
								{
									if( ((*ptrActiveDim)[1] != -1) || (he!=0) )
									{	
										if( ((*ptrActiveDim)[1] != 1) || (he!=height) )
										{
											if( ((*ptrActiveDim)[0] != -1) || (wi!=0) )
											{	
												if( ((*ptrActiveDim)[0] != 1) || (wi!=width) )
												{
													voteSum += *(ptrHoughSpace + *ptrDispl);
												}
											}
										}
									}
								}
								ptrDispl++;
								ptrActiveDim++;
							}
							if(voteSum > maxSum)
							{
								maxSum = voteSum;
								maxWidth = wi;
								maxHeight = he;
								maxDepth = deOut;
							}
						}
						ptrHoughSpace++;
					}
				}

				deOut=height;
				ptrHoughSpace = hough->m_HoughSpace + depth*(height+1)*(width+1);
				for(int he=0; he<(height+1); he++)
				{
					for(int wi=0; wi<(width+1); wi++)
					{
						if(*ptrHoughSpace != 0)
						{
							//sum all the 27 neighs
							ptrDispl = &vDispl[0];
							ptrActiveDim = &vvActiveDim[0];
							CodomainType voteSum = 0;
							for(size_t si=0; si<vectorSize; si++)
							{
								if( (*ptrActiveDim)[2] != 1)
								{
									if( ((*ptrActiveDim)[1] != -1) || (he!=0) )
									{	
										if( ((*ptrActiveDim)[1] != 1) || (he!=height) )
										{
											if( ((*ptrActiveDim)[0] != -1) || (wi!=0) )
											{	
												if( ((*ptrActiveDim)[0] != 1) || (wi!=width) )
												{
													voteSum += *(ptrHoughSpace + *ptrDispl);
												}
											}
										}
									}
								}
								ptrDispl++;
								ptrActiveDim++;
							}
							if(voteSum > maxSum)
							{
								maxSum = voteSum;
								maxWidth = wi;
								maxHeight = he;
								maxDepth = deOut;
							}
						}
						ptrHoughSpace++;
					}
				}


				int heOut=0;
				ptrHoughSpace = hough->m_HoughSpace + (height+1)*(width+1);
				for(int de=1; de<(depth); de++)
				{
					for(int wi=0; wi<(width+1); wi++)
					{
						if(*ptrHoughSpace != 0)
						{
							//sum all the 27 neighs
							ptrDispl = &vDispl[0];
							ptrActiveDim = &vvActiveDim[0];
							CodomainType voteSum = 0;
							for(size_t si=0; si<vectorSize; si++)
							{
								if( (*ptrActiveDim)[1] != -1)
								{
									if( ((*ptrActiveDim)[0] != -1) || (wi!=0) )
									{	
										if( ((*ptrActiveDim)[0] != 1) || (wi!=width) )
										{
											voteSum += *(ptrHoughSpace + *ptrDispl);
										}
									}
								}
								ptrDispl++;
								ptrActiveDim++;
							}
							if(voteSum > maxSum)
							{
								maxSum = voteSum;
								maxWidth = wi;
								maxHeight = heOut;
								maxDepth = de;
							}
						}
						ptrHoughSpace++;
					}
					ptrHoughSpace += height*(width+1);
				}

				heOut=height;
				ptrHoughSpace = hough->m_HoughSpace + (height+1)*(width+1) + (width+1)*height;
				for(int de=1; de<(depth); de++)
				{
					for(int wi=0; wi<(width+1); wi++)
					{
						if(*ptrHoughSpace != 0)
						{
							//sum all the 27 neighs
							ptrDispl = &vDispl[0];
							ptrActiveDim = &vvActiveDim[0];
							CodomainType voteSum = 0;
							for(size_t si=0; si<vectorSize; si++)
							{
								if( (*ptrActiveDim)[1] != +1)
								{
									if( ((*ptrActiveDim)[0] != -1) || (wi!=0) )
									{	
										if( ((*ptrActiveDim)[0] != 1) || (wi!=width) )
										{
											voteSum += *(ptrHoughSpace + *ptrDispl);
										}
									}
								}
								ptrDispl++;
								ptrActiveDim++;
							}
							if(voteSum > maxSum)
							{
								maxSum = voteSum;
								maxWidth = wi;
								maxHeight = heOut;
								maxDepth = de;
							}
						}
						ptrHoughSpace++;
					}
					ptrHoughSpace += height*(width+1);
				}


				int wiOut=0;
				ptrHoughSpace = hough->m_HoughSpace + (height+1)*(width+1) + (width+1);
				for(int de=1; de<depth; de++)
				{
					for(int he=1; he<height; he++)
					{
						if(*ptrHoughSpace != 0)
						{
							//sum all the 27 neighs
							ptrDispl = &vDispl[0];
							ptrActiveDim = &vvActiveDim[0];
							CodomainType voteSum = 0;
							for(size_t si=0; si<vectorSize; si++)
							{
								if( (*ptrActiveDim)[0] != -1)
								{
									voteSum += *(ptrHoughSpace + *ptrDispl);
								}
								ptrDispl++;
								ptrActiveDim++;
							}
							if(voteSum > maxSum)
							{
								maxSum = voteSum;
								maxWidth = wiOut;
								maxHeight = he;
								maxDepth = de;
							}
						}
						ptrHoughSpace += (width+1);
					}
					ptrHoughSpace += (width+1);
				}

				wiOut=width;
				ptrHoughSpace = hough->m_HoughSpace + (height+1)*(width+1) + (width+1) + width;
				for(int de=1; de<depth; de++)
				{
					for(int he=1; he<height; he++)
					{
						if(*ptrHoughSpace != 0)
						{
							//sum all the 27 neighs
							ptrDispl = &vDispl[0];
							ptrActiveDim = &vvActiveDim[0];
							CodomainType voteSum = 0;
							for(size_t si=0; si<vectorSize; si++)
							{
								if( (*ptrActiveDim)[0] != 1)
								{
									voteSum += *(ptrHoughSpace + *ptrDispl);
								}
								ptrDispl++;
								ptrActiveDim++;
							}
							if(voteSum > maxSum)
							{
								maxSum = voteSum;
								maxWidth = wiOut;
								maxHeight = he;
								maxDepth = de;
							}
						}
						ptrHoughSpace += (width+1);
					}
					ptrHoughSpace += (width+1);
				}


				width++;
				height++;
				depth++;

				int maxIndex = height*width*maxDepth + width*maxHeight + maxWidth;

				maximum.vote = maxSum;

				//compute coordinates and neighbourhoodContinuousVoteSpace
				maximum.coordinate.clear();
				maximum.coordinate.resize(hough->m_nDim, 0.0);
				maximum.continuousVote.clear();

				ptrDispl = &vDispl[0];
				int currIndex;
				CodomainType voteSum = 0;
				for(size_t si=0; si<vectorSize; si++)
				{
					currIndex = maxIndex + *ptrDispl;

					//find average coordinates
					voteSum = hough->m_ContinuousVoteSpace[currIndex].size();
					//for (unsigned int j = 0; j < hough->m_ContinuousVoteSpace[currIndex].size(); j++)
					//{
					//	for(int k=0; k < hough->m_nDim; k++)

					//		maximum.coordinate[k] += hough->m_ContinuousVoteSpace[currIndex][j].vote * hough->m_ContinuousVoteSpace[currIndex][j].coordinate[k];

					//	voteSum += hough->m_ContinuousVoteSpace[currIndex][j].vote;
					//}

					maximum.continuousVote.insert(maximum.continuousVote.end(), hough->m_ContinuousVoteSpace[currIndex].begin(), hough->m_ContinuousVoteSpace[currIndex].end());


					ptrDispl++;
				}
				for(int k=0; k < hough->m_nDim; k++)
					maximum.coordinate[k] /= voteSum;


			}
			else
			{
				cout<<"ERROR: findMaximum for dimension " << hough->m_nDim << " is not implemented yet. Sorry!";
				getchar();
			}

		};


	private:

		template <typename CodomainType, typename AdditionalInfoType, typename VoteStrategy>
		CodomainType sumVotes(pstHough<CodomainType, AdditionalInfoType, VoteStrategy, pstHoughNeighborhoodFindMaximaStrategy> *hough, std::vector<int> indexesMax, std::vector < ContinuousVote<CodomainType> > & neighbourhoodContinuousVoteSpace)
		{
			CodomainType sum = 0;
			int nNeigh = 1;	// Total number of neighbours = 3^nDim
			for(int n = 0; n < hough->m_nDim; n++)
				nNeigh *= 3;

			for(int n = 0; n < nNeigh; n++)
			{
				int indexTot = 0;
				int exp = 1;
				int currNeighIndex = 0;

				for(int d = 0; d < hough->m_nDim; d++)
				{
					currNeighIndex = indexesMax[d] + ( n % (exp*3) ) / exp - 1; // (n % 3^(nDim+1) / 3^nDim) - 1
					if( currNeighIndex >= 0 && currNeighIndex <= hough->m_Bins[d]-1 )
					{
						indexTot += currNeighIndex * hough->m_partialBinProducts[d];
					}
					else
					{
						indexTot = -1;
						break;
					}
					exp *= 3;
				}

				if(indexTot != -1)
				{
					sum += hough->m_HoughSpace[indexTot];
					neighbourhoodContinuousVoteSpace.insert(neighbourhoodContinuousVoteSpace.end(), hough->m_ContinuousVoteSpace[indexTot].begin(), hough->m_ContinuousVoteSpace[indexTot].end());
				}
			}

			return sum;
		};

		template <typename CodomainType, typename AdditionalInfoType, typename VoteStrategy>
		CodomainType sumVotes(pstHough<CodomainType, AdditionalInfoType, VoteStrategy, pstHoughNeighborhoodFindMaximaStrategy> *hough, std::vector<int> indexesMax, std::vector < ContinuousVote<CodomainType> > & neighbourhoodContinuousVoteSpace, std::vector < AdditionalInfoType > &neighbourhoodAdditionalInfoSpace)
		{
			CodomainType sum = 0;
			int nNeigh = 1;	// Total number of neighbours = 3^nDim
			for(int n = 0; n < hough->m_nDim; n++)
				nNeigh *= 3;

			for(int n = 0; n < nNeigh; n++)
			{
				int indexTot = 0;
				int exp = 1;
				int currNeighIndex = 0;

				for(int d = 0; d < hough->m_nDim; d++)
				{
					currNeighIndex = indexesMax[d] + ( n % (exp*3) ) / exp - 1; // (n % 3^(nDim+1) / 3^nDim) - 1
					if( currNeighIndex >= 0 && currNeighIndex <= hough->m_Bins[d]-1 )
					{
						indexTot += currNeighIndex * hough->m_partialBinProducts[d];
					}
					else {
						indexTot = -1;
						break;
					}
					exp *= 3;
				}

				if(indexTot != -1)
				{
					sum += hough->m_HoughSpace[indexTot];
					neighbourhoodContinuousVoteSpace.insert(neighbourhoodContinuousVoteSpace.end(), hough->m_ContinuousVoteSpace[indexTot].begin(), hough->m_ContinuousVoteSpace[indexTot].end());
					neighbourhoodAdditionalInfoSpace.insert(neighbourhoodAdditionalInfoSpace.end(), hough->m_AdditionalInfoSpace[indexTot].begin(), hough->m_AdditionalInfoSpace[indexTot].end());
				}
			}

			return sum;
		};
	};


	//given the coordinates of the centroid of mesh_ref represented w.r.t. LRFs of mesh_ref (vCentroidWrtLRF), represents the coordinates w.r.t the LRFs of mesh_trg (vTransformedCentroidSpace)
	template<typename T> static void CreateTransformedCentroidSpace(const FeatureMatch* matches, const int numMatches, const Feature3D* feat3D_trg, const std::vector<float> &vCentroidWrtLRF, std::vector<float> &vTransformedCentroidSpace, std::vector<T> &vWeights, const bool useMatchScoreAsWeight )
	{

		vTransformedCentroidSpace.resize(numMatches*3);
		float* ptrTransformedCentroidSpace = &vTransformedCentroidSpace[0];

		vWeights.resize(numMatches);
		T* ptrWeights = &vWeights[0];

		int idx_ref;
		int idx_trg;
		FeatureMatch* ptrMatches = (FeatureMatch*)matches;
		const float* ptrCentroid;
		for(int ma = 0; ma < numMatches; ma++)
		{
			idx_ref = ptrMatches->refIndex;
			idx_trg = ptrMatches->trgIndex;

			ptrCentroid = &vCentroidWrtLRF[3*idx_ref];

			//change the basis of the centroid
			*ptrTransformedCentroidSpace++ = feat3D_trg[idx_trg].rf[0]*ptrCentroid[0] + feat3D_trg[idx_trg].rf[3]*ptrCentroid[1] + feat3D_trg[idx_trg].rf[6]*ptrCentroid[2] + feat3D_trg[idx_trg].x;
			*ptrTransformedCentroidSpace++ = feat3D_trg[idx_trg].rf[1]*ptrCentroid[0] + feat3D_trg[idx_trg].rf[4]*ptrCentroid[1] + feat3D_trg[idx_trg].rf[7]*ptrCentroid[2] + feat3D_trg[idx_trg].y;
			*ptrTransformedCentroidSpace++ = feat3D_trg[idx_trg].rf[2]*ptrCentroid[0] + feat3D_trg[idx_trg].rf[5]*ptrCentroid[1] + feat3D_trg[idx_trg].rf[8]*ptrCentroid[2] + feat3D_trg[idx_trg].z;


			//get weigth
			*ptrWeights = (useMatchScoreAsWeight) ? ((T)(ptrMatches->matchScore)) : 1;
			ptrWeights++;
			ptrMatches++;
		}
	}

	template<typename CodomainType, class VoteStrategy, class MaximaStrategy> static void PopulateHoughSpaceAndFindMaxima(pstHough<CodomainType, bool, VoteStrategy, MaximaStrategy >* houghSpace, const std::vector<float> &vVotes, const std::vector<CodomainType> &vWeights, std::vector< MaximumContainer<CodomainType> > &vMaxima )
	{
		houghSpace->Reset();

		//populate Hough space
		size_t nVotes = vWeights.size();
		float const* ptrVotes = &vVotes[0];
		CodomainType const* ptrWeights = &vWeights[0];
		double vote_double[3];
		for(size_t vo = 0; vo < nVotes; vo++)
		{
			vote_double[0] = ptrVotes[0];
			vote_double[1] = ptrVotes[1];
			vote_double[2] = ptrVotes[2];
			houghSpace->Vote(vote_double, *ptrWeights, (int)vo);
			ptrWeights++;
			ptrVotes+=3;
		}

		//find maximum
		vMaxima.resize(1);
		houghSpace->findMaximum(vMaxima[0]);
	}

	template<typename T> static void ObtainHoughMatches(std::vector< MaximumContainer<T> > & vMaxima, const FeatureMatch* matches, std::vector<FeatureMatch> &houghMatches, const int iMaxima = 0)
	{
		if(vMaxima.size() == 0)
		{
			houghMatches.clear();
			return;
		}
		int nHoughMatches = (int)vMaxima[iMaxima].continuousVote.size();
		houghMatches.resize(nHoughMatches);

		for(int ma=0; ma<nHoughMatches; ma++)
		{
			houghMatches[ma] = matches[ vMaxima[iMaxima].continuousVote[ma].id ];
		}	
	}


}


#endif



