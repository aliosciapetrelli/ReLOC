#include "ReLOC_Hough.h"


using namespace std;



void ReLOC::GetCentroidsWrtLRFs(const std::vector<float> &vCentroids, const std::vector<Feature3D*> &vFeatures3D, const std::vector<int> &vNumFeatures3D, std::vector< std::vector<float> > &vvCentroidWrtLRF)
{
	int nClouds = (int)vFeatures3D.size();
	vvCentroidWrtLRF.resize(nClouds);

	const float* ptrCentroid = &(vCentroids[0]);

	float basisChange4x4[16];
	for(int cl=0; cl<nClouds; cl++)
	{
		const Feature3D* feat3D = vFeatures3D[cl];
		const int nFeat3D = vNumFeatures3D[cl];

		vvCentroidWrtLRF[cl].resize(nFeat3D*3);
		float* ptrCentroidWrtLRF = &(vvCentroidWrtLRF[cl][0]);

		for(int fe=0; fe<nFeat3D; fe++)
		{
			if(feat3D->rf[0] != feat3D->rf[0])
			{
				*ptrCentroidWrtLRF++ = std::numeric_limits<float>::quiet_NaN();
				*ptrCentroidWrtLRF++ = std::numeric_limits<float>::quiet_NaN();
				*ptrCentroidWrtLRF++ = std::numeric_limits<float>::quiet_NaN();
			}	
			else
			{
				//find change-of-basis matrix
				memcpy(basisChange4x4, feat3D->rf, sizeof(float)*3);
				memcpy(basisChange4x4 + 4, feat3D->rf+3, sizeof(float)*3);
				memcpy(basisChange4x4 + 8, feat3D->rf+6, sizeof(float)*3);
				basisChange4x4[3] = - basisChange4x4[0]*feat3D->x - basisChange4x4[1]*feat3D->y - basisChange4x4[2]*feat3D->z;
				basisChange4x4[7] = - basisChange4x4[4]*feat3D->x - basisChange4x4[5]*feat3D->y - basisChange4x4[6]*feat3D->z;
				basisChange4x4[11] = - basisChange4x4[8]*feat3D->x - basisChange4x4[9]*feat3D->y - basisChange4x4[10]*feat3D->z;
				
				//change the basis of the centroid
				*ptrCentroidWrtLRF++ = basisChange4x4[0]*ptrCentroid[0] + basisChange4x4[1]*ptrCentroid[1] + basisChange4x4[2]*ptrCentroid[2] + basisChange4x4[3];
				*ptrCentroidWrtLRF++ = basisChange4x4[4]*ptrCentroid[0] + basisChange4x4[5]*ptrCentroid[1] + basisChange4x4[6]*ptrCentroid[2] + basisChange4x4[7];
				*ptrCentroidWrtLRF++ = basisChange4x4[8]*ptrCentroid[0] + basisChange4x4[9]*ptrCentroid[1] + basisChange4x4[10]*ptrCentroid[2] + basisChange4x4[11];
			}

			feat3D++;
		}
		
		ptrCentroid+=3;
	}
}

