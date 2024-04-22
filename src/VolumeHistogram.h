#ifndef VOLUME_HISTOGRAM_H
#define VOLUME_HISTOGRAM_H

#include "Volume.h"

namespace vpt{
	class VolumeHistogram{
	protected: 
		int* _distr; 
		int _res; 
		int _max; 

	public: 
		VolumeHistogram(const Volume* v, int resolution=80){
			_distr = new int[resolution]; 
			_res = resolution; 

			for(int j=0;j<resolution;j++){
				_distr[j] = 0; 
			}

			int dimx = v->getDimX(); 
			int dimy = v->getDimY(); 
			int dimz = v->getDimZ(); 

			int resm1 = _res-1; 
			_max = 0; 

            #pragma omp parallel for
			for(int k=0;k<dimz;k++){
				for(int j=0;j<dimy;j++){
					for(int i=0;i<dimx;i++){
						int d = resm1*vval_inverse_scale((*v)(i,j,k)); 
						_distr[d]++; 
					}
				}
			}
			for(int j=0;j<resolution;j++)
				_max = max(_distr[j],_max);
		}
		~VolumeHistogram() { delete [] _distr; }

		inline int getResolution() { return _res; }
		inline int getDistrAt(int i) { 
			if(_distr[i]<0) return 0; 
			return _distr[i]; } 
		inline int getMaxCount() { return _max; }
	}; 
}

#endif