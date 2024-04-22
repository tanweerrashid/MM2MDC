#ifndef VPT_COMBINED_VOLUME_H
#define VPT_COMBINED_VOLUME_H

#include <vector>

#include "Volume.h"

namespace vpt{
	template<typename T>
	class CombinedVolume{
	protected: 
		T* _data; 
		VolumeData<T>** _vols; 
		uint _numVols; 
		int _actualSize; 

	public: 
		static const uint BAD_STRIDE = 5; 

		CombinedVolume(const std::vector<VolumeData<T>*>& vols){
			_numVols = vols.size(); 
			if(_numVols>4 || _numVols==0){
				_data = NULL; 
				_vols = NULL; 
				return ; 
			}

			_vols = new VolumeData<T>*[_numVols]; 
			int nx = (int) vols[0]->getDimX(); 
			int ny = (int) vols[0]->getDimY();
			int nz = (int) vols[0]->getDimZ(); 

			for(uint j=0;j<_numVols;j++){
				_vols[j] = vols[j]; 

				int mx = (int) vols[j]->getDimX(); 
				int my = (int) vols[j]->getDimY(); 			
				int mz = (int) vols[j]->getDimZ(); 

				if(mx!=nx || my!=ny || mz!=nz){
					cerr<<"cannot have volumes of different dimensions"<<endl;
					delete [] _vols; 
					_data = NULL; 
					_vols = NULL; 
					return; 
				}

			}

			_data = new T[nx*ny*nz*_numVols];
			for(uint j=0;j<nx*ny*nz*_numVols;j++)
				_data[j] = 0;

			int flag = 0; 
			for(uint j=0;j<_numVols;j++)
				flag = flag | (1<<j);

			fillAll(flag);
		}


		~CombinedVolume(){
			delete [] _data; 
			delete [] _vols; 
		}
		inline void fillAll(int flag){
			uint dimx = getDimX(); 
			uint dimy = getDimY(); 
			uint dimz = getDimZ();

       //     #pragma omp parallel for
			for(uint p=0;p<_numVols;p++){
				if((1<<p)&flag){
					for(uint k=0;k<dimz;k++){
						for(uint j=0;j<dimy;j++){
							for(uint i=0;i<dimx;i++){
								// better for cache locality in case of parallel execution
								_data[((k*dimy+j)*dimx + i)*_numVols+p]
								= (*_vols[p])(i,j,k);
/*
								if(i+j+k>200 && i+j+k<225 && p==1){
									cout<<(int)(*_vols[p])(i,j,k)<<" "; 
								}
*/
							}
						}
					}
				}
			}
		}

		inline void fillRestAll(int flag){
			fillAll(~flag);
		}

		inline T* data() { return _data; }

		// use this sparingly
		inline uint getStride(VolumeData<T>* v) { 
			for(int j=0;j<_numVols;j++)
				if(_vols[j]==v) return j; 
			return BAD_STRIDE; 
		}

		inline int getNumVolumes() { return _numVols; }
		inline VolumeData<T>* getVolume(int i) { return _vols[i]; }

		inline T& operator()(
			const uint& x, const uint& y, const uint& z, 
			const uint& dimx, const uint& dimy, const uint& dimz, 
			const uint& st){
				return _data[((z*dimy + y)*dimx + x)*_numVols + st]; 
		}

		inline T operator()(
			const uint& x, const uint& y, const uint& z, 
			const uint& dimx, const uint& dimy, const uint& dimz, 
			const uint& st) const {
				return _data[((z*dimy + y)*dimx + x)*_numVols + st]; 
		}

		inline int getFlag(VolumeData<T>* v){
			for(int j=0;j<_numVols;j++){
				if(_vols[j]==v){
					return 1<<j;
				}
			}
			return 1<<5;
		}

		inline uint getDimX() { return _vols[0]->getDimX(); }
		inline uint getDimY() { return _vols[0]->getDimY(); }
		inline uint getDimZ() { return _vols[0]->getDimZ(); }
	};

}


#endif
