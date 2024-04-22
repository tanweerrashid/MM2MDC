#ifndef VPT_BOUNDARY_MAP_H
#define VPT_BOUNDARY_MAP_H

#include "Volume.h"

namespace vpt{
/**/
	typedef unsigned char bmvval; // BM Volume type
	#define BMV_TEXTURE_TYPE GL_UNSIGNED_BYTE
	#define BMV_INTERNAL_TEXTURE_TYPE GL_RGBA
    #define bmv_scale(a) (bmvval)(255*a)

/*
	typedef float bmvval; // BM Volume type
	#define BMV_TEXTURE_TYPE GL_FLOAT
	#define BMV_INTERNAL_TEXTURE_TYPE GL_RGBA_FLOAT32
    #define bmv_scale(a) (bmvval)(a)
/**/

	class Project; 

	class BoundaryMap{
	protected: 
		Volume* _boundary; 
		Volume* _auxBoundary; 
		Project* _project; 

		FVolume* _function; 
		FVolume* _auxFunction; 

		VolumeData<vval>* _corners;
		VolumeData<bmvval>* _packedVol; 

		FVolume* _perSegFunction[NUM_SEGMENTS]; 

		int _dimx,_dimy,_dimz; 
	public: 

		BoundaryMap(){
			_boundary = NULL; 
			_auxBoundary = NULL; 
			_function = NULL; 
			_auxFunction = NULL; 
			_corners = NULL; 
			_packedVol = NULL; 
			_dimx = 0; 
			_dimy = 0; 
			_dimz = 0; 
			_project = NULL; 

			for(int j=0;j<NUM_SEGMENTS;j++)
				_perSegFunction[j] = NULL; 
		}

		BoundaryMap(Volume* vol){
			_boundary = NULL; 
			_auxBoundary = NULL; 
			_function = NULL; 
			_auxFunction = NULL; 
			_corners = NULL; 
			_packedVol = NULL; 
			_dimx = 0; 
			_dimy = 0; 
			_dimz = 0; 
			_project = NULL; 

			int x = vol->getDimX();
			int y = vol->getDimY();
			int z = vol->getDimZ();
			resize(x, y, z); 

			
			for(int j=0;j<NUM_SEGMENTS;j++)
				_perSegFunction[j] = NULL; 
		}

		void resize(int dimx, int dimy, int dimz){
			if(_boundary) delete _boundary; 
			_boundary = new Volume(); 
			_boundary->resize(dimx,dimy,dimz); 

			if(_auxBoundary) delete _auxBoundary; 
			_auxBoundary = new Volume(); 
			_auxBoundary->resize(dimx,dimy,dimz); 

			if(_function) delete _function; 
			_function = new VolumeData<float>(); 
			_function->resize(dimx, dimy, dimz); 

			if(_auxFunction) delete _auxFunction;
			_auxFunction = new VolumeData<float>(); 
			_auxFunction->resize(dimx, dimy, dimz); 

			if(_corners) delete _corners; 
			_corners = new VolumeData<unsigned char>(); 
			_corners->resize((dimx-1)*8,dimy-1,dimz-1);

			if(_packedVol) delete _packedVol; 
			_packedVol = new VolumeData<bmvval>(); 
			_packedVol->resize((dimx-1)*8,dimy-1,dimz-1);

			_dimx = _corners->getDimX()/8; 
			_dimy = _corners->getDimY(); 
			_dimz = _corners->getDimZ(); 

			//getFunction()->setAll(1.5f);
                        //getAuxFunction()->setAll(1.5f);
                        getFunction()->setAll(1.0f);
                        getAuxFunction()->setAll(1.0f);			
		}

		~BoundaryMap() { 
			delete _boundary; 
			delete _auxBoundary; 
			delete _corners; 
			delete _packedVol;
			delete _function; 
			delete _auxFunction; 
		}

		template<typename T>
		inline void writeVolume(VolumeData<T>* v){
			int dimx = v->getDimX(); 
			int dimy = v->getDimY(); 
			int dimz = v->getDimZ(); 

            #pragma omp parallel for
			for(int k=0;k<dimz;k++)
				for(int j=0;j<dimy;j++)
					for(int i=0;i<dimx;i++)
						writeVolumeVal(i,j,k,(*v)(i,j,k)); 
		}

		template<typename T>
		inline void writeVolumeVal(int x, int y, int z, T v){
			//x+=20; 
			for(int j=0;j<8;j++){
				int nx=x-1,ny=y-1,nz=z-1;
				if(j&1)	nx++; 
				if(j&2)	ny++; 
				if(j&4) nz++; 
				if(nx<0 || ny<0 || nz<0 || nx>=_dimx || ny>=_dimy || nz>=_dimz) continue; 

				int nj = (~j) & 7; 
				(*_packedVol)(nx*8+nj,ny,nz) = (bmvval) v; 
			}
		}

		template<typename T>
		inline void writeVolumeVal(int x, int y, int z, int dx, int dy, int dz, T v){
			//x+=20; 
			for(int j=0;j<8;j++){
				int nx=x-1,ny=y-1,nz=z-1;
				if(j&1)	nx++; 
				if(j&2)	ny++; 
				if(j&4) nz++; 
				if(nx<0 || ny<0 || nz<0 || nx>=dx|| ny>=dy || nz>=dz) continue; 

				int nj = (~j) & 7; 
				int ind = nx*8+nj + (ny + (nz*dy))*dx*8; 
				(*_packedVol)(ind) = (bmvval) v; 
			}
		}

		inline void writeCornerMask(int x, int y, int z, vval v){
			for(int j=0;j<8;j++){
				int nx=x-1,ny=y-1,nz=z-1;
				if(j&1)	nx++; 
				if(j&2)	ny++; 
				if(j&4) nz++; 
				if(nx<0 || ny<0 || nz<0 || nx>=_dimx || ny>=_dimy || nz>=_dimz) continue; 

				int nj = (~j) & 7; 
				(*_corners)(nx*8+nj,ny,nz) = v; 
			}

			for(int j=0;j<8;j++){
				int nx=x-1,ny=y-1,nz=z-1;
				if(j&1)	nx++; 
				if(j&2)	ny++; 
				if(j&4) nz++; 
				if(nx<0 || ny<0 || nz<0 || nx>=_dimx || ny>=_dimy || nz>=_dimz) continue; 

				vval c = (*_corners)(nx*8,ny,nz)%100; 
				if(isHomogenous(nx,ny,nz))
					for(int k=0;k<8;k++)
						(*_corners)(nx*8+k,ny,nz)=c+100; 
				else
					for(int k=0;k<8;k++)
						(*_corners)(nx*8+k,ny,nz)%=100;
			}
		}

		inline void writeCornerMask(int x, int y, int z, int dx, int dy, int dz, vval v){
			for(int j=0;j<8;j++){
				int nx=x-1,ny=y-1,nz=z-1;
				if(j&1)	nx++; 
				if(j&2)	ny++; 
				if(j&4) nz++; 
				if(nx<0 || ny<0 || nz<0 || nx>=dx || ny>=dy || nz>=dz) continue; 

				int nj = (~j) & 7; 
				int ind = nx*8+nj + (ny + (nz*dy))*dx*8; 
				(*_corners)(ind) = v; 
				//for(int h=0;h<8;h++)
				//(*_corners)(nx*8+h,ny,nz) = ((*_corners)(nx*8+h,ny,nz)%100);
			}

			for(int j=0;j<8;j++){
				int nx=x-1,ny=y-1,nz=z-1;
				if(j&1)	nx++; 
				if(j&2)	ny++; 
				if(j&4) nz++; 
				if(nx<0 || ny<0 || nz<0 || nx>=dx || ny>=dy || nz>=dz) continue; 

				int ind = nx*8 + (ny + (nz*dy))*dx*8; 
				vval c = (*_corners)(ind)%100; 
				if(isHomogenous(ind))
					for(int k=0;k<8;k++)
						(*_corners)(ind+k)=c+100; 
				else
					for(int k=0;k<8;k++)
						(*_corners)(ind+k)%=100;
			}
		}

		inline bool isHomogenous(int ind){
			vval c = (*_corners)(ind)%100; 
			for(int k=1;k<8;k++){
				if(c!=(*_corners)(ind+k)%100){
					return false; 
				}
			}
			return true; 
		}

		inline bool isHomogenous(int nx, int ny, int nz){
			vval c = (*_corners)(nx*8,ny,nz)%100; 
			for(int k=1;k<8;k++){
				if(c!=(*_corners)(nx*8+k,ny,nz)%100){
					return false; 
				}
			}
			return true; 
		}

		Volume* getBoundary() { return _boundary; }
		Volume* getAuxBoundary() { return _auxBoundary; }

		void setFunction(VolumeData<float>* func) { _function = func; }
		VolumeData<float>* getFunction() { return _function; }
		VolumeData<float>* getAuxFunction() { return _auxFunction; }

		VolumeData<vval>* getCorners() { return _corners; }
		VolumeData<bmvval>* getPackedVolume() { return _packedVol; }

		inline void setProject(Project* p) { _project = p; }
		inline Project* getProject() { return _project; }

		inline FVolume* getFunction(int segment){
			return _perSegFunction[segment]; 
		}

		inline void setFunction(int segment, FVolume* fvol){
			_perSegFunction[segment] = fvol; 
		}
	};
}

#endif
