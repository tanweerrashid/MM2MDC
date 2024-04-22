#ifndef VPT_VOLUME_UTILS_H
#define VPT_VOLUME_UTILS_H

#include <fstream> 
#include <iostream> 
#include <stdio.h>
#include <string> 
#include "Volume.h"
#include "VolumeHistogram.h"
#include "mrc.h"
#include "Common.h"
//#include "Data/VolumeHistogram.h"

#include <osg/Image> 
#include <osg/Texture3D>

namespace vpt{

#define VOLUME_X_AXIS 0
#define VOLUME_Y_AXIS 1
#define VOLUME_Z_AXIS 2

	class BoundaryMap; 

	class VolumeUtils{
	protected: 
		static bool readVolume(const std::string& filename, Volume& vol, bool scale); 
		static bool readVolume(const std::string& filename, FVolume& vol, bool scale); 
	public: 
		// returns true if volume was successfully read
		static Volume* readByteVolumeFromFile(const std::string& filename, bool scale=true){
			Volume* ret = new Volume();
			bool res = VolumeUtils::readVolume(filename,*ret, scale);
			ret->setFileName(getFileName(filename));
			ret->setFilePath(getFilePath(filename));
			ret->setSavableState(vpt::SAVABLE_OPENED_FILE);
			if(!res){
				delete ret; 
				return NULL; 
			}
			return ret;
		} 

		static FVolume* readFloatVolumeFromFile(const std::string& filename, bool scale=true){
			FVolume* ret = new FVolume();
			bool res = VolumeUtils::readVolume(filename,*ret, scale);
			ret->setFileName(getFileName(filename));
			ret->setFilePath(getFilePath(filename));
			ret->setSavableState(vpt::SAVABLE_OPENED_FILE);
			if(!res){
				delete ret; 
				return NULL; 
			}
			return ret;
		}


		static void writeToFile(const std::string& filename, Volume* vol, VolumeIOType type=MRC_TYPE);
		static void writeToFile(const std::string& filename, FVolume* vol, VolumeIOType type=MRC_TYPE);

		static osg::Image* convertToOSGImage(Volume& vol); 
		static osg::Image* convertToOSGImage(FVolume& vol); 
		static void updateTexture(Volume* v, uint xs, uint ys, uint zs, uint xdim, uint ydim, uint zdim, vval* ndata); 
		static osg::Vec3 smoothGradient(Volume* v, const osg::Vec3& pos); 

		static inline VolumeHistogram* buildHistogram(Volume* vin){
			return new VolumeHistogram(vin); 
		}

		static void writeSegmentedVolToMRCFile(const std::string& filename, FVolume* vol, Volume* mask, SegmentMap* segmap, VolumeIOType type=MRC_TYPE); 

		// grow the mask border by pixels to make it look better
		static void growMask(Volume* mask, Volume* vol, TransferFunction* tf, vval segment, vval baseSegment); 

		template<typename T> 
		static VolumeData<T>* make2DSlice(VolumeData<T>* vol, int slice, int axis); 

		template<typename T> 
		inline static void copyVolume(VolumeData<T>* from, VolumeData<T>* to); 

		template<typename T> 
		static mrcH createMRCHeader(VolumeData<T>* vol, int modes=MODE_char); 

		static void subsample(FVolume* vol, Volume* mask, FVolume* nvol, Volume* nmask, int tx, int ty, int tz);

		static void buildPerMaterialFunction(FVolume* vol, Volume* mask, vval mat, FVolume* nvol, int ntx, int nty, int ntz); 
		static void buildMaterialFunctions(SegmentMap* smap, BoundaryMap* bm, FVolume* vol, Volume* mask, int ntx, int nty, int ntz); 
	};

	// This probably should not be here
	template<typename T>
	class VolumeSubloadCallback : public osg::Texture3D::SubloadCallback{
	protected: 
		bool _subload; 
		uint _xs, _ys, _zs, _xdim, _ydim, _zdim;
		T* _ndata; 
		GLenum _iformat, _format, _type; 
	public: 
		VolumeSubloadCallback() { 
			_subload = false; 
			_iformat = VV_ITEXTURE_FORMAT; 
			_format = VV_PIXEL_FORMAT;
			_type = VV_DATA_FORMAT; 
		}

		VolumeSubloadCallback(GLenum iformat, GLenum format, GLenum type){
			_subload = false; 
			_iformat = iformat; 
			_format = format; 
			_type = type; 
		}

		inline void initSubload(uint xs, uint ys, uint zs, uint xdim, uint ydim, uint zdim, T* data){
			_xs = xs; 
			_ys = ys; 
			_zs = zs; 
			_xdim = xdim; 
			_ydim = ydim; 
			_zdim = zdim; 
			_ndata = data; 
			_subload = true; 
		}

		inline void finishedSubload(){
			_subload = false; 
			_ndata = NULL; 
		}
		void load(const osg::Texture3D& tex, osg::State& state) const{
			const osg::Image* img = tex.getImage(); 
			tex.getExtensions(state.getContextID(),true)->glTexImage3D(GL_TEXTURE_3D,0,_iformat,img->s(),img->t(),img->r(),0,_format,_type,img->data());
		}

		void subload(const osg::Texture3D& tex, osg::State& state) const{
			if(_subload && _ndata){
				tex.getExtensions(state.getContextID(),true)->glTexSubImage3D(GL_TEXTURE_3D, state.getContextID(), _xs, _ys, _zs, _xdim, _ydim, _zdim, _format, _type, _ndata);
				bool* sb = const_cast<bool*>(&_subload); 
				(*sb) = false; 
			}
		}
	};

	template<typename T> 
	VolumeData<T>* VolumeUtils::make2DSlice(VolumeData<T>* vol, int slice, int axis){
		VolumeData<T>* ret = new VolumeData<T>(); 

		int dim = 0; 
		int rbound=0, cbound=0; 
		int xOrigDim = vol->getDimX(); 
		int yOrigDim = vol->getDimY(); 
		int zOrigDim = vol->getDimZ(); 

		int rstep[3] = {0,0,0};
		int cstep[3] = {0,0,0};

		int x=0,y=0,z=0;

		switch(axis){
		case VOLUME_X_AXIS: 
			dim = vol->getDimX(); 
			rbound = vol->getDimZ(); 
			cbound = vol->getDimY(); 
			rstep[2]++; 
			cstep[1]++; 
			x=slice;
			break;
		case VOLUME_Y_AXIS: 
			dim = vol->getDimY(); 
			rbound = vol->getDimX(); 
			cbound = vol->getDimZ(); 
			rstep[0]++; 
			cstep[2]++; 
			y=slice;
			break;
		case VOLUME_Z_AXIS: 
			dim = vol->getDimZ(); 
			rbound = vol->getDimX(); 
			cbound = vol->getDimY(); 
			rstep[0]++; 
			cstep[1]++; 
			z=slice;
			break; 
		}

		if(slice>=dim) return NULL; 

		ret->resize(rbound,cbound,1);

		for(int k=0;k<cbound;k++){
			for(int j=0;j<rbound;j++){
				(*ret)(j,k,0) = (*vol)(x,y,z); 
				x+=rstep[0]; 
				y+=rstep[1]; 
				z+=rstep[2]; 
			}
			if(x>=xOrigDim) x=0;
			if(y>=yOrigDim) y=0;
			if(z>=zOrigDim) z=0;

			x+=cstep[0]; 
			y+=cstep[1]; 
			z+=cstep[2]; 
		}

		return ret; 
	}

	// need to think of a better way to handle this 
	template<typename T> 
	mrcH VolumeUtils::createMRCHeader(VolumeData<T>* vol, int modes){
		mrcH ret; 

		ret.nc = (int) vol->getDimX(); 
		ret.nr = (int) vol->getDimY(); 
		ret.ns = (int) vol->getDimZ(); 

		ret.mode = modes; 
		ret.ncstart = 0; 
		ret.nrstart = 0; 
		ret.nsstart = 0; 

		ret.mx = ret.nc; 
		ret.my = ret.nr; 
		ret.mz = ret.ns; 

		ret.xlen = 1;
		ret.ylen = 1;
		ret.zlen = 1;

		ret.alpha = 90; 
		ret.beta = 90; 
		ret.gamma = 90; 

		ret.mapc = 1; 
		ret.mapr = 2; 
		ret.maps = 3; 

		float tmin=1e9, tmax=-1e9;
		int size = vol->getSize(); 

        #pragma omp parallel for 
		for(int j=0;j<size;j++){
			tmin = min(tmin,(float)(*vol)(j));
			tmax = max(tmax,(float)(*vol)(j)); 
		}

		ret.amin = tmin; 
		ret.amax = tmax; 
		ret.amean = (tmax-tmin)/2; 

		ret.ispg = 0; 
		ret.nsymbt = 0; 

		ret.xorigin = 0; 
		ret.yorigin = 0; 
		ret.zorigin = 0; 

		ret.map[0] = 'M'; 
		ret.map[1] = 'A'; 
		ret.map[2] = 'P'; 
		ret.map[3] = ' '; 

		ret.machinestamp = 0; 
		ret.rms = 0; 
		ret.nlabl = 0; 
		for(int j=0;j<MRC_NUM_LABELS;j++)
			for(int k=0;k<MRC_LABEL_SIZE;k++)
				ret.labels[j][k] = 0; 

		return ret; 
	}

	template<typename T>
	void VolumeUtils::copyVolume(VolumeData<T>* from, VolumeData<T>* to){
		if(from->getDimX()!=to->getDimX() ||
			from->getDimY()!=to->getDimY() || 
			from->getDimZ()!=to->getDimZ()){
				cerr<<"trying to copy volumes of different sizes"<<endl;
				return; 
		}

		for(uint k=0;k<from->getDimZ();k++)
			for(uint j=0;j<from->getDimY();j++)
				for(uint i=0;i<from->getDimX();i++){
					(*to)(i,j,k) = (*from)(i,j,k);
				}
	}
}

#endif
