#ifndef VPT_NORMAL_MAP
#define VPT_NORMAL_MAP

#include "Volume.h"
#include "TransferFunction.h"
#include "SegmentMap.h"

#define NORMAL_ITEXTURE_FORMAT GL_LUMINANCE
#define NORMAL_PIXEL_FORMAT GL_LUMINANCE
#define NORMAL_DATA_FORMAT GL_UNSIGNED_BYTE
#define normal_type_scale(x) vval_scale(x)

namespace vpt{

	typedef vval normal_type; 
	typedef VolumeData<normal_type> NormalVolume; 

	class Project; 

	class NormalMap{
	protected: 
		NormalVolume* _norm; // this is the current normal buffer
		NormalVolume* _original; // this is the original normal buffer
		NormalVolume* _buffer; // this is a buffer for working
		FVolume* _fbuffer; 

		FVolume* _volume; 
		Project* _project; 

	public: 
		NormalMap(FVolume* vol){
			_volume = vol; 

			uint dimx = vol->getDimX(); 
			uint dimy = vol->getDimY();
			uint dimz = vol->getDimZ();

			_norm = new NormalVolume(); 
			_norm->resize((uint)dimx,(uint)dimy,(uint)dimz); 

			_buffer = new NormalVolume(); 
			_buffer->resize((uint)dimx*2,(uint)dimy,(uint)dimz); 

			_original = new NormalVolume(); 
			_original->resize((uint)dimx,(uint)dimy,(uint)dimz); 

			_fbuffer = new FVolume(); 
			_fbuffer->resize((uint)dimx,(uint)dimy,(uint)(dimz));
		}

		~NormalMap(){
			delete _norm; 
			delete _buffer; 
			delete _original; 
			delete _fbuffer; 
		}

		inline NormalVolume* getNormal() { return _norm; }
		inline NormalVolume* getOriginalNormal() { return _original; }
		inline FVolume* getVolume() { return _volume; }
		inline NormalVolume* getBuffer() { return _buffer; }
		inline FVolume* getFBuffer() { return _fbuffer; }

		inline void setProject(Project* p){ _project = p; }
		inline Project* getProject() { return _project; }
	};
}

#endif