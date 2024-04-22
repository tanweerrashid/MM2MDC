#ifndef VPT_COMBINED_VOLUME_UTILS_H
#define VPT_COMBINED_VOLUME_UTILS_H

#include "CombinedVolume.h"

namespace vpt{
	class CombinedVolumeUtils{
	protected: 

	public: 
		template<typename T> 
		inline static osg::Image* convertToOSGImage(CombinedVolume<T>* cv){
			T t=0;

			GLenum type = GL_UNSIGNED_BYTE; 
			switch(sizeof(t)){
				case 1: type = GL_UNSIGNED_BYTE; break; 
				case 4: type = GL_FLOAT; break; 
				default: type = GL_UNSIGNED_BYTE; break; 
			}

			GLenum iformat = GL_LUMINANCE; 
			GLenum format = GL_LUMINANCE; 
			switch(cv->getNumVolumes()){
				case 2: 
					format = GL_LUMINANCE_ALPHA; 
					iformat = GL_LUMINANCE_ALPHA; 
					break;
				case 3: 
					format = GL_RGB; 
					iformat = GL_RGB; 
					break;
				case 4: 
					format = GL_RGBA; 
					iformat = GL_RGBA; 
					break; 
				case 1: 
				default: break;
			}
//			format = GL_LUMINANCE_ALPHA; 
//			iformat = GL_LUMINANCE_ALPHA; 
			osg::Image* ret = new osg::Image(); 
			ret->setImage(cv->getDimX(),cv->getDimY(),cv->getDimZ(),iformat,format,type,cv->data(),osg::Image::NO_DELETE);

			return ret;
		}
	};
}

#endif