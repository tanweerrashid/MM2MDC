#ifndef VPT_TRANSFER_FUNCTION_UTILS_H
#define VPT_TRANSFER_FUNCTION_UTILS_H

#include "TransferFunction.h"
#include "TextureManager.h"
#include "Session.h"

#include <osg/Image>
#include <iostream>
#include <fstream> 
#include <string> 

namespace vpt{

	class TransferFunctionUtils{
	private: 
		inline static vval* discretizeTransFunc(TransferFunction* tf, int resolution){
			vval* data = new vval[resolution]; 

			for(int j=0;j<resolution;j++)
				data[j] = vval_scale(tf->sample(j/1.f/(resolution-1))); 
			return data; 
		}

	public: 
		static void writeToFile(const std::string& filename, TransferFunction* transFunc){
			std::fstream file(filename.c_str(), std::fstream::out); 

			int num = transFunc->getNumDataPoints(); 

			file<<num<<endl;
			for(int j=0;j<num;j++){
				TFDPoint* pt = transFunc->getDataPoint(j); 
				file<<pt->getIsoValue()<<" "<<pt->getAlpha()<<endl;
			}

			file.close(); 
		}

		static TransferFunction* readFromFile(const std::string& filename){
			std::fstream file(filename.c_str(), std::fstream::in); 
			if(file.fail()){
				cerr<<"cannot read: "<<filename<<endl;
				return NULL; 
			}

			TransferFunction* ret = new TransferFunction(false); 
			int num; 
			file>>num; 

			for(int j=0;j<num;j++){
				float alpha, isovalue; 
				file>>isovalue>>alpha; 
				ret->createDataPoint(isovalue,alpha);
			}

			ret->setFileName(getFileName(filename));
			ret->setFilePath(getFilePath(filename));
			ret->setSavableState(SAVABLE_OPENED_FILE);

			file.close();

			return ret; 
		}

		inline static osg::Image* convertToOSGImage(TransferFunction* tf){
			osg::Image* ret = new osg::Image(); 
			int resolution = TF_RESOLUTION; 
			vval* data = discretizeTransFunc(tf, resolution);
			ret->setImage(resolution,1,1,VV_ITEXTURE_FORMAT,VV_PIXEL_FORMAT,VV_DATA_FORMAT,data,osg::Image::USE_NEW_DELETE); 

			return ret; 
		}

		inline static void updateTexture(TransferFunction* tf){
			if(TextureManager::inst()->hasTransFuncTex1D(tf)){
				osg::ref_ptr<osg::Texture1D> tex = TextureManager::inst()->getOrCreateTransFuncTex1D(tf); 
				osg::ref_ptr<osg::Image> img = tex->getImage(); 

				int resolution = TF_RESOLUTION; 
				vval* data = discretizeTransFunc(tf, resolution);
				GLenum iformat = img->getInternalTextureFormat(); 
				GLenum pformat = img->getPixelFormat(); 
				GLenum type = img->getDataType(); 

				img->setImage(resolution,1,1,iformat,pformat,type,data,osg::Image::USE_NEW_DELETE); 
				img->dirty();
			}
		}
	}; 

}

#endif