#ifndef VPT_SEGMENT_MAP_UTILS_H
#define VPT_SEGMENT_MAP_UTILS_H

#include "SegmentMap.h"
#include "TextureManager.h"
#include "Session.h"
#include <osg/Image>
#include <string>
#include <fstream> 
#include <sstream> 

namespace vpt{
	class SegmentMapUtils{
	public: 
		static void writeToFile(const std::string& filename, SegmentMap* sm){
			std::fstream file(filename.c_str(), std::fstream::out); 

			int num = sm->getNumSegments(); 
			file<<num<<std::endl; 
			for(int j=0;j<num;j++){
				SegmentMapEntry* ent = sm->getEntry(j); 
				osg::Vec4 color = ent->getColor(); 
				float state = ent->getTransparency();

				file<<ent->getName()<<";"<<ent->getSegmentCode()<<";"; 
				file<<color.x()<<";"<<color.y()<<";"<<color.z()<<";"<<state<<endl;
			}
			file.close();

			std::cout<<"saved "<<filename<<endl;
		}

		static SegmentMap* readFromFile(const std::string& filename){
			std::fstream file(filename.c_str(), std::fstream::in);

			if(file.fail()){
				cerr<<"cannot open file: "<<filename<<" as segmentmap"<<endl;
				return NULL; 
			}

			SegmentMap* ret = new SegmentMap(false); 

			int num; 
			file>>num; 
			std::string tmp; 
			getline(file,tmp); // get the newline at the end

			for(int j=0;j<num;j++){
				getline(file,tmp); 
				std::vector<std::string> ss = strSplit(tmp,';'); 

				SegmentMapEntry* ent = ret->createSegmentMapEntry(); 

				std::string name = ss[0]; 
				int code=0; 
				float r=0,g=0,b=0,state=0;
				code = strToNum(ss[1],code); 
				r = strToNum(ss[2],r); 
				g = strToNum(ss[3],g); 
				b = strToNum(ss[4],b); 
				state = strToNum(ss[5],state); 

				ent->setName(name); 
				ent->setSegmentCode(code); 
				ent->setState(state); 
				ent->setColor(osg::Vec4(r,g,b,state)); 
			}

			ret->setFileName(getFileName(filename));
			ret->setFilePath(getFilePath(filename));
			ret->setSavableState(SAVABLE_OPENED_FILE);

			file.close(); 

			return ret; 
		}

		static osg::Image* convertToOSGImage(SegmentMap* sm){
			osg::Image* ret = new osg::Image(); 

			GLfloat* data = new GLfloat[NUM_SEGMENTS*4]; 

			for(int j=0;j<NUM_SEGMENTS;j++){
				SegmentMapEntry* ent = sm->getEntry(j); 
				if(ent!=NULL){
					osg::Vec4 color = ent->getColor(); 

					data[j*4] = color.x(); 
					data[j*4+1] = color.y(); 
					data[j*4+2] = color.z(); 
					if(ent->getState()==vpt::HIDE_SEGMENT)
						data[j*4+3] = 0;
					else
						data[j*4+3] = ent->getTransparency(); 
				}
				else{
					data[j*4] = 0; 
					data[j*4+1] = 0; 
					data[j*4+2] = 0; 
					data[j*4+3] = 0; 
				}
			}

			ret->setImage(NUM_SEGMENTS,1,1,GL_RGBA,GL_RGBA,GL_FLOAT,(unsigned char*) data,osg::Image::USE_NEW_DELETE); 
			return ret; 
		}

		inline static void updateTexture(SegmentMap* sm){
			if(TextureManager::inst()->hasSegmentMapTex1D(sm)){
				osg::ref_ptr<osg::Texture1D> tex = TextureManager::inst()->getOrCreateSegmentMapTex1D(sm); 
				osg::ref_ptr<osg::Image> img = tex->getImage(); 

				GLfloat* data = new GLfloat[NUM_SEGMENTS*4]; 

				for(int j=0;j<NUM_SEGMENTS;j++){
					SegmentMapEntry* ent = sm->getEntry(j); 
					if(ent!=NULL){
						osg::Vec4 color = ent->getColor(); 

						data[j*4] = color.x(); 
						data[j*4+1] = color.y(); 
						data[j*4+2] = color.z(); 
						if(ent->getState()==vpt::HIDE_SEGMENT)
							data[j*4+3] = 0;
						else
							data[j*4+3] = ent->getTransparency(); 
					}
					else{
						data[j*4] = 0; 
						data[j*4+1] = 0; 
						data[j*4+2] = 0; 
						data[j*4+3] = 0; 
					}
				}

				img->setImage(NUM_SEGMENTS,1,1,GL_RGBA,GL_RGBA,GL_FLOAT,(unsigned char*) data,osg::Image::USE_NEW_DELETE); 
				img->dirty();
			}
		}
	}; 
}

#endif
