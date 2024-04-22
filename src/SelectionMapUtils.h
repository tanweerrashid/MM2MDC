#ifndef VPT_SELECTION_MAP_UTILS_H
#define VPT_SELECTION_MAP_UTILS_H

#include <osg/Shape>
#include <osg/Texture3D>
#include "SelectionMap.h"
#include "TextureManager.h"
#include "VolumeUtils.h"
#include "Common.h"

#include <vector>

namespace vpt{
	class SelectionMapUtils{
	protected: 
	public: 
		static inline osg::Image* convertToOSGImage(SelectionMap* sm){
			osg::Image* ret = new osg::Image(); 
			Volume* data = sm->getSelection();
			ret->setImage(data->getDimX(), data->getDimY(), data->getDimZ(), VV_ITEXTURE_FORMAT, VV_PIXEL_FORMAT, VV_DATA_FORMAT, (unsigned char*) data->data(),osg::Image::NO_DELETE);
			return ret;
		}

		static inline void updateTexture(SelectionMap* sm, uint xs, uint ys, uint zs, uint xdim, uint ydim, uint zdim, vval* ndata){
			if(TextureManager::inst()->hasSelectionMapTex3D(sm)){
				osg::ref_ptr<osg::Texture3D> tex = TextureManager::inst()->getOrCreateSelectionMapTex3D(sm); 
				if(tex->getSubloadCallback()==NULL)
					tex->setSubloadCallback( new VolumeSubloadCallback<vval>(VV_ITEXTURE_FORMAT, VV_PIXEL_FORMAT, VV_DATA_FORMAT)); 

				osg::ref_ptr<VolumeSubloadCallback<vval> > cb = dynamic_cast<VolumeSubloadCallback<vval>*>(tex->getSubloadCallback());
				cb->initSubload(xs,ys,zs,xdim,ydim,zdim,ndata); 
				//cout<<"finished "<<v<<endl;
			}
		}

		static void selectWithPlanes(SelectionMap* selectMap){
			Volume* select = selectMap->getSelection(); 
			//Volume* buffer = selectMap->getBuffer(); 
			int dimx = (int) select->getDimX(); 
			int dimy = (int) select->getDimY(); 
			int dimz = (int) select->getDimZ(); 
			float fdimx = dimx-1; 
			float fdimy = dimy-1; 
			float fdimz = dimz-1;

			int nplanes = 0;
			std::vector<osg::Vec3> centers; 
			std::vector<osg::Vec3> normals; 

			for(int j=0;j<selectMap->getNumPlanes();j++){
				if(selectMap->getPlane(j)->isActive()){
					centers.push_back(selectMap->getPlane(j)->getPoint()); 
					normals.push_back(selectMap->getPlane(j)->getNormal()); 
					nplanes++; 
				}
			}

			//float d = - (center*normal);
			if(nplanes>0){
                #pragma omp parallel for
				for(int k=0;k<dimz;k++){
					for(int j=0;j<dimy;j++){
						for(int i=0;i<dimx;i++){
							for(int pi=0;pi<nplanes;pi++){
								float vx = i/fdimx-centers[pi].x()-.5; 
								float vy = j/fdimy-centers[pi].y()-.5; 
								float ydot = vy*normals[pi].y(); 
								float vz = k/fdimz-centers[pi].z()-.5;
								float zdot = vz*normals[pi].z();

								float dot = vx*normals[pi].x() + ydot + zdot; 
								if(dot>EPS32)
									(*select)(i,j,k) = 255; 
								else{
									(*select)(i,j,k) = 0; 
									break;
								}
							}
						}
					}
				}
			}
			else{
				uint size = select->getSize();
				for(uint j=0;j<size;j++)
					(*select)(j)=255;
			}

			updateTexture(selectMap,0,0,0,(uint)dimx,(uint)dimy,(uint)dimz,select->data()); 
		}
	};

}

#endif