#ifndef VPT_NORMAL_MAP_UTILS_H
#define VPT_NORMAL_MAP_UTILS_H

#include "Volume.h"
#include "TextureManager.h"
#include "NormalMap.h"
#include "VolumeUtils.h"
#include "BoundaryMapUtils.h"
#include "Session.h"

#include <iostream>

namespace vpt{
#define BAD_VAL -1e9

	class NormalMapUtils{
	protected: 

		static inline float getVal(FVolume* v, Volume* m, SegmentMap* sg, TransferFunction* tf, 
			int x, int y, int z, int maxx, int maxy, int maxz){
			if(x>maxx || x<0 || y>maxy || y<0 || z>maxz || z<0) return 0; 


			//if(sg->getEntryFast((*m)(x,y,z))->getState()!=vpt::HIDE_SEGMENT){
				//return tf->sampleCoarse((*v)(x,y,z)); 
				return (*v)(x,y,z); 
//			}
//			return 0; 
		}

		static inline float getValWrtMask(FVolume* v, Volume* m, TransferFunction* tf, bool* show,
			int x, int y, int z, int maxx, int maxy, int maxz, FVolume* fbuffer){
			if(x>maxx || x<0 || y>maxy || y<0 || z>maxz || z<0) return 0; 
/*
			if(!show[(*m)(x,y,z)]){
				return -(*v)(x,y,z); 
			}
			else if(true)
				return (*v)(x,y,z); 
*/


			if(fbuffer && (*fbuffer)(x,y,z)>BAD_VAL+1)
				return (*fbuffer)(x,y,z); 

			float normal_mask[125]={
				1,1,1,1,1,
				1,1,3,1,1,
				1,3,3,3,1,
				1,1,3,1,1,
				1,1,1,1,1,

				1,1,3,1,1,
				1,3,3,3,1,
				3,3,4,3,3,
				1,3,3,3,1,
				1,1,3,1,1,

				1,3,3,3,1,
				3,3,4,3,3,
				3,4,5,4,3,
				3,3,4,3,3,
				1,3,3,3,1,

				1,1,3,1,1,
				1,3,3,3,1,
				3,3,4,3,3,
				1,3,3,3,1,
				1,1,3,1,1,

				1,1,1,1,1,
				1,1,3,1,1,
				1,3,3,3,1,
				1,1,3,1,1,
				1,1,1,1,1
			}; 

			float total = 247;
			int dd = 2; 
			int mdim = (2*dd+1); 
			int mdim2 = mdim*mdim;

			float samp = (float) (*v)(x,y,z); 
			if(samp==0)
				return 0; 

			float acc=0; 

			int ist = x>=dd ? -dd : 0; 
			int ied = x<=maxx-dd ? dd : 0; 
			int jst = y>=dd ? -dd : 0; 
			int jed = y<=maxy-dd ? dd : 0; 
			int kst = z>=dd ? -dd: 0; 
			int ked = z<=maxz-dd ? dd : 0; 


			for(int k=kst;k<=ked;k++){
				for(int j=jst;j<=jed;j++){
					for(int i=ist;i<=ied;i++){
						int nx = x+i, ny = y+j, nz = z+k; 
						if(show[(*m)(nx,ny,nz)])
							acc+=normal_mask[(k+dd)*mdim2 + (j+dd)*mdim + (i+dd)]/total;
					}
				}					
			}
//			cout<<"acc: "<<acc<<endl;
			if(fbuffer)
				return (*fbuffer)(x,y,z) = samp * acc; 
			return samp * acc; 
		}

		static inline float getValWrtDFunction(FVolume* v, Volume* m, TransferFunction* tf, bool* show,
			int x, int y, int z, int maxx, int maxy, int maxz, FVolume* fbuffer){
			if(x>maxx || x<0 || y>maxy || y<0 || z>maxz || z<0) return 0; 

			//if(fbuffer && (*fbuffer)(x,y,z)>BAD_VAL+1)
			//return (*fbuffer)(x,y,z); 

			if(show[(*m)(x,y,z)])
				return (*v)(x,y,z); 
			return -(*v)(x,y,z);

/*
			if(fbuffer)
				return (*fbuffer)(x,y,z) = samp * acc; 
			return samp * acc; 
*/
		}

	public: 
		static inline osg::Image* convertToOSGImage(NormalMap* nm){
			osg::Image* ret = new osg::Image(); 
			NormalVolume* data = nm->getNormal(); 
			ret->setImage(data->getDimX(), data->getDimY(), data->getDimZ(), NORMAL_ITEXTURE_FORMAT, NORMAL_PIXEL_FORMAT, NORMAL_DATA_FORMAT, (unsigned char*) data->data(),osg::Image::NO_DELETE);
			return ret; 
		}

		static inline void updateTexture(NormalMap* nm, uint xs, uint ys, uint zs, uint xdim, uint ydim, uint zdim, normal_type* ndata){
			if(TextureManager::inst()->hasNormalMapTex3D(nm)){
				osg::ref_ptr<osg::Texture3D> tex = TextureManager::inst()->getOrCreateNormalMapTex3D(nm); 
				if(tex->getSubloadCallback()==NULL){
					GLenum iformat = tex->getImage()->getInternalTextureFormat(); 
					GLenum pformat = tex->getImage()->getPixelFormat(); 
					GLenum type = tex->getImage()->getDataType(); 
					tex->setSubloadCallback( new VolumeSubloadCallback<normal_type>(iformat, pformat, type)); 
				}

				osg::ref_ptr<VolumeSubloadCallback<normal_type> > cb = dynamic_cast<VolumeSubloadCallback<normal_type>*>(tex->getSubloadCallback());
				cb->initSubload(xs,ys,zs,xdim,ydim,zdim,ndata); 
				//cout<<"finished "<<v<<endl;
			}
		}

		template<typename T>
		static inline void blur(VolumeData<T>* v, Volume* border, int xst, int yst, int zst, int xed, int yed, int zed){
			VolumeData<T>* tmp = new VolumeData<T>(); 
			tmp->resize(v->getDimX(),v->getDimY(),v->getDimZ());
/**/
			float normal_mask[125]={
				1,1,1,1,1,
				1,1,3,1,1,
				1,3,3,3,1,
				1,1,3,1,1,
				1,1,1,1,1,

				1,1,3,1,1,
				1,3,3,3,1,
				3,3,4,3,3,
				1,3,3,3,1,
				1,1,3,1,1,

				1,3,3,3,1,
				3,3,4,3,3,
				3,4,5,4,3,
				3,3,4,3,3,
				1,3,3,3,1,

				1,1,3,1,1,
				1,3,3,3,1,
				3,3,4,3,3,
				1,3,3,3,1,
				1,1,3,1,1,

				1,1,1,1,1,
				1,1,3,1,1,
				1,3,3,3,1,
				1,1,3,1,1,
				1,1,1,1,1
			}; 
			float total = 247;
			int dd = 2; 
			int mdim = (2*dd+1); 
			int mdim2 = mdim*mdim;
/**/
			int maxx = v->getDimX()-1;
			int maxy = v->getDimY()-1;
			int maxz = v->getDimZ()-1;

            #pragma omp parallel for
			for(int z=zst;z<=zed;z++){
				for(int y=yst;y<=yed;y++){
					for(int x=xst;x<=xed;x++){
						if((*border)(x,y,z)){
							int ist = x>=dd ? -dd : 0; 
							int ied = x<=maxx-dd ? dd : 0; 
							int jst = y>=dd ? -dd : 0; 
							int jed = y<=maxy-dd ? dd : 0; 
							int kst = z>=dd ? -dd: 0; 
							int ked = z<=maxz-dd ? dd : 0; 
							float acc=0,t2=0; 

							for(int k=kst;k<=ked;k++){
								for(int j=jst;j<=jed;j++){
									for(int i=ist;i<=ied;i++){
										int nx = x+i, ny = y+j, nz = z+k; 
										if((*border)(nx,ny,nz) && (*v)(nx,ny,nz)>0){
											acc+=normal_mask[(k+dd)*mdim2 + (j+dd)*mdim + (i+dd)]*(*v)(nx,ny,nz); 
											t2+=normal_mask[(k+dd)*mdim2 + (j+dd)*mdim + (i+dd)]; 
										}
									}
								}					
							}
							if(t2>0)
								(*tmp)(x,y,z) = (T) (acc/t2); 
						}
					}
				}
			}

            #pragma omp parallel for
			for(int z=zst;z<=zed;z++){
				for(int y=yst;y<=yed;y++)
					for(int x=xst;x<=xed;x++)
						if((*border)(x,y,z)){
							(*v)(x,y,z) = (*tmp)(x,y,z); 
							//(*v)(x,y,z) = 255;
						}
			}

			delete tmp; 
		}

		static inline void updateNormals(NormalMap* nm, const osg::Vec3& light,
			int xst, int yst, int zst, int xed, int yed, int zed, vval segcode, bool updateOrig=false, bool updateBySegment=false){

	//			if(!updateOrig) return;

				NormalVolume* norm = nm->getNormal();  
				NormalVolume* orig = nm->getOriginalNormal();
				FVolume* v = nm->getVolume(); // this is the un-truncated Volume
				Volume* m = nm->getProject()->getMask();
				Volume* tv = nm->getProject()->getVolume();
				normal_type* buffer = nm->getBuffer()->data();
				SegmentMap* segMap = nm->getProject()->getSegmentMap();
				TransferFunction* tf = nm->getProject()->getTransFunc();
				FVolume* fbuffer = nm->getFBuffer(); 
				fbuffer->setAll(xst,yst,zst,xed,yed,zed,BAD_VAL);

				FVolume* dfunc = nm->getProject()->getBoundaryMap()->getFunction();
				//FVolume* dfunc = nm->getVolume();

				if(!updateOrig)
					BoundaryMapUtils::updateBoundary(nm->getProject()->getBoundaryMap(), xst, yst, zst, xed, yed, zed);

				Volume* border = nm->getProject()->getBoundaryMap()->getBoundary(); 

				osg::Vec3 light1(100,20,0); 
				osg::Vec3 light2(20,100,0); 
				osg::Vec3 light3(0,5,100); 


				float maxx = (int) (v->getDimX()-1); 
				float maxy = (int) (v->getDimY()-1); 
				float maxz = (int) (v->getDimZ()-1); 

				int imaxx =  (int) (v->getDimX()-1); 
				int imaxy =  (int) (v->getDimY()-1); 
				int imaxz =  (int) (v->getDimZ()-1); 

				int xdim = xed-xst+1; 
				int ydim = yed-yst+1; 
				int zdim = zed-zst+1;

				bool show[300]; // table for indexing as to whether to show a segment or not

				int num = segMap->getNumSegments(); 
				for(int r=0;r<num;r++){
					SegmentMapEntry* ent = segMap->getEntry(r); 
					show[ent->getSegmentCode()] = ent->getState()!=vpt::HIDE_SEGMENT;
				}

                #pragma omp parallel for
				for(int k=zst;k<=zed;k++){
					for(int j=yst;j<=yed;j++){
						for(int i=xst;i<=xed;i++){
							float samp = tf->sampleCoarse(vval_scale((*v)(i,j,k))); 
							{
								float px = i/maxx-.5f;
								float py = j/maxy-.5f;
								float pz = k/maxz-.5f;
								osg::Vec3 L1(light1.x()-px,light1.y()-py,light1.z()-pz); 
								L1.normalize(); 
								osg::Vec3 L2(light2.x()-px,light2.y()-py,light2.z()-pz); 
								L2.normalize(); 
								osg::Vec3 L3(light3.x()-px,light3.y()-py,light3.z()-pz); 
								L3.normalize(); 

								float vx,vy,vz;
								bool skipBySegment = false; 
/*

								if(updateBySegment){
									skipBySegment = 
										segcode!=(*m)(i,j,k) && 
										(i<imaxx && segcode!=(*m)(i+1,j,k)) &&
										(i>0 && segcode!=(*m)(i-1,j,k)) &&
										(j<imaxy && segcode!=(*m)(i,j+1,k)) &&
										(j>0 && segcode!=(*m)(i,j-1,k)) &&
										(k<imaxz && segcode!=(*m)(i,j,k+1)) &&
										(k>0 && segcode!=(*m)(i,j,k-1)); 
								}

								skipBySegment = false; 
*/

								if(!skipBySegment && (updateOrig || (*border)(i,j,k))){
									if((*border)(i,j,k) ){
										vx = getValWrtDFunction(dfunc,m,tf,show,i+1,j,k,imaxx,imaxy,imaxz,fbuffer) - 
											getValWrtDFunction(dfunc,m,tf,show,i,j,k,imaxx,imaxy,imaxz,fbuffer); 

										vy = getValWrtDFunction(dfunc,m,tf,show,i,j+1,k,imaxx,imaxy,imaxz,fbuffer) - 
											getValWrtDFunction(dfunc,m,tf,show,i,j,k,imaxx,imaxy,imaxz,fbuffer); 

										vz = getValWrtDFunction(dfunc,m,tf,show,i,j,k+1,imaxx,imaxy,imaxz,fbuffer) - 
											getValWrtDFunction(dfunc,m,tf,show,i,j,k,imaxx,imaxy,imaxz,fbuffer); 										
/**/

									}
									else{
										vx = getVal(v,m,segMap,tf,i+1,j,k,imaxx,imaxy,imaxz) - 
											getVal(v,m,segMap,tf,i-1,j,k,imaxx,imaxy,imaxz); 
											//getVal(v,m,segMap,tf,i,j,k,imaxx,imaxy,imaxz); 

										vy = getVal(v,m,segMap,tf,i,j+1,k,imaxx,imaxy,imaxz) - 
											getVal(v,m,segMap,tf,i,j-1,k,imaxx,imaxy,imaxz); 
											//getVal(v,m,segMap,tf,i,j,k,imaxx,imaxy,imaxz); 

										vz = getVal(v,m,segMap,tf,i,j,k+1,imaxx,imaxy,imaxz) - 
											getVal(v,m,segMap,tf,i,j,k-1,imaxx,imaxy,imaxz); 
											//getVal(v,m,segMap,tf,i,j,k,imaxx,imaxy,imaxz); 
									}							 

									float vnorm = sqrt(vx*vx + vy*vy + vz*vz); 
									if(vnorm>=.0001){
										vx/=vnorm; 
										vy/=vnorm; 
										vz/=vnorm; 

										double res = abs(L1.x()*vx + L1.y()*vy + L1.z()*vz)*.3; 
										res += abs(L2.x()*vx + L2.y()*vy + L2.z()*vz)*.3; 
										res += abs(L3.x()*vx + L3.y()*vy + L3.z()*vz)*.3; 
										res = max(0.,min(1.,res));
										(*norm)(i,j,k) = normal_type_scale(res); 
									}
									else
										(*norm)(i,j,k) = 0; 

									if(updateOrig)
										(*orig)(i,j,k) = (*norm)(i,j,k); 
								}
								else
									(*norm)(i,j,k) = (*orig)(i,j,k); 
							}
						}
					}
				}	

//				blur(norm,border,xst,yst,zst,xed,yed,zed);

                #pragma omp parallel for
				for(int k=zst;k<=zed;k++){
					for(int j=yst;j<=yed;j++){
						for(int i=xst;i<=xed;i++){
							buffer[(k-zst)*xdim*ydim*2 + (j-yst)*xdim*2 + (i-xst)*2 +1] = (*norm)(i,j,k); 
							buffer[(k-zst)*xdim*ydim*2 + (j-yst)*xdim*2 + (i-xst)*2] = (*tv)(i,j,k);
						}
					}
				}

				updateTexture(nm,xst,yst,zst,xdim,ydim,zdim,buffer); 
		}

		static inline void updateAllNormals(NormalMap* nm, bool updateOrig=false){
			updateNormals(nm,osg::Vec3(0,0,-100), 0, 0, 0, 
				nm->getNormal()->getDimX()-1, 
				nm->getNormal()->getDimY()-1, 
				nm->getNormal()->getDimZ()-1,0,updateOrig); 
		}
	};
}

#endif