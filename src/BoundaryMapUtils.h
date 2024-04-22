#ifndef VPT_BOUNDARY_MAP_UTILS_H
#define VPT_BOUNDARY_MAP_UTILS_H

#include "BoundaryMap.h"
#include "VolumeUtils.h"
#include <vector> 
#include <algorithm>
#include <minc2.h>
#include <minc.h>

namespace vpt{
	class BoundaryMapUtils{
	public: 

		static BoundaryMap* readFromFile(const std::string& fstr){

		}

		static void writeToFile(const std::string& str){

		}

		static osg::Image* convertCornersToOSGImage(BoundaryMap* bm){
			osg::Image* ret = new osg::Image(); 
			VolumeData<unsigned char>* c = bm->getCorners(); 
			ret->setImage(c->getDimX()/4, c->getDimY(), c->getDimZ(), GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, (unsigned char*) c->data(),osg::Image::NO_DELETE);
			return ret; 
		}

		template<typename T>
		static inline void updateCornersTexture(BoundaryMap* bm, uint xs, uint ys, uint zs, uint xdim, uint ydim, uint zdim, T* ndata){
			if(TextureManager::inst()->hasBoundaryMapCornersTex3D(bm)){
				osg::ref_ptr<osg::Texture3D> tex = TextureManager::inst()->getOrCreateBoundaryMapCornersTex3D(bm); 
				if(tex->getSubloadCallback()==NULL){
					GLenum iformat = tex->getImage()->getInternalTextureFormat(); 
					GLenum pformat = tex->getImage()->getPixelFormat(); 
					GLenum type = tex->getImage()->getDataType(); 
					tex->setSubloadCallback( new VolumeSubloadCallback<T>(iformat, pformat, type)); 
				}

				osg::ref_ptr<VolumeSubloadCallback<T> > cb = dynamic_cast<VolumeSubloadCallback<T>*>(tex->getSubloadCallback());
				cb->initSubload(xs,ys,zs,xdim,ydim,zdim,ndata); 
			}
		}


		static inline void updateCornersTexture(BoundaryMap* bm){
			VolumeData<vval>* corners = bm->getCorners(); 
			updateCornersTexture(bm,0,0,0,corners->getDimX()/4,corners->getDimY(),corners->getDimZ(),corners->data()); 
		}

		static inline void updateCorners(BoundaryMap* bm, int xst, int yst, int zst, int xed, int yed, int zed, bool updateTexture=false){

			Volume* m = bm->getProject()->getMask();
			int xdim = xed-xst+1; 
			int ydim = yed-yst+1; 
			int zdim = zed-zst+1;

			for(int k=zst;k<=zed;k++){
				for(int j=yst;j<=yed;j++){
					for(int i=xst;i<=xed;i++){
						int ind = (k-zst)*xdim*ydim + (j-yst)*xdim + (i-xst);
						bm->writeCornerMask(i,j,k,(*m)(i,j,k)); 
					}
				}
			}

			if(updateTexture)
				updateCornersTexture(bm);
/*

			Volume* m = bm->getProject()->getMask();
			int xdim = xed-xst+1; 
			int ydim = yed-yst+1; 
			int zdim = zed-zst+1;

			int nxst = max(xst-1,0); 
			int nyst = max(yst-1,0);
			int nzst = max(zst-1,0);
			int nxed = min((int)m->getDimX()-2,xed);
			int nyed = min((int)m->getDimY()-2,yed);
			int nzed = min((int)m->getDimZ()-2,zed);

			int nxdim = nxed-nxst+1; 
			int nydim = nyed-nyst+1; 
			int nzdim = nzed-nzst+1; 

			for(int k=0;k<zdim;k++){
				for(int j=0;j<ydim;j++){
					for(int i=0;i<xdim;i++){
						bm->writeCornerMask(i,j,k,nxdim,nydim,nzdim,(*m)(i+xst,j+yst,k+zst)); 
					}
				}
			}

			updateCornersTexture(bm,nxst*2,nyst,nzst,nxdim*2,nydim,nzdim,bm->getCorners()->data());
*/
		}

		static inline void updateCorners(BoundaryMap* bm, int xst, int yst, int zst, int xed, int yed, int zed, vval* mask, bool updateTexture=false){

			int xdim = xed-xst+1; 
			int ydim = yed-yst+1; 
			int zdim = zed-zst+1;

			for(int k=zst;k<=zed;k++){
				for(int j=yst;j<=yed;j++){
					for(int i=xst;i<=xed;i++){
						int ind = (k-zst)*xdim*ydim + (j-yst)*xdim + (i-xst);
						bm->writeCornerMask(i,j,k,mask[ind]);
					}
				}
			}

			if(updateTexture)
				updateCornersTexture(bm);
/*
			Volume* m = bm->getProject()->getMask();
			int xdim = xed-xst+1; 
			int ydim = yed-yst+1; 
			int zdim = zed-zst+1;

			int nxst = max(xst-1,0); 
			int nyst = max(yst-1,0);
			int nzst = max(zst-1,0);
			int nxed = min((int)m->getDimX()-2,xed);
			int nyed = min((int)m->getDimY()-2,yed);
			int nzed = min((int)m->getDimZ()-2,zed);

			int nxdim = nxed-nxst+1; 
			int nydim = nyed-nyst+1; 
			int nzdim = nzed-nzst+1; 

			for(int k=zst;k<=zed;k++){
				for(int j=yst;j<=yed;j++){
					for(int i=xst;i<=xed;i++){
						int ind = (k-zst)*xdim*ydim + (j-yst)*xdim + i-xst;
						cout<<((int)mask[ind]); 
						bm->writeCornerMask(i-xst,j-yst,k-zst,nxdim,nydim,nzdim,mask[ind]);
					}
					cout<<endl;
				}
				cout<<endl;
			}

			updateCornersTexture(bm,nxst*2,nyst,nzst,nxdim*2,nydim,nzdim,bm->getCorners()->data());
*/
		}

		static inline void updateCorners(BoundaryMap* bm, bool updateTexture=true){
			Volume* m = bm->getProject()->getMask();
			updateCorners(bm,0,0,0,m->getDimX()-1,m->getDimY()-1,m->getDimZ()-1,updateTexture); 
		}

		template<typename T>
		static inline void updateVolumeTexture(BoundaryMap* bm, uint xs, uint ys, uint zs, uint xdim, uint ydim, uint zdim, T* ndata){
			if(TextureManager::inst()->hasBoundaryMapVolumeTex3D(bm)){
				osg::ref_ptr<osg::Texture3D> tex = TextureManager::inst()->getOrCreateBoundaryMapVolumeTex3D(bm); 
				if(tex->getSubloadCallback()==NULL){
					GLenum iformat = tex->getImage()->getInternalTextureFormat(); 
					GLenum pformat = tex->getImage()->getPixelFormat(); 
					GLenum type = tex->getImage()->getDataType(); 
					tex->setSubloadCallback( new VolumeSubloadCallback<T>(iformat, pformat, type)); 
				}

				osg::ref_ptr<VolumeSubloadCallback<T> > cb = dynamic_cast<VolumeSubloadCallback<T>*>(tex->getSubloadCallback());
				cb->initSubload(xs,ys,zs,xdim,ydim,zdim,ndata); 
			}
		}

		static inline void updateVolumeTexture(BoundaryMap* bm){
			VolumeData<bmvval>* volume = bm->getPackedVolume();
			updateVolumeTexture(bm,0,0,0,volume->getDimX()/4,volume->getDimY(),volume->getDimZ(),volume->data()); 
		}


    static inline void updateVolume(BoundaryMap* bm, int xst, int yst, int zst, int xed, int yed, int zed){
        VolumeData<float>* function = bm->getFunction();

        int xdim = xed-xst+1; 
        int ydim = yed-yst+1; 
        int zdim = zed-zst+1;

        for(int k=zst;k<=zed;k++){
            for(int j=yst;j<=yed;j++){
                for(int i=xst;i<=xed;i++){
                    bm->writeVolumeVal(i,j,k,bmv_scale((*function)(i,j,k)));
                }
            }
        }

        updateVolumeTexture(bm);
    }

    static inline void updateVolume(BoundaryMap* bm){
        VolumeData<float>* function = bm->getFunction();
        updateVolume(bm,0,0,0,function->getDimX()-1,function->getDimY()-1,function->getDimZ()-1); 
    }

		static osg::Image* convertVolumeToOSGImage(BoundaryMap* bm){
			osg::Image* ret = new osg::Image(); 
/*
			VolumeData<unsigned char>* c = bm->getPackedVolume(); 
			ret->setImage(c->getDimX()/4, c->getDimY(), c->getDimZ(), GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, (unsigned char*) c->data(),osg::Image::NO_DELETE);
*/

			VolumeData<bmvval>* c = bm->getPackedVolume(); 
			ret->setImage(c->getDimX()/4, c->getDimY(), c->getDimZ(), BMV_INTERNAL_TEXTURE_TYPE, GL_RGBA, BMV_TEXTURE_TYPE, (unsigned char*) c->data(),osg::Image::NO_DELETE);
			return ret; 
		}

		static inline void updateBoundary(BoundaryMap* bm, int xst, int yst, int zst, int xed, int yed, int zed){
			Volume* m = bm->getProject()->getMask();
			Volume* border1 = bm->getAuxBoundary(); 
			Volume* border2 = bm->getBoundary(); 

			//Volume* border1 = bm->getBoundary();
			Volume* v = bm->getProject()->getVolume(); 

			float maxx = (int) (v->getDimX()-1); 
			float maxy = (int) (v->getDimY()-1); 
			float maxz = (int) (v->getDimZ()-1); 

			int imaxx =  (int) (v->getDimX()-1); 
			int imaxy =  (int) (v->getDimY()-1); 
			int imaxz =  (int) (v->getDimZ()-1); 

			int xdim = xed-xst+1; 
			int ydim = yed-yst+1; 
			int zdim = zed-zst+1;

			int di[] = {1,-1, 1,-1, 1,-1, 1,-1, 1, 0, 1, 0, 1, 0};
			int dj[] = {1, 1,-1,-1, 1, 1,-1,-1, 0, 1, 1, 0, 0, 1};
			int dk[] = {1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1}; 

            #pragma omp parallel for
			for(int k=zst;k<=zed;k++){
				for(int j=yst;j<=yed;j++){
					for(int i=xst;i<=xed;i++){
						vval v = (*m)(i,j,k);
						for(int p=0;p<14;p++){
							int ni = i+di[p], nj = j+dj[p], nk = k+dk[p];
							if(border1->valid(ni,nj,nk) && v!=(*m)(ni,nj,nk)){
								(*border1)(i,j,k) = 255;
								break;
							}
						}
					}
				}
			}

            #pragma omp parallel for
			for(int k=zst;k<=zed;k++){
				for(int j=yst;j<=yed;j++){
					for(int i=xst;i<=xed;i++){
						(*border2)(i,j,k) = ((*border1)(i,j,k) ||
							(i<imaxx && (*border1)(i+1,j,k)) ||
							(i>0 && (*border1)(i-1,j,k)) ||
							(j<imaxy && (*border1)(i,j+1,k)) ||
							(j>0 && (*border1)(i,j-1,k)) ||
							(k<imaxz && (*border1)(i,j,k+1)) ||
							(k>0 && (*border1)(i,j,k-1)))? 255 : 0; 
					}
				}
			}
		}

		static inline void updateBoundaryWithVolume2(BoundaryMap* bm, FVolume* fvol, Volume* mask, float iso, vval seg){
			VolumeData<float>* function = bm->getFunction();
			VolumeData<float>* afunction = bm->getAuxFunction();

			int dimx = (int) fvol->getDimX(); 
			int dimy = (int) fvol->getDimY(); 
			int dimz = (int) fvol->getDimZ(); 

			int dx[] = {-1, 0, 0, 1, 0, 0}; 
			int dy[] = { 0,-1, 0, 0, 1, 0}; 
			int dz[] = { 0, 0,-1, 0, 0, 1}; 

			float mmin = 1e9;
			float mmax = -1e9; 

			for(int k=0;k<dimz;k++){
				for(int j=0;j<dimy;j++){
					for(int i=0;i<dimx;i++){
						float val = (*fvol)(i,j,k);
						if(val >= iso){
							float minval = val;
							float xval = 0;
							int count=0;
							for(int n=0;n<6;n++){
								int ni = dx[n] + i;
								int nj = dy[n] + j; 
								int nk = dz[n] + k; 

								if(ni>=0 && ni<dimx && nj>=0 && nj<dimy && nk>=0 && nk<dimz){
									if((*fvol)(ni,nj,nk)<iso){
										minval = min(minval,(*fvol)(ni,nj,nk)); 
										xval += (*fvol)(ni,nj,nk);
										count++; 
									}
								}
							}
							// all six edges do not intersect the iso-surface
							xval = minval; 
							if(minval>iso)
								(*function)(i,j,k) = 1; 
							else
								(*function)(i,j,k) = 1.f-(iso-xval)/(val-minval); 
						}
					}
				}
			}

			for(int k=0;k<dimz;k++){
				for(int j=0;j<dimy;j++){
					for(int i=0;i<dimx;i++){
						float val = (*fvol)(i,j,k);
						if(val >= iso)
							(*mask)(i,j,k) = seg;
					}
				}
			}

			for(int k=0;k<dimz;k++)
				for(int j=0;j<dimy;j++)
					for(int i=0;i<dimx;i++)
						(*afunction)(i,j,k) = (*function)(i,j,k);

			VolumeUtils::updateTexture(mask, 0, 0, 0, dimx, dimy, dimz, mask->data());	

			updateCorners(bm); 
			updateVolume(bm); 
		}

		static inline void updateBoundaryWithVolume(BoundaryMap* bm, FVolume* fvol, Volume* mask, float iso, vval seg){
			VolumeData<float>* function = bm->getFunction();
			VolumeData<float>* afunction = bm->getAuxFunction();
			VolumeData<bool>* use = new VolumeData<bool>(); 
			std::vector<float> sarray; // for sorting

			int dimx = (int) fvol->getDimX(); 
			int dimy = (int) fvol->getDimY(); 
			int dimz = (int) fvol->getDimZ(); 

			use->resize(dimx,dimy,dimz); 

			//int dx[] = {-1, 0, 0, 1, 0, 0}; 
			//int dy[] = { 0,-1, 0, 0, 1, 0}; 
			//int dz[] = { 0, 0,-1, 0, 0, 1}; 

			int dx[] = {1,-1, 1,-1, 1,-1, 1,-1, 1, 0, 1, 0, 1, 0,-1, 0,-1, 0,-1, 0, 1, 1,-1,-1, 0, 0};
			int dy[] = {1, 1,-1,-1, 1, 1,-1,-1, 0, 1, 1, 0, 0, 1, 0,-1,-1, 0, 0,-1,-1, 0, 1, 0,-1, 1};
			int dz[] = {1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 0, 0, 0,-1,-1,-1, 0,-1, 0, 1, 1,-1}; 

			for(int k=0;k<dimz;k++){
				for(int j=0;j<dimy;j++){
					for(int i=0;i<dimx;i++){
						float val = (*fvol)(i,j,k);
						bool good = false; 
						for(int n=0;n<26;n++){
							int ni = dx[n] + i;
							int nj = dy[n] + j; 
							int nk = dz[n] + k; 

							if(ni>=0 && ni<dimx && nj>=0 && nj<dimy && nk>=0 && nk<dimz){
								if(((*fvol)(ni,nj,nk)<iso && val>=iso) || 
									((*fvol)(ni,nj,nk)>=iso && val<iso)){
										good = true; 
										break;
								}
							}
						}
						(*use)(i,j,k) = good; 

						if(good){
							sarray.push_back((*fvol)(i,j,k)-iso); 
							//(*function)(i,j,k) = abs(val-iso);
						}
					}
				}
			}
			std::sort(sarray.begin(),sarray.end()); 
			float frac=.80f;
			int offset = (sarray.size() - (frac* sarray.size()))/2; 
			float mmin = sarray[offset]; 
			float mmax = sarray[sarray.size()-1-offset];
			float mrange = mmax-mmin;

			for(int k=0;k<dimz;k++){
				for(int j=0;j<dimy;j++){
					for(int i=0;i<dimx;i++){
						float val = (*fvol)(i,j,k)-iso; 
						if((*use)(i,j,k)){
							if(val>=mmin && val<=mmax)
								(*function)(i,j,k) = abs((val-mmin)/(mrange)*2-1.f); 
							else
								(*function)(i,j,k) = 1; 
						}
						else{
						//	(*function)(i,j,k) = 1; 
						}
					}
				}
			}



			for(int k=0;k<dimz;k++){
				for(int j=0;j<dimy;j++){
					for(int i=0;i<dimx;i++){
						float val = (*fvol)(i,j,k);
						if(val >= iso)
							(*mask)(i,j,k) = seg;
					}
				}
			}

			afunction->copy(function); 
			VolumeUtils::updateTexture(mask, 0, 0, 0, dimx, dimy, dimz, mask->data());	

			updateCorners(bm); 
			updateVolume(bm); 

			delete use; 
		}


		static inline void updateBoundaryWithVolume3(BoundaryMap* bm, FVolume* fvol, Volume* mask, float iso, vval seg){
			VolumeData<float>* function = bm->getFunction();
			VolumeData<float>* afunction = bm->getAuxFunction();
			VolumeData<float>* nlen = new VolumeData<float>(); 
			VolumeData<vval>* count = new VolumeData<vval>(); 

			int dimx = (int) fvol->getDimX(); 
			int dimy = (int) fvol->getDimY(); 
			int dimz = (int) fvol->getDimZ(); 
			nlen->resize(dimx,dimy,dimz);
			count->resize(dimx,dimy,dimz);

			int dx[] = {-1, 0, 0, 1, 0, 0}; 
			int dy[] = { 0,-1, 0, 0, 1, 0}; 
			int dz[] = { 0, 0,-1, 0, 0, 1}; 

			float mmin = 1e9;
			float mmax = -1e9; 

			float p[8]; 
/*
			for(int k=0;k<dimz;k++){
				for(int j=0;j<dimy;j++){
					for(int i=0;i<dimx;i++){
						cout<<(*fvol)(i,j,k)<<endl;
					}
					cout<<endl;
				}
				cout<<"**************"<<endl;
			}
*/
			function->setAll(0);
			nlen->setAll(0);
			count->setAll(0);

			for(int k=0;k<dimz-1;k++){
				for(int j=0;j<dimy-1;j++){
					for(int i=0;i<dimx-1;i++){
						float val = (*fvol)(i,j,k);
						bool good = false; 

						p[0] = (*fvol)(i,j,k); 
						p[1] = (*fvol)(i+1,j,k); 
						p[2] = (*fvol)(i,j+1,k); 
						p[3] = (*fvol)(i+1,j+1,k); 
						p[4] = (*fvol)(i,j,k+1); 
						p[5] = (*fvol)(i+1,j,k+1); 
						p[6] = (*fvol)(i,j+1,k+1); 
						p[7] = (*fvol)(i+1,j+1,k+1); 

						for(int jj=0;jj<8;jj++){
							if(p[jj]>=iso){
								for(int kk=0;kk<8;kk++){
									if(p[kk]<iso){
										good = true;
										break;
									}
								}
								break;
							}
						}

						if(!good)
							continue;

						float dx = .25 * (p[1]-p[0]+p[3]-p[2]+p[5]-p[4]+p[7]-p[6])*(dimx-1); 
						float dy = .25 * (p[2]-p[0]+p[3]-p[1]+p[6]-p[4]+p[7]-p[5])*(dimy-1); 
						float dz = .25 * (p[4]-p[0]+p[5]-p[1]+p[6]-p[2]+p[7]-p[3])*(dimz-1);

						float len = sqrt(dx*dx + dy*dy + dz*dz); 

						for(int jj=0;jj<8;jj++){
							int nx=i,ny=j,nz=k;
							if(jj&1) nx++; 
							if(jj&2) ny++; 
							if(jj&4) nz++; 

							(*nlen)(nx,ny,nz) += len;
							(*count)(nx,ny,nz)++; 

							(*function)(nx,ny,nz) += ((*fvol)(nx,ny,nz)/len);
						}
					}
				}
			}

			for(int k=0;k<dimz;k++){
				for(int j=0;j<dimy;j++){
					for(int i=0;i<dimx;i++){
						if((*count)(i,j,k)>0){
							(*function)(i,j,k)/=(*count)(i,j,k); 
							(*nlen)(i,j,k)/=(*count)(i,j,k);
							//cout<<"nlen: "<<(*nlen)(i,j,k)<<endl;
						}
						//cout<<"nf: "<<(*function)(i,j,k)<<" of:"<<(*fvol)(i,j,k)<<endl;
						if((*count)(i,j,k)>0)
							(*function)(i,j,k) = abs(min(1.f,(*function)(i,j,k)-(iso/(*nlen)(i,j,k))));
						else
							(*function)(i,j,k) = 1;

						(*afunction)(i,j,k) = (*function)(i,j,k);
//						cout<<"f: "<<(*function)(i,j,k)<<endl;
					}
				}
			}


			for(int k=0;k<dimz;k++){
				for(int j=0;j<dimy;j++){
					for(int i=0;i<dimx;i++){
						float val = (*fvol)(i,j,k);
						if(val >= iso)
							(*mask)(i,j,k) = seg;
					}
				}
			}
	

			VolumeUtils::updateTexture(mask, 0, 0, 0, dimx, dimy, dimz, mask->data());	

			updateCorners(bm); 
			updateVolume(bm); 

			delete nlen; 
			delete count;
		}



    template<typename T>
    static inline void blur(VolumeData<T>* v, int xst, int yst, int zst, int xed, int yed, int zed){
        VolumeData<T>* tmp = new VolumeData<T>(); 
        tmp->resize(v->getDimX(),v->getDimY(),v->getDimZ());
        
        /*float normal_mask[27] = {
            1,3,1,
            3,3,3,
            1,3,1,

            1,3,1,
            3,5,3,
            1,3,1,

            1,3,1,
            3,3,3,
            1,3,1
        };*/
        

        float normal_mask[27]= {
            0.030838, 0.0571701, 0.030838, 
            0.0571701, 0.105987, 0.0571701, 
            0.030838, 0.0571701, 0.030838, 

            0.0571701, 0.105987, 0.0571701, 
            0.105987, 0.196488, 0.105987, 
            0.0571701, 0.105987, 0.0571701, 

            0.030838, 0.0571701, 0.030838, 
            0.0571701, 0.105987, 0.0571701, 
            0.030838, 0.0571701, 0.030838
        };
        
        float total = 1.76515;
        int dd = 1; 
        int mdim = (2*dd+1); 
        int mdim2 = mdim*mdim;
        int maxx = v->getDimX()-1;
        int maxy = v->getDimY()-1;
        int maxz = v->getDimZ()-1;

        #pragma omp parallel for
        for(int z=zst;z<=zed;z++){
            for(int y=yst;y<=yed;y++){
                for(int x=xst;x<=xed;x++){
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
                                float mval = normal_mask[(k+dd)*mdim2 + (j+dd)*mdim + (i+dd)]; 
                                acc+=mval*(*v)(nx,ny,nz); 
                                t2+=mval; 
                            }
                        }					
                    }
                    if(t2>0)
                        (*tmp)(x,y,z) = (T) (acc/t2); 
                    else
                        (*tmp)(x,y,z) = (*v)(x,y,z); 
                }
            }
        }

        #pragma omp parallel for
        for(int z=zst;z<=zed;z++){
            for(int y=yst;y<=yed;y++)
                for(int x=xst;x<=xed;x++)
                    (*v)(x,y,z) = (*tmp)(x,y,z); 
        }

        delete tmp; 
    }



// Update Boundary with mask
// Given a discrete set of masks, how can we update the boundary? 
static inline void updateBoundaryWithMask(BoundaryMap* bm, SegmentMap* sm, Volume* mask){
    // First create binary volumes for 1
    FVolume* fvol = bm->getFunction();         
    
    vector<FVolume*> nbounds; 

    int dimx = mask->getDimX(); 
    int dimy = mask->getDimY(); 
    int dimz = mask->getDimZ(); 

    nbounds.assign(sm->getNumSegments(),NULL); 
    //cout << "nbounds size: " << nbounds.size() << endl;
    
    for(int s=0;s<sm->getNumSegments();s++){
        SegmentMapEntry* ent = sm->getEntry(s); 
        vval code = ent->getSegmentCode(); 
        //cout << "Code: " << (int)code << endl;
        if(code==SegmentMap::SELECTION_SEGMENT)
        continue; 

        nbounds[s] = new FVolume(); 
        nbounds[s]->resize(dimx,dimy,dimz); 

        for(int k=0;k<dimz;k++){
            for(int j=0;j<dimy;j++){
                for(int i=0;i<dimx;i++){
                    if((*mask)(i,j,k)==code)                        
                        (*nbounds[s])(i,j,k) = (*fvol)(i,j,k); 
                    else
                        (*nbounds[s])(i,j,k) = 0; 
                }
            }
        }
    }

    
    /*ofstream f1;
    f1.open("/home/trash001/Desktop/PoweiSmoothing_nbounds[2].txt");
    f1 << "nbounds" << "\t" << "fvol" << "\t" << "mask" << endl;    
    for (int i = 0; i < dimx; i++) {
        for (int j = 0; j < dimy; j++) {
            for (int k = 0; k < dimz; k++) {
                //if ((int)(*mask)(i,j,k) == 2) {
                    f1 << (*nbounds[2])(i, j, k) << "\t" << (*fvol)(i, j, k) << "\t" << (int)(*mask)(i, j, k) << endl;
                //}    
            }
        }
    }
    f1.close();*/
    
    
    for(uint j=0;j<nbounds.size();j++){
        if(nbounds[j]){
            blur(nbounds[j],0,0,0,dimx-1,dimy-1,dimz-1); 
            //blur(nbounds[j],0,0,0,dimx-1,dimy-1,dimz-1); 
        }
    }

    
    
    
    /*ofstream f2;
    f2.open("/home/trash001/Desktop/PoweiSmoothing_AfterBlurring.txt");
    f2 << "nbounds" << "\t" << "fvol" << "\t" << "mask" << endl;    
    for (int i = 0; i < dimx; i++) {
        for (int j = 0; j < dimy; j++) {
            for (int k = 0; k < dimz; k++) {
                //if ((int)(*mask)(i,j,k) == 2) {
                    f2 << (*nbounds[0])(i, j, k) << "\t" << (*fvol)(i, j, k) << "\t" << (int)(*mask)(i, j, k) << endl;
                //}    
            }
        }
    }
    f2.close();*/
    
    
    int nboundsSize = (int) nbounds.size(); 
    
    for(int i=0;i<dimx;i++){
        for(int j=0;j<dimy;j++){
            for(int k=0;k<dimz;k++){

                float best1 = -1000; 
                float best2 = -1000; 
                for(int s=0;s<nboundsSize;s++){
                    if(nbounds[s]){
                        if((*nbounds[s])(i,j,k) > best1){
                            best2 = best1; 
                            best1 = (*nbounds[s])(i,j,k); 
                        }
                        else if((*nbounds[s])(i,j,k) > best2){
                            best2 = (*nbounds[s])(i,j,k); 
                        }
                    }
                }
                (*fvol)(i,j,k) = best1-best2; 
                //(*fvol)(i,j,k) = best1; 
            }
        }
    }

    /*ofstream f3;
    f3.open("/home/trash001/Desktop/PoweiSmoothing_AfterBest1Best2.txt");
    f3 << "nbounds" << "\t" << "fvol" << "\t" << "mask" << endl;    
    for (int i = 0; i < dimx; i++) {
        for (int j = 0; j < dimy; j++) {
            for (int k = 0; k < dimz; k++) {
                //if ((int)(*mask)(i,j,k) == 2) {
                    f3 << (*nbounds[0])(i, j, k) << "\t" << (*fvol)(i, j, k) << "\t" << (int)(*mask)(i, j, k) << endl;
                //}    
            }
        }
    }
    f3.close();*/
    
    for(uint j=0;j<nbounds.size();j++){
        if(nbounds[j])
            delete nbounds[j]; 
    }
    updateVolume(bm); 
}




/**
static void readBlurredMincVolume(std::string filename, FVolume* vol) {
    // Read MINC volume. 
    mihandle_t minc_volume;
    double voxel;
    int result;
    misize_t location[3];

    // Open MINC file
    //std::string filename = "filename;
    //std::string filename = "data/mincblur/Volume.mnc";
    cout << "Opening mincblurred file " << filename << endl;
    result = miopen_volume(filename.c_str(), MI2_OPEN_READ, &minc_volume);	
    if (result == MI_ERROR) cout << "Error opening input file: " << result << endl;

    // Set order of MINC x, y and z spaces
    char *dimorder[3];
    dimorder[0] = "xspace";
    dimorder[1] = "yspace";
    dimorder[2] = "zspace";
    result = miset_apparent_dimension_order_by_name(minc_volume, 3, dimorder); 
    if (result == MI_ERROR) cout << "Error in setting dimension order" << endl;

    // Get Volume dimension and step sizes
    midimhandle_t dimensions[3];
    misize_t sizes[3];
    result = miget_volume_dimensions(minc_volume, MI_DIMCLASS_SPATIAL, MI_DIMATTR_ALL, MI_DIMORDER_APPARENT, 3, dimensions);
    if (result == MI_ERROR) cout << "Error in retrieving volume dimensions" << endl;	
    result = miget_dimension_sizes(dimensions, 3, sizes);
    if (result == MI_ERROR) cout << "Error in retrieving dimension sizes" << endl;

    double stepSize[3];
    result = miget_dimension_separations(dimensions, MI_ORDER_APPARENT, 3, stepSize);
    if (result == MI_ERROR) cout << "Error in retrieving step sizes" << endl;

    double start[3];
    result = miget_dimension_starts(dimensions, MI_ORDER_FILE, 3, start);
    if (result == MI_ERROR) cout << "Error in retrieving volume origin" << endl;
   
    // w is along x axis;  h is along y axis; d is along z axis. 
    int w = sizes[0], h = sizes[1], d = sizes[2];
    cout << "Volume dimensions: [" << w << ", " << h << ", " << d << "]" << endl;
    cout << "Volume stepsize: [" << stepSize[0] << ", " << stepSize[1] << ", " << stepSize[2] << "]" << endl;
    cout << "Volume start: [" <<  start[0] << ", " << start[1] << ", " << start[2] << "]" << endl;       
    
    cout << "Reading mincblurred voxels" << endl;
    #pragma omp parallel for
    for (int wi = 0; wi < w; wi = wi + 1) {
        for (int hi = 0; hi < h; hi = hi + 1) {
            for (int di = 0; di < d; di = di + 1) {                
                location[0] = wi;
                location[1] = hi;
                location[2] = di;                                

                // Get voxel value at specified voxel location. 
                result = miget_voxel_value(minc_volume, location, 3, &voxel);
                if (result == MI_ERROR) cout << "Error reading voxels." << endl;
                
                (*vol)(wi, hi, di) = voxel; 
            }
        }
    }

    result = miclose_volume(minc_volume);
    if (result == MI_ERROR) cout << "Error in closing MINC volume" << endl;

    cout << "Finished reading mincblurred voxels" << endl;
}


// Using mincblurred volumes for smoothing
/*static inline void updateBoundaryWithMincBlurredVolume(BoundaryMap* bm, SegmentMap* sm, Volume* mask, std::string blurredFilenames[]){
    // First create binary volumes for 1
    FVolume* fvol = bm->getFunction(); 

    vector<FVolume*> nbounds; 

    int dimx = mask->getDimX(); 
    int dimy = mask->getDimY(); 
    int dimz = mask->getDimZ(); 

    nbounds.assign(sm->getNumSegments(), NULL);     
    
    /*std::string filenames[] = {
      "data/mincblur labels_on_colin_Nov2010/Object26_material0_stddev0.2_blur.mnc",    // Material index 0, i.e. background
      "",                                                                               // Material index 1, i.e. selection
      "data/mincblur labels_on_colin_Nov2010/Object26_material1_stddev4.0_blur.mnc"     // Material index 2, material of interest
    };*/
    
    // Populate nbounds[] for remaining all segments except SELECTION_SEGMENT
    /*for (int s = 0; s < sm->getNumSegments(); s++) {
        if (s == SegmentMap::SELECTION_SEGMENT)
            continue;
        
        nbounds[s] = new FVolume(); 
        nbounds[s]->resize(dimx, dimy, dimz); 
        
        readBlurredMincVolume(blurredFilenames[s], nbounds[s]);
    }
    
    int nboundsSize = (int) nbounds.size(); 
    
    for(int i=0;i<dimx;i++){
        for(int j=0;j<dimy;j++){
            for(int k=0;k<dimz;k++){

                float best1 = -1000; 
                float best2 = -1000; 
                for(int s=0;s<nboundsSize;s++){
                    if(nbounds[s]){
                        if((*nbounds[s])(i,j,k) > best1){
                            best2 = best1; 
                            best1 = (*nbounds[s])(i,j,k); 
                        }
                        else if((*nbounds[s])(i,j,k) > best2){
                            best2 = (*nbounds[s])(i,j,k); 
                        }
                    }
                }
                (*fvol)(i,j,k) = best1-best2; 
                //(*fvol)(i,j,k) = best1; 
            }
        }
    }

    /*ofstream f1;
    f1.open("/home/trash001/Desktop/MincBlurringAfter.txt");
    f1 << "nbounds" << "\t" << "fvol" << "\t" << "mask" << endl;    
    for (int i = 0; i < dimx; i++) {
        for (int j = 0; j < dimy; j++) {
            for (int k = 0; k < dimz; k++) {
                //if ((int)(*mask)(i,j,k) == 2) {
                    f1 << (*nbounds[0])(i, j, k) << "\t" << (*fvol)(i, j, k) << "\t" << (int)(*mask)(i, j, k) << endl;
                //}    
            }
        }
    }*/
    
    /*for(uint j=0;j<nbounds.size();j++){
        if(nbounds[j])
            delete nbounds[j]; 
    }
    updateVolume(bm); 
}*/















	};
}


#endif
