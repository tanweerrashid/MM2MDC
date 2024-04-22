#include "Volume.h"
#include "SelectionMap.h"
#include <cmath>
#include <iostream>

#include <osg/Vec4>
#include <osg/Plane>

using namespace std; 
using namespace vpt;

void VolumeOctree::build(){
	NodeRange r; 
	r.lx=r.ly=r.lz = 0; 
	r.ux = _volume->getDimX()-1; 
	r.uy = _volume->getDimY()-1; 
	r.uz = _volume->getDimZ()-1; 

	//_root = buildR(r,_depth);
	_root = NULL; 
}

// returns a node that contains the min/max of a subvolume from l[x,y,z] to u[x,y,z] (upper-inclusive)
VNode* VolumeOctree::buildR(const NodeRange& r, int depth){
	VNode* ret = new VNode(); 
	vval vmin = VV_MAX; 
	vval vmax = VV_MIN; 

	if(depth==0){
/*
		for(uint k=r.lz;k<=r.uz;k++){
			for(uint j=r.ly;j<=r.uy;j++){
				for(uint i=r.lx;i<=r.ux;i++){
					vmin = min(vmin,(*_volume)(i,j,k));
					vmax = max(vmax,(*_volume)(i,j,k));
				}
			}
		}
*/
	//	ret->setMin(vmin); 
//		ret->setMax(vmax);

		return ret; 
	}

	for(int j=0;j<8;j++){
		NodeRange nr = NodeRange::newRange(r, j); 

		VNode* child = buildR(nr, depth-1);
		//vmin = min(vmin, child->vmin());
		//vmax = max(vmax, child->vmax());
		ret->setChild(j, child);
	}

	return ret; 
}

std::pair<bool,float> VolumeOctree::intersectR(const NodeRange& range, VNode* cnode, const IntersectData& idata, int depth){
	std::pair<bool,float> failRet = make_pair(false,(float)1e9);

	//if(cnode==NULL){
		//cout<<"fail 1"<<endl;
	//return failRet; 
//}

	// intersect skip
	if(!isBoxIntersecting(range,idata,depth)) {
		//cout<<"fail 2"<<endl;
		return failRet; 
	}

	//cout<<depth<<endl;
	if(depth==0){
		// early skip 
		// Test each entry whether it is within the iso-range, if so, find the 
		// closest point that satisfies that condition.
		float bestd = 1e9; 
		osg::Vec3 best; 
		for(uint k=range.lz;k<=range.uz;k++){
			float zp = k/idata.dimZm1*idata.sclz-idata.offsetz; 
			for(uint j=range.ly;j<=range.uy;j++){
				float yp = j/idata.dimYm1*idata.scly-idata.offsety;
				for(uint i=range.lx;i<=range.ux;i++){
					//						if((*_volume)(i,j,k) > idata.vmax || (*_volume)(i,j,k)<idata.vmin)
					//if(idata.transFunc->sample((*_volume)(i,j,k)/1.f/VV_MAX)<.05 ||
					if(idata.transFunc->sampleCoarse((*_volume)(i,j,k))<.05 ||
						(*(idata.selectMap->getSelection()))(i,j,k)==0 ||
						idata.segMap->getEntryFast((int)((*idata.mask)(i,j,k)))->getState()==vpt::HIDE_SEGMENT || 
						idata.segMap->getEntryFast((int)((*idata.mask)(i,j,k)))->getSegmentCode() == idata.curId)
						continue; 

					osg::Vec3 C(i/idata.dimXm1*idata.sclx-idata.offsetx, yp, zp);
					if(idata.dimXm1==0) C[0]=0;
					if(idata.dimYm1==0) C[1]=0;
					if(idata.dimZm1==0) C[2]=0;

					osg::Vec3 PC = C - idata.pt; 
					float t = PC*idata.dir; 
					osg::Vec3 w = PC - (idata.dir * t);

					if(w.length2() < idata.r2){
						osg::Vec3 Q = C + w; 
						if(t < bestd){
							bestd = t; 
							best = Q; 
						}
					}
				}
			}
		}

		//cout<<"pt: "<<best.x()<<" "<<best.y()<<" "<<best.z()<<endl;

		// projection of the best point onto the intersection ray
		osg::Vec3 pt2isect = best - idata.pt; 
		float ret = pt2isect*idata.dir / idata.dirlen;

		//if(bestd>1e9-1)
		//	cout<<"fail4"<<endl;

		return make_pair(bestd<1e9-1, ret);
		//	return make_pair(bestd<1e9-1, best);
	}


	bool resb = false; 
	float resd = 1e9; 

	//if(cnode->hasChildren()){
	for(int j=0;j<8;j++){
		NodeRange nr = NodeRange::newRange(range, j); 

		pair<bool,float> cres = intersectR(nr, NULL, idata, depth-1);

		resb = cres.first || resb; 
		if(cres.first)
			resd = min(resd, cres.second);
	}
//}

	return make_pair(resb, resd);
}

