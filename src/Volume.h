#ifndef VPT_VOLUME_H
#define VPT_VOLUME_H

#include "Data.h"
#include "Common.h"
#include "TransferFunction.h" 
#include "SegmentMap.h"
//#include "Rendering/SelectionMap.h"

#include "Savable.h"

#include <osg/Vec3>

#include <algorithm> 
#include <cmath>
#include "mrc.h"

using namespace std; 

namespace vpt{

	enum VolumeIOType {MRC_TYPE};
	enum VolumeDataType {VDT_BYTE,VDT_FLOAT};
	class SelectionMap;

	class MRCData{
	public: 
		mrcH actualHeader; 
		double stepx,stepy,stepz; // length of a cell in each dimension
		double scalex, scaley, scalez; // the max the scales should be 1
		double offsetx, offsety, offsetz; // offset[x,y,z] = dim[x,y,z]/dimmax/2
		double vmax, vmin; 
	}; 

	// a three dimension grid of data
	template <typename T> 
	class VolumeData : public Savable{
	protected: 
		T* _data; 
		unsigned int _dimX,_dimY,_dimZ; 
		int _type;
		// additional header data
		MRCData* _mrcData; 

	public: 
		VolumeData() 
			: _data(NULL), _type(sizeof(T)), 
			_mrcData(NULL),
			_dimX(0),_dimY(0),_dimZ(0)
		{}

        VolumeData(uint x, uint y, uint z)
			: _data(NULL), _type(sizeof(T)), 
			_mrcData(NULL){
                resize(x,y,z);
        }

		~VolumeData(){ 
			clear(); 
			if(_mrcData)
				delete _mrcData; 
		}

		// create a volume of this size
		inline void resize(uint x, uint y, uint z){
			clear();
			_data = new T[x*y*z];
			_dimX = x; 
			_dimY = y; 
			_dimZ = z;
			zero();
		}


		/****************/
		inline T getData(int i) {return _data[i];}
		/****************/


		inline void zero(){	setAll(0); }
		inline void setAll(T v){ for(uint j=0;j<_dimX*_dimY*_dimZ;j++) _data[j]=v; }
		inline void setAll(int xst, int yst, int zst, int xed, int yed, int zed, T v){
			for(int k=zst;k<=zed;k++)
				for(int j=yst;j<=yed;j++)
					for(int i=xst;i<=xed;i++)
						(*this)(i,j,k) = v; 
		}

		inline bool copy(VolumeData<T>* src){
			if(src->getDimX()!=getDimX() || src->getDimY()!=getDimY() || 
				src->getDimZ()!=getDimZ())
				return false; 
			int size = (int) getSize(); 

            #pragma omp parallel for
			for(int j=0;j<size;j++)
				_data[j] = src->_data[j]; 			
			return true; 
		}

		inline void clear(){
			if(_data!=NULL)
				delete [] _data; 
			_data = NULL;
		}

		inline MRCData* getMRCData() { return _mrcData; }
		inline void setMRCData(MRCData* m) { _mrcData = m; }

		inline unsigned int getDimX() const { return _dimX; }
		inline unsigned int getDimY() const { return _dimY; }
		inline unsigned int getDimZ() const { return _dimZ; }
		inline unsigned int getSize() const { return _dimX*_dimY*_dimZ; }

		inline bool ready() const { return _data!=NULL; };
		inline virtual T* data(){ return _data; }

		inline double sample(double x, double y, double z){
			double stepx = (getDimX()-1) * x; 
			double stepy = (getDimY()-1) * y; 
			double stepz = (getDimZ()-1) * z; 

			int minx = (int) stepx; 
			int miny = (int) stepy; 
			int minz = (int) stepz; 

			if(minx<0 || minx>=(int)getDimX()-1 ||
				miny<0 || miny>=(int)getDimY()-1 ||
				minz<0 || minz>=(int)getDimZ()-1)
				return 0; 

			int	maxx = minx+1; 
			int maxy = miny+1; 
			int maxz = minz+1; 

			double res=0; 

			res += (*this)(minx,miny,minz) * (1-x) * (1-y) * (1-z); 
			res += (*this)(maxx,miny,minz) * x * (1-y) * (1-z); 
			res += (*this)(minx,maxy,minz) * (1-x) * y * (1-z); 
			res += (*this)(maxx,maxy,minz) * x * y * (1-z); 
			res += (*this)(minx,miny,maxz) * (1-x) * (1-y) * z; 
			res += (*this)(maxx,miny,maxz) * x * (1-y) * z; 
			res += (*this)(minx,maxy,maxz) * (1-x) * y * z; 
			res += (*this)(maxx,maxy,maxz) * x * y * z; 

			return res; 
		}


		inline double sampleSmooth(double x, double y, double z, int ksize, double* mask){
			double stepx = 1./(getDimX()-1); 
			double stepy = 1./(getDimY()-1); 
			double stepz = 1./(getDimZ()-1); 

			int c = (ksize-1)/2; 

			double res = 0; 
			double norm = 0; 

			for(int k=0;k<ksize;k++){
				for(int j=0;j<ksize;j++){
					for(int i=0;i<ksize;i++){
						res+=sample(x-stepx*(i-c),y-stepy*(j-c),z-stepz*(k-c))*mask[k*ksize*ksize+j*ksize+i]; 
						norm+=mask[k*ksize*ksize+j*ksize+i]; 
					}
				}
			}

			return res/norm; 
		}

		inline T& operator()(const uint& j){return _data[j]; }
		inline T operator()(const uint& j) const{ return _data[j]; }

		// need to be careful not to exceed the bounds
		inline T& operator()(const uint& x, const uint& y, const uint& z){
			return _data[(z*_dimY + y)*_dimX + x];
		}

		inline T operator() (const uint& x, const uint& y, const uint& z) const{
			return _data[(z*_dimY + y)*_dimX + x];
		}

		inline bool valid(int x, int y, int z){
			return x>=0 && x<(int)_dimX && y>=0 && y<(int)_dimY && z>=0 && z<(int)_dimZ; 
		}
	};
	typedef VolumeData<vval> Volume; 
	typedef VolumeData<float> FVolume; 


	template <typename T>
	class VolumeOctreeNode{
	protected: 
		VolumeOctreeNode** _child; 
		T _min, _max; 
	public: 
		inline VolumeOctreeNode(){
			_child = NULL; 
			_min = VV_MAX; 
			_max = VV_MIN; 
		}

		~VolumeOctreeNode(){
			if(_child){
				for(int j=0;j<8;j++)
					delete _child[j]; 
			}
			delete [] _child; 
		}

		inline T vmin() const { return _min; }
		inline T vmax() const { return _max; }
		inline VolumeOctreeNode* child(int i) { 
			if(_child)
				return _child[i]; 
			return NULL;
		}

		inline void setMin(T val) { _min = val; }
		inline void setMax(T val) { _max = val; }

		inline void setChild(int i, VolumeOctreeNode* v){
			if(_child) 
				_child[i] = v; 
			else{
				_child = new VolumeOctreeNode*[8]; 
				for(int j=0;j<8;j++)
					_child[j] = NULL; 
				_child[i] = v; 
			}
		}
		inline bool hasChildren() { return _child!=NULL; }

	}; 
	typedef VolumeOctreeNode<vval> VNode; 

	typedef struct VolumeOctreeIntersectBundle{
		osg::Vec3 pt; 
		osg::Vec3 dir; 
		vval vmin; 
		vval vmax; 
		TransferFunction* transFunc; 
		SegmentMap* segMap; 
		SelectionMap* selectMap; 
		Volume* mask; 
		int curId; 
		float r2; 
		float dimXm1; 
		float dimYm1; 
		float dimZm1; 
		float dimmaxm1; 
		float offsetx; 
		float offsety; 
		float offsetz; 
		float sclx; 
		float scly; 
		float sclz; 

		float dirlen;
	} IntersectData; 

	typedef struct VolumeOctreeNodeRange{
		uint lx, ly, lz; 
		uint ux, uy, uz; 

		static inline VolumeOctreeNodeRange newRange(const VolumeOctreeNodeRange& r, int j){
			VolumeOctreeNodeRange nr; 
			nr.lx = r.lx; 
			nr.ly = r.ly; 
			nr.lz = r.lz; 
			nr.ux = r.lx + (r.ux-r.lx+1)/2; 
			nr.uy = r.ly + (r.uy-r.ly+1)/2; 
			nr.uz = r.lz + (r.uz-r.lz+1)/2;

			if(j&(1<<0)){
				nr.lx = nr.ux; 
				nr.ux = r.ux; 
			}
			if(j&(1<<1)){
				nr.ly = nr.uy; 
				nr.uy = r.uy; 
			}
			if(j&(1<<2)){
				nr.lz = nr.uz; 
				nr.uz = r.uz; 
			}
			return nr; 
		}
	} NodeRange; 
	typedef NodeRange VolumeRange; 

	class VolumeOctree {
	protected: 
		VNode* _root; 
		Volume* _volume; 
		int _depth; 

		void build(); 
		VNode* buildR(const NodeRange& nrange, int depth);

		inline bool isXYPlaneIntersecting(const NodeRange& range, const IntersectData& idata, int depth);
		inline bool isYZPlaneIntersecting(const NodeRange& range, const IntersectData& idata, int depth);
		inline bool isXZPlaneIntersecting(const NodeRange& range, const IntersectData& idata, int depth);

		inline bool isBoxIntersecting(const NodeRange& range, const IntersectData& idata, int depth);
		std::pair<bool,float> intersectR(const NodeRange& range, VNode* cnode, const IntersectData& idata, int depth); 

	public: 
		VolumeOctree(Volume* v, int depth){
			_root = NULL; 
			_depth = depth;
			_volume = v; 
			build();
		}
		~VolumeOctree() { if(_root) delete _root; }

		inline std::pair<bool,float> intersect(const osg::Vec3& pt, const osg::Vec3& dir, float radius, TransferFunction* tf, SegmentMap* sm, SelectionMap* selectm, Volume* mask, int color); 
	};

	typedef VolumeOctree VolumeIntersecter; 


	bool VolumeOctree::isBoxIntersecting(const NodeRange& range, const IntersectData& idata, int depth){
		float tol = .01/(1<<(depth+1));

		osg::Vec3 bmin(
			range.lx/idata.dimXm1*idata.sclx-idata.offsetx, 
			range.ly/idata.dimYm1*idata.scly-idata.offsety, 
			range.lz/idata.dimZm1*idata.sclz-idata.offsetz);
		osg::Vec3 bmax(
			range.ux/idata.dimXm1*idata.sclx-idata.offsetx, 
			range.uy/idata.dimYm1*idata.scly-idata.offsety, 
			range.uz/idata.dimZm1*idata.sclz-idata.offsetz); 
		// always assume axis algined box

		if(idata.dimXm1==0) {bmin[0]=0; bmax[0]=0;}
		if(idata.dimYm1==0) {bmin[1]=0; bmax[1]=0;}
		if(idata.dimZm1==0) {bmin[2]=0; bmax[2]=0;}

		if(idata.dir.z()!=0){
			// z=-.5 plane
			{
				float t = (bmin.z()- idata.pt.z())/idata.dir.z(); 
				float nx = idata.pt.x() + idata.dir.x()*t; 
				float ny = idata.pt.y() + idata.dir.y()*t; 

				if(nx<=bmax.x()+tol && nx>=bmin.x()-tol && ny<=bmax.y()+tol && ny>=bmin.y()-tol)
					return true; 
			}

			// z=.5 plane
			{
				float t = (bmax.z()-idata.pt.z())/idata.dir.z(); 
				float nx = idata.pt.x() + idata.dir.x()*t; 
				float ny = idata.pt.y() + idata.dir.y()*t; 

				if(nx<=bmax.x()+tol && nx>=bmin.x()-tol && ny<=bmax.y()+tol && ny>=bmin.y()-tol)
					return true; 
			}
		}

		if(idata.dir.x()!=0){
			// x=-.5 plane
			{
				float t = (bmin.x()-idata.pt.x())/idata.dir.x(); 
				float nz = idata.pt.z() + idata.dir.z()*t; 
				float ny = idata.pt.y() + idata.dir.y()*t; 

				if(nz<=bmax.z()+tol && nz>=bmin.z()-tol && ny<=bmax.y()+tol && ny>=bmin.y()-tol)
					return true; 
			}

			// x=.5 plane
			{
				float t = (bmax.x()-idata.pt.x())/idata.dir.x(); 
				float nz = idata.pt.z() + idata.dir.z()*t; 
				float ny = idata.pt.y() + idata.dir.y()*t; 

				if(nz<=bmax.z()+tol && nz>=bmin.z()-tol && ny<=bmax.y()+tol && ny>=bmin.y()-tol)
					return true; 
			}
		}

		if(idata.dir.y()!=0){
			// y=-.5 plane
			{
				float t = (bmin.y()-idata.pt.y())/idata.dir.y(); 
				float nz = idata.pt.z() + idata.dir.z()*t; 
				float nx = idata.pt.x() + idata.dir.x()*t; 

				if(nz<=bmax.z()+tol && nz>=bmin.z()-tol && nx<=bmax.x()+tol && nx>=bmin.x()-tol)
					return true; 
			}

			// y=.5 plane
			{
				float t = (bmax.y()-idata.pt.y())/idata.dir.y(); 
				float nz = idata.pt.z() + idata.dir.z()*t; 
				float nx = idata.pt.x() + idata.dir.x()*t; 

				if(nz<=bmax.z()+tol && nz>=bmin.z()-tol && nx<=bmax.x()+tol && nx>=bmin.x()-tol)
					return true; 
			}
		}

		return false;
	}

	std::pair<bool,float> VolumeOctree::intersect(const osg::Vec3& pt, const osg::Vec3& dir, float radius, TransferFunction* tf, SegmentMap* sm, SelectionMap* selectm, Volume* mask, int curColor){
		IntersectData idata; 
		idata.pt = pt; 
		idata.dir = dir; 
		idata.vmin = 0;  
		idata.vmax = 255; 
		idata.transFunc = tf; 
		idata.segMap = sm; 
		idata.selectMap = selectm; 
		idata.mask = mask; 
		idata.curId = curColor; 
		idata.r2 = radius * radius; 
		idata.dimXm1 = _volume->getDimX()-1; 
		idata.dimYm1 = _volume->getDimY()-1; 
		idata.dimZm1 = _volume->getDimZ()-1; 

		// we need the following to account for non-cube dimensions
		idata.dimmaxm1 = max(idata.dimXm1,max(idata.dimYm1,idata.dimZm1));
		//idata.offsetx = idata.dimXm1/idata.dimmaxm1/2; 
		//idata.offsety = idata.dimYm1/idata.dimmaxm1/2; 
		//idata.offsetz = idata.dimZm1/idata.dimmaxm1/2; 

		idata.offsetx = _volume->getMRCData()->offsetx; 
		idata.offsety = _volume->getMRCData()->offsety; 
		idata.offsetz = _volume->getMRCData()->offsetz; 

		idata.sclx = _volume->getMRCData()->scalex; 
		idata.scly = _volume->getMRCData()->scaley;
		idata.sclz = _volume->getMRCData()->scalez; 

		idata.dirlen = dir.length();

		//cout<<"pt: "<<pt.x()<<" "<<pt.y()<<" "<<pt.z()<<endl;
		//cout<<"dir: "<<dir.x()<<" "<<dir.y()<<" "<<pt.z()<<endl;

		NodeRange r; 
		r.lx=r.ly=r.lz = 0; 
		r.ux = _volume->getDimX()-1; 
		r.uy = _volume->getDimY()-1; 
		r.uz = _volume->getDimZ()-1; 

		return intersectR(r, _root, idata, _depth);
	}
}
#endif