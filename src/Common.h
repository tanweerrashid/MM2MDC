#ifndef VPT_COMMON_H
#define VPT_COMMON_H

// Common includes

#include <sstream>
#include <string> 
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include <osg/Matrix>
#include <osg/Vec3>
#include <osg/Vec4>

#include <omp.h>

// we typedef the voxel type here, hopefully changing it should be 
//  as easy as modifying this typedef
typedef unsigned int uint; 
typedef unsigned char vval; 

namespace vpt{
#define VV_ITEXTURE_FORMAT GL_LUMINANCE
#define VV_PIXEL_FORMAT GL_LUMINANCE
#define VV_DATA_FORMAT GL_UNSIGNED_BYTE
#define VV_MIN 0
#define VV_MAX 255
// since this is a special function, use all lower case
inline vval vval_scale(double v){ return (vval) (v*VV_MAX); }
inline double vval_inverse_scale(vval t) { return (t-VV_MIN)/((double)(VV_MAX-VV_MIN)); }

// floating point comparison
#define EPS32 1e-6
#define EPS64 1e-9

#define eq32(x,y) y-EPS32<x && y+EPS32>x


inline bool isSameFloat(float a, float b) {
    if (fabs(a - b) < 0.001) {
        //std::cout << "\n\nFloat difference between " << a << " and " << b << " is: " << fabs(a - b) << std::endl << std::endl;
        return true;
    }
    else {
        //std::cout << "\n\nFloat difference between " << a << " and " << b << " is: " << fabs(a - b) << std::endl << std::endl;
        return false;
    }
}

// string functions
inline std::string toLower(const std::string& in){
    std::string ret = ""; 
    for(uint j=0;j<in.size();j++){
        if('A'<=in[j] && 'Z'>=in[j])
            ret+=(in[j]-'A')+'a'; 
        else
            ret+=in[j]; 
    }
    return ret; 
}

inline std::string getFileExtension(const std::string& file){
    if(file.rfind(".")!=std::string::npos)
        return toLower(file.substr(file.rfind(".")+1));
    return ""; 
}

inline std::string getFilePath(const std::string& file){
	if(file.rfind("/")!=std::string::npos)
		return file.substr(0,file.rfind("/")+1);
	return ""; 
}

inline std::string getFileName(const std::string& file){
	if(file.rfind("/")!=std::string::npos)
		return file.substr(file.rfind("/")+1);
	return file; 
}

inline std::string getFileNameNoExt(const std::string& file){
	std::string ret = getFileName(file); 

    if(ret.rfind(".")!=std::string::npos)
        ret = ret.substr(0,ret.rfind(".")); 

	return ret;
}

inline std::vector<std::string> strSplit(const std::string& str, char delim){
	std::vector<std::string> ret; 
	std::string buffer =""; 

	for(uint j=0;j<str.length();j++){
		if(str[j]==delim){
			ret.push_back(buffer); 
			buffer =""; 
		}
		else
			buffer+=str[j]; 
	}
	if(buffer.length()>0)
		ret.push_back(buffer); 

	return ret; 
}

template<typename T>
void printMatrix(T& oldProj){
	std::cout<<"*************"<<std::endl;
	std::cout<<oldProj(0,0)<<" "<<oldProj(0,1)<<" "<<oldProj(0,2)<<" "<<oldProj(0,3)<<std::endl
		<<oldProj(1,0)<<" "<<oldProj(1,1)<<" "<<oldProj(1,2)<<" "<<oldProj(1,3)<<std::endl
		<<oldProj(2,0)<<" "<<oldProj(2,1)<<" "<<oldProj(2,2)<<" "<<oldProj(2,3)<<std::endl
		<<oldProj(3,0)<<" "<<oldProj(3,1)<<" "<<oldProj(3,2)<<" "<<oldProj(3,3)<<std::endl;
	std::cout<<"*************"<<std::endl;
}


template<typename T>
inline std::string numToStr(T num){
	std::stringstream ss; 
	ss<<std::fixed<<std::setprecision(1)<<num; 
	std::string ret; 
	ss>>ret; 
	return ret; 
}

template<typename T>
inline T strToNum(const std::string& str, T num){
	std::stringstream ss(str);
	T ret; 
	ss>>ret; 
	return ret; 
}

struct VertexData{
	VertexData(){
		verts = NULL; 
		inds = NULL; 
		isQuad = false;
		isTri = false;
		restart(); 
	}

	~VertexData(){	clear();}

	inline void clear(){
		if(verts!=NULL)
			delete [] verts; 
		if(inds!=NULL)
			delete [] inds; 
		verts = NULL;
		inds = NULL;
		isQuad = false;
		isTri = false;
	}

	inline void resize(int vs, int is){
		if(verts!=NULL)
			delete [] verts; 
		if(inds!=NULL)
			delete [] inds; 
		m_actualVSize = vs; 
		m_actualISize = is; 
		verts = new float[vs*3]; 
		inds = new unsigned int[is]; 
	}

	inline void restart(){
		isize = 0; 
		vsize = 0; 
	}

	inline int addVertex(float a, float b, float c){
		if(vsize>=m_actualVSize){
			std::cerr<<"exceeding vertex limit"<<std::endl;
			return -1; 
		}

		verts[vsize*3] = a; 
		verts[vsize*3+1] = b; 
		verts[vsize*3+2] = c; 
		return vsize++;
	}

	inline void addTriangle(int a, int b, int c){
		if(isize+2>=m_actualISize){
			std::cerr<<"exceeding index limit"<<std::endl;
			return; 
		}
		isTri = true;

		addIndex(a);
		addIndex(b);
		addIndex(c);
	}

	inline void addQuad(int a, int b, int c, int d){                                
            if(isize+3>=m_actualISize){
                std::cerr<<"exceeding index limit"<<std::endl;
                return; 
            }

            isQuad=true;

            addIndex(a);
            addIndex(b);
            addIndex(c);
            addIndex(d);                
        }
	

	inline void addIndex(int a){
            if(isize>=m_actualISize){
                std::cerr<<"exceeding index limit"<<std::endl;
                return; 
            }

            inds[isize++] = a; 
	}

	float* verts; 
	unsigned int* inds; 
	int isize,vsize; 
	bool isTri, isQuad;
private: 
	int m_actualVSize; 
	int m_actualISize; 
};

#define LARGE_FLOAT 10000000.f;
#define SMALL_FLOAT -10000000.f;

class MeshGeometry{
protected: 
	float* _verts; 
	float* _norms; // assume per vertex normal
	unsigned int* _inds; 
	
	int _icap, _vcap;

	float _minExtents[3]; 
	float _maxExtents[3]; 
	float _center[3]; 

public: 
	MeshGeometry(){
		_verts = NULL; 
		_inds = NULL; 
		_norms = NULL;

		_vsize = _vcap = 0; 
		_isize = _icap = 0; 

	}

	~MeshGeometry(){
		clear();
	}

        int _isize,_vsize; 
        
	void clear(){
		if(_verts!=NULL)
			delete [] _verts; 
		if(_norms!=NULL)
			delete [] _norms;
		if(_inds!=NULL)
			delete [] _inds; 
		_verts = NULL;
		_norms = NULL;
		_inds = NULL;
		_vsize = _vcap = 0; 
		_isize = _icap = 0; 
	}

        
	inline unsigned int addVertex(float* in){
		return addVertex(in[0],in[1],in[2]);
	}

	inline unsigned int addVertex(float a, float b, float c){
            if(_vsize>=_vcap)
                resizeVerts(2*_vcap+1); 

            _verts[_vsize*3] = a; 
            _verts[_vsize*3+1] = b; 
            _verts[_vsize*3+2] = c; 
            return _vsize++;
	}

	inline void addTriangle(int a, int b, int c){
        if(_isize+3>=_icap)
			resizeInds(2*_icap+1); 

		addIndex(a);
		addIndex(b);
		addIndex(c);
	}


	inline void addQuad(int a, int b, int c, int d){            
            if(_isize+4>=_icap)
                resizeInds(2*_icap+1); 

            addIndex(a);
            addIndex(b);
            addIndex(c);
            addIndex(d);
        }

	inline void getCenter(float* in) const { 
		in[0] = _center[0];
		in[1] = _center[1];
		in[2] = _center[2];
	}

	inline void getMin(float* in) const { 
		in[0] = _minExtents[0];
		in[1] = _minExtents[1];
		in[2] = _minExtents[2];
	}

	inline void getMax(float* in) const { 
		in[0] = _maxExtents[0];
		in[1] = _maxExtents[1];
		in[2] = _maxExtents[2];
	}


	inline float* getVerts()  { return _verts; }
	inline unsigned int* getInds()  { return _inds; }
	virtual float* getVert(unsigned int i) const { return &_verts[i*3]; }
	inline unsigned int getInd(int i) const { return _inds[i]; }
	virtual unsigned int* getTri(int i) const { return &_inds[i*3]; }
	virtual unsigned int* getQuad(int i) const { return &_inds[i*4]; }
	inline float* getNormal(unsigned int i) const  { return &_norms[i*3]; }
	inline float* getNormals() { return _norms; }

	inline int getNumVerts() const { return _vsize; }
	inline int getNumTris() const { return _isize/3; }
	inline int getNumInds() const { return _isize; }
	inline int getNumQuads() const { return _isize/4; }


	void cross(float* a, float* b, float* out){
		out[0] += a[1]*b[2] - a[2]*b[1]; 
		out[1] += -a[0]*b[2] + a[2]*b[0]; 
		out[2] += a[0]*b[1] - a[1]*b[0];
	}

        // Computes normals for triangles. 
	void computeNormals(){
		for(int j=0;j<getNumVerts()*3;j++)
			_norms[j] = 0; 

		for(int j=0;j<getNumTris();j++){
			unsigned int* tri = getTri(j);
			float* p0 = getVert(tri[0]);
			float* p1 = getVert(tri[1]);
			float* p2 = getVert(tri[2]);
			float a[3] = {p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]};
			float b[3] = {p2[0]-p0[0],p2[1]-p0[1],p2[2]-p0[2]};

			for(int k=0;k<3;k++)
				cross(a,b,&_norms[tri[k]*3]);
		}

		// normalize
		for(int j=0;j<getNumVerts();j++){
			float norm = sqrt(
				_norms[j*3]*_norms[j*3] +
				_norms[j*3+1]*_norms[j*3+1] + 
				_norms[j*3+2]*_norms[j*3+2]);

			_norms[j*3]/=norm;
			_norms[j*3+1]/=norm;
			_norms[j*3+2]/=norm;
		}
	}

        
        // Computes normals for quads. 
	void computeNormalsQ(){
		for(int j=0;j<getNumVerts()*3;j++)
			_norms[j] = 0; 

		for(int j=0;j<getNumQuads();j++){
			unsigned int* quad = getQuad(j);
			float* p0 = getVert(quad[0]);
			float* p1 = getVert(quad[1]);
			float* p2 = getVert(quad[2]);
			float a[3] = {p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]};
			float b[3] = {p2[0]-p0[0],p2[1]-p0[1],p2[2]-p0[2]};

			for(int k=0;k<4;k++)
				cross(a,b,&_norms[quad[k]*3]);
		}

		// normalize
		for(int j=0;j<getNumVerts();j++){
			float norm = sqrt(
				_norms[j*3]*_norms[j*3] +
				_norms[j*3+1]*_norms[j*3+1] + 
				_norms[j*3+2]*_norms[j*3+2]);

			_norms[j*3]/=norm;
			_norms[j*3+1]/=norm;
			_norms[j*3+2]/=norm;
                        
                        //std::cout << "Vertex " << j << " has normal value: " << _norms[j*3] << ", " << _norms[j*3+1] << ", " << _norms[j*3+2] << std::endl;
		}
	}

	void computeCenter(){
		_center[0] = _center[1] = _center[2] = 0; 
		int nverts = getNumVerts();

		for(int j=0;j<nverts;j++){
			float* p = getVert(j);

			_center[0]+=p[0];
			_center[1]+=p[1];
			_center[2]+=p[2];
		}
		_center[0]/=nverts;
		_center[1]/=nverts;
		_center[2]/=nverts;
	}

	void computeExtents(){
		for(int j=0;j<3;j++){
			_minExtents[j] = LARGE_FLOAT;
			_maxExtents[j] = SMALL_FLOAT;
		}

		for(int j=0;j<getNumVerts();j++){
			float* p = getVert(j);
			for(int k=0;k<3;k++){
//				_minExtents[k] = min(_minExtents[k],p[k]);
//				_maxExtents[k] = max(_maxExtents[k],p[k]);
			}
		}
	}

	// don't use this function
	void normalize(){
		computeExtents();

		float dim[3] = {
			_maxExtents[0] - _minExtents[0],
			_maxExtents[1] - _minExtents[1],
			_maxExtents[2] - _minExtents[2]}; 

//		float maxDim = max(dim[0],dim[1],dim[2]); 
//		float maxDim = max(dim[0],max(dim[1],dim[2])); 
		float maxDim; 

		for(int j=0;j<getNumVerts();j++){
			float* vert = getVert(j); 
			for(int k=0;k<3;k++)
				vert[k]=(vert[k]-_minExtents[k])/(maxDim)*2-1.f; 
		}

		computeCenter(); 

		for(int j=0;j<getNumVerts();j++){
			float* vert = getVert(j); 
			for(int k=0;k<3;k++)
				vert[k]-=_center[k];
		}
	}

//	void computeNormals(); 
//	void computeCenter(); 
//	void computeExtents();
//	void normalize(); 

protected: 


	void resizeVerts(int vs){
		float* oldv = NULL;
		if(_verts!=NULL)
			oldv = _verts; 

		_verts = new float[vs*3];
		if(oldv){
			for(int j=0;j<_vcap*3;j++)
				_verts[j] = oldv[j]; 
			delete [] oldv; 
		}

		float* oldn = NULL;
		if(_norms!=NULL)
			oldn = _norms; 

		_norms = new float[vs*3];
		if(oldn){
			for(int j=0;j<_vcap*3;j++)
				_norms[j] = oldn[j]; 
			delete [] oldn; 
		}

		_vcap = vs; 
	}

	void resizeInds(int is){
		unsigned int* oldi = NULL;
		if(_inds!=NULL)
			oldi = _inds; 

		_inds = new unsigned int[is];
		if(oldi){
			for(int j=0;j<_icap;j++)
				_inds[j] = oldi[j]; 
			delete [] oldi; 
		}

		_icap = is; 
	}

    inline void addIndex(int a){
        if(_isize>=_icap)
			resizeInds(2*_icap+1); 

        _inds[_isize++] = a; 
    }
};


#define NUM_PRESET_COLORS 42
const float PRESET_COLORS[NUM_PRESET_COLORS][3]=
{
    {0.90222, 0.101808, 0.198306}, {0.480598, 0.391471, 0.999921}, 
    {0.272896, 0.77160, 1.000000}, {0.414583, 0.979547, 0.278776}, 
    {0.986072, 0.957383, 0.394318}, {1.00000, 0.651215, 0.303137}, 
    {0.97894, 0.494649, 0.63177}, {0.926088, 0.224025, 0.333161}, 
    {0.569159, 0.49042, 0.999833}, {0.394252, 0.82298, 1.00000},
    {0.526328, 0.986979, 0.380287}, {0.992211, 0.965904, 0.484183}, 
    {1.000000, 0.706396, 0.414043}, {0.411014, 0.313725, 0.99999}, 
    {0.949957, 0.346242, 0.468017}, {0.657721, 0.589369, 0.999745}, 
    {0.515608, 0.87436, 1.00000}, {0.638074, 0.994411, 0.481799}, 
    {0.998349, 0.974424, 0.574048}, {0.907335, 0.127997, 0.227204}, 
    {0.499576, 0.412674, 0.999902}, {0.973825, 0.46846, 0.602872}, 
    {0.746282, 0.688319, 0.999657}, {0.636963, 0.925741, 1.00000}, 
    {0.973796, 0.940343, 0.214587}, {1., 0.540853, 0.0813252}, 
    {0.931203, 0.250215, 0.362059}, {0.588137, 0.511623, 0.999814}, 
    {0.997694, 0.590677, 0.737727}, {0.15154, 0.72022, 1.000000}, 
    {0.302837, 0.972116, 0.177265}, {0.979934, 0.948863, 0.304452}, 
    {1.000000, 0.596034, 0.192231}, {0.955071, 0.372432, 0.496914}, 
    {0.676698, 0.610573, 0.999726}, {0.480598, 0.391471, 0.999921}, 
    {0.272896, 0.771600, 1.00000}, {0.414583, 0.979547, 0.278776}, 
    {0.986072, 0.957383, 0.394318}, {1.00000, 0.651215, 0.303137}, 
    {0.97894, 0.494649, 0.63177}, {0.765259, 0.709522, 0.999638}}; 

}

#define GL_RGBA_FLOAT32 0x8814
#define GL_RGB_FLOAT32 0x8815

#endif