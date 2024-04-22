#ifndef VPT_SELECTION_MAP_H
#define VPT_SELECTION_MAP_H

#include "Volume.h"
#include "Common.h"

#include <vector>

namespace vpt{

	class SelectionPlane{
	protected: 
		osg::Vec3 _normal; 
		osg::Vec3 _pt; 
		osg::Vec4 _color; 
		bool _active; 
		bool _visible; 
	public: 
		SelectionPlane(const osg::Vec3& norm, const osg::Vec3& p, bool act=true, bool visible=true) {
			_normal = norm; _pt = p; _active = act; _visible = visible; 
		}

		inline osg::Vec3 getNormal() const { return _normal; }
		inline osg::Vec3 getPoint() const { return _pt; }
		inline osg::Vec4 getColor() const { return _color; }
		inline const void setNormal(const osg::Vec3& n) { _normal = n; }
		inline const void setPoint(const osg::Vec3& p) { _pt = p; }
		inline const void setColor(const osg::Vec4& c) { _color = c; }
		inline bool isVisible() { return _visible; }
		inline bool isActive() { return _active; }
		inline void setActive(bool f) { _active = f; }
		inline void setVisible(bool v) { _visible = v; }
	};


	class SelectionMap{
	protected: 
		Volume* _selection; 
		Volume* _volume; 
		std::vector<SelectionPlane*> _planes; 

		int _curEditPlaneInd; 

	public: 
		SelectionMap(Volume* vol);
		~SelectionMap(); 
		inline Volume* getSelection() { return _selection; }
		inline void addPlane(SelectionPlane* plane) {_planes.push_back(plane); }
		inline void addPlanes(const std::vector<SelectionPlane*>& ps) { _planes = ps; }
		inline SelectionPlane* getPlane(int ind) { return _planes[ind]; }
		inline int getNumPlanes() { return (int) _planes.size(); }
		inline void removeAllPlanes() { _planes.clear(); }

		inline std::vector<SelectionPlane*> getPlanes() { return _planes; }

		inline void setCurEditPlane(SelectionPlane* sp){
			for(uint j=0;j<_planes.size();j++)
				if(_planes[j]==sp){
					_curEditPlaneInd = j; 
					break;
				}
		}

		inline SelectionPlane* getCurEditPlane() { 
			if(_curEditPlaneInd>=0 && _curEditPlaneInd<(int)_planes.size())
				return _planes[_curEditPlaneInd]; 
			return NULL; 
		}

	//	inline SelectionPlane* createSelectionPlane(); 
	};
}


#endif