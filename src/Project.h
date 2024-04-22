#ifndef VPT_PROJECT_H
#define VPT_PROJECT_H

#include "Volume.h"
#include "TransferFunction.h"
#include "SegmentMap.h"
#include "NormalMap.h"
#include "UndoVolume.h"
#include "BoundaryMap.h"

#include <string> 

namespace vpt{
	class Project : public Savable {
	protected: 
		Volume* _volume; 
		Volume* _mask; 
		TransferFunction* _transFunc; 
		SegmentMap* _segmentMap; 

		// The following are not savable 
		VolumeData<float>* _fvolume; 
		SelectionMap* _selectionMap;
		NormalMap* _normalMap; 
		Volume* _skeleton; 
		UndoVolume* _undoVolume; 
		BoundaryMap* _boundaryMap; 

	public: 

		Project(){
			_volume = NULL; 
			_fvolume = NULL; 
			_mask = NULL; 
			_transFunc = NULL; 
			_segmentMap = NULL; 
			_selectionMap = NULL; 
			_normalMap = NULL; 
			_undoVolume = NULL; 
			_boundaryMap = NULL;
		}

		// don't really have to do anything
		~Project() { }

		inline void setSegmentMap(SegmentMap* sg){	_segmentMap = sg; }
		inline void setMask(Volume* m) { _mask = m; }
		inline void setVolume(Volume* v) { _volume = v; }
		inline void setFVolume(VolumeData<float>* v) { _fvolume = v; }
		inline void setTransFunc(TransferFunction* tf) { _transFunc = tf; }
		inline void setSelectionMap(SelectionMap* select) { _selectionMap = select; }
		inline void setSkeleton(Volume* v) { _skeleton = v; }
		inline void setNormalMap(NormalMap* nm) { _normalMap = nm; }
		inline void setUndoVolume(UndoVolume* uv) { _undoVolume = uv; }
		inline void setBoundaryMap(BoundaryMap* bm) { _boundaryMap = bm; }

		inline Volume* getMask() { return _mask; }
		inline Volume* getVolume() { return _volume; }
		inline VolumeData<float>* getFVolume() { return _fvolume; }
		inline Volume* getSkeleton() { return _skeleton; }
		inline TransferFunction* getTransFunc() { return _transFunc; }
		inline SegmentMap* getSegmentMap() { return _segmentMap; }
		inline SelectionMap* getSelectionMap() { return _selectionMap; }
		inline NormalMap* getNormalMap() { return _normalMap; }
		inline UndoVolume* getUndoVolume() { return _undoVolume; }
		inline BoundaryMap* getBoundaryMap() { return _boundaryMap; }
	}; 
}


#endif