#ifndef VPT_SEGMENT_MAP_H
#define VPT_SEGMENT_MAP_H

#include <string>
#include <osg/Vec3>
#include <map>
#include <osg/Texture1D>

#include <iostream>

#include "Savable.h"

#include "QtCommon.h"
#include "Common.h"

#define NUM_SEGMENTS 256

namespace vpt{
	const int SHOW_SEGMENT = 1; 
	const int HIDE_SEGMENT = 0; 
	const int SEMI_HIDE_SEGMENT = 2; 

	class SegmentMapEntry {
	protected: 
		// these three are mutable
		std::string _name; 
		osg::Vec4 _color; 
		int _state; // either on or off
		float _trans; 

		// id should not be mutable
		int _id; 
	public:

		SegmentMapEntry(int id) 
			: _id(id), _state(1){
				_name = std::string("color-") + numToStr(id); 
				_trans = 1; 
		}

		inline std::string getName() { return _name; }
		inline void setName(const std::string& n) { _name = n; }

		inline osg::Vec4 getColor() { return _color; }
		inline void setColor(const osg::Vec4& col) {
			//std::cout<<"setting: "<<col.x()<<" "<<col.y()<<" "<<col.z()<<std::endl;
			_color = col; 
		}

		inline int getState() { return _state; }
		inline void setState(int st) { _state = st; }

		inline float getTransparency() { return _trans; }
		inline void setTransparency(float v) { _trans = v; }

		inline int getSegmentCode() { return _id; }
		inline void setSegmentCode(int id) { _id = id; }
	};

	class SegmentMap : public Savable{
	protected: 
		SegmentMapEntry* _segmentMapEntries[NUM_SEGMENTS]; 
		bool _reserved[NUM_SEGMENTS]; 
		int _num; 

	public: 
		static const int SELECTION_SEGMENT = 1; 
		static const int BASE_SEGMENT = 0; 

	public: 
		SegmentMap(bool initBaseAndSelection=true) : _num(0){
			for(int j=0;j<NUM_SEGMENTS;j++)
				_segmentMapEntries[j] = NULL; 

			for(int j=0;j<NUM_SEGMENTS;j++)
				_reserved[j] = false; 

			if(initBaseAndSelection){
				SegmentMapEntry* zeroth = createSegmentMapEntry(false); 
				zeroth->setColor(osg::Vec4(.5,.5,.5,.1));
				zeroth->setName(std::string("Base Segment")); 
				zeroth->setState(vpt::SHOW_SEGMENT); 

				SegmentMapEntry* first = createSegmentMapEntry(false); 
				first->setColor(osg::Vec4(.6,0,0,1)); 
				first->setName(std::string("Current Selection")); 
				first->setState(vpt::HIDE_SEGMENT);
			}
		}

		~SegmentMap(){
			for(int j=0;j<_num;j++){
				if(_segmentMapEntries[j]){
					delete _segmentMapEntries[j]; 
					_segmentMapEntries[j] = NULL; 
				}
			}
		}

		inline void clearEntries(){
			for(int j=0;j<_num;j++){
				if(j!=BASE_SEGMENT && j!=SELECTION_SEGMENT){
					delete _segmentMapEntries[j]; 
					_segmentMapEntries[j] = NULL; 
				}
			}
			_num = 2; 
		}

		SegmentMapEntry* createSegmentMapEntry(bool reserved=false) { 
			_segmentMapEntries[_num] = new SegmentMapEntry(_num); // the two numbers HAVE TO MATCH UP
			_reserved[_num] = reserved; 
			return _segmentMapEntries[_num++]; 
		}

		inline int getNumSegments() { return _num; }

		// return null for reserved entries
		inline SegmentMapEntry* getEntry(int id) {
			if(id>=NUM_SEGMENTS || _reserved[id]) return NULL; 
			return _segmentMapEntries[id]; 
		}

		// for now, the segment code and it's entry index is the same
		inline SegmentMapEntry* getEntryFromCode(int id){
			if(_reserved[id]) return NULL; 
			return _segmentMapEntries[id]; 
		}
		// use this if you know that it is not a reserved entry
		inline SegmentMapEntry* getEntryFast(int id){
			return _segmentMapEntries[id]; 
		}
	}; 
}

#endif