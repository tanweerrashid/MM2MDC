#ifndef VPT_UNDO_VOLUME_H
#define VPT_UNDO_VOLUME_H

#include "Volume.h"
#include <list> 

#define NUM_MAX_UNDOS 1

namespace vpt{
	struct UndoRange{
		int px,py,pz;
		int dx,dy,dz; 
	}; 

	class UndoVolume{
	protected: 
		int DATA_SIZE; 

		std::list<UndoRange> mUndoRangeHistory; 
		std::list<Volume*> mUndoHistory; 
		std::list<Volume*> mUnusedBuffer; 

	public: 
		UndoVolume(){}

		~UndoVolume(){
			std::vector<Volume*> l(mUndoHistory.begin(),mUndoHistory.end()); 
			for(uint j=0;j<l.size();j++)
				delete l[j]; 
		}


		inline UndoRange getUndoRange() { return mUndoRangeHistory.back(); }

		inline void resize(uint dx, uint dy, uint dz){
			Volume* vol= NULL; 
			for(int j=0;j<NUM_MAX_UNDOS;j++){
				 vol = new Volume(); 
				 vol->resize(dx,dy,dz); 
				 mUnusedBuffer.push_back(vol);
			}
		}

		// write the last edit into vol without actually popping 
		inline void peekUndo(Volume* vol){
/*
			if(getNumUndos()>0){
				Volume* undovol = mUndoHistory.back(); 
				UndoRange r = mUndoRangeHistory.back(); 
                #pragma omp parallel for
				for(int k=r.pz;k<r.pz+r.dz;k++)
					for(int j=r.py;j<r.py+r.dy;j++)
						for(int i=r.px;i<r.px+r.dx;i++){
							(*vol)(i,j,k) = (*undovol)(i,j,k); 
						}
			}
*/
		}

		inline void pushUndo(Volume* vol, const UndoRange& r){
/*
			uint dimx = vol->getDimX(); 
			uint dimy = vol->getDimY(); 
			uint dimz = vol->getDimZ(); 

			mUndoRangeHistory.push_back(r); 
			if(mUndoRangeHistory.size()>NUM_MAX_UNDOS)
				mUndoRangeHistory.pop_front(); 

			Volume* undovol = NULL; 
			if(mUnusedBuffer.size()==0){
				undovol = mUndoHistory.front(); 
				mUndoHistory.pop_front(); 
			}
			else{
				undovol = mUnusedBuffer.front(); 
				mUnusedBuffer.pop_front(); 
			}
			//undovol->clear(); 

            #pragma omp parallel for
			for(int k=r.pz;k<r.pz+r.dz;k++)
				for(int j=r.py;j<r.py+r.dy;j++)
					for(int i=r.px;i<r.px+r.dx;i++)
						(*undovol)(i,j,k) = (*vol)(i,j,k); 

			mUndoHistory.push_back(undovol); 
*/
		}

		inline void popUndo(Volume* vol, UndoRange& r){
/*
			if(getNumUndos()==0){
				r.dx = r.dy = r.dz = 0; 
				return; 
			}

			r = mUndoRangeHistory.back(); 
			mUndoRangeHistory.pop_back(); 

			Volume* undovol = mUndoHistory.back(); 
			mUndoHistory.pop_back(); 

            #pragma omp parallel for
			for(int k=r.pz;k<r.pz+r.dz;k++)
				for(int j=r.py;j<r.py+r.dy;j++)
					for(int i=r.px;i<r.px+r.dx;i++){
						(*vol)(i,j,k) = (*undovol)(i,j,k); 
					//	(*vol)(i,j,k) = 0;
					}

			mUnusedBuffer.push_back(undovol); 
*/
		}

		inline int getNumUndos() { return (int) mUndoHistory.size(); }
	};
}

#endif