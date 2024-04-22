#ifndef TRANSFER_FUNCTION_H
#define TRANSFER_FUNCTION_H

#include <set> 
#include <cmath>
#include <algorithm>
#include <vector>

#include "Common.h"
#include "Savable.h"

#define TF_RESOLUTION 256

using namespace std; 

namespace vpt{
	// the isovalue is always between 0 and 1.  

	class TransferFunctionDataPoint{
	protected: 
		float _isovalue, _alpha; 
	public: 
		TransferFunctionDataPoint(float isovalue, float alpha)
			:_isovalue(isovalue), _alpha(alpha){}
		inline float getIsoValue() { return _isovalue; }
		inline float getAlpha() { return _alpha; }

		inline void setIsoValue(float nval) { _isovalue = nval; }
		inline void setAlpha(float nval) { _alpha = nval; }
	};

	typedef TransferFunctionDataPoint TFDPoint; 
	typedef std::pair<float,TFDPoint*> SortedPoint; 
	typedef std::vector<SortedPoint> SortedPoints; 

	class TransferFunction : public Savable {
	protected: 
		SortedPoints _sortedPts; 
		float _coarseCache[TF_RESOLUTION]; 
	public: 

		// pin two data points at the beginning and the end of the function
		TransferFunction(bool initBegEnd=true){
			//_sortedPts.insert(std::make_pair(0,new TFDPoint(0,0))); 
			//_sortedPts.insert(std::make_pair(1,new TFDPoint(1,1)));
			if(initBegEnd){
				createDataPoint(0,0); 
				createDataPoint(1,1); 
			}
		} 

		~TransferFunction(){
			for(SortedPoints::iterator it = _sortedPts.begin();it!=_sortedPts.end();it++)
				delete (*it).second; 
			_sortedPts.clear();
		}

		// samples the transfer function at a particular isovalue
		inline float sample(float iso){
			if(_sortedPts.size() < 2) return 0; 

			int st = 0; 
			int ed = _sortedPts.size()-1; 
			int mid = (st+ed)/2; 

			while(st+1<ed){
				int mid = (st+ed)/2; 
				if(_sortedPts[mid].first+EPS32 > iso){
					ed = mid; 
				}
				else if(_sortedPts[mid].first-EPS32 < iso){
					st = mid; 
				}
				else{
					st = mid; 
					break;
				}
			}

			mid = st; 

			if(mid==_sortedPts.size()-1) return 0; 

			float range = _sortedPts[mid+1].first - _sortedPts[mid].first; 
			float t = (iso - _sortedPts[mid].first) / range; // normalize; 

			float prevAlpha = _sortedPts[mid].second->getAlpha(); 
			float nextAlpha = _sortedPts[mid+1].second->getAlpha(); 

			// linear interpolation
			float res =  prevAlpha + (nextAlpha-prevAlpha)*t; 

			res = max(0.f,min(1.f,res)); 
			int ndiv = 20; 
			float div = 1.f/ndiv; 
			res = (pow(ndiv,res)/ndiv-div)/(1-div);
			return res; 
		}

		inline float sampleCoarse(int iso){
			if(iso<0 || iso>=TF_RESOLUTION) return 0;
			return _coarseCache[iso]; 
		}

		inline void updateFunction(){
			// resort the points
			sort(_sortedPts.begin(),_sortedPts.end()); 

			// build the coarse cache
			for(int j=0;j<TF_RESOLUTION;j++){
				float iso = j/1.f/(TF_RESOLUTION-1); 
				_coarseCache[j] = sample(iso); 
			}
		}

		inline TFDPoint* getDataPoint(int ind){
			if((uint) ind < _sortedPts.size())
				return _sortedPts[ind].second; 
			return NULL; 
		}

		inline int getNumDataPoints() { return (int) _sortedPts.size(); }

		inline TFDPoint* createDataPoint(float iso,float alpha){
			if(iso< -EPS32 || iso> 1+EPS32)
				return NULL; 

			TFDPoint* ret =  new TFDPoint(iso,alpha); 
			_sortedPts.push_back(std::make_pair(iso,ret)); 

			updateFunction(); 

			return ret; 
		}	   

		inline void removeDataPoint(TFDPoint* p){
			for(std::vector<std::pair<float,TFDPoint*> >::iterator i = _sortedPts.begin(); 
				i!=_sortedPts.end(); i++){
					if((*i).second == p){
						_sortedPts.erase(i); 
						break;
					}
			}

			updateFunction(); 
		}

		inline void moveTFDPoint(TFDPoint* p, float nval, float nalpha){
			for(uint j=0;j<_sortedPts.size();j++){
				if(_sortedPts[j].second==p){
					_sortedPts[j] = std::make_pair(nval,p); 
					break;
				}
			}

			p->setIsoValue(nval); 
			p->setAlpha(nalpha); 

			updateFunction(); 
		}
	}; 
}

#endif