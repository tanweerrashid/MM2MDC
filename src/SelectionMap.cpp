#include "SelectionMap.h"
#include "Volume.h"

using namespace vpt; 

SelectionMap::SelectionMap(Volume* vol){
	_volume = vol;

	_selection = new Volume(); 
	_selection->resize(vol->getDimX(), vol->getDimY(), vol->getDimZ()); 

	for(uint j=0;j<_selection->getSize();j++)
		(*_selection)(j)=1;

	_curEditPlaneInd = -1; 
}

SelectionMap::~SelectionMap(){
	delete _selection; 
}
