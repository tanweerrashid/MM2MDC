#include "Session.h"
#include "VolumeUtils.h"
#include "TransferFunctionUtils.h"
#include "SegmentMapUtils.h"
#include "NormalMapUtils.h"

using namespace vpt; 
using namespace std; 

Session* Session::_instance = NULL;

Session* Session::inst() { 
	if(_instance==NULL)
		return _instance = new Session(); 
	return _instance; 
}

Volume* Session::createVolumeMask(Volume* vol, const std::string& name){
	Volume* ret = new Volume(); 
	ret->resize(vol->getDimX(), vol->getDimY(), vol->getDimZ());
	return ret; 
}


TransferFunction* Session::createTransferFunction(const std::string& name){
	TransferFunction* ret = new TransferFunction(); 
	return ret; 
}

SegmentMap* Session::createSegmentMap(const std::string& name){
	SegmentMap* ret = new SegmentMap(); 
	return ret; 
}

NormalMap* Session::createNormalMap(FVolume* vol, Volume* mask, SegmentMap* segmap, TransferFunction* tf){
/*
	NormalMap* ret = new NormalMap(vol,mask,segmap,tf); 
	return ret; 
*/
	return NULL;
}

