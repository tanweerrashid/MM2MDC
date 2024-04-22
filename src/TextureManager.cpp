#include "TextureManager.h"

#include "VolumeUtils.h"
#include "TransferFunctionUtils.h"
#include "SegmentMapUtils.h"
#include "NormalMapUtils.h"
#include "SelectionMapUtils.h"
#include "CombinedVolumeUtils.h"
#include "BoundaryMapUtils.h"

#include <vector>

using namespace vpt; 
using namespace std; 

TextureManager* TextureManager::_instance = NULL;

TextureManager* TextureManager::inst() { 
	if(_instance==NULL)
		return _instance = new TextureManager(); 
	return _instance; 
}

osg::ref_ptr<osg::Texture3D> TextureManager::createInitialVolumeNormalCombinedTexture(Volume* vol, NormalMap* nm){
	osg::ref_ptr<osg::Texture3D> ret = new osg::Texture3D();
	ret->setResizeNonPowerOfTwoHint(false); 

	std::vector<Volume*> tv; 
	tv.push_back(vol); 
	tv.push_back(nm->getNormal());

	CombinedVolume<vval>* cv = new CombinedVolume<vval>(tv);
	osg::ref_ptr<osg::Image> img = CombinedVolumeUtils::convertToOSGImage(cv);
	ret->setImage(img.get());

	_volTex3DManager[vol] = ret; 
	_normalMapTex3DManager[nm] = ret; 

	//_volCombinedVolManager[vol] = cv; 
	//_volCombinedVolManager[nm->getNormal()] = cv; 

	return ret; 
}

osg::ref_ptr<osg::Texture3D> TextureManager::getOrCreateVolumeTex3D(Volume* vol){
	if(hasVolumeTex3D(vol))
		return _volTex3DManager[vol];

	osg::ref_ptr<osg::Texture3D> ret = new osg::Texture3D(); 
	ret->setResizeNonPowerOfTwoHint(false);
	osg::ref_ptr<osg::Image> img = VolumeUtils::convertToOSGImage(*vol); 
	ret->setImage(img.get());

	return _volTex3DManager[vol] = ret; 
}

osg::ref_ptr<osg::Texture3D> TextureManager::getOrCreateVolumeTex3D(FVolume* fvol){
	if(hasVolumeTex3D(fvol))
		return _fvolTex3DManager[fvol];

	osg::ref_ptr<osg::Texture3D> ret = new osg::Texture3D(); 
	ret->setResizeNonPowerOfTwoHint(false);
	osg::ref_ptr<osg::Image> img = VolumeUtils::convertToOSGImage(*fvol); 
	ret->setImage(img.get());

	return _fvolTex3DManager[fvol] = ret; 
}

void TextureManager::removeVolumeTex3D(Volume* vol){
	if(hasVolumeTex3D(vol)){
		_volTex3DManager[vol] = NULL; // for ref pointers?
		_volTex3DManager.erase(_volTex3DManager.find(vol)); 
		if(_volCombinedVolManager.find(vol)!=_volCombinedVolManager.end()){
			CombinedVolume<vval>* cv = _volCombinedVolManager[vol] ; 
			for(int j=0;j<cv->getNumVolumes();j++){
				Volume* v = cv->getVolume(j); 
				_volCombinedVolManager.erase(_volCombinedVolManager.find(v)); 
			}
			delete cv; 
		}
	}
}

osg::ref_ptr<osg::Texture1D> TextureManager::getOrCreateTransFuncTex1D(TransferFunction* tf){
	if(hasTransFuncTex1D(tf))
		return _tfTex1DManager[tf];

	osg::ref_ptr<osg::Texture1D> ret = new osg::Texture1D(); 
	ret->setResizeNonPowerOfTwoHint(false);
	osg::ref_ptr<osg::Image> img = TransferFunctionUtils::convertToOSGImage(tf); 
	ret->setImage(img.get());

	return _tfTex1DManager[tf] = ret; 
}


void TextureManager::removeTransFuncTex1D(TransferFunction* tf){
	if(hasTransFuncTex1D(tf)){
		_tfTex1DManager[tf] = NULL; 
		_tfTex1DManager.erase(_tfTex1DManager.find(tf)); 
	}
}

osg::ref_ptr<osg::Texture1D> TextureManager::getOrCreateSegmentMapTex1D(SegmentMap* sg){
	if(hasSegmentMapTex1D(sg))
		return _segMapTex1DManager[sg];

	osg::ref_ptr<osg::Texture1D> ret = new osg::Texture1D(); 
	ret->setResizeNonPowerOfTwoHint(false);
	osg::ref_ptr<osg::Image> img = SegmentMapUtils::convertToOSGImage(sg); 
	ret->setImage(img.get());

	return _segMapTex1DManager[sg] = ret; 
}

void TextureManager::removeSegmentMaptex1D(SegmentMap* smp){
	if(hasSegmentMapTex1D(smp)){
		_segMapTex1DManager[smp] = NULL; 
		_segMapTex1DManager.erase(_segMapTex1DManager.find(smp)); 
	}
}

osg::ref_ptr<osg::Texture3D> TextureManager::getOrCreateNormalMapTex3D(NormalMap* nm){
	if(hasNormalMapTex3D(nm))
		return _normalMapTex3DManager[nm]; 

	osg::ref_ptr<osg::Texture3D> ret = new osg::Texture3D(); 
	ret->setResizeNonPowerOfTwoHint(false);
	osg::ref_ptr<osg::Image> img = NormalMapUtils::convertToOSGImage(nm); 
	ret->setImage(img.get()); 

	return _normalMapTex3DManager[nm] = ret; 
}

void TextureManager::removeNormalMapTex3D(NormalMap* nm){
	if(hasNormalMapTex3D(nm)){
		_normalMapTex3DManager[nm] = NULL; 
		_normalMapTex3DManager.erase(_normalMapTex3DManager.find(nm)); 
	}
}

osg::ref_ptr<osg::Texture3D> TextureManager::getOrCreateSelectionMapTex3D(SelectionMap* sm){
	if(hasSelectionMapTex3D(sm))
		return _selectionMapTex3DManager[sm]; 

	osg::ref_ptr<osg::Texture3D> ret = new osg::Texture3D(); 
	ret->setResizeNonPowerOfTwoHint(false);
	osg::ref_ptr<osg::Image> img = SelectionMapUtils::convertToOSGImage(sm); 
	ret->setImage(img.get()); 

	return _selectionMapTex3DManager[sm] = ret; 
}

void TextureManager::removeSelectionMapTex3D(SelectionMap* sm){
	if(hasSelectionMapTex3D(sm)){
		_selectionMapTex3DManager[sm] = NULL; 
		_selectionMapTex3DManager.erase(_selectionMapTex3DManager.find(sm)); 
	}
}

osg::ref_ptr<osg::Texture3D> TextureManager::getOrCreateBoundaryMapCornersTex3D(BoundaryMap* bm){
	if(hasBoundaryMapCornersTex3D(bm))
		return _boundaryMapCornersTex3DManager[bm]; 

	osg::ref_ptr<osg::Texture3D> ret = new osg::Texture3D(); 
	ret->setResizeNonPowerOfTwoHint(false); 
	osg::ref_ptr<osg::Image> img = BoundaryMapUtils::convertCornersToOSGImage(bm); 
	ret->setImage(img.get()); 

	return _boundaryMapCornersTex3DManager[bm] = ret; 
}

void TextureManager::removeBoundaryMapCornersTex3D(BoundaryMap* bm){
	if(hasBoundaryMapCornersTex3D(bm)){
		_boundaryMapCornersTex3DManager[bm] = NULL; 
		_boundaryMapCornersTex3DManager.erase(_boundaryMapCornersTex3DManager.find(bm)); 
	}
}

osg::ref_ptr<osg::Texture3D> TextureManager::getOrCreateBoundaryMapVolumeTex3D(BoundaryMap* bm){
	if(hasBoundaryMapVolumeTex3D(bm))
		return _boundaryMapVolumeTex3DManager[bm]; 

	osg::ref_ptr<osg::Texture3D> ret = new osg::Texture3D(); 
	ret->setResizeNonPowerOfTwoHint(false); 
	osg::ref_ptr<osg::Image> img = BoundaryMapUtils::convertVolumeToOSGImage(bm); 
	ret->setImage(img.get()); 

	return _boundaryMapVolumeTex3DManager[bm] = ret; 
}

void TextureManager::removeBoundaryMapVolumeTex3D(BoundaryMap* bm){
	if(hasBoundaryMapVolumeTex3D(bm)){
		_boundaryMapVolumeTex3DManager[bm] = NULL; 
		_boundaryMapVolumeTex3DManager.erase(_boundaryMapVolumeTex3DManager.find(bm)); 
	}
}
