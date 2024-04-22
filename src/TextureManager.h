#ifndef VPT_TEXTURE_MANAGER_H
#define VPT_TEXTURE_MANAGER_H

#include <map>
#include <list>
#include <string>

#include <osg/Texture3D> 
#include <osg/Texture1D>

#include "Volume.h"

#include "TransferFunction.h"
#include "SegmentMap.h"
#include "NormalMap.h"
#include "SelectionMap.h"
#include "CombinedVolume.h"
#include "BoundaryMap.h"

namespace vpt{
	//typedef std::pair<vpt::Volume*,vpt::NormalMap*> VolumeAndNormal;
	//typedef std::map<VolumeAndNormal,osg::ref_ptr<osg::Texture3D> > VolumeAndNormalTexture3DManager;
	typedef std::map<vpt::Volume*,osg::ref_ptr<osg::Texture3D> > VolumeTexture3DManager; 
	typedef std::map<vpt::FVolume*,osg::ref_ptr<osg::Texture3D> > FVolumeTexture3DManager; 
	typedef std::map<vpt::NormalMap*,osg::ref_ptr<osg::Texture3D> > NormalMapTex3DManager; 
	typedef std::map<vpt::SegmentMap*,osg::ref_ptr<osg::Texture1D> > SegmentMapTex1DManager; 
	typedef std::map<vpt::Volume*,osg::ref_ptr<osg::Texture3D> > Texture3DManager; 
	typedef std::map<vpt::TransferFunction*,osg::ref_ptr<osg::Texture1D> > TransFuncTex1DManager; 
	typedef std::map<vpt::SelectionMap*,osg::ref_ptr<osg::Texture3D> > SelectionMapTex3DManager;
	typedef std::map<vpt::BoundaryMap*,osg::ref_ptr<osg::Texture3D> > BoundaryMapCornersTex3DManager; 
	typedef std::map<vpt::BoundaryMap*,osg::ref_ptr<osg::Texture3D> > BoundaryMapVolumeTex3DManager; 

	// a global manager for all textures
	class TextureManager{
	public: 
		static TextureManager* _instance; 
	protected: 
		VolumeTexture3DManager _volTex3DManager;
		FVolumeTexture3DManager _fvolTex3DManager;
		TransFuncTex1DManager _tfTex1DManager; 
		SegmentMapTex1DManager _segMapTex1DManager; 
		NormalMapTex3DManager _normalMapTex3DManager; 
		SelectionMapTex3DManager _selectionMapTex3DManager; 
		BoundaryMapCornersTex3DManager _boundaryMapCornersTex3DManager; 
		BoundaryMapVolumeTex3DManager _boundaryMapVolumeTex3DManager; 
//		VolumeAndNormalTexture3DManager _volAndNormTex3DManager;
		std::map<vpt::Volume*,vpt::CombinedVolume<vval>*> _volCombinedVolManager;

	public: 
		static TextureManager* inst(); 

		// this must be called before any other getOrCreate calls 
		osg::ref_ptr<osg::Texture3D> createInitialVolumeNormalCombinedTexture(Volume* vol, NormalMap* nm); 

		osg::ref_ptr<osg::Texture3D> getOrCreateVolumeTex3D(Volume* vol);
		osg::ref_ptr<osg::Texture3D> getOrCreateVolumeTex3D(FVolume* vol);
		inline bool hasVolumeTex3D(Volume* vol) { return _volTex3DManager.find(vol)!=_volTex3DManager.end(); }
		inline bool hasVolumeTex3D(FVolume* fvol) { return _fvolTex3DManager.find(fvol)!=_fvolTex3DManager.end(); }
		void removeVolumeTex3D(Volume* vol); 

		osg::ref_ptr<osg::Texture1D> getOrCreateTransFuncTex1D(TransferFunction* tf); 
		inline bool hasTransFuncTex1D(TransferFunction* tf){ return _tfTex1DManager.find(tf)!=_tfTex1DManager.end(); }
		void removeTransFuncTex1D(TransferFunction* tf); 

		osg::ref_ptr<osg::Texture1D> getOrCreateSegmentMapTex1D(SegmentMap* smp); 
		inline bool hasSegmentMapTex1D(SegmentMap* smp) { return _segMapTex1DManager.find(smp)!=_segMapTex1DManager.end(); }
		void removeSegmentMaptex1D(SegmentMap* smp); 

		osg::ref_ptr<osg::Texture3D> getOrCreateNormalMapTex3D(NormalMap* nm); 
		inline bool hasNormalMapTex3D(NormalMap* nm) { return _normalMapTex3DManager.find(nm)!=_normalMapTex3DManager.end(); }
		void removeNormalMapTex3D(NormalMap* nm); 

		osg::ref_ptr<osg::Texture3D> getOrCreateSelectionMapTex3D(SelectionMap* sm); 
		inline bool hasSelectionMapTex3D(SelectionMap* sm) { return _selectionMapTex3DManager.find(sm)!=_selectionMapTex3DManager.end(); }
		void removeSelectionMapTex3D(SelectionMap* sm); 

		osg::ref_ptr<osg::Texture3D> getOrCreateBoundaryMapCornersTex3D(BoundaryMap* bm); 
		inline bool hasBoundaryMapCornersTex3D(BoundaryMap* bm) { return _boundaryMapCornersTex3DManager.find(bm)!=_boundaryMapCornersTex3DManager.end(); }
		void removeBoundaryMapCornersTex3D(BoundaryMap* bm); 

		osg::ref_ptr<osg::Texture3D> getOrCreateBoundaryMapVolumeTex3D(BoundaryMap* bm); 
		inline bool hasBoundaryMapVolumeTex3D(BoundaryMap* bm) { return _boundaryMapVolumeTex3DManager.find(bm)!=_boundaryMapVolumeTex3DManager.end(); }
		void removeBoundaryMapVolumeTex3D(BoundaryMap* bm); 

		inline CombinedVolume<vval>* getCombinedVolume(Volume* v){
			if(_volCombinedVolManager.find(v)!=_volCombinedVolManager.end())
				return _volCombinedVolManager[v]; 
			return NULL; 
		}


	protected: 
		TextureManager() {}
	}; 

}

#endif