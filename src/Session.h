#ifndef VPT_SESSION_H
#define VPT_SESSION_H

#include <map>
#include <list>
#include <string>

#include "Project.h"

#include "TransferFunction.h"
#include "SegmentMap.h"
#include "NormalMap.h"

#include <osg/Texture3D>
#include <osg/Texture1D>

#define make_normal_key(a,b,c,d) std::make_pair(a,std::make_pair(b,std::make_pair(c,d)))


namespace vpt{

	typedef std::vector<vpt::Project*> ProjectManager; 

    class Session{
	private: 
		static Session* _instance; 

    protected: 
		ProjectManager _projectManager; 

    public: 
		static Session* inst(); 

		Volume* createVolumeMask(Volume* vol, const std::string& name); 
		TransferFunction* createTransferFunction(const std::string& name); 
		SegmentMap* createSegmentMap(const std::string& name); 
		NormalMap* createNormalMap(FVolume* vol, Volume* mask, SegmentMap* segmap, TransferFunction* tf); 

		Project* createNewProject() { return NULL; }
		Project* openProjectFromFile(const std::string& filename) { return NULL; }

	protected: 
		Session() {}
    }; 
}

#endif