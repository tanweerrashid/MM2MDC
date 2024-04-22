#ifndef VPT_PROJECTION_UTILS_H
#define VPT_PROJECTION_UTILS_H

#include "Project.h"

#include "TransferFunctionUtils.h"
#include "SegmentMapUtils.h"
#include "VolumeUtils.h"
#include "BoundaryMap.h"
#include "BoundaryMapUtils.h"

#include <fstream> 
#include <iostream> 
#include <string>

#include <QDir>

namespace vpt{

	class ProjectUtils{
	protected: 

	public: 

		static std::pair<std::string,std::string> getFileNameAndPath(const std::string& line, const QDir& projDir){
			QString file(getFileName(line).c_str());
			QString path(getFilePath(line).c_str());

			if(QDir::isRelativePath(path))
				return make_pair(file.toStdString(), path.toStdString());

			QString rel = projDir.relativeFilePath(path); 
			if(rel.length()!=0 && rel[rel.length()-1]!='/') rel+="/";

			return make_pair(file.toStdString(), rel.toStdString());
		}

		static void writeToFile(Project* p){
			std::string SECRET_FILE_CODE = "NULL_1234234";
			std::string path = p->getFilePath() + p->getFileName(); 
			std::fstream file(path.c_str(), std::fstream::out); 

			Volume* vol = p->getVolume();
			if(vol){
				file<<vol->getFilePath()<<vol->getFileName()<<endl;
				cout<<"vol: "<<vol->getFilePath()<<vol->getFileName()<<endl;
			}				
			else
				file<<SECRET_FILE_CODE<<endl;

			Volume* mask = p->getMask();
			if(mask){
				file<<mask->getFilePath()<<mask->getFileName()<<endl;
				cout<<"mask: "<<mask->getFilePath()<<mask->getFileName()<<endl;
			}
			else
				file<<SECRET_FILE_CODE<<endl;

			TransferFunction* tf = p->getTransFunc();
			if(tf)
				file<<tf->getFilePath()<<tf->getFileName()<<endl;
			else
				file<<SECRET_FILE_CODE<<endl;

			SegmentMap* sm = p->getSegmentMap();
			if(sm)
				file<<sm->getFilePath()<<sm->getFileName()<<endl;
			else
				file<<SECRET_FILE_CODE<<endl;

			file.close();
		}


		static Project* newProjectFromVolume(const std::string& vfilename){
			Volume* vol = VolumeUtils::readByteVolumeFromFile(vfilename);
			VolumeData<float>* fvol = VolumeUtils::readFloatVolumeFromFile(vfilename);

			return newProjectFromVolume(vol,fvol); 
		}


		static Project* newProjectFromVolume(Volume* vol, VolumeData<float>* fvol){
			Project* ret = new Project(); 
			ret->setFileName("Untitled-Project.vpj");

			//Volume* vol = VolumeUtils::readByteVolumeFromFile(vfilename);
			ret->setVolume(vol); 
			//VolumeData<float>* fvol = VolumeUtils::readFloatVolumeFromFile(vfilename);
			ret->setFVolume(fvol);

			Volume* mask = new Volume(); 
			mask->resize(vol->getDimX(),vol->getDimY(),vol->getDimZ());
			ret->setMask(mask); 
			mask->setFileName("Untitled-mask.vmask");
			mask->setFilePath("");

			TransferFunction* tf = new TransferFunction();
			ret->setTransFunc(tf); 
			tf->setFileName("Untitled-transfunc.vtrans");
			tf->setFilePath("");

			SegmentMap* smap = new SegmentMap();
			ret->setSegmentMap(smap); 
			smap->setFileName("Untitled-segmap.vsmap");
			smap->setFilePath("");

			UndoVolume* undo = new UndoVolume(); 
			undo->resize(vol->getDimX(),vol->getDimY(),vol->getDimZ()); 
			ret->setUndoVolume(undo);

			BoundaryMap* bmap = new BoundaryMap(vol); 
			bmap->setProject(ret); 
			ret->setBoundaryMap(bmap);

			// set the project path to the path of the volume
			ret->setFilePath(vol->getFilePath());

			vol->setFilePath(""); 
			fvol->setFilePath(""); 

			return ret; 
		}


		// filename here should be absolute?
		static Project* readFromFile(const std::string& filename, bool autoGen=false){
			std::fstream file(filename.c_str(), std::fstream::in); 

			if(file.fail()) return NULL; 

			Project* ret = new Project(); 
			std::string projPath = getFilePath(filename);
			ret->setFileName(getFileName(filename)); 
			ret->setFilePath(projPath);
			ret->setSavableState(SAVABLE_OPENED_FILE);

			QString prevDir = QDir::currentPath();
			QDir::setCurrent(QString(projPath.c_str()));

			QDir projDir(QString(projPath.c_str()));

			std::string line; 
			std::pair<std::string,std::string> name_path;

			getline(file,line); 
			Volume* vol = VolumeUtils::readByteVolumeFromFile(line);
			if(!vol && autoGen)
				vol = new Volume(); 

			name_path = getFileNameAndPath(line, projDir);
			vol->setFileName(name_path.first);
			vol->setFilePath(name_path.second);
			//cout<<"volume: "<<vol->getFileName()<<endl;
			ret->setVolume(vol); 

			VolumeData<float>* fvol = VolumeUtils::readFloatVolumeFromFile(line);
			if(!fvol && autoGen){
				fvol = new VolumeData<float>(); 
			}
			ret->setFVolume(fvol);


			getline(file,line); 
			Volume* mask = VolumeUtils::readByteVolumeFromFile(line,false);
			if(!mask && autoGen){
				mask = new Volume(); 
				mask->resize(vol->getDimX(),vol->getDimY(),vol->getDimZ());
			}
			name_path = getFileNameAndPath(line, projDir);
			mask->setFileName(name_path.first);
			mask->setFilePath(name_path.second);
			ret->setMask(mask); 

			BoundaryMap* bm  = new BoundaryMap(); 
			bm->resize(mask->getDimX(),mask->getDimY(),mask->getDimZ()); 
			std::string distFname = getFilePath(line) + getFileNameNoExt(line) + ".vdist"; 
			FVolume* func = VolumeUtils::readFloatVolumeFromFile(distFname,false);
			if(func){
				bm->getFunction()->copy(func); 
				bm->getAuxFunction()->copy(func); 

				delete func;
				func = NULL; 
			}
			bm->setProject(ret);
			BoundaryMapUtils::updateCorners(bm); 
			BoundaryMapUtils::updateVolume(bm); 

			ret->setBoundaryMap(bm); 

			//cout<<"Mask: "<<mask->getFileName()<<endl;

			getline(file,line);
			TransferFunction* tf = TransferFunctionUtils::readFromFile(line);
			if(!tf && autoGen)
				tf = new TransferFunction();

			name_path = getFileNameAndPath(line, projDir);
			tf->setFileName(name_path.first);
			tf->setFilePath(name_path.second);
			ret->setTransFunc(tf); 

			getline(file,line);
			SegmentMap* smap = SegmentMapUtils::readFromFile(line);
			if(!smap && autoGen)
				smap = new SegmentMap();
				//smap->setFileName(line);

			name_path = getFileNameAndPath(line, projDir);
			smap->setFileName(name_path.first);
			smap->setFilePath(name_path.second);
			ret->setSegmentMap(smap); 

			UndoVolume* undo = new UndoVolume(); 
			undo->resize(vol->getDimX(),vol->getDimY(),vol->getDimZ()); 
			ret->setUndoVolume(undo);

			file.close();

			QDir::setCurrent(prevDir);

			return ret; 
		}
	}; 
}

#endif
