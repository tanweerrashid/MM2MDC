#include "VolumeUtils.h"
#include <stdio.h>
#include "mrc.h"
#include "Session.h"
#include "Common.h"
#include "TextureManager.h"
#include "ddsbase.h"
#include "BoundaryMap.h"

#include <minc2.h>
#include <omp.h>
#include <fstream>

using namespace vpt; 
using namespace std; 

// the input is max, then min, but really we want to pass in the minimum possible
//   value and then the maximum possible value. 
// If scale is true, then scale to within the range of vval.  Else just read
template <typename T> 
inline void readMRCInputByte(FILE* fin, T mmin, T mmax, int size, vval* data, bool scale=true){
	T* store = new T[size]; 

	fseek(fin,1024,SEEK_SET);
	fread(store,sizeof(T),size,fin);

	cout<<"read byte"<<endl;

	cout<<"min: "<<mmin<<" mmax: "<<mmax<<endl;


	if(scale){
/*
		mmin = 100000;
		mmax = -1000000;
        #pragma omp parallel for shared(mmax,mmin)
		for(int j=0;j<size;j++){
			mmin = min(store[j],mmin); 
			mmax = max(store[j],mmax); 
			//if(store[j]>mmax) store[j]=mmax; 
			//data[j] = vval_scale((store[j]-mmin)/(mmax-mmin));
		}
*/

		cout<<"min: "<<mmin<<" mmax: "<<mmax<<endl;

        #pragma omp parallel for shared(mmax,mmin)
		for(int j=0;j<size;j++){
			if(store[j]<mmin) store[j]=mmin; 
			if(store[j]>mmax) store[j]=mmax; 
			data[j] = vval_scale((store[j]-mmin)/(float)(mmax-mmin));
		}
	}
	else{
        #pragma omp parallel 
		for(int j=0;j<size;j++)
			data[j] = store[j];
	}

    delete [] store;
}


template <typename T> 
inline void readMRCInputFloat(FILE* fin, T mmin, T mmax, int size, float* data, bool scale=true){
	T* store = new T[size]; 

	fseek(fin,1024,SEEK_SET);
	fread(store,sizeof(T),size,fin);

	cout<<"read float"<<endl;

	cout<<"min: "<<mmin<<" mmax: "<<mmax<<endl;


	if(scale){
/*
		mmin = 100000;
		mmax = -1000000;

        #pragma omp parallel for shared(mmax,mmin)
		for(int j=0;j<size;j++){
			mmin = min(store[j],mmin); 
			mmax = max(store[j],mmax); 
			//if(store[j]>mmax) store[j]=mmax; 
			//data[j] = vval_scale((store[j]-mmin)/(mmax-mmin));
		}
*/

		cout<<"min: "<<mmin<<" mmax: "<<mmax<<endl;

        #pragma omp parallel for
		for(int j=0;j<size;j++){
			if(store[j]<mmin) store[j]=mmin; 
			if(store[j]>mmax) store[j]=mmax; 
			data[j] = (store[j]-mmin)/(float)(mmax-mmin);
		}
	}
	else{
        #pragma omp parallel for
		for(int j=0;j<size;j++)
			data[j] = store[j];
	}

    delete [] store;
}

template <typename T>
bool readMINCVolume(const string& filename, VolumeData<T>* vol, VolumeDataType type, bool scale){
	mihandle_t minc_volume;
	double voxel;
	int result;
	misize_t location[3];
 
	// Open MINC file
	cout << "Opening MINC file " << filename << endl;
	result = miopen_volume(filename.c_str(), MI2_OPEN_READ, &minc_volume);	
	if (result == MI_ERROR) cout << "Error opening input file: " << result << endl;
	
	// Set order of MINC x, y and z spaces
	char *dimorder[3];
	dimorder[0] = "xspace";
	dimorder[1] = "yspace";
	dimorder[2] = "zspace";
	result = miset_apparent_dimension_order_by_name(minc_volume, 3, dimorder); 
	if (result == MI_ERROR) cout << "Error in setting dimension order" << endl;

	// Get Volume dimension and step sizes
	midimhandle_t dimensions[3];
	misize_t sizes[3];
	result = miget_volume_dimensions(minc_volume, MI_DIMCLASS_SPATIAL, MI_DIMATTR_ALL, MI_DIMORDER_APPARENT, 3, dimensions);
	if (result == MI_ERROR) cout << "Error in retrieving volume dimensions" << endl;	
	result = miget_dimension_sizes(dimensions, 3, sizes);
	if (result == MI_ERROR) cout << "Error in retrieving dimension sizes" << endl;

	double stepSize[3];
	result = miget_dimension_separations(dimensions, MI_ORDER_APPARENT, 3, stepSize);
	if (result == MI_ERROR) cout << "Error in retrieving step sizes" << endl;
	
	// w -> x;  h -> y; d -> z
	int w = sizes[0], h = sizes[1], d = sizes[2];
	cout << "Volume dimensions:" << endl << "X: " << w << endl << "Y: " << h << endl << "Z: " << d << endl;

	double *v = new double[w * h * d];
	uint k = 0;	
	float mmax = -99999999999, mmin = 999999999; 
	

	cout << "Reading MINC voxels" << endl;
        
	#pragma omp parallel for
	for (int wi = 0; wi < w; wi = wi + 1) {
		for (int hi = 0; hi < h; hi = hi + 1) {
			for (int di = 0; di < d; di = di + 1) {				
				location[0] = wi;
				location[1] = hi;
				location[2] = di;                                
                                
				result = miget_voxel_value(minc_volume, location, 3, &voxel);
				if (result == MI_ERROR) cout << "Error reading voxels." << endl;
				
				v[k] = voxel;					
				
				if (v[k] >= mmax) mmax = v[k];
				if (v[k] <= mmin) mmin = v[k];
				
				//mmax = max(v[k], mmax);
				//mmin = min(v[k], mmin);				

				k = k + 1;
			}
		}
	}

	result = miclose_volume(minc_volume);
	if (result == MI_ERROR) cout << "Error in closing MINC volume" << endl;

	cout << "Finished reading voxels" << endl;
	cout << "mmin = " << mmin << endl;
	cout << "mmax = " << mmax << endl;
			
	vol->resize(w, h, d);
		
	#pragma omp parallel for
	for (int j = 0; j < w * h * d; j = j + 1) {
		if (type == VDT_BYTE) {
			//cout << "VDT_BYTE" << endl;
			(*vol)(j) = vval_scale((v[j] - mmin) / (mmax - mmin));			
		}
		else {
			//cout << "VDT_FLOAT" << endl;
			(*vol)(j) = (T)(v[j] - mmin) / (mmax - mmin);			
		}

	}

	/*ofstream file;
	file.open("C:\\Users\\Tanweer Rashid\\Desktop\\ZZZ.txt");

	for (int j = 0; j < w * h * d; j++) {
		file << vol->getData(j) << "\n";
	}
	file.close();*/

	MRCData *mrcdata = new MRCData();
	mrcdata->scalex = 1;
	mrcdata->scaley = 1;
	mrcdata->scalez = 1;

	mrcdata->offsetx = -45; 
	mrcdata->offsety = -51; 
	mrcdata->offsetz = -42; 

	mrcdata->stepx = stepSize[0];
	mrcdata->stepy = stepSize[1];
	mrcdata->stepz = stepSize[2];

	mrcdata->vmax = 1;
	mrcdata->vmin = 0;

	mrcdata->actualHeader = VolumeUtils::createMRCHeader(vol, MODE_float);

	vol->setMRCData(mrcdata);

	return true; 
}

template <typename T>
bool readMRCVolume(const string& filename, VolumeData<T>* vol, VolumeDataType type, bool scale){
		mrcH header; 
		FILE* in = fopen(filename.c_str(), "r");

		if(!in) return false; 
 

		if(fread(&header, sizeof(mrcH),1,in)==3){
			cerr<<"cannot read MRC header"<<endl;
			fclose(in);
			return false;
		}

		vol->resize(header.nc,header.nr,header.ns);
		T* data = vol->data();

		MRCData* mrcdata = new MRCData(); 
		mrcdata->stepx = header.xlen/header.mx; 
		mrcdata->stepy = header.ylen/header.my; 
		mrcdata->stepz = header.zlen/header.mz; 

		float maxdim = max(header.mx,max(header.my,header.mz));
		double maxstep = max(mrcdata->stepx,max(mrcdata->stepy,mrcdata->stepz));
		mrcdata->scalex = mrcdata->stepx/maxstep * header.mx/maxdim;
		mrcdata->scaley = mrcdata->stepy/maxstep * header.my/maxdim;
		mrcdata->scalez = mrcdata->stepz/maxstep * header.mz/maxdim;

		mrcdata->offsetx = header.mx/maxdim/2;
		mrcdata->offsety = header.my/maxdim/2;
		mrcdata->offsetz = header.mz/maxdim/2;

		mrcdata->vmax = header.amax; 
		mrcdata->vmin = header.amin; 
		mrcdata->actualHeader = header;

		vol->setMRCData(mrcdata);
		if(type==VDT_BYTE){
			switch(header.mode){
		    case MODE_char: readMRCInputByte(in,(unsigned char)0,(unsigned char)255,vol->getSize(),(vval*)data,scale); break; 
			case MODE_float: readMRCInputByte(in,(float) header.amin,(float)header.amax,vol->getSize(),(vval*)data,scale); break; 
			case MODE_short: readMRCInputByte(in,(short) header.amin,(short)header.amax,vol->getSize(),(vval*)data,scale); break;
			case 6: readMRCInputByte(in,(unsigned short) header.amin,(unsigned short)header.amax,vol->getSize(),(vval*)data,scale); break;
			//case 6: readMRCInputByte(in,(unsigned long) header.amin,(unsigned long)header.amax,vol->getSize(),(vval*)data,scale); break;
			default: cout<<"Unsupported MRC format "<<header.mode<<endl;
			}
		}
		else if(type==VDT_FLOAT){
			switch(header.mode){
		    case MODE_char: readMRCInputFloat(in,(unsigned char)header.amin,(unsigned char)header.amax,vol->getSize(),(float*)data,scale); break; 
			case MODE_float: readMRCInputFloat(in,(float) header.amin,(float)header.amax,vol->getSize(),(float*)data,scale); break; 
			case MODE_short: readMRCInputFloat(in,(short)header.amin,(short)header.amax,vol->getSize(),(float*)data,scale); break;
			case 6: readMRCInputFloat(in,(unsigned short) header.amin,(unsigned short)header.amax,vol->getSize(),(float*)data,scale); break;
			default: break;
		}

		cout<<"reading volume of size: "<<vol->getDimX()<<" "<<vol->getDimY()<<" "<<vol->getDimZ()<<endl;
	}

	fclose(in);
	return true; 
}

template <typename T>
bool readPVMVolume(const string& filename, VolumeData<T>* vol, VolumeDataType type, bool scale){
	uint w,h,d;
	float sx, sy, sz; 
	uint comp; 
	char* name = const_cast<char*>(filename.c_str()); 
	unsigned char* v = readPVMvolume(name,&w,&h,&d,&comp,&sx,&sy,&sz); // this function is from ddsbase

	if(!v) {
		cerr<<"failed in reading PVM file "<<filename<<endl;		
		return false; 
	}

	float mmax=0, mmin=255; 
	if(comp==2)
		mmin = 1<<16; 

	for(uint j=0;j<w*h*d*comp;j+=comp){
		int val = 0; 
		if(comp==2){
			val= v[j]<<16; 
			val= val | v[j+1]; 
		}
		else if(comp==1)
			val = v[j]; 

		mmax = max(val,mmax); 
		mmin = min(val,mmin);
	}

	vol->resize(w,h,d); 

    #pragma omp parallel for
	for(int j=0;j<(int)(w*h*d*comp);j+=((int)comp)){
		int val = 0; 
		if(comp==2){
			val= v[j]<<16; 
			val= val | v[j+1]; 
		}
		else if(comp==1)
			val = v[j]; 
		if(type==VDT_BYTE)
			(*vol)(j/comp) = vval_scale((val-mmin)/(mmax-mmin));
		else
			(*vol)(j/comp) = (T)(val-mmin)/(mmax-mmin);
	}

	MRCData* mrcdata = new MRCData(); 

	float maxdim = max(w,max(h,d));

	cout<<"size: "<<w<<" "<<h<<" "<<d<<endl;

	float maxstep = max(sx,max(sy,sz));

	mrcdata->scalex = sx/maxstep * w/maxdim; 
	mrcdata->scaley = sy/maxstep * h/maxdim; 
	mrcdata->scalez = sz/maxstep * d/maxdim; 

	float maxscale = max(mrcdata->scalex,max(mrcdata->scaley,mrcdata->scalez)); 
	mrcdata->scalex/=maxscale; 
	mrcdata->scaley/=maxscale; 
	mrcdata->scalez/=maxscale; 

	cout<<"scale: "<<sx<<" "<<sy<<" "<<sz<<" "<<w<<" "<<h<<" "<<d<<endl;
	cout<<"scalx: "<<mrcdata->scalex<<" "<<mrcdata->scaley<<" "<<mrcdata->scalez<<endl;

	mrcdata->offsetx = mrcdata->scalex/2;
	mrcdata->offsety = mrcdata->scaley/2;
	mrcdata->offsetz = mrcdata->scalez/2;
	mrcdata->vmax = 1; 
	mrcdata->vmin = 0; 

	switch(sizeof(T)){
       case 1: 
		   mrcdata->actualHeader = VolumeUtils::createMRCHeader(vol,MODE_char); 
		   break; 
	   case 4: 
		   mrcdata->actualHeader = VolumeUtils::createMRCHeader(vol,MODE_float); 
		   break; 
	   default: 
		   mrcdata->actualHeader = VolumeUtils::createMRCHeader(vol,MODE_char); 
		   break; 
	}

	vol->setMRCData(mrcdata);

	// the other side (v^3) uses malloc, must use free here.
	free(v);

	return true; 
}

template <typename T>
bool readIMGVolume(const string& filename, VolumeData<T>* vol, VolumeDataType type, bool scale){
	uint w,h,d;
	float sx, sy, sz; 
	uint comp; 
	char* name = const_cast<char*>(filename.c_str()); 
	//unsigned char* v = readPVMvolume(name,&w,&h,&d,&comp,&sx,&sy,&sz); // this function is from ddsbase

	FILE* file = fopen(name,"r");

	w = 256; 
	h = 256; 
	d = 104; 

	short* v = (new short[w*h*d]);
	fread(v,sizeof(short),w*h*d,file);

	sx=sy=sz=1;
	comp=1;

	if(!v) {
		cerr<<"failed in reading PVM file "<<filename<<endl;		
		return false; 
	}

	float mmax=0, mmin=255; 
	if(comp==2)
		mmin = 1<<16; 

	for(uint j=0;j<w*h*d*comp;j+=comp){
		int val = 0; 
		if(comp==2){
			val= v[j]<<16; 
			val= val | v[j+1]; 
		}
		else if(comp==1)
			val = v[j]; 

		mmax = max(val,mmax); 
		mmin = min(val,mmin);
	}

	vol->resize(w,h,d); 

    #pragma omp parallel for
	for(int j=0;j<(int)(w*h*d*comp);j+=((int)comp)){
		int val = 0; 
		if(comp==2){
			val= v[j]<<16; 
			val= val | v[j+1]; 
		}
		else if(comp==1)
			val = v[j]; 
		if(type==VDT_BYTE)
			(*vol)(j/comp) = vval_scale((val-mmin)/(mmax-mmin));
		else
			(*vol)(j/comp) = (T)(val-mmin)/(mmax-mmin);
	}

	MRCData* mrcdata = new MRCData(); 

	float maxdim = max(w,max(h,d));

	cout<<"size: "<<w<<" "<<h<<" "<<d<<endl;

	float maxstep = max(sx,max(sy,sz));

	mrcdata->scalex = sx/maxstep * w/maxdim; 
	mrcdata->scaley = sy/maxstep * h/maxdim; 
	mrcdata->scalez = sz/maxstep * d/maxdim; 

	float maxscale = max(mrcdata->scalex,max(mrcdata->scaley,mrcdata->scalez)); 
	mrcdata->scalex/=maxscale; 
	mrcdata->scaley/=maxscale; 
	mrcdata->scalez/=maxscale; 

	cout<<"scale: "<<sx<<" "<<sy<<" "<<sz<<" "<<w<<" "<<h<<" "<<d<<endl;
	cout<<"scalx: "<<mrcdata->scalex<<" "<<mrcdata->scaley<<" "<<mrcdata->scalez<<endl;

	mrcdata->offsetx = mrcdata->scalex/2;
	mrcdata->offsety = mrcdata->scaley/2;
	mrcdata->offsetz = mrcdata->scalez/2;
	mrcdata->vmax = 1; 
	mrcdata->vmin = 0; 

	switch(sizeof(T)){
       case 1: 
		   mrcdata->actualHeader = VolumeUtils::createMRCHeader(vol,MODE_char); 
		   break; 
	   case 4: 
		   mrcdata->actualHeader = VolumeUtils::createMRCHeader(vol,MODE_float); 
		   break; 
	   default: 
		   mrcdata->actualHeader = VolumeUtils::createMRCHeader(vol,MODE_char); 
		   break; 
	}

	vol->setMRCData(mrcdata);

	// the other side (v^3) uses malloc, must use free here.
	free(v);

	return true; 
}

template<typename T> 
bool readVolumePrime(const string& filename, VolumeData<T>& vol, VolumeDataType type, bool scale){
	std::string ext = getFileExtension(filename);
	if(ext=="mrc" || ext=="dat" || ext=="map" || ext=="vmask" || ext=="vdist")
		return readMRCVolume(filename,&vol,type,scale);
	else if(ext=="pvm")
		return readPVMVolume(filename,&vol,type,scale);
	else if(ext=="img")
		return readIMGVolume(filename,&vol,type,scale);
	else if (ext == "mnc") 
		return readMINCVolume(filename, &vol, type, scale);
	return false; 
}


bool  VolumeUtils::readVolume(const string& filename, Volume& vol, bool scale){
	return readVolumePrime(filename,vol,VDT_BYTE,scale);
}

bool  VolumeUtils::readVolume(const string& filename, FVolume& vol, bool scale){
	return readVolumePrime(filename,vol,VDT_FLOAT,scale); 
}


void VolumeUtils::writeToFile(const std::string& filename, Volume* vol, VolumeIOType type){
	std::fstream file(filename.c_str(),std::ios::out | std::ios::binary); 

	if(type==MRC_TYPE){
		mrcH mrcHeader;
		mrcHeader = createMRCHeader(vol,MODE_char);
		file.write((char*)&mrcHeader, sizeof(mrcH)); 
		file.write((char*)vol->data(), sizeof(vval)*vol->getSize()); 
	}

	file.close();
}

void VolumeUtils::writeToFile(const std::string& filename, FVolume* vol, VolumeIOType type){
	std::fstream file(filename.c_str(),std::ios::out | std::ios::binary); 

	if(type==MRC_TYPE){
		mrcH mrcHeader = createMRCHeader(vol,MODE_float);
		file.write((char*)&mrcHeader, sizeof(mrcH)); 
		file.write((char*)vol->data(), sizeof(float)*vol->getSize()); 
	}

	file.close();
}

void VolumeUtils::writeSegmentedVolToMRCFile(const std::string& filename, FVolume* vol, Volume* mask, SegmentMap* segmap, VolumeIOType type){
	bool show[NUM_SEGMENTS]; 
	int nsegs = segmap->getNumSegments(); 
	for(int j=0;j<nsegs;j++)
		show[j] = segmap->getEntry(j)->getState()!=vpt::HIDE_SEGMENT; 

	FVolume* output = new FVolume(); 
	output->resize(vol->getDimX(),vol->getDimY(),vol->getDimZ()); 

	int dz = (int) output->getDimZ(); 
	int dx = (int) output->getDimX(); 
	int dy = (int) output->getDimY(); 

    #pragma omp parallel for
	for(int k=0;k<dz;k++){
		for(int j=0;j<dy;j++)
			for(int i=0;i<dx;i++){
				if(show[(*mask)(i,j,k)])
					(*output)(i,j,k) = (*vol)(i,j,k); 
				else
					(*output)(i,j,k) = vol->getMRCData()->vmin; 
			}
	}

	writeToFile(filename, output, MRC_TYPE); 

	delete output; 
}

// I need to think about a better way to handle this later
osg::Image* VolumeUtils::convertToOSGImage(Volume& vol){
	osg::Image* ret = new osg::Image(); 
	ret->setImage(vol.getDimX(),vol.getDimY(),vol.getDimZ(),VV_ITEXTURE_FORMAT,VV_PIXEL_FORMAT,VV_DATA_FORMAT,vol.data(),osg::Image::NO_DELETE);
	return ret; 
}

osg::Image* VolumeUtils::convertToOSGImage(FVolume& vol){
	osg::Image* ret = new osg::Image(); 
	ret->setImage(vol.getDimX(),vol.getDimY(),vol.getDimZ(),VV_ITEXTURE_FORMAT,VV_PIXEL_FORMAT,GL_FLOAT,(unsigned char*)vol.data(),osg::Image::NO_DELETE);
	return ret; 
}


void VolumeUtils::updateTexture(Volume* v, uint xs, uint ys, uint zs, uint xdim, uint ydim, uint zdim, vval* ndata){
	if(TextureManager::inst()->hasVolumeTex3D(v)){
		osg::ref_ptr<osg::Texture3D> tex = TextureManager::inst()->getOrCreateVolumeTex3D(v); 
		if(tex->getSubloadCallback()==NULL){
			GLenum iformat = tex->getImage()->getInternalTextureFormat(); 
			GLenum pformat = tex->getImage()->getPixelFormat(); 
			GLenum type = tex->getImage()->getDataType(); 

			tex->setSubloadCallback( new VolumeSubloadCallback<vval>(iformat,pformat,type)); 
		}

		osg::ref_ptr<VolumeSubloadCallback<vval> > cb = dynamic_cast<VolumeSubloadCallback<vval>*>(tex->getSubloadCallback());
		cb->initSubload(xs,ys,zs,xdim,ydim,zdim,ndata); 
		//cout<<"finished "<<v<<endl;
	}
}

osg::Vec3 VolumeUtils::smoothGradient(Volume* v, const osg::Vec3& pos){
	//cout<<"p: "<<pos.x()<<" "<<pos.y()<<" "<<pos.z()<<endl;

	double qx = 1./(v->getDimX()-1);
	double qy = 1./(v->getDimY()-1);
	double qz = 1./(v->getDimZ()-1);

//	cout<<"v: "<<v->sample(pos.x(),pos.y(),pos.z())<<endl;

	double mask[125]={
		1,1,1,1,1,
		1,1,1,1,1,
		1,1,1,1,1,
		1,1,1,1,1,
		1,1,1,1,1,

		1,1,1,1,1,
		1,3,3,3,1,
		1,3,3,3,1,
		1,3,3,3,1,
		1,1,1,1,1,

		1,1,1,1,1,
		1,3,3,3,1,
		1,3,5,3,1,
		1,3,3,3,1,
		1,1,1,1,1,

		1,1,1,1,1,
		1,3,3,3,1,
		1,3,3,3,1,
		1,3,3,3,1,
		1,1,1,1,1,

		1,1,1,1,1,
		1,1,1,1,1,
		1,1,1,1,1,
		1,1,1,1,1,
		1,1,1,1,1,
	}; 

	double dx = 
		v->sampleSmooth(pos.x()+qx,pos.y(),pos.z(),5,mask) - 
		v->sampleSmooth(pos.x()-qx,pos.y(),pos.z(),5,mask); 

	double dy = 
		v->sampleSmooth(pos.x(),pos.y()+qy,pos.z(),5,mask) - 
		v->sampleSmooth(pos.x(),pos.y()-qy,pos.z(),5,mask); 

	double dz = 
		v->sampleSmooth(pos.x(),pos.y(),pos.z()+qz,5,mask) - 
		v->sampleSmooth(pos.x(),pos.y(),pos.z()-qz,5,mask); 

	return osg::Vec3(dx,dy,dz);
}

void VolumeUtils::growMask(Volume* mask, Volume* vol, TransferFunction* tf, vval segment, vval baseSegment){
	int dimx = (int) mask->getDimX(); 
	int dimy = (int) mask->getDimY(); 
	int dimz = (int) mask->getDimZ(); 

	int dx[6]={-1,0,0,0,0,1}; 
	int dy[6]={0,0,1,0,-1,0}; 
	int dz[6]={0,1,0,-1,0,0}; 

	Volume* tmpMask = new Volume(dimx,dimy,dimz); 

	float visible = .3f; 

    #pragma omp parallel for
	for(int k=0;k<dimz;k++){
		for(int j=0;j<dimy;j++){
			for(int i=0;i<dimx;i++){
				if((*mask)(i,j,k)==baseSegment){
					for(int p=0;p<6;p++){
						int ni=dx[p]+i;
						int nj=dy[p]+j;
						int nk=dz[p]+k;

						if(ni>=0 && ni<dimx && 
							nj>=0 && nj<dimy && 
							nk>=0 && nk<dimz &&
							(*mask)(ni,nj,nk)==segment && 
							tf->sampleCoarse((*vol)(i,j,k))>visible){
							(*tmpMask)(i,j,k)=segment; 
							break;
						}
					}
				}
			}
		}
	}

    #pragma omp parallel for
	for(int k=0;k<dimz;k++)
		for(int j=0;j<dimy;j++)
			for(int i=0;i<dimx;i++)
				if((*tmpMask)(i,j,k)==segment)
					(*mask)(i,j,k)=segment;

	delete tmpMask; 
}

float vuTrilinear(float* vals, float px, float py, float pz){
	return vals[0]*(1.f-px)*(1.f-py)*(1.f-pz) + 
		vals[1]*(px)*(1.f-py)*(1.f-pz) + 
		vals[2]*(1.f-px)*(py)*(1.f-pz) + 
		vals[3]*(px)*(py)*(1.f-pz) + 
		vals[4]*(1.f-px)*(1.f-py)*(pz) + 
		vals[5]*(px)*(1.f-py)*(pz) + 
		vals[6]*(1.f-px)*(py)*(pz) + 
		vals[7]*(px)*(py)*(pz); 
}


std::pair<float,vval> vuGetValAndMask(FVolume* vol, Volume* mask, int qi, int qj, int qk, float fi, float fj, float fk){
	float vals[8];
	vals[0] = (*vol)(qi,qj,qk); 
	vals[1] = (*vol)(qi+1,qj,qk); 
	vals[2] = (*vol)(qi,qj+1,qk); 
	vals[3] = (*vol)(qi+1,qj+1,qk); 
	vals[4] = (*vol)(qi,qj,qk+1); 
	vals[5] = (*vol)(qi+1,qj,qk+1); 
	vals[6] = (*vol)(qi,qj+1,qk+1); 
	vals[7] = (*vol)(qi+1,qj+1,qk+1); 

	int maskvals[8]; 
	maskvals[0] = (*mask)(qi,qj,qk); 
	maskvals[1] = (*mask)(qi+1,qj,qk); 
	maskvals[2] = (*mask)(qi,qj+1,qk); 
	maskvals[3] = (*mask)(qi+1,qj+1,qk); 
	maskvals[4] = (*mask)(qi,qj,qk+1); 
	maskvals[5] = (*mask)(qi+1,qj,qk+1); 
	maskvals[6] = (*mask)(qi,qj+1,qk+1); 
	maskvals[7] = (*mask)(qi+1,qj+1,qk+1); 

	int dmap[8]; 
	int numMats = 0; 
	float func[8][8]; 

	for(int j=0;j<8;j++){
		int mat = maskvals[j]; 
		if(mat>=0){
			for(int k=0;k<8;k++){
				if(mat==maskvals[k]){
					func[numMats][k] = vals[k];
					maskvals[k] = -1; 
				}
				else{
					//func[numMats][k] = -vals[k];
					func[numMats][k] = 0;
				}
			}
			dmap[numMats] = mat; 
			numMats++;
		}
	}

	if(numMats==1)
		return std::make_pair(1.5f,/*vuTrilinear(func[0], fi, fj, fk),*/ (vval)dmap[0]);

	float firstBestVal = -1000000000; 
	int firstBestMat = -1; 
	float secondBestVal = -1000000000; 
	int secondBestMat = -1; 
 
	for(int j=0;j<numMats;j++){
		int mat = dmap[j]; 
		float tval = vuTrilinear(func[j],fi,fj,fk);
		if(tval > firstBestVal){
			firstBestVal = tval;
			firstBestMat = mat;
		}
	}

	return std::make_pair(abs(firstBestVal),firstBestMat);
}

float vuGetVal(FVolume* vol, Volume* mask, vval mat, int qi, int qj, int qk, float fi, float fj, float fk){
	float vals[8];
	vals[0] = (*vol)(qi,qj,qk); 
	vals[1] = (*vol)(qi+1,qj,qk); 
	vals[2] = (*vol)(qi,qj+1,qk); 
	vals[3] = (*vol)(qi+1,qj+1,qk); 
	vals[4] = (*vol)(qi,qj,qk+1); 
	vals[5] = (*vol)(qi+1,qj,qk+1); 
	vals[6] = (*vol)(qi,qj+1,qk+1); 
	vals[7] = (*vol)(qi+1,qj+1,qk+1); 

	int maskvals[8]; 
	maskvals[0] = (*mask)(qi,qj,qk); 
	maskvals[1] = (*mask)(qi+1,qj,qk); 
	maskvals[2] = (*mask)(qi,qj+1,qk); 
	maskvals[3] = (*mask)(qi+1,qj+1,qk); 
	maskvals[4] = (*mask)(qi,qj,qk+1); 
	maskvals[5] = (*mask)(qi+1,qj,qk+1); 
	maskvals[6] = (*mask)(qi,qj+1,qk+1); 
	maskvals[7] = (*mask)(qi+1,qj+1,qk+1); 

//	int dmap[8]; 
	int numMats = 0; 
	float func[8]; 

	for(int k=0;k<8;k++){
		if(mat==maskvals[k])
			func[k] = vals[k];
		else
			func[k] = 0;
	}

	return vuTrilinear(func,fi,fj,fk);

	//return firstBestVal;
}

void VolumeUtils::buildPerMaterialFunction(FVolume* vol, Volume* mask, vval mat, FVolume* nvol, int ntx, int nty, int ntz){
	int dimx = vol->getDimX(); 
	int dimy = vol->getDimY(); 
	int dimz = vol->getDimZ(); 

	nvol->resize(ntx,nty,ntz);

	for(int k=0;k<ntz;k++){
		float pz = (k/(float)(ntz-1)) * (dimz-1);
		int qz = (int) pz; 
		if(qz==dimz-1) { pz=1.f; qz--; }
		else pz-=qz;

		for(int j=0;j<nty;j++){
			float py = (j/(float)(nty-1)) * (dimy-1);
			int qy = (int) py; 
			if(qy==dimy-1) { py=1.f; qy--; }
			else py-=qy;

			for(int i=0;i<ntx;i++){
				float px = (i/(float)(ntx-1)) * (dimx-1);
				int qx = (int) px; 
				if(qx==dimx-1) { px=1.f; qx--; }
				else px-=qx;

				(*nvol)(i,j,k) = vuGetVal(vol,mask,mat,qx,qy,qz,px,py,pz);
				//(*nvol)(i,j,k) = p; 
			}
		}
	}
}

void VolumeUtils::buildMaterialFunctions(SegmentMap* smap, BoundaryMap* bm, FVolume* vol, Volume* mask, int ntx, int nty, int ntz){
	int dimx = vol->getDimX(); 
	int dimy = vol->getDimY(); 
	int dimz = vol->getDimZ(); 

	for(int j=0;j<smap->getNumSegments();j++){
		SegmentMapEntry* ent = smap->getEntry(j); 
		vval mat = (vval) ent->getSegmentCode(); 
		FVolume* fvol = bm->getFunction(j); 
		if(fvol==NULL)
			fvol = new FVolume(); 

		bm->setFunction(j,fvol);
		buildPerMaterialFunction(vol,mask,mat,fvol,ntx,nty,ntz);
	}
}

void VolumeUtils::subsample(FVolume* vol, Volume* mask, FVolume* nvol, Volume* nmask, int ntx, int nty, int ntz){
	int dimx = vol->getDimX(); 
	int dimy = vol->getDimY(); 
	int dimz = vol->getDimZ(); 

	nvol->resize(ntx,nty,ntz);
	nmask->resize(ntx,nty,ntz);

	float fmin = 10000000; 
	float fmax = -10000000; 

	for(int k=0;k<ntz;k++){
		float pz = (k/(float)(ntz-1)) * (dimz-1);
		int qz = (int) pz; 
		if(qz==dimz-1) { pz=1.f; qz--; }
		else pz-=qz;

		for(int j=0;j<nty;j++){
			float py = (j/(float)(nty-1)) * (dimy-1);
			int qy = (int) py; 
			if(qy==dimy-1) { py=1.f; qy--; }
			else py-=qy;

			for(int i=0;i<ntx;i++){
				float px = (i/(float)(ntx-1)) * (dimx-1);
				int qx = (int) px; 
				if(qx==dimx-1) { px=1.f; qx--; }
				else px-=qx;

				std::pair<float,vval> p = vuGetValAndMask(vol,mask,qx,qy,qz,px,py,pz);
				(*nvol)(i,j,k) = p.first; 
				(*nmask)(i,j,k) = (vval) p.second; 

				fmin = min(fmin,p.first); 
				fmax = max(fmax,p.first); 
			}
		}
	}

	float length = fmax - fmin;  
	for(int j=0;j<nvol->getSize();j++)
		(*nvol)(j) = ((*nvol)(j)-fmin) / length;
}
