//#include <QString>
//#include <QFileDialog>
#include "Project.h"
#include "ProjectUtils.h"
#include "DualContouring.h"
#include "NormalMapUtils.h"
#include "Common.h"
#include "eigen.h"

#include <iostream>
#include <minc2.h>
#include <minc.h>
#include <omp.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkFillHolesFilter.h>

using namespace vpt; 

/**
 *  Saves the mesh as a VTK PolyData file 
 */
void saveMeshAsPolyData(CorrectOctree* _dctree, SegmentMap* smap, std::string outputFile, std::string mincFile, bool inWorldCoordinates) {
    std::fstream file(outputFile.c_str(), std::fstream::out);
    
    MeshGeometry* geomData = _dctree->m_Geom; 
    MeshGeometry* geomColors = _dctree->m_GeomColors;

    int numberOfVertices = geomData->getNumVerts(); 
    int numberOfCells = geomData->getNumTris(); 
    
    cout << "Number of vertices: " << numberOfVertices << endl;
    cout << "Number of cells: " << numberOfCells << endl;
    
    file << "# vtk DataFile Version 3.0\n";
    file << "Zhang Hybrid DC Output\n"; 
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    file << "POINTS " << numberOfVertices << " float\n";
    
    if (inWorldCoordinates) { // Write the vertices in world coordinates
        mihandle_t minc_volume;
        int result;

        // Open MINC file
        cout << "Opening MINC file " << mincFile << " for coordinate conversion" << endl;
        result = miopen_volume(mincFile.c_str(), MI2_OPEN_READ, &minc_volume);	
        if (result == MI_ERROR) cout << "Error opening input file: " << result << endl;

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

        double start[3];
        result = miget_dimension_starts(dimensions, MI_ORDER_FILE, 3, start);
        if (result == MI_ERROR) cout << "Error in retrieving volume origin" << endl;
                
        for (int i = 0; i < numberOfVertices; i++) {
            float *vert = geomData->getVert(i);
            
            double world[3];
            double voxel[3];
            voxel[0] = vert[0];
            voxel[1] = vert[1];
            voxel[2] = vert[2];

            result = miconvert_voxel_to_world(minc_volume, voxel, world);
            if (result == MI_ERROR) cout << "Error converting voxel to world coordinates: " << result << endl;
            
            // This extra translation is needed to ensure that the converted mesh
            // will be properly aligned with a .nii version of the MINC file in ParaView
            // For whatever reasons, a simple conversion from voxel to world coordinates is 
            // always offset by an amount of <minc volume's start values>
            double world_translated[3];

            start[0] = 0.0; start[1] = 0.0; start[2] = 0.0;
            world_translated[0] = world[0] - start[0];
            world_translated[1] = world[1] - start[1];
            world_translated[2] = world[2] - start[2];
            
            file << world_translated[0] << " " << world_translated[1] << " " << world_translated[2] << "\n";
        }
        file << "\n";
    }
    else { // Otherwise write the vertices in voxel coordinates. 
        for (int i = 0; i < numberOfVertices; i++) {
            float *vert = geomData->getVert(i);
            file << vert[0] << " " << vert[1] << " " << vert[2] << "\n";
        }
        file << "\n";
    }
//    file << "VERTICES " << numberOfVertices << " " << (numberOfVertices * 2) << "\n";
//    for (int i = 0; i < numberOfVertices; i++) {
//        file << "1 " << i << "\n";
//    }
//    file << "\n";
        
    file << "POLYGONS " << numberOfCells << " " << (numberOfCells * 4) << "\n";
    for (int i = 0; i < numberOfCells; i++) {
        unsigned int *cell = geomData->getTri(i);
        file << "3 " << cell[0] << " " << cell[1] << " " << cell[2] << "\n";
    }
    file << "\n";

    // Array 0. Used for coloring the mesh in Paraview.
//    file << "CELL_DATA " << numberOfCells << "\n";
//    file << "SCALARS paraview_color int 1"; 
//    file << "\nLOOKUP_TABLE paraview_color_table\n";
    file << "CELL_DATA " << numberOfCells << "\n";
    file << "SCALARS scalars float\n"; 
    file << "\nLOOKUP_TABLE default\n";
    for (int i = 0; i < numberOfCells; i++) {
        unsigned int *color = geomColors->getTri(i);
        file << color[0] << "\n";
    }
    file << "\n";
    
    // Array 1. Contains information regarding the materials of the faces/cells. The information is stored as an array of 2 elements [m0, m1]. 
    // For example, [0, 3] means that the cell lies between materials 0 and 3. 
    // For example, [4, 5] means that the cell lies between materials 4 and 5. 
    // 0 represents background. 1 is not used by the program. Materials are labeled starting from 2 and onwards.     
//    file << "SCALARS material_id int 2"; 
//    file << "\nLOOKUP_TABLE material_id_table\n";
    file << "FIELD Materials 1\n";
    file << "MaterialIndices 2 " << numberOfCells << " int\n";
    for (int i = 0; i < numberOfCells; i++) {
        unsigned int *color = geomColors->getTri(i);
        file << color[1] << " " << color[2] << " \n";        
    }
    file << "\n";

    file.close();
}

void saveMesh(CorrectOctree* _dctree, SegmentMap* smap, std::string outputFile){

	////////////////////////////////////////////
	// This is for writing out the object with 
	// cells being quads. 
	////////////////////////////////////////////
	//SegmentMap* smap = proj->getSegmentMap(); 	
	/*
	MeshGeometry* geomData = _dctree->m_Geom; 
//	VertexData* normalData = _dctree->m_GeomNormals;
	MeshGeometry* geomColors = _dctree->m_GeomColors;

	int nverts = geomData->getNumVerts(); 
	int nquads = geomData->getNumQuads(); 
        
	unsigned int* inds = geomData->getInds(); 
	float* verts = geomData->getVerts();

	std::fstream file(outputFile.c_str(), std::fstream::out);

	int* indMap = new int[nverts]; 

	for(int j=0;j<nverts;j++)
		indMap[j] = -1; 

	for(int j=0;j<nquads;j++){
		for(int k=0;k<4;k++){
			int i =inds[j*4+k]; 
			if(indMap[i]<0)
				indMap[i]=0;
		}
	}

	cout<<"nverts: "<<nverts<<endl;
	cout<<"nfaces: "<<nquads<<endl;

	int count=0;

	for(int j=0;j<nverts;j++){
		if(indMap[j]>=0){
			file<<"v "<<verts[j*3]<<" "<<verts[j*3+1]<<" "<<verts[j*3+2]<<endl;
			indMap[j]=count++;
			indMap[j]++;
		}
	}

	file<<endl;
	file<<"g all"<<endl;
	for(int j=0;j<nquads;j++){
		file << "f " << indMap[inds[j*4]] << " " <<indMap[inds[j*4+1]] << " " << indMap[inds[j*4+2]] << " " <<indMap[inds[j*4+3]] << endl;
	}

	file<<endl;

	file<<"bind m face"<<endl;

	float* colors = geomColors->getVerts(); 
        
        //cout << "Size of *colors: " << (sizeof(colors)/sizeof(float)) << endl;
        
	for(int j=0;j<nquads;j++){
		file<<"m ";
		int color = (int) colors[j*3]; 
		int rc=0;

		for(int j=0;j<32;j++){
			if(color & (1<<j)){
				rc = j; 
				if(j==0 || !smap->getEntryFromCode(rc)->getState()==vpt::SHOW_SEGMENT)
					rc = -1; 
				file<<rc<<" "; 
			}
		}
		file<<endl;
	}

	file.close();

	delete [] indMap; 
	*/









	////////////////////////////////////////////
	// This is for writing out the object with 
	// cells being triangles. 
	////////////////////////////////////////////

	//SegmentMap* smap = proj->getSegmentMap(); 	
	
	MeshGeometry* geomData = _dctree->m_Geom; 
//	VertexData* normalData = _dctree->m_GeomNormals;
	MeshGeometry* geomColors = _dctree->m_GeomColors;

	int nverts = geomData->getNumVerts(); 
	int ntri = geomData->getNumTris(); 
        
	unsigned int* inds = geomData->getInds(); 
	float* verts = geomData->getVerts();

	std::fstream file(outputFile.c_str(), std::fstream::out);

	int* indMap = new int[nverts]; 

	for(int j=0;j<nverts;j++)
		indMap[j] = -1; 

	for(int j=0;j<ntri;j++){
		for(int k=0;k<3;k++){
			int i =inds[j*3+k]; 
			if(indMap[i]<0)
				indMap[i]=0;
		}
	}


        cout<<"nverts: "<<nverts<<endl;
	cout<<"nfaces: "<<ntri<<endl;

	int count=0;

	for(int j=0;j<nverts;j++){
		if(indMap[j]>=0){
			file<<"v "<<verts[j*3]<<" "<<verts[j*3+1]<<" "<<verts[j*3+2]<<endl;
			indMap[j]=count++;
			indMap[j]++;
		}
	
        }
        
	file<<endl;
	file<<"g all"<<endl;
	for(int j=0;j<ntri;j++){
		file << "f " << indMap[inds[j*3]] << " " <<indMap[inds[j*3+1]] << " " << indMap[inds[j*3+2]] << endl;// " " <<indMap[inds[j*4+3]] << endl;
	}

	file<<endl;

	file<<"bind m face"<<endl;

	unsigned int* colors = geomColors->getTri(0);
        
        //cout << "Size of *colors: " << (sizeof(colors)/sizeof(float)) << endl;
        
	for(int j=0;j<ntri;j++){
		file<<"m ";
		int color = (int) colors[j*3]; 
		int rc=0;

		for(int j=0;j<32;j++){
			if(color & (1<<j)){
				rc = j; 
				if(j==0)// || !smap->getEntryFromCode(rc)->getState()==vpt::SHOW_SEGMENT)
					rc = -1; 
				file<<rc<<" "; 
			}
		}
		file<<endl;
	}

	file.close();

	delete [] indMap; 

}

/**
 * Initializes the DC process. Powei's function. 
 * 
 */
void buildDualContour(Volume* mask, FVolume* vol, MRCData* mrcdata, SegmentMap* smap, std::string filePath, std::string outputFile, bool inWorldCoordinates, double segmentValues[], int numSegments, BoundaryMap* bmap=NULL) {
    CorrectOctree* _dctree; 

    if (vol && mask){
        int twoTo = 0; 
        int maxDim = max(vol->getDimX(),max(vol->getDimY(),vol->getDimZ())); 

        cout << "DimX: " << vol->getDimX() << "\tDimY: " << vol->getDimY() << "\tDimZ: " << vol->getDimZ() << endl;
        
        // Bitwise left shift
        while((1<<twoTo)+1<maxDim)
                twoTo++; 

        cout << "twoTo: " << twoTo << endl;
        int actualSize = (1<<twoTo) + 1; 

        cout << "ActualSize: " << actualSize << endl;

        //cout<<"actualSize: "<<actualSize<<" "<<maxDim<<endl;
                
        // nvol is an array of double type that probably contains the volume data in linear order. 
        double* nvol = new double[actualSize * actualSize * actualSize]; 
        
        FVolume* tdvol = new FVolume();
        Volume* tmask = new Volume();
        tdvol->resize(actualSize, actualSize, actualSize);
        tmask->resize(actualSize, actualSize, actualSize);

        //SegmentMap* smap = proj->getSegmentMap(); 

        double nsx = (actualSize) / ((double)vol->getDimX()) * mrcdata->scalex; 
        double nsy = (actualSize) / ((double)vol->getDimY()) * mrcdata->scaley;
        double nsz = (actualSize) / ((double)vol->getDimZ()) * mrcdata->scalez;

        cout << "nsx: " << nsx << "\tnsy: " << nsy << "\tnsz: " << nsz << endl;

        //#pragma omp parallel for
        for (int k = 0; k < actualSize; k++) {
            for (int j = 0; j < actualSize; j++) {
                for (int i = 0; i < actualSize; i++) {
                    if (vol->valid(i, j, k)) { // Valid voxel coordinates or not. 
                        vval ent = (*mask)(i, j, k);
                        if (ent != SegmentMap::SELECTION_SEGMENT && ent != SegmentMap::BASE_SEGMENT && smap->getEntryFromCode(ent)->getState() == vpt::SHOW_SEGMENT) {
                            nvol[i + j * actualSize + k * actualSize * actualSize] = (*vol)(i, j, k); 
                            (*tmask)(i, j, k) = (*mask)(i, j, k);
                        }
                        else {
                            nvol[i + j * actualSize + k * actualSize * actualSize] = -(*vol)(i, j, k); 
                            (*tmask)(i, j, k) = SegmentMap::BASE_SEGMENT;
                        }
                        (*tdvol)(i, j, k) = (*vol)(i, j, k);
                    }
                    else {
                        nvol[i + j * actualSize + k * actualSize * actualSize] = 0;
                        (*tdvol)(i, j, k) = 0;
                    }
                }
            }
        }

        vector<osg::Vec4> colorTable; 

        for(int j=0;j<smap->getNumSegments();j++)
                colorTable.push_back(smap->getEntry(j)->getColor());

        //if(_dctree==NULL)
                _dctree = new CorrectOctree();
        //else{
        //	delete _dctree; 
        //	_dctree = new CorrectOctree();
        //}

        CorrectOctree* co = _dctree; 
                
        // *******************************************
        // cmin and cmax controls the amount of scaling and translation. 
        // Basically, these two parameters define the size of tdvol and tdmask, 
        // cmin determines the starting corner.         
        // *******************************************
        
        // **************************************************************************************************************
        // For computing the mesh in Voxel coordinates: 
        //double cmin[3] = {-0.5, -0.5, -0.5};
        //double cmax[3] = {0.5, 0.5, 0.5};
        
        double cmin[3] = {0, 0, 0};
        //double cmax[3] = {actualSize * mrcdata->stepx, actualSize * mrcdata->stepy, actualSize * mrcdata->stepz}; // This results in zero non-manifold edges
        double cmax[3] = {actualSize, actualSize, actualSize}; 
        //double cmax[3] = {nsx, nsy, nsz};
        //double cmax[3] = {mask->getDimX(), mask->getDimY(), mask->getDimZ()}; // This results in a few non-manifold edges. 
        //double cmax[3] = {10, 10, 20};
        // **************************************************************************************************************
        
        // **************************************************************************************************************
        // For computing the mesh in World coordinates: 
        //double cmin[3] = {mrcdata->offsetx, mrcdata->offsety, mrcdata->offsetz};
        //double cmax[3] = {actualSize * mrcdata->stepx + mrcdata->offsetx, actualSize * mrcdata->stepy + mrcdata->offsety, actualSize * mrcdata->stepz + mrcdata->offsetz};
        //double cmax[3] = {mask->getDimX() * mrcdata->stepx + mrcdata->offsetx, mask->getDimY() * mrcdata->stepy + mrcdata->offsety, mask->getDimZ() * mrcdata->stepz + mrcdata->offsetz};
        // **************************************************************************************************************
        
        
        double fx = ((double)vol->getDimX())/actualSize/2+cmin[0];
        double fy = ((double)vol->getDimY())/actualSize/2+cmin[1]; 
        double fz = ((double)vol->getDimZ())/actualSize/2+cmin[2]; 
        
        cout << "fx = " << fx << "\tfy = " << fy << "\tfz = " << fz << endl;        
        cout << "twoTo: " << twoTo << endl;
        cout << "actualSize: " << actualSize << endl;
        
        MeshGeometry* geomData = co->m_Geom; 
        
        //CorrectOctree::initialize(double coarseIso, uint coarseLevel, double* coarseMin, double* coarseMax, double* coarseFuncVals, int segment, Volume* mask, FVolume* dvol, BoundaryMap* bmap){
        
        
//        float *a = vol->data();
//        double *volData = new double[vol->getDimX() * vol->getDimY() * vol->getDimZ()];
//        for (int i = 0; i < vol->getDimX() * vol->getDimY() * vol->getDimZ(); i++) {
//            volData[i] = a[i];
//        }        
//        co->initialize(0, twoTo, cmin, cmax, volData, 0, mask, vol, bmap);
        
        co->initialize(0, twoTo, cmin, cmax, nvol, numSegments, tmask, tdvol, bmap);
        co->build(); 
        co->contour();

        cout << "Saving mesh..." << endl;        
        
        //saveMesh(co, smap, outputFile);
        saveMeshAsPolyData(co, smap, outputFile, filePath, inWorldCoordinates);        
        
        delete co; 
        delete [] nvol; 
        delete tdvol; 
        delete tmask;
    }	
}

/**
 * Read volume data. 
 * 
 * @param filename - Location of the minc volume
 * @param segmentValues - A list of volume labels for which surfaces will be generated. 
 * @param numSegments - Number of volume labels. 
 * 
 * @return 
 */
Project* readVolumeData(std::string filename, double segmentValues[], int numSegments) {
    // Initialize Vol and FVol variables. 
    Volume* vvol = new Volume();
    VolumeData<float>* fvol = new FVolume();
    Volume* mvol = new Volume(); 
    
    // Read MINC volume. 
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

    double start[3];
    result = miget_dimension_starts(dimensions, MI_ORDER_FILE, 3, start);
    if (result == MI_ERROR) cout << "Error in retrieving volume origin" << endl;
   
    // w is along x axis;  h is along y axis; d is along z axis. 
    int w = sizes[0], h = sizes[1], d = sizes[2];
    cout << "Volume dimensions: [" << w << ", " << h << ", " << d << "]" << endl;
    cout << "Volume stepsize: [" << stepSize[0] << ", " << stepSize[1] << ", " << stepSize[2] << "]" << endl;
    cout << "Volume start: [" <<  start[0] << ", " << start[1] << ", " << start[2] << "]" << endl;
    
    vvol->resize(w, h, d);
    fvol->resize(w, h, d);
    mvol->resize(w, h, d);
    
    double *v = new double[w * h * d];
    uint k = 0;	
    float mmax = -FLT_MAX, mmin = FLT_MAX; 

    cout << "Reading MINC voxels" << endl;

    #pragma omp parallel for
    for (int wi = 0; wi < w; wi = wi + 1) {
        for (int hi = 0; hi < h; hi = hi + 1) {
            for (int di = 0; di < d; di = di + 1) {                
                location[0] = wi;
                location[1] = hi;
                location[2] = di;                                

                // Get voxel value at specified voxel location. 
                result = miget_voxel_value(minc_volume, location, 3, &voxel);
                if (result == MI_ERROR) cout << "Error reading voxels." << endl;

                v[k] = voxel; // Populate the v[] array. 
                
                // Begin mask modification
                for (int q = 0; q < numSegments; q++) {
                    if (voxel == segmentValues[q]) {
                        (*mvol)(wi, hi, di) = q + 2; // q + 2 because 0 and 1 values are already used. 0 is background and 1 is current_selection in Powei's implementation. 
                        break;
                    }
                }
                // End mask modification

                // Compute min and max voxels. 
                if (v[k] >= mmax) mmax = v[k];
                if (v[k] <= mmin) mmin = v[k];

                k = k + 1;
            }
        }
    }

    result = miclose_volume(minc_volume);
    if (result == MI_ERROR) cout << "Error in closing MINC volume" << endl;

    cout << "Finished reading voxels" << endl;
    cout << "mmin = " << mmin << endl;
    cout << "mmax = " << mmax << endl;

    
    
    // Apply scaling for vvol and fvol. 
    #pragma omp parallel for
    for (int j = 0; j < w * h * d; j = j + 1) {
        (*vvol)(j) = vval_scale((v[j] - mmin) / (mmax - mmin));	
        (*fvol)(j) = (v[j] - mmin) / (mmax - mmin);
    }

    MRCData *vmrcdata = new MRCData();
    vmrcdata->scalex = 1;
    vmrcdata->scaley = 1;
    vmrcdata->scalez = 1;

    vmrcdata->offsetx = start[0]; 
    vmrcdata->offsety = start[1]; 
    vmrcdata->offsetz = start[2]; 

    vmrcdata->stepx = stepSize[0];
    vmrcdata->stepy = stepSize[1];
    vmrcdata->stepz = stepSize[2];

    vmrcdata->vmax = 1;
    vmrcdata->vmin = 0;

    vmrcdata->actualHeader = VolumeUtils::createMRCHeader(vvol, MODE_char);

    MRCData *fmrcdata = new MRCData();
    fmrcdata->scalex = 1;
    fmrcdata->scaley = 1;
    fmrcdata->scalez = 1;

    fmrcdata->offsetx = start[0]; 
    fmrcdata->offsety = start[1]; 
    fmrcdata->offsetz = start[2]; 

    fmrcdata->stepx = stepSize[0];
    fmrcdata->stepy = stepSize[1];
    fmrcdata->stepz = stepSize[2];

    fmrcdata->vmax = 1;
    fmrcdata->vmin = 0;

    fmrcdata->actualHeader = VolumeUtils::createMRCHeader(fvol, MODE_float);
    
    vvol->setMRCData(vmrcdata);
    fvol->setMRCData(fmrcdata);
    
    
    vvol->setFileName(getFileName(filename));
    vvol->setFilePath(getFilePath(filename));
    vvol->setSavableState(vpt::SAVABLE_OPENED_FILE);
    
    fvol->setFileName(getFileName(filename));
    fvol->setFilePath(getFilePath(filename));
    fvol->setSavableState(vpt::SAVABLE_OPENED_FILE);
    
    Project *proj = ProjectUtils::newProjectFromVolume(vvol,fvol); 
    proj->setMask(mvol); 
    
    BoundaryMap *bmap = new BoundaryMap(vvol);
    bmap->setProject(proj);
    proj->setBoundaryMap(bmap);
    
    return proj;
}

/**
 * Initializes some parameters for DC. 
 * 
 * @param filename - Location of input minc volume to be read. 
 * @param outputFile - Location of where the output .vtk file will be written
 * @param segmentValues - The volume labels for which surfaces will be generated
 * @param numSegments - The number of volume labels
 * @param smooth - Whether to apply a smoothing operation or not. 
 * @param blurredFilenames - A list of mincblur-ed files, to be used if smoothing is mincblur-ed based. 
 */
void init(std::string filename, std::string outputFile, bool inWorldCoordinates, double segmentValues[], int numSegments, bool smooth, std::string blurredFilenames[]){
    Project* proj = readVolumeData(filename, segmentValues, numSegments);    
    SegmentMap* ret = new SegmentMap(false); 
    
    // Base segment is background and has value 0. 
    SegmentMapEntry *baseSegment = ret->createSegmentMapEntry();
    baseSegment->setName("BaseSegment");
    baseSegment->setSegmentCode(0);
    baseSegment->setState(1);
    baseSegment->setColor(osg::Vec4(0.5, 0.5, 0.5, 1));
   
    SegmentMapEntry *currentSelection = ret->createSegmentMapEntry();
    currentSelection->setName("CurrentSelection");
    currentSelection->setSegmentCode(1);
    currentSelection->setState(1);
    currentSelection->setColor(osg::Vec4(0.6, 0, 0, 1));

    // Generate segment maps for specified number of materials. 
    for(int j = 2; j < numSegments + 2; j++) {		
        SegmentMapEntry* ent = ret->createSegmentMapEntry(); 

        int code = j; // j or 0. 
        float r = ((double) rand() / (RAND_MAX));
        float g = ((double) rand() / (RAND_MAX)); 
        float b = ((double) rand() / (RAND_MAX));
        int state = 1;

        std::stringstream s;
        s << "color-";
        s << j;
        ent->setName(s.str()); 
        ent->setSegmentCode(code); 
        ent->setState(state); 
        ent->setColor(osg::Vec4(r, g, b, state)); 
    }
    proj->setSegmentMap(ret); 
    
    // Turn smoothing on or off. 
    // updateBoundaryWithMask() is Powei's smoothing, based on Gaussian curve
    // updateBoundaryWithMincBlurredVolume() is my smoothing using mincblur-ed images
    if (smooth) {
        BoundaryMapUtils::updateBoundaryWithMask(proj->getBoundaryMap(), proj->getSegmentMap(), proj->getMask());
        
        //BoundaryMapUtils::updateBoundaryWithMincBlurredVolume(proj->getBoundaryMap(), proj->getSegmentMap(), proj->getMask(), blurredFilenames);
        
        NormalMap *normalMap = proj->getNormalMap();
        if (normalMap) {
            NormalMapUtils::updateAllNormals(normalMap);
        }
    }
    buildDualContour(proj->getMask(), proj->getBoundaryMap()->getFunction(), proj->getVolume()->getMRCData(), proj->getSegmentMap(), filename, outputFile, inWorldCoordinates, segmentValues, numSegments); 

    return;
}



int main() {
    
    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/MultiMaterialCase4.mnc";
    //std::string filename = "data/TwoBoxesEqualSideBySide.mnc";
    //std::string filename = "data/FourBoxesEqualSideBySide.mnc";
    //std::string filename = "data\\labels_on_colin_Nov2010_minc2_mincreshape_xyz.mnc";
    //std::string filename = "data/Object26_cropped.mnc";
    //std::string filename = "data/Object26_part2.mnc";
    std::string filename = "data/LargeNonManifoldCube.mnc";
    //std::string filename = "data/Case4.mnc";
    //std::string filename = "data/Object26_part2.mnc";
    //std::string filename = "data/Object26_part1.mnc";
    //std::string filename = "data/Object26_part3_padded.mnc";
    //std::string filename = "/home/trash001/Desktop/asd.mnc"; 
    //std::string filename = "/home/trash001/NetBeansProjects/Zhang_HybridDC/Zhang_Hybrid_DualContouring/data/FourBoxesEqualSideBySide.mnc";
    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/R(F(S))/Object26.mnc";
    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/R(F(S))/Object27.mnc";
    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/R(F(S))/Object28.mnc";
    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/R(F(S))/Object26_close_open_reshape.mnc";
    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/R(F(S))/Object27_close_open_reshape.mnc";
    //std::string filename = "/home/trash001/Desktop/mincmorphing/Object26-27-28_morphed_combined.mnc";
    //std::string filename = "/home/trash001/Desktop/Cases/cropped1.mnc";
    //std::string filename = "/home/trash001/Desktop/Object27_28_part1_padded_11_16_6_padded.mnc";
    //std::string filename = "data/Object27_28_part1_padded.mnc";
    //std::string filename = "/home/trash001/Desktop/MultimaterialVolumeOverlappingExample.mnc";
    //std::string filename = "/home/trash001/Desktop/BrainWeb_Atlases/subject04_crisp_v_minc2_mincreshape_xyz.mnc"; 
    //std::string filename = "/home/trash001/Desktop/MalarAtlas_SeparatedObjects/Object27_close_open.mnc";
    //std::string filename = "/home/trash001/Desktop/Cases/Case4-1-1.mnc";
    //std::string filename = "/home/trash001/Desktop/Case7_padded_modified.mnc";
    //std::string filename = "/media/trash001/SAMSUNG USB/MallarAtlasLeftRight/labels_on_colin_Nov2010_minc2_mincreshape_xyz_RightSideCropped.mnc";
    
    //const std::string outputFile = "data/TwoBoxesEqualSideBySide_MM2M_VoxelCoordinates.vtk";
    //const std::string outputFile = "data/Object_26_MM2M_PoweiUnsmoothed_Minimizers_NotClamped.vtk";
    //const std::string outputFile = "/media/trash001/SAMSUNG USB/MallarAtlasLeftRight/Objects_1-4-5-11-12_RightSide_World.vtk";
    const std::string outputFile = "data/IntersectionFree_Manifold_DC_Output.vtk";
    
    std::string blurredFilenames[] = {""
        //"data/mincblur labels_on_colin_Nov2010/Object26_material0_stddev7.0_blur.mnc",
        //"",
        //"data/mincblur labels_on_colin_Nov2010/Object26_material1_stddev0.2_blur.mnc",
    };
    

    int numSeg = 2;
    double segmentValues[] = {0.3, 0.5};
    bool wordlCoord = true;
    bool smooth = false;
    
    
    init(filename, outputFile, wordlCoord, segmentValues, numSeg, smooth, blurredFilenames);
    cout << "Contouring completed" << endl;
    
    
	std::cout << "Execution completed." << std::endl;
	getchar();
	
            
    return 0;
}
