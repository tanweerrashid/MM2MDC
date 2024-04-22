#include "Common.h"
#include "OCell.h"
#include "TCell.h"
#include <fstream>
#include <minc2.h>
#include <iostream>

using namespace std;

/************ begin OCell class *************/
OCell::OCell(unsigned char level, uint ind[3], OCell* parent, double wpc[3]){
    m_Children = NULL; 
    m_Level = level; 

    for (int j = 0; j < 3; j++) {
        m_Ind[j] = ind[j]; 
    }

    m_Parent = parent; 
    m_BiPoint = -1;
    isAmbiguous = false;
    hasTetrahedra = false;
    
    int pi = (int)m_Ind[0];
    int pj = (int)m_Ind[1];
    int pk = (int)m_Ind[2];
    
    corners[0][0] = (pi * wpc[0]);                corners[0][1] = (pj * wpc[1]);            corners[0][2] = (pk * wpc[2]);
    corners[1][0] = (pi * wpc[0]) + wpc[0];       corners[1][1] = (pj * wpc[1]);            corners[1][2] = (pk * wpc[2]);
    corners[2][0] = (pi * wpc[0]);                corners[2][1] = (pj * wpc[1]) + wpc[1];   corners[2][2] = (pk * wpc[2]);
    corners[3][0] = (pi * wpc[0]) + wpc[0];       corners[3][1] = (pj * wpc[1]) + wpc[1];   corners[3][2] = (pk * wpc[2]);
    corners[4][0] = (pi * wpc[0]);                corners[4][1] = (pj * wpc[1]);            corners[4][2] = (pk * wpc[2]) + wpc[2];
    corners[5][0] = (pi * wpc[0]) + wpc[0];       corners[5][1] = (pj * wpc[1]);            corners[5][2] = (pk * wpc[2]) + wpc[2];
    corners[6][0] = (pi * wpc[0]);                corners[6][1] = (pj * wpc[1]) + wpc[1];   corners[6][2] = (pk * wpc[2]) + wpc[2];
    corners[7][0] = (pi * wpc[0]) + wpc[0];       corners[7][1] = (pj * wpc[1]) + wpc[1];   corners[7][2] = (pk * wpc[2]) + wpc[2];
    
    corners[8][0] = 0; corners[8][1] = 0; corners[8][2] = 0;
    for (int i = 0; i < 8; i++) {
        corners[8][0] = corners[8][0] + corners[i][0];
        corners[8][1] = corners[8][1] + corners[i][1];
        corners[8][2] = corners[8][2] + corners[i][2];
    }
    
    corners[8][0] = corners[8][0] / 8; 
    corners[8][1] = corners[8][1] / 8; 
    corners[8][2] = corners[8][2] / 8; 
    
    for (int i = 0; i < 8; i++) {
        interiorEdgeProcessed[i] = false;
    }
}

OCell::OCell(unsigned char level, uint ix, uint iy, uint iz, OCell* parent, double wpc[3]){
    m_Children = NULL; 
    m_Level = level; 

    m_Ind[0] = ix; 	
    m_Ind[1] = iy; 
    m_Ind[2] = iz; 
    
    m_Parent = parent; 
    m_BiPoint = -1; 
    isAmbiguous = false;
    hasTetrahedra = false;
    
    int pi = (int)m_Ind[0];
    int pj = (int)m_Ind[1];
    int pk = (int)m_Ind[2];
    
    corners[0][0] = (pi * wpc[0]);                corners[0][1] = (pj * wpc[1]);            corners[0][2] = (pk * wpc[2]);
    corners[1][0] = (pi * wpc[0]) + wpc[0];       corners[1][1] = (pj * wpc[1]);            corners[1][2] = (pk * wpc[2]);
    corners[2][0] = (pi * wpc[0]);                corners[2][1] = (pj * wpc[1]) + wpc[1];   corners[2][2] = (pk * wpc[2]);
    corners[3][0] = (pi * wpc[0]) + wpc[0];       corners[3][1] = (pj * wpc[1]) + wpc[1];   corners[3][2] = (pk * wpc[2]);
    corners[4][0] = (pi * wpc[0]);                corners[4][1] = (pj * wpc[1]);            corners[4][2] = (pk * wpc[2]) + wpc[2];
    corners[5][0] = (pi * wpc[0]) + wpc[0];       corners[5][1] = (pj * wpc[1]);            corners[5][2] = (pk * wpc[2]) + wpc[2];
    corners[6][0] = (pi * wpc[0]);                corners[6][1] = (pj * wpc[1]) + wpc[1];   corners[6][2] = (pk * wpc[2]) + wpc[2];
    corners[7][0] = (pi * wpc[0]) + wpc[0];       corners[7][1] = (pj * wpc[1]) + wpc[1];   corners[7][2] = (pk * wpc[2]) + wpc[2];
    
    corners[8][0] = 0; corners[8][1] = 0; corners[8][2] = 0;
    for (int i = 0; i < 8; i++) {
        corners[8][0] = corners[8][0] + corners[i][0];
        corners[8][1] = corners[8][1] + corners[i][1];
        corners[8][2] = corners[8][2] + corners[i][2];
    }
    
    corners[8][0] = corners[8][0] / 8; 
    corners[8][1] = corners[8][1] / 8; 
    corners[8][2] = corners[8][2] / 8; 
    
    for (int i = 0; i < 8; i++) {
        interiorEdgeProcessed[i] = false;
    }
}

OCell::OCell() { 
    m_Children = NULL; 
    m_Parent = NULL; 
    m_BiPoint = -1;
    tetrahedra = NULL;
    isAmbiguous = false;
    hasTetrahedra = false;
    
    for (int i = 0; i < 8; i++) {
        interiorEdgeProcessed[i] = false;
    }
}

OCell::~OCell() {
    if (m_Children != NULL) {
        for (int j = 0; j < 8; j++) {
            if (m_Children[j] != NULL) 
                delete m_Children[j]; 
        }
        delete [] m_Children; 
    }
}

void OCell::print(){
    cout << "Cell: (" << m_Ind[0] << "," << m_Ind[1] << "," << m_Ind[2] << ")" 
            << " Level: " << (int)m_Level << " Parent: " << m_Parent 
            << " Has_Children: " << this->hasChildren() << " bpoint: " << m_BiPoint;
    
    if (m_Parent != NULL) {
        cout << " childInd: "; 
        for (int j = 0; j < 8; j++)
            if (m_Parent->m_Children[j] == this) {
                cout << j;
                break;
        }
    }
cout << endl;
}
    
void OCell::setBPoint(int index) { 
    m_BiPoint = index; 
}

int OCell::getBPoint() { 
    return m_BiPoint; 
}

bool OCell::hasChildren() { 
    return m_Children != NULL; 
}

void OCell::setHasChildren() { 
    m_Children = new OCell*[8]; 
}

int getFaceCenterMaterialValues_Direct_MINC_Query(float faceCenter[3]) {
    mihandle_t minc_volume;
    double voxel;
    int result;
    double world_location[3];
    double dvoxel_location[3];
    misize_t voxel_location[3];

    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/MultiMaterialCase4.mnc";
    //std::string filename = "data/TwoBoxesEqualSideBySide.mnc";
    //std::string filename = "data/FourBoxesEqualSideBySide.mnc";
    //std::string filename = "data\\labels_on_colin_Nov2010_minc2_mincreshape_xyz.mnc";
    //std::string filename = "data/Object26_part1.mnc";
    //std::string filename = "data/Object26_part2.mnc";
    //std::string filename = "data/Object26_part3_padded.mnc";
    //std::string filename = "data/Object26_cropped.mnc";
    //std::string filename = "/home/trash001/Desktop/Object26_part2.mnc";
    std::string filename = "data/LargeNonManifoldCube.mnc";
    //std::string filename = "data/Case4.mnc";
    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/R(F(S))/Object26.mnc";
    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/R(F(S))/Object27.mnc";
    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/R(F(S))/Object28.mnc";
    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/R(F(S))/Object26_close_open_reshape.mnc";
    //std::string filename = "/home/trash001/NetBeansProjects/MINC/data/R(F(S))/Object27_close_open_reshape.mnc";
    //std::string filename = "/home/trash001/Desktop/mincmorphing/Object28_open.mnc";
    //std::string filename = "/home/trash001/Desktop/Cases/cropped1.mnc";
    //std::string filename = "/home/trash001/Desktop/Object27_28_part1_padded_11_16_6_padded.mnc";
    //std::string filename = "data/Object27_28_part1_padded.mnc";
    //std::string filename = "/home/trash001/Desktop/MultimaterialVolumeOverlappingExample.mnc";
    //std::string filename = "/home/trash001/Desktop/BrainWeb_Atlases/subject04_crisp_v_minc2_mincreshape_xyz.mnc";
    //std::string filename = "/home/trash001/Desktop/MalarAtlas_SeparatedObjects/Object27_close_open.mnc";
    //std::string filename = "/home/trash001/Desktop/Cases/Case4-1-1.mnc";
    //std::string filename = "/home/trash001/Desktop/Case7_padded_modified.mnc";
    //std::string filename = "/media/trash001/SAMSUNG USB/MallarAtlasLeftRight/labels_on_colin_Nov2010_minc2_mincreshape_xyz_RightSideCropped.mnc";
    
    result = miopen_volume(filename.c_str(), MI2_OPEN_READ, &minc_volume);
    if (result != MI_NOERROR) {
        cout << "Error opening input file: " << result << endl;
    }
    
    /** Convert the voxel coordinates into world coordinates. */
    midimhandle_t dimension[3];
    miget_volume_dimensions(minc_volume, MI_DIMCLASS_SPATIAL, MI_DIMATTR_ALL, MI_DIMORDER_FILE, 3, dimension);
    
    double xcosine[3], ycosine[3], zcosine[3];
    miget_dimension_cosines(dimension[0], xcosine);
    miget_dimension_cosines(dimension[1], ycosine);
    miget_dimension_cosines(dimension[2], zcosine);
    
    double stepSize[3];
    result = miget_dimension_separations(dimension, MI_ORDER_FILE, 3, stepSize);
    if (result == MI_ERROR) cout << "Error in retrieving step sizes" << endl;
    
    double start[3];
    result = miget_dimension_starts(dimension, MI_ORDER_FILE, 3, start);
    if (result == MI_ERROR) cout << "Error in retrieving volume origin" << endl;
    
    double vx = (xcosine[0] * stepSize[0] * faceCenter[0]) + (ycosine[0] * stepSize[1] * faceCenter[1]) + (zcosine[0] * stepSize[2] * faceCenter[2]) + start[0];
    double vy = (xcosine[1] * stepSize[0] * faceCenter[0]) + (ycosine[1] * stepSize[1] * faceCenter[1]) + (zcosine[1] * stepSize[2] * faceCenter[2]) + start[1];
    double vz = (xcosine[2] * stepSize[0] * faceCenter[0]) + (ycosine[2] * stepSize[1] * faceCenter[1]) + (zcosine[2] * stepSize[2] * faceCenter[2]) + start[2];
    
    world_location[0] = vx;
    world_location[1] = vy;
    world_location[2] = vz;

    result = miconvert_world_to_voxel(minc_volume, world_location, dvoxel_location);
    if (result != MI_NOERROR) cout << "Error in converting world to voxel coordinates: " << result << endl;
    
    for (int i = 0; i < 3; i++) { 
        voxel_location[i] = (misize_t)dvoxel_location[i];
    }

    //cout << "Voxel location of world coordinate (" << world_location[0] << ", " << world_location[1] << ", " << world_location[2] << ") is [" << voxel_location[0] << ", " << voxel_location[1] << ", " << voxel_location[2] << "]" << endl;
    
    result = miget_voxel_value(minc_volume, voxel_location, 3, &voxel);
    if (result != MI_NOERROR) cout << "Error in retrieving voxel value: " << result << endl;
    
    miclose_volume(minc_volume);
    
    if (voxel == 0.3) {
        voxel = 2;
    }
    else if (voxel == 0.5) {
        voxel = 3;
    }
    /*else if (voxel == 28) {
        voxel == 4;
    }
    else if (voxel == 39) {
        voxel == 5;
    }
    else if (voxel == 47) {
        voxel == 6;
    }
    /*else if (voxel == 27) {
        voxel == 7;
    }
    else if (voxel == 28) {
        voxel == 8;
    }
    else if (voxel == 37) {
        voxel == 9;
    }
    else if (voxel == 39) {
        voxel == 10;
    }
    else if (voxel == 47) {
        voxel == 11;
    }
    else if (voxel == 48) {
        voxel == 12;
    }
    else if (voxel == 49) {
        voxel == 13;
    }*/
    else {
        voxel = 0;
    }
    
    //cout << "voxel value: " << voxel << endl;
    return ((int)voxel);
}

int getFaceCenterMaterialValues_TrilinearInterpolation(float faceCenter[3], int *corners) {
    float x = faceCenter[0], y = faceCenter[1], z = faceCenter[2];
    int V000 = corners[0];
    int V100 = corners[1];
    int V010 = corners[2];
    int V110 = corners[3];
    int V001 = corners[4];
    int V101 = corners[5];
    int V011 = corners[6];
    int V111 = corners[7];
    //                       0       1       2       3       4       5       6       7
    //const int ac[8][3]={{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}}; 
    
    int Vxyz = V000 * (1 - x) * (1 - y) * (1 - z) + 
        V100 * x * (1 - y) * (1 - z) + 
        V010 * (1 - x) * y * (1 - z) + 
        V001 * (1 - x) * (1 - y) * z +
        V101 * x * (1 - y) * z + 
        V011 * (1 - x) * y * z + 
        V110 * x * y * (1 - z) + 
        V111 * x * y * z;

    return Vxyz;
}

void OCell::setTetrahedralCells() {
    tetrahedra = new TCell*[12];
    
    // Compute center of face coordinates. 
    float f0123[3], f4567[3], f0145[3], f2367[3], f0246[3], f1357[3];
    int mv0123, mv4567, mv0145, mv2367, mv0246, mv1357; // Material values of the face centers 
    
    f0123[0] = (corners[0][0] + corners[1][0] + corners[2][0] + corners[3][0]) / 4;
    f0123[1] = (corners[0][1] + corners[1][1] + corners[2][1] + corners[3][1]) / 4;
    f0123[2] = (corners[0][2] + corners[1][2] + corners[2][2] + corners[3][2]) / 4;
    
    f4567[0] = (corners[4][0] + corners[5][0] + corners[6][0] + corners[7][0]) / 4;
    f4567[1] = (corners[4][1] + corners[5][1] + corners[6][1] + corners[7][1]) / 4;
    f4567[2] = (corners[4][2] + corners[5][2] + corners[6][2] + corners[7][2]) / 4;
    
    f0145[0] = (corners[0][0] + corners[1][0] + corners[4][0] + corners[5][0]) / 4;
    f0145[1] = (corners[0][1] + corners[1][1] + corners[4][1] + corners[5][1]) / 4;
    f0145[2] = (corners[0][2] + corners[1][2] + corners[4][2] + corners[5][2]) / 4;
    
    f2367[0] = (corners[2][0] + corners[3][0] + corners[6][0] + corners[7][0]) / 4;
    f2367[1] = (corners[2][1] + corners[3][1] + corners[6][1] + corners[7][1]) / 4;
    f2367[2] = (corners[2][2] + corners[3][2] + corners[6][2] + corners[7][2]) / 4;
    
    f0246[0] = (corners[0][0] + corners[2][0] + corners[4][0] + corners[6][0]) / 4;
    f0246[1] = (corners[0][1] + corners[2][1] + corners[4][1] + corners[6][1]) / 4;
    f0246[2] = (corners[0][2] + corners[2][2] + corners[4][2] + corners[6][2]) / 4;
    
    f1357[0] = (corners[1][0] + corners[3][0] + corners[5][0] + corners[7][0]) / 4;
    f1357[1] = (corners[1][1] + corners[3][1] + corners[5][1] + corners[7][1]) / 4;
    f1357[2] = (corners[1][2] + corners[3][2] + corners[5][2] + corners[7][2]) / 4;
   
    mv0123 = getFaceCenterMaterialValues_Direct_MINC_Query(f0123);
    mv4567 = getFaceCenterMaterialValues_Direct_MINC_Query(f4567);
    mv0145 = getFaceCenterMaterialValues_Direct_MINC_Query(f0145);
    mv2367 = getFaceCenterMaterialValues_Direct_MINC_Query(f2367);
    mv0246 = getFaceCenterMaterialValues_Direct_MINC_Query(f0246);
    mv1357 = getFaceCenterMaterialValues_Direct_MINC_Query(f1357);
       

//    int v = 2;
//    mv0123 = v;
//    mv4567 = v;
//    mv0145 = v;
//    mv2367 = v;
//    mv0246 = v;
//    mv1357 = v;
    
//    int *corners = this->getMaterialValues();
//    mv0123 = getFaceCenterMaterialValues_TrilinearInterpolation(f0123, corners);
//    mv4567 = getFaceCenterMaterialValues_TrilinearInterpolation(f4567, corners);
//    mv0145 = getFaceCenterMaterialValues_TrilinearInterpolation(f0145, corners);
//    mv2367 = getFaceCenterMaterialValues_TrilinearInterpolation(f2367, corners);
//    mv0246 = getFaceCenterMaterialValues_TrilinearInterpolation(f0246, corners);
//    mv1357 = getFaceCenterMaterialValues_TrilinearInterpolation(f1357, corners);
    
    // cout << "For OCell [" << m_Ind[0] << ", " << m_Ind[1] << ", " << m_Ind[2] << "]" << endl;
    
//    cout << "Center of Face0123  {" << f0123[0] << ", " << f0123[1] << ", " << f0123[2] << "}" << endl;
//    cout << "Center of Face4567  {" << f4567[0] << ", " << f4567[1] << ", " << f4567[2] << "}" << endl;
//    cout << "Center of Face0145  {" << f0145[0] << ", " << f0145[1] << ", " << f0145[2] << "}" << endl;
//    cout << "Center of Face2367  {" << f2367[0] << ", " << f2367[1] << ", " << f2367[2] << "}" << endl;
//    cout << "Center of Face0246  {" << f0246[0] << ", " << f0246[1] << ", " << f0246[2] << "}" << endl;
//    cout << "Center of Face1357  {" << f1357[0] << ", " << f1357[1] << ", " << f1357[2] << "}" << endl;
    
//    cout << "Center of Face0123  has material index: " << mv0123 << endl;
//    cout << "Center of Face4567  has material index: " << mv4567 << endl;
//    cout << "Center of Face0145  has material index: " << mv0145 << endl;
//    cout << "Center of Face2367  has material index: " << mv2367 << endl;
//    cout << "Center of Face0246  has material index: " << mv0246 << endl;
//    cout << "Center of Face1357  has material index: " << mv1357 << endl;    
    
    
//    if ((int)m_Ind[0] == 4 && (int)m_Ind[1] == 4 && ((int)m_Ind[2] == 4)) { // || (int)m_Ind[2] == 1 || (int)m_Ind[2] == 2 || (int)m_Ind[2] == 3)) {
//        materialVals[8] = 0;
//    }
    
    // Face 0, 1, 2, 3
    if ((materialVals[0] == materialVals[3]) && (materialVals[0] != materialVals[2]) && (materialVals[0] != materialVals[1]) && (materialVals[0] > 1)) {
        if (materialVals[0] == mv0123) {
            tetrahedra[0] = new TCell(this, 0, 0, 1, 3); 
            tetrahedra[1] = new TCell(this, 1, 0, 2, 3); 
            //cout << "TCells 0 and 1 created" << endl;
        }
        else {
            tetrahedra[0] = new TCell(this, 0, 0, 1, 2); 
            tetrahedra[1] = new TCell(this, 1, 1, 2, 3); 
        }
    }
    else if ((materialVals[1] == materialVals[2]) && (materialVals[1] != materialVals[0]) && (materialVals[1] != materialVals[3]) && (materialVals[1] > 1)) {
        if (materialVals[1] == mv0123) {
            tetrahedra[0] = new TCell(this, 0, 0, 1, 2); 
            tetrahedra[1] = new TCell(this, 1, 1, 2, 3); 
            //cout << "TCells 0 and 1 created" << endl;
        }
        else {
            tetrahedra[0] = new TCell(this, 0, 0, 1, 3); 
            tetrahedra[1] = new TCell(this, 1, 0, 2, 3); 
            //cout << "TCells 0 and 1 created" << endl;
        }
    }
    else {
        tetrahedra[0] = new TCell(this, 0, 0, 1, 3); 
        tetrahedra[1] = new TCell(this, 1, 0, 2, 3); 
        //cout << "TCells 0 and 1 created" << endl;
    }

   // Face 0, 1, 4, 5
    if ((materialVals[0] == materialVals[5]) && (materialVals[0] != materialVals[4]) && (materialVals[0] != materialVals[1]) && (materialVals[0] > 1)) {
        if (materialVals[0] == mv0145) {
            tetrahedra[2] = new TCell(this, 2, 0, 1, 5); 
            tetrahedra[3] = new TCell(this, 3, 0, 4, 5); 
            //cout << "TCells 2 and 3 created" << endl;
        }
        else {
            tetrahedra[2] = new TCell(this, 2, 1, 5, 4); 
            tetrahedra[3] = new TCell(this, 3, 1, 0, 4); 
            //cout << "TCells 2 and 3 created" << endl;
        }
    }
    else if ((materialVals[4] == materialVals[1]) && (materialVals[4] != materialVals[0]) && (materialVals[4] != materialVals[5]) && (materialVals[4] > 1)) {
        if (materialVals[4] == mv0145) {
            tetrahedra[2] = new TCell(this, 2, 1, 5, 4); 
            tetrahedra[3] = new TCell(this, 3, 1, 0, 4); 
            //cout << "TCells 2 and 3 created" << endl;
        }
        else {
            tetrahedra[2] = new TCell(this, 2, 0, 1, 5); 
            tetrahedra[3] = new TCell(this, 3, 0, 4, 5); 
            //cout << "TCells 2 and 3 created" << endl;
        }
    }
    else {
        tetrahedra[2] = new TCell(this, 2, 0, 1, 5); 
        tetrahedra[3] = new TCell(this, 3, 0, 4, 5); 
        //cout << "TCells 2 and 3 created" << endl;
    }
    
    // Face 2, 3, 6, 7
    if ((materialVals[2] == materialVals[7]) && (materialVals[2] != materialVals[3]) && (materialVals[2] != materialVals[6]) && (materialVals[2] > 1)) {
        if (materialVals[2] == mv2367) {
            tetrahedra[6] = new TCell(this, 6, 2, 3, 7); 
            tetrahedra[7] = new TCell(this, 7, 2, 6, 7);
            //cout << "TCells 6 and 7 created" << endl;
        }
        else {
            tetrahedra[6] = new TCell(this, 6, 2, 3, 6); 
            tetrahedra[7] = new TCell(this, 7, 3, 6, 7);
            //cout << "TCells 6 and 7 created" << endl;
        }
    }
    else if ((materialVals[3] == materialVals[6]) && (materialVals[3] != materialVals[2]) && (materialVals[3] != materialVals[7]) && (materialVals[3] > 1)) {
        if (materialVals[3] == mv2367) {
            tetrahedra[6] = new TCell(this, 6, 2, 3, 6); 
            tetrahedra[7] = new TCell(this, 7, 3, 6, 7);
            //cout << "TCells 6 and 7 created" << endl;
        }
        else {
            tetrahedra[6] = new TCell(this, 6, 2, 3, 7); 
            tetrahedra[7] = new TCell(this, 7, 2, 6, 7);
            //cout << "TCells 6 and 7 created" << endl;
        }
    }
    else {
        tetrahedra[6] = new TCell(this, 6, 2, 3, 7); 
        tetrahedra[7] = new TCell(this, 7, 2, 6, 7);
        //cout << "TCells 6 and 7 created" << endl;
    }
    
    // Face 0, 2, 4, 6
    if ((materialVals[0] == materialVals[6]) && (materialVals[0] != materialVals[2]) && (materialVals[0] != materialVals[4]) && (materialVals[0] > 1)) {
        if (materialVals[0] == mv0246) {
            tetrahedra[4] = new TCell(this, 4, 0, 2, 6); 
            tetrahedra[5] = new TCell(this, 5, 0, 4, 6); 
            //cout << "TCells 4 and 5 created" << endl;
        }
        else {
            tetrahedra[4] = new TCell(this, 4, 0, 2, 4); 
            tetrahedra[5] = new TCell(this, 5, 2, 4, 6); 
            //cout << "TCells 4 and 5 created" << endl;
        }
    }
    else if ((materialVals[2] == materialVals[4]) && (materialVals[2] != materialVals[0]) && (materialVals[2] != materialVals[6]) && (materialVals[2] > 1)) {
        if (materialVals[2] == mv0246) {
            tetrahedra[4] = new TCell(this, 4, 0, 2, 4); 
            tetrahedra[5] = new TCell(this, 5, 2, 4, 6); 
            //cout << "TCells 4 and 5 created" << endl;
        }
        else {
            tetrahedra[4] = new TCell(this, 4, 0, 2, 6); 
            tetrahedra[5] = new TCell(this, 5, 0, 4, 6); 
            //cout << "TCells 4 and 5 created" << endl;
        }
    }
    else {
        tetrahedra[4] = new TCell(this, 4, 0, 2, 6); 
        tetrahedra[5] = new TCell(this, 5, 0, 4, 6); 
        //cout << "TCells 4 and 5 created" << endl;
    }
    
    // Face 4, 5, 6, 7
    if ((materialVals[4] == materialVals[7]) && (materialVals[4] != materialVals[5]) && (materialVals[4] != materialVals[6]) && (materialVals[4] > 1)) {
        if (materialVals[4] == mv4567) {
            tetrahedra[10] = new TCell(this, 10, 4, 6, 7); 
            tetrahedra[11] = new TCell(this, 11, 4, 5, 7);
            //cout << "TCells 10 and 11 created" << endl;
        }
        else {
            tetrahedra[10] = new TCell(this, 10, 5, 6, 7); 
            tetrahedra[11] = new TCell(this, 11, 5, 4, 6);
            //cout << "TCells 10 and 11 created" << endl;
        }
    }
    else if ((materialVals[5] == materialVals[6]) && (materialVals[5] != materialVals[4]) && (materialVals[5] != materialVals[7]) && (materialVals[5] > 1)) {
        if (materialVals[5] == mv4567) {
            tetrahedra[10] = new TCell(this, 10, 5, 6, 7); 
            tetrahedra[11] = new TCell(this, 11, 5, 4, 6);
            //cout << "TCells 10 and 11 created" << endl;
        }
        else {
            tetrahedra[10] = new TCell(this, 10, 4, 6, 7); 
            tetrahedra[11] = new TCell(this, 11, 4, 5, 7);
            //cout << "TCells 10 and 11 created" << endl;
        }
    }
    else {
        tetrahedra[10] = new TCell(this, 10, 4, 6, 7); 
        tetrahedra[11] = new TCell(this, 11, 4, 5, 7);
        //cout << "TCells 10 and 11 created" << endl;
    }

    // Face 1, 3, 5, 7
    if ((materialVals[1] == materialVals[7]) && (materialVals[1] != materialVals[5]) && (materialVals[1] != materialVals[3]) && (materialVals[1] > 1)) {
        if (materialVals[1] == mv1357) {
            tetrahedra[8] = new TCell(this, 8, 1, 3, 7); 
            tetrahedra[9] = new TCell(this, 9, 1, 5, 7); 
            //cout << "TCells 8 and 9 created" << endl;
        }
        else {
            tetrahedra[8] = new TCell(this, 8, 1, 3, 5); 
            tetrahedra[9] = new TCell(this, 9, 3, 5, 7); 
            //cout << "TCells 8 and 9 created" << endl;
        }
    }
    else if ((materialVals[3] == materialVals[5]) && (materialVals[3] != materialVals[1]) && (materialVals[3] != materialVals[7]) && (materialVals[3] > 1)) {
        if (materialVals[3] == mv1357) {
            tetrahedra[8] = new TCell(this, 8, 1, 3, 5); 
            tetrahedra[9] = new TCell(this, 9, 3, 5, 7); 
            //cout << "TCells 8 and 9 created" << endl;
        }
        else {
            tetrahedra[8] = new TCell(this, 8, 1, 3, 7); 
            tetrahedra[9] = new TCell(this, 9, 1, 5, 7); 
            //cout << "TCells 8 and 9 created" << endl;
        }
    }
    else {
        tetrahedra[8] = new TCell(this, 8, 1, 3, 7); 
        tetrahedra[9] = new TCell(this, 9, 1, 5, 7); 
        //cout << "TCells 8 and 9 created" << endl;
    }
  
    setHasTetrahedra(true);
};

TCell** OCell::getTetrahedralCells() { 
    return tetrahedra; 
}

void OCell::setIsAmbiguous(bool b) { 
    isAmbiguous = b; 
}

bool OCell::getIsAmbiguous() { 
    return isAmbiguous; 
}

void OCell::setHasTetrahedra(bool b) { 
    hasTetrahedra = b; 
}

bool OCell::getHasTetrahedra() { 
    return hasTetrahedra; 
}

void OCell::setMaterialValues(int mat[9]) {
    for (int i = 0; i < 9; i++) {
        materialVals[i] = mat[i];
    }
}

int* OCell::getMaterialValues() {
    int *mat = new int[9];
    for (int i = 0; i < 9; i++) {
        mat[i] = materialVals[i];
    }
    return mat;
}

float** OCell::getCorners() {
    float **ret = new float*[9];
    for (int i = 0; i < 9; i++) {
        ret[i] = new float[3];
        ret[i][0] = corners[i][0];
        ret[i][1] = corners[i][1];
        ret[i][2] = corners[i][2];
    }
    return ret;
}

bool OCell::hasEdge(float p[3], float q[3]) {
    int edges[12][2] = {
        {0, 1}, {1, 3}, {3, 2}, {2, 0},
        {4, 5}, {5, 7}, {7, 6}, {6, 4}, 
        {1, 5}, {0, 4}, {2, 6}, {3, 7}, 
    };
    
    bool ret = false;
    
    for (int i = 0; i < 12; i++) {
        float r[3], s[3];
        r[0] = corners[edges[i][0]][0];
        r[1] = corners[edges[i][0]][1];
        r[2] = corners[edges[i][0]][2];
        
        s[0] = corners[edges[i][1]][0];
        s[1] = corners[edges[i][1]][1];
        s[2] = corners[edges[i][1]][2];
        
        if (((vpt::isSameFloat(p[0], r[0]) == true && vpt::isSameFloat(p[1], r[1]) == true && vpt::isSameFloat(p[2], r[2]) == true) && 
                (vpt::isSameFloat(q[0], s[0]) == true && vpt::isSameFloat(q[1], s[1]) == true && vpt::isSameFloat(q[2], s[2]) == true)) || 
                ((vpt::isSameFloat(p[0], s[0]) == true && vpt::isSameFloat(p[1], s[1]) == true && vpt::isSameFloat(p[2], s[2]) == true) && 
                (vpt::isSameFloat(q[0], r[0]) == true && vpt::isSameFloat(q[1], r[1]) == true && vpt::isSameFloat(q[2], r[2]) == true))) {
            ret = true;
            break;
        }
    }
    return ret;
}

int OCell::getCornerMaterialValue(float p[3]) {
    int ret = -1;
    
    for (int i = 0; i < 9; i++) {
        //if ((compareFloats(p[0], corners[i][0]) == true) && (compareFloats(p[1], corners[i][1]) == true) && (compareFloats(p[2], corners[i][2]) == true)) {
        if (vpt::isSameFloat(p[0], corners[i][0]) == true && vpt::isSameFloat(p[1], corners[i][1]) == true && vpt::isSameFloat(p[2], corners[i][2]) == true) {
        //if (p[0] == corners[i][0] && p[1] == corners[i][1] && p[2] == corners[i][2]) {
            ret = materialVals[i];
            break;
        }
    }
    
    return ret;
}

/**
 * Determines whether two points are diagonally opposite points of an 
 * OCell making up a face diagonal. 
 * p and q are assumed to be points of the TCell/OCell. 
 * 
 * @param p - First point of an edge
 * @param q - Second point of an edge
 * @return - True if p and q are diagonally opposite points of the parent OCell, false otherwise. 
 */
bool OCell::isFaceDiagonal(float p[3], float q[3], int faceDiagonalMaterial[2]) {
    bool ret = false;
    /*
    int diagonalArray[12][2] = { // Two possible diagonals for each face
        {0, 3}, {1, 2}, // Face 0, 1, 3, 2
        {3, 6}, {2, 7}, // Face 2, 3, 7, 6
        {6, 5}, {4, 7}, // Face 6, 7, 5, 4
        {0, 5}, {1, 4}, // Face 0, 4, 5, 1
        {0, 6}, {4, 2}, // Face 0, 2, 6, 4
        {1, 7}, {3, 5}, // Face 1, 3, 7, 5
    };
    
    //OCell *oc = this->getParent();
    float **corners = getCorners();
    int pIndex = -1, qIndex = -1;
    
    // Assign indices for points p and q 
    for (int i = 0; i < 8; i++) {
        if ((isSameFloat(p[0], corners[i][0]) == true) && (isSameFloat(p[1], corners[i][1]) == true) && (isSameFloat(p[2], corners[i][2]) == true)) {
            pIndex = i;
        }
        if ((isSameFloat(q[0], corners[i][0]) == true) && (isSameFloat(q[1], corners[i][1]) == true) && (isSameFloat(q[2], corners[i][2]) == true)) {
            qIndex = i;
        }
    }
    
    if (pIndex >= 0 && qIndex >= 0) {
        for (int j = 0; j < 12; j++) {
            if ((pIndex == diagonalArray[j][0] && qIndex == diagonalArray[j][1]) || (pIndex == diagonalArray[j][1] && qIndex == diagonalArray[j][0])) {
                ret = true;
                break;
            }
        }
    }
    */
    return ret;
}

/**
 * Identifies the specific ambiguity of OCells. 
 * 
 * @return 3, 4, 6, 7, 10, 12 and 13 for cases 3, 4, 6, 7, 10, 12 and 13, respectively. 
 * -1 for unambiguous cases. 
 */
int OCell::identityAmbiguousCase() {
    int ret = -1;
    
    
    
    return ret;
}

void OCell::setInteriorEdgeProcessed(int interiorEdgeIndex, bool b) {
    if (interiorEdgeIndex >= 0 && interiorEdgeIndex <= 7) {
        interiorEdgeProcessed[interiorEdgeIndex] = b;
    }
    else {
        cout << "ERROR \t\t InteriorEdgeIndex outside array index." << endl;
    }
}

bool OCell::getInteriorEdgeProcessed(int interiorEdgeIndex) {
    if (interiorEdgeIndex >= 0 && interiorEdgeIndex <= 7) {
        return interiorEdgeProcessed[interiorEdgeIndex];
    }
    else {
        cout << "ERROR \t\t InteriorEdgeIndex outside array index." << endl;
        return false;
    }
}


void OCell::setFaceProcessed(double face[4][3], bool b) {
    int faceIndices[6][4] = {
        {0, 4, 6, 2}, 
        {1, 5, 7, 3}, 
        {0, 2, 3, 1}, 
        {4, 6, 7, 5}, 
        {2, 3, 7, 6},
        {0, 1, 5, 4}
    };
    
}

//bool OCell::getFaceProcessed(double face[4][3]) {
//    
//}

void OCell::vtkFile(std::string path) {
//    // With center
//    ofstream file;
//    file.open(path.c_str());
//    
//    file << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS " << 9 << " float\n";
//    for (int i = 0; i < 9; i++) {        
//        file << corners[i][0] << " " << corners[i][1] << " " << corners[i][2] << "\n";
//    }
//    
//    file << "LINES " << 6 << " " << (6 * 6) << "\n";
//    file << "5 0 1 3 2 0\n";
//    file << "5 0 4 6 2 0\n";
//    file << "5 4 5 7 6 4\n";
//    file << "5 1 3 7 5 1\n";
//    file << "5 0 1 5 4 0\n";
//    file << "5 2 3 7 6 2\n";
//    
//    file << "VERTICES " << 1 << " " << (2) << "\n";
//    file << "1 8\n";
//    
//    file.close();
    
    
    // Without center
    ofstream file;
    file.open(path.c_str());
    
    file << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS " << 9 << " float\n";
    for (int i = 0; i < 9; i++) {        
        file << corners[i][0] << " " << corners[i][1] << " " << corners[i][2] << "\n";
    }
    
    file << "VERTICES " << 1 << " " << 2 << "\n";
    file << "1 8\n";
    
    file << "LINES " << 6 << " " << (6 * 6) << "\n";
    file << "5 0 1 3 2 0\n";
    file << "5 0 4 6 2 0\n";
    file << "5 4 5 7 6 4\n";
    file << "5 1 3 7 5 1\n";
    file << "5 0 1 5 4 0\n";
    file << "5 2 3 7 6 2\n";
    
    file << "\n";
    file << "POINT_DATA 9\n";
    file << "COLOR_SCALARS Colors 1\n";
    for (int i = 0; i < 9; i++) {
        file << materialVals[i] << "\n";
    }
    
    file.close();
}

/************** end OCell class ***************/
