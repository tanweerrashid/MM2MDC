
#include "TCell.h"
#include <float.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "Common.h"
#include "DualContouring.h"

using namespace std;
using namespace vpt;

TCell::TCell(OCell *pc, int lastIndex, int corner1, int corner2, int corner3) {
    m_BiPoint = -1;
    
    parentOCell = pc;
    index[0] = (int)parentOCell->m_Ind[0];
    index[1] = (int)parentOCell->m_Ind[1];
    index[2] = (int)parentOCell->m_Ind[2];
    index[3] = lastIndex;

    // First point of the tetrahedra, i.e., the center of the grid cell. 
    points[0][0] = parentOCell->corners[8][0];
    points[0][1] = parentOCell->corners[8][1];
    points[0][2] = parentOCell->corners[8][2];

    // Second point of the tetrahedra. 
    points[1][0] = parentOCell->corners[corner1][0];
    points[1][1] = parentOCell->corners[corner1][1];
    points[1][2] = parentOCell->corners[corner1][2];

    // Third point of the tetrahedra. 
    points[2][0] = parentOCell->corners[corner2][0];
    points[2][1] = parentOCell->corners[corner2][1];
    points[2][2] = parentOCell->corners[corner2][2];

    // Fourth point of the tetrahedra. 
    points[3][0] = parentOCell->corners[corner3][0];
    points[3][1] = parentOCell->corners[corner3][1];
    points[3][2] = parentOCell->corners[corner3][2];

    // Center of the tetrahedra. 
    points[4][0] = (points[0][0] + points[1][0] + points[2][0] + points[3][0]) / 4;
    points[4][1] = (points[0][1] + points[1][1] + points[2][1] + points[3][1]) / 4;
    points[4][2] = (points[0][2] + points[1][2] + points[2][2] + points[3][2]) / 4;
    
    int mat[4] = {
        parentOCell->materialVals[8], 
        parentOCell->materialVals[corner1], 
        parentOCell->materialVals[corner2], 
        parentOCell->materialVals[corner3]
    };
    setMaterialValues(mat);
};
        
TCell::TCell() {
}

TCell::~TCell() {
    //points = NULL;
    //index = NULL;
    parentOCell = NULL;
    m_BiPoint = -1;
};

void TCell::setBPoint(int index) { 
    m_BiPoint = index; 
}

int TCell::getBPoint() { 
    return m_BiPoint; 
}

OCell* TCell::getParent() {
    return parentOCell;
}

int* TCell::getIndex() {
    int *ind = new int[4];
    ind[0] = index[0];
    ind[1] = index[1];
    ind[2] = index[2];
    ind[3] = index[3];
    
    return ind;
}

void TCell::setMaterialValues(int mat[4]) {
    for (int i = 0; i < 4; i++) {
        materialVals[i] = mat[i];
    }
}

int* TCell::getMaterialValues() {
    int *mat = new int[4];
    for (int i = 0; i < 4; i++) {
        mat[i] = materialVals[i];
    }
    return mat;
}

bool TCell::hasSignChangeEdge() {
    
    // Checks whether any edge is a sign change edge
    if ((materialVals[0] != materialVals[1]) || (materialVals[0] != materialVals[2]) || (materialVals[0] != materialVals[3]) || 
    (materialVals[1] != materialVals[2]) || (materialVals[2] != materialVals[3]) || (materialVals[3] != materialVals[1])) {
        return true;
    }
    else {
        return false;
    }

    
    
//    // Checks only the edges that are edges of the cube as well as the face diagonal
//    if ((materialVals[1] != materialVals[2]) || (materialVals[1] != materialVals[3]) || (materialVals[2] != materialVals[3])) {
//        return true;
//    }
//    else {
//        return false;
//    }
}

float* TCell::getCenter() {
    float *c = new float[3];
    c[0] = points[4][0];
    c[1] = points[4][1];
    c[2] = points[4][2];
    return c;
}

bool TCell::hasEdge(float p[3], float q[3]) {
    int count = 0;
    
    // Make sure p and q are not the same point. 
    if (isSameFloat(p[0], q[0]) == false || isSameFloat(p[1], q[1]) == false || isSameFloat(p[2], q[2]) == false) {
        for (int i = 0; i < 4; i++) {
            if ((isSameFloat(p[0], points[i][0]) == true) && (isSameFloat(p[1], points[i][1]) ==  true) && (isSameFloat(p[2], points[i][2]) == true)) {
                count = count + 1;
            }
            if ((isSameFloat(q[0], points[i][0]) == true) && (isSameFloat(q[1], points[i][1]) ==  true) && (isSameFloat(q[2], points[i][2]) == true)) {
                count = count + 1;
            }
        }
    }
    
    if (count == 2) {
        return true;
    }
    else {
        return false;
    }
}

bool TCell::hasVertex(float p[3]) {
    bool ret = false;
    for (int i = 0; i < 4; i++) {
        //if (p[0] == points[i][0] && p[1] == points[i][1] && p[2] == points[i][2]) {
        if ((isSameFloat(p[0], points[i][0]) == true) && (isSameFloat(p[1], points[i][1]) == true) && (isSameFloat(p[2], points[i][2]) == true)) {
            ret = true;
            break;
        }
    }
    return ret;
}

bool TCell::isSubsetOfFace(float p0[3], float p1[3], float p2[3], float p3[3]) {
    int count = 0;
    for (int i = 1; i < 4; i++) {
        if ((isSameFloat(p0[0], points[i][0]) == true) && (isSameFloat(p0[1], points[i][1]) == true) && (isSameFloat(p0[2], points[i][2]) == true)) {
            count = count + 1;
        }
        else if ((isSameFloat(p1[0], points[i][0]) == true) && (isSameFloat(p1[1], points[i][1]) == true) && (isSameFloat(p1[2], points[i][2]) == true)) {
            count = count + 1;
        }
        else if ((isSameFloat(p2[0], points[i][0]) == true) && (isSameFloat(p2[1], points[i][1]) == true) && (isSameFloat(p2[2], points[i][2]) == true)) {
            count = count + 1;
        }
        else if ((isSameFloat(p3[0], points[i][0]) == true) && (isSameFloat(p3[1], points[i][1]) == true) && (isSameFloat(p3[2], points[i][2]) == true)) {
            count = count + 1;
        }
    }
    
    if (count == 3) return true;
    else return false;
}

float** TCell::getCorners() {
    float **ret = new float*[4];
    //cout << parentOCell->getHasTetrahedra() << endl;
    //cout << "index: " << getIndex()[0] << ", " << getIndex()[1] << ", " << getIndex()[2] << ", " << getIndex()[3] << endl;
    
    for (int i = 0; i < 4; i++) {
        ret[i] = new float[3];
        ret[i][0] = points[i][0];
        ret[i][1] = points[i][1];
        ret[i][2] = points[i][2];
    }
    
    return ret;
}


/**
 * Determines whether two points are diagonally opposite points of the parent 
 * OCell making up a face diagonal. 
 * p and q are assumed to be points of the TCell/OCell. 
 * 
 * @param p - First point of an edge
 * @param q - Second point of an edge
 * @return - True if p and q are diagonally opposite points of the parent OCell, false otherwise. 
 */
bool TCell::isFaceDiagonal(float p[3], float q[3], int pq_mat[2]) {
    bool ret = false;
    
    int diagonalArray[12][2] = { // Two possible diagonals for each face
        {0, 3}, {1, 2}, // Face 0, 1, 3, 2
        {3, 6}, {2, 7}, // Face 2, 3, 7, 6
        {6, 5}, {4, 7}, // Face 6, 7, 5, 4
        {0, 5}, {1, 4}, // Face 0, 4, 5, 1
        {0, 6}, {4, 2}, // Face 0, 2, 6, 4
        {1, 7}, {3, 5}, // Face 1, 3, 7, 5
    };
    
    OCell *oc = this->getParent();
    float **corners = oc->getCorners();
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
    
    // If both pIndex and qIndex have values >= 0
    if (pIndex >= 0 && qIndex >= 0) {
        for (int j = 0; j < 12; j++) {
            if ((pIndex == diagonalArray[j][0] && qIndex == diagonalArray[j][1]) || (pIndex == diagonalArray[j][1] && qIndex == diagonalArray[j][0])) {
                ret = true;
                break;
            }
        }
    }
    
    int *mats = oc->getMaterialValues();
    pq_mat[0] = mats[pIndex];
    pq_mat[1] = mats[qIndex];
    
    return ret;
}

void TCell::vtkFile(std::string path) {
    ofstream file;
    file.open(path.c_str());
    
    file << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS " << 5 << " float\n";
    for (int i = 0; i < 5; i++) {        
        file << points[i][0] << " " << points[i][1] << " " << points[i][2] << "\n";
    }
    
    file << "POLYGONS " << 4 << " " << (4 * 4) << "\n";
    file << "3 0 1 2\n";
    file << "3 0 2 3\n";
    file << "3 0 3 1\n";
    file << "3 1 2 3\n";
    
    file << "VERTICES " << 5 << " " << (5 * 2) << "\n";
    file << "1 0\n";
    file << "1 1\n";
    file << "1 2\n";
    file << "1 3\n";
    file << "1 4\n";
    file << "\n";
    
    file << "POINT_DATA 5\n";
    file << "COLOR_SCALARS Colors 1\n";
    file << materialVals[0] << "\n";
    file << materialVals[1] << "\n";
    file << materialVals[2] << "\n";
    file << materialVals[3] << "\n";
    file << "0\n";
    
    
//    file << "POINT_DATA " << 5 << "\nSCALARS scalars float\nLOOKUP_TABLE default\n";
//    file << materialVals[0] << "\n";
//    file << materialVals[1] << "\n";
//    file << materialVals[2] << "\n";
//    file << materialVals[3] << "\n";
//    file << "0\n" << endl;
    
    file.close();
}
