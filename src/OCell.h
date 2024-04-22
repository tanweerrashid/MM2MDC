/* 
 * File:   OCell.h
 * Author: trash001
 *
 * Created on February 26, 2014, 12:51 PM
 */

#ifndef OCELL_H
#define	OCELL_H

#include "Common.h"
#include "TCell.h"
#include "DualContouring.h"

/* This class represents each cell in the octree */
class OCell{
    //friend class OctreeContouring; 
    
    public:
        OCell(unsigned char level, uint ind[3], OCell* parent, double wpc[3]); 
        OCell(unsigned char level, uint ix, uint iy, uint iz, OCell* parent, double wpc[3]); 
        OCell();
        ~OCell();

        void setBPoint(int index); 
        int getBPoint(); 
        bool hasChildren(); 
        void setHasChildren();

        void setTetrahedralCells() ;
        TCell** getTetrahedralCells();
        void setIsAmbiguous(bool b);
        bool getIsAmbiguous();
        void setHasTetrahedra(bool b);
        bool getHasTetrahedra();
        void vtkFile(std::string path); 
        void print(); 
        int* getMaterialValues();
        void setMaterialValues(int mat[9]); 
        float** getCorners();
        bool hasEdge(float p[3], float q[3]);
        int getCornerMaterialValue(float p[3]);
        bool isFaceDiagonal(float p[3], float q[3], int faceDiagonalMaterial[2]);
        int identityAmbiguousCase();
        
        void setInteriorEdgeProcessed(int interiorEdgeIndex, bool b);
        bool getInteriorEdgeProcessed(int interiorEdgeIndex);
        
        void setFaceProcessed(double face[4][3], bool b);
        bool getFaceProcessed(double face[4][3]);
        
        
        OCell** m_Children; 
        OCell* m_Parent; 
        unsigned char m_Level ; 
        uint m_Ind[3];
        int m_BiPoint;  // an index into the vertex table 	  

        float corners[9][3]; // First 8 are the coordinates of cell corners, the last one is the coordinate of the cell center. 
        int materialVals[9]; // The material values of the corners and center of the cube. First 8 elements are for corners, last one is for center. 
        
        
    private: 
        bool isAmbiguous; // Whether the cell is ambiguous or not. 
        bool hasTetrahedra; // Whether the tetrahedral cells have been created or not. 
        TCell **tetrahedra; 
        bool interiorEdgeProcessed[8]; // To determine whether an interior edge has been triangulated using the InteriorEdgeRule. To avoid overlapping triangles. 
        bool faceProcessed[6]; // To determine whether a face has been triangulated using the FaceDiagonalRules. To avoid overlapping triangles. 
        
    /* OCell does not keep its function values due to the redundancy and 
    extrage storage it incurs.  Instead, every cell keeps its cell
    level and index that can be used to resolve its function values. */
};

#endif	/* OCELL_H */

