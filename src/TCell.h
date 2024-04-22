/* 
 * File:   TCell.h
 * Author: trash001
 *
 * Created on February 24, 2014, 2:43 PM
 */

#ifndef TCELL_H
#define	TCELL_H

#include <string>

class OCell;

/** 
 * This class represents a tetrahedral cell. 
 * An ambiguous OCell (octree cell) is divided into 12 tetrahedral cells.  
 */
class TCell {
    //friend class OctreeContouring;
    
    public:
        TCell(OCell *pc, int lastIndex, int corner1, int corner2, int corner3);
        TCell();
        ~TCell();
        
        void setBPoint(int index);
        int getBPoint();
        OCell* getParent();
        void vtkFile(std::string path); 
        int* getIndex();
        int* getMaterialValues();
        void setMaterialValues(int mat[5]); 
        bool hasSignChangeEdge();
        float* getCenter();
        bool hasEdge(float p[3], float q[3]);
        bool hasVertex(float p[3]);
        bool isSubsetOfFace(float p0[3], float p1[3], float p2[3], float p3[3]);
        float** getCorners();
        bool isFaceDiagonal(float p[3], float q[3], int pq_mat[2]);
        
        
        
    private:
        int m_BiPoint; // BiPoint, i.e. center of the tetrahedra
        OCell *parentOCell;
        float points[5][3]; // The first four points are the four coordinates of the corners of the tetrahedra, the last point is the coordinate of the centroid of the tetrahedra. 
        int index[4];
        int materialVals[4]; // Material values of the vertices of the tetrahedra. The first element represents the material value of the parent OCell's center. 
            
};

#endif	/* TCELL_H */

