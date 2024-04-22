#ifndef _VPT_NEW_ISOCONTOUR_H_
#define _VPT_NEW_ISOCONTOUR_H_

#include <vector> 
#include <iostream>

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <malloc.h>
#include <string.h>
#include <sys/types.h>

#include "Common.h"

#include "Volume.h"

#include "TCell.h"
#include "OCell.h"


#define XPLUS 1
#define YPLUS 2
#define ZPLUS 4

typedef unsigned int uint; 

namespace vpt{

typedef struct _vtx { 
  float x; 
  float y; 
  float z; 
} Vtx; 

//-----------------------------------------------------------------------
//
// cellQueue.h - queue of cell identifiers.  The circular queue dyanmically
//               resizes itself when full.  Elements in the queue are of
//               type and size specified by the user.  Not even template'd
//               yet, just using a void *.
//
// Copyright (c) 1997 Dan Schikore
//------------------------------------------------------------------------
class CellQueue {
   public:
      // constructor/destructor
      inline CellQueue(int size=100);
      ~CellQueue();

      // add item to the queue
      inline void Add(unsigned int cell);

      // remove and return the first item in queue
      inline int  Get(int &cell);

      // return the first item in queue
      inline int  Peek(int &cell);

      // remove the first item in queue
      inline void Pop();

      // reset to empty
      void Reset(void) { nel = 0; }

      // check if queue is empty
      int  Empty(void) { return(nel == 0); }

   protected:

   private:
      int nel;
      int cellsize;  /* # of elements in cell array   */
      int start;
      unsigned int *cells;
};

//------------------------------------------------------------------------
//
// CellQueue() - create a new cell queue with elements of specified size
//
//------------------------------------------------------------------------
inline CellQueue::CellQueue(int size)
{
   nel    = 0;
   start  = 0;
   cellsize = size;
   cells    = (unsigned int *)malloc(sizeof(unsigned int) * cellsize);
}


//------------------------------------------------------------------------
//
// ~CellQueue() - free storage
//
//------------------------------------------------------------------------
inline CellQueue::~CellQueue()
{
   if (cells != NULL)
      free(cells);
}


//------------------------------------------------------------------------
//
// Add() - add an item to the queue
//
//------------------------------------------------------------------------
inline void
CellQueue::Add(unsigned int c)
{
   int n;
   int oldsize;
   int atend;

   n = nel++;

   // resize the queue if needed
   if (nel > cellsize) {
      oldsize = cellsize;
      cellsize *= 2;
      cells = (unsigned int *)realloc(cells, sizeof(int) * cellsize);

      // move everything from 'start' to the end
      if (start != 0) {
         atend = oldsize - start;
         memmove(&cells[cellsize-atend], &cells[start], sizeof(unsigned int)*atend);
         start = cellsize-atend;
      }
   }

   n += start;
   if (n >= cellsize)
      n-=cellsize;

   cells[n] = c;
}

//------------------------------------------------------------------------
//
// Get() - return the top item from the queue
//
//------------------------------------------------------------------------
inline int
CellQueue::Get(int &c)
{
   if (Peek(c) == -1)
      return(-1);
   Pop();
   return(1);
}

//------------------------------------------------------------------------
//
// Peek() - return the top item, but don't remove it
//
//------------------------------------------------------------------------
inline int
CellQueue::Peek(int &c)
{
   if (nel == 0)
      return(-1);

   c = cells[start];

   return(1);
}

//------------------------------------------------------------------------
//
// Pop() - delete the top item in the queue
//
//------------------------------------------------------------------------
inline void
CellQueue::Pop(void)
{
   start++;
   if (start == cellsize)
      start=0;
   nel--;
}
// cellQueue stuff

// Begin bishoulder stuff


#define X_SLICE 0
#define Y_SLICE 1
#define Z_SLICE 2

#define min2(x,y) (x<y) ? (x) : (y)
#define max2(x,y) (x>y) ? (x) : (y)

class BoundaryMap; 
class OctreeContouring; 

enum QuadType  {QUAD_NOT_RENDERED, QUAD_FINE, QUAD_COARSE}; 

/* This class keeps the meta-data (and more) about the input.  
 *  In particular, it stores the size of the input, the 
 *  dimension and scales of the input, etc.
 */
class ContourData{
 public: 
  double decay;  // the rate of decay
  double isovalue;  
  int level;  // level of subdivision
  int gridSize; 
  double* funcVals; // for storage reason, this is the only copy of the function values
  double bdMin[3]; 
  double bdMax[3]; 
  Volume* mask; 
  FVolume* dvol; 
  vval segment;
  BoundaryMap* bmap; 
  /*
  inline getVal(int i){
	  if((*mask)(i)==segment)
		  return funcVals[i]; 
	  return -funcVals[i];
  }
*/

  double widthPerCell[3]; 
  std::vector<double> dimTable; 
};



///* This class represents each cell in the octree */
//class OCell{
//  friend class OctreeContouring; 
// public:
//  OCell(unsigned char level, uint ind[3], OCell* parent); 
//  OCell(unsigned char level, uint ix, uint iy, uint iz, OCell* parent); 
////  OCell(){ 
////	  m_Children = NULL; 
////	  m_Parent = NULL; 
////	  m_BiPoint = -1;
////          tetrahedra = NULL;
////  }
////
////  inline ~OCell(){
////	if(m_Children!=NULL){
////		for(int j=0; j<8; j++){
////			if(m_Children[j] != NULL) 
////				delete m_Children[j]; 
////		}
////		delete [] m_Children; 
////	}
////  }
////	
////  inline void setBPoint(int index){ m_BiPoint = index; }
////  inline int getBPoint(){ return m_BiPoint; }
////  inline bool hasChildren(){ return m_Children!=NULL; }
////  inline void setHasChildren() { m_Children = new OCell*[8]; }
//  
//  OCell();
//  ~OCell();
//	
//  void setBPoint(int index); 
//  int getBPoint(); 
//  bool hasChildren(); 
//  void setHasChildren();
//  
//  void setTetrahedralCells();
//  TCell** getTetrahedralCells();
//  void setIsAmbiguous(bool b);
//  bool getIsAmbiguous();
//  void setHasTetrahedra(bool b);
//  bool getHasTetrahedra();
//  
//  void print(); 
//
//  OCell** m_Children; 
//  OCell* m_Parent; 
//  unsigned char m_Level ; 
//  uint m_Ind[3];
//  
//  float corners[9][3]; // First 8 are the coordinates of cell corners, the last one is the coordinate of the cell center. 
//
// private: 
//  int m_BiPoint;  // an index into the vertex table 	  
//  
//  bool isAmbiguous; // Whether the cell is ambiguous or not. 
//  bool hasTetrahedra; // Whether the tetrahedral cells have been created or not. 
//  TCell **tetrahedra; 
//  
//  /* OCell does not keep its function values due to the redundancy and 
//     extrage storage it incurs.  Instead, every cell keeps its cell
//     level and index that can be used to resolve its function values. */
//};

// builds an adaptive tree given two sets of function values
class CorrectOctree{
  friend class OctreeContouring; 

 public: 
  CorrectOctree();
  ~CorrectOctree(); 

  void testFunction(); 

  bool initialize(double coarseIso, uint coarseLevel, double* coarseMin, double* coarseMax, double* coarseFuncVals, int segment, Volume* mask, FVolume* dvol, BoundaryMap* bmap); 

  bool build(); 
  bool contour();
  void testPrinting1(); 

  uint putMinimizer(ContourData* cd, uint indVec[3], bool clamp); 
  uint putMinimizer(OCell *oc, ContourData* cd, uint indVec[3], bool clamp); 
  uint putBishoulderPoint(ContourData* cd, uint indVec[3]); 

  void printGrid(OCell *oc);    
  
  //void printGridInfo(double edgeLength, int i, int j, int k, double x, double y, double z);
  

  // int corner probably indicates the corner of the octree, i.e. the leaves of the octree. ????? 
  static inline vval getVolumeVal(Volume* mask, int pi, int pj, int pk, int corner){
	  const int d[8][3]={{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}}; 
          //const int d[8][3] = { {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1} }; 
	  _runtimes++; 
	  return (*mask)(pi+d[corner][0],pj+d[corner][1],pk+d[corner][2]);
  }

  static inline float getVolumeVal(FVolume* mask, int pi, int pj, int pk, int corner){
	  const int d[8][3]={{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}}; 
          //const int d[8][3] = { {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1} }; 
	  return (*mask)(pi+d[corner][0],pj+d[corner][1],pk+d[corner][2]);
  }

  MeshGeometry* m_Geom; 
  MeshGeometry* m_GeomColors; 
  OctreeContouring* m_Oc; 

  double m_InnerColor[4]; 
  double m_OuterColor[4];
  ContourData m_CoarseData; 
  
  static int _runtimes; 

  vector<OCell*> _toDel; 
  list<OCell*> _cellq; 

  double _s0; 
  double _s1; 
  double _s2; 
  double _s3; 
  double _s4; 

  OCell* m_Root; 

//  double length; // This is the length of a single edge of a grid cell. 
  
 private: 

  OCell* newOCell(unsigned char level, uint i0, uint i1, uint i2, OCell* parent); 

  void computeNormals();
  void buildTree(OCell* cell, uint level, uint accuInd[3]);
//  OCell* buildTreeBetter(OCell* parent, uint level, uint accuInd[3],bool& same);
  OCell* buildTreeBetter(OCell* parent, uint level, uint ax, uint ay, uint az, bool& same);
  OCell* buildTreeBetterNR(OCell* parent, uint level, uint ax, uint ay, uint az);
  
  //OCell* m_Root; 
};


class OctreeContouring{
 public: 
  OctreeContouring(CorrectOctree* oct); 
  void initializeAndRun(); 

private: 
  void contourTree(OCell* cell, uint cellInd, int level); 
  void contourTree(OCell* cell); 
  void findEdgeValues(uint cellInd, double values[24], int level);        

  void cellProc(OCell* q); 
  void faceProc(OCell* q1, OCell* q2); 
  void edgeProc(OCell* q1, OCell* q2, OCell* q3, OCell* q4); 
  
  /* Given two cells find their relative positions */
  /*  i.e. qB is x+, y+, or z+ in relation to qA */
  int findRelPositionOf2Cells(OCell* qA, OCell* qB); 
  void findEdgeValue(OCell* qs[4], double val[2], vval mat[2]); 
  void faceProcMap1(OCell* qA, OCell* qB, OCell** out8, OCell** out16);
  CorrectOctree* m_Oct; 
  std::vector<double*> m_FuncValArray;
  std::vector<int*> m_SpecialTriangles; // Contains triangles for checking for duplicates. 
  //std::vector<int> m_FaceRuleTriangles; // Contains indices of triangles created by the Face Rule. 
  
  
  
  float** getSharedFace(OCell *c1, OCell *c2);
  float** getSharedFace(TCell *t1, TCell *t2);
  void getSharedEdge(OCell *qs[4], float p[3], float q[3], int pq_material[2]);  
  bool hasSharedFace(OCell *oc1, OCell *oc2);
  bool hasSharedFace(TCell *t1, TCell *t2);
  bool hasSharedFace(OCell *c1, TCell *t1);  

  void materialChangeEdge_WithBackgoundZhang(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[2], float p[3], float q[3], int pq_mat[2]);
  void materialChangeEdge_WithOtherMaterialZhang(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[2], float p[3], float q[3], int pq_mat[2]);
  void interiorEdgeZhang();
  
  void materialChangeEdge_WithBackgound(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[2], float p[3], float q[3], int pq_mat[2]);
  void materialChangeEdge_WithOtherMaterial();
  void interiorEdge();
  
  vector<int> FaceDiagonalRule_Ambiguous_Ambiguous(OCell* q1, OCell* q2, OCell* q3, OCell* q4);
  vector<int> individualAmbiguousCellRule(OCell* q1, OCell* q2, OCell* q3, OCell* q4, float p[3], float q[3], int pq_mat[2]);
  
  double computeRadiusRatio(int pointId1, int pointId2, int pointId3);
  
  vector<int> minimalEdgeTriangulation_WithAmbiguousOCell(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[2], float p[3], float q[3], int pq_mat[2]);
  vector<int> FaceDiagonalRule_Ambiguous_Unambiguous(OCell* q1, OCell* q2, OCell* q3, OCell* q4);
  vector<int> InteriorEdgeRule(OCell* q1, OCell* q2, OCell* q3, OCell* q4);
  void minimalTetraEdge(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[2], float p[3], float q[3], int pq_mat[2]);
  void testPolygonizationFunction(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[2], float p[3], float q[3], int pq_mat[2]); 
  
  // These counters are used to count the number of low quality triangles created in the polygon generation rules. 
  int regularDCPolygonizerCounter;
  int minimalEdgeTriangulation_WithAmbiguousOCellCounter;
  int FaceDiagonalRule_Ambiguous_UnambiguousCounter;
  int FaceDiagonalRule_Ambiguous_AmbiguousCounter;
  int InteriorEdgeRuleCounter;
  
  int processHoles();
  
  bool isFlatish(std::vector<int*> tris, std::string outputFile);
  void cleanup();
  void printTriangles(std::vector<int*> tris, std::string outputFile);
  void printTriangle(int *triangle, std::string outputFile);
  void printOCellGroup(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[], float p[], float q[], int pq_mat[], std::string outputFile, vector<int> t0, vector<int> t1, vector<int> t2, vector<int> t3);
  bool isDuplicateTriangle(int a, int b, int c);
  std::vector<int*> getTriangles(std::vector<int> pointIndex, std::string filename);
  std::vector<int*> getTriangles(int *pointIndex, float **pts, int numberOfPoints, std::string filename);
  std::vector<int*> getCGALTriangles(std::vector<int> pointIndex, std::string filename);
  std::vector<int*> tetrahedralizeAndGetTriangles(std::vector<int> pointIndex, std::string filename);
  
  
}; 


}
#endif
