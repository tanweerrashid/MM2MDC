#include "BoundaryMap.h"
//#include "DualContouring.h"
#include "DualContouring.h"

#include <iostream>
#include <map>
#include <vector>
#include <algorithm> 
#include <ctime>

#include <omp.h>
#include <fstream>

#include "eigen.h"
#include "TCell.h"

#include "OCell.h"

#include <minc2.h>
#include "Dense"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkDelaunay2D.h>
#include <vtkDelaunay3D.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFeatureEdges.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>

//#define USE_TRADITIONAL_DC

using namespace std; 
using namespace vpt;




//#define signChange(a, b, c) ( ((a >= c && b < c) || (a < c && b >=c)) ? true : false)
//#define DEBUG 0
//#define EPS 1E-9
//
//Vtx cube_pos[8] = {
//  {0,0,0} , {1,0,0} , {1,0,1} , {0,0,1} , 
//  {0,1,0} , {1,1,0} , {1,1,1} , {0,1,1} };

int face_edge[4][2] = { {0, 1}, {1, 3}, {3, 2}, {2, 0} };

int CorrectOctree::_runtimes = 0; 

void recurseAndPrint(OCell* c, int level, bool printAllLevel);


//bool compareFloats(float a, float b) {
//    if (fabs(a - b) < 0.0001) {
//        return true;
//    }
//    else return false;
//}

//int computeFaceShoulder(float* face_val, Vtx& face_shoulder, 
//			int slice_direction, float displacement,float isovalue);
//void getSliceValue(float* slice_val , float* val , int slice_direction , 
//		   float displacement);
//
//int face_vtx[6][4] = {
//  {0,1,4,5},
//  {0,1,3,2},
//  {1,5,2,6},
//  {4,5,7,6},
//  {0,4,3,7},
//  {3,2,7,6}
//};
//
//
//float dist2(Vtx a, Vtx b)
//{
//  float x,y; 
//  x= a.x-b.x;
//  y= a.y-b.y;
//  
//  return sqrt(x*x + y*y );
//}
//
//
//// val[8] order
//// val[0] : (0,0,0)
//// val[1] : (1,0,0)
//// val[2] : (1,0,1)
//// val[3] : (0,0,1)
//// val[4] : (0,1,0)
//// val[5] : (1,1,0)
//// val[6] : (1,1,1)
//// val[7] : (0,1,1)
//void GetBishoulder(float* val , Vtx& bishoulder , float isovalue)
//{
//  float face_val[4] , slice_val[4];
//  //Vtx intersection[4];
//  int f,j,e;
//  int edge_idx=0;
//  Vtx face_shoulder1;
//  float min,max;
//  
//  
//  for (f = 0 ; f < 6 ; f++) {
//    
//    for (j = 0 ; j < 4 ; j++) face_val[j] = val[face_vtx[f][j]];
//    
//    for (e = 0 ; e < 4 ; e++) {
//      
//      min=min2(face_val[face_edge[e][0]] , face_val[face_edge[e][1]]);
//      max=max2(face_val[face_edge[e][0]] , face_val[face_edge[e][1]]);
//      
//      if (min < isovalue && max >= isovalue) {
//	
//	edge_idx++;
//	
//      }
//      
//    }
//    if (edge_idx > 0) break;
//    
//  }
//  
//  switch (f) {
//  case 0 :
//    computeFaceShoulder(face_val,face_shoulder1,Z_SLICE,0,isovalue);
//    getSliceValue(slice_val,val,X_SLICE,face_shoulder1.x);
//    computeFaceShoulder(slice_val,bishoulder,X_SLICE,face_shoulder1.x,isovalue);
//    break;
//    
//  case 1 :
//    //	slope=(intersection[1].z-intersection[0].z)/(intersection[1].x-intersection[0].x);
//    
//    computeFaceShoulder(face_val,face_shoulder1,Y_SLICE,0,isovalue);
//    getSliceValue(slice_val,val,X_SLICE,face_shoulder1.x);
//    computeFaceShoulder(slice_val,bishoulder,X_SLICE,face_shoulder1.x,isovalue);
//    break;
//    
//  case 2 :
//    //	slope=(intersection[1].z-intersection[0].z)/(intersection[1].y-intersection[0].y);
//    
//    computeFaceShoulder(face_val,face_shoulder1,X_SLICE,1,isovalue);
//    getSliceValue(slice_val,val,Y_SLICE,face_shoulder1.y);
//    computeFaceShoulder(slice_val,bishoulder,Y_SLICE,face_shoulder1.y,isovalue);
//    break;
//    
//  case 3 :
//    //	slope=(intersection[1].z-intersection[0].z)/(intersection[1].x-intersection[0].x);
//    
//    computeFaceShoulder(face_val,face_shoulder1,Y_SLICE,1,isovalue);
//    getSliceValue(slice_val,val,X_SLICE,face_shoulder1.x);
//    computeFaceShoulder(slice_val,bishoulder,X_SLICE,face_shoulder1.x,isovalue);
//    break;
//    
//  case 4 :
//    //	slope=(intersection[1].z-intersection[0].z)/(intersection[1].y-intersection[0].y);
//    
//    computeFaceShoulder(face_val,face_shoulder1,X_SLICE,0,isovalue);
//    getSliceValue(slice_val,val,Y_SLICE,face_shoulder1.y);
//    computeFaceShoulder(slice_val,bishoulder,Y_SLICE,face_shoulder1.y,isovalue);
//    break;
//    
//  case 5 :
//    //	slope=(intersection[1].y-intersection[0].y)/(intersection[1].x-intersection[0].x);
//    
//    computeFaceShoulder(face_val,face_shoulder1,Z_SLICE,1,isovalue);
//    getSliceValue(slice_val,val,X_SLICE,face_shoulder1.x);
//    computeFaceShoulder(slice_val,bishoulder,X_SLICE,face_shoulder1.x,isovalue);
//    break;
//    
//  }
//}
//
//int computeFaceShoulder(float* face_val, Vtx& face_shoulder, int slice_direction, float displacement,float isovalue)
//{
//  int e , edge_idx;
//  float min,max,ratio;
//  
//  Vtx face_pos[4]={{0,0,0},{1,0,0},{0,1,0},{1,1,0}};
//  Vtx intersection[4];
//  
//  edge_idx=0;	
//  
//  for (e = 0 ; e < 4 ; e++) {
//    
//    min=min2(face_val[face_edge[e][0]] , face_val[face_edge[e][1]]);
//    max=max2(face_val[face_edge[e][0]] , face_val[face_edge[e][1]]);
//    
//    if (min < isovalue && max >= isovalue) {
//      
//      // intersect_edge[edge_idx] = e;
//      ratio = (isovalue - min) / (max - min);
//      
//      if (face_val[face_edge[e][0]] < face_val[face_edge[e][1]]) {
//	
//	intersection[edge_idx].x = face_pos[face_edge[e][0]].x * (1 - ratio) + face_pos[face_edge[e][1]].x * ratio;
//	intersection[edge_idx].y = face_pos[face_edge[e][0]].y * (1 - ratio) + face_pos[face_edge[e][1]].y * ratio;
//	
//      } else {
//	
//	intersection[edge_idx].x = face_pos[face_edge[e][1]].x * (1 - ratio) + face_pos[face_edge[e][0]].x * ratio;
//	intersection[edge_idx].y = face_pos[face_edge[e][1]].y * (1 - ratio) + face_pos[face_edge[e][0]].y * ratio;
//	
//      }
//      edge_idx++;
//      
//    }
//    
//  }
//  
//  if (edge_idx!=2) {
//    face_shoulder.x=0.5;
//    face_shoulder.y=0.5;
//    face_shoulder.z=0.5;
//    return false;
//  }
//  
//  double p,q,r,s;
//  
//  float tface_val[4];
//  tface_val[0]=face_val[0]-isovalue;
//  tface_val[1]=face_val[1]-isovalue;
//  tface_val[2]=face_val[2]-isovalue;
//  tface_val[3]=face_val[3]-isovalue;
//  
//  double slope_sign;
//  
//  
//  p=tface_val[0];
//  q=-tface_val[0]+tface_val[1];
//  r=-tface_val[0]+tface_val[2];
//  s=tface_val[0]-tface_val[1]-tface_val[2]+tface_val[3];
//  
//  double slope;
//  Vtx middle_point,face_candidate1,face_candidate2,face_candidate;
//  
//  slope=(intersection[1].y-intersection[0].y)/(intersection[1].x-intersection[0].x);
//  if (slope>0) slope_sign=1;
//  else slope_sign=-1;
//  
//  middle_point.x = (float)(intersection[0].x + intersection[1].x)/2.f;
//  middle_point.y = (float)(intersection[0].y + intersection[1].y)/2.f;
//  
//  if ((s<0.01 && s>-0.01) || 
//      (slope<0.01 && slope>-0.01) ||
//      (slope<-1000 || slope>1000)) {
//    face_candidate=middle_point;
//  } else {				
//    
//    face_candidate1.x=(float)((-r+sqrt(-(q*r-s*p)/slope))/s);	
//    face_candidate1.y=(float)((-q-slope_sign*sqrt(-(q*r-s*p)*slope))/s);	
//    face_candidate2.x=(float)((-r-sqrt(-(q*r-s*p)/slope))/s);
//    face_candidate2.y=(float)((-q+slope_sign*sqrt(-(q*r-s*p)*slope))/s);
//    
//    
//    if (face_candidate1.x>=0 && face_candidate1.x<=1 && face_candidate1.y>=0 && face_candidate1.y<=1) {
//      face_candidate=face_candidate1;
//    } else if (face_candidate2.x>=0 && face_candidate2.x<=1 && face_candidate2.y>=0 && face_candidate2.y<=1) {
//      face_candidate=face_candidate2;
//    } else {
//      if (dist2(middle_point,face_candidate1)<dist2(middle_point,face_candidate2))
//	face_candidate=face_candidate1;
//      else face_candidate=face_candidate2;
//    }
//  }
//  
//  
//  if (face_candidate.x<0) {
//    face_candidate.x=0;
//    face_candidate=middle_point;
//  }
//  if (face_candidate.x>1) {
//    face_candidate.x=1;
//    face_candidate=middle_point;
//  }
//  if (face_candidate.y<0) {
//    face_candidate.y=0;
//    face_candidate=middle_point;
//  }
//  if (face_candidate.y>1) {
//    face_candidate.y=1;
//    face_candidate=middle_point;
//  }
//  
//  
//  switch (slice_direction) {
//  case X_SLICE :
//    face_shoulder.x=displacement;
//    face_shoulder.y=face_candidate.x;
//    face_shoulder.z=face_candidate.y;
//    break;
//    
//  case Y_SLICE :			
//    face_shoulder.x=face_candidate.x;
//    face_shoulder.y=displacement;
//    face_shoulder.z=face_candidate.y;
//    break;
//    
//  case Z_SLICE :
//    face_shoulder.x=face_candidate.x;
//    face_shoulder.y=face_candidate.y;
//    face_shoulder.z=displacement;
//    break;
//  }
//  
//  return true;
//}
//
//void getSliceValue(float* slice_val , float* val , int slice_direction , float displacement)
//{
//  int x_edge[4][2] = { {0,1},{4,5},{3,2},{7,6} };
//  int y_edge[4][2] = { {0,4},{1,5},{3,7},{2,6} };
//  int z_edge[4][2] = { {0,3},{1,2},{4,7},{5,6} };
//  int i;
//  
//  switch (slice_direction) {
//  case X_SLICE :
//    for (i=0;i<4;i++)
//      slice_val[i] = val[x_edge[i][0]]*(1-displacement) + val[x_edge[i][1]]*displacement;
//    break;
//  case Y_SLICE :
//    for (i=0;i<4;i++)
//      slice_val[i] = val[y_edge[i][0]]*(1-displacement) + val[y_edge[i][1]]*displacement;
//    break;
//  case Z_SLICE :
//    for (i=0;i<4;i++)
//      slice_val[i] = val[z_edge[i][0]]*(1-displacement) + val[z_edge[i][1]]*displacement;
//    break;
//  }
//}
//
//// input : normalized_bishoulder - the location of the computed bishoulder point
////         which is assumed to be in a cell (0,0,0)-(1,1,1)
////         x,y,z - minimum coordinate of a cell which has the bishoulder point
//// output: bishoulder - actual coordinate of the bishoulder point 
//
//void Norm2Read_Bishoulder(Vtx normalized_bishoulder,float* bishoulder,float x,float y,float z,float cell_size)
//{
//  bishoulder[0] = x + cell_size*normalized_bishoulder.x;
//  bishoulder[1] = y + cell_size*normalized_bishoulder.y;
//  bishoulder[2] = z + cell_size*normalized_bishoulder.z;
//}
//
//void Norm2Read_Bishoulder2(Vtx normalized_bishoulder,float* bishoulder, float xyz[3] , float cell_size[3]) 
//{ 
//  bishoulder[0] = xyz[0] + cell_size[0]*normalized_bishoulder.x; 
//  bishoulder[1] = xyz[1] + cell_size[1]*normalized_bishoulder.y; 
//  bishoulder[2] = xyz[2] + cell_size[2]*normalized_bishoulder.z; 
//} 
//
//
//
//bool is_ambiguous_case(float* val,float iso_val,int Case)
//{
//  CellQueue vtx_queue;
//  int vtx_connectivity[8][3]={{1,3,4} , {0,2,5} , {1,3,6} , {0,2,7} , {0,5,7} , {1,4,6} , {2,5,7} , {3,4,6}};
//  int vtx_idx=-1;
//  
//  int code[8] = {0,0,0,0,0,0,0,0} ;
//  int i;
//  for (i = 0 ; i < 8 ; i++) {
//    if (Case==0) {
//      if (val[i] >= iso_val) {
//	code[i] = 1;
//	vtx_idx=i;
//      }
//    } else {
//      if (val[i] < iso_val) {
//	code[i] = 1;
//	vtx_idx=i;
//      }
//    }
//  }
//  
//  if (vtx_idx==-1) return false;
//  else {
//    vtx_queue.Add(vtx_idx);
//    code[vtx_idx]=0;
//    
//    while (vtx_queue.Get(vtx_idx)>=0) {
//      for (i = 0 ; i < 3 ; i++) {
//	if (code[vtx_connectivity[vtx_idx][i]]==1) {
//	  vtx_queue.Add(vtx_connectivity[vtx_idx][i]);
//	  code[vtx_connectivity[vtx_idx][i]]=0;
//	}
//      }
//    }
//  }
//  
//  for (i = 0 ; i < 8 ; i++) {
//    if (code[i]==1) return true;
//  }
//  return false;
//}
//
//// this function decides whether an isosurface in the cell specified by
//// its value is ambiguous or not.
//// if a cell has an ambiguous cell, you must refine the cell until
//// isosurfaces in the cells are not ambiguous
//bool is_ambiguous(float* val,float iso_val)
//{
//  bool flag_val_case1,flag_val_case2;
//  flag_val_case1=is_ambiguous_case(val,iso_val,0);
//  flag_val_case2=is_ambiguous_case(val,iso_val,1);
//  if (flag_val_case1==true || flag_val_case2==true) return true;
//  else return false;
//}

/* the face adjacent children in a cell */
int c_FaceAdj[12][2] = {
    {0,1}, {0,2}, {0,4},
    {1,3}, {1,5},
    {2,3}, {2,6},
    {3,7},
    {4,5}, {5,7},
    {4,6}, 
    {6,7}
}; 

/* the edge adjacent children in a cell */
int c_EdgeAdj[6][4] = {
    {0,4,6,2}, 
    {1,5,7,3}, 
    {0,2,3,1}, 
    {4,6,7,5}, 
    {2,3,7,6},
    {0,1,5,4}
}; 

// the edge and vertex conversion table, each entry is two vertex indices
int EV[12][2] = {
    {0,4}, {4,6}, {4,5}, 
    {1,5}, {5,7}, {0,1},
    {3,7}, {1,3}, {2,3},
    {2,6}, {0,2}, {6,7}
};

//template<typename T>
//void setArray3(T* v, T x, T y, T z){
//	v[0] = x; 
//	v[1] = y;
//	v[2] = z; 
//}
//
template<typename T>
void copyVect(T* src, T* dest, int n) {
	for(int j=0;j<n;j++)
		dest[j] = src[j]; 
}

///************ begin OCell class *************/
//OCell::OCell(unsigned char level, uint ind[3], OCell* parent){
//    int j; 
////    for(j=0; j<8; j++)
//  //      m_Children[j] = NULL; 
//	m_Children = NULL; 
//    m_Level = level; 
//
//    for(j=0; j<3; j++)
//        m_Ind[j] = ind[j]; 
//  
//    m_Parent = parent; 
//    m_BiPoint = -1;     
//}
//
//OCell::OCell(unsigned char level, uint ix, uint iy, uint iz, OCell* parent){
//	m_Children = NULL; 
//    m_Level = level; 
//
//	m_Ind[0] = ix; 	
//	m_Ind[1] = iy; 
//	m_Ind[2] = iz; 
//
//	m_Parent = parent; 
//    m_BiPoint = -1; 
//}
//
//void OCell::print(){
//    int j; 
//    cout<<"Cell: ("<<m_Ind[0]<<","<<m_Ind[1]<<","<<m_Ind[2]<<")"
//        <<" Level: "<<(int)m_Level<<" Parent: "<<m_Parent
//        <<" Has_Children: "<<this->hasChildren()
//        <<" bpoint: "<<m_BiPoint;
//    if(m_Parent != NULL){
//        cout<<" childInd: "; 
//        for(j=0; j<8; j++)
//            if(m_Parent->m_Children[j] == this){
//                cout<<j;
//                break;
//            }  
//    }
//    cout<<endl;
//}
//    
//OCell::OCell() { 
//    m_Children = NULL; 
//    m_Parent = NULL; 
//    m_BiPoint = -1;
//    tetrahedra = NULL;
//}
//
//OCell::~OCell() {
//    if(m_Children!=NULL){
//        for(int j=0; j<8; j++){
//            if(m_Children[j] != NULL) 
//            delete m_Children[j]; 
//        }
//        delete [] m_Children; 
//    }
//}
//
//void OCell::setBPoint(int index) { 
//    m_BiPoint = index; 
//}
//
//int OCell::getBPoint() { 
//    return m_BiPoint; 
//}
//
//bool OCell::hasChildren() { 
//    return m_Children != NULL; 
//}
//
//void OCell::setHasChildren() { 
//    m_Children = new OCell*[8]; 
//}
//
//void OCell::setTetrahedralCells() {
//    tetrahedra = new TCell*[12];
//    tetrahedra[0] = new TCell(0, 0, 1, 3); 
//    tetrahedra[1] = new TCell(1, 0, 2, 3); 
//    tetrahedra[2] = new TCell(2, 0, 1, 5); 
//    tetrahedra[3] = new TCell(3, 0, 4, 5); 
//    tetrahedra[4] = new TCell(4, 0, 2, 6); 
//    tetrahedra[5] = new TCell(5, 0, 4, 6); 
//    tetrahedra[6] = new TCell(6, 2, 3, 7); 
//    tetrahedra[7] = new TCell(7, 2, 6, 7); 
//    tetrahedra[8] = new TCell(8, 1, 3, 7); 
//    tetrahedra[9] = new TCell(9, 1, 5, 7); 
//    tetrahedra[10] = new TCell(10, 4, 6, 7); 
//    tetrahedra[11] = new TCell(11, 4, 5, 7); 
//};
//
//TCell** OCell::getTetrahedralCells() { 
//    return tetrahedra; 
//}
//
//void OCell::setIsAmbiguous(bool b) { 
//    isAmbiguous = b; 
//}
//bool OCell::getIsAmbiguous() { 
//    return isAmbiguous; 
//}
//void OCell::setHasTetrahedra(bool b) { 
//    hasTetrahedra = b; 
//}
//bool OCell::getHasTetrahedra() { 
//    return hasTetrahedra; 
//}
//
//
///************** end OCell class ***************/



/*********** begin CorrectOctree **************/ 
CorrectOctree::CorrectOctree(){
    m_Root = NULL; 

    m_InnerColor[0] = 0.30; 
    m_InnerColor[1] = 0.65; 
    m_InnerColor[2] = 0.20; 
    m_InnerColor[3] = 1.0; 

    m_OuterColor[0] = 0.40; 
    m_OuterColor[1] = 0.65; 
    m_OuterColor[2] = 0.90; 
    m_OuterColor[3] = 1;
	//m_Geom = new DynamicGeometry(); 
	m_Geom = new MeshGeometry(); 
	m_GeomColors = new MeshGeometry(); 
    m_Oc = new OctreeContouring(this); 
}

CorrectOctree::~CorrectOctree(){
    delete m_Root; 
	delete m_Geom; 
	delete m_GeomColors; 
//	delete m_GeomNormals; 
	delete m_Oc; 
}


void CorrectOctree::testFunction(){
   cout<<"\n\n\n**************************\ntest function works\n*****************************\n\n\n"<<endl;
}
	
/* coarseDecay, fineDecay, coarseIso, fineIso are the standard particle model parameters
 * coarseLevel and fineLevel are the level of subdivision for the octree.  So 6 would yield a 64 + 1 sized grid.  
 */
bool CorrectOctree::initialize(double coarseIso, uint coarseLevel, double* coarseMin, double* coarseMax, double* coarseFuncVals, int segment, Volume* mask, FVolume* dvol, BoundaryMap* bmap){
   // string header = "CorrectOctree::initialize";
    int j, k; 

    m_CoarseData.gridSize = (1<<coarseLevel) + 1; 
    m_CoarseData.level = coarseLevel; 
    m_CoarseData.isovalue = coarseIso; 

    copyVect(coarseMin, m_CoarseData.bdMin, 3); 
    copyVect(coarseMax, m_CoarseData.bdMax, 3); 

    m_CoarseData.dimTable.clear(); 

    for(j=0; j<3; j++)
        for(k=0; k<m_CoarseData.gridSize; k++)
            m_CoarseData.dimTable.push_back( (k / (double) m_CoarseData.gridSize) * (coarseMax[j] - coarseMin[j]) + coarseMin[j]); 

    
    
    /*ofstream df;
    df.open("data/dimtable.vtk");
    df << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS " <<  (m_CoarseData.gridSize) << "  float\n\n";

    for(k=0; k<m_CoarseData.gridSize; k++) {
        for(j=0; j<3; j++) {
            df << ((k / (double) m_CoarseData.gridSize) * (coarseMax[j] - coarseMin[j]) + coarseMin[j]) << " "; 
        }
        df << "\n";
    }
    df << "\nVERTICES " << (m_CoarseData.gridSize) << " " << (2 * m_CoarseData.gridSize) << "\n";
    for (k = 0; k < m_CoarseData.gridSize; k++) {
        df << "1 " << k << "\n";
    }
    df << "\n";
    df.close();*/
    
    
    
    for(j=0; j<3; j++)
        m_CoarseData.widthPerCell[j] = (coarseMax[j] - coarseMin[j]) / (m_CoarseData.gridSize -1); 

    m_CoarseData.funcVals = coarseFuncVals; 
	int ngs = m_CoarseData.gridSize -1; 
/*
	m_Geom->resize(ngs*ngs*ngs,ngs*ngs*ngs);
	m_GeomColors->resize(ngs*ngs*ngs,0);
	m_GeomNormals->resize(ngs*ngs*ngs,ngs*ngs*ngs);
*/

	m_CoarseData.segment = (vval) segment; 
	m_CoarseData.mask = mask; 
	m_CoarseData.dvol = dvol; 
	m_CoarseData.bmap = bmap;         
        
//        // Computing the length of an edge of the grid. 
//        // Assuming that all edges are of the same length. 
//        double d = pow(((0 / (double) m_CoarseData.gridSize) * (coarseMax[0] - coarseMin[0]) + coarseMin[0]) - ((1 / (double) m_CoarseData.gridSize) * (coarseMax[0] - coarseMin[0]) + coarseMin[0]), 2) + 
//        pow(((0 / (double) m_CoarseData.gridSize) * (coarseMax[1] - coarseMin[1]) + coarseMin[1]) - ((1 / (double) m_CoarseData.gridSize) * (coarseMax[1] - coarseMin[1]) + coarseMin[1]), 2) + 
//        pow(((0 / (double) m_CoarseData.gridSize) * (coarseMax[2] - coarseMin[2]) + coarseMin[2]) - ((1 / (double) m_CoarseData.gridSize) * (coarseMax[2] - coarseMin[2]) + coarseMin[2]), 2);
//        length = sqrt(d) / sqrt(3);
        
//        cout << "Length of a grid cell's edge: " << length << endl;        
        cout << "Width per cell: {" << m_CoarseData.widthPerCell[0] << ", " << m_CoarseData.widthPerCell[1] << ", " << m_CoarseData.widthPerCell[2] << "}" << endl;
        cout << "Grid Size: " << m_CoarseData.gridSize << endl;
        cout << "Isovalue: " << m_CoarseData.isovalue << endl;
    return true; 
}

bool CorrectOctree::build(){
   // string header="CorrectOctree::build";
//	cout<<"deleting root"<<endl;
    uint ind[3] = {0, 0, 0}; 
    if(m_Root != NULL)
        delete m_Root; 
//	cout<<"done deleting root"<<endl;
  
//	cout<<"building root"<<endl;
	bool same; 
	_runtimes = 0; 
	double start = omp_get_wtime(); 
	_s0 = _s1 = _s2 = _s3 = 0; 
	//m_Root = buildTreeBetter(NULL,0,ind[0],ind[1],ind[2],same); 
	m_Root = buildTreeBetterNR(NULL,0,ind[0],ind[1],ind[2]); 
	cout<<"time: "<<(omp_get_wtime()-start)<<" "<<_runtimes<<endl;
	cout<<_s0<<" "<<_s1<<" "<<_s2<<endl;
        
    return true; 
}

OCell* CorrectOctree::newOCell(unsigned char level, uint i0, uint i1, uint i2, OCell* parent){
	if(_cellq.size()==0){
		for(int j=0;j<1000000;j++){
			_cellq.push_back(new OCell()); 
		}
	}

	OCell* ret = _cellq.front(); 
	ret->m_Level = level; 
	ret->m_Ind[0] = i0; 
	ret->m_Ind[1] = i1; 
	ret->m_Ind[2] = i2; 
	ret->m_Parent = parent; 
        
	_cellq.pop_front(); 
	return ret; 
}

class NRStack{
public: 
	OCell* parent; 
	uint level; 
	uint ax; 
	uint ay; 
	uint az; 
	OCell* lChildren[8]; 
	bool rSame;
	OCell* rCell; 
	int cCPoint; 

	NRStack(OCell* parent, uint level, uint ax, uint ay, uint az){
		this->parent = parent; 
		this->level = level; 
		this->ax = ax; 
		this->ay = ay; 
		this->az = az; 
	}
	NRStack(){
//		parent = NULL; 
	}
}; 

void setNRStack(OCell* p, uint level, uint a, uint b, uint c, NRStack& s){
	s.parent = p; 
	s.level = level; 
	s.ax = a; 
	s.ay = b; 
	s.az = c; 
}

// not recursive version
OCell* CorrectOctree::buildTreeBetterNR(OCell* parent, uint level, uint ax, uint ay, uint az){
	NRStack* lastRet = new NRStack(NULL,0,0,0,0); 
	const int ac[8][3]={{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}}; 
        //const int ac[8][3] = { {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1} }; 
	const int NR_BEGIN = 1; 
	const int NR_END = 2; 
	const int NR_LOOPIND0 = 10; 

	NRStack callStack[20]; 
	int si = 0; 
	setNRStack(parent,level,ax,ay,az,callStack[si]); 
	callStack[si].cCPoint = NR_BEGIN; 

	while(si>=0){
		NRStack* call = &callStack[si]; 

		if(call->cCPoint==NR_BEGIN){
			double start = omp_get_wtime(); 
			call->rCell = new OCell(call->level, call->ax, call->ay, call->az, call->parent, this->m_CoarseData.widthPerCell); 
                        
                        int materialVals[9]; 
                        for (int j = 0; j < 8; j++) {
                            materialVals[j] = (int)getVolumeVal(m_CoarseData.mask, call->ax, call->ay, call->az, j);
                        }
                        materialVals[8] = (*this->m_CoarseData.mask)(call->ax, call->ay, call->az);
                        
                        //cout << "Material value of index [" << call->ax << ", " << call->ay << ", " << call->az << "] " << materialVals[8] << endl;
                        call->rCell->setMaterialValues(materialVals);
                        
                       
			if(call->level < (uint)m_CoarseData.level){
				call->rSame = true; 
				setNRStack(call->rCell,call->level+1,
					call->ax*2+ac[0][0],call->ay*2+ac[0][1],call->az*2+ac[0][2],callStack[si+1]); 
				callStack[si+1].cCPoint = NR_BEGIN; 
				call->cCPoint = NR_LOOPIND0; 
				si++; 
			}
			else{
				vval maskvals[8]; 
				for(int j=0;j<8;j++){
					maskvals[j] = getVolumeVal(m_CoarseData.mask,call->ax,call->ay,call->az,j);
					_runtimes++; 
				}

				bool diffmask = false; 
				for(int j=1;j<8;j++){
					if(maskvals[j]!=maskvals[0]){
						_runtimes++; 
						diffmask = true; 
						break; 
					}
				}
				call->rSame = !diffmask; 
				call->cCPoint = NR_END; 
			}
			_s0+=(omp_get_wtime()-start); 
		}
		else if(call->cCPoint>=NR_LOOPIND0){
			double start = omp_get_wtime(); 
			int lj = call->cCPoint-NR_LOOPIND0; 
			if(lj<8){
				call->lChildren[lj] = lastRet->rCell; 
				if(!lastRet->rSame)
					call->rSame = false; 

				lj++; 
				if(lj<8){
					setNRStack(call->rCell,call->level+1,
						call->ax*2+ac[lj][0],call->ay*2+ac[lj][1],call->az*2+ac[lj][2],callStack[si+1]); 
					callStack[si+1].cCPoint = NR_BEGIN; 
					si++; 
				}
				call->cCPoint = NR_LOOPIND0+lj; 
			}
			else{
				if(!call->rSame){
					call->rCell->setHasChildren(); 
					for(int j=0;j<8;j++){
						call->rCell->m_Children[j] = call->lChildren[j]; 
						_runtimes++; 
					}
				}
				else{
					for(int j=0;j<8;j++){
						delete call->lChildren[j]; 
						_runtimes++; 
					}
				}
				call->cCPoint = NR_END; 
			}
			_s1+=(omp_get_wtime()-start); 
		}
		else{ // call->cCPoint==NR_END
			double start = omp_get_wtime(); 
			lastRet->rCell = call->rCell; 
			lastRet->rSame = call->rSame; 
			si--; 
			_s2+=(omp_get_wtime()-start); 
		}
	}

	OCell* ret = lastRet->rCell; 
	delete lastRet; 

	return ret; 
}


//OCell* CorrectOctree::buildTreeBetter(OCell* parent, uint level, uint aInd[3],bool& same){
OCell* CorrectOctree::buildTreeBetter(OCell* parent, uint level, uint ax, uint ay, uint az, bool& same){

	OCell* ret = newOCell(level, ax, ay, az, parent);			
	const int ac[8][3]={{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}}; 


        // First part
	if(level < (uint) m_CoarseData.level){

		OCell* children[8]; 
		same = true; 
		for(int j=0; j<8; j++){
			bool tsame; 
			children[j] = buildTreeBetter(ret, level+1, ax*2+ac[j][0], ay*2+ac[j][1], az*2+ac[j][2], tsame); 
			if(!tsame)
				same = false; 
		}

                // If the leaves/children have the same values, then they are 
                // deleted from the tree. If they have different values, they
                // are kept in the tree. 
		if(!same){
			ret->setHasChildren(); 
			for(int j=0;j<8;j++)
				ret->m_Children[j] = children[j]; 
		}
		else{
			for(int j=0;j<8;j++){
				delete children[j]; 
				_runtimes++; 
			}
		}
	}
        
        // Second part
	else{
		double start = omp_get_wtime(); 
                
                // Retrieve voxel values, or at least the eight corners of a 
                // voxel ?????
		vval maskvals[8];                 
		for(int j=0;j<8;j++){
			maskvals[j] = getVolumeVal(m_CoarseData.mask,ax,ay,az,j);
			_runtimes++; 
		}

                // Check if the 8 leaves/children have the same values or not. 
		bool diffmask = false; 
		for(int j=1;j<8;j++)
			if(maskvals[j]!=maskvals[0]){
				_runtimes++; 
				diffmask = true; 
				break; 
			}
		same = !diffmask; 
	}

	return ret; 
}

/**
 * Determines whether a point is inside a box. 
 * 
 * @param A - Defines the 8 corners of the box
 * @param p - A point
 * @return - True if point is within the box, false otherwise
 */
bool isPointInsideBox(double A[8][3], float p[3]) {
    bool ret = false;
    double xMax = -DBL_MAX, xMin = DBL_MAX, yMax = -DBL_MAX, yMin = DBL_MAX, zMax = -DBL_MAX, zMin = DBL_MAX;
    for (int i = 0; i < 8; i++) {
        if (A[i][0] > xMax) xMax = A[i][0];
        if (A[i][0] < xMin) xMin = A[i][0];
        
        if (A[i][1] > yMax) yMax = A[i][1];
        if (A[i][1] < yMin) yMin = A[i][1];
        
        if (A[i][2] > zMax) zMax = A[i][2];
        if (A[i][2] < zMin) zMin = A[i][2];
    }
    
    if (p[0] <= xMax && p[0] >= xMin && p[1] <= yMax && p[1] >= yMin && p[2] <= zMax && p[2] >= zMin) {
        ret = true;
    }
    return ret;
}

/**
 * Determines whether a point is inside a box. 
 * 
 * @param A - Defines the 8 corners of the box
 * @param p - A point
 * @return - True if point is within the box, false otherwise
 */
bool isPointInsideBox(double A[8][3], double p[3]) {
    bool ret = false;
    double xMax = -DBL_MAX, xMin = DBL_MAX, yMax = -DBL_MAX, yMin = DBL_MAX, zMax = -DBL_MAX, zMin = DBL_MAX;
    for (int i = 0; i < 8; i++) {
        if (A[i][0] > xMax) xMax = A[i][0];
        if (A[i][0] < xMin) xMin = A[i][0];
        
        if (A[i][1] > yMax) yMax = A[i][1];
        if (A[i][1] < yMin) yMin = A[i][1];
        
        if (A[i][2] > zMax) zMax = A[i][2];
        if (A[i][2] < zMin) zMin = A[i][2];
    }
    
    if (p[0] <= xMax && p[0] >= xMin && p[1] <= yMax && p[1] >= yMin && p[2] <= zMax && p[2] >= zMin) {
        ret = true;
    }
    return ret;
}

void printGridInfo(double wpc[3], int i, int j, int k, double x, double y, double z, int pointNo) {
    double A[3] = {(i * wpc[0]),           (j * wpc[1]),           (k * wpc[2])};
    double B[3] = {(i * wpc[0]),           (j * wpc[1]),           (k * wpc[2]) + wpc[2]};
    double C[3] = {(i * wpc[0]),           (j * wpc[1]) + wpc[1],  (k * wpc[2])};
    double D[3] = {(i * wpc[0]),           (j * wpc[1]) + wpc[1],  (k * wpc[2]) + wpc[2]};
    double E[3] = {(i * wpc[0]) + wpc[0],  (j * wpc[1]),           (k * wpc[2])};
    double F[3] = {(i * wpc[0]) + wpc[0],  (j * wpc[1]),           (k * wpc[2]) + wpc[2]};
    double G[3] = {(i * wpc[0]) + wpc[0],  (j * wpc[1]) + wpc[1],  (k * wpc[2])};
    double H[3] = {(i * wpc[0]) + wpc[0],  (j * wpc[1]) + wpc[1],  (k * wpc[2]) + wpc[2]};
    
    
//    cout << "{" << x << ", " << y << ", " << z << "}, {";
//    cout << A[0] << ", " << A[1] << ", " << A[2] << "}, {";
//    cout << B[0] << ", " << B[1] << ", " << B[2] << "}, {";
//    cout << C[0] << ", " << C[1] << ", " << C[2] << "}, {";
//    cout << D[0] << ", " << D[1] << ", " << D[2] << "}, {";
//    cout << E[0] << ", " << E[1] << ", " << E[2] << "}, {";
//    cout << F[0] << ", " << F[1] << ", " << F[2] << "}, {";
//    cout << G[0] << ", " << G[1] << ", " << G[2] << "}, {";
//    cout << H[0] << ", " << H[1] << ", " << H[2] << "}, {";
//    cout << "\n";
    
    double grid[][3] = {
        {A[0], A[1], A[2]}, 
        {B[0], B[1], B[2]}, 
        {C[0], C[1], C[2]}, 
        {D[0], D[1], D[2]}, 
        {E[0], E[1], E[2]}, 
        {F[0], F[1], F[2]}, 
        {G[0], G[1], G[2]}, 
        {H[0], H[1], H[2]}
    };
    
    double p[] = {x, y, z};
    
    //if (isPointInsideBox(grid, p) == false) {    
        ofstream f("data/gridInfo.txt", ios::app);
        f << pointNo << " " << x << " " << y << " " << z << " ";
        f << A[0] << " " << A[1] << " " << A[2] << " ";
        f << B[0] << " " << B[1] << " " << B[2] << " ";
        f << C[0] << " " << C[1] << " " << C[2] << " ";
        f << D[0] << " " << D[1] << " " << D[2] << " ";
        f << E[0] << " " << E[1] << " " << E[2] << " ";
        f << F[0] << " " << F[1] << " " << F[2] << " ";
        f << G[0] << " " << G[1] << " " << G[2] << " ";
        f << H[0] << " " << H[1] << " " << H[2] << " ";
        f << "\n";    
        f.close();
    //}
}

void CorrectOctree::printGrid(OCell *oc) {
    if (oc != NULL) {
        if (oc->hasChildren()) {
            for (int i = 0; i < 8; i++) {
                printGrid(oc->m_Children[i]);
            }
        }
        else {
            if ((int)oc->m_Level == m_CoarseData.level) {
                if (oc->getBPoint() != -1) {
                    printGridInfo(this->m_CoarseData.widthPerCell, (int)oc->m_Ind[0], (int)oc->m_Ind[1], (int)oc->m_Ind[2], this->m_Geom->getVert(oc->getBPoint())[0], this->m_Geom->getVert(oc->getBPoint())[1], this->m_Geom->getVert(oc->getBPoint())[2], oc->getBPoint());
                }
            }
        }
    }
}

bool CorrectOctree::contour(){
    //  string header = "CorrectOctree::contour()"; 
    m_Oc->initializeAndRun(); 

    // build normals
    //	computeNormals();
    m_Geom->computeNormalsQ(); 	
    
    
    /***********************************************/
    
    //recurseAndPrint(m_Root, 5, false);
    //m_CoarseData.dimTable.
    //cout << "\n\n\nPrinting Octree\n\n\n" << endl;
    //printOCellInfo(m_Root);
    //OCell **rootChildren = m_Root->m_Children;
    
    cout << "Max Level/Octree Depth: " << m_CoarseData.level << endl;
    cout << "Grid size: " << m_CoarseData.gridSize << endl;
    
    //traverseTree(m_Root);
    
    
    
    
    
    
    /***********************************************/
        
        
    return true; 
}

void CorrectOctree::computeNormals(){}

/* 
 * cell : the current cell of subdivision
 * accuInd : accumulates the index from above, 
 * level : the current level of subdivision 
 * flag : for passing on boolean flags   
 */
void CorrectOctree::buildTree(OCell* cell, uint level, uint aInd[3]){
    int j, k;  
    uint aCopy[3]; 

	if(level < (uint) m_CoarseData.level){
		cell->setHasChildren();
		for(j=0; j<8; j++){

			for(k=0; k<3; k++)
				aCopy[k] = aInd[k] * 2; 
      
			if(j & ZPLUS) aCopy[2]++; 
			if(j & YPLUS) aCopy[1]++; 
			if(j & XPLUS) aCopy[0]++; 
      
			cell->m_Children[j] = new OCell(level+1, aCopy, cell, this->m_CoarseData.widthPerCell);
			buildTree(cell->m_Children[j], level+1, aCopy); 
		}
    }
    else{
		//cell->setBPoint(putMinimizer(&m_CoarseData, cell->m_Ind)); 
    }
}

int edges[12][2]={
	{0,1},{0,2},{0,4},
	{1,3},{1,5},
	{2,3},{2,6},
	{3,7},
	{4,5},{4,6},
	{5,7},
	{6,7}
};
	

float corners[8][3]={
	{0,0,0},
	{1,0,0},
	{0,1,0},
	{1,1,0},
	{0,0,1},
	{1,0,1},
	{0,1,1},
	{1,1,1}
};


void computeNormal(vval maskvals[8], float dvals[8], vval m1, vval m2, const osg::Vec3& pos, float* out){
	float v[8];
	for(int j=0;j<8;j++){
		if(maskvals[j]==m1)
			v[j] = dvals[j];
		else if(maskvals[j]==m2)
			v[j] = -dvals[j];
		else
			v[j] = 0; 
	}

	float s = pos.x(); 
	float t = pos.y(); 
	float r = pos.z(); 

	float nx = 
		- (1-r)*(1-t)*v[0] - r*(1-t)*v[4] 
	    - (1-r)*t*v[2] - r*t*v[6] 
        + (1-r)*(1-t)*v[1] + r*(1-t)*v[5] 
        + (1-r)*t*v[3] + r*t*v[7];

	float ny =  
		- (1-r)*(1-s)*v[0] - r*(1-s)*v[4] 
	    + (1-r)*(1-s)*v[2] + r*(1-s)*v[6] 
        - (1-r)*s*v[1] - r*s*v[5] 
        + (1-r)*s*v[3] + r*s*v[7]; 

	float nz =  
        - (1-s)*(1-t)*v[0] + (1-s)*(1-t)*v[4] 
        - (1-s)*t*v[2] + (1-s)*t*v[6] 
        - s*(1-t)*v[1] + s*(1-t)*v[5] 
        - s*t*v[3] + s*t*v[7]; 

	float norm = sqrt(nx*nx + ny*ny + nz*nz);

	out[0] = nx/norm;
	out[1] = ny/norm;
	out[2] = nz/norm;
}

//void computeSubdivisionNormal(float a[8], float b[8], const osg::Vec3& pos, float* out){
//	float s = pos.x(); 
//	float t = pos.y(); 
//	float r = pos.z(); 
//
//	float nx = a[1]+b[0]-b[1] 
//        - (a[1]+a[4]-a[5]+b[0]-b[1]-b[4]+b[5])*r 
//		+ (-a[1]-a[2]+a[3]-b[0]+b[1]+b[2]-b[3] + 
//		(a[1]+a[2]-a[3]+a[4]-a[5]-a[6]+a[7]+b[0]-b[1]-b[2]+b[3]-b[4]+b[5]+b[6]-b[7]) * r) * t 
//		+ a[0]*(-1+r+t-r*t);
//
//	float ny = a[2]+b[0]-b[2] 
//	    - (a[2]+a[4]-a[6]+b[0]-b[2]-b[4]+b[6])*r 
//        + (-a[1]-a[2]+a[3]-b[0]+b[1]+b[2]-b[3] + 
//		(a[1]+a[2]-a[3]+a[4]-a[5]-a[6]+a[7]+b[0]-b[1]-b[2]+b[3]-b[4]+b[5]+b[6]-b[7]) * r) * s 
//		+ a[0]*(-1+r+s-r*s);
//
//	float nz = a[4]+b[0]-b[4]
//	    - (a[1]+a[4]-a[5]+b[0]-b[1]-b[4]+b[5])*s 
//		+ (-a[2]-a[4]+a[6]-b[0]+b[2]+b[4]-b[6] + 
//		(a[1]+a[2]-a[3]+a[4]-a[5]-a[6]+a[7]+b[0]-b[1]-b[2]+b[3]-b[4]+b[5]+b[6]-b[7]) * s) * t 
//		+ a[0]*(-1+s+t-s*t);
//
//	float norm = sqrt(nx*nx + ny*ny + nz*nz);
//
//	out[0] = nx/norm;
//	out[1] = ny/norm;
//	out[2] = nz/norm;
//}

void qefMinimizer(ContourData* cd, int pi, int pj, int pk, float vert[3], float interceptAvg[3]){
    Volume* mask = cd->mask; 
    FVolume* dvol = cd->dvol; 
    int dirx = 1<<0; 
    int diry = 1<<1; 
    int dirz = 1<<2; 

    vval maskvals[8]; 
    float dvals[8]; 

    // Retrieve volume values at specified corner j. 
    // maskvals[] contains corner sign values. 
    for(int j=0;j<8;j++){
        maskvals[j] = CorrectOctree::getVolumeVal(mask,pi,pj,pk,j);
        dvals[j] = CorrectOctree::getVolumeVal(dvol,pi,pj,pk,j);
    }

    float qrmat[10];
    for(int j=0;j<10;j++) qrmat[j]=0;

    float cellPt[3] = {0.f,0.f,0.f};
    int count=0;

    for(int eg=0;eg<12;eg++){
        int e1 = edges[eg][0]; 
        int e2 = edges[eg][1]; 
        //int dir = e2-e1; 
        vval m1=0,m2=0;
        float v1=0,v2=0;

        m1 = maskvals[e1]; 
        m2 = maskvals[e2];

        //if ((int)m1 == (int)m2) {
        //   cout << "\tSame Corner (" << (int)m1 << " and " << (int)m2 << ")" << endl;
        //}
        
        if(m1!=m2){
            v1 = dvals[e1];
            v2 = dvals[e2];

            //v2 = -v2; 
            float frac = v1/(v1+v2);
            osg::Vec3 p1(corners[e1][0],corners[e1][1],corners[e1][2]);
            osg::Vec3 p2(corners[e2][0],corners[e2][1],corners[e2][2]);

            // intpt is intercept point (probably???)
            osg::Vec3 intpt = p1 * (1-frac) + p2 * frac; 

            // cellPt is minimizer
            cellPt[0] += intpt.x();
            cellPt[1] += intpt.y();
            cellPt[2] += intpt.z();
            count++;

            float eq[4]; 
            computeNormal(maskvals,dvals,m1,m2,intpt,eq);

            eq[3] = eq[0]*intpt.x() + eq[1]*intpt.y() + eq[2]*intpt.z();            
            qr(qrmat,eq);
        }
    }
    
    cellPt[0]/=count;
    cellPt[1]/=count;
    cellPt[2]/=count;
    //float res[3];

    interceptAvg[0] = cellPt[0];
    interceptAvg[1] = cellPt[1];
    interceptAvg[2] = cellPt[2];
    
    float error = calcPoint(cellPt,vert,qrmat);

    // Tao Ju's CLAMPING
    /*double A[3] = {(pi * wpc[0]), (pj * wpc[1]), (pk * wpc[2])};
    double B[3] = {(pi * wpc[0]) + wpc[0], (pj * wpc[1]), (pk * wpc[2])};
    double C[3] = {(pi * wpc[0]) + wpc[0], (pj * wpc[1]), (pk * wpc[2]) + wpc[2]};
    double D[3] = {(pi * wpc[0]), (pj * wpc[1]), (pk * wpc[2]) + wpc[2]};

    double E[3] = {(pi * wpc[0]), (pj * wpc[1]) + wpc[1], (pk * wpc[2])};
    double F[3] = {(pi * wpc[0]) + wpc[0], (pj * wpc[1]) + wpc[1], (pk * wpc[2])};
    double G[3] = {(pi * wpc[0]) + wpc[0], (pj * wpc[1]) + wpc[1], (pk * wpc[2]) + wpc[2]};
    double H[3] = {(pi * wpc[0]), (pj * wpc[1]) + wpc[1], (pk * wpc[2]) + wpc[2]};    

    double grid[][3] = {
        {A[0], A[1], A[2]}, 
        {B[0], B[1], B[2]}, 
        {C[0], C[1], C[2]}, 
        {D[0], D[1], D[2]}, 
        {E[0], E[1], E[2]}, 
        {F[0], F[1], F[2]}, 
        {G[0], G[1], G[2]}, 
        {H[0], H[1], H[2]}
    };

    if (isPointInsideBox(grid, vert) == false) {
        vert[0] = cellPt[0];
        vert[1] = cellPt[1];
        vert[2] = cellPt[2];
    }*/
    
    //return error;
}

//void qefSubdivisionMinimizer(ContourData* cd, int pi, int pj, int pk, float vert[3]){
//    cout << "qefSubdivisionMinimizer" << endl;
//	Volume* mask = cd->mask; 
//	FVolume* dvol = cd->dvol; 
//	BoundaryMap* bm = cd->bmap; 
//	int dirx = 1<<0; 
//	int diry = 1<<1; 
//	int dirz = 1<<2; 
//
//	vval maskvals[8]; 
//	float dvals[8][8]; 
//	int dmap[NUM_SEGMENTS]; 
//	for(int j=0;j<NUM_SEGMENTS;j++) dmap[j] = -1; 
//	int numMats = 0; 
//
//	for(int j=0;j<8;j++){
//		maskvals[j] = CorrectOctree::getVolumeVal(mask,pi,pj,pk,j);
//		int mat = maskvals[j]; 
//		if(dmap[mat]<0){
//			//cout<<numMats<<": ";
//			//cout<<"mat: "<<mat<<endl;
//			for(int k=0;k<8;k++){
//				dvals[numMats][k] = CorrectOctree::getVolumeVal(bm->getFunction(mat),pi,pj,pk,k);
//				//cout<<dvals[numMats][k]<<" "; 
//			}
//			//cout<<endl;
//			dmap[mat] = numMats++; 
//		}
//	}
//
//	float qrmat[10];
//	for(int j=0;j<10;j++) qrmat[j]=0;
//	float cellPt[3] = {0.f,0.f,0.f};
//	int count=0;
//
//	for(int eg=0;eg<12;eg++){
//		int e1 = edges[eg][0]; 
//		int e2 = edges[eg][1]; 
//		vval m1 = maskvals[e1]; 
//		vval m2 = maskvals[e2];
//		if(m1!=m2){
//			int d1 = dmap[m1]; 
//			int d2 = dmap[m2]; 
//			float a0 = dvals[d1][e1]; 
//			float a1 = dvals[d1][e2]; 
//			float b0 = dvals[d2][e1]; 
//			float b1 = dvals[d2][e2]; 
//
//			float denom = a0-a1-b0+b1; 
//
//			if(denom>1e-6 || denom<-1e-6){
//				float frac = (a0-b0)/denom; 
//				//cout<<a0<<" "<<a1<<" "<<b0<<" "<<b1<<" "<<frac<<endl;
//				//cout<<"f: "<<frac<<endl;
//				//if(frac>0)
//				//cout<<"f: "<<(a0-b0)<<"/"<<denom<<" "<<frac<<endl;
//				if(frac>=0 && frac<=1){
//					osg::Vec3 p1(corners[e1][0],corners[e1][1],corners[e1][2]);
//					osg::Vec3 p2(corners[e2][0],corners[e2][1],corners[e2][2]);
//
//					osg::Vec3 intpt = p1 * (1-frac) + p2 * frac; 
//
//					cellPt[0] += intpt.x();
//					cellPt[1] += intpt.y();
//					cellPt[2] += intpt.z();
//					count++;
//
//					float eq[4]; 
//					computeSubdivisionNormal(dvals[d1],dvals[d2],intpt,eq);
//					eq[3] = eq[0]*intpt.x() + eq[1]*intpt.y() + eq[2]*intpt.z();
//					qr(qrmat,eq);
//				}
//			}
//		}
//	}
//	if(count>0){
//		cellPt[0]/=count;
//		cellPt[1]/=count;
//		cellPt[2]/=count;
//	}
//	else{   
//		cellPt[0] = .5f; 
//		cellPt[1] = .5f; 
//		cellPt[2] = .5f; 
//	}
//
//	float error = calcPoint(cellPt,vert,qrmat);
//	//cout<<"vert: "<<vert[0]<<" "<<vert[1]<<" "<<vert[2]<<endl;
//        //return error;
//}

//uint CorrectOctree::putBishoulderPoint(ContourData* cd, uint indVec[3]){
//    uint ind; 
//    float tmpVal[8]; 
//    uint plus_x, plus_y, plus_z; 
//    bool allGreater, allLess; 
//    uint grid_size; 
//    double bivert[3]; 
//
//    double tmp[3]; 
//    Vtx bPoint; 
//    int j; 
//  
//    grid_size = cd->gridSize; 
//    ind = (indVec[2] * grid_size * grid_size) + (indVec[1] * grid_size) + indVec[0]; 
//  
//    plus_x = 1; 
//    plus_y = grid_size; 
//    plus_z = grid_size * grid_size; 
//
//    tmpVal[0] = (float) cd->funcVals[ind];
//
//    ind += plus_x;   // 100
//    tmpVal[1] = (float) cd->funcVals[ind];
//
//    ind += plus_z;   // 101
//    tmpVal[2] = (float) cd->funcVals[ind];
//
//    ind -= plus_x;   // 001
//    tmpVal[3] = (float) cd->funcVals[ind];
//
//    ind = ind - plus_z + plus_y;   // 010
//    tmpVal[4] = (float) cd->funcVals[ind];
//
//    ind += plus_x;   // 110
//    tmpVal[5] = (float) cd->funcVals[ind];
//
//    ind += plus_z;   // 111
//    tmpVal[6] = (float) cd->funcVals[ind];
//
//    ind -= plus_x;   // 011
//    tmpVal[7] = (float) cd->funcVals[ind];
//
//	int best = 0; 
//    allGreater = allLess = true; 
//  
//    for(j=0; j<8; j++){
//        if(tmpVal[j] > cd->isovalue)
//            allLess = false; 
//        if(tmpVal[j] < cd->isovalue)
//            allGreater = false;           
//    }
//
//	float cellLower[3];   
//
//    for(int j=0; j<3; j++)
//        cellLower[j] = cd->dimTable[indVec[j]+(j * cd->gridSize)]; 
//
//
//    if((allLess || allGreater))
//        for(j=0;j<3; j++)
//            bivert[j] = cellLower[j] + tmp[j]; 
//    else{
//        GetBishoulder(tmpVal, bPoint, cd->isovalue); 
//        bivert[0] = ((double) bPoint.x) * cd->widthPerCell[0] + cellLower[0] ; 
//        bivert[1] = ((double) bPoint.y) * cd->widthPerCell[1] + cellLower[1] ; 
//        bivert[2] = ((double) bPoint.z) * cd->widthPerCell[2] + cellLower[2] ; 
//    }    
//    
//    return m_Geom->addVertex(bivert[0],bivert[1],bivert[2]);;
//}


/**
 * Outputs the cube as a VTK file. 
 * 
 * @param path - Location of the file to be written. 
 * @param corners - The corners of the cube. 
 * @param materials - The material values of the cube's corners
 * @param p - The center or minimizer of the cube. 
 */
void printCube(std::string path, float **corners, int materials[8], float p[3]) {
    ofstream file;
    file.open(path.c_str());

    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output PointsOutsideGridCube\n";
    file << "ASCII\nDATASET POLYDATA\n";
    file << "POINTS " << 9 << " float\n";
    for (int i = 0; i < 8; i++) {        
        file << corners[i][0] << " " << corners[i][1] << " " << corners[i][2] << "\n";
    }
    file << p[0] << " " << p[1] << " " << p[2] << "\n";

    file << "LINES " << 6 << " " << (6 * 6) << "\n";
    file << "5 0 1 3 2 0\n";
    file << "5 0 4 6 2 0\n";
    file << "5 4 5 7 6 4\n";
    file << "5 1 3 7 5 1\n";
    file << "5 0 1 5 4 0\n";
    file << "5 2 3 7 6 2\n";

    file << "VERTICES " << 1 << " " << (2) << "\n";
    file << "1 8\n";

    file << "POINT_DATA " << 9 << "\nSCALARS scalars float\nLOOKUP_TABLE default\n";
    file << materials[0] << "\n";
    file << materials[1] << "\n";
    file << materials[2] << "\n";
    file << materials[3] << "\n";
    file << materials[4] << "\n";
    file << materials[5] << "\n";
    file << materials[6] << "\n";
    file << materials[7] << "\n";
    file << "10\n" << endl;
    
    file.close();
}

/**
 * Adjusts the position of the minimizer if the minimizer is outside its cube. 
 * 
 * @param p - Minimizer coordinates
 * @param A - Minimizer's cube. 
 * @param wpc - Lengths of the cube's edges. 
 */
void correctPoint(float p[3], double A[8][3], double wpc[3]) {
    while (isPointInsideBox(A, p) == false) {
        double xMax = -DBL_MAX, xMin = DBL_MAX, yMax = -DBL_MAX, yMin = DBL_MAX, zMax = -DBL_MAX, zMin = DBL_MAX;
        for (int i = 0; i < 8; i++) {
            if (A[i][0] > xMax) xMax = A[i][0];
            if (A[i][0] < xMin) xMin = A[i][0];

            if (A[i][1] > yMax) yMax = A[i][1];
            if (A[i][1] < yMin) yMin = A[i][1];

            if (A[i][2] > zMax) zMax = A[i][2];
            if (A[i][2] < zMin) zMin = A[i][2];
        }
        
        if (p[0] > xMax) p[0] = p[0] - (wpc[0] / 100); 
        if (p[0] < xMin) p[0] = p[0] + (wpc[0] / 100);
        
        if (p[1] > yMax) p[1] = p[1] - (wpc[1] / 100); 
        if (p[1] < yMin) p[1] = p[1] + (wpc[1] / 100);
        
        if (p[2] > zMax) p[2] = p[2] - (wpc[2] / 100); 
        if (p[2] < zMin) p[2] = p[2] + (wpc[2] / 100);
        
        //correctPoint(p, A, wpc);
    }
}

uint CorrectOctree::putMinimizer(ContourData* cd, uint indVec[3], bool clamp){
    float cellLower[3];   

    for(int j=0; j<3; j++) {
        cellLower[j] = cd->dimTable[indVec[j]+(j * cd->gridSize)];
    }
    
    //if (((int)indVec[0] != cellLower[0]) || ((int)indVec[1] != cellLower[1]) || ((int)indVec[2] != cellLower[2])) {
    //    cout << "Mismatch indVec [" << indVec[0] << ", " << indVec[1] << ", " << indVec[2] << "] and cellLower [" << cellLower[0] << ", " << cellLower[1] << ", " << cellLower[2] << "]" << endl;
    //}
    
    float bivert[3], interceptAvg[3], point[3]; // bivert is the QEF computed point, interceptAvg is the average of the intercepts. Point is the final value. 
    if(cd->bmap == NULL)
            qefMinimizer(cd, indVec[0], indVec[1], indVec[2], bivert, interceptAvg);
    //else
            //qefSubdivisionMinimizer(cd,indVec[0],indVec[1],indVec[2],bivert);

    point[0] = (bivert[0]) * cd->widthPerCell[0] + cellLower[0]; 
    point[1] = (bivert[1]) * cd->widthPerCell[1] + cellLower[1]; 
    point[2] = (bivert[2]) * cd->widthPerCell[2] + cellLower[2]; 
    
    // Tao Ju's CLAMPING
    if (clamp == true) {
        int pi = (int)indVec[0]; 
        int pj = (int)indVec[1]; 
        int pk = (int)indVec[2];
        
        double wpc[3] = {this->m_CoarseData.widthPerCell[0], this->m_CoarseData.widthPerCell[1], this->m_CoarseData.widthPerCell[2]};

        double A[3], B[3], C[3], D[3], E[3], F[3], G[3], H[3];
        A[0] = (pi * wpc[0]);                A[1] = (pj * wpc[1]);            A[2] = (pk * wpc[2]);
        B[0] = (pi * wpc[0]) + wpc[0];       B[1] = (pj * wpc[1]);            B[2] = (pk * wpc[2]);
        C[0] = (pi * wpc[0]);                C[1] = (pj * wpc[1]) + wpc[1];   C[2] = (pk * wpc[2]);
        D[0] = (pi * wpc[0]) + wpc[0];       D[1] = (pj * wpc[1]) + wpc[1];   D[2] = (pk * wpc[2]);
        E[0] = (pi * wpc[0]);                E[1] = (pj * wpc[1]);            E[2] = (pk * wpc[2]) + wpc[2];
        F[0] = (pi * wpc[0]) + wpc[0];       F[1] = (pj * wpc[1]);            F[2] = (pk * wpc[2]) + wpc[2];
        G[0] = (pi * wpc[0]);                G[1] = (pj * wpc[1]) + wpc[1];   G[2] = (pk * wpc[2]) + wpc[2];
        H[0] = (pi * wpc[0]) + wpc[0];       H[1] = (pj * wpc[1]) + wpc[1];   H[2] = (pk * wpc[2]) + wpc[2];
//        double A[3] = {(pi * wpc[0]), (pj * wpc[1]), (pk * wpc[2])};
//        double B[3] = {(pi * wpc[0]) + wpc[0], (pj * wpc[1]), (pk * wpc[2])};
//        double C[3] = {(pi * wpc[0]) + wpc[0], (pj * wpc[1]), (pk * wpc[2]) + wpc[2]};
//        double D[3] = {(pi * wpc[0]), (pj * wpc[1]), (pk * wpc[2]) + wpc[2]};
//
//        double E[3] = {(pi * wpc[0]), (pj * wpc[1]) + wpc[1], (pk * wpc[2])};
//        double F[3] = {(pi * wpc[0]) + wpc[0], (pj * wpc[1]) + wpc[1], (pk * wpc[2])};
//        double G[3] = {(pi * wpc[0]) + wpc[0], (pj * wpc[1]) + wpc[1], (pk * wpc[2]) + wpc[2]};
//        double H[3] = {(pi * wpc[0]), (pj * wpc[1]) + wpc[1], (pk * wpc[2]) + wpc[2]};    

        double grid[][3] = {
            {A[0], A[1], A[2]}, 
            {B[0], B[1], B[2]}, 
            {C[0], C[1], C[2]}, 
            {D[0], D[1], D[2]}, 
            {E[0], E[1], E[2]}, 
            {F[0], F[1], F[2]}, 
            {G[0], G[1], G[2]}, 
            {H[0], H[1], H[2]}
        };

        if (isPointInsideBox(grid, point) == false) {
            point[0] = (interceptAvg[0]) * cd->widthPerCell[0] + cellLower[0]; 
            point[1] = (interceptAvg[1]) * cd->widthPerCell[1] + cellLower[1]; 
            point[2] = (interceptAvg[2]) * cd->widthPerCell[2] + cellLower[2];
        }
    }
    
    int pt = m_Geom->addVertex(point[0], point[1], point[2]);
    //if (pt == 776 || pt == 825 || pt == 829) cout << point[0] << ", " << point[1] << ", " << point[2] << endl;
    return pt;
}

uint CorrectOctree::putMinimizer(OCell *oc, ContourData* cd, uint indVec[3], bool clamp){
    float cellLower[3];   

    for(int j=0; j<3; j++) {
        cellLower[j] = cd->dimTable[indVec[j]+(j * cd->gridSize)];
    }
    
    //if (((int)indVec[0] != cellLower[0]) || ((int)indVec[1] != cellLower[1]) || ((int)indVec[2] != cellLower[2])) {
    //    cout << "Mismatch indVec [" << indVec[0] << ", " << indVec[1] << ", " << indVec[2] << "] and cellLower [" << cellLower[0] << ", " << cellLower[1] << ", " << cellLower[2] << "]" << endl;
    //}
    
    float bivert[3], interceptAvg[3], point[3]; // bivert is the QEF computed point, interceptAvg is the average of the intercepts. Point is the final value. 
    if(cd->bmap == NULL)
            qefMinimizer(cd, indVec[0], indVec[1], indVec[2], bivert, interceptAvg);
    //else
            //qefSubdivisionMinimizer(cd,indVec[0],indVec[1],indVec[2],bivert);

    point[0] = (bivert[0]) * cd->widthPerCell[0] + cellLower[0]; 
    point[1] = (bivert[1]) * cd->widthPerCell[1] + cellLower[1]; 
    point[2] = (bivert[2]) * cd->widthPerCell[2] + cellLower[2]; 
    
    // Tao Ju's CLAMPING
    if (clamp == true) {
        float **corners = oc->getCorners();
        int *materials = oc->getMaterialValues();

        double grid[8][3];
        grid[0][0] = corners[0][0]; grid[0][1] = corners[0][1]; grid[0][2] = corners[0][2];
        grid[1][0] = corners[1][0]; grid[1][1] = corners[1][1]; grid[1][2] = corners[1][2];
        grid[2][0] = corners[2][0]; grid[2][1] = corners[2][1]; grid[2][2] = corners[2][2];
        grid[3][0] = corners[3][0]; grid[3][1] = corners[3][1]; grid[3][2] = corners[3][2];
        grid[4][0] = corners[4][0]; grid[4][1] = corners[4][1]; grid[4][2] = corners[4][2];
        grid[5][0] = corners[5][0]; grid[5][1] = corners[5][1]; grid[5][2] = corners[5][2];
        grid[6][0] = corners[6][0]; grid[6][1] = corners[6][1]; grid[6][2] = corners[6][2];
        grid[7][0] = corners[7][0]; grid[7][1] = corners[7][1]; grid[7][2] = corners[7][2];
        
        if (isPointInsideBox(grid, point) == false) {
            point[0] = (interceptAvg[0]) * cd->widthPerCell[0] + cellLower[0]; 
            point[1] = (interceptAvg[1]) * cd->widthPerCell[1] + cellLower[1]; 
            point[2] = (interceptAvg[2]) * cd->widthPerCell[2] + cellLower[2];
          
            correctPoint(point, grid, cd->widthPerCell);
        }
    }
    
    int pt = m_Geom->addVertex(point[0], point[1], point[2]);
    //if (pt == 9090) 
    //    cout << point[0] << ", " << point[1] << ", " << point[2] << endl;
    return pt;
}

double retrieveMaterialValueAtCoordinates(float faceCenter[3]) {
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
    //std::string filename = "data/Object26_part2.mnc";
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
        voxel_location[i] = (unsigned long)dvoxel_location[i];
    }

    //cout << "Voxel location of world coordinate (" << world_location[0] << ", " << world_location[1] << ", " << world_location[2] << ") is [" << voxel_location[0] << ", " << voxel_location[1] << ", " << voxel_location[2] << "]" << endl;
    
    result = miget_voxel_value(minc_volume, voxel_location, 3, &voxel);
    if (result != MI_NOERROR) cout << "Error in retrieving voxel value: " << result << endl;
    
    miclose_volume(minc_volume);
    //cout << "voxel value: " << voxel << endl;
    return (voxel);
}

void recurseAndCheck(OCell *c) {
    if (c != NULL) {
        if (c->hasChildren()) {
            for (int i = 0; i < 8; i++) {
                recurseAndCheck(c->m_Children[i]);
            }
        }
        else {
            if (c->getHasTetrahedra() && c->getIsAmbiguous()) {
                float **cornerCoord = c->getCorners();
                int *cornerMat = c->getMaterialValues();

                cout << "\t\t\tOCell: [" << c->m_Ind[0] << ", " << c->m_Ind[1] << ", " << c->m_Ind[2] << "]  Corners are: " << cornerMat[0] << ", " << cornerMat[1] << ", " << cornerMat[2] << ", " << cornerMat[3] << ", " << cornerMat[4] << ", " << cornerMat[5] << ", " << cornerMat[6] << ", " << cornerMat[7] << ", " << cornerMat[8] << endl;

                for (int j = 0; j < 9; j++) {
                    double MINCmat = retrieveMaterialValueAtCoordinates(cornerCoord[j]);

                    if (cornerMat[j] == 0) {
                        if (MINCmat != 0) {
                            cout << "\t\t\tERROR: Background material does not match. Corner: " << j << ", cornerMat = " << cornerMat[j] << " and MINCmat = " << MINCmat << " for cube index " << c->m_Ind[0] << ", " << c->m_Ind[1] << ", " << c->m_Ind[2] << endl;
                        }
                    }
                    else if (cornerMat[j] == 2) {
                        if (MINCmat != 26) {
                            cout << "\t\t\tERROR: Foreground material does not match. Corner: " << j << ", cornerMat = " << cornerMat[j] << " and MINCmat = " << MINCmat << " for cube index " << c->m_Ind[0] << ", " << c->m_Ind[1] << ", " << c->m_Ind[2] << endl;
                        }
                    }
                }
            }
        }
    }
    
}

void recurseAndPrint(OCell* c, int level, bool printAllLevel){
    int j; 
    if (c != NULL) {
        if (c->hasChildren()) {
            for (j = 0; j < 8; j++) {                
                recurseAndPrint(c->m_Children[j], level, printAllLevel); 
            }
        }
        if (!printAllLevel) {
            if ((int)c->m_Level == level)
                c->print(); 
        }
        else {
            c->print();
        }
    }
}


void CorrectOctree::testPrinting1(){}
/************ end CorrectOctree ***************/ 



/************ begin OctreeContouring *************/ 
OctreeContouring::OctreeContouring(CorrectOctree* oct){
    this->m_Oct = oct;   
    
}

void OctreeContouring::initializeAndRun(){
    
    
    
    contourTree(m_Oct->m_Root);   
    
    //recurseAndCheck(m_Oct->m_Root);
}

/*void OctreeContouring::cleanup() {
    // First create a vtk PolyData
    vtkSmartPointer<vtkPolyData> pdata_whole = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> pdata_points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < this->m_Oct->m_Geom->getNumVerts(); i++) {
        pdata_points->InsertNextPoint(this->m_Oct->m_Geom->getVert(i));
    }
    
    vtkSmartPointer<vtkCellArray> pdata_cells = vtkSmartPointer<vtkCellArray>::New();
    //vtkSmartPointer<vtkCellData> pdata_celldata = vtkSmartPointer<vtkCellData>::New();
    
    for (int i = 0; i < this->m_Oct->m_Geom->getNumTris(); i++) {
        vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
        
        unsigned int *ids = this->m_Oct->m_Geom->getTri(i);
        cellList->InsertNextId(ids[0]);
        cellList->InsertNextId(ids[1]);
        cellList->InsertNextId(ids[2]);
        
        pdata_cells->InsertNextCell(cellList);
    }
    
    pdata_whole->SetPoints(pdata_points);
    pdata_whole->SetPolys(pdata_cells);
    
    pdata_whole->BuildCells();
    pdata_whole->BuildLinks();
    
    // Cleanup. 
    for (int i = 0; i < pdata_whole->GetNumberOfCells(); i++) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        pdata_whole->GetCellPoints(i, cellPoints);
        
        vtkSmartPointer<vtkIdList> cellEdgeNeighbors0 = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> cellEdgeNeighbors1 = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> cellEdgeNeighbors2 = vtkSmartPointer<vtkIdList>::New();
        
        pdata_whole->GetCellEdgeNeighbors(i, cellPoints->GetId(0), cellPoints->GetId(1), cellEdgeNeighbors0);
        pdata_whole->GetCellEdgeNeighbors(i, cellPoints->GetId(1), cellPoints->GetId(2), cellEdgeNeighbors1);
        pdata_whole->GetCellEdgeNeighbors(i, cellPoints->GetId(2), cellPoints->GetId(0), cellEdgeNeighbors2);
        
        //cout << "Cell: " << i << " has edges: " << cellPoints->GetId(0) << ", " << cellPoints->GetId(1) << ", " << cellPoints->GetId(2) << ": has ";
        //cout << "neighbors: " << cellEdgeNeighbors0->GetNumberOfIds() << ", " << cellEdgeNeighbors1->GetNumberOfIds() << ", " << cellEdgeNeighbors2->GetNumberOfIds() << endl;
        if (cellEdgeNeighbors0->GetNumberOfIds() == 0 || cellEdgeNeighbors1->GetNumberOfIds() == 0 || cellEdgeNeighbors2->GetNumberOfIds() == 0) {
            cout << "Cell to remove: " << i << endl;
            pdata_whole->DeleteCell(i); 
        }
    }
    pdata_whole->RemoveDeletedCells();
    
}*/

//void OctreeContouring::processHoles() {
//    // First create a vtk PolyData
//    vtkSmartPointer<vtkPolyData> pdata_whole = vtkSmartPointer<vtkPolyData>::New();
//    vtkSmartPointer<vtkPoints> pdata_points = vtkSmartPointer<vtkPoints>::New();
//    for (int i = 0; i < this->m_Oct->m_Geom->getNumVerts(); i++) {
//        pdata_points->InsertNextPoint(this->m_Oct->m_Geom->getVert(i));
//    }
//    
//    vtkSmartPointer<vtkCellArray> pdata_cells = vtkSmartPointer<vtkCellArray>::New();
//    //vtkSmartPointer<vtkCellData> pdata_celldata = vtkSmartPointer<vtkCellData>::New();
//    
//    for (int i = 0; i < this->m_Oct->m_Geom->getNumTris(); i++) {
//        vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
//        
//        unsigned int *ids = this->m_Oct->m_Geom->getTri(i);
//        cellList->InsertNextId(ids[0]);
//        cellList->InsertNextId(ids[1]);
//        cellList->InsertNextId(ids[2]);
//        
//        pdata_cells->InsertNextCell(cellList);
//    }
//    
//    pdata_whole->SetPoints(pdata_points);
//    pdata_whole->SetPolys(pdata_cells);
//    
//    pdata_whole->BuildCells();
//    pdata_whole->BuildLinks();
//    
//    // Apply FeatureEdges filter
//    vtkSmartPointer<vtkFeatureEdges> fedges = vtkSmartPointer<vtkFeatureEdges>::New();
//    fedges->SetInput(pdata_whole);
//    fedges->BoundaryEdgesOn();
//    fedges->FeatureEdgesOff();
//    fedges->NonManifoldEdgesOff();
//    fedges->ManifoldEdgesOff();
//    fedges->ColoringOff();
//    fedges->Update();
//    
//    vtkSmartPointer<vtkPolyData> pdata = fedges->GetOutput();
//    pdata->BuildCells();
//    pdata->BuildLinks();
//    
//    // Remove isolated pieces by checking whether cells, i.e. edges have neighbors. 
//    bool stop = false;
//    while (stop == false) {
//        stop = true;
//        cout << "Pass..." << endl;
//        for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
//            vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
//            pdata->GetCellPoints(i, cellPoints);
//
//            vtkSmartPointer<vtkIdList> cellNeighborLeft = vtkSmartPointer<vtkIdList>::New();
//            vtkSmartPointer<vtkIdList> cellNeighborRight = vtkSmartPointer<vtkIdList>::New();
//            
//            pdata->GetCellEdgeNeighbors(i, cellPoints->GetId(0), cellPoints->GetId(0), cellNeighborLeft);
//            pdata->GetCellEdgeNeighbors(i, cellPoints->GetId(1), cellPoints->GetId(1), cellNeighborRight);
//
//            if ((cellNeighborLeft->GetNumberOfIds() == 0) || (cellNeighborRight->GetNumberOfIds() == 0)) {
//                pdata->DeleteCell(i);
//                stop = false;
//            }
//        }
//        pdata->RemoveDeletedCells();
//        pdata->BuildCells();
//        pdata->BuildLinks();
//    }
//    
//    int visited[pdata->GetNumberOfCells()];
//    for (int i = 0; i < pdata->GetNumberOfCells(); i++) visited[i] = 0;
//    
//    // First loop: locate neighboring lines of each line in one direction
//    for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
//        if (visited[i] == 0) { // If the ith line has not been visited. 
//            int currentCell = i;
//            vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
//            pdata->GetCellPoints(currentCell, cellPoints);      // The ith line. 
//
//            // Get the two points of the line
//            int first = cellPoints->GetId(0);
//            int second = cellPoints->GetId(1);
//            
//            std::vector<int> loop;
//            bool loop_found = true;
//            vtkSmartPointer<vtkIdList> cellNeighbor = vtkSmartPointer<vtkIdList>::New();
//            while (first != second) {
//                loop.push_back(currentCell);
//                pdata->GetCellEdgeNeighbors(currentCell, second, second, cellNeighbor); // Get the neighboring line using one point
//                
//                vtkSmartPointer<vtkIdList> cellNeighborPoints = vtkSmartPointer<vtkIdList>::New();
//                pdata->GetCellPoints(cellNeighbor->GetId(0), cellNeighborPoints);
//
//                // Assume that a line will have only one neighboring line for each point. 
//                if (cellNeighbor->GetNumberOfIds() == 1) {
//                    if (cellNeighborPoints->GetId(0) == second) {
//                        second = cellNeighborPoints->GetId(1);
//                        currentCell = cellNeighbor->GetId(0);
//                    }
//                    else if (cellNeighborPoints->GetId(1) == second) {
//                        second = cellNeighborPoints->GetId(0);
//                        currentCell = cellNeighbor->GetId(0);
//                    }
//                }
//                // Break out of loop if the line has more than one neighboring line for one point. 
//                else if (cellNeighbor->GetNumberOfIds() > 1) {
//                    loop_found = false;
//                    break;
//                }
//            }
//            
//            if (loop_found) {
//                loop.push_back(currentCell);
//                for (int k = 0; k < loop.size(); k++) {
//                    visited[loop.at(k)] = 1;
//                }
//                
//                // Store all point indices from this hole 
//                std::vector<int> pointsIndices;
//                for (int i = 0; i < loop.size(); i++) {
//                    vtkSmartPointer<vtkIdList> pts_ids = pdata->GetCell(loop.at(i))->GetPointIds();
//
//                    pointsIndices.push_back(pts_ids->GetId(0));
//                    pointsIndices.push_back(pts_ids->GetId(1));
//                }
//                
//                // Remove duplicate point indices. 
//                for (int p = 0; p < pointsIndices.size(); p++) {
//                    int pi = pointsIndices.at(p);
//                    for (int q = 0; q < pointsIndices.size(); q++) {
//                        int qi = pointsIndices.at(q);
//
//                        if (p != q) {
//                            if (pi == qi) {
//                                pointsIndices.erase(pointsIndices.begin() + q);
//                                break;
//                            }
//                        }
//                    }
//                }
//                
//                // Retrieve the point coordinates. 
//                std::vector<double*> points;
//                for (int i = 0; i < pointsIndices.size(); i++) {
//                    double *d0 = new double[3];
//                    pdata->GetPoint(pointsIndices.at(i), d0);
//                    points.push_back(d0);
//                }
//                
//                for (int i = 0; i < points.size(); i++) {
//                    double *d = points.at(i);
//
//                    for (int j = 0; j < this->m_Oct->m_Geom->getNumVerts(); j++) {
//                        float *p = this->m_Oct->m_Geom->getVert(j);
//
//                        if (compareFloats((float)d[0], p[0]) == true && compareFloats((float)d[1], p[1]) == true && compareFloats((float)d[2], p[2]) ==  true) {
//                            pointsIndices.at(i) = j;
//                        }
//                    }
//                }
//                
//                
//                // Triangulation by N-polygoning, and then triangulating
//                // First create a vtk PolyData
//                cout << "\tHole is made up of: ";
//                vtkSmartPointer<vtkPolyData> pdata_whole = vtkSmartPointer<vtkPolyData>::New();
//                vtkSmartPointer<vtkPoints> pdata_points = vtkSmartPointer<vtkPoints>::New();
//                for (int i = 0; i < pointsIndices.size(); i++) {
//                    float *point = this->m_Oct->m_Geom->getVert(pointsIndices.at(i));
//                    cout << pointsIndices.at(i) << ", ";
//                    pdata_points->InsertNextPoint(point[0], point[1], point[2]);
//                }
//                cout << endl;
//                
//                vtkSmartPointer<vtkCellArray> pdata_cells = vtkSmartPointer<vtkCellArray>::New();
//
//                vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
//                for (int i = 0; i < pointsIndices.size(); i++) {
//                    cellList->InsertNextId(i);
//                }
//                pdata_cells->InsertNextCell(cellList);
//
//                pdata_whole->SetPoints(pdata_points);
//                pdata_whole->SetPolys(pdata_cells);
//
//                pdata_whole->BuildCells();
//                pdata_whole->BuildLinks();
//
//                // Triangulate the polydata. 
//                vtkSmartPointer<vtkTriangleFilter> tri = vtkSmartPointer<vtkTriangleFilter>::New();
//                tri->SetInput(pdata_whole);
//                tri->Update();
//
//                vtkSmartPointer<vtkPolyData> pdata = tri->GetOutput();
//                pdata->BuildCells();
//                pdata->BuildLinks();
//                
//                for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
//                    vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
//                    pdata->GetCellPoints(i, cellPoints);
//
//                    int *tri = new int[3];
//                    tri[0] = pointsIndices[cellPoints->GetId(0)]; 
//                    tri[1] = pointsIndices[cellPoints->GetId(1)]; 
//                    tri[2] = pointsIndices[cellPoints->GetId(2)];
//
//                    m_Oct->m_Geom->addTriangle(tri[0], tri[1], tri[2]);
//                    m_Oct->m_GeomColors->addTriangle(999, 9999, 99999);
//                    cout << "\t\tHOLES vtkTriangle " << (this->m_Oct->m_Geom->getNumTris() - 1) <<  ": " << tri[0] << ", " << tri[1] << ", " << tri[2] << endl;
//                }
//                
//                
//                
//                
//                
//                
////                // Triangulation using CGAL
////                std::stringstream ss;
////                ss << "data/Tets/HOLES_";
////                for (int st = 0; st < pointsIndices.size(); st++) ss << pointsIndices.at(st) << "_";
////                ss << ".vtk";
////                
////                std::vector<int*> tris = getCGALTriangles(pointsIndices, ss.str());
////                for (int t = 0; t < tris.size(); t++) {
////                    int *triangle = tris.at(t);
////
////                    if (isDuplicateTriangle(triangle[0], triangle[1], triangle[2]) == false) {
////                        m_Oct->m_Geom->addTriangle(triangle[0], triangle[1], triangle[2]);
////                        m_Oct->m_GeomColors->addTriangle(999, 9999, 99999);
////
////                        m_SpecialTriangles.push_back(triangle);
////                    }
////                }
//                
//            }
//        }
//    }
//
//    // Second loop: locate neighboring lines of each line in the other direction
//    for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
//        if (visited[i] == 0) {
//            int currentCell = i;
//            vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
//            pdata->GetCellPoints(currentCell, cellPoints);
//
//            int first = cellPoints->GetId(1);
//            int second = cellPoints->GetId(0);
//            
//            std::vector<int> loop;
//            bool loop_found = true;
//            vtkSmartPointer<vtkIdList> cellNeighbor = vtkSmartPointer<vtkIdList>::New();
//            while (first != second) {
//                loop.push_back(currentCell);
//                pdata->GetCellEdgeNeighbors(currentCell, second, second, cellNeighbor);
//                
//                vtkSmartPointer<vtkIdList> cellNeighborPoints = vtkSmartPointer<vtkIdList>::New();
//                pdata->GetCellPoints(cellNeighbor->GetId(0), cellNeighborPoints);
//
//                if (cellNeighbor->GetNumberOfIds() == 1) {
//                    if (cellNeighborPoints->GetId(0) == second) {
//                        second = cellNeighborPoints->GetId(1);
//                        currentCell = cellNeighbor->GetId(0);
//                    }
//                    else if (cellNeighborPoints->GetId(1) == second) {
//                        second = cellNeighborPoints->GetId(0);
//                        currentCell = cellNeighbor->GetId(0);
//                    }
//                }
//                
//                else if (cellNeighbor->GetNumberOfIds() > 1) {
//                    loop_found = false;
//                    break;
//                }
//            }
//            
//            if (loop_found) {
//                loop.push_back(currentCell);
//                for (int k = 0; k < loop.size(); k++) {
//                    visited[loop.at(k)] = 1;
//                }
//                
//                // Store all point indices from this hole 
//                std::vector<int> pointsIndices;
//                for (int i = 0; i < loop.size(); i++) {
//                    vtkSmartPointer<vtkIdList> pts_ids = pdata->GetCell(loop.at(i))->GetPointIds();
//
//                    pointsIndices.push_back(pts_ids->GetId(0));
//                    pointsIndices.push_back(pts_ids->GetId(1));
//                }
//                
//                // Remove duplicate point indices. 
//                for (int p = 0; p < pointsIndices.size(); p++) {
//                    int pi = pointsIndices.at(p);
//                    for (int q = 0; q < pointsIndices.size(); q++) {
//                        int qi = pointsIndices.at(q);
//
//                        if (p != q) {
//                            if (pi == qi) {
//                                pointsIndices.erase(pointsIndices.begin() + q);
//                                break;
//                            }
//                        }
//                    }
//                }
//                
//                // Retrieve the point coordinates. 
//                std::vector<double*> points;
//                for (int i = 0; i < pointsIndices.size(); i++) {
//                    double *d0 = new double[3];
//                    pdata->GetPoint(pointsIndices.at(i), d0);
//                    points.push_back(d0);
//                }
//                
//                for (int i = 0; i < points.size(); i++) {
//                    double *d = points.at(i);
//
//                    for (int j = 0; j < this->m_Oct->m_Geom->getNumVerts(); j++) {
//                        float *p = this->m_Oct->m_Geom->getVert(j);
//
//                        if (compareFloats((float)d[0], p[0]) == true && compareFloats((float)d[1], p[1]) == true && compareFloats((float)d[2], p[2]) ==  true) {
//                            pointsIndices.at(i) = j;
//                        }
//                    }
//                }
//                
//                
//                // Triangulation by N-polygoning, and then triangulating
//                // First create a vtk PolyData
//                cout << "\tHole is made up of: ";
//                vtkSmartPointer<vtkPolyData> pdata_whole = vtkSmartPointer<vtkPolyData>::New();
//                vtkSmartPointer<vtkPoints> pdata_points = vtkSmartPointer<vtkPoints>::New();
//                for (int i = 0; i < pointsIndices.size(); i++) {
//                    float *point = this->m_Oct->m_Geom->getVert(pointsIndices.at(i));
//                    cout << pointsIndices.at(i) << ", ";
//                    pdata_points->InsertNextPoint(point[0], point[1], point[2]);
//                }
//                cout << endl;
//                
//                vtkSmartPointer<vtkCellArray> pdata_cells = vtkSmartPointer<vtkCellArray>::New();
//
//                vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
//                for (int i = 0; i < pointsIndices.size(); i++) {
//                    cellList->InsertNextId(i);
//                }
//                pdata_cells->InsertNextCell(cellList);
//
//                pdata_whole->SetPoints(pdata_points);
//                pdata_whole->SetPolys(pdata_cells);
//
//                pdata_whole->BuildCells();
//                pdata_whole->BuildLinks();
//
//                // Triangulate the polydata. 
//                vtkSmartPointer<vtkTriangleFilter> tri = vtkSmartPointer<vtkTriangleFilter>::New();
//                tri->SetInput(pdata_whole);
//                tri->Update();
//
//                vtkSmartPointer<vtkPolyData> pdata = tri->GetOutput();
//                pdata->BuildCells();
//                pdata->BuildLinks();
//                
//                for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
//                    vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
//                    pdata->GetCellPoints(i, cellPoints);
//
//                    int *tri = new int[3];
//                    tri[0] = pointsIndices[cellPoints->GetId(0)]; 
//                    tri[1] = pointsIndices[cellPoints->GetId(1)]; 
//                    tri[2] = pointsIndices[cellPoints->GetId(2)];
//
//                    m_Oct->m_Geom->addTriangle(tri[0], tri[1], tri[2]);
//                    m_Oct->m_GeomColors->addTriangle(999, 9999, 99999);
//                    cout << "\t\tHOLES vtkTriangle " << (this->m_Oct->m_Geom->getNumTris() - 1) <<  ": " << tri[0] << ", " << tri[1] << ", " << tri[2] << endl;
//                }
//                
////                // Triangulation using CGAL
////                std::stringstream ss;
////                ss << "data/Tets/HOLES_";
////                for (int st = 0; st < pointsIndices.size(); st++) ss << pointsIndices.at(st) << "_";
////                ss << ".vtk";
////                
////                std::vector<int*> tris = getCGALTriangles(pointsIndices, ss.str());
////                for (int t = 0; t < tris.size(); t++) {
////                    int *triangle = tris.at(t);
////
////                    if (isDuplicateTriangle(triangle[0], triangle[1], triangle[2]) == false) {
////                        m_Oct->m_Geom->addTriangle(triangle[0], triangle[1], triangle[2]);
////                        m_Oct->m_GeomColors->addTriangle(999, 9999, 99999);
////
////                        m_SpecialTriangles.push_back(triangle);
////                    }
////                }
//                
//            }
//        }
//    }
//}

/**
 * This function locates holes and attempts to patch the holes. 
 */
int OctreeContouring::processHoles() {
    int numberOfHoles = 0; // For counting the number of holes. 
    
    // First create a vtk PolyData using point and cell information in m_Geom
    vtkSmartPointer<vtkPolyData> pdata_whole = vtkSmartPointer<vtkPolyData>::New(); // The surface with holes. 
    vtkSmartPointer<vtkPoints> pdata_points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < this->m_Oct->m_Geom->getNumVerts(); i++) {
        pdata_points->InsertNextPoint(this->m_Oct->m_Geom->getVert(i));
    }
    
    vtkSmartPointer<vtkCellArray> pdata_cells = vtkSmartPointer<vtkCellArray>::New();
    //vtkSmartPointer<vtkCellData> pdata_celldata = vtkSmartPointer<vtkCellData>::New();
    
    for (int i = 0; i < this->m_Oct->m_Geom->getNumTris(); i++) {
        vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
        
        unsigned int *ids = this->m_Oct->m_Geom->getTri(i);
        cellList->InsertNextId(ids[0]);
        cellList->InsertNextId(ids[1]);
        cellList->InsertNextId(ids[2]);
        
        pdata_cells->InsertNextCell(cellList);
    }
    
    pdata_whole->SetPoints(pdata_points);
    pdata_whole->SetPolys(pdata_cells);
    
    pdata_whole->BuildCells();
    pdata_whole->BuildLinks();
    
    // Apply vtkFeatureEdges filter to locate boundary edges
    vtkSmartPointer<vtkFeatureEdges> fedges = vtkSmartPointer<vtkFeatureEdges>::New();
    fedges->SetInput(pdata_whole);
    fedges->BoundaryEdgesOn();
    fedges->FeatureEdgesOff();
    fedges->NonManifoldEdgesOff();
    fedges->ManifoldEdgesOff();
    fedges->ColoringOff();
    fedges->Update();
    
    vtkSmartPointer<vtkPolyData> pdata = fedges->GetOutput(); // Only the boundary edges of the surface
    pdata->BuildCells();
    pdata->BuildLinks();
    
    // Remove isolated pieces/edges by checking whether cells, i.e. edges have neighbors. 
    bool stop = false;
    while (stop == false) {
        stop = true;
        cout << "Pass..." << endl;
        for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
            vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
            pdata->GetCellPoints(i, cellPoints);

            vtkSmartPointer<vtkIdList> cellNeighborLeft = vtkSmartPointer<vtkIdList>::New();
            vtkSmartPointer<vtkIdList> cellNeighborRight = vtkSmartPointer<vtkIdList>::New();
            
            pdata->GetCellEdgeNeighbors(i, cellPoints->GetId(0), cellPoints->GetId(0), cellNeighborLeft);
            pdata->GetCellEdgeNeighbors(i, cellPoints->GetId(1), cellPoints->GetId(1), cellNeighborRight);

            if ((cellNeighborLeft->GetNumberOfIds() == 0) || (cellNeighborRight->GetNumberOfIds() == 0)) {
                pdata->DeleteCell(i);
                stop = false;
            }
        }
        pdata->RemoveDeletedCells();
        pdata->BuildCells();
        pdata->BuildLinks();
    }
    
    
    // First loop: locate neighboring lines of each line in one direction
    for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
        int currentCell = i;
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        pdata->GetCellPoints(currentCell, cellPoints);      // The ith line. 

        // Get the two points of the line
        int first = cellPoints->GetId(0);
        int second = cellPoints->GetId(1);

        std::vector<int> loop; // To store the points of a loop/cycle of boundary edges. 
        bool loop_found = true;
        vtkSmartPointer<vtkIdList> cellNeighbor = vtkSmartPointer<vtkIdList>::New();
        while (first != second) {
            loop.push_back(currentCell);
            pdata->GetCellEdgeNeighbors(currentCell, second, second, cellNeighbor); // Get the neighboring line using one point

            vtkSmartPointer<vtkIdList> cellNeighborPoints = vtkSmartPointer<vtkIdList>::New();
            pdata->GetCellPoints(cellNeighbor->GetId(0), cellNeighborPoints);

            // Assume that a line will have only one neighboring line for each point. 
            if (cellNeighbor->GetNumberOfIds() == 1) {
                if (cellNeighborPoints->GetId(0) == second) {
                    second = cellNeighborPoints->GetId(1);
                    currentCell = cellNeighbor->GetId(0);
                }
                else if (cellNeighborPoints->GetId(1) == second) {
                    second = cellNeighborPoints->GetId(0);
                    currentCell = cellNeighbor->GetId(0);
                }
            }
            // Break out of loop if the line has more than one neighboring line for one point. 
            else if (cellNeighbor->GetNumberOfIds() > 1) {
                loop_found = false;
                break;
            }
        }

        if (loop_found) {
            loop.push_back(currentCell);
            
            // Store all point indices from this hole 
            std::vector<int> pointsIndices;
            for (int ai = 0; ai < loop.size(); ai++) {
                vtkSmartPointer<vtkIdList> pts_ids = pdata->GetCell(loop.at(ai))->GetPointIds();

                pointsIndices.push_back(pts_ids->GetId(0));
                pointsIndices.push_back(pts_ids->GetId(1));
            }

            // Remove duplicate point indices. 
            for (int p = 0; p < pointsIndices.size(); p++) {
                int pi = pointsIndices.at(p);
                for (int q = 0; q < pointsIndices.size(); q++) {
                    int qi = pointsIndices.at(q);

                    if (p != q) {
                        if (pi == qi) {
                            pointsIndices.erase(pointsIndices.begin() + q);
                            break;
                        }
                    }
                }
            }
            
            
            // Point indices from vtkFeatureEdges are incorrect. 
            // So, retrieve the correct point coordinates and indices from m_Geom
            std::vector<double*> points;
            for (int bi = 0; bi < pointsIndices.size(); bi++) {
                double *d0 = new double[3];
                pdata->GetPoint(pointsIndices.at(bi), d0);
                points.push_back(d0);
            }

            for (int ci = 0; ci < points.size(); ci++) {
                double *d = points.at(ci);

                for (int j = 0; j < this->m_Oct->m_Geom->getNumVerts(); j++) {
                    float *p = this->m_Oct->m_Geom->getVert(j);

                    if (isSameFloat((float)d[0], p[0]) == true && isSameFloat((float)d[1], p[1]) == true && isSameFloat((float)d[2], p[2]) ==  true) {
                        pointsIndices.at(ci) = j;
                    }
                }
            }
            
            // Make sure that all points are in the proper clockwise/anti-clockwise sequence. Otherwise, overlapping polygons will occur. 
            bool rearranged = false;
            pdata_whole->BuildLinks();
            for (int qi = 0; qi < pointsIndices.size() - 1; qi++) {
                int pt0 = pointsIndices.at(qi);
                int pt1 = pointsIndices.at(qi + 1);
                
                int res = pdata_whole->IsEdge(pt0, pt1);
                if (res == 0) {
                    cout << "Loop points are out of sequence. Rearranging points for hole..." << endl;
                    int pt_popped = pointsIndices.at(qi + 1);
                    pointsIndices.erase(pointsIndices.begin() + qi + 1);
                    
                    pointsIndices.push_back(pt_popped);
                    qi = qi - 1;
                    rearranged = true;
                }
            }

            if (rearranged) {
                cout << "\tRearranged hole is ";
                for (int qwe = 0; qwe < pointsIndices.size(); qwe++) {
                    cout << pointsIndices.at(qwe) << ", ";
                }cout << endl;
            }
            
            // Triangulation by N-polygoning, and then triangulating
            // First create a vtk PolyData
            numberOfHoles = numberOfHoles + 1;
            cout << "\tHole is made up of: ";
            vtkSmartPointer<vtkPolyData> pdata_loop = vtkSmartPointer<vtkPolyData>::New();
            vtkSmartPointer<vtkPoints> pdata_points = vtkSmartPointer<vtkPoints>::New();
            for (int di = 0; di < pointsIndices.size(); di++) {
                float *point = this->m_Oct->m_Geom->getVert(pointsIndices.at(di));
                cout << pointsIndices.at(di) << ", ";
                pdata_points->InsertNextPoint(point[0], point[1], point[2]);
            }
            cout << endl;

            vtkSmartPointer<vtkCellArray> pdata_cells = vtkSmartPointer<vtkCellArray>::New();

            vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
            for (int ei = 0; ei < pointsIndices.size(); ei++) {
                cellList->InsertNextId(ei);
            }
            pdata_cells->InsertNextCell(cellList);

            pdata_loop->SetPoints(pdata_points);
            pdata_loop->SetPolys(pdata_cells);

            pdata_loop->BuildCells();
            pdata_loop->BuildLinks();

            // Triangulate the polydata. 
            vtkSmartPointer<vtkTriangleFilter> tri = vtkSmartPointer<vtkTriangleFilter>::New();
            tri->SetInput(pdata_loop);
            tri->Update();

            vtkSmartPointer<vtkPolyData> pdata_tri = tri->GetOutput();
            pdata_tri->BuildCells();
            pdata_tri->BuildLinks();

            if (pdata_tri->GetNumberOfCells() < 1) {
                cout << "\n\nNO POLYGONS GENERATED\n\n" << endl;
            }
            
            for (int fi = 0; fi < pdata_tri->GetNumberOfCells(); fi++) {
                vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
                pdata_tri->GetCellPoints(fi, cellPoints);

                int *tri = new int[3];
                tri[0] = pointsIndices[cellPoints->GetId(0)]; 
                tri[1] = pointsIndices[cellPoints->GetId(1)]; 
                tri[2] = pointsIndices[cellPoints->GetId(2)];

                m_Oct->m_Geom->addTriangle(tri[0], tri[1], tri[2]);
                m_Oct->m_GeomColors->addTriangle(999, 9999, 99999);
                cout << "\t\tHOLES vtkTriangle " << (this->m_Oct->m_Geom->getNumTris() - 1) <<  ": " << tri[0] << ", " << tri[1] << ", " << tri[2] << endl;
            }

            
            i = 0; // Reset the for loop
            // Mark and delete the detected loops. 
            for (int k = 0; k < loop.size(); k++) {
                pdata->DeleteCell(loop.at(k));
            }

            pdata->RemoveDeletedCells();
            pdata->BuildCells();
            pdata->BuildLinks();
        }
    }

    // Second loop: locate neighboring lines of each line in the other direction
    for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
        int currentCell = i;
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        pdata->GetCellPoints(currentCell, cellPoints);      // The ith line. 

        // Get the two points of the line
        int first = cellPoints->GetId(1);
        int second = cellPoints->GetId(0);

        std::vector<int> loop;
        bool loop_found = true;
        vtkSmartPointer<vtkIdList> cellNeighbor = vtkSmartPointer<vtkIdList>::New();
        while (first != second) {
            loop.push_back(currentCell);
            pdata->GetCellEdgeNeighbors(currentCell, second, second, cellNeighbor); // Get the neighboring line using one point

            vtkSmartPointer<vtkIdList> cellNeighborPoints = vtkSmartPointer<vtkIdList>::New();
            pdata->GetCellPoints(cellNeighbor->GetId(0), cellNeighborPoints);

            // Assume that a line will have only one neighboring line for each point. 
            if (cellNeighbor->GetNumberOfIds() == 1) {
                if (cellNeighborPoints->GetId(0) == second) {
                    second = cellNeighborPoints->GetId(1);
                    currentCell = cellNeighbor->GetId(0);
                }
                else if (cellNeighborPoints->GetId(1) == second) {
                    second = cellNeighborPoints->GetId(0);
                    currentCell = cellNeighbor->GetId(0);
                }
            }
            // Break out of loop if the line has more than one neighboring line for one point. 
            else if (cellNeighbor->GetNumberOfIds() > 1) {
                loop_found = false;
                break;
            }
        }

        if (loop_found) {
            loop.push_back(currentCell);
            
            // Store all point indices from this hole 
            std::vector<int> pointsIndices;
            for (int ai = 0; ai < loop.size(); ai++) {
                vtkSmartPointer<vtkIdList> pts_ids = pdata->GetCell(loop.at(ai))->GetPointIds();

                pointsIndices.push_back(pts_ids->GetId(0));
                pointsIndices.push_back(pts_ids->GetId(1));
            }

            // Remove duplicate point indices. 
            for (int p = 0; p < pointsIndices.size(); p++) {
                int pi = pointsIndices.at(p);
                for (int q = 0; q < pointsIndices.size(); q++) {
                    int qi = pointsIndices.at(q);

                    if (p != q) {
                        if (pi == qi) {
                            pointsIndices.erase(pointsIndices.begin() + q);
                            break;
                        }
                    }
                }
            }

            
            // Retrieve the point coordinates. 
            std::vector<double*> points;
            for (int bi = 0; bi < pointsIndices.size(); bi++) {
                double *d0 = new double[3];
                pdata->GetPoint(pointsIndices.at(bi), d0);
                points.push_back(d0);
            }

            for (int ci = 0; ci < points.size(); ci++) {
                double *d = points.at(ci);

                for (int j = 0; j < this->m_Oct->m_Geom->getNumVerts(); j++) {
                    float *p = this->m_Oct->m_Geom->getVert(j);

                    if (isSameFloat((float)d[0], p[0]) == true && isSameFloat((float)d[1], p[1]) == true && isSameFloat((float)d[2], p[2]) ==  true) {
                        pointsIndices.at(ci) = j;
                    }
                }
            }


            // Make sure that all points are in the proper clockwise/anti-clockwise sequence
            bool rearranged = false;
            pdata_whole->BuildLinks();
            for (int qi = 0; qi < pointsIndices.size() - 1; qi++) {
                int pt0 = pointsIndices.at(qi);
                int pt1 = pointsIndices.at(qi + 1);
                
                int res = pdata_whole->IsEdge(pt0, pt1);
                if (res == 0) {
                    cout << "Loop points are out of sequence. Rearranging points for hole..." << endl;
                    int pt_popped = pointsIndices.at(qi + 1);
                    pointsIndices.erase(pointsIndices.begin() + qi + 1);
                    
                    pointsIndices.push_back(pt_popped);
                    qi = qi - 1;
                    rearranged = true;
                }
            }

            if (rearranged) {
                cout << "\tRearranged hole is ";
                for (int qwe = 0; qwe < pointsIndices.size(); qwe++) {
                    cout << pointsIndices.at(qwe) << ", ";
                }cout << endl;
            }
            
            
            // Triangulation by N-polygoning, and then triangulating
            // First create a vtk PolyData
            numberOfHoles = numberOfHoles + 1;
            cout << "\tHole is made up of: ";
            vtkSmartPointer<vtkPolyData> pdata_loop = vtkSmartPointer<vtkPolyData>::New();
            vtkSmartPointer<vtkPoints> pdata_points = vtkSmartPointer<vtkPoints>::New();
            for (int di = 0; di < pointsIndices.size(); di++) {
                float *point = this->m_Oct->m_Geom->getVert(pointsIndices.at(di));
                cout << pointsIndices.at(di) << ", ";
                pdata_points->InsertNextPoint(point[0], point[1], point[2]);
            }
            cout << endl;

            vtkSmartPointer<vtkCellArray> pdata_cells = vtkSmartPointer<vtkCellArray>::New();

            vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
            for (int ei = 0; ei < pointsIndices.size(); ei++) {
                cellList->InsertNextId(ei);
            }
            
            pdata_cells->InsertNextCell(cellList);

            pdata_loop->SetPoints(pdata_points);
            pdata_loop->SetPolys(pdata_cells);

            pdata_loop->BuildCells();
            pdata_loop->BuildLinks();

            // Triangulate the polydata. 
            vtkSmartPointer<vtkTriangleFilter> tri = vtkSmartPointer<vtkTriangleFilter>::New();
            tri->SetInput(pdata_loop);
            tri->Update();

            vtkSmartPointer<vtkPolyData> pdata_tri = tri->GetOutput();
            pdata_tri->BuildCells();
            pdata_tri->BuildLinks();

            for (int fi = 0; fi < pdata_tri->GetNumberOfCells(); fi++) {
                vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
                pdata_tri->GetCellPoints(fi, cellPoints);

                int *tri = new int[3];
                tri[0] = pointsIndices[cellPoints->GetId(0)]; 
                tri[1] = pointsIndices[cellPoints->GetId(1)]; 
                tri[2] = pointsIndices[cellPoints->GetId(2)];

                m_Oct->m_Geom->addTriangle(tri[0], tri[1], tri[2]);
                m_Oct->m_GeomColors->addTriangle(999, 9999, 99999);
                cout << "\t\tHOLES vtkTriangle " << (this->m_Oct->m_Geom->getNumTris() - 1) <<  ": " << tri[0] << ", " << tri[1] << ", " << tri[2] << endl;
            }

            
            i = 0; // Reset the for loop
            // Mark and delete the detected loops. 
            for (int k = 0; k < loop.size(); k++) {
                pdata->DeleteCell(loop.at(k));
            }
            
            pdata->RemoveDeletedCells();
            pdata->BuildCells();
            pdata->BuildLinks();

        }
    }
    return numberOfHoles;
}


void OctreeContouring::contourTree(OCell* cell){
    
    //recurseAndCheck(m_Oct->m_Root);
    
    regularDCPolygonizerCounter = 0;
    minimalEdgeTriangulation_WithAmbiguousOCellCounter = 0;
    FaceDiagonalRule_Ambiguous_UnambiguousCounter = 0;
    FaceDiagonalRule_Ambiguous_AmbiguousCounter = 0;
    InteriorEdgeRuleCounter = 0;
    
    time_t t = time(0);
    struct tm *now = localtime(&t);
    cout << "Parsing octree started at: " << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << endl;
    
    
    
    cellProc(cell);
    
    cout << "\n\n**** POST PROCESSING ****" << endl;
    
//    std::vector<int> ambiguousCellBPoints; // To list all BiPoints in ambiguous OCells. 
//    
//    std::vector<OCell*> octreeCells; // Temporary vector to store and OCells during octree traversal. 
//    octreeCells.push_back(this->m_Oct->m_Root);
//    
//    int ambCount = 0;
//    while (octreeCells.size() != 0) {
//        OCell *c = octreeCells.at(octreeCells.size() - 1);
//        octreeCells.pop_back(); // Retrieve the last OCell from vector. 
//        
//        if (c->hasChildren() == true) {
//            for (int j = 0; j < 8; j++) {
//                octreeCells.push_back(c->m_Children[j]);
//            }
//        }
//        
//        // Check if OCell is a leaf cell 
//        if ((int)c->m_Level == this->m_Oct->m_CoarseData.level) {
//            if ((c->getIsAmbiguous() == true) && (c->getHasTetrahedra() == true)) {
//                TCell **tcells = c->getTetrahedralCells();
//                for (int t = 0; t < 12; t++) {
//                    if (tcells[t]->hasSignChangeEdge() == true) { // Make sure TCell has BiPoint. 
//                        ambiguousCellBPoints.push_back(tcells[t]->getBPoint());
//                    }
//                }
//                ambCount = ambCount + 1;
//            } 
//        }
//    }
//    
//    cout << "Number of ambiguous cells in octree: " << ambCount << endl;
//    //this->m_Oct->m_Geom->getTri()
    
    
    //int h = processHoles();
    //cout << "Total number of holes: " << h << endl;
    
    //cleanup();
    //processHoles();
    //cleanup();
    //processHoles();
    
    //vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    //writer->SetFileName("data/TEMP.vtk");
    //writer->SetInput(pdata);
    //writer->Write();
    
    cout << "\n\n";
    cout << "Number of poor quality triangles created by RegularDCPolygonizing: " << regularDCPolygonizerCounter << endl;
    cout << "Number of poor quality triangles created by minimalEdgeTriangulation_WithAmbiguousOCell(): " << minimalEdgeTriangulation_WithAmbiguousOCellCounter << endl;
    cout << "Number of poor quality triangles created by FaceDiagonalRule_Ambiguous_Unambiguous(): " << FaceDiagonalRule_Ambiguous_UnambiguousCounter << endl;
    cout << "Number of poor quality triangles created by FaceDiagonalRule_Ambiguous_Ambiguous(): " << FaceDiagonalRule_Ambiguous_AmbiguousCounter << endl;
    cout << "Number of poor quality triangles created by InteriorEdgeRule(): " << InteriorEdgeRuleCounter << endl;
    
    cout << "**** FINISHED POST PROCESSING ****\n\n" << endl;
}

/* algorithm from "Dual Contouring with Hermite Data" */
void OctreeContouring::cellProc(OCell* q){
    int j; 
  
    if(q == NULL)
        return;   

    if(q->hasChildren()){
        /* eight calls to cellProc */
        for(j=0; j<8;j++)
            cellProc(q->m_Children[j]); 

        /* twelve calls to faceProc */
        for(j=0; j<12; j++)
            faceProc(q->m_Children[c_FaceAdj[j][0]], q->m_Children[c_FaceAdj[j][1]]); 

        /* six calls to edgeProc */
        for(j=0; j<6; j++)
            edgeProc(q->m_Children[c_EdgeAdj[j][0]], q->m_Children[c_EdgeAdj[j][1]],
                     q->m_Children[c_EdgeAdj[j][2]], q->m_Children[c_EdgeAdj[j][3]]);        
    }
}

/* Given two cells find their relative positions */
/*  i.e. qB is x+, y+, or z+ in relation to qA */
int OctreeContouring::findRelPositionOf2Cells(OCell* qA, OCell* qB){
    int j; 
  
    while(qB->m_Level > qA->m_Level)
        qB = qB->m_Parent; 
  
    while(qA->m_Level > qB->m_Level)
        qA = qA->m_Parent;   

    for(j=0; j<3; j++)
        if(qB->m_Ind[j] > qA->m_Ind[j])
            switch(j){
            case 0: return XPLUS; 
            case 1: return YPLUS; 	
            case 2: return ZPLUS; 	
            }
	
    //printf("findRePositionOf2Cells: should not be here\n");   
    cout<<"findRePositionOf2Cells: should not be here"<<endl; 
    qA->print(); 
    qB->print();

    return 0; 
}

void OctreeContouring::faceProcMap1(OCell* qA, OCell* qB, OCell** out8, OCell** out16){
	OCell* qac[8]; 	// Children of qA
	OCell* qbc[8];  // Children of qB

        // If qA has children, assign them to qac[]. Otherwise assign NULL. 
	if(qA->hasChildren())
		for(int j=0;j<8;j++)
			qac[j] = qA->m_Children[j]; 
	else{
		for(int j=0;j<8;j++)
			qac[j] = NULL; 
	}

        // If qB has children, assign them to qbc[]. Otherwise assign NULL. 
	if(qB->hasChildren())
		for(int j=0;j<8;j++)
			qbc[j] = qB->m_Children[j]; 
	else{
		for(int j=0;j<8;j++)
			qbc[j] = NULL; 
	}

        
    switch (findRelPositionOf2Cells(qA, qB)){
    case XPLUS: 
        out8[0] = qac[1]; 
        out8[1] = qbc[0]; 
        out8[2] = qac[5]; 
        out8[3] = qbc[4]; 
        out8[4] = qac[7]; 
        out8[5] = qbc[6]; 
        out8[6] = qac[3]; 
        out8[7] = qbc[2]; 	

        out16[0] = qac[3]; 
        out16[1] = qbc[2]; 
        out16[2] = qbc[6]; 
        out16[3] = qac[7]; 

        out16[4] = qac[1]; 
        out16[5] = qac[3]; 
        out16[6] = qbc[2]; 
        out16[7] = qbc[0]; 

        out16[8] = qac[1]; 
        out16[9] = qbc[0]; 
        out16[10] = qbc[4]; 
        out16[11] = qac[5];

        out16[12] = qac[5]; 
        out16[13] = qac[7]; 
        out16[14] = qbc[6]; 
        out16[15] = qbc[4]; 	 		
	  
        return; 

    case YPLUS: 
        out8[0] = qac[2]; 
        out8[1] = qbc[0]; 
        out8[2] = qac[3]; 
        out8[3] = qbc[1]; 
        out8[4] = qac[7]; 
        out8[5] = qbc[5]; 
        out8[6] = qac[6]; 
        out8[7] = qbc[4]; 

        out16[0] = qac[2]; 
        out16[1] = qbc[0]; 
        out16[2] = qbc[1]; 
        out16[3] = qac[3]; 

        out16[4] = qac[3]; 
        out16[5] = qac[7]; 
        out16[6] = qbc[5]; 
        out16[7] = qbc[1]; 	

        out16[8] = qac[6]; 
        out16[9] = qbc[4]; 
        out16[10] = qbc[5]; 
        out16[11] = qac[7];
 
        out16[12] = qac[2]; 
        out16[13] = qac[6]; 
        out16[14] = qbc[4]; 
        out16[15] = qbc[0]; 	

        return; 
    case ZPLUS: 
        out8[0] = qac[4]; 
        out8[1] = qbc[0]; 
        out8[2] = qac[6]; 
        out8[3] = qbc[2]; 
        out8[4] = qac[7]; 
        out8[5] = qbc[3]; 
        out8[6] = qac[5]; 
        out8[7] = qbc[1]; 

        out16[0] = qac[5]; 
        out16[1] = qbc[1]; 
        out16[2] = qbc[3]; 
        out16[3] = qac[7]; 
	
        out16[4] = qac[6]; 
        out16[5] = qac[7]; 
        out16[6] = qbc[3]; 
        out16[7] = qbc[2]; 
 
        out16[8] = qac[4]; 
        out16[9] = qbc[0]; 
        out16[10] = qbc[2]; 
        out16[11] = qac[6]; 
 
        out16[12] = qac[4]; 
        out16[13] = qac[5]; 
        out16[14] = qbc[1]; 
        out16[15] = qbc[0]; 	
    
        return; 
    }	

    printf("faceProcMap1: should not be here\n"); 
}

void OctreeContouring::faceProc(OCell* qA, OCell* qB){
    int j; 
    OCell* faceRecCells[8]; 
    OCell* edgeRecCells[16]; 
//    OCell* tmpC; 
    bool qAchild = qA->hasChildren(); 
    bool qBchild = qB->hasChildren(); 

    
    if(!qAchild && !qBchild) // Do nothing if qA and qB both have no children. 
        return; 
  
	faceProcMap1(qA, qB, faceRecCells, edgeRecCells); 

    if(!qAchild){
        for(j=0; j<4; j++)
            faceRecCells[j*2] = qA; 

        edgeRecCells[0] = qA; 
        edgeRecCells[3] = qA;

        edgeRecCells[4] = qA;
        edgeRecCells[5] = qA;

        edgeRecCells[8] = qA;
        edgeRecCells[11] = qA;

        edgeRecCells[12] = qA;
        edgeRecCells[13] = qA;
    }

    if(!qBchild){
        for(j=0; j<4; j++)
            faceRecCells[j*2+1] = qB; 

        edgeRecCells[1] = qB;
        edgeRecCells[2] = qB;

        edgeRecCells[6] = qB;
        edgeRecCells[7] = qB;

        edgeRecCells[9] = qB;
        edgeRecCells[10] = qB;

        edgeRecCells[14] = qB;
        edgeRecCells[15] = qB;

    }

    for(j=0; j<4; j++)
        faceProc(faceRecCells[j*2], faceRecCells[j*2+1]); 
        
  
    for(j=0; j<4; j++)
        edgeProc(edgeRecCells[j*4], edgeRecCells[j*4+1], edgeRecCells[j*4+2], edgeRecCells[j*4+3]);      
}

/**
 * Tests for Cases 3, 6, 7, 10, 12, 13 ambiguity. 
 * All these cases have one or more face ambiguities. 
 * 
 * @param vals - Material values of the corners of the grid cube. 
 * @return True if cube has face ambiguity. False otherwise. 
 */
bool isFaceAmbiguous(int vals[8]) {
    if (
            // Face 0, 1, 2, 3
            ((vals[0] == vals[3]) && (vals[0] != vals[2]) && (vals[0] != vals[1]) && (vals[0] > 1)) || 
            ((vals[1] == vals[2]) && (vals[1] != vals[0]) && (vals[1] != vals[3]) && (vals[1] > 1)) ||
            
            // Face 4, 5, 6, 7
            ((vals[4] == vals[7]) && (vals[4] != vals[5]) && (vals[4] != vals[6]) && (vals[4] > 1)) || 
            ((vals[5] == vals[6]) && (vals[5] != vals[4]) && (vals[5] != vals[7]) && (vals[5] > 1)) || 
            
            // Face 0, 1, 4, 5
            ((vals[0] == vals[5]) && (vals[0] != vals[4]) && (vals[0] != vals[1]) && (vals[0] > 1)) || 
            ((vals[4] == vals[1]) && (vals[4] != vals[0]) && (vals[4] != vals[5]) && (vals[4] > 1)) || 
            
            // Face 2, 3, 6, 7
            ((vals[2] == vals[7]) && (vals[2] != vals[3]) && (vals[2] != vals[6]) && (vals[2] > 1)) || 
            ((vals[3] == vals[6]) && (vals[3] != vals[2]) && (vals[3] != vals[7]) && (vals[3] > 1)) || 
            
            // Face 0, 2, 4, 6
            ((vals[0] == vals[6]) && (vals[0] != vals[2]) && (vals[0] != vals[4]) && (vals[0] > 1)) || 
            ((vals[2] == vals[4]) && (vals[2] != vals[0]) && (vals[2] != vals[6]) && (vals[2] > 1)) || 
            
            // Face 1, 3, 5, 7
            ((vals[1] == vals[7]) && (vals[1] != vals[5]) && (vals[1] != vals[3]) && (vals[1] > 1)) || 
            ((vals[3] == vals[5]) && (vals[3] != vals[1]) && (vals[3] != vals[7]) && (vals[3] > 1))            
            ) {
        return true;
    }
    else {
        return false; 
    }
}

/**
 * Tests for Case 4 ambiguity. 
 * 
 * @param vals - Material values of the corners of the grid cube. 
 * @return True if cube has Case 4 ambiguity. False otherwise. 
 */
bool isCase4Ambiguous(int vals[8]) {
//    // This is for only normal Case 4 ambiguity where only foreground material values
//    // are considered. Background values are not considered. 
//    if (
//            ((vals[0] > 1) && 
//            (vals[0] == vals[7]) && 
//            (vals[0] != vals[1]) && 
//            (vals[0] != vals[2]) && 
//            (vals[0] != vals[3]) && 
//            (vals[0] != vals[4]) && 
//            (vals[0] != vals[5]) && 
//            (vals[0] != vals[6])) || 
//            
//            ((vals[1] > 1) && 
//            (vals[1] == vals[6]) && 
//            (vals[1] != vals[0]) && 
//            (vals[1] != vals[2]) && 
//            (vals[1] != vals[3]) && 
//            (vals[1] != vals[4]) && 
//            (vals[1] != vals[5]) && 
//            (vals[1] != vals[7])) ||     
//    
//            ((vals[2] > 1) && 
//            (vals[2] == vals[5]) && 
//            (vals[2] != vals[0]) && 
//            (vals[2] != vals[1]) && 
//            (vals[2] != vals[3]) && 
//            (vals[2] != vals[4]) && 
//            (vals[2] != vals[6]) && 
//            (vals[2] != vals[7])) || 
//            
//            ((vals[3] > 1) && 
//            (vals[3] == vals[4]) && 
//            (vals[3] != vals[0]) && 
//            (vals[3] != vals[1]) && 
//            (vals[3] != vals[2]) && 
//            (vals[3] != vals[5]) && 
//            (vals[3] != vals[6]) && 
//            (vals[3] != vals[7]))
//            ) {     
//        return true; 
//    }
//    else return false;
    
    
    // This is for normal Case 4 ambiguity, as well as complimentary Case 4 
    // ambiguity where diagonal corner values are background values. 
    // Works only for single material 2M_DC
    if (
            ((vals[0] >= 0) && 
            (vals[0] == vals[7]) && 
            (vals[0] != vals[1]) && 
            (vals[0] != vals[2]) && 
            (vals[0] != vals[3]) && 
            (vals[0] != vals[4]) && 
            (vals[0] != vals[5]) && 
            (vals[0] != vals[6])) || 
            
            ((vals[1] >= 0) && 
            (vals[1] == vals[6]) && 
            (vals[1] != vals[0]) && 
            (vals[1] != vals[2]) && 
            (vals[1] != vals[3]) && 
            (vals[1] != vals[4]) && 
            (vals[1] != vals[5]) && 
            (vals[1] != vals[7])) ||     
    
            ((vals[2] >= 0) && 
            (vals[2] == vals[5]) && 
            (vals[2] != vals[0]) && 
            (vals[2] != vals[1]) && 
            (vals[2] != vals[3]) && 
            (vals[2] != vals[4]) && 
            (vals[2] != vals[6]) && 
            (vals[2] != vals[7])) || 
            
            ((vals[3] >= 0) && 
            (vals[3] == vals[4]) && 
            (vals[3] != vals[0]) && 
            (vals[3] != vals[1]) && 
            (vals[3] != vals[2]) && 
            (vals[3] != vals[5]) && 
            (vals[3] != vals[6]) && 
            (vals[3] != vals[7]))
            ) {
        return true; 
    }
    //else return false;
    else if (
            // 0 and 7 are diagonally opposite
            ((vals[1] == vals[2]) && (vals[1] == vals[3]) && (vals[1] == vals[4]) && (vals[1] == vals[5]) && (vals[1] == vals[6]) && 
            (((vals[0] == vals[7]) && (vals[0] != vals[1])) || ((vals[0] != vals[7]) && vals[0] != vals[1] && vals[7] != vals[1]))) || 
            
            // 2 and 5 are diagonally opposite
            ((vals[0] == vals[1]) && (vals[0] == vals[3]) && (vals[0] == vals[4]) && (vals[0] == vals[6]) && (vals[0] == vals[7]) && 
            (((vals[2] == vals[5]) && (vals[2] != vals[0])) || ((vals[2] != vals[5]) && vals[2] != vals[0] && vals[5] != vals[0]))) || 
            
            // 1 and 6 are diagonally opposite 
            ((vals[0] == vals[2]) && (vals[0] == vals[3]) && (vals[0] == vals[4]) && (vals[0] == vals[5]) && (vals[0] == vals[7]) && 
            (((vals[1] == vals[6]) && (vals[1] != vals[0])) || ((vals[1] != vals[6]) && vals[1] != vals[0] && vals[6] != vals[0]))) || 
            
            // 3 and 4 are diagonally opposite
            ((vals[0] == vals[1]) && (vals[0] == vals[2]) && (vals[0] == vals[5]) && (vals[0] == vals[6]) && (vals[0] == vals[7]) && 
            (((vals[3] == vals[4]) && (vals[3] != vals[0])) || ((vals[3] != vals[4]) && vals[3] != vals[0] && vals[4] != vals[0])))
            ) {
        return true; 
    }
    else return false;
}

bool isCase13Ambiguous(vval maskvals[8]) {
    bool ret = false;
    int c[8];
    c[0] = (int)maskvals[0];
    c[1] = (int)maskvals[1];
    c[2] = (int)maskvals[2];
    c[3] = (int)maskvals[3];
    c[4] = (int)maskvals[4];
    c[5] = (int)maskvals[5];
    c[6] = (int)maskvals[6];
    c[7] = (int)maskvals[7];    
    
    if (
    
    (c[1] == c[2] && c[1] == c[4] && c[1] == c[7] && c[1] != c[0] && c[0] == c[3] && c[0] == c[5] && c[0] == c[6])    ||    (c[0] == c[3] && c[0] == c[5] && c[0] == c[6] && c[0] != c[1] && c[1] == c[2] && c[1] == c[4] && c[1] == c[7]) || // Face 0132 and Face 4576
    (c[3] == c[6] && c[3] == c[0] && c[3] == c[5] && c[3] != c[2] && c[2] == c[7] && c[2] == c[1] && c[2] == c[4])    ||    (c[2] == c[7] && c[2] == c[1] && c[2] == c[4] && c[2] != c[3] && c[3] == c[6] && c[3] == c[0] && c[3] == c[5]) || // Face 2376 and Face 0154
    (c[3] == c[5] && c[3] == c[0] && c[3] == c[6] && c[3] != c[1] && c[1] == c[7] && c[1] == c[2] && c[1] == c[4])    ||    (c[1] == c[7] && c[1] == c[2] && c[1] == c[4] && c[1] != c[3] && c[3] == c[5] && c[3] == c[0] && c[3] == c[6])   // Face 1375 and Face 0264
    ) {
        ret = true;
    }
    return ret;
}

/**
 * Determines if the grid cell is ambiguous (either face or interior ambiguous). 
 * 
 * @param maskvals - The material values at the corner of the grid cell. 
 * @return - True if the grid cell has ambiguity. False otherwise. 
 */
bool isAmbiguous(vval maskvals[8]) {
    int corners[8];
    corners[0] = (int)maskvals[0];
    corners[1] = (int)maskvals[1];
    corners[2] = (int)maskvals[2];
    corners[3] = (int)maskvals[3];
    corners[4] = (int)maskvals[4];
    corners[5] = (int)maskvals[5];
    corners[6] = (int)maskvals[6];
    corners[7] = (int)maskvals[7];    
    
    
    return (isFaceAmbiguous(corners) || isCase4Ambiguous(corners));    
}



/**
 * Computes and returns the euclidean distance between 2 points. 
 * 
 * @param p1 - First point in 3D. 
 * @param p2 - Second point in 3D. 
 * @return - The euclidean distance. 
 */
double euclideanDistance(float p1[3], float p2[3]) {
    return sqrt(((p2[0] - p1[0]) * (p2[0] - p1[0])) + ((p2[1] - p1[1]) * (p2[1] - p1[1])) + ((p2[2] - p1[2]) * (p2[2] - p1[2])));
}

/**
 * Computes and returns the euclidean distance between 2 points. 
 * 
 * @param x1 - X coordinate of first point. 
 * @param y1 - Y coordinate of first point. 
 * @param z1 - Z coordinate of first point. 
 * @param x2 - X coordinate of second point. 
 * @param y2 - Y coordinate of second point. 
 * @param z2 - Z coordinate of second point. 
 * 
 * @return - The euclidean distance. 
 */
double euclideanDistance(float x1, float y1, float z1, float x2, float y2, float z2) {
    return sqrt(((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) + ((z2 - z1) * (z2 - z1)));
}

/** 
 * Returns a double[4][3] array containing the coordinates of the shared face. 
 * Assumes that the OCells c1 and c2 are side-by-side.  
 */
float** OctreeContouring::getSharedFace(OCell *c1, OCell *c2) {
    float **face; 
    face = new float*[4];
    
    float **firstCorners = c1->getCorners();
    float **secondCorners = c2->getCorners();
    
    int count = 0;
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            if ((isSameFloat(firstCorners[i][0], secondCorners[j][0]) == true) && (isSameFloat(firstCorners[i][1], secondCorners[j][1]) == true) && (isSameFloat(firstCorners[i][2], secondCorners[j][2]) == true)) {
                face[count] = new float[3];
                face[count][0] = firstCorners[i][0]; face[count][1] = firstCorners[i][1]; face[count][2] = firstCorners[i][2];
                count = count + 1;
            }
        }
    }
    
    return face;
}

float** OctreeContouring::getSharedFace(TCell *t1, TCell *t2) {
    float **face; 
    face = new float*[3];
    
    float **corners1 = t1->getCorners();
    float **corners2 = t2->getCorners();
    
    int count = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if ((isSameFloat(corners1[i][0], corners2[j][0]) == true) && (isSameFloat(corners1[i][1], corners2[j][1]) == true) && (isSameFloat(corners1[i][2], corners2[j][2]) == true)) {
                face[count] = new float[3];
                face[count][0] = corners1[i][0];
                face[count][1] = corners1[i][1];
                face[count][2] = corners1[i][2];
                count = count + 1;
            }
        }
    }
    
    return face;
}

/**
 * Determines whether 2 TCells have a shared face. Basically check if 
 * the two TCells have 3 common points. 
 * 
 * @param t1 - First TCell
 * @param t2 - Second TCell
 * @return - True if the 2 TCells have a shared face, false otherwise. 
 */
bool OctreeContouring::hasSharedFace(TCell *t1, TCell *t2) {
    float **corners1 = t1->getCorners();
    float **corners2 = t2->getCorners();
    
    int count = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if ((isSameFloat(corners1[i][0], corners2[j][0]) == true) && (isSameFloat(corners1[i][1], corners2[j][1]) == true) && (isSameFloat(corners1[i][2], corners2[j][2]) == true)) {
                count = count + 1;
            }
        }
    }
    
    if (count == 3) return true;
    else return false;
}

/**
 * Determines whether a TCell and an OCell have a shared face. Basically check if 
 * the TCell and OCell have 3 common points. 
 * 
 * @param oc1 - The OCell
 * @param t1 - The TCell
 * @return - True if the OCell and TCell have a shared face, false otherwise. 
 */
bool OctreeContouring::hasSharedFace(OCell *oc1, TCell *t1) {
    float **corners1 = t1->getCorners();
    float **corners2 = oc1->getCorners();
    
    int count = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            if ((isSameFloat(corners1[i][0], corners2[j][0]) == true) && (isSameFloat(corners1[i][1], corners2[j][1]) == true) && (isSameFloat(corners1[i][2], corners2[j][2]) == true)) {
                count = count + 1;
            }
        }
    }
    
    if (count == 3) return true;
    else return false;
}

/**
 * Determines whether two OCells have a shared face. Basically check if 
 * the OCells have 4 common points. 
 * 
 * @param oc1 - First OCell
 * @param oc2 - Second OCell
 * @return - True if the OCells have a shared face, false otherwise. 
 */
bool OctreeContouring::hasSharedFace(OCell *oc1, OCell *oc2) {
    float **corners1 = oc1->getCorners();
    float **corners2 = oc2->getCorners();
    
    int count = 0;
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            if ((isSameFloat(corners1[i][0], corners2[j][0]) == true) && (isSameFloat(corners1[i][1], corners2[j][1]) == true) && (isSameFloat(corners1[i][2], corners2[j][2]) == true)) {
                count = count + 1;
            }
        }
    }
    
    if (count == 4) return true;
    else return false;
}

/**
 * Returns the coordinates of the common edge shared by the 4 OCells as well as the material indices of the edge points
 * @param qs - Array of the 4 neighboring OCells
 * 
 * @return - p and q are the coordinates of the shared edge, pq_material contains the material indices of p and q, respectively. 
 */
void OctreeContouring::getSharedEdge(OCell *qs[4], float p[3], float q[3], int pq_material[2]) {
    int edges[12][2] = {
        {0, 1}, {1, 3}, {3, 2}, {2, 0},
        {4, 5}, {5, 7}, {7, 6}, {6, 4}, 
        {1, 5}, {0, 4}, {2, 6}, {3, 7}, 
    };
    
    float **corners = qs[0]->getCorners();
    
    for (int i = 0; i < 12; i++) {
        float a[3] = {corners[edges[i][0]][0], corners[edges[i][0]][1], corners[edges[i][0]][2]};
        float b[3] = {corners[edges[i][1]][0], corners[edges[i][1]][1], corners[edges[i][1]][2]};
        
        if (qs[1]->hasEdge(a, b) && qs[2]->hasEdge(a, b) && qs[3]->hasEdge(a, b)) {
            p[0] = a[0]; p[1] = a[1]; p[2] = a[2];
            q[0] = b[0]; q[1] = b[1]; q[2] = b[2];
        }
    }
    
    pq_material[0] = qs[0]->getCornerMaterialValue(p);
    pq_material[1] = qs[0]->getCornerMaterialValue(q);
    
    
    //cout << "\t\t\t\t\t\t\t\t\t\t\t\tMaterial value of p: " << pq_material[0] << endl;
    //cout << "\t\t\t\t\t\t\t\t\t\t\t\tMaterial value of q: " << pq_material[1] << endl;
}

std::vector<int*> OctreeContouring::getTriangles(std::vector<int> pointIndex, std::string filename) {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkIdList> pointList = vtkSmartPointer<vtkIdList>::New();
    
    for (int i = 0; i < pointIndex.size(); i++) {
        float *pt = this->m_Oct->m_Geom->getVert(pointIndex.at(i));
        points->InsertPoint(i, pt);
        pointList->InsertNextId(i);
    }
    
    vtkSmartPointer<vtkCellArray> vertexArray = vtkSmartPointer<vtkCellArray>::New();
    vertexArray->InsertNextCell(pointList);
    
    vtkSmartPointer<vtkPolyData> pdata = vtkSmartPointer<vtkPolyData>::New();
    pdata->SetPoints(points);
    pdata->SetVerts(vertexArray);
    
    vtkSmartPointer<vtkDelaunay2D> del = vtkSmartPointer<vtkDelaunay2D>::New();
    del->SetProjectionPlaneMode(VTK_BEST_FITTING_PLANE);
    del->SetInput(pdata);
    del->SetOffset(100);
    del->Update();
    
    vtkSmartPointer<vtkPolyData> output = del->GetOutput();
    
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(output);
    writer->Update();
    
    std::vector<int*> ret;
    for (int i = 0; i < output->GetNumberOfCells(); i++) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        output->GetCellPoints(i, cellPoints);
        
        int *tri = new int[3];
        tri[0] = pointIndex[cellPoints->GetId(0)]; 
        tri[1] = pointIndex[cellPoints->GetId(1)]; 
        tri[2] = pointIndex[cellPoints->GetId(2)];
        
        ret.push_back(tri);
    }
    
    return ret;
}

std::vector<int*> OctreeContouring::getTriangles(int *pointIndex, float **pts, int numberOfPoints, std::string filename) {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkIdList> pointList = vtkSmartPointer<vtkIdList>::New();
    
    for (int i = 0; i < numberOfPoints; i++) {
        points->InsertPoint(i, pts[i][0], pts[i][1], pts[i][2]);
        pointList->InsertNextId(i);
    }
    
    vtkSmartPointer<vtkCellArray> vertexArray = vtkSmartPointer<vtkCellArray>::New();
    vertexArray->InsertNextCell(pointList);
    
    vtkSmartPointer<vtkPolyData> pdata = vtkSmartPointer<vtkPolyData>::New();
    pdata->SetPoints(points);
    pdata->SetVerts(vertexArray);
    
    vtkSmartPointer<vtkDelaunay2D> del = vtkSmartPointer<vtkDelaunay2D>::New();
    del->SetProjectionPlaneMode(VTK_BEST_FITTING_PLANE);
    del->SetInput(pdata);
    del->Update();
    
    vtkSmartPointer<vtkPolyData> output = del->GetOutput();
    
//    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
//    writer->SetFileName(filename.c_str());
//    writer->SetInput(output);
//    writer->Update();
    
    std::vector<int*> ret;
    for (int i = 0; i < output->GetNumberOfCells(); i++) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        output->GetCellPoints(i, cellPoints);
        
        int *tri = new int[3];
        tri[0] = pointIndex[cellPoints->GetId(0)]; 
        tri[1] = pointIndex[cellPoints->GetId(1)]; 
        tri[2] = pointIndex[cellPoints->GetId(2)];
        
        ret.push_back(tri);
    }
    
    return ret;
}

void OctreeContouring::printTriangle(int *triangle, std::string outputFile) {
    std::fstream file(outputFile.c_str(), std::fstream::out);
    int numberOfCells = 1; 
    int numberOfVertices = 3; 
    
    file << "# vtk DataFile Version 3.0\n";
    file << "Zhang Hybrid DC Triangulations\n"; 
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    file << "POINTS " << numberOfVertices << " float\n";
    
    float *vert0 = this->m_Oct->m_Geom->getVert(triangle[0]);
    float *vert1 = this->m_Oct->m_Geom->getVert(triangle[1]);
    float *vert2 = this->m_Oct->m_Geom->getVert(triangle[2]);
        
    file << vert0[0] << " " << vert0[1] << " " << vert0[2] << "\n";
    file << vert1[0] << " " << vert1[1] << " " << vert1[2] << "\n";
    file << vert2[0] << " " << vert2[1] << " " << vert2[2] << "\n";
        

    file << "\n";
    
    file << "POLYGONS " << numberOfCells << " " << (numberOfCells * 4) << "\n";
    for (int i = 0; i < numberOfVertices; i = i + 3) {
        file << "3 " << i << " " << (i + 1) << " " << (i + 2) << "\n";
    }
    file << "\n";

    file.close();
}

void OctreeContouring::printTriangles(std::vector<int*> tris, std::string outputFile) {
    
    std::fstream file(outputFile.c_str(), std::fstream::out);
    int numberOfCells = tris.size(); 
    int numberOfVertices = numberOfCells * 3; 
    
    file << "# vtk DataFile Version 3.0\n";
    file << "Zhang Hybrid DC Triangulations\n"; 
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    file << "POINTS " << numberOfVertices << " float\n";
    
    for (int i = 0; i < numberOfCells; i++) {
        int * t = tris.at(i);
        float *vert0 = this->m_Oct->m_Geom->getVert(t[0]);
        float *vert1 = this->m_Oct->m_Geom->getVert(t[1]);
        float *vert2 = this->m_Oct->m_Geom->getVert(t[2]);
        
        file << vert0[0] << " " << vert0[1] << " " << vert0[2] << "\n";
        file << vert1[0] << " " << vert1[1] << " " << vert1[2] << "\n";
        file << vert2[0] << " " << vert2[1] << " " << vert2[2] << "\n";
        
    }
    file << "\n";
    
    file << "POLYGONS " << numberOfCells << " " << (numberOfCells * 4) << "\n";
    for (int i = 0; i < numberOfVertices; i = i + 3) {
        file << "3 " << i << " " << (i + 1) << " " << (i + 2) << "\n";
    }
    file << "\n";

    file.close();
}

int getBestProjectionPlane(float **points, int numPoints) {
    double xPlane = 0, yPlane = 0, zPlane = 0;
    int count = 0;
    for (int i = 0; i < numPoints; i++) {
        for (int j = 0; j < numPoints; j++) {
            if (i != j) {
                xPlane = xPlane + euclideanDistance(0, points[i][1], points[i][2], 0, points[j][1], points[j][2]);
                yPlane = yPlane + euclideanDistance(points[i][0], 0, points[i][2], points[j][0], 0, points[j][2]);
                zPlane = zPlane + euclideanDistance(points[i][0], points[i][1], 0, points[j][0], points[j][1], 0);
                
                count = count + 1;
            }
        }
    }
    
    xPlane = xPlane / (double)count;
    yPlane = yPlane / (double)count;
    zPlane = zPlane / (double)count;
    
    //cout << "Avg distances: [" << xPlane << ", " << yPlane << ", " << zPlane << "]" << endl;
    
    if (xPlane >= yPlane && xPlane >= zPlane) {
        return 0;
    }
    else if (yPlane >= xPlane && yPlane >= zPlane) {
        return 1;
    }
    else if (zPlane >= xPlane && zPlane >= yPlane) {
        return 2;
    }
    else {
        return -1;
    }
}

std::vector<int*> OctreeContouring::getCGALTriangles(std::vector<int> pointIndex, std::string filename) {
    //cout << "\t\t\t\tCGAL Triangulating" << endl;
    ofstream f;
    f.open("data/CGAL/CGAL_TRIANGLES.cin");
    
    std::vector<float*> points;
    float **pointArray = new float*[pointIndex.size()];
    
    //cout << "pointIndex.size() = " << pointIndex.size() << endl;
    for (int i = 0; i < pointIndex.size(); i++) {
        
        float *p = this->m_Oct->m_Geom->getVert(pointIndex.at(i));
        points.push_back(p);        
        f << p[0] << " " << p[1] << " " << p[2] << "\n";
        
        pointArray[i] = new float[3];
        pointArray[i][0] = p[0];
        pointArray[i][1] = p[1];
        pointArray[i][2] = p[2];
    }
    f.close();
    
    int result = -1;
    int plane = getBestProjectionPlane(pointArray, pointIndex.size());
    if (plane == 0) {
        std::stringstream ss;
        ss << "data/CGAL/Delaunay2D_yzPlane " << filename.c_str();
        //ss << "data\\CGAL\\Delaunay2D_yzPlane.exe " << filename.c_str();
        result = system(ss.str().c_str());
        
        //cout << "\t\t\tyz plane" << endl;
    }
    else if (plane == 1) {
        std::stringstream ss;
        ss << "data/CGAL/Delaunay2D_xzPlane " << filename.c_str();
        //ss << "data\\CGAL\\Delaunay2D_xzPlane.exe " << filename.c_str();
        result = system(ss.str().c_str());
        
        //cout << "\t\t\txz plane" << endl;
    }
    else if (plane == 2) {
        std::stringstream ss;
        ss << "data/CGAL/Delaunay2D_xyPlane " << filename.c_str();
        //ss << "data\\CGAL\\Delaunay2D_xyPlane.exe " << filename.c_str();
        result = system(ss.str().c_str());
        
        //cout << "\t\t\txy plane" << endl;
    }
    else {
        cout << "CGAL Triangulation should not be here..." << endl;
    }
    
    std::vector<int*> tris;
    if (result == 0) {
        std::vector<int> ret;
        ifstream fi;
        fi.open("data/CGAL/CGAL_TRIANGLES.cout");
        int num, count = 0;
        fi >> num;
        
        while (count < num) {
            float x0,  y0, z0;
            fi >> x0 >> y0 >> z0;
            
            for (int i = 0; i < points.size(); i++) {
                float *tempPoint = points.at(i);
                if ((isSameFloat(tempPoint[0], x0) == true) && (isSameFloat(tempPoint[1], y0) == true) && (isSameFloat(tempPoint[2], z0) == true)) {
                    ret.push_back(pointIndex.at(i));

                    break;
                }
            } 
            
            count = count + 1;
        }
        fi.close();

        
        for (int i = 0; i < ret.size(); i = i + 3) {
            int *t = new int[3];
            t[0] = ret.at(i + 0);
            t[1] = ret.at(i + 1);
            t[2] = ret.at(i + 2);
            
            tris.push_back(t);
        }
    }
    
//    cout << "Number of polygons: " << tris.size() << "\t[";
//    for (int p = 0; p < tris.size(); p++) {
//        int *t = tris.at(p);
//        cout << t[0] << ", " << t[1] << ", " << t[2] << "], [";
//    }cout << endl;

    return tris;
}

/**
 * Computes a quality metric for triangles. 
 * Quality metric = 2 * (r_in) / (r_circ) 
 * where r_in is the radius of the circle inscribed in a triangle
 * and r_circ is the radius of the circumscribed circle
 * 
 * @return (2 * (r_in) / (r_circ))
 */
double OctreeContouring::computeRadiusRatio(int pointId1, int pointId2, int pointId3) {
    float *p1 = this->m_Oct->m_Geom->getVert(pointId1);
    float *p2 = this->m_Oct->m_Geom->getVert(pointId2);
    float *p3 = this->m_Oct->m_Geom->getVert(pointId3);
    
    Eigen::Vector3d a(p1[0], p1[1], p1[2]);
    Eigen::Vector3d b(p2[0], p2[1], p2[2]);
    Eigen::Vector3d c(p3[0], p3[1], p3[2]);
    
    Eigen::Vector3d r = a - c;
    Eigen::Vector3d s = b - c;
    Eigen::Vector3d t = a - b;
    
    Eigen::Vector3d cp = r.cross(s);
    double A = cp.norm() * 0.5;
    
    double norm_r2 = r.norm() * r.norm();
    double norm_s2 = s.norm() * s.norm();
    
    Eigen::Vector3d r2s_s2r = (norm_r2 * s) - (norm_s2 * r);
    Eigen::Vector3d cp_rs = r.cross(s);
    
    Eigen::Vector3d numerator = r2s_s2r.cross(cp_rs);
    double r_circ = numerator.norm() / (8 * A * A);
    
    double l1 = r.norm();
    double l2 = s.norm();
    double l3 = t.norm();
    
    double r_in = 2 * A / (l1 + l2 + l3);
    
    return (2 * r_in / r_circ) ;
}



/**
 * This method is used to generate a polygon using minimizers of TCells and OCells that share the minimal edge. 
 * Resulting polygon does not need to be convex. 
 * 
 * @param q1 - First OCell sharing the minimal edge. 
 * @param q2 - Second OCell sharing the minimal edge. 
 * @param q3 - Third OCell sharing the minimal edge. 
 * @param q4 - Fourth OCell sharing the minimal edge. 
 * @param v
 * @param mat
 * @param p - Coordinates of one point of the minimal edge. 
 * @param q - Coordinates of the other point of the minimal edge. 
 * @param pq_mat - The material values of the points of the minimal edge. First index is p's material value, second index is q's material value. 
 */
vector<int> OctreeContouring::minimalEdgeTriangulation_WithAmbiguousOCell(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[2], float p[3], float q[3], int pq_mat[2]) {
    std::vector<int> ret;
    
    std::vector<OCell*> olist; // List of all the OCells. 
    olist.push_back(q1);
    olist.push_back(q2);
    olist.push_back(q3);
    olist.push_back(q4);
    
    std::vector<int> pointIndices; // For storing all the point indices in order. 
    
    // Since this triangulation method is used when ambiguous OCells are present
    // there is always at least one ambiguous OCell. We start with that. 
    // However, the first 2 point indices from this OCell may not be in the same 
    // order as the rest of the pointIndices, so put the first OCell at the back 
    // of olist. Putting this OCell at the very end of olist will ensure that
    // its point indices will be added to pointIndices in the correct order. 
    
    // Remove the first two point indices from pointIndices at the end. 
    OCell *first; // Current OCell
    for (int i = 0; i < 4; i++) {
        OCell *temp = olist.at(i);
        if (temp->getIsAmbiguous() == true && temp->getHasTetrahedra() == true) {
            first = olist.at(i);
            
            olist.erase(olist.begin() + i);
            olist.push_back(first); 
            break;
        }
    }
    
    // Insert the bipoints of the 2 TCells sharing the minimal edge into pointIndices. 
    TCell **tcells = first->getTetrahedralCells();
    for (int i = 0; i < 12; i++) {
        if (tcells[i]->hasEdge(p, q) == true && tcells[i]->getBPoint() >= 0) {
            pointIndices.push_back(tcells[i]->getBPoint());
        }
    }
    
    
    // For all remaining OCells in olist, 
    int i = 0;
    while (olist.size() != 0) {
        OCell *oc = olist.at(i);
        
        // If oc and first have a shared face, 
        if (hasSharedFace(first, oc) == true) {
            // If oc is a normal OCell, then insert oc's bipoint into pointIndices
            // and remove oc from olist. Also check if oc shares a face with tc01
            if (oc->getIsAmbiguous() == false) {
                pointIndices.push_back(oc->getBPoint());
                
//                cout << "cube " << (int)oc->m_Ind[0] << ", " << (int)oc->m_Ind[1] << ", " << (int)oc->m_Ind[2] << " has BP: " << oc->getBPoint() << endl;
                
                first = oc; // Assign oc as the current OCell
                olist.erase(olist.begin() + i); // remove oc from the list of OCells
                i = -1; // Reset the counter. 
            }

            // If oc is an ambiguous OCell, 
            else if (oc->getIsAmbiguous() == true && oc->getHasTetrahedra() ==  true) {
                TCell **tempTs = oc->getTetrahedralCells(); // Retrieve the TCells of the ambiguous OCell
                
                TCell *ts[2];// = {t1, t2}; // An array to store the 2 TCells sharing the minimal edge. 
                
                int tcounter = 0;
                for (int j = 0; j < 12; j++) {
                    if (tempTs[j]->hasEdge(p, q) == true && tempTs[j]->getBPoint() >= 0) {
                        ts[tcounter] = tempTs[j];
                        tcounter = tcounter + 1;
                    }
                }
                
                
                // Once we have the 2 TCells sharing the minimal edge, check 
                // which one shares a face with first. 
                if (hasSharedFace(first, ts[0]) == true) {
                    pointIndices.push_back(ts[0]->getBPoint());
                    pointIndices.push_back(ts[1]->getBPoint());
                    
//                    cout << "BPt: " << ts[0]->getBPoint() << endl;
//                    cout << "BPt: " << ts[1]->getBPoint() << endl;
                    
                    first = oc; // Assign oc as the current OCell
                    olist.erase(olist.begin() + i); // remove oc from the list of OCells
                    i = -1; // Reset the counter. 
                }
                else if (hasSharedFace(first, ts[1]) == true) {
                    pointIndices.push_back(ts[1]->getBPoint());
                    pointIndices.push_back(ts[0]->getBPoint());
                    
//                    cout << "BPt: " << ts[1]->getBPoint() << endl;
//                    cout << "BPt: " << ts[0]->getBPoint() << endl;
                    
                    first = oc; // Assign oc as the current OCell
                    olist.erase(olist.begin() + i); // remove oc from the list of OCells
                    i = -1; // Reset the counter. 
                }
            }
        }
        
        i = i + 1;
    }
    
    // Remove the first 2 point indices from pointIndices. 
    // These 2 point indices are from the very first OCell, and may not 
    // be in the correct order. 
    pointIndices.erase(pointIndices.begin());
    pointIndices.erase(pointIndices.begin());
    
//    cout << "\n\n\tPolygon is made of points: ";
//    for (int i = 0; i < pointIndices.size(); i++) {
//        cout << pointIndices.at(i) << ", ";
//    }cout << "\n\n" << endl;
    
    
    
    
    
    
    std::stringstream ss;
    ss << "data/Tets/MANUAL_TRIS_";
    
    // First create a vtkPolyData
    vtkSmartPointer<vtkPolyData> pdata_whole = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> pdata_points = vtkSmartPointer<vtkPoints>::New();
    
//    cout << "\t\t\t\t\t" << pointIndices.size() << endl;
//    for (int q = 0; q < pointIndices.size(); q++) {
//        cout << pointIndices.at(q) << ", ";
//    }cout << endl;
    
    for (int i = 0; i < pointIndices.size(); i++) {
        float *point = this->m_Oct->m_Geom->getVert(pointIndices.at(i));
        //cout << point[0] << ", " << point[1] << ", " << point[2] << endl;
        
        pdata_points->InsertNextPoint(point[0], point[1], point[2]);
        ss << pointIndices.at(i) << "-";
    }
    ss << ".vtk";
    
    vtkSmartPointer<vtkCellArray> pdata_cells = vtkSmartPointer<vtkCellArray>::New();
    
    vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
    for (int i = 0; i < pointIndices.size(); i++) {
        cellList->InsertNextId(i);
    }
    pdata_cells->InsertNextCell(cellList);
    
    pdata_whole->SetPoints(pdata_points);
    pdata_whole->SetPolys(pdata_cells);
    
    pdata_whole->BuildCells();
    pdata_whole->BuildLinks();
    
    //cout << "\n\n" << pdata_whole->GetNumberOfPoints() << endl;
    
    // Triangulate the polydata. 
    vtkSmartPointer<vtkTriangleFilter> tri = vtkSmartPointer<vtkTriangleFilter>::New();
    tri->SetInput(pdata_whole);
    tri->Update();
    
    vtkSmartPointer<vtkPolyData> pdata = tri->GetOutput();
    pdata->BuildCells();
    pdata->BuildLinks();
    
    //cout << pdata->GetNumberOfPoints() << "\n\n" << endl;
    
    
    // Triangulation of simple polygon with n vertices will result in (n - 2) triangles
    if (pdata->GetNumberOfCells() == (pdata_whole->GetNumberOfPoints() - 2)) {
    //if (pdata->GetNumberOfCells() > 1) {
        for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
            vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
            pdata->GetCellPoints(i, cellPoints);

            int *tri = new int[3];
            tri[0] = pointIndices[cellPoints->GetId(0)]; 
            tri[1] = pointIndices[cellPoints->GetId(1)]; 
            tri[2] = pointIndices[cellPoints->GetId(2)];
            
            double ratio0 = computeRadiusRatio(tri[0], tri[1], tri[2]);
            if (ratio0 < 0.2) {
                //cout << "\t\tWarning: Poor quality triangle created by minimalEdgeTriangulation_WithAmbiguousOCell(): Triangle made of points " << tri[0] << ", " << tri[1] << " and " << tri[2] << endl;
                minimalEdgeTriangulation_WithAmbiguousOCellCounter++;
            }
                        
            m_Oct->m_Geom->addTriangle(tri[0], tri[1], tri[2]);
            m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
//            cout << "\t\tManual vtkTriangle " << (this->m_Oct->m_Geom->getNumTris() - 1) <<  ": " << tri[0] << ", " << tri[1] << ", " << tri[2] << endl;
            
            ret.push_back((this->m_Oct->m_Geom->getNumTris() - 1));
        }
        
//        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
//        writer->SetInput(pdata);
//        writer->SetFileName(ss.str().c_str());
//        writer->Update();
    }
    
    
    // vtkTriangleFilter failed to produce a triangulation of (n - 2) triangles
    else {
    //else if (pdata->GetNumberOfCells() <= 1) {
        
        //cout << "\t\tvtkTriangleFilter failed to produce a triangulation of (n - 2) triangles. Resorting to manual triangulation with centroid." << endl;
        // Generate a new point, the centroid of all the minimizers surrounding the minimal edge. 
        float centroid[3] = {0.0, 0.0, 0.0};
        for (int c = 0; c < pointIndices.size(); c++) {
            float *p = this->m_Oct->m_Geom->getVert(pointIndices.at(c));
            centroid[0] = centroid[0] + p[0];
            centroid[1] = centroid[1] + p[1];
            centroid[2] = centroid[2] + p[2];
        }

        centroid[0] = centroid[0] / pointIndices.size();
        centroid[1] = centroid[1] / pointIndices.size();
        centroid[2] = centroid[2] / pointIndices.size();

        int cnt = this->m_Oct->m_Geom->addVertex(centroid);
        pointIndices.push_back(cnt);


        // Manually generate triangles
        for (int t = 0; t < pointIndices.size() - 2; t++) {
            this->m_Oct->m_Geom->addTriangle(pointIndices.at(t), pointIndices.at(t + 1), pointIndices.at(pointIndices.size() - 1));
            //cout << "\t\tManual centroidTriangle " << (this->m_Oct->m_Geom->getNumTris() - 1) <<  ": " << pointIndices.at(t) << ", " << pointIndices.at(t + 1) << ", " << pointIndices.at(pointIndices.size() - 1) << endl;
            this->m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
           
            double ratio1 = computeRadiusRatio(pointIndices.at(t), pointIndices.at(t + 1), pointIndices.at(pointIndices.size() - 1));
            if (ratio1 < 0.2) {
                //cout << "\t\tWarning: Poor quality triangle created by minimalEdgeTriangulation_WithAmbiguousOCell()_manualtriangulation: Triangle made of points " << pointIndices.at(t) << " ," << pointIndices.at(t + 1) << " and " << pointIndices.at(pointIndices.size() - 1) << endl;
                minimalEdgeTriangulation_WithAmbiguousOCellCounter++;
            }
            
            ret.push_back((this->m_Oct->m_Geom->getNumTris() - 1));
        }
        
        this->m_Oct->m_Geom->addTriangle(pointIndices.at(0), pointIndices.at(pointIndices.size() - 2), pointIndices.at(pointIndices.size() - 1));
       
        double ratio2 = computeRadiusRatio(pointIndices.at(0), pointIndices.at(pointIndices.size() - 2), pointIndices.at(pointIndices.size() - 1));
        if (ratio2 < 0.2) {
            //cout << "\t\tWarning: Poor quality triangle created by minimalEdgeTriangulation_WithAmbiguousOCell()_manualtriangulation: Triangle made of points " << pointIndices.at(0) << "," << pointIndices.at(pointIndices.size() - 2) << " and " << pointIndices.at(pointIndices.size() - 1) << endl;
            minimalEdgeTriangulation_WithAmbiguousOCellCounter++;
        }
        
        //cout << "\t\tManual centroidTriangle " << (this->m_Oct->m_Geom->getNumTris() - 1) <<  ": " << pointIndices.at(0) << ", " << pointIndices.at(pointIndices.size() - 2) << ", " << pointIndices.at(pointIndices.size() - 1) << endl;
        this->m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
        
        ret.push_back((this->m_Oct->m_Geom->getNumTris() - 1));
    }
    
    return ret;
}

std::vector<int*> OctreeContouring::tetrahedralizeAndGetTriangles(std::vector<int> pointIndex, std::string filename) {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkIdList> pointList = vtkSmartPointer<vtkIdList>::New();
    
    for (int i = 0; i < pointIndex.size(); i++) {
        float *pt = this->m_Oct->m_Geom->getVert(pointIndex.at(i));
        points->InsertPoint(i, pt);
        pointList->InsertNextId(i);
    }
    
    vtkSmartPointer<vtkCellArray> vertexArray = vtkSmartPointer<vtkCellArray>::New();
    vertexArray->InsertNextCell(pointList);
    
    vtkSmartPointer<vtkPolyData> pdata = vtkSmartPointer<vtkPolyData>::New();
    pdata->SetPoints(points);
    pdata->SetVerts(vertexArray);
    
    vtkSmartPointer<vtkDelaunay3D> del3 = vtkSmartPointer<vtkDelaunay3D>::New();
    del3->SetInput(pdata);
    del3->SetAlpha(0);
    del3->SetOffset(0);
    del3->SetTolerance(0.001);
    del3->Update();
    
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = del3->GetOutput();
    //cout << "Cells in unstructured grid: " << unstructuredGrid->GetNumberOfCells() << endl;
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInput(unstructuredGrid);
    surfaceFilter->Update(); 
 
    vtkSmartPointer<vtkPolyData> output3D = surfaceFilter->GetOutput();
    //cout << "Cells in polydata: " << output3D->GetNumberOfCells() << endl;
    
    std::vector<int*> ret;
    for (int i = 0; i < output3D->GetNumberOfCells(); i++) {
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        output3D->GetCellPoints(i, cellPoints);
        
        int *tri = new int[3];
        tri[0] = pointIndex[cellPoints->GetId(0)]; 
        tri[1] = pointIndex[cellPoints->GetId(1)]; 
        tri[2] = pointIndex[cellPoints->GetId(2)];
        
        ret.push_back(tri);
    }
    
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(output3D);
    writer->Update();
    
    return ret;
}

bool OctreeContouring::isFlatish(std::vector<int*> tris, std::string outputFile) {
    // First create a vtk PolyData
    vtkSmartPointer<vtkPolyData> pdata_whole = vtkSmartPointer<vtkPolyData>::New();
    
    vtkSmartPointer<vtkPoints> pdata_points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < tris.size(); i++) {
        int *triangle = tris.at(i);
        pdata_points->InsertNextPoint(this->m_Oct->m_Geom->getVert(triangle[0]));
        pdata_points->InsertNextPoint(this->m_Oct->m_Geom->getVert(triangle[1]));
        pdata_points->InsertNextPoint(this->m_Oct->m_Geom->getVert(triangle[2]));
    }
    
    vtkSmartPointer<vtkCellArray> pdata_cells = vtkSmartPointer<vtkCellArray>::New();
    
    for (int i = 0; i < tris.size() * 3; i = i + 3) {
        vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
        
        cellList->InsertNextId(i);
        cellList->InsertNextId(i + 1);
        cellList->InsertNextId(i + 2);
        
        pdata_cells->InsertNextCell(cellList);
    }
    
    pdata_whole->SetPoints(pdata_points);
    pdata_whole->SetPolys(pdata_cells);
    
    pdata_whole->BuildCells();
    pdata_whole->BuildLinks();
    
    vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInput(pdata_whole);
    clean->Update();
            
    
    // Apply FeatureEdges filter
    vtkSmartPointer<vtkFeatureEdges> fedges = vtkSmartPointer<vtkFeatureEdges>::New();
    fedges->SetInput(clean->GetOutput());
    fedges->BoundaryEdgesOff();
    fedges->FeatureEdgesOn();
    fedges->NonManifoldEdgesOff();
    fedges->ManifoldEdgesOff();
    fedges->ColoringOff();
    fedges->SetFeatureAngle(90);
    fedges->Update();
    
    vtkSmartPointer<vtkPolyData> pdata = fedges->GetOutput();
    pdata->BuildCells();
    pdata->BuildLinks();
    
    bool ret = true;
    if (pdata->GetNumberOfCells() > 0) {
        ret = false;
        std::stringstream ss;
        ss << outputFile.replace(0, 10, "data/Tets/NOT_FLATISH_");
        
        printTriangles(tris, ss.str());
    }
    return ret;
}

void OctreeContouring::materialChangeEdge_WithBackgoundZhang(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[2], float p[3], float q[3], int pq_mat[2]) {
    OCell* qs[4] = {q1, q2, q3, q4}; 
    
    float t[3];
    int ept;
    if (pq_mat[0] > 0) {
        t[0] = p[0]; t[1] = p[1]; t[2] = p[2];
        //ept = this->m_Oct->m_Geom->addVertex(t[0], t[1], t[2]);
    }
    else if (pq_mat[1] > 0) {
        t[0] = q[0]; t[1] = q[1]; t[2] = q[2];
        //ept = this->m_Oct->m_Geom->addVertex(t[0], t[1], t[2]);
    }

    // Get all the minimizers from all the OCells/TCells that share the edge. 
    std::vector<int> minimizers; // Storing all the bi_Points. 
    for (int e = 0; e < 4; e++) {
        if (qs[e]->getIsAmbiguous() == true) {
//            cout << qs[e]->m_Ind[0] << ", " << qs[e]->m_Ind[1] << ", " << qs[e]->m_Ind[2] << endl;
//            int *matVals = qs[e]->getMaterialValues();
//            cout << "\tCorners: [";
//            for (int mt = 0; mt < 9; mt++) { cout << matVals[mt] << ", "; } cout << "]" << endl;
//
//            if (qs[e]->m_Ind[0] == 119 && qs[e]->m_Ind[1] == 142 && qs[e]->m_Ind[2] == 136) {
//                std::stringstream ss;
//                int indoc[] = {(int)qs[e]->m_Ind[0], (int)qs[e]->m_Ind[1], (int)qs[e]->m_Ind[2]};
//                ss << "/home/trash001/NetBeansProjects/Zhang_HybridDC/Zhang_Hybrid_DualContouring/data/Tets/XXXCube_" << indoc[0] << "_" << indoc[1] << "_" << indoc[2] << ".vtk";
//                qs[e]->vtkFile(ss.str());
//
//                cout << "getHasTetrahedra: " << qs[e]->getHasTetrahedra() << endl;
//                cout << "getIsAmbiguous: " << qs[e]->getIsAmbiguous() << endl;
//            }


            TCell **qsTetras = qs[e]->getTetrahedralCells();
            for (int tets = 0; tets < 12; tets++) {
                //int *ind = qsTetras[tets]->getIndex();
//                cout << "Tetrahedra Index: " << ind[0] << ", " << ind[1] << ", " << ind[2] << ", " << ind[3] << endl;
//                cout << "Edge: [" << p[0] << ", " << p[1] << ", " << p[2] << "] and [" << q[0] << ", " << q[1] << ", " << q[2] << "]" << endl;
                if ((qsTetras[tets]->hasEdge(p, q) == true) && (qsTetras[tets]->hasSignChangeEdge() == true) && qsTetras[tets]->getBPoint() >= 0) {
                //if ((qsTetras[tets]->hasVertex(t) == true) && (qsTetras[tets]->hasSignChangeEdge() == true) && qsTetras[tets]->getBPoint() >= 0) {
                //if ((qsTetras[tets]->hasEdge(p, q) == true)) {
                //if ((qsTetras[tets]->hasEdge(p, q) == true) || (qsTetras[tets]->hasVertex(t))) {
                //if ((qsTetras[tets]->hasVertex(t))) {
                    minimizers.push_back(qsTetras[tets]->getBPoint());
                    if (qsTetras[tets]->getBPoint() < 0) {
                        cout << endl;
                    }
                }
            }
        }
        else if (qs[e]->getIsAmbiguous() == false) {
            minimizers.push_back(qs[e]->getBPoint());
        }
    }

    ept = this->m_Oct->m_Geom->addVertex(t[0], t[1], t[2]);
    minimizers.push_back(ept);

    std::stringstream ss;
    ss << "data/Tets/CGAL_TRIS_";
    for (int pti = 0; pti < minimizers.size(); pti++) {
        ss << minimizers.at(pti) << "-";
    }
    ss << ".vtk";

    std::vector<int*> tris = tetrahedralizeAndGetTriangles(minimizers, ss.str());
    //std::vector<int*> tris = getTriangles(minimizers, ss.str());
    //std::vector<int*> tris = getCGALTriangles(minimizers, ss.str());
    for (int t = 0; t < tris.size(); t++) {
//        int *triangle = tris.at(t);
//        //cout << "CGAL Triangle: " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << endl;
//        if (isDuplicateTriangle(triangle[0], triangle[1], triangle[2]) == false) {
//            m_Oct->m_Geom->addTriangle(triangle[0], triangle[1], triangle[2]);
//            m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
//            
//            m_SpecialTriangles.push_back(triangle);
//        }
//        m_Oct->m_Geom->addTriangle(ept, triangle[1], triangle[2]);
//        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
//
//        m_Oct->m_Geom->addTriangle(triangle[0], ept, triangle[2]);
//        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
//
//        m_Oct->m_Geom->addTriangle(triangle[0], triangle[1], ept);
//        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
    }
}

void OctreeContouring::materialChangeEdge_WithBackgound(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[2], float p[3], float q[3], int pq_mat[2]) {
    OCell* qs[4] = {q1, q2, q3, q4}; 
    
    float t[3];
    int ept;
    if (pq_mat[0] > 0) {
        t[0] = p[0]; t[1] = p[1]; t[2] = p[2];
        //ept = this->m_Oct->m_Geom->addVertex(t[0], t[1], t[2]);
    }
    else if (pq_mat[1] > 0) {
        t[0] = q[0]; t[1] = q[1]; t[2] = q[2];
        //ept = this->m_Oct->m_Geom->addVertex(t[0], t[1], t[2]);
    }

    // Get all the minimizers from all the OCells/TCells that share the edge. 
    std::vector<int> minimizers; // Storing all the bi_Points. 
    for (int e = 0; e < 4; e++) {
        if (qs[e]->getIsAmbiguous() == true) {
//            cout << qs[e]->m_Ind[0] << ", " << qs[e]->m_Ind[1] << ", " << qs[e]->m_Ind[2] << endl;
//            int *matVals = qs[e]->getMaterialValues();
//            cout << "\tCorners: [";
//            for (int mt = 0; mt < 9; mt++) { cout << matVals[mt] << ", "; } cout << "]" << endl;
//
//            if (qs[e]->m_Ind[0] == 119 && qs[e]->m_Ind[1] == 142 && qs[e]->m_Ind[2] == 136) {
//                std::stringstream ss;
//                int indoc[] = {(int)qs[e]->m_Ind[0], (int)qs[e]->m_Ind[1], (int)qs[e]->m_Ind[2]};
//                ss << "/home/trash001/NetBeansProjects/Zhang_HybridDC/Zhang_Hybrid_DualContouring/data/Tets/XXXCube_" << indoc[0] << "_" << indoc[1] << "_" << indoc[2] << ".vtk";
//                qs[e]->vtkFile(ss.str());
//
//                cout << "getHasTetrahedra: " << qs[e]->getHasTetrahedra() << endl;
//                cout << "getIsAmbiguous: " << qs[e]->getIsAmbiguous() << endl;
//            }


            TCell **qsTetras = qs[e]->getTetrahedralCells();
            for (int tets = 0; tets < 12; tets++) {
                //int *ind = qsTetras[tets]->getIndex();
//                cout << "Tetrahedra Index: " << ind[0] << ", " << ind[1] << ", " << ind[2] << ", " << ind[3] << endl;
//                cout << "Edge: [" << p[0] << ", " << p[1] << ", " << p[2] << "] and [" << q[0] << ", " << q[1] << ", " << q[2] << "]" << endl;
                if ((qsTetras[tets]->hasEdge(p, q) == true) && (qsTetras[tets]->hasSignChangeEdge() == true) && qsTetras[tets]->getBPoint() >= 0) {
                //if ((qsTetras[tets]->hasVertex(t) == true) && (qsTetras[tets]->hasSignChangeEdge() == true) && qsTetras[tets]->getBPoint() >= 0) {
                //if ((qsTetras[tets]->hasEdge(p, q) == true)) {
                //if ((qsTetras[tets]->hasEdge(p, q) == true) || (qsTetras[tets]->hasVertex(t))) {
                //if ((qsTetras[tets]->hasVertex(t))) {
                    minimizers.push_back(qsTetras[tets]->getBPoint());
//                    if (qsTetras[tets]->getBPoint() < 0) {
//                        cout << endl;
//                    }
                }
            }
        }
        else if (qs[e]->getIsAmbiguous() == false) {
            minimizers.push_back(qs[e]->getBPoint());
        }
    }

//    ept = this->m_Oct->m_Geom->addVertex(t[0], t[1], t[2]);
//    minimizers.push_back(ept);

    std::stringstream ss;
    ss << "data/Tets/CGAL_TRIS_";
    for (int pti = 0; pti < minimizers.size(); pti++) {
        ss << minimizers.at(pti) << "-";
    }
    ss << ".vtk";

    //std::vector<int*> tris = tetrahedralizeAndGetTriangles(minimizers, ss.str());
    //std::vector<int*> tris = getTriangles(minimizers, ss.str());
    std::vector<int*> tris = getCGALTriangles(minimizers, ss.str());
    
    printTriangles(tris, ss.str());
    for (int t = 0; t < tris.size(); t++) {
        int *triangle = tris.at(t);
        //cout << "CGAL Triangle: " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << endl;
//        if (isDuplicateTriangle(triangle[0], triangle[1], triangle[2]) == false) {
            m_Oct->m_Geom->addTriangle(triangle[0], triangle[1], triangle[2]);
            m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);

            //m_SpecialTriangles.push_back(triangle);
            cout << "\t\tCGAL Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << endl;
//        }
    }
}

void OctreeContouring::materialChangeEdge_WithOtherMaterial() {
    
}

void OctreeContouring::interiorEdge() {
    
}

/**
 * One half of the FaceDiagonal Rule, for two ambiguous cubes having a shared face. 
 * For two ambiguous cubes having a shared face, if the face diagonal is a sign
 * change edge, then there are 4 valid TCells sharing the face diagonal. 
 * This function generates a polygon using the minimizers of the 4 TCells. 
 * 
 * @param q1, q2, q3, q4 - 4 OCells (ambiguous and unambiguous) sharing the minimal edge. 
 * @return - A vector containing all the generated triangles' indices. 
 */
vector<int> OctreeContouring::FaceDiagonalRule_Ambiguous_Ambiguous(OCell* q1, OCell* q2, OCell* q3, OCell* q4) {
    vector<int> ret;
    
    OCell *qs[4] = {q1, q2, q3, q4};
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i != j) {
                if ((qs[i]->getHasTetrahedra() == true && qs[i]->getIsAmbiguous() == true) && (qs[j]->getHasTetrahedra() == true && qs[j]->getIsAmbiguous() == true)) {
                    vector<int> pointIndices;
                    
                    // If both ambiguous cubes have a shared face
                    if (hasSharedFace(qs[i], qs[j]) == true) {
                        float **face = getSharedFace(qs[i], qs[j]);
                        
                        int ti[2], tj[2];
                        int icounter = 0, jcounter = 0;
                        
                        // Locate the 2 valid TCells of one ambiguous cube that share the face
                        TCell **iTCells = qs[i]->getTetrahedralCells();
                        for (int p = 0; p < 12; p++) {
                            if (iTCells[p]->isSubsetOfFace(face[0], face[1], face[2], face[3]) ==  true) {
                                if (iTCells[p]->hasSignChangeEdge() && iTCells[p]->getBPoint() >= 0) {
                                    ti[icounter] = p;
                                    icounter = icounter + 1;
                                }
                            }
                        }
                        
                        // Locate the 2 valid TCells of the other ambiguous cube that share the face
                        TCell **jTCells = qs[j]->getTetrahedralCells();
                        for (int k = 0; k < 12; k++) {
                            if (jTCells[k]->isSubsetOfFace(face[0], face[1], face[2], face[3]) ==  true) {
                                if (jTCells[k]->hasSignChangeEdge() && jTCells[k]->getBPoint() >= 0) {
                                    tj[jcounter] = k;
                                    jcounter = jcounter + 1;
                                }
                            }
                        }
                        
                        
                        
                        // Make sure there are 4 valid TCells sharing the face
                        if (icounter == 2 && jcounter == 2) {
                            bool makeTri = false; // Flag for determining whether or not to proceed with triangulation
                            int m1, m2;
                            
                            
                            // Locate the face diagonal. Since, at the moment there is no way to explicitly identify the face diagonal
                            // We are looking for the 1 edge that is shared by all 4 TCells.
                            // It will be the face diagonal. 
                            // If the material values of the two points of the 
                            // face diagonal are different, then makeTri = true, and proceed with triangulation
                            float **pts = iTCells[ti[0]]->getCorners(); // Ignore pts[0] because it is the center of the cube. 
                            if ((iTCells[ti[1]]->hasEdge(pts[1], pts[2]) == true) && (jTCells[tj[0]]->hasEdge(pts[1], pts[2]) == true) && (jTCells[tj[1]]->hasEdge(pts[1], pts[2]) == true)) {
                                m1 = qs[i]->getCornerMaterialValue(pts[1]);
                                m2 = qs[i]->getCornerMaterialValue(pts[2]);
                                
                                if (m1 != m2) {
                                    makeTri = true;
                                }
                            }
                            else if ((iTCells[ti[1]]->hasEdge(pts[2], pts[3]) == true) && (jTCells[tj[0]]->hasEdge(pts[2], pts[3]) == true) && (jTCells[tj[1]]->hasEdge(pts[2], pts[3]) == true)) {
                                m1 = qs[i]->getCornerMaterialValue(pts[2]);
                                m2 = qs[i]->getCornerMaterialValue(pts[3]);
                                
                                if (m1 != m2) {
                                    makeTri = true;
                                }
                            }
                            else if ((iTCells[ti[1]]->hasEdge(pts[1], pts[3]) == true) && (jTCells[tj[0]]->hasEdge(pts[1], pts[3]) == true) && (jTCells[tj[1]]->hasEdge(pts[1], pts[3]) == true)) {
                                m1 = qs[i]->getCornerMaterialValue(pts[1]);
                                m2 = qs[i]->getCornerMaterialValue(pts[3]);
                                
                                if (m1 != m2) {
                                    makeTri = true;
                                }
                            }
                            
                            // If the material values of the two points of the face diagonal are different, then triangulate
                            if (makeTri == true) {
                                int newv = 0;
                                vval newmat[2];

                                newmat[0] = m1;
                                newmat[1] = m2;

                                newv |= (1 << newmat[0]);
                                newv |= (1 << newmat[1]);
                                
                                
                                vector<int> cell; 
                                cell.push_back(iTCells[ti[0]]->getBPoint());
                                cell.push_back(iTCells[ti[1]]->getBPoint());
                                cell.push_back(jTCells[tj[0]]->getBPoint());
                                cell.push_back(jTCells[tj[1]]->getBPoint());
                                
                                std::stringstream ss;
                                ss << "data/Tets/FaceDiagonalRuleTri_" << iTCells[ti[0]]->getBPoint() << "_" << iTCells[ti[1]]->getBPoint() << "_" << jTCells[tj[0]]->getBPoint() << "_" << jTCells[tj[1]]->getBPoint() << ".vtk";
                                
                                std::vector<int*> tris = getCGALTriangles(cell, ss.str());
                               
                                //std::stringstream qsi;
                                //qsi << "data/Tets/FDR_" << (int)qs[i]->m_Ind[0] << "-" << (int)qs[i]->m_Ind[1] << "-" << (int)qs[i]->m_Ind[2] << ".vtk";
                                //qs[i]->vtkFile(qsi.str());

                                //std::stringstream qsj;
                                //qsj << "data/Tets/FDR_" << (int)qs[j]->m_Ind[0] << "-" << (int)qs[j]->m_Ind[1] << "-" << (int)qs[j]->m_Ind[2] << ".vtk";
                                //qs[j]->vtkFile(qsj.str());
                                
                                for (int a = 0; a < tris.size(); a++) {
                                    int *tri0 = tris.at(a);
                                    if (isDuplicateTriangle(tri0[0], tri0[1], tri0[2]) == false) {
                                        m_Oct->m_Geom->addTriangle(tri0[0], tri0[1], tri0[2]);
                                        m_Oct->m_GeomColors->addTriangle(newv, newmat[0], newmat[1]);
                                        //cout << "\t\tFaceDiagonal Rule Tri " << (this->m_Oct->m_Geom->getNumTris() - 1) <<  ": " << tri0[0] << ", " << tri0[1] << ", " << tri0[2] << endl;
                                        
                                        double ratio = computeRadiusRatio(tri0[0], tri0[1], tri0[2]);
                                        if (ratio < 0.2) {
                                            //cout << "\t\tWarning: Poor quality triangle created by InteriorEdgeRule(): Triangle made of points " << tri0[0] << ", " << tri0[1] << " and " << tri0[2] << endl;
                                            FaceDiagonalRule_Ambiguous_AmbiguousCounter++;
                                        }
                                        
                                        ret.push_back((this->m_Oct->m_Geom->getNumTris() - 1));

                                        m_SpecialTriangles.push_back(tri0);
                                    }
                                    else {
                                        //cout << "\t\tDuplicates triangles occured for Face Diagonal Rule tri: triangle " << tri0[0] << ", " << tri0[1] << ", " << tri0[2] << endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
    }
    
    return ret;
}

/**
 * Other half of the FaceDiagonal Rule, for an ambiguous cube and unambiguous cube having a shared face. 
 * For an ambiguous cube and an unambiguous cube having a shared face, if the 
 * face diagonal is a sign change edge, then there are two TCells for the 
 * ambiguous cube that make up the shared face. 
 * This function generate a triangle using the 2 minimizer of the two TCells 
 * as well as the minimizer of the unambiguous cube. 
 * 
 * @param q1, q2, q3, q4 - 4 OCells (ambiguous and unambiguous) sharing the minimal edge. 
 * @return - A vector containing all the generated triangles' indices. 
 */
vector<int> OctreeContouring::FaceDiagonalRule_Ambiguous_Unambiguous(OCell* q1, OCell* q2, OCell* q3, OCell* q4) {
    vector<int> ret;
    
    OCell* qs[4] = {q1, q2, q3, q4}; 
    
    for (int oc1 = 0; oc1 < 4; oc1++) {
        for (int oc2 = 0; oc2 < 4; oc2++) {
            if (oc1 != oc2) {
                // First locate ambiguous and unambiguous OCells. 
                if (hasSharedFace(qs[oc1], qs[oc2]) == true) {
                    float **face = getSharedFace(qs[oc1], qs[oc2]);

                    // Only if one OCell is ambiguous and the other is unambiguous. 
                    if (((qs[oc1]->getIsAmbiguous() == true) && (qs[oc2]->getIsAmbiguous() == false)) || ((qs[oc1]->getIsAmbiguous() == false) && (qs[oc2]->getIsAmbiguous() == true))) {
                        int unambiguous, ambiguous;
                        if (qs[oc1]->getIsAmbiguous() == true) {
                            ambiguous = oc1;
                            unambiguous = oc2;
                        }
                        else {
                            ambiguous = oc2;
                            unambiguous = oc1;
                        }

//                        //if (qs[ambiguous]->m_Ind[0] == 51 && qs[ambiguous]->m_Ind[1] == 62 && qs[ambiguous]->m_Ind[2] == 57) {
//                            std::stringstream ss1;
//                            ss1 << "data/Tets/AMB_" << qs[ambiguous]->m_Ind[0] << "_" << qs[ambiguous]->m_Ind[1] << "_" << qs[ambiguous]->m_Ind[2] << ".vtk";
//                            qs[ambiguous]->vtkFile(ss1.str());
//
//                            TCell **amT1 = qs[ambiguous]->getTetrahedralCells();
//
//                            for (int tc = 0; tc < 12; tc++) {
//                                std::stringstream ss;
//                                ss << "data/Tets/TET_" << amT1[tc]->getIndex()[0] << "_" << amT1[tc]->getIndex()[1] << "_" << amT1[tc]->getIndex()[2] << "_" << amT1[tc]->getIndex()[3] << ".vtk";
//                                amT1[tc]->vtkFile(ss.str());
//                            }
//                        //}

                        // Retrieve the two TCells that make up the face. 
                        int tets[2];
                        int tetcount = 0;
                        TCell **amT = qs[ambiguous]->getTetrahedralCells();
                        for (int tc = 0; tc < 12; tc++) {
                            if ((amT[tc]->isSubsetOfFace(face[0], face[1], face[2], face[3]) == true)) {
                                tets[tetcount] = tc;
                                tetcount++;
                            }
                        }
                        
                        
                        
                        if (hasSharedFace(amT[tets[0]], amT[tets[1]]) == true) { // Make sure that the two TCells have a shared face. 
                            float **face = getSharedFace(amT[tets[0]], amT[tets[1]]); // Get the shared face. 
                            
                            // The shared face between the two TCells would have 3 points.
                            // One of the points is the center of the ambiguous cube and 
                            // the other two points would make up the face diagonal. 
                            // Identify the center, and therefore identify the face diagonal. 
                            float **pts = qs[ambiguous]->getCorners();
                            float center[3];
                            center[0] = pts[8][0];
                            center[1] = pts[8][1];
                            center[2] = pts[8][2];
                            
                            float p1[3], p2[3];
                            if (isSameFloat(face[0][0], center[0]) == true && isSameFloat(face[0][1], center[1]) == true && isSameFloat(face[0][2], center[2]) == true) {
                                p1[0] = face[1][0];
                                p1[1] = face[1][1];
                                p1[2] = face[1][2];
                                
                                p2[0] = face[2][0];
                                p2[1] = face[2][1];
                                p2[2] = face[2][2];
                            }
                            else if (isSameFloat(face[1][0], center[0]) == true && isSameFloat(face[1][1], center[1]) == true && isSameFloat(face[1][2], center[2]) == true) {
                                p1[0] = face[0][0];
                                p1[1] = face[0][1];
                                p1[2] = face[0][2];
                                
                                p2[0] = face[2][0];
                                p2[1] = face[2][1];
                                p2[2] = face[2][2];
                            }
                            else if (isSameFloat(face[2][0], center[0]) == true && isSameFloat(face[2][1], center[1]) == true && isSameFloat(face[2][2], center[2]) == true) {
                                p1[0] = face[0][0];
                                p1[1] = face[0][1];
                                p1[2] = face[0][2];
                                
                                p2[0] = face[1][0];
                                p2[1] = face[1][1];
                                p2[2] = face[1][2];
                            }
                            
                            // Once the points of the face diagonal are identified, 
                            // retrieve their respective material values. 
                            int mat1 = qs[ambiguous]->getCornerMaterialValue(p1);
                            int mat2 = qs[ambiguous]->getCornerMaterialValue(p2);
                            
                            // If the material values are different, meaning the 
                            // face diagonal is a sign change edge, then create 
                            // a triangle.                             
                            if (mat1 != mat2) {
                                if (amT[tets[0]]->getBPoint() > 0 && amT[tets[1]]->getBPoint() > 0) {
                                    // Create a new triangle using the minimizers of the two TCells and the unambiguous OCell. 
                                    // First make sure this triangle does not already exist. 
                                    if (isDuplicateTriangle(qs[unambiguous]->getBPoint(), amT[tets[0]]->getBPoint(), amT[tets[1]]->getBPoint()) == false) {
                                        int newv = 0;
                                        newv |= (1 << mat1);
                                        newv |= (1 << mat2);
                                        
                                        m_Oct->m_Geom->addTriangle(qs[unambiguous]->getBPoint(), amT[tets[0]]->getBPoint(), amT[tets[1]]->getBPoint());
                                        m_Oct->m_GeomColors->addTriangle(newv, mat1, mat2);

                                        int *tri = new int[3];
                                        tri[0] = qs[unambiguous]->getBPoint(); 
                                        tri[1] = amT[tets[0]]->getBPoint(); 
                                        tri[2] = amT[tets[1]]->getBPoint();
                                        ret.push_back((m_Oct->m_Geom->getNumTris() - 1));

                                        double ratio = computeRadiusRatio(tri[0], tri[1], tri[2]);
                                        if (ratio < 0.2) {
                                            //cout << "\t\tWarning: Poor quality triangle created by FaceDiagonalRule_Ambiguous_Unambiguous(): Triangle made of points " << tri[0] << ", " << tri[1] << " and " << tri[2] << endl;
                                            FaceDiagonalRule_Ambiguous_UnambiguousCounter++;
                                        }
                                        
                                        //cout << "\t\tFace Rule Triangle for cubes ambiguous[" << (int)qs[ambiguous]->m_Ind[0] << ", " << (int)qs[ambiguous]->m_Ind[1] << ", " << (int)qs[ambiguous]->m_Ind[2] << "] and unambiguous[" << (int)qs[unambiguous]->m_Ind[0] << ", " << (int)qs[unambiguous]->m_Ind[1] << ", " << (int)qs[unambiguous]->m_Ind[2] << "]: " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << tri[0] << ", " << tri[1] << ", " << tri[2] << endl;

                                        m_SpecialTriangles.push_back(tri);                                
                                        //m_FaceRuleTriangles.push_back(m_Oct->m_Geom->getNumTris() - 1);

                                        //std::vector<int*> tris;
                                        //tris.push_back(tri);
                                        //std::stringstream ss;
                                        //ss << "data/Tets/FACE_TRIS-" << tri[0] << "-" << tri[1] << "-" << tri[2] << ".vtk";
                                        //printTriangles(tris, ss.str());
                                    }
                                    else {
                                        //cout << "\t\tDuplicates triangles occured for FaceRule: triangle " << qs[unambiguous]->getBPoint() << ", " << amT[tets[0]]->getBPoint() << ", " <<  amT[tets[1]]->getBPoint() << endl;
                                    }
                                }
                                else {
                                    //cout << "\n\nNegative bipoints detected\n\n" << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return ret;
}

/**
 * Individual ambiguous OCell rule implementation
 * For each ambiguous OCell sharing the minimal edge,
 * For each end point of the minimal edge not in the background
 * Retrieve all the TCells that contain the end point. 
 * Create triangle(s) using the minimizers of those TCells. 
    
 * @param q1 - OCell
 * @param q2 - OCell
 * @param q3 - OCell
 * @param q4 - OCell
 * @param v - 
 * @param mat - 
 * @param p - First endpoint of shared (minimal) edge. 
 * @param q - Second endpoint of shared (minimal) edge. 
 * @param pq_mat - The material indices of shared (minimal) edge. pq_mat[0] is material index for p, and  pq_mat[1] is material index for q. 
 */
//vector<int> OctreeContouring::individualAmbiguousCellRule(OCell* q1, OCell* q2, OCell* q3, OCell* q4, float p[3], float q[3], int pq_mat[2]) {
//    vector<int> ret;
//
//    OCell* qs[4] = {q1, q2, q3, q4}; 
//    
//    ofstream fmissedpoints("data/Tets/MissedPoints.txt", ios::app);
//    
//    for (int oc = 0; oc < 4; oc++) {
//        if (qs[oc]->getHasTetrahedra() == true && qs[oc]->getIsAmbiguous() == true) {
//            
//            TCell **TList = qs[oc]->getTetrahedralCells();
//            //int pArray[4], qArray[4], pCounter = 0, qCounter = 0;
//            std::vector<int> pArray, qArray;
//            for (int tc = 0; tc < 12; tc++) {
//                // Make sure that the first end point is not in the background
//                if (pq_mat[0] > 1) {
//                    if (TList[tc]->hasVertex(p) == true) {
//                        //pArray[pCounter] = tc;
//                        //pCounter = pCounter + 1;
//                        if (TList[tc]->getBPoint() >= 0) {
//                            pArray.push_back(TList[tc]->getBPoint());
//                        }
//                    }
//                }
//                // Make sure that the second end point is not in the background
//                if (pq_mat[1] > 1) {
//                    if (TList[tc]->hasVertex(q) == true) {
//                        //qArray[qCounter] = tc;
//                        //qCounter = qCounter + 1;
//                        if (TList[tc]->getBPoint() >= 0) {
//                            qArray.push_back(TList[tc]->getBPoint());
//                        }
//                    }
//                }
//            }
//
//            // Once all vertices are retrieved, triangulate. 
//            if (pArray.size() > 0) {
//                int newv = 0;
//                vval newmat[2];
//                
//                int *tm = qs[oc]->getMaterialValues();
//                
//                newmat[0] = tm[8];
//                newmat[1] = pq_mat[0];
//                
//                newv |= (1 << newmat[0]);
//                newv |= (1 << newmat[1]);
//                
//                cout << "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t" << (int)newmat[0] << ", " << (int)newmat[1] << ", " << newv << endl;
//
//                std::stringstream ss;
//                ss << "data/Tets/IAORule_p_";
//                for (int st = 0; st < pArray.size(); st++) ss << pArray[st] << "_";
//                ss << ".vtk";
//                std::vector<int*> tris = getCGALTriangles(pArray, ss.str());
//                
//                if (pArray.size() < 3) {
//                    cout << "pArray has less than 3 elements: ";
//                    for (int qw = 0; qw < pArray.size(); qw++) {
//                        cout << pArray[qw] << ", ";
//                        fmissedpoints << pArray[qw] << "\n";
//                    }
//                    cout << endl;
//                    
//                    std::stringstream fn;
//                    fn << "data/Tets/HasHoles_" << qs[oc]->m_Ind[0] << "-" << qs[oc]->m_Ind[1] << "-" << qs[oc]->m_Ind[2] << ".vtk";
//                    qs[oc]->vtkFile(fn.str());
//                }
//                
//                
//                for (int t = 0; t < tris.size(); t++) {
//                    int *triangle = tris.at(t);
//
//                    if (isDuplicateTriangle(triangle[0], triangle[1], triangle[2]) == false) {
//                        m_Oct->m_Geom->addTriangle(triangle[0], triangle[1], triangle[2]);
//                        m_Oct->m_GeomColors->addTriangle(newv, (int)newmat[0], (int)newmat[1]);
//                        
//                        if ((int)newmat[0] == (int)newmat[1]) {
//                            cout << "pSAME " << (int)newmat[0] << ", "<< (int)newmat[1] << endl;
//                        }
//                        
//                        cout << "\t\tIndividualAmbiguousCellRule_p Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << endl;
//                        
//                        ret.push_back((m_Oct->m_Geom->getNumTris() - 1));
//                        m_SpecialTriangles.push_back(triangle);
//                        
//                        std::stringstream pss;
//                        pss << "data/Tets/IAORule_p_" << triangle[0] << "_" << triangle[1] << "_" << triangle[2] << ".vtk";
//                        printTriangle(triangle, pss.str());
//                    }
//                    else {
//                        cout << "\t\tDuplicates triangles occured for individualAmbiguousCellRule_p: triangle " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << endl;
//                    }
//                }
//            }
//            if (qArray.size() > 0) {
//                int newv = 0;
//                vval newmat[2];
//                
//                int *tm = qs[oc]->getMaterialValues();
//                newmat[0] = tm[8];
//                newmat[1] = pq_mat[1];
//                
//                newv |= (1 << newmat[0]);
//                newv |= (1 << newmat[1]);
//                
//                cout << "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t" << (int)newmat[0] << ", " << (int)newmat[1] << ", " << newv << endl;
//                
//                std::stringstream ss;
//                ss << "data/Tets/IAORule_q_";
//                for (int st = 0; st < qArray.size(); st++) ss << qArray[st] << "_";
//                ss << ".vtk";
//                std::vector<int*> tris = getCGALTriangles(qArray, ss.str());
//                
//                if (qArray.size() < 3) {
//                    cout << "qArray has less than 3 elements: ";
//                    for (int qw = 0; qw < qArray.size(); qw++) {
//                        cout << qArray[qw] << ", ";
//                        fmissedpoints << qArray[qw] << "\n";
//                    }
//                    cout << endl;
//                    
//                    std::stringstream fn;
//                    fn << "data/Tets/HasHoles_" << qs[oc]->m_Ind[0] << "-" << qs[oc]->m_Ind[1] << "-" << qs[oc]->m_Ind[2] << ".vtk";
//                    qs[oc]->vtkFile(fn.str());
//                }
//                
//                
//                for (int t = 0; t < tris.size(); t++) {
//                    int *triangle = tris.at(t);
//
//                    if (isDuplicateTriangle(triangle[0], triangle[1], triangle[2]) == false) {
//                        m_Oct->m_Geom->addTriangle(triangle[0], triangle[1], triangle[2]);
//                        m_Oct->m_GeomColors->addTriangle(newv, (int)newmat[0], (int)newmat[1]);
//                        
//                        if ((int)newmat[0] == (int)newmat[1]) {
//                            cout << "qSAME " << (int)newmat[0] << ", "<< (int)newmat[1] << endl;
//                        }
//                        
//                        cout << "\t\tIndividualAmbiguousCellRule_q Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << endl;
//                        
//                        ret.push_back((m_Oct->m_Geom->getNumTris() - 1));
//                        m_SpecialTriangles.push_back(triangle);
//                        
//                        std::stringstream qss;
//                        qss << "data/Tets/IAORule_q_" << triangle[0] << "_" << triangle[1] << "_" << triangle[2] << ".vtk";
//                        printTriangle(triangle, qss.str());
//                    }
//                    cout << "\t\tDuplicates triangles occured for individualAmbiguousCellRule_p: triangle " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << endl;
//                }
//            }
//        }
//    }
//    
//    fmissedpoints.close();
//    
//    return ret;
//}

/**
 * Originally called individualAmbiguousCellRule(). 
 * This is modified to check whether the interior edges (an edge that is made 
 * up of a corner of the cube and the cube's center) of an AMBIGUOUS cube  is 
 * a sign change edges. For all 8 interior edges, if the interior edge is a 
 * sign change edges, then the function locates and identifies all valid TCells 
 * sharing the interior edge and generates a polygon using the minimizers 
 * of those TCells. 
 * 
 * @param q1, q2, q3, q4 - 4 OCells (ambiguous and unambiguous) sharing the minimal edge. 
 * @return - A vector containing all the generated triangles' indices. 
 */
vector<int> OctreeContouring::InteriorEdgeRule(OCell* q1, OCell* q2, OCell* q3, OCell* q4) {
    vector<int> ret;

    OCell* qs[4] = {q1, q2, q3, q4}; 
    
    //ofstream fmissedpoints("data/Tets/MissedPoints.txt", ios::app);
    
    // For all OCells
    for (int oc = 0; oc < 4; oc++) {
        if (qs[oc]->getHasTetrahedra() == true && qs[oc]->getIsAmbiguous() == true) { // If the OCell is ambiguous
            
            // Get the corners and their respective material values
            int *cornerMats = qs[oc]->getMaterialValues();
            float **corners = qs[oc]->getCorners();
            
            // Get the center and its material value. 
            int centerMat = cornerMats[8];
            float center[3] = {corners[8][0], corners[8][1], corners[8][2]};
            
            TCell **TList = qs[oc]->getTetrahedralCells(); // A list of all TCells
            
            // For all 8 interior edges
            for (int i = 0; i < 8; i++) {
                if (qs[oc]->getInteriorEdgeProcessed(i) == false) {
                    // Coordinates and material value of a corner
                    int tempCornerMat = cornerMats[i];
                    float tempCorner[3] = {corners[i][0], corners[i][1], corners[i][2]};

                    vector<int> pArray; // An array to contain the TCells' minimizers

                    // If the material value of the corner and center are different, 
                    // i.e., a sign change edge
                    // Add the TCell's minimizer into the array. 
                    if (centerMat != tempCornerMat) {
                        for (int j = 0; j < 12; j++) {
                            if (TList[j]->hasEdge(tempCorner, center) == true && TList[j]->hasSignChangeEdge()) {
                                pArray.push_back(TList[j]->getBPoint()); 
                            }
                        }
                    }

                    // Triangulation
                    if (pArray.size() > 0) {
                        int newv = 0;
                        vval newmat[2];

                        newmat[0] = centerMat;
                        newmat[1] = tempCornerMat;

                        newv |= (1 << newmat[0]);
                        newv |= (1 << newmat[1]);

                        std::stringstream ss;
                        ss << "data/Tets/IAORule_p_";
                        for (int st = 0; st < pArray.size(); st++) ss << pArray[st] << "_";
                        ss << ".vtk";
                        std::vector<int*> tris = getCGALTriangles(pArray, ss.str());

//                        if (pArray.size() < 3) {
//                            cout << "pArray has less than 3 elements: ";
//                            for (int qw = 0; qw < pArray.size(); qw++) {
//                                cout << pArray[qw] << ", ";
//                                fmissedpoints << pArray[qw] << "\n";
//                            }
//                            cout << endl;
//
//                            std::stringstream fn;
//                            fn << "data/Tets/HasHoles_" << qs[oc]->m_Ind[0] << "-" << qs[oc]->m_Ind[1] << "-" << qs[oc]->m_Ind[2] << ".vtk";
//                            qs[oc]->vtkFile(fn.str());
//                        }


                        for (int t = 0; t < tris.size(); t++) {
                            int *triangle = tris.at(t);

//                            if (isDuplicateTriangle(triangle[0], triangle[1], triangle[2]) == false) {
                                m_Oct->m_Geom->addTriangle(triangle[0], triangle[1], triangle[2]);
                                m_Oct->m_GeomColors->addTriangle(newv, (int)newmat[0], (int)newmat[1]);

                                double ratio = computeRadiusRatio(triangle[0], triangle[1], triangle[2]);
                                if (ratio < 0.2) {
                                    //cout << "\t\tWarning: Poor quality triangle created by InteriorEdgeRule(): Triangle made of points " << triangle[0] << ", " << triangle[1] << " and " << triangle[2] << endl;
                                    InteriorEdgeRuleCounter++;
                                }
                                
                                //cout << "\t\tIndividualAmbiguousCellRule_p Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << endl;

                                ret.push_back((m_Oct->m_Geom->getNumTris() - 1));
//                                m_SpecialTriangles.push_back(triangle);

//                                std::stringstream pss;
//                                pss << "data/Tets/IAORule_p_" << triangle[0] << "_" << triangle[1] << "_" << triangle[2] << ".vtk";
//                                printTriangle(triangle, pss.str());
//                            }
//                            else {
//                                cout << "\t\tDuplicates triangles occured for individualAmbiguousCellRule_p: triangle " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << endl;
//                            }
                        }
                        
                        qs[oc]->setInteriorEdgeProcessed(i, true);
                    }

                }
            }
        }
    }
    
    //fmissedpoints.close();
    
    return ret;
}

bool OctreeContouring::isDuplicateTriangle(int a, int b, int c) {
    bool ret = false;
    for (int i = 0; i < m_SpecialTriangles.size(); i++) {
        int *tri = m_SpecialTriangles.at(i);
        
        if (
                ((a == tri[0]) && (b == tri[1]) && (c == tri[2])) || 
                ((a == tri[0]) && (c == tri[1]) && (b == tri[2])) || 
                ((b == tri[0]) && (a == tri[1]) && (c == tri[2])) || 
                ((b == tri[0]) && (c == tri[1]) && (a == tri[2])) || 
                ((c == tri[0]) && (a == tri[1]) && (b == tri[2])) || 
                ((c == tri[0]) && (b == tri[1]) && (a == tri[2]))
                ) {
            ret = true;
        }
    }
    return ret;
}


/**
 * This process works in a similar fashion to the minimal edge rule. 
 * 
 * 
 * @param q1
 * @param q2
 * @param q3
 * @param q4
 * @param v
 * @param mat
 * @param p
 * @param q
 * @param pq_mat
 */
void OctreeContouring::minimalTetraEdge(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[2], float p[3], float q[3], int pq_mat[2]) {
    OCell* qs[4] = {q1, q2, q3, q4}; 
    
    for (int i = 0; i < 4; i++) {
        if ((qs[i]->getIsAmbiguous() == true) && (qs[i]->getHasTetrahedra() == true)) {
            
            // Retrieve the ambiguous OCell's center the and center's material value. 
            float **corners = qs[i]->getCorners();
            float center[3];
            center[0] = corners[8][0];
            center[1] = corners[8][1];
            center[2] = corners[8][2];
            
            int centerMaterial2 = qs[i]->getCornerMaterialValue(center); // Works for single material contouring
             
            //cout << "\t\t\t\t\t\t\t\tCenter is: " << center[0] << ", " << center[1] << ", " << center[2] << endl;
            //cout << "\t\t\t\t\t\t\tCorners Center is: " << corners[8][0] << ", " << corners[8][1] << ", " << corners[8][2] << endl;
            
            
            int *matVals = qs[i]->getMaterialValues(); 
            int centerMaterial = matVals[8]; 
            
            
            cout << "\t\t\t\t\t\t\t\t\t\t\t\t\tcenterMaterial: " << centerMaterial << endl;
            cout << "\t\t\t\t\t\t\t\t\t\t\t\t\tcenterMaterial2: " << centerMaterial2 << endl;
            
            // If the first endpoint of the minimal edge and the center has different material values
            // Also make sure that the first endpoint of the minimal edge is not in background. 
            if (centerMaterial != pq_mat[0] && pq_mat[0] != 0) {
                //if (centerMaterial == 0 || pq_mat[0] == 0) {
                //if (centerMaterial != pq_mat[0] && pq_mat[0] != 0 && centerMaterial != pq_mat[0]) {
                //if (centerMaterial != pq_mat[0] && pq_mat[0] != 0 && (centerMaterial == 0 || pq_mat[0] == 0)) {
                    std::vector<int> pointIndices;

                    TCell **tcells = qs[i]->getTetrahedralCells();
                    for (int t = 0; t < 12; t++) {
                        if ((tcells[t]->hasEdge(center, p) == true) && (tcells[t]->getBPoint() >= 0)) {
                            pointIndices.push_back(tcells[t]->getBPoint());
                        }
                    }

                    if (pointIndices.size() < 3) {
                        cout << "\n\t\t\t\tNO TRIS for points ";
                        for (int a = 0; a < pointIndices.size(); a++) {
                            cout << pointIndices.at(a) << ", ";
                        }cout << endl << endl;
                    }

                    std::stringstream ss;
                    ss << "data/Tets/MinimalTetraEdge_";
                    for (int st = 0; st < pointIndices.size(); st++) ss << pointIndices[st] << "_";
                    ss << ".vtk";
                    std::vector<int*> tris = getCGALTriangles(pointIndices, ss.str());
                    printTriangles(tris, ss.str());
                    for (int t = 0; t < tris.size(); t++) {
                        int *triangle = tris.at(t);

                        if (isDuplicateTriangle(triangle[0], triangle[1], triangle[2]) == false) {
                            m_Oct->m_Geom->addTriangle(triangle[0], triangle[1], triangle[2]);
                            m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);

                            m_SpecialTriangles.push_back(triangle);
//                            cout << "\t\tMinimal TetraEdge Rule Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << endl;
                        }
                    }
                //}
            }
            
            // If the second endpoint of the minimal edge and the center has different material values
            // Also make sure that the second endpoint of the minimal edge is not in background. 
            else if (centerMaterial != pq_mat[1] && pq_mat[1] != 0) {
            //else if (centerMaterial != pq_mat[1] && pq_mat[1] != 0 && centerMaterial != pq_mat[1]) {
            //if (centerMaterial != pq_mat[1] && pq_mat[1] != 0 && (centerMaterial == 0 || pq_mat[1] == 0)) {
                //if (centerMaterial == 0 || pq_mat[1] == 0)  {
                    std::vector<int> pointIndices;

                    TCell **tcells = qs[i]->getTetrahedralCells();
                    for (int t = 0; t < 12; t++) {
                        if ((tcells[t]->hasEdge(center, q) == true) && (tcells[t]->getBPoint() >= 0)) {
                            pointIndices.push_back(tcells[t]->getBPoint());
                        }
                    }

                    if (pointIndices.size() < 3) {
                        cout << "\n\t\t\t\tNO TRIS for points ";
                        for (int a = 0; a < pointIndices.size(); a++) {
                            cout << pointIndices.at(a) << ", ";
                        }cout << endl << endl;
                    }

                    std::stringstream ss;
                    ss << "data/Tets/MinimalTetraEdge_";
                    for (int st = 0; st < pointIndices.size(); st++) ss << pointIndices[st] << "_";
                    ss << ".vtk";

                    std::vector<int*> tris = getCGALTriangles(pointIndices, ss.str());
                    printTriangles(tris, ss.str());
                    for (int t = 0; t < tris.size(); t++) {
                        int *triangle = tris.at(t);

                        if (isDuplicateTriangle(triangle[0], triangle[1], triangle[2]) == false) {
                            m_Oct->m_Geom->addTriangle(triangle[0], triangle[1], triangle[2]);
                            m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);

                            m_SpecialTriangles.push_back(triangle);
                            cout << "\t\tMinimal TetraEdge Rule Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << endl;
                        }
                    }
                //}
            }
        }
    }
}

void OctreeContouring::testPolygonizationFunction(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[2], float p[3], float q[3], int pq_mat[2]) {
    
}

void OctreeContouring::printOCellGroup(OCell* q1, OCell* q2, OCell* q3, OCell* q4, int v, vval mat[], float p[], float q[], int pq_mat[], std::string outputFile, vector<int> t0, vector<int> t1, vector<int> t2, vector<int> t3) {
    OCell* qs[4] = {q1, q2, q3, q4};
    
    vtkSmartPointer<vtkPolyData> pdata = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> cubesCellArray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkCellArray> trianglesCellArray = vtkSmartPointer<vtkCellArray>::New();
    
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    
    int pointCounter = 0;
    
    // Insert points and then lines for all 4 OCells. 
    for (int j = 0; j < 4; j++) {
        float **qPoints = qs[j]->getCorners();
        for (int i = 0; i < 9; i++) {
            points->InsertNextPoint(qPoints[i][0], qPoints[i][1], qPoints[i][2]);
        }

        vtkSmartPointer<vtkIdList> cellList0 = vtkSmartPointer<vtkIdList>::New();
        cellList0->InsertNextId(pointCounter);
        cellList0->InsertNextId(pointCounter + 1);
        cellList0->InsertNextId(pointCounter + 3);
        cellList0->InsertNextId(pointCounter + 2);
        cellList0->InsertNextId(pointCounter);
        cubesCellArray->InsertNextCell(cellList0);

        vtkSmartPointer<vtkIdList> cellList1 = vtkSmartPointer<vtkIdList>::New();
        cellList1->InsertNextId(pointCounter);
        cellList1->InsertNextId(pointCounter + 1);
        cellList1->InsertNextId(pointCounter + 5);
        cellList1->InsertNextId(pointCounter + 4);
        cellList1->InsertNextId(pointCounter);
        cubesCellArray->InsertNextCell(cellList1);

        vtkSmartPointer<vtkIdList> cellList2 = vtkSmartPointer<vtkIdList>::New();
        cellList2->InsertNextId(pointCounter + 6);
        cellList2->InsertNextId(pointCounter + 7);
        cellList2->InsertNextId(pointCounter + 5);
        cellList2->InsertNextId(pointCounter + 4);
        cellList2->InsertNextId(pointCounter + 6);
        cubesCellArray->InsertNextCell(cellList2);

        vtkSmartPointer<vtkIdList> cellList3 = vtkSmartPointer<vtkIdList>::New();
        cellList3->InsertNextId(pointCounter + 2);
        cellList3->InsertNextId(pointCounter + 3);
        cellList3->InsertNextId(pointCounter + 7);
        cellList3->InsertNextId(pointCounter + 6);
        cellList3->InsertNextId(pointCounter + 2);
        cubesCellArray->InsertNextCell(cellList3);

        vtkSmartPointer<vtkIdList> cellList4 = vtkSmartPointer<vtkIdList>::New();
        cellList4->InsertNextId(pointCounter);
        cellList4->InsertNextId(pointCounter + 2);
        cellList4->InsertNextId(pointCounter + 6);
        cellList4->InsertNextId(pointCounter + 4);
        cellList4->InsertNextId(pointCounter);
        cubesCellArray->InsertNextCell(cellList4);

        vtkSmartPointer<vtkIdList> cellList5 = vtkSmartPointer<vtkIdList>::New();
        cellList5->InsertNextId(pointCounter + 1);
        cellList5->InsertNextId(pointCounter + 3);
        cellList5->InsertNextId(pointCounter + 7);
        cellList5->InsertNextId(pointCounter + 5);
        cellList5->InsertNextId(pointCounter + 1);
        cubesCellArray->InsertNextCell(cellList5);
        
        pointCounter = pointCounter + 9;
    }
    
    // Insert points and lines for each vectors containing triangles. Make sure the vectors are not empty.     
    while (t0.size() > 0) {
        unsigned int *pt_ind = this->m_Oct->m_Geom->getTri(t0.at(t0.size() - 1));
        t0.pop_back();
        
        float *p0 = this->m_Oct->m_Geom->getVert(pt_ind[0]);
        float *p1 = this->m_Oct->m_Geom->getVert(pt_ind[1]);
        float *p2 = this->m_Oct->m_Geom->getVert(pt_ind[2]);
        
        points->InsertNextPoint(p0);
        points->InsertNextPoint(p1);
        points->InsertNextPoint(p2);
        
        vtkSmartPointer<vtkIdList> cList = vtkSmartPointer<vtkIdList>::New();
        cList->InsertNextId(pointCounter);
        cList->InsertNextId(pointCounter + 1);
        cList->InsertNextId(pointCounter + 2);
        
        trianglesCellArray->InsertNextCell(cList);
        
        pointCounter = pointCounter + 3;
    }
    
    while (t1.size() > 0) {
        unsigned int *pt_ind = this->m_Oct->m_Geom->getTri(t1.at(t1.size() - 1));
        t1.pop_back();
        
        float *p0 = this->m_Oct->m_Geom->getVert(pt_ind[0]);
        float *p1 = this->m_Oct->m_Geom->getVert(pt_ind[1]);
        float *p2 = this->m_Oct->m_Geom->getVert(pt_ind[2]);
        
        points->InsertNextPoint(p0);
        points->InsertNextPoint(p1);
        points->InsertNextPoint(p2);
        
        vtkSmartPointer<vtkIdList> cList = vtkSmartPointer<vtkIdList>::New();
        cList->InsertNextId(pointCounter);
        cList->InsertNextId(pointCounter + 1);
        cList->InsertNextId(pointCounter + 2);
        
        trianglesCellArray->InsertNextCell(cList);
        
        pointCounter = pointCounter + 3;
    }
    
    while (t2.size() > 0) {
        unsigned int *pt_ind = this->m_Oct->m_Geom->getTri(t2.at(t2.size() - 1));
        t2.pop_back();
        
        float *p0 = this->m_Oct->m_Geom->getVert(pt_ind[0]);
        float *p1 = this->m_Oct->m_Geom->getVert(pt_ind[1]);
        float *p2 = this->m_Oct->m_Geom->getVert(pt_ind[2]);
        
        points->InsertNextPoint(p0);
        points->InsertNextPoint(p1);
        points->InsertNextPoint(p2);
        
        vtkSmartPointer<vtkIdList> cList = vtkSmartPointer<vtkIdList>::New();
        cList->InsertNextId(pointCounter);
        cList->InsertNextId(pointCounter + 1);
        cList->InsertNextId(pointCounter + 2);
        
        trianglesCellArray->InsertNextCell(cList);
        
        pointCounter = pointCounter + 3;
    }
    
    while (t3.size() > 0) {
        unsigned int *pt_ind = this->m_Oct->m_Geom->getTri(t3.at(t3.size() - 1));
        t3.pop_back();
        
        float *p0 = this->m_Oct->m_Geom->getVert(pt_ind[0]);
        float *p1 = this->m_Oct->m_Geom->getVert(pt_ind[1]);
        float *p2 = this->m_Oct->m_Geom->getVert(pt_ind[2]);
        
        points->InsertNextPoint(p0);
        points->InsertNextPoint(p1);
        points->InsertNextPoint(p2);
        
        vtkSmartPointer<vtkIdList> cList = vtkSmartPointer<vtkIdList>::New();
        cList->InsertNextId(pointCounter);
        cList->InsertNextId(pointCounter + 1);
        cList->InsertNextId(pointCounter + 2);
        
        trianglesCellArray->InsertNextCell(cList);
        
        pointCounter = pointCounter + 3;
    }
    
    
    pdata->SetPoints(points);
    pdata->SetLines(cubesCellArray);
    pdata->SetPolys(trianglesCellArray);
    
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInput(pdata);
    writer->SetFileTypeToASCII();
    writer->SetFileName(outputFile.c_str());
    writer->Write();
}


#ifdef USE_TRADITIONAL_DC
// This copy of edgeProc() is for traditional DC. 
void OctreeContouring::edgeProc(OCell* q1, OCell* q2, OCell* q3, OCell* q4) {
    OCell* qs[4] = {q1, q2, q3, q4}; 
    OCell* param[4];
    
    int edgeArr[4]; 
    int relPos; 
    int j;
        
    if (!qs[0]->hasChildren()  && !qs[1]->hasChildren() && !qs[2]->hasChildren() && !qs[3]->hasChildren()) {
        double val[2]; 
        vval mat[2]; // Probably the material indices. 
        findEdgeValue(qs, val, mat);         
        
        if (mat[0] != mat[1]) {
            int v = 0; 
            v |= (1 << mat[0]);
            v |= (1 << mat[1]);
            
            for (int j = 0; j < 4; j++) {
                vval maskvals[8];
                for (int corner = 0; corner < 8; corner++) {
                    maskvals[corner] = CorrectOctree::getVolumeVal(this->m_Oct->m_CoarseData.mask, (int)qs[j]->m_Ind[0], (int)qs[j]->m_Ind[1], (int)qs[j]->m_Ind[2], corner);
                }
                
                
                
//                // If the OCell is ambiguous. 
//                if ((isAmbiguous(maskvals) == true) && (qs[j]->getHasTetrahedra() == false)) {
//                    qs[j]->setTetrahedralCells(); // Set the tetrahedral cells for the OCell
//                    qs[j]->setIsAmbiguous(true);
////                    cout << "\tAmbiguous cell has index [" << (int)qs[j]->m_Ind[0] << ", " << (int)qs[j]->m_Ind[1] << ", " << (int)qs[j]->m_Ind[2] << "]" << endl;
////                    int *matVals = qs[j]->getMaterialValues();
////                    cout << "\tCorners: [";
////                    for (int mt = 0; mt < 9; mt++) {
////                        cout << matVals[mt] << ", ";
////                    }cout << "]" << endl;
//
//                    std::stringstream ss;
//                    int indoc[] = {(int)qs[j]->m_Ind[0], (int)qs[j]->m_Ind[1], (int)qs[j]->m_Ind[2]};
//                    ss << "/home/trash001/NetBeansProjects/Zhang_HybridDC/Zhang_Hybrid_DualContouring/data/Tets/Amb_Cube_" << indoc[0] << "_" << indoc[1] << "_" << indoc[2] << ".vtk";
//                    qs[j]->vtkFile(ss.str());
//                    
//                    // Check whether the tetrahedral cells have material change edges. 
//                    TCell ** tcs = qs[j]->getTetrahedralCells();
//                    for (int tc = 0; tc < 12; tc++) {
//                        if (tcs[tc]->hasSignChangeEdge()) {
//                            tcs[tc]->setBPoint(this->m_Oct->m_Geom->addVertex(tcs[tc]->getCenter()));
//                            
//                            std::stringstream ss;
//                            int* ind = tcs[tc]->getIndex();
//                            ss << "/home/trash001/NetBeansProjects/Zhang_HybridDC/Zhang_Hybrid_DualContouring/data/Tets/Tetra_" << ind[0] << "_" << ind[1] << "_" << ind[2] << "_" << ind[3] << ".vtk";
//                            tcs[tc]->vtkFile(ss.str());
//                        }
//                    }
//                }
                // If the OCell is unambiguous. 
//                else if (isAmbiguous(maskvals) == false) {
                    qs[j]->setIsAmbiguous(false);
                    
                    if (qs[j]->getBPoint() < 0 && (int)qs[j]->m_Level == m_Oct->m_CoarseData.level) {
                        qs[j]->setBPoint(m_Oct->putMinimizer(qs[j], &m_Oct->m_CoarseData, qs[j]->m_Ind, true));
                        
                                                
                        std::stringstream ss;
                        int indoc[] = {(int)qs[j]->m_Ind[0], (int)qs[j]->m_Ind[1], (int)qs[j]->m_Ind[2]};
                        ss << "/home/trash001/NetBeansProjects/Zhang_HybridDC/Zhang_Hybrid_DualContouring/data/Tets/NORMAL_Cube_" << indoc[0] << "_" << indoc[1] << "_" << indoc[2] << ".vtk";
                        qs[j]->vtkFile(ss.str());
//                        int indoc[] = {(int)qs[j]->m_Ind[0], (int)qs[j]->m_Ind[1], (int)qs[j]->m_Ind[2]};
//                        if (indoc[0] == 6 && indoc[1] == 8 && indoc[2] == 6) {
//                            int pointNo = qs[j]->getBPoint();
//                            float *point = this->m_Oct->m_Geom->getVert(pointNo);
//                            cout << point[0] << ", " << point[1] << ", " << point[2] << endl;
//                        }
                    }
                    
//                    if ((int)qs[j]->m_Ind[0] == 11 && (int)qs[j]->m_Ind[1] == 10 && (int)qs[j]->m_Ind[2] == 5) {
//                        qs[j]->setTetrahedralCells(); // Set the tetrahedral cells for the OCell
//                        qs[j]->setIsAmbiguous(true);
//                    }
//                }
            }

//            // Triangulation when no tetrahedra are present in the 4 OCells. 
//            if ((qs[0]->getIsAmbiguous() == false) && (qs[1]->getIsAmbiguous() == false) && (qs[2]->getIsAmbiguous() == false) && (qs[3]->getIsAmbiguous() == false)) {
                // Normal (Powei's) cell/face generation
                if ((int)qs[0]->m_Level == m_Oct->m_CoarseData.level && (int)qs[1]->m_Level == m_Oct->m_CoarseData.level && (int)qs[2]->m_Level == m_Oct->m_CoarseData.level && (int)qs[3]->m_Level == m_Oct->m_CoarseData.level) {
                    if(val[0] <= val[1]) {
                        // Clockwise
                        // For quads
                        //m_Oct->m_Geom->addQuad(qs[0]->m_BiPoint, qs[1]->m_BiPoint, qs[2]->m_BiPoint, qs[3]->m_BiPoint);
                        //cout << qs[0]->m_BiPoint << ", " << qs[1]->m_BiPoint << ", " << qs[2]->m_BiPoint << ", " << qs[3]->m_BiPoint << endl;

                        // For triangles
                        m_Oct->m_Geom->addTriangle(qs[0]->getBPoint(), qs[1]->getBPoint(), qs[3]->getBPoint());
                        //cout << "\t\tREGULAR DC Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << qs[0]->getBPoint() << ", " << qs[1]->getBPoint() << ", " << qs[3]->getBPoint() << endl;
                        m_Oct->m_Geom->addTriangle(qs[1]->getBPoint(), qs[2]->getBPoint(), qs[3]->getBPoint());
                        //cout << "\t\tREGULAR DC Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << qs[1]->getBPoint() << ", " << qs[2]->getBPoint() << ", " << qs[3]->getBPoint() << endl;
                        
                        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
                        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
                    }
                    else {
                        // Counter-clockwise
                        // For quads
                        //m_Oct->m_Geom->addQuad(qs[3]->m_BiPoint, qs[2]->m_BiPoint, qs[1]->m_BiPoint, qs[0]->m_BiPoint);
                        //cout << "    " << qs[0]->m_BiPoint << ", " << qs[1]->m_BiPoint << ", " << qs[2]->m_BiPoint << ", " << qs[3]->m_BiPoint << endl;

                        // For triangles
                        m_Oct->m_Geom->addTriangle(qs[3]->getBPoint(), qs[2]->getBPoint(), qs[1]->getBPoint());
                        //cout << "\t\tREGULAR DC Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << qs[3]->getBPoint() << ", " << qs[2]->getBPoint() << ", " << qs[1]->getBPoint() << endl;
                        m_Oct->m_Geom->addTriangle(qs[3]->getBPoint(), qs[1]->getBPoint(), qs[0]->getBPoint());
                        //cout << "\t\tREGULAR DC Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << qs[3]->getBPoint() << ", " << qs[1]->getBPoint() << ", " << qs[0]->getBPoint() << endl;
                        
                        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
                        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
                    }
                }
//            }
            
            // When tetrahedra are present amongst the 4 OCells. 
//            else {
//                cout << "Cubes: ";
//                for (int as = 0; as < 4; as++) {
//                    //std::stringstream ss;
//                    int indoc[] = {(int)qs[as]->m_Ind[0], (int)qs[as]->m_Ind[1], (int)qs[as]->m_Ind[2]};
//                    //ss << "/home/trash001/NetBeansProjects/Zhang_HybridDC/Zhang_Hybrid_DualContouring/data/Tets/XXXCube_" << indoc[0] << "_" << indoc[1] << "_" << indoc[2] << ".vtk";
//                    //qs[as]->vtkFile(ss.str());
//                    
//                    cout << "[" << indoc[0] << ", " << indoc[1] << ", " << indoc[2] << "]   ";
//                    
////                    if (indoc[0] == 4 && indoc[1] == 10 & indoc[2] == 14) {
////                        cout << "\n\nBiPoint: " << qs[as]->getBPoint() << endl << endl;
////                    }
//                } cout << endl;
//                
//          
//                
//                // Determine the coordinates and material indices of the shared (minimal) edge. 
//                float p[3], q[3]; // The two ends of the shared (minimal) edge. 
//                int pq_mat[2]; // The material indices of the two edge end points. 
//                getSharedEdge(qs, p, q, pq_mat);
//                
////                cout << "p: " << p[0] << ", " << p[1] << ", " << p[2] << endl;
////                cout << "q: " << q[0] << ", " << q[1] << ", " << q[2] << endl;
////                cout << "p has material index " << pq_mat[0] << endl;
////                cout << "q has material index " << pq_mat[1] << endl;
//                
//                
//                // If one end point of the shared edge is in the background. 
//                if ((pq_mat[0] == 0 && pq_mat[1] != 0) || (pq_mat[0] != 0 && pq_mat[1] == 0)) {
//                    //cout << "\t\tMaterial change with background" << endl;
//                    //materialChangeEdge_WithBackgound(q1, q2, q3, q4, v, mat, p, q, pq_mat);
//                    //materialChangeEdge_WithBackgoundZhang(q1, q2, q3, q4, v, mat, p, q, pq_mat);
//                    //minimalEdgeTriangulation_WithAmbiguousOCell(q1, q2, q3, q4, v, mat, p, q, pq_mat);
//                }
//                if (pq_mat[0] > 1 && pq_mat[1] > 1 && pq_mat[0] != pq_mat[1]) {
//                    //materialChangeEdge_WithOtherMaterial();
//                    //cout << "\t\t\tMaterial change with other materials" << endl;
//                }
//                if (pq_mat[0] > 1 && pq_mat[1] > 1 && pq_mat[0] == pq_mat[1]) {
//                    //cout << "\t\t\t\tInterior edge" << endl;
//                    //interiorEdge();
//                }
//                
//                vector<int> tris0, tris1, tris2;
//                
//                tris0 = minimalEdgeTriangulation_WithAmbiguousOCell(q1, q2, q3, q4, v, mat, p, q, pq_mat);                
//                tris1 = faceRuleImplementation(q1, q2, q3, q4, v, mat, p, q, pq_mat);
//                //minimalTetraEdge(q1, q2, q3, q4, v, mat, p, q, pq_mat);
//                tris2 = individualAmbiguousCellRule(q1, q2, q3, q4, p, q, pq_mat);
//                
//                                
//                //if (tris0.size() > 0) {
//                    cout << "minimalEdgeTriangulation_WithAmbiguousOCell Rule created triangles: ";
//                    for (int i = 0; i < tris0.size(); i++) {
//                        cout << tris0.at(i) << ", ";
//                    }cout << endl;
//                //}
//                
//                //if (tris1.size() > 0) {
//                    cout << "faceRule Rule created triangles: ";
//                    for (int i = 0; i < tris1.size(); i++) {
//                        cout << tris1.at(i) << ", ";
//                    }cout << endl;
//                //}
//                
//                //if (tris2.size() > 0) {
//
//                    cout << "individualAmbiguousOCell Rule created triangles: ";
//                    for (int i = 0; i < tris2.size(); i++) {
//                        cout << tris2.at(i) << ", ";
//                    }cout << endl;
//                //}
//                
//                std::stringstream groupfilename;
//                groupfilename << "data/Tets/Group_" << 
//                        q1->m_Ind[0] << "-" << q1->m_Ind[1] << "-" << q1->m_Ind[2] << "_" << 
//                        q2->m_Ind[0] << "-" << q2->m_Ind[1] << "-" << q2->m_Ind[2] << "_" << 
//                        q3->m_Ind[0] << "-" << q3->m_Ind[1] << "-" << q3->m_Ind[2] << "_" << 
//                        q4->m_Ind[0] << "-" << q4->m_Ind[1] << "-" << q4->m_Ind[2] << ".vtk";
//                
//                
//                printOCellGroup(q1, q2, q3, q4, v, mat, p, q, pq_mat, groupfilename.str(), tris0, tris1, tris2);
//                
//                cout << "-------------------------------------------------------------" << endl << endl;
//            }
        }
    }
    else {
        if(qs[0] == qs[1]){
            relPos = findRelPositionOf2Cells(qs[0], qs[2]); 
            /* pay attention here */
            switch(relPos){
                case XPLUS: relPos = YPLUS; break; 
                case YPLUS: relPos = ZPLUS; break; 
                case ZPLUS: relPos = XPLUS; break;
            }
        }
        else
            relPos = findRelPositionOf2Cells(qs[0], qs[1]);
        
        switch(relPos){
        case XPLUS: 
            edgeArr[0] = 4;
            edgeArr[1] = 1;
            edgeArr[2] = 10;
            edgeArr[3] = 7;            
            break; 
        case YPLUS: 
            edgeArr[0] = 6; 
            edgeArr[1] = 3; 
            edgeArr[2] = 0; 
            edgeArr[3] = 9;
            break; 
        case ZPLUS: 
            edgeArr[0] = 11;
            edgeArr[1] = 8;
            edgeArr[2] = 5;
            edgeArr[3] = 2;
            break; 
        } 
    
        for(j=0; j<4; j++){
            if(qs[j]->hasChildren())
                param[j] = qs[j]->m_Children[ EV[edgeArr[j]][0] ]; 
            else
                param[j] = qs[j]; 
        }
            
        edgeProc(param[0], param[1], param[2], param[3]); 

        for(j=0; j<4; j++){
            if(qs[j]->hasChildren())
                param[j] = qs[j]->m_Children[ EV[edgeArr[j]][1] ]; 
            else
                param[j] = qs[j]; 
        }
    
        //second recursive call to edgeProc        
        edgeProc(param[0], param[1], param[2], param[3]);     
    }
}

#else
 // This copy of edgeProc() is for modified 2-manifold DC
void OctreeContouring::edgeProc(OCell* q1, OCell* q2, OCell* q3, OCell* q4) {
    OCell* qs[4] = {q1, q2, q3, q4}; 
    OCell* param[4];
    
    int edgeArr[4]; 
    int relPos; 
    int j;
    
    
    if (!qs[0]->hasChildren()  && !qs[1]->hasChildren() && !qs[2]->hasChildren() && !qs[3]->hasChildren()) {
        double val[2]; 
        vval mat[2]; // Probably the material indices. 
        findEdgeValue(qs, val, mat);         
        
        if (mat[0] != mat[1]) {
            int v = 0; 
            v |= (1 << mat[0]);
            v |= (1 << mat[1]);
            
            for (int j = 0; j < 4; j++) {
                vval maskvals[8];
                for (int corner = 0; corner < 8; corner++) {
                    maskvals[corner] = CorrectOctree::getVolumeVal(this->m_Oct->m_CoarseData.mask, (int)qs[j]->m_Ind[0], (int)qs[j]->m_Ind[1], (int)qs[j]->m_Ind[2], corner);
                }
                
                
//                if ((int)qs[j]->m_Ind[0] == 4 && (int)qs[j]->m_Ind[1] == 4 && (int)qs[j]->m_Ind[2] == 4) { // || (int)m_Ind[2] == 1 || (int)m_Ind[2] == 2 || (int)m_Ind[2] == 3)) {
//                    maskvals[6] = 2;
//                }
                
                // If the OCell is ambiguous. 
                if ((isAmbiguous(maskvals) == true) && (qs[j]->getHasTetrahedra() == false)) {
                    qs[j]->setTetrahedralCells(); // Set the tetrahedral cells for the OCell
                    qs[j]->setIsAmbiguous(true);
//                    cout << "\tAmbiguous cell has index [" << (int)qs[j]->m_Ind[0] << ", " << (int)qs[j]->m_Ind[1] << ", " << (int)qs[j]->m_Ind[2] << "]" << endl;
//                    int *matVals = qs[j]->getMaterialValues();
//                    cout << "\tCorners: [";
//                    for (int mt = 0; mt < 9; mt++) {
//                        cout << matVals[mt] << ", ";
//                    }cout << "]" << endl;

//                    if (isCase13Ambiguous(maskvals)) {
//                        cout << "Case 13: " << (int)qs[j]->m_Ind[0] << ", " << (int)qs[j]->m_Ind[1] << ", " << (int)qs[j]->m_Ind[2] << endl;
//                        std::stringstream ss;
//                        ss << "/home/trash001/NetBeansProjects/Zhang_HybridDC/Zhang_Hybrid_DualContouring/data/Tets/Case13_" << (int)qs[j]->m_Ind[0] << "_" << (int)qs[j]->m_Ind[1] << "_" << (int)qs[j]->m_Ind[2] << ".vtk";
//                        qs[j]->vtkFile(ss.str());
//                    }
                    
                    std::stringstream ss;
                    int indoc[] = {(int)qs[j]->m_Ind[0], (int)qs[j]->m_Ind[1], (int)qs[j]->m_Ind[2]};
                    ss << "data\\Tets\\Amb_Cube_" << indoc[0] << "_" << indoc[1] << "_" << indoc[2] << ".vtk";
                    qs[j]->vtkFile(ss.str());
                    
                    // Check whether the tetrahedral cells have material change edges. 
                    TCell ** tcs = qs[j]->getTetrahedralCells();
                    for (int tc = 0; tc < 12; tc++) {
                        if (tcs[tc]->hasSignChangeEdge()) {
                            tcs[tc]->setBPoint(this->m_Oct->m_Geom->addVertex(tcs[tc]->getCenter()));
                            
                            std::stringstream ss;
                            int* ind = tcs[tc]->getIndex();
                            ss << "data\\Tets\\Tetra_" << ind[0] << "_" << ind[1] << "_" << ind[2] << "_" << ind[3] << ".vtk";
                            tcs[tc]->vtkFile(ss.str());
                        }
                    }
                    
//                    if ((int)qs[j]->m_Ind[0] == 4 && (int)qs[j]->m_Ind[1] == 4 && (int)qs[j]->m_Ind[2] == 4) {
//                        int *mv = qs[j]->getMaterialValues();
//                        cout << mv[0] << ", " << mv[1] << ", " << mv[2] << ", " << mv[3] << ", " << mv[4] << ", " << mv[5] << ", " << mv[6] << ", " << mv[7] << ", " << mv[8] << endl;
//                        
//                        mv[8] = 0;
//                        qs[j]->setMaterialValues(mv);
//                    }
                }
                // If the OCell is unambiguous. 
                else if (isAmbiguous(maskvals) == false) {
//                    if ((int)qs[j]->m_Ind[0] == 3 && (int)qs[j]->m_Ind[1] == 2 && (int)qs[j]->m_Ind[2] == 2) {
//                        cout << "here ";
//                        cout << "bpt = " << qs[j]->getBPoint() << endl;
//                    }
                    
                    qs[j]->setIsAmbiguous(false);
                    
                    if (qs[j]->getBPoint() < 0 && (int)qs[j]->m_Level == m_Oct->m_CoarseData.level) {
                        int bpt = m_Oct->putMinimizer(qs[j], &m_Oct->m_CoarseData, qs[j]->m_Ind, true);
                        qs[j]->setBPoint(bpt);
                        
//                        std::stringstream ss;
//                        int indoc[] = {(int)qs[j]->m_Ind[0], (int)qs[j]->m_Ind[1], (int)qs[j]->m_Ind[2]};
//                        ss << "/home/trash001/NetBeansProjects/Zhang_HybridDC/Zhang_Hybrid_DualContouring/data/Tets/NORMAL_Cube_" << indoc[0] << "_" << indoc[1] << "_" << indoc[2] << ".vtk";
//                        qs[j]->vtkFile(ss.str());
//                        int indoc[] = {(int)qs[j]->m_Ind[0], (int)qs[j]->m_Ind[1], (int)qs[j]->m_Ind[2]};
//                        if (indoc[0] == 6 && indoc[1] == 8 && indoc[2] == 6) {
//                            int pointNo = qs[j]->getBPoint();
//                            float *point = this->m_Oct->m_Geom->getVert(pointNo);
//                            cout << point[0] << ", " << point[1] << ", " << point[2] << endl;
//                        }
                    }
                }
            }

//            // Triangulation when no tetrahedra are present in the 4 OCells. 
            if ((qs[0]->getIsAmbiguous() == false) && (qs[1]->getIsAmbiguous() == false) && (qs[2]->getIsAmbiguous() == false) && (qs[3]->getIsAmbiguous() == false)) {
                // Normal (Powei's) cell/face generation
                if ((int)qs[0]->m_Level == m_Oct->m_CoarseData.level && (int)qs[1]->m_Level == m_Oct->m_CoarseData.level && (int)qs[2]->m_Level == m_Oct->m_CoarseData.level && (int)qs[3]->m_Level == m_Oct->m_CoarseData.level) {
                    if(val[0] <= val[1]) {
                        // Clockwise
                        // For quads
                        //m_Oct->m_Geom->addQuad(qs[0]->m_BiPoint, qs[1]->m_BiPoint, qs[2]->m_BiPoint, qs[3]->m_BiPoint);
                        //cout << qs[0]->m_BiPoint << ", " << qs[1]->m_BiPoint << ", " << qs[2]->m_BiPoint << ", " << qs[3]->m_BiPoint << endl;

                        // For triangles
                        m_Oct->m_Geom->addTriangle(qs[0]->getBPoint(), qs[1]->getBPoint(), qs[3]->getBPoint());
                        //cout << "\t\tREGULAR DC Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << qs[0]->getBPoint() << ", " << qs[1]->getBPoint() << ", " << qs[3]->getBPoint() << endl;
                        m_Oct->m_Geom->addTriangle(qs[1]->getBPoint(), qs[2]->getBPoint(), qs[3]->getBPoint());
                        //cout << "\t\tREGULAR DC Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << qs[1]->getBPoint() << ", " << qs[2]->getBPoint() << ", " << qs[3]->getBPoint() << endl;
                        
                        double ratio1 = computeRadiusRatio(qs[0]->getBPoint(), qs[1]->getBPoint(), qs[3]->getBPoint());
                        double ratio2 = computeRadiusRatio(qs[1]->getBPoint(), qs[2]->getBPoint(), qs[3]->getBPoint());
                        if (ratio1 < 0.2) {
                            //cout << "\t\tWarning: Poor quality triangle created by RegularDCPolygonizing: Triangle made of points " << qs[0]->getBPoint() << ", " << qs[1]->getBPoint() << " and " << qs[3]->getBPoint() << endl;
                            regularDCPolygonizerCounter++;
                        }
                        if (ratio2 < 0.2) {
                            //cout << "\t\tWarning: Poor quality triangle created by RegularDCPolygonizing: Triangle made of points " << qs[1]->getBPoint() << ", " << qs[2]->getBPoint() << " and " << qs[3]->getBPoint() << endl;
                            regularDCPolygonizerCounter++;
                        }
                        
                        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
                        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
                    }
                    else {
                        // Counter-clockwise
                        // For quads
                        //m_Oct->m_Geom->addQuad(qs[3]->m_BiPoint, qs[2]->m_BiPoint, qs[1]->m_BiPoint, qs[0]->m_BiPoint);
                        //cout << "    " << qs[0]->m_BiPoint << ", " << qs[1]->m_BiPoint << ", " << qs[2]->m_BiPoint << ", " << qs[3]->m_BiPoint << endl;

                        // For triangles
                        m_Oct->m_Geom->addTriangle(qs[3]->getBPoint(), qs[2]->getBPoint(), qs[1]->getBPoint());
                        //cout << "\t\tREGULAR DC Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << qs[3]->getBPoint() << ", " << qs[2]->getBPoint() << ", " << qs[1]->getBPoint() << endl;
                        m_Oct->m_Geom->addTriangle(qs[3]->getBPoint(), qs[1]->getBPoint(), qs[0]->getBPoint());
                        //cout << "\t\tREGULAR DC Triangle " << (m_Oct->m_Geom->getNumTris() - 1) <<  ": " << qs[3]->getBPoint() << ", " << qs[1]->getBPoint() << ", " << qs[0]->getBPoint() << endl;
                        
                        double ratio1 = computeRadiusRatio(qs[3]->getBPoint(), qs[2]->getBPoint(), qs[1]->getBPoint());
                        double ratio2 = computeRadiusRatio(qs[3]->getBPoint(), qs[1]->getBPoint(), qs[0]->getBPoint());
                        if (ratio1 < 0.2) {
                            //cout << "\t\tWarning: Poor quality triangle created by RegularDCPolygonizing: Triangle made of points " << qs[3]->getBPoint() << ", " << qs[2]->getBPoint() << " and " << qs[1]->getBPoint() << endl;
                            regularDCPolygonizerCounter++;
                        }
                        if (ratio2 < 0.2) {
                            //cout << "\t\tWarning: Poor quality triangle created by RegularDCPolygonizing: Triangle made of points " << qs[3]->getBPoint() << ", " << qs[1]->getBPoint() << " and " << qs[0]->getBPoint() << endl;
                            regularDCPolygonizerCounter++;
                        }
                        
                        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
                        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
                    }
                    
//                    vector<int> cell; 
//                    cell.push_back(qs[0]->getBPoint());
//                    cell.push_back(qs[1]->getBPoint());
//                    cell.push_back(qs[2]->getBPoint());
//                    cell.push_back(qs[3]->getBPoint());
//                    
//                    std::stringstream ss1;
//                    ss1 << "data/Tets/FaceDiagonalRuleTri_" << qs[0]->getBPoint() << "_" << qs[1]->getBPoint() << "_" << qs[2]->getBPoint() << "_" << qs[3]->getBPoint() << ".vtk";
//                    std::vector<int*> tris = getCGALTriangles(cell, ss1.str());
//                    
//                    for (int a = 0; a < tris.size(); a++) {
//                        int *tri0 = tris.at(a);
//                        
//                        m_Oct->m_Geom->addTriangle(tri0[0], tri0[1], tri0[2]);
//                        m_Oct->m_GeomColors->addTriangle(v, (int)mat[0], (int)mat[1]);
//                        
//                        double ratio = computeRadiusRatio(tri0[0], tri0[1], tri0[2]);
//                        if (ratio < 0.2) {
//                            cout << "\t\tWarning: Poor quality triangle created by RegularDCPolygonizing: Triangle made of points " << tri0[0] << ", " << tri0[1] << " and " << tri0[2] << endl;
//                            regularDCPolygonizerCounter++;
//                        }
//                    }
                }
            }
            
            // When tetrahedra are present amongst the 4 OCells. 
            else {
//                cout << "Cubes: ";
//                for (int as = 0; as < 4; as++) {
//                    //std::stringstream ss;
//                    int indoc[] = {(int)qs[as]->m_Ind[0], (int)qs[as]->m_Ind[1], (int)qs[as]->m_Ind[2]};
//                    //ss << "/home/trash001/NetBeansProjects/Zhang_HybridDC/Zhang_Hybrid_DualContouring/data/Tets/XXXCube_" << indoc[0] << "_" << indoc[1] << "_" << indoc[2] << ".vtk";
//                    //qs[as]->vtkFile(ss.str());
//                    
//                    cout << "[" << indoc[0] << ", " << indoc[1] << ", " << indoc[2] << "]: " << qs[as]->getBPoint() << "   ";
//                    
////                    if (indoc[0] == 4 && indoc[1] == 10 & indoc[2] == 14) {
////                        cout << "\n\nBiPoint: " << qs[as]->getBPoint() << endl << endl;
////                    }
//                } cout << endl;
                
          
                
                // Determine the coordinates and material indices of the shared (minimal) edge. 
                float p[3], q[3]; // The two ends of the shared (minimal) edge. 
                int pq_mat[2]; // The material indices of the two edge end points. 
                getSharedEdge(qs, p, q, pq_mat);
                
//                cout << "p: " << p[0] << ", " << p[1] << ", " << p[2] << endl;
//                cout << "q: " << q[0] << ", " << q[1] << ", " << q[2] << endl;
//                cout << "p has material index " << pq_mat[0] << endl;
//                cout << "q has material index " << pq_mat[1] << endl;
                
                
                // If one end point of the shared edge is in the background. 
//                if ((pq_mat[0] == 0 && pq_mat[1] != 0) || (pq_mat[0] != 0 && pq_mat[1] == 0)) {
//                    //cout << "\t\tMaterial change with background" << endl;
//                    //materialChangeEdge_WithBackgound(q1, q2, q3, q4, v, mat, p, q, pq_mat);
//                    //materialChangeEdge_WithBackgoundZhang(q1, q2, q3, q4, v, mat, p, q, pq_mat);
//                    //minimalEdgeTriangulation_WithAmbiguousOCell(q1, q2, q3, q4, v, mat, p, q, pq_mat);
//                }
//                if (pq_mat[0] > 1 && pq_mat[1] > 1 && pq_mat[0] != pq_mat[1]) {
//                    //materialChangeEdge_WithOtherMaterial();
//                    //cout << "\t\t\tMaterial change with other materials" << endl;
//                }
//                if (pq_mat[0] > 1 && pq_mat[1] > 1 && pq_mat[0] == pq_mat[1]) {
//                    //cout << "\t\t\t\tInterior edge" << endl;
//                    //interiorEdge();
//                }
                
                vector<int> tris0, tris1, tris2, tris3;
                
                tris0 = minimalEdgeTriangulation_WithAmbiguousOCell(q1, q2, q3, q4, v, mat, p, q, pq_mat);                
                tris1 = FaceDiagonalRule_Ambiguous_Unambiguous(q1, q2, q3, q4);
                //minimalTetraEdge(q1, q2, q3, q4, v, mat, p, q, pq_mat);
                tris2 = InteriorEdgeRule(q1, q2, q3, q4);
                tris3 = FaceDiagonalRule_Ambiguous_Ambiguous(q1, q2, q3, q4);
                                
//                //if (tris0.size() > 0) {
//                    cout << "minimalEdgeTriangulation_WithAmbiguousOCell Rule created triangles: ";
//                    for (int i = 0; i < tris0.size(); i++) {
//                        cout << tris0.at(i) << ", ";
//                    }cout << endl;
//                //}
//                
//                //if (tris1.size() > 0) {
//                    cout << "faceRule Rule created triangles: ";
//                    for (int i = 0; i < tris1.size(); i++) {
//                        cout << tris1.at(i) << ", ";
//                    }cout << endl;
//                //}
//                
//                //if (tris2.size() > 0) {
//
//                    cout << "individualAmbiguousOCell Rule created triangles: ";
//                    for (int i = 0; i < tris2.size(); i++) {
//                        cout << tris2.at(i) << ", ";
//                    }cout << endl;
//                //}
//                
//                    cout << "FaceDiagonalRule created triangles: ";
//                    for (int i = 0; i < tris3.size(); i++) {
//                        cout << tris3.at(i) << ", ";
//                    }cout << endl;
//                    
                std::stringstream groupfilename;
                groupfilename << "data/Tets/Group_" << 
                        q1->m_Ind[0] << "-" << q1->m_Ind[1] << "-" << q1->m_Ind[2] << "_" << 
                        q2->m_Ind[0] << "-" << q2->m_Ind[1] << "-" << q2->m_Ind[2] << "_" << 
                        q3->m_Ind[0] << "-" << q3->m_Ind[1] << "-" << q3->m_Ind[2] << "_" << 
                        q4->m_Ind[0] << "-" << q4->m_Ind[1] << "-" << q4->m_Ind[2] << ".vtk";
//                
//                
                printOCellGroup(q1, q2, q3, q4, v, mat, p, q, pq_mat, groupfilename.str(), tris0, tris1, tris2, tris3);
//                
//                cout << "-----------------------------------------------------------------------------------------------" << endl << endl;
//                
            }
        }
    }
    else {
        if(qs[0] == qs[1]){
            relPos = findRelPositionOf2Cells(qs[0], qs[2]); 
            /* pay attention here */
            switch(relPos){
                case XPLUS: relPos = YPLUS; break; 
                case YPLUS: relPos = ZPLUS; break; 
                case ZPLUS: relPos = XPLUS; break;
            }
        }
        else
            relPos = findRelPositionOf2Cells(qs[0], qs[1]);
        
        switch(relPos){
        case XPLUS: 
            edgeArr[0] = 4;
            edgeArr[1] = 1;
            edgeArr[2] = 10;
            edgeArr[3] = 7;            
            break; 
        case YPLUS: 
            edgeArr[0] = 6; 
            edgeArr[1] = 3; 
            edgeArr[2] = 0; 
            edgeArr[3] = 9;
            break; 
        case ZPLUS: 
            edgeArr[0] = 11;
            edgeArr[1] = 8;
            edgeArr[2] = 5;
            edgeArr[3] = 2;
            break; 
        } 
    
        for(j=0; j<4; j++){
            if(qs[j]->hasChildren())
                param[j] = qs[j]->m_Children[ EV[edgeArr[j]][0] ]; 
            else
                param[j] = qs[j]; 
        }
            
        edgeProc(param[0], param[1], param[2], param[3]); 

        for(j=0; j<4; j++){
            if(qs[j]->hasChildren())
                param[j] = qs[j]->m_Children[ EV[edgeArr[j]][1] ]; 
            else
                param[j] = qs[j]; 
        }
    
        //second recursive call to edgeProc        
        edgeProc(param[0], param[1], param[2], param[3]);     
    }
}
#endif


/* Given four neighboring cells, find the value of the intersecting smallest edge */
void OctreeContouring::findEdgeValue(OCell* qs[4], double val[2], vval mat[2]){
    int j; 

    int sm; /* index of the smallest node */
    int relPos;
    uint ind1, ind2; 
    uint gridSize, gridSize_2; 
    uint edgeNum = 0; 

    OCell* tmp; 
  
    sm  = 0; 
    ind1 = ind2 = 0; 

    for(j=1; j<4; j++)
        if(qs[sm]->m_Level < qs[j]->m_Level)
            sm = j; 

    if(qs[0] == qs[1]){
        relPos = findRelPositionOf2Cells(qs[0], qs[2]); 
        /* pay attention here */
        switch(relPos){
        case XPLUS: relPos = YPLUS; break; 
        case YPLUS: relPos = ZPLUS; break; 
        case ZPLUS: relPos = XPLUS; break; 		  
        }
    }
    else
        relPos = findRelPositionOf2Cells(qs[0], qs[1]);

    tmp = qs[sm]; 


    gridSize = m_Oct->m_CoarseData.gridSize; 
    gridSize_2 = gridSize * gridSize; 

    switch(relPos){
    case XPLUS: 
        switch(sm){
        case 0: edgeNum = 4; break; 
        case 1: edgeNum = 1; break; 
        case 2: edgeNum = 10; break; 
        case 3: edgeNum = 7; break; 
        }
        break;
    case YPLUS:
        switch(sm){
        case 0: edgeNum = 6; break; 
        case 1: edgeNum = 3; break; 
        case 2: edgeNum = 0; break; 
        case 3: edgeNum = 9; break; 
        }    
        break ;
    case ZPLUS: 
        switch(sm){
        case 0: edgeNum = 11; break; 
        case 1: edgeNum = 8; break; 
        case 2: edgeNum = 5; break; 
        case 3: edgeNum = 2; break; 
        }
        break;    
    }
    
    ind1 += (EV[edgeNum][0] & XPLUS) ? (tmp->m_Ind[0]+1) : tmp->m_Ind[0]; 
    ind1 += (EV[edgeNum][0] & YPLUS) ? (tmp->m_Ind[1]+1) * gridSize : tmp->m_Ind[1] * gridSize; 
    ind1 += (EV[edgeNum][0] & ZPLUS) ? (tmp->m_Ind[2]+1) * gridSize_2 : tmp->m_Ind[2] * gridSize_2; 

    ind2 += (EV[edgeNum][1] & XPLUS) ? (tmp->m_Ind[0]+1) : tmp->m_Ind[0]; 
    ind2 += (EV[edgeNum][1] & YPLUS) ? (tmp->m_Ind[1]+1) * gridSize : tmp->m_Ind[1] * gridSize; 
    ind2 += (EV[edgeNum][1] & ZPLUS) ? (tmp->m_Ind[2]+1) * gridSize_2 : tmp->m_Ind[2] * gridSize_2; 
  
	val[0] = m_Oct->m_CoarseData.funcVals[ind1];
	val[1] = m_Oct->m_CoarseData.funcVals[ind2];   

	mat[0] = (*(m_Oct->m_CoarseData.mask))(ind1);
	mat[1] = (*(m_Oct->m_CoarseData.mask))(ind2);
}


void OctreeContouring::findEdgeValues(uint cellInd, double values[24], int level){
    uint xplus, yplus, zplus; 
    int j, k; 
    uint tmpInd; 

	xplus = 1; 
	yplus = m_Oct->m_CoarseData.gridSize; 
	zplus = yplus  * yplus; 
 
	for(j=0; j<12; j++){
		for(k=0; k<2; k++){
			tmpInd = cellInd; 
			if(EV[j][k] & 1)
				tmpInd += xplus; 
			if(EV[j][k] & 2)
				tmpInd += yplus; 
			if(EV[j][k] & 4)
				tmpInd += zplus;       
			values[j*2+k] = m_Oct->m_CoarseData.funcVals[tmpInd]; 
		}
	}
}



/************ end octreeContouring **************/ 


