// header for MRC 

#ifndef MRC_H
#define MRC_H

/* for MRC header */
#define MRC_LABEL_SIZE      80   
#define MRC_USER            25
#define MRC_NUM_LABELS      10

/*  different modes supported by the MRC format. */
#define MODE_char           0
#define MODE_short          1
#define MODE_float          2
#define MODE_short_COMPLEX  3  
#define MODE_float_COMPLEX  4

struct mrcH
{
	// num of columns, rows, sections in map 
	// columns are fastest changing, sections slowest changing
	int	    nc;
	int     nr;
	int     ns;

	int     mode;       // modes supported by MRC format

    int     ncstart;    // first column, row, section in map
    int     nrstart;    // default 0,0,0
    int     nsstart;            

    int     mx;                 /* Number of intervals along X. */
    int     my;                 /* Number of intervals along Y. */
    int     mz;                 /* Number of intervals along Z. */

    float   xlen;               /* Cell dimensions (Angstroms). */
    float   ylen;               /* Cell dimensions (Angstroms). */
    float   zlen;               /* Cell dimensions (Angstroms). */

    float   alpha;              /* Cell angles (Degrees). */
    float   beta;               /* Cell angles (Degrees). */
    float   gamma;              /* Cell angles (Degrees). */

    int     mapc;               /* Which axis corresponds to Columns.  */
                                /* (1,2,3 for X,Y,Z.                   */
    int     mapr;               /* Which axis corresponds to Rows.     */
                                /* (1,2,3 for X,Y,Z.                   */
    int     maps;               /* Which axis corresponds to Sections. */
                                /* (1,2,3 for X,Y,Z.                   */

    float   amin;               /* Minimum density value. */
    float   amax;               /* Maximum density value. */
    float   amean;              /* Mean density value.    */

    int     ispg;               /* Space group number (0 for images). */

    int     nsymbt;             /* Number of chars used for storing   */
                                /* symmetry operators. */


    int     user[MRC_USER];  
	
    float   xorigin;            /* X origin. */
    float   yorigin;            /* Y origin. */
    float   zorigin;            /* Z origin. */

    char    map[4];				/* constant string "MAP "  */
    int     machinestamp;	    /* machine stamp in CCP4 convention: big endian=0x11110000 little endian=0x4444000 */

    float   rms;                /* 	rms deviation of map from mean density */

    int     nlabl;              /* Number of labels being used. */
                                /* 10 text labels of 80 characters each. */
    char    labels[MRC_NUM_LABELS][MRC_LABEL_SIZE];
    
};
#endif
