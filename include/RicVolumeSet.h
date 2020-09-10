// ------------------------ RicVolumeSet.h --------------------------------
/*! \file
Header file for the RicVolumeSet class. The class can open a variety
of different volume file formats. The data is read as a set of
RicVolume structures.

Copyright (C) 2007 by Bill Rogers - Research Imaging Center - UTHSCSA
rogers@uthscsa.edu     
 */
 
#ifndef _RICVOLUMESET_H
#define _RICVOLUMESET_H

using namespace std;

#include <cstdlib>
#include <iostream>
#include <string>
#include "RicUtil.h"
#include "RicVolume.h"
#include "nifti1_io.h"

/// File type
#define RIC_NO_FILE 0
#define RIC_NEMA_FILE 1
#define RIC_NIFTI_FILE 2
#define RIC_ANALYZE_FILE 3
#define RIC_GIS_FILE 4

// Data type
#define RIC_INTEGER 1	// 16 bit signed integer
#define RIC_FLOAT	2	// 32 bit float

/// Units 
#define RIC_UNITS_UNKNOWN 0
#define RIC_UNITS_METER   1
#define RIC_UNITS_MM      2
#define RIC_UNITS_MICRON  3
#define RIC_UNITS_SEC     8
#define RIC_UNITS_MSEC   16
#define RIC_UNITS_USEC   24

/*!
Header file for the RicVolumeSet class. The class can open a variety
of different volume file formats. The data is read as a set of
RicVolume structures. The orientation of the source file is converted to
the internal class format which is the NEMA equivalent of XYZ---.

Copyright (C) 2007 by Bill Rogers - Research Imaging Center - UTHSCSA
rogers@uthscsa.edu     
 */

//Contains Quateronion data for orientation for NIFTI_1 files Edit: Brian Donohue
typedef struct NIFTI_Orient {

	//Quateronion parameters
	float b;
	float c;
	float d;
	//Quateronion offsets
	 float qx;
	 float qy;
	 float qz;
	 //Voxel dimensions offsets
	  float dx;
	  float dy;
	  float dz;
	  //Quateronion scale factor
	   float qfac;
           mat44 qto_xyz;
           mat44 qto_ijk;
           mat44 sto_xyz;
           mat44 sto_ijk;
	   bool emptyStatus;
          NIFTI_Orient& operator=(const NIFTI_Orient& other){
                this->b = other.b;
                 this->c = other.c;
                 this->d = other.d;
                 this->qx = other.qx;
                 this->qy = other.qy;
                 this->qz = other.qz;
                 this->dx = other.dx;
                  this->dy = other.dy;
                  this->dz = other.dz;
                   this->qfac = other.qfac;
                 this->emptyStatus  = other.emptyStatus;
                  for(int i =0; i < 4; i++) for(int r=0; r< 4; r++){

                    this->qto_xyz.m[i][r] = other.qto_xyz.m[i][r];
                    this->sto_xyz.m[i][r] = other.sto_xyz.m[i][r];
                     this->qto_ijk.m[i][r] = other.qto_ijk.m[i][r];
                    this->sto_ijk.m[i][r] = other.sto_ijk.m[i][r];
                   }
            return *this;
	 }
};
class RicVolumeSet
{
	friend class RicVolume;
	
	public:
		
	RicVolume	*VolSet;///< volume structure set
	
	int		nvol;		///< number of volumes in the file
	int		nz;			///< number of data slices in volume
	int		ny;			///< number of rows per slice
	int		nx;			///< number of columns per row
	
	public:
	// stuff read in from volume files
	NIFTI_Orient NIFTIorientation; //Edit Brian Donohue
	int		filetype;	///< type of file (NEMA, Nifti, etc)
	string	filename;	///< name of volume file
	string	orientation;///< data orientation string
	int		s_units;	///< spatial units
	int		t_units;	///< time units
	int		dtype;		///< data type
	float	dx;			///< spacing along x axis
	float	dy;			///< spacing along y axis
	float	dz;			///< spacing along z axis
	float	dt;			///< time increment between volumes
	long	totvox;		///< total number of voxels in a volume
	float	scale;		///< scale factor for voxel values
	float	xoffset;	///< Talairach x offset
	float	yoffset;	///< Talairach y offset
	float	zoffset;	///< Talairach z offset

	// file info
	float	TE1;		///< echo time
	float	TE2;		///< echo time
	float	TR;			///< relaxation time
	string	scan_date;	///< date of scan
	string	pname;		///< patient name
	string	pnum;		///< patient number

	// constructors
	RicVolumeSet();
	RicVolumeSet(int nx, int ny, int nz, int nvol);
	RicVolumeSet(string filename);
	~RicVolumeSet();

	// methods
	int Init(int nx, int ny, int nz, int nvol);
	int Write_NEMA(string filename);
	int Write_NEMA(string filename, string orient_str);
	int Write_NIFTI(string filename);
	int Write_GIS(string filename);
	int	Write(string filename);
	
	/// get number of volumes in set
	int get_numvols(){return nvol;};


				
	private:
	int Read_NEMA(string filename);
	int Read_NIFTI(string filename);
	int Read_GIS(string filename);
};

// other function defines

/*!
 * function to convert orientation between nifti and nema
 */
string convertNiftiSFormToNEMA(mat44 R);

#endif //_RICVOLUMESET_H
