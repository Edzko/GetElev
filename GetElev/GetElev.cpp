// GetElev.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <windows.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "mex.h"

#define FSIZE 1000000
#define GSIZE 5000
#define NOBJECTS 500
#define MAX_FILENAME 250
#define PI 3.1415

char wbtdata[FSIZE];

struct {
	int xdim, zdim;
	double x0, z0;
	double xspacing, zspacing;
	double elevation[GSIZE];
} Elevation;

double get_height(double x, double z)
{
	if (x > Elevation.x0) x = Elevation.x0;
	if (x < -Elevation.x0) x = -Elevation.x0;
	if (z < Elevation.z0) z = Elevation.z0;
	if (z > -Elevation.z0) z = -Elevation.z0;

	double fx = (Elevation.x0 - x) / Elevation.xspacing;
	double fz = (-(Elevation.z0 - z) / Elevation.zspacing);
	int mx = (int)fx;
	int mz = (int)fz;
	double y0 = (fx - mx)*Elevation.elevation[mz + mx * Elevation.zdim] +
		(mx + 1 - fx)*Elevation.elevation[mz + (mx + 1)*Elevation.zdim];
	double y1 = (fx - mx)*Elevation.elevation[(mz + 1) + mx * Elevation.zdim] +
		(mx + 1 - fx)*Elevation.elevation[(mz + 1) + (mx + 1)*Elevation.zdim];
	return (fz - mz)*y0 + (mz + 1 - fz)*y1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	FILE *fid;
	double *outh, *outr = NULL, *outp = NULL, *inx, *iny, *inh = NULL;
	char ifn[500];
	int n,m,in,im,k;

	if (nrhs == 0)
	{
		mexPrintf("GetElev.mex64 - Get Elevation and gradients from WebBots elevation surface.\r");
		mexPrintf("GetElev(<filename.wbt>)   - Open WebBots world model with elevation surface.\r");
		mexPrintf("[H] = GetElev(X,Y)	     - Get height at location (X,Y).\r");
		mexPrintf("[H,R,P] = GetElev(X,Y,H)  - Get height, roll and pitch at location (X,Y) in direction H (rad).\r");
		return;
	}
	if (nrhs == 1)
	{
		if (mxIsChar(prhs[0]))
		{
			// load Webots project file
			mxGetString(prhs[0], ifn, 500);
			if (fopen_s(&fid, ifn, "r")) {
				mexPrintf("Cannot open Webots Project file.\r");
				return;
			}
			size_t nwbt = fread(wbtdata, 1, FSIZE, fid);
			fclose(fid);
			mexPrintf("Wbt file = %i bytes.\r", nwbt);
			

			// Collect ground plane data
			char *eg = strstr(wbtdata, "FLOOR Solid");
			while ((*eg < '+') || (*eg > '9')) eg++;
			sscanf_s(eg, "%lf 0 %lf", &Elevation.x0, &Elevation.z0);
			eg = strstr(wbtdata, "ElevationGrid");
			while (*eg != '[') eg++; eg++;
			k = 0;
			while (true) {
				sscanf_s(eg, "%lf", &Elevation.elevation[k++]);
				while (*eg > ' ') eg++;
				while (*eg <= ' ') eg++;
				if (*eg == ']')
					break;
				if (k >= GSIZE)
					break;
			}
			while ((*eg < '+') || (*eg > '9')) eg++;
			sscanf_s(eg, "%i", &Elevation.zdim);
			while (*eg != 10) eg++;
			while ((*eg < '+') || (*eg > '9')) eg++;
			sscanf_s(eg, "%lf", &Elevation.zspacing);
			while (*eg != 10) eg++;
			while ((*eg < '+') || (*eg > '9')) eg++;
			sscanf_s(eg, "%i", &Elevation.xdim);
			while (*eg != 10) eg++;
			while ((*eg < '+') || (*eg > '9')) eg++;
			sscanf_s(eg, "%lf", &Elevation.xspacing);
			while (*eg != 10) eg++;

			mexPrintf("Grid = (%i x %i) points at (%1.2f, %1.2f).\r",
					Elevation.xdim, Elevation.zdim, Elevation.x0, Elevation.z0);
			mexPrintf("X dim = %i, spacing = %1.2f\r", Elevation.xdim, Elevation.xspacing);
			mexPrintf("Z dim = %i, spacing = %1.2f\r", Elevation.zdim, Elevation.zspacing);

			return;
		}
		else {
			mexPrintf("Argument must be a string with the full pathname of the WebBots world model file.\r");
			return;
		}
	}
	if (nrhs >= 2)
	{
		n = (int)mxGetN(prhs[0]);
		m = (int)mxGetM(prhs[0]);
		if (((int)mxGetN(prhs[1]) != n) || ((int)mxGetM(prhs[1]) != m)) {
			mexPrintf("Arguments X and Y must have same dimensions.\r");
			return;
		}

		inx = mxGetPr(prhs[0]);
		iny = mxGetPr(prhs[1]);
		if (nrhs > 2)
			inh = mxGetPr(prhs[2]);

		plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
		outh = mxGetPr(plhs[0]);
		for (in = 0; in < n; in++) {
			for (im = 0; im < m; im++) {
				k = in*m + im;
				outh[k] = get_height(inx[k], iny[k]);
			}
		}

		if ((nlhs > 1) && (nrhs > 2)) {
			plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
			outp = mxGetPr(plhs[0]);
		}
		if ((nlhs > 2) && (nrhs > 2)) {
			plhs[2] = mxCreateDoubleMatrix(m, n, mxREAL);
			outr = mxGetPr(plhs[0]);
		}

		for (in = 0; in < n; in++) {
			for (im = 0; im < m; im++) {
				k = in * m + im;
				double h = get_height(inx[k], iny[k]);
				double h1 = get_height(inx[k] + cos(inh[k]), iny[k]+sin(inh[k]));
				outp[k] = h1 - h;
				if (nlhs > 2) {
					h1 = get_height(inx[k] + cos(inh[k]+PI/2.0), iny[k] + sin(inh[k]+PI/2.0));
					outr[k] = h1 - h;
				}
			}
		}
	}
}
