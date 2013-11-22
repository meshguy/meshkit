#ifndef TFIBLEND_H
#define TFIBLEND_H

void blend_from_corners ( double *x, int m );
void blend_from_corners ( double *x, int nx, int ny );
void blend_from_corners ( double *x, int nx, int ny, int nz );

void blend_from_edges   ( double *x, int nx, int ny );
void blend_from_edges   ( double *x, int nx, int ny, int nz );

void blend_from_faces   ( double *x, int nx, int ny, int nz );

#endif
