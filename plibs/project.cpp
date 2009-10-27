  int ITAPSurface :: project (Point<3> & p, double tolerance) const
  {
     double x[3], xold[3], n[3];
     double dx, dy, dz;
       
     iGeom_getEntClosestPt( geometry, faceHandle, p[0], p[1], p[2],
                            &x[0], &x[1], &x[2], &err);

     dx = p[0] - x[0];
     dy = p[1] - x[1];
     dz = p[2] - x[2];

     double sqr_tol = tolerance*tolerance;

     if( dx*dx + dy*dy + dz*dz  < sqr_tol ) return 0;

     xold[0] = x[0];
     xold[1] = x[1];
     xold[2] = x[2];

     double u,v;
     SearchUV seachUV;
     searchUV.getUV( geometry, faceHandle, x[0], x[1], x[2], u, v);

     SimpleArray<double> du, dv;
     iGeom_getEnt1stDrvt( geometry, faceHandle, u, v, ARRAY_INOUT(du), ARRAY_INOUT(dv), &err);

     int count = 0;
     double det, lambda, mu;

    while(1)
    {
      count++;

      n = du^dv;

      det = Det3 (n[0], du[0], dv[0],
		  n[1], du[1], dv[1],
		  n[2], du[2], dv[2] );

      if (det < 1e-15) return false;

      lambda = Det3 (n[0], p[0] - x[0], dv[0],
		     n[1], p[1] - x[1], dv[1],
		     n[2], p[2] - x[2], dv[2] )/det;

      mu     = Det3 (n[0], du[0], p[0] - x[0],
		     n[1], du[1], p[1] - x[1],
		     n[2], du[2], p[2] - x[2] )/det;

      u += lambda;
      v += mu;

      iGeom_getEntUVtoXYZ( geometry, faceHandle, u, v, &x[0], &x[1], &x[2], &err);

      dx = xold[0] - x[0];
      dy = xold[1] - x[1];
      dz = xold[2] - x[2];

      if( dx*dx + dy*dy + dz*dz < sqr_tol || count >= 100) break;

      xold[0] = x[0];
      xold[1] = x[1];
      xold[2] = x[2];

      du.clear();
      dv.clear();
      iGeom_getEnt1stDrvt( geometry, faceHandle, u, v, ARRAY_INOUT(du), ARRAY_INOUT(dv), &err);
    } 

    if( count >=100) return 2;

    p[0] = x[0];
    p[1] = x[1];
    p[2] = x[2];

    return 0;
  }

////////////////////////////////////////////////////////////////////////////////

