#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#include "allvars.h"
#include "proto.h"


/* returns the acceleration at coordinate pos[] */
double get_timestep(double *pos, double *vel, double *acc, int icell)
{
  // double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
  double v = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
  double aa = sqrt(acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]);

  double torbit = All.V200 / aa;
  double tcross = DG_CellSize[icell] / v;

  return dmin(All.TimeStepFactorOrbit * torbit, All.TimeStepFactorCellCross * tcross);
}



/* calculate the density response of a single particle starting from pos[]/vel[], averaged over time 'timespan'. If timespan=0, the routine
 * determines an appropriate time itself.
 */
double produce_orbit_response_field( double *pos, double *vel, int id, double *mfield, double mass,
												 double timespan, int *orbitstaken ) {
	
	int i, norbit, icell, flag = 0, iR, iz;
	double x[3], v[3], a[3], dt, tall, radsign_previous = 0, radsign, fR, fz;

	for(i = 0; i < 3; i++) {
		x[i] = pos[i];
		v[i] = vel[i];
	}

	for(i = 0; i < DG_Ngrid; i++)
		mfield[i] = 0;

	norbit = 0;
	tall = 0;


	forcegrid_get_acceleration(x, a);

	densitygrid_get_cell(x, &iR, &iz, &fR, &fz);
	icell = iz * DG_Nbin + iR;

	int Norbits = 100000000;

	double E0 = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) + forcegrid_get_potential(x);
	int steps = 0;

	do {
		
		dt = get_timestep(x, v, a, icell);

		if (0 < timespan)
		if (timespan < dt + tall) {
			dt = timespan - tall;
			flag = 1;
		}

		mfield[iz * DG_Nbin + iR] += 0.5 * dt * (1 - fR) * (1 - fz);
		mfield[iz * DG_Nbin + (iR + 1)] += 0.5 * dt * (fR) * (1 - fz);
		mfield[(iz + 1) * DG_Nbin + iR] += 0.5 * dt * (1 - fR) * (fz);
		mfield[(iz + 1) * DG_Nbin + (iR + 1)] += 0.5 * dt * (fR) * (fz);

		/*
		
		insertion place
		
		*/
		
		for(i = 0; i < 3; i++)
			v[i] += 0.5 * dt * a[i];

		for(i = 0; i < 3; i++)
			x[i] += dt * v[i];

		forcegrid_get_acceleration(x, a);

		for(i = 0; i < 3; i++)
			v[i] += 0.5 * dt * a[i];

		densitygrid_get_cell(x, &iR, &iz, &fR, &fz);
		icell = iz * DG_Nbin + iR;

		mfield[iz * DG_Nbin + iR] += 0.5 * dt * (1 - fR) * (1 - fz);
		mfield[iz * DG_Nbin + (iR + 1)] += 0.5 * dt * (fR) * (1 - fz);
		mfield[(iz + 1) * DG_Nbin + iR] += 0.5 * dt * (1 - fR) * (fz);
		mfield[(iz + 1) * DG_Nbin + (iR + 1)] += 0.5 * dt * (fR) * (fz);

		
		/*
		
		insertion place
		
		*/
		
		tall += dt;

		radsign = v[0] * x[0] + v[1] * x[1] + v[2] * x[2];

		if(radsign > 0 && radsign_previous < 0)
			norbit++;

		radsign_previous = radsign;

		steps++;
		if(steps > 100000000) {
			printf("too many steps...  pos=(%g|%g|%g)  vel=(%g|%g|%g)  dt=%g\n",
			pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], dt);
			double E1 = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) + forcegrid_get_potential(x);
			printf("steps=%d:  rel error = %g\n", steps, fabs(E1 - E0) / fabs(E0));
			exit(1);
		}
		
	} while ((timespan == 0 && norbit < Norbits) || (timespan != 0 && flag == 0));

	double E1 = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) + forcegrid_get_potential(x);

	double rel_egy_error = fabs((E1 - E0) / E0);

	if(rel_egy_error > 0.5) {
      mpi_printf("relative energy error= %g  orbits=%d   steps=%d  pos(=%g|%g|%g) vel=(%g|%g|%g)\n", rel_egy_error, norbit, steps,
		 pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]);
      /*
		terminate("error seems large, we better stop:  pos=(%g|%g|%g)  vel=(%g|%g|%g) id=%d  v=%g  vesc=%g",
		pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], id, 
		sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]), 
		forcegrid_get_escape_speed(pos));
       */
	}

	double fac = mass / tall;

	for(i = 0; i < DG_Ngrid; i++)
		mfield[i] *= fac;

	*orbitstaken = norbit;

	return tall;
}


#ifdef VER_1_1
double produce_orbit_response_field_mod( double *pos, double *vel, int id, 
													  double *mfield, double *egyfield_r, double *egyfield_t, double *egyfield_q, double *egyfield_p,
													  double mass, double timespan, int *orbitstaken, int type ) {
	
	int typeOfVelocityStructure = 0;

  if(type == 1)               /* a halo particle */
    typeOfVelocityStructure = All.TypeOfHaloVelocityStructure;
  else if(type == 2)          /* disk */
    typeOfVelocityStructure = All.TypeOfDiskVelocityStructure;
  else if(type == 3)          /* bulge */
    typeOfVelocityStructure = All.TypeOfBulgeVelocityStructure;
  else
    terminate("unknown type");
	
	int i, norbit, icell, flag = 0, iR, iz;
	double x[3], v[3], a[3], dt, tall, radsign_previous = 0, radsign, fR, fz;
	double r2, v_dot_x, vr2;
	double Z[] = {0,0,-1};
	double T[3], Q[3];
	double q, q2, vq, vq2, v_dot_Q;
	double t2, vt2, v_dot_T;
	double vstr, vp2;
	int irz[2][2];
	double m[2][2];
		
	for(i = 0; i < 3; i++) {
		x[i] = pos[i];
		v[i] = vel[i];
	}

	for(i = 0; i < DG_Ngrid; i++)
		mfield[i] = 0;

	for(i = 0; i < EG_Ngrid; i++) {
		egyfield_r[i] = 0;
		egyfield_t[i] = 0;
		egyfield_q[i] = 0;
		egyfield_p[i] = 0;		
	}
	
	
	norbit = 0;
	tall = 0;


	forcegrid_get_acceleration(x, a);

	densitygrid_get_cell(x, &iR, &iz, &fR, &fz);
	icell = iz * DG_Nbin + iR;

	int Norbits = 100000000;

	double E0 = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) + forcegrid_get_potential(x);
	int steps = 0;

	do {
		
		dt = get_timestep(x, v, a, icell);

		if (0 < timespan)
		if (timespan < dt + tall) {
			dt = timespan - tall;
			flag = 1;
		}

		
		if(typeOfVelocityStructure == 2) {

			// radial
			r2 = x[0]*x[0] + x[1]*x[1];
			v_dot_x = v[0]*x[0] + v[1]*x[1];
			vr2 = v_dot_x * v_dot_x / r2;

			// phi
			Q[0] = -x[1];
			Q[1] =  x[0];
			q2 = Q[0]*Q[0] + Q[1]*Q[1];
			v_dot_Q = v[0]*Q[0] + v[1]*Q[1];
			q = sqrt(q2);
			vq = v_dot_Q / q;
			vq2 = vq*vq;
			
			// phi - vstr
			vstr = get_vstream(x, type);
			vp2 = (vq-vstr)*(vq-vstr);	
					
			// theta
			vt2 = v[2]*v[2];
				
		} else {
		
			// radial
			r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
			v_dot_x = v[0] * x[0] + v[1] * x[1] + v[2] * x[2];
			vr2 = v_dot_x * v_dot_x / r2;

			// phi
			Q[0] = x[1]*Z[2] - x[2]*Z[1];
			Q[1] = x[2]*Z[0] - x[0]*Z[2];
			Q[2] = x[0]*Z[1] - x[1]*Z[0];
			q2 = Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2];
			v_dot_Q = v[0]*Q[0] + v[1]*Q[1] + v[2]*Q[2];
			q = sqrt(q2);
			vq = v_dot_Q / q;
			vq2 = vq*vq;
			
			// phi - vstr
			vstr = get_vstream(x, type);
			vp2 = (vq-vstr)*(vq-vstr);
			
			// theta
			T[0] = x[1]*Q[2] - x[2]*Q[1];
			T[1] = x[2]*Q[0] - x[0]*Q[2];
			T[2] = x[0]*Q[1] - x[1]*Q[0];
			t2 = T[0]*T[0] + T[1]*T[1] + T[2]*T[2];
			v_dot_T = v[0]*T[0] + v[1]*T[1] + v[2]*T[2];
			vt2 = v_dot_T * v_dot_T / t2;
		}
		
		
		// mass
		m[0][0] = 0.5 * dt * (1 - fR) * (1 - fz);
		m[1][0] = 0.5 * dt * (fR) * (1 - fz);
		m[0][1] = 0.5 * dt * (1 - fR) * (fz);
		m[1][1] = 0.5 * dt * (fR) * (fz);

		
		irz[0][0] = iz * DG_Nbin + iR;
		irz[1][0] = iz * DG_Nbin + (iR + 1);
		irz[0][1] = (iz + 1) * DG_Nbin + iR;
		irz[1][1] = (iz + 1) * DG_Nbin + (iR + 1);
		
		
		// m
		mfield[irz[0][0]] += m[0][0];
		mfield[irz[1][0]] += m[1][0];
		mfield[irz[0][1]] += m[0][1];
		mfield[irz[1][1]] += m[1][1];
		
		// mvr2
		egyfield_r[irz[0][0]] += m[0][0] * vr2;
		egyfield_r[irz[1][0]] += m[1][0] * vr2;
		egyfield_r[irz[0][1]] += m[0][1] * vr2;
		egyfield_r[irz[1][1]] += m[1][1] * vr2;

		// mvt2
		egyfield_t[irz[0][0]] += m[0][0] * vt2;
		egyfield_t[irz[1][0]] += m[1][0] * vt2;
		egyfield_t[irz[0][1]] += m[0][1] * vt2;
		egyfield_t[irz[1][1]] += m[1][1] * vt2;

		
		// mvq2 (2nd-moment)
		egyfield_q[irz[0][0]] += m[0][0] * vq2;
		egyfield_q[irz[1][0]] += m[1][0] * vq2;
		egyfield_q[irz[0][1]] += m[0][1] * vq2;
		egyfield_q[irz[1][1]] += m[1][1] * vq2;		
		
		
		// mvp2 (dispersion)
		egyfield_p[irz[0][0]] += m[0][0] * vp2;
		egyfield_p[irz[1][0]] += m[1][0] * vp2;
		egyfield_p[irz[0][1]] += m[0][1] * vp2;
		egyfield_p[irz[1][1]] += m[1][1] * vp2;		
		
		
		
		for(i = 0; i < 3; i++)
			v[i] += 0.5 * dt * a[i];

		for(i = 0; i < 3; i++)
			x[i] += dt * v[i];

		forcegrid_get_acceleration(x, a);

		for(i = 0; i < 3; i++)
			v[i] += 0.5 * dt * a[i];

		densitygrid_get_cell(x, &iR, &iz, &fR, &fz);
		icell = iz * DG_Nbin + iR;

		
		if(typeOfVelocityStructure == 2) {

			// radial
			r2 = x[0]*x[0] + x[1]*x[1];
			v_dot_x = v[0]*x[0] + v[1]*x[1];
			vr2 = v_dot_x * v_dot_x / r2;

			// phi
			Q[0] = -x[1];
			Q[1] =  x[0];
			q2 = Q[0]*Q[0] + Q[1]*Q[1];
			v_dot_Q = v[0]*Q[0] + v[1]*Q[1];
			q = sqrt(q2);
			vq = v_dot_Q / q;
			vq2 = vq*vq;
			
			// phi - vstr
			vstr = get_vstream(x, type);
			vp2 = (vq-vstr)*(vq-vstr);	
					
			// theta
			vt2 = v[2]*v[2];
				
		} else {
		
			r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
			v_dot_x = v[0] * x[0] + v[1] * x[1] + v[2] * x[2];
			vr2 = v_dot_x * v_dot_x / r2;
			
			
			// phi
			Q[0] = x[1]*Z[2] - x[2]*Z[1];
			Q[1] = x[2]*Z[0] - x[0]*Z[2];
			Q[2] = x[0]*Z[1] - x[1]*Z[0];
			q2 = Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2];
			v_dot_Q = v[0]*Q[0] + v[1]*Q[1] + v[2]*Q[2];
			q = sqrt(q2);
			vq = v_dot_Q / q;
			vq2 = vq*vq;
			
			// phi - vstr
			vstr = get_vstream(x, type);
			vp2 = (vq-vstr)*(vq-vstr);
			
			
			// theta
			T[0] = x[1]*Q[2] - x[2]*Q[1];
			T[1] = x[2]*Q[0] - x[0]*Q[2];
			T[2] = x[0]*Q[1] - x[1]*Q[0];
			t2 = T[0]*T[0] + T[1]*T[1] + T[2]*T[2];
			v_dot_T = v[0]*T[0] + v[1]*T[1] + v[2]*T[2];
			vt2 = v_dot_T * v_dot_T / t2;
		
		}
		
		// mass
		m[0][0] = 0.5 * dt * (1 - fR) * (1 - fz);
		m[1][0] = 0.5 * dt * (fR) * (1 - fz);
		m[0][1] = 0.5 * dt * (1 - fR) * (fz);
		m[1][1] = 0.5 * dt * (fR) * (fz);

		
		irz[0][0] = iz * DG_Nbin + iR;
		irz[1][0] = iz * DG_Nbin + (iR + 1);
		irz[0][1] = (iz + 1) * DG_Nbin + iR;
		irz[1][1] = (iz + 1) * DG_Nbin + (iR + 1);
		
		
		// m
		mfield[irz[0][0]] += m[0][0];
		mfield[irz[1][0]] += m[1][0];
		mfield[irz[0][1]] += m[0][1];
		mfield[irz[1][1]] += m[1][1];
		
		// mvr2
		egyfield_r[irz[0][0]] += m[0][0] * vr2;
		egyfield_r[irz[1][0]] += m[1][0] * vr2;
		egyfield_r[irz[0][1]] += m[0][1] * vr2;
		egyfield_r[irz[1][1]] += m[1][1] * vr2;
		
		// mvt2
		egyfield_t[irz[0][0]] += m[0][0] * vt2;
		egyfield_t[irz[1][0]] += m[1][0] * vt2;
		egyfield_t[irz[0][1]] += m[0][1] * vt2;
		egyfield_t[irz[1][1]] += m[1][1] * vt2;
		

		// mvp2
		egyfield_q[irz[0][0]] += m[0][0] * vq2;
		egyfield_q[irz[1][0]] += m[1][0] * vq2;
		egyfield_q[irz[0][1]] += m[0][1] * vq2;
		egyfield_q[irz[1][1]] += m[1][1] * vq2;
		
		
		// mvp2 (sigmap^2)
		egyfield_p[irz[0][0]] += m[0][0] * vp2;
		egyfield_p[irz[1][0]] += m[1][0] * vp2;
		egyfield_p[irz[0][1]] += m[0][1] * vp2;
		egyfield_p[irz[1][1]] += m[1][1] * vp2;		
		
		
		
		tall += dt;

		radsign = v[0] * x[0] + v[1] * x[1] + v[2] * x[2];

		if(radsign > 0 && radsign_previous < 0)
			norbit++;

		radsign_previous = radsign;

		steps++;
		if(steps > 100000000) {
			printf("too many steps...  pos=(%g|%g|%g)  vel=(%g|%g|%g)  dt=%g\n",
			pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], dt);
			double E1 = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) + forcegrid_get_potential(x);
			printf("steps=%d:  rel error = %g\n", steps, fabs(E1 - E0) / fabs(E0));
			exit(1);
		}
		
	} while ((timespan == 0 && norbit < Norbits) || (timespan != 0 && flag == 0));

	double E1 = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) + forcegrid_get_potential(x);

	double rel_egy_error = fabs((E1 - E0) / E0);

	if(rel_egy_error > 0.5) {
      mpi_printf("relative energy error= %g  orbits=%d   steps=%d  pos(=%g|%g|%g) vel=(%g|%g|%g)\n", rel_egy_error, norbit, steps,
		 pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]);
      /*
		terminate("error seems large, we better stop:  pos=(%g|%g|%g)  vel=(%g|%g|%g) id=%d  v=%g  vesc=%g",
		pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], id, 
		sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]), 
		forcegrid_get_escape_speed(pos));
       */
	}

	double fac = mass / tall;

	for(i = 0; i < DG_Ngrid; i++)
		mfield[i] *= fac;
	
	for(i = 0; i < EG_Ngrid; i++) {

		egyfield_r[i] *= fac;
		egyfield_t[i] *= fac;
		egyfield_q[i] *= fac;
		egyfield_p[i] *= fac;		
	}

	*orbitstaken = norbit;

	return tall;
}

#endif
