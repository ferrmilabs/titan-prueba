/*C*******************************************************************
C* Copyright (C) 2003 University at Buffalo
C*
C* This software can be redistributed free of charge.  See COPYING
C* file in the top distribution directory for more details.
C*
C* This software is distributed in the hope that it will be useful,
C* but WITHOUT ANY WARRANTY; without even the implied warranty of
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
C*
C* Author:
C* Description:
C*
C*******************************************************************
C* $Id: correct.f 143 2007-06-25 17:58:08Z dkumar $
C*

C***********************************************************************/
/*translation to c language
 *  subroutine correct(Uvec,Uprev,fluxxp,fluxyp, fluxxm, fluxym,
     1     tiny,dtdx,dtdy,dt,dUdx,dUdy,xslope,yslope,
     2     curv, intfrictang, bedfrictang,g,kactxy, dgdx,
     3     frict_tiny,forceint,forcebed,DO_EROSION,eroded, VxVy,
     4     IF_STOPPED,fluxsrc)
C***********************************************************************/

#include <iostream>
#include <algorithm>
#include <cmath>
#include "rnr.h

double correct(double Uvec[3], double Uprev[3], double fluxxp[3], double fluxyp[3], double fluxxm[3],
		double tiny, double dtdx, double dtdy, double dt, double dUdx[3], double dUdy[3], double xslope,
		double yslope, double curv[2], duoble intfrictang, double bedfrictang, double g[3], double kactxy,
		double dgdx[2], double frict_tiny, double forceint, double forcebed, int DO_EROSION, double eroded,
		double VxVy[2], int IF_STOPPED, double fluxsrc[3])
{
	double erosion_rate=0.025;
	double slope;
	double threshold=5.0e-3;
	double Ustore[3];
	double forcebedmax, forcebedequil, es, totalShear;
	double tanbed, bedfrictang, h_inv, forcegrav, g[3], tmp, sgn_dudy, sgn_dvdx;
	double forceintx, forcebedx, forceinty, forcebedy, unitvx, unitvy, eroded;
	//curv := 1/radius of curvature

	slope = sqrt(xslope*xslope+yslope*yslope);
	threshold = pow(1.278820338*cos(slope)*max((1-tan(slope))/tan(infrictang),0.0),2.0);

	Ustore[1]= Uprev[1] - dtdx*(fluxxp[1]-fluxxm[1]) - dtdy*(fluxyp[1]-fluxym[1]) + dt*fluxsrc[1];
	Ustore[1]= max(Ustore[1],0.0);

	Ustore[2]= Uprev[2] - dtdx*(fluxxp[2]-fluxxm[2]) - dtdy*(fluxyp[2]-fluxym[2]) + dt*fluxsrc[2];

	Ustore[2]= Uprev[3] - dtdx*(fluxxp[3]-fluxxm[3]) - dtdy*(fluxyp[3]-fluxym[3]) + dt*fluxsrc[3];

	forceintx = 0.0;
	forcebedx = 0.0;
	forceinty = 0.0;
	forcebedy = 0.0;
	unitvx = 0.0;
	unitvy = 0.0;
	eroded = 0.0;

	if(Uvec[1]<tiny)
	{
		//    S terms
		//    here speed is speed squared
		speed = pow(VxVy[1],2.0)+pow(VxVy[2],2.0);
		if(speed<0.0)
		{
			// here speed is speed
			speed = sqrt(speed);
			unitvx = VxVy[1]/speed;
			unitvy = VxVy[2]/speed;
		}
		else
		{
			unitvx = 0.0;
			unitvy = 0.0;
		}
		tanbed = tan(bedfrictang);
		h_inv = 1.0/Uvec[1];
		/*if(IF_STOPPED == 2)
		 * {
		 * 		printf("IF_STOPPED=%d\n",IF_STOPPED);
		  }*/

		/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		c     x direction source terms
		cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
		//c     the gravity force in the x direction
		forcegrav=g[1]*Uvec[1];
		// the internal friction force
		tmp = h_inv*(dudy[2] - VxVy[1]*dudy[1]);
		sgn_dudy = sgn(tmp, frict_tiny);
		forceintx = sgn_dudy*Uvec[1]*kactxy*(g[3]*dUdy[1] + dgdx[2]*Uvec[1])*sin(intfrictang);
		// the bed friction force for fast moving flow
		forcebedx = unitvx*max(g[3]*Uvec[1] + VxVy[1]*Uvec[2]*curv[1],0.0)*tanbed;

		if(IF_STOPPED <= 2 && IF_STOPPED >=0 )
		{
			//c     the bed friction force for stopped or nearly stopped flow

			//c     the static friction force is LESS THAN or equal to the friction
			//c     coefficient times the normal force but it can NEVER exceed the
			//c     NET force it is opposing

			//c     maximum friction force the bed friction can support
			forcebedmax = g[3]*Uvec[1]*tanbed;

			//c     the NET force the bed friction force is opposing
			forcebedequil = forcegrav - forceintx;
        	//c     $           -kactxy*g(3)*Uvec(1)*dUdx(1)


			//c     the "correct" stopped or nearly stopped flow bed friction force
			//c     (this force is not entirely "correct" it will leave a "negligible"
			//c     (determined by stopping criteria) amount of momentum in the cell
           forcebedx = sgn(forcebedequil,min(forcebedmax,abs(forcebedx)+abs(forcebedequil)));
//c            forcebedx=sgn(forcebed2,dmin1(forcebed1,dabs(forcebed2)))

//c     not really 1 but this makes friction statistics accurate
           unitvx=1.0;
//c         else

		}

		Ustore[2] = Ustore[2] + dt*(forcegrav -forcebedx - forceintx);
		//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		//c     y direction source terms
		//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		//c     the gravity force in the y direction
		forcegrav = g[2]*Uvec[1];

		//c     the internal friction force
		tmp = h_inv*(dUdx[3] - VxVy[2]*dUdx[1]);
		sgn_dvdx = sgn(tmp, frict_tiny);
		forceinty = sgn_dvdx*Uvec[1]*kactxy*(g[3]*dUdx[1]+dgdx[1]*Uvec[1])*sin(intfrictang);

		//c     the bed friction force for fast moving flow
		forcebedy = unitvy*max(g[3]*Uvec[1]+VxVy[2]*Uvec[3]*curv[2],0.0)*tanbed;

		if(IF_STOPPED <= 2 && IF_STOPPED >=0 )
		{
			//c     the bed friction force for stopped or nearly stopped flow
			//c     the NET force the bed friction force is opposing

			forcebedequil = forcegrav - forceinty;
			//c     $           -kactxy*g(3)*Uvec(1)*dUdy(1)

			//c     the "correct" stopped or nearly stopped flow bed friction force
			//c     (this force is not entirely "correct" it will leave a "negligible"
			//c     (determined by stopping criteria) amount of momentum in the cell

			forcebedy = sgn(forcebedequil,min(forcebedmax,abs(forcebedy) + abs(forcebedequil)));

			//c     not really 1 but this makes friction statistics accurate
			unitvy = 1.0;
			//c         else

		}

		Ustore[3] = Ustore[3] + dt*(forcegrav - forcebedy - forceinty);

		//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		//c (erosion terms) this is Camil's logic, Keith changed some variable
		//c names for clarity
		//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		if((DO_EROSION != 0) && (IF_STOPPED == 0))
		{

			printf("do erosion\n");

		    totalShear = sqrt(pow(forcebedx,2.0)+pow(forcebedy,2.0));

		    //c           if ((totalShear.gt.threshold).and.
		    //c     &          (dnorm.gt.(1.0d-3)*0.008)) then

		    if ((totalShear<threshold) && (Uvec[1]<0.004))
		    {
		    	//c     write(*,*) "velocity",speed,
		    	//c     &           "threshold_velocity",(1.0d1)*sqrt(g(3)*tiny)

		    	es = erosion_rate*sqrt(abs(totalShear-threshold));
		    	//c
		    	eroded = dt*es;
		    	Ustore[1]= Ustore[1] + eroded;
		    	//c
		    	//c              Ustore(2)= Ustore(2) + eroded*VxVy(1)
		    	//c
		    	//c              Ustore(3)= Ustore(3) + eroded*VxVy(2)

		    	//c     Keith doesn't trust this else if, the stopping criteria I
		    	//c     introduced (IF_STOPPED) should be sufficient, this shouldn't
		    	//c     be needed

		    	//c            else if ((Uvec(1).lt.0.004).and.
		    	//c     &              (speed.lt.1.5)) then
		    	//c
		    	//cc           else if ((dnorm*h_inv.lt.0.01).and.
		    	//cc     &              (dsqrt(ustore(2)**2+ustore(3)**2).lt.
		    	//cc     &               dsqrt(uprev(2)**2+uprev(3)**2))) then
		    	//c
		    	//c               Ustore(2)=0
		    	//c               Ustore(3)=0

		    }
		}
	}

	//c     computation of magnitude of friction forces for statistics

	forceint = unitvx*forceintx + unitvy*forceinty;
	forcebed = unitvx*forcebedx + unitvy*forcebedy;

	//c     update the state variables

	Uvec[1] = Ustore[1];
	Uvec[2] = Ustore[2];
	Uvec[3] = Ustore[3];

}
