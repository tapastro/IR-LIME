/*
 *  photon.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 15/11/06.
 *  Copyright 2006-2017, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *  All rights reserved.
 *
 */

#include "lime.h"

int
sortangles(double *inidir, int id, struct grid *g, const gsl_rng *ran) {
  int i,n[2];
  double angle,exitdir[2];

  exitdir[0]=1e30;
  exitdir[1]=1e31;
  n[0]=-1;
  n[1]=-1;
  for(i=0;i<g[id].numNeigh;i++){
    angle=( inidir[0]*g[id].dir[i].xn[0]
           +inidir[1]*g[id].dir[i].xn[1]
           +inidir[2]*g[id].dir[i].xn[2]);
    if(angle<exitdir[0]){
      exitdir[1]=exitdir[0];
      n[1]=n[0];
      exitdir[0]=angle;
      n[0]=i;
    } else if(angle<exitdir[1]) {
      exitdir[1]=angle;
      n[1]=i;
    }
  }
  if(gsl_rng_uniform(ran)<1./((1-exitdir[0])/(1-exitdir[1])+1) ) {
    if(n[0]==-1){
      if(!silent) bail_out("Photon propagation error");
      exit(1);
    }
    return n[0];
  } else {
    if(n[1]==-1){
      if(!silent) bail_out("Photon propagation error");
      exit(1);
    }
    return n[1];
  }
}



void
velocityspline(struct grid *g, int id, int k, double binv, double deltav, double *vfac){
  int nspline,ispline,naver,iaver;
  double v1,v2,s1,s2,sd,v,vfacsub,d;

  v1=deltav-veloproject(g[id].dir[k].xn,g[id].vel);
  v2=deltav-veloproject(g[id].dir[k].xn,g[id].neigh[k]->vel);

  nspline=(fabs(v1-v2)*binv < 1) ? 1 : (int)(fabs(v1-v2)*binv);
  *vfac=0.;
  s2=0;
  v2=v1;

  for(ispline=0;ispline<nspline;ispline++){
    s1=s2;
    s2=((double)(ispline+1))/(double)nspline;
    v1=v2;
    d=s2*g[id].ds[k];
    v2=deltav-(g[id].a4[k]*pow(d,4)+g[id].a3[k]*pow(d,3)+g[id].a2[k]*pow(d,2)+g[id].a1[k]*d+g[id].a0[k]);
    naver=(1 > fabs(v1-v2)*binv) ? 1 : (int)(fabs(v1-v2)*binv);
    for(iaver=0;iaver<naver;iaver++){
      sd=s1+(s2-s1)*((double)iaver-0.5)/(double)naver;
      d=sd*g[id].ds[k];
      v=deltav-(g[id].a4[k]*pow(d,4)+g[id].a3[k]*pow(d,3)+g[id].a2[k]*pow(d,2)+g[id].a1[k]*d+g[id].a0[k]);
      vfacsub=gaussline(v,binv);
      *vfac+=vfacsub/(double)naver;
    }
  }
  *vfac= *vfac/(double)nspline;
  return;
}


void
velocityspline_lin(struct grid *g, int id, int k, double binv, double deltav, double *vfac){
  int nspline,ispline,naver,iaver;
  double v1,v2,s1,s2,sd,v,vfacsub,d;

  v1=deltav-veloproject(g[id].dir[k].xn,g[id].vel);
  v2=deltav-veloproject(g[id].dir[k].xn,g[id].neigh[k]->vel);

  nspline=(fabs(v1-v2)*binv < 1) ? 1 : (int)(fabs(v1-v2)*binv);
  *vfac=0.;
  s2=0;
  v2=v1;

  for(ispline=0;ispline<nspline;ispline++){
    s1=s2;
    s2=((double)(ispline+1))/(double)nspline;
    v1=v2;
    d=s2*g[id].ds[k];
    v2=deltav-(g[id].a1[k]*d+g[id].a0[k]);
    naver=(1 > fabs(v1-v2)*binv) ? 1 : (int)(fabs(v1-v2)*binv);
    for(iaver=0;iaver<naver;iaver++){
      sd=s1+(s2-s1)*((double)iaver-0.5)/(double)naver;
      d=sd*g[id].ds[k];
      v=deltav-(g[id].a1[k]*d+g[id].a0[k]);
      vfacsub=gaussline(v,binv);
      *vfac+=vfacsub/(double)naver;
    }
  }
  *vfac= *vfac/(double)nspline;
  return;
}

double veloproject(double dx[3], double *vel){
  return dx[0]*vel[0]+dx[1]*vel[1]+dx[2]*vel[2];
}


double gaussline(double v, double sigma){
  int maxgau=101,maxsig=4,ival;
  double fac,val;

  fac=(maxgau-1)/maxsig;
  ival=(int)(fac*(fabs(v)*sigma))+1;
  if((ival-1)>=maxgau) return 0.;
  val=(ival*ival)/(fac*fac);
  return exp(-val);
}

void
photon(int id, struct grid *g, molData *m, int iter, const gsl_rng *ran,inputPars *par,blend *matrix,double *phot, double *vfac, double *ds){
  int here,there,firststep,dir,np_per_line,ip_at_line,l,iphot,iline,jline;
  double deltav,segment,vblend,snu,dtau,jnu,alpha,pt_theta,pt_z, lds,*lvfac;
  double *tau,x[3], inidir[3],vel[3];
  int *counta, *countb,nlinetot;


  lineCount(par->nSpecies, m, &counta, &countb, &nlinetot);
  velocity(g[id].x[0],g[id].x[1],g[id].x[2],vel);

  tau=malloc(sizeof(double)*nlinetot);
  lvfac=malloc(sizeof(double)*par->nSpecies);
  
  
  for(iphot=0;iphot<g[id].nphot;iphot++){
    firststep=1;
    for(iline=0;iline<nlinetot;iline++){
      phot[iline+iphot*m[0].nline]=0.;
      tau[iline]=0.;
    }

    /* Initial velocity, direction and frequency offset  */
    pt_theta=gsl_rng_uniform(ran)*2*PI;
    pt_z=2*gsl_rng_uniform(ran)-1;

    inidir[0]=sqrt(1-pow(pt_z,2))*cos(pt_theta);
    inidir[1]=sqrt(1-pow(pt_z,2))*sin(pt_theta);
    inidir[2]=pt_z;

    iter=(int) (gsl_rng_uniform(ran)*3.);
    np_per_line=(int) g[id].nphot/g[id].numNeigh;
    ip_at_line=(int) iphot/g[id].numNeigh;
    segment=1/(2.*np_per_line)*(2*ip_at_line-np_per_line+iter);

    dir=sortangles(inidir,id,g,ran);
    here=g[id].id;
    there=g[here].neigh[dir]->id;
    deltav=segment*4.3*g[id].dopb+veloproject(g[id].dir[dir].xn,vel);


    /* Photon propagation loop */
    do{
      if(firststep){
        firststep=0;
        lds=g[here].ds[dir]/2.;
        
   
        for(l=0;l<par->nSpecies;l++){
          if(!par->pregrid && !par->discreteVelocities) velocityspline(g,here,dir,g[id].mol[l].binv,deltav,&lvfac[l]);
          else velocityspline_lin(g,here,dir,g[id].mol[l].binv,deltav,&lvfac[l]);
          vfac[iphot]=lvfac[l];
          ds[iphot]=lds;
        }
        for(l=0;l<3;l++) x[l]=g[here].x[l]+(g[here].dir[dir].xn[l] * g[id].ds[dir]/2.);
      } else {
        lds=g[here].ds[dir];
        for(l=0;l<3;l++) x[l]=g[here].x[l];
      }

      for(l=0;l<par->nSpecies;l++){
        if(!par->pregrid && !par->discreteVelocities) velocityspline(g,here,dir,g[id].mol[l].binv,deltav,&lvfac[l]);
        else velocityspline_lin(g,here,dir,g[id].mol[l].binv,deltav,&lvfac[l]);
      }

      for(iline=0;iline<nlinetot;iline++){
        jnu=0.;
        alpha=0.;
        snu=0.;
        dtau=0.;

        sourceFunc_line(&jnu,&alpha,m,lvfac[counta[iline]],g,here,counta[iline],countb[iline]);
        sourceFunc_cont(&jnu,&alpha,g,here,counta[iline],countb[iline]);
        if(fabs(alpha)>0.){
          snu=(jnu/alpha)*m[0].norminv;
          dtau=alpha*lds;
          if(dtau < -30) dtau = -30;
        }

        phot[iline+iphot*m[0].nline]+=exp(-tau[iline])*(1.-exp(-dtau))*snu;
        tau[iline]+=dtau;
        if(tau[iline] < -30.){
          if(!silent) warning("Maser warning: optical depth has dropped below -30");
          tau[iline]= -30.;
        }

        /* Line blending part */
        if(par->blend){
          jnu=0.;
          alpha=0.;
          for(jline=0;jline<sizeof(matrix)/sizeof(blend);jline++){
            if(matrix[jline].line1 == jline || matrix[jline].line2 == jline){
              if(!par->pregrid && !par->discreteVelocities) velocityspline(g,here,dir,g[id].mol[counta[jline]].binv,deltav-matrix[jline].deltav,&vblend);
              else velocityspline_lin(g,here,dir,g[id].mol[counta[jline]].binv,deltav-matrix[jline].deltav,&vblend);
              sourceFunc_line(&jnu,&alpha,m,vblend,g,here,counta[jline],countb[jline]);
              if(fabs(alpha)>0.){
                snu=(jnu/alpha)*m[0].norminv;
                dtau=alpha*lds;
                if(dtau < -30) dtau = -30;
              }
              phot[jline+iphot*m[0].nline]+=exp(-tau[jline])*(1.-exp(-dtau))*snu;
              tau[jline]+=dtau;
              if(tau[jline] < -30.){
                if(!silent) warning("Optical depth has dropped below -30");
                tau[jline]= -30.;
              }
            }
          }
        }
        /* End of line blending part */
      }

      dir=sortangles(inidir,there,g,ran);
      here=there;
      there=g[here].neigh[dir]->id;
    } while(!g[there].sink);

    /* Add cmb contribution */
    if(m[0].cmb[0]>0.){
      for(iline=0;iline<nlinetot;iline++){
        phot[iline+iphot*m[0].nline]+=exp(-tau[iline])*m[counta[iline]].cmb[countb[iline]];
      }
    }
  }

  free(lvfac);
  free(tau);
  free(counta);
  free(countb);
}

void
getjbar(int posn, molData *m, struct grid *g, inputPars *par, double *jbar, double *phot, double *vfac, double *ds){
 int iline,iphot;
 double tau, snu, vsum=0., jnu, alpha;
 int *counta, *countb,nlinetot;

 lineCount(par->nSpecies, m, &counta, &countb, &nlinetot);

 for(iline=0;iline<m[0].nline;iline++) jbar[iline]=0.;
 for(iphot=0;iphot<g[posn].nphot;iphot++){
   if(vfac[iphot]>0){
     for(iline=0;iline<m[0].nline;iline++){
       jnu=0.;
       alpha=0.;
       snu=0.;
       tau=0.;

       sourceFunc_line(&jnu,&alpha,m,vfac[iphot],g,posn,counta[iline],countb[iline]);
       sourceFunc_cont(&jnu,&alpha,g,posn,counta[iline],countb[iline]);
       if(fabs(alpha)>0.){
         snu=(jnu/alpha)*m[0].norminv;
         tau=alpha*ds[iphot];
       }
       jbar[iline]+=vfac[iphot]*(exp(-tau)*phot[iline+iphot*m[0].nline]+(1.-exp(-tau))*snu);
     }
     vsum+=vfac[iphot];
   }
 }
 for(iline=0;iline<m[0].nline;iline++) jbar[iline]=m[0].norm*jbar[iline]/vsum;
 free(counta);
 free(countb);
}
