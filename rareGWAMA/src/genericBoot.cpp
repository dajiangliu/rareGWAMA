#include <algorithm> 
#include <iostream>
#include <ctime>
#include <vector>
#include <cmath>
#include <R.h>
#include <Rmath.h>
//#include <Ext/print.h>
extern "C"{
  void qfc1(double* lambda, double* delta, int* h, int *r, double *sigma, double *q, int *lim, double *acc, double* trace, int* ifault, double *res, int* ql);
  double quadForm(std::vector<double>& x,std::vector<std::vector<double> >& M,std::vector<double>& y);
  double rarecoverSumstatFunc(std::vector<double >& X_T_times_Y, std::vector<std::vector<double> >& X_T_times_X,double varY, int alternative,double*);
  double skatoSumstatFunc(std::vector<double >& X_T_times_Y, std::vector<std::vector<double> >& X_T_times_X,double varY, int alternative,std::vector<double>&,std::vector<std::vector<double> >&);
  double vtSumstatFunc(std::vector<double >& X_T_times_Y, std::vector<std::vector<double> >& X_T_times_X,double varY, int alternative,double* extraPar);
  double matrixByCol(double* mat,int m,int n, int nrow,int ncol);
  void rmvnorm(std::vector<std::vector<double> >& res,std::vector<std::vector<double> >& tmat);
  double pMixChisq(double,std::vector<double> lambda);
  double corXY(double* x, double*y, std::vector<int>& ixPerm, double& denominator)
  {
    double xMean=0.0,yMean=0.0,xyMean=0.0,noPerson=((double) ixPerm.size());
    for(int ii=0;ii<ixPerm.size();ii++)
      {
	xMean+=x[ii];
	yMean+=y[ii];
	xyMean+=(x[ii]*y[ixPerm[ii]]);
      }
    xMean/=((double) noPerson);
    yMean/=((double) noPerson);
    xyMean/=((double) noPerson);
    return (xyMean-xMean*yMean)/denominator;
  }
  
  double vectorCor(std::vector<double>& x, std::vector<double>& y,std::vector<int>& ixPerm)
  {
    double xMean=0.0,yMean=0.0,xyMean=0.0,xxMean=0.0,yyMean=0.0,noPerson=((double) ixPerm.size());
    for(int ii=0;ii<ixPerm.size();ii++)
      {
	xMean+=x[ii];
	xxMean+=(x[ii])*(x[ii]);
	yyMean+=(y[ii])*(y[ii]);
	yMean+=y[ii];
	xyMean+=(x[ii]*y[ixPerm[ii]]);
      }
    xMean/=noPerson;
    yMean/=noPerson;
    xxMean/=noPerson;
    yyMean/=noPerson;
    xyMean/=noPerson;    
    return (xyMean-xMean*yMean)/sqrt((xxMean-xMean*xMean)*(yyMean-yMean*yMean));
  }
  
  double vectorMax(std::vector<double>& x)
  {
    return(*max_element(x.begin(),x.end()));
  }
  
  void rmvnorm(std::vector<std::vector<double> >& res,std::vector<std::vector<double> >& tmat)
  {
    int ii=0,jj=0,kk=0;
    int n=res.size(),d=(res[0]).size();
    std::vector<std::vector<double> > res0(n,std::vector<double>(d,0.0));
    for(ii=0;ii<n;ii++)
      {
	for(jj=0;jj<d;jj++)
	  {
	    res0[ii][jj]=rnorm(0,1);
	    //Rprintf("res0[%d,%d]=%1.10f\t",ii,jj,res0[ii][jj]);
	  }
      }
    for(ii=0;ii<n;ii++)
      {
	for(jj=0;jj<d;jj++)
	  {
	    for(kk=0;kk<d;kk++)
	      {
		res[ii][jj]+=res0[ii][kk]*tmat[kk][jj];
	      }
	  }
      }
  }  
  void genericBoot(double* U,double* V,double* D,int* noMarkerPt, double* X_T_times_XPt, int* Npt, int* noBootPt, double* varYPt, double* statisticPt,int* alternativePt,int* testPt,double* extraPar)
  {
    double statData=extraPar[0];
    double alpha=extraPar[1];
    int d=noMarkerPt[0],ii=0,jj=0,kk=0;
    //Rprintf("no of marker%d",d);
    std::vector<std::vector<double> > tmat(d,std::vector<double>(d,0.0));
    std::vector<std::vector<double> > UDhalf(d,std::vector<double> (d,0.0));
    for(ii=0;ii<d;ii++)
      {
	for(jj=0;jj<d;jj++)
	  {
	    UDhalf[ii][jj]+=(matrixByCol(U,ii,jj,d,d)*sqrt(D[jj]));
	  }
      }
    for(ii=0;ii<d;ii++)
      {
	for(jj=0;jj<d;jj++)
	  {
	    for(kk=0;kk<d;kk++)
	      {
		tmat[ii][jj]+=(UDhalf[ii][kk]*matrixByCol(V,jj,kk,d,d));
	      }
	    //Rprintf("%d,%d,%8.3f ",ii,jj,tmat[ii][jj]);
	  }
      }
    std::vector<std::vector<double> > X_T_times_X(d,std::vector<double>(d,0.0));    
    for(jj=0;jj< noMarkerPt[0];jj++)
      {
        for(ii=0;ii< noMarkerPt[0];ii++)
          {
            X_T_times_X[ii][jj]=X_T_times_XPt[ii*noMarkerPt[0]+jj];
          }
      }
    int iiMax=noBootPt[0]/1000;
    
    if(testPt[0]==1)
      {
	double (*statFunc)(std::vector<double>&, std::vector<std::vector<double> >&, double, int,double* );
	statFunc=&rarecoverSumstatFunc;
	int noBootCompleted=0;
	double pValueCount=0.0,pValue=0.0,statPerm=0.0,pValueMin=0.0;
	GetRNGstate();	
	for(ii=0;ii<iiMax;ii++)
	  {
	    std::vector<std::vector<double> > X_T_times_Y_perm(1000,std::vector<double>(noMarkerPt[0],0.0));
	    rmvnorm(X_T_times_Y_perm,tmat);
	    for(jj=0;jj<1000;jj++)
	      {
		statPerm=(*statFunc)(X_T_times_Y_perm[jj],X_T_times_X,varYPt[0],alternativePt[0],extraPar);
		extraPar[3+noBootCompleted]=statPerm;
		noBootCompleted++;
		pValueCount+=(statPerm>statData ? 1.0 : 0.0);
	      }
	    pValue=pValueCount/((double) noBootCompleted);
	    pValueMin=pValue-1.96*sqrt(pValue*(1-pValue)/noBootCompleted);
	    if(pValueMin>alpha) ii=iiMax+1;
	  }
	extraPar[2]=pValue;
	PutRNGstate();
      }
    if(testPt[0]==3)
      {
	double (*statFunc)(std::vector<double>&, std::vector<std::vector<double> >&, double, int,double* );
	statFunc=&vtSumstatFunc;
	int noBootCompleted=0;
	double pValueCount=0.0,pValue=0.0,statPerm=0.0,pValueMin=0.0;
	GetRNGstate();
	for(ii=0;ii<iiMax;ii++)
	  {
	    std::vector<std::vector<double> > X_T_times_Y_perm(1000,std::vector<double>(noMarkerPt[0],0.0));
	    rmvnorm(X_T_times_Y_perm,tmat);
	    for(jj=0;jj<1000;jj++)
	      {
		statPerm=(*statFunc)(X_T_times_Y_perm[jj],X_T_times_X,varYPt[0],alternativePt[0],extraPar);
		//extraPar[3+noBootCompleted]=statPerm;
		//Rprintf("%d statperm %8.5f ",jj,statPerm);
		noBootCompleted++;
		pValueCount+=(statPerm>statData ? 1.0 : 0.0);
	      }
	    pValue=pValueCount/((double) noBootCompleted);
	    pValueMin=pValue-1.96*sqrt(pValue*(1-pValue)/noBootCompleted);
	    if(pValueMin>alpha) ii=iiMax+1;
	  }
	extraPar[2]=pValue;
	PutRNGstate();
      }
    if(testPt[0]==2)
      {
	double (*statFunc)(std::vector<double>&, std::vector<std::vector<double> >&, double, int,std::vector<double>&, std::vector<std::vector<double> >& );
	std::vector<double> weight(noMarkerPt[0],0.0);
	std::vector<std::vector<double> > lambdaMat(11,std::vector<double>(noMarkerPt[0],0.0));
	statFunc=&skatoSumstatFunc;

	for(ii=0;ii<noMarkerPt[0];ii++)
	  {
	    weight[ii]=extraPar[3+ii];
	  }
	for(ii=0;ii<11;ii++)
	  {
	    for(jj=0;jj<noMarkerPt[0];jj++)
	      {
		lambdaMat[ii][jj]=matrixByCol(&extraPar[noMarkerPt[0]+3],ii,jj,11,noMarkerPt[0]);
	      }
	  }
	int noBootCompleted=0;
	double pValueCount=0.0,pValue=0.0,statPerm=0.0,pValueMin=0.0;
	GetRNGstate();
	for(ii=0;ii<iiMax;ii++)
	  {
	    std::vector<std::vector<double> > X_T_times_Y_perm(1000,std::vector<double>(noMarkerPt[0],0.0));
	    rmvnorm(X_T_times_Y_perm,tmat);
	    for(jj=0;jj<1000;jj++)
	      {		
		statPerm=(*statFunc)(X_T_times_Y_perm[jj],X_T_times_X,varYPt[0],alternativePt[0],weight,lambdaMat);
		//Rprintf("%d,%8.3f,%8.3f\n",ii*1000+jj,statData,statPerm);
		statisticPt[jj]=statPerm;
		noBootCompleted++;
		pValueCount+=(statPerm<statData ? 1.0 : 0.0);
	      }
	    pValue=pValueCount/((double) noBootCompleted);
	    pValueMin=pValue-1.96*sqrt(pValue*(1-pValue)/noBootCompleted);
	    Rprintf("pValue,pValueMin:%8.3f,%8.3f\n",pValue,pValueMin);
	    if(pValueMin>alpha) ii=iiMax+1;
	  }
	extraPar[2]=pValue;
	PutRNGstate();
      }    
  }
  
  void genericStat(double* X_T_times_Y_permPt, double* X_T_times_XPt, int* noMarkerPt, int* Npt, int* noBootPt, double* varYPt, double* statisticPt, int* alternativePt,int* testPt,double* extraPar)
  {
    int ii=0,jj=0,kk=0;
    std::vector<std::vector<double> > X_T_times_Y_perm(noBootPt[0],std::vector<double>(noMarkerPt[0]));
    std::vector<std::vector<double> > X_T_times_X(noMarkerPt[0],std::vector<double>(noMarkerPt[0]));
    int alternative=alternativePt[0];
    for(jj=0;jj< noMarkerPt[0];jj++)
      {
	for(ii=0;ii< noMarkerPt[0];ii++)
	  {
	    X_T_times_X[ii][jj]=X_T_times_XPt[ii*noMarkerPt[0]+jj];
	  }
	for(kk=0;kk<noBootPt[0];kk++)
	  {
	    X_T_times_Y_perm[kk][jj]=X_T_times_Y_permPt[jj*noBootPt[0]+kk];
	  }
      }
    if(testPt[0]==1)
      {
	statisticPt[0]=rarecoverSumstatFunc(X_T_times_Y_perm[0],X_T_times_X,varYPt[0],alternative,extraPar);	
      }
    if(testPt[0]==3)
      {
	statisticPt[0]=vtSumstatFunc(X_T_times_Y_perm[0],X_T_times_X,varYPt[0],alternative,extraPar);	
      }
    if(testPt[0]==2)
      {
	std::vector<double> weight(noMarkerPt[0],0.0);
	std::vector<std::vector<double> > lambdaMat(11,std::vector<double>(noMarkerPt[0],0.0));
	for(ii=0;ii<noMarkerPt[0];ii++)
	  {
	    weight[ii]=extraPar[ii];
	  }
	for(ii=0;ii<11;ii++)
	  {
	    for(jj=0;jj<noMarkerPt[0];jj++)
	      {
		lambdaMat[ii][jj]=matrixByCol(&(extraPar[noMarkerPt[0]]),ii,jj,11,noMarkerPt[0]);
	      }
	  }
	statisticPt[0]=skatoSumstatFunc(X_T_times_Y_perm[0],X_T_times_X,varYPt[0],alternative,weight,lambdaMat);
	//Rprintf("Rprintf%8.3f",statisticPt[0]);
      }       
  }
  
  double quadForm(std::vector<double>& x,std::vector<std::vector<double> >& M,std::vector<double>& y)
  {
    double res=0.0;
    for(int ii=0;ii<x.size();ii++)
      {
	for(int jj=0;jj<y.size();jj++)
	  {
	    res+=(x[ii])*(M[ii][jj])*(y[jj]);
	  }
      }
    return res;
  }
  //good for R input;
  double matrixByCol(double* mat,int m,int n,int nrow,int ncol)
  {
    if(m <nrow & n<ncol)
      {
	return mat[m+n*nrow];
      }
    return(0.0);
  }
  /*  double vtSumstatFunc(std::vector<double >& X_T_times_Y, std::vector<std::vector<double> >& X_T_times_X,double varY, int alternative,double* extraPar)
  {
    int ii=0,noMarker=X_T_times_Y.size();
    std::vector<double> statVec(noMarker,0.0);
    for(ii=0;ii<X_T_times_Y.size();ii++)
      {
	statVec[ii]=(X_T_times_Y[ii])*(X_T_times_Y[ii])/(X_T_times_X[ii][ii])/varY;
      }
    return *max_element(statVec.begin(),statVec.end());
  }
  */
  double skatoSumstatFunc(std::vector<double >& X_T_times_Y, std::vector<std::vector<double> >& X_T_times_X,double varY, int alternative, std::vector<double>& weight,std::vector<std::vector<double> >& lambdaMat)
  {
    int ii=0,jj=0,kk=0,noMarker=X_T_times_Y.size();
    std::vector<double> rho(11,0.0),stat(11,0.0),pValue(11,0.0);    
    for(ii=0;ii<11;ii++) rho[ii]=0.1*ii;
    for(ii=0;ii<11;ii++)
      {
	for(jj=0;jj<noMarker;jj++)
	  {
	    for(kk=0;kk<noMarker;kk++)
	      {
		stat[ii]+=(1-rho[ii])*(X_T_times_Y[jj])*(X_T_times_Y[kk]);
		if(jj==kk) stat[ii]+=(rho[ii])*(X_T_times_Y[jj])*weight[kk]*(X_T_times_Y[kk]);
	      }
	  }
	pValue[ii]=pMixChisq(stat[ii],lambdaMat[ii]);
      }
    return (*min_element(pValue.begin(),pValue.end()));
  }
  //double (*statFunc)(std::vector<double>&, std::vector<std::vector<double> >&, double, int,double* );
  double vtSumstatFunc(std::vector<double >& X_T_times_Y, std::vector<std::vector<double> >& X_T_times_X,double varY, int alternative,double* extraPar)
  {
    int ii=0,jj=0,kk=0,noMarker=X_T_times_Y.size(),noTH=extraPar[3];
    std::vector<double> statVT(noTH,0.0);
    for(ii=0;ii<noTH;ii++)
      {
	double varStatVT=0.0;
	for(jj=0;jj<noMarker;jj++)
	  {
	    statVT[ii]+=(X_T_times_Y[jj]*matrixByCol(&(extraPar[4]),ii,jj,noTH,noMarker));
	  }
	statVT[ii]*=(statVT[ii]);
	for(jj=0;jj<noMarker;jj++)
	  {
	    for(kk=0;kk<noMarker;kk++)
	      {
		varStatVT+=(matrixByCol(&(extraPar[4]),ii,jj,noTH,noMarker)*(X_T_times_X[jj][kk])*matrixByCol(&(extraPar[4]),ii,kk,noTH,noMarker));
	      }
	  }
	statVT[ii]/=varStatVT;
      }
    return(*max_element(statVT.begin(),statVT.end()));
  }

  double rarecoverSumstatFunc(std::vector<double >& X_T_times_Y, std::vector<std::vector<double> >& X_T_times_X,double varY, int alternative,double* extraPar)
  {
    int cont=1,ii=0,jj=0;
    double Q=0.0;
    std::vector<double> ixIn(X_T_times_X.size(),0.0);
    std::vector<double> ixOut(X_T_times_X.size(),1.0);
    double statOld=0.0,statNewTemp=0.0;
    int noOut=X_T_times_X.size(), noIn=0,ixMax=0;
    std::vector<double> statVec(X_T_times_X.size(),0.0);
    double numer=0.0,denom=0.0;
    while(noOut>0 & cont==1)
      {
	for(ii=0;ii<ixOut.size();ii++)
	  {
	    if(ixOut[ii]==1.0)
	      {
		ixIn[ii]=1.0;
		denom=sqrt(quadForm(ixIn,X_T_times_X,ixIn)*varY);
		if(denom==0.0) denom=10000000.0;
		for(jj=0;jj<X_T_times_X.size();jj++)
		  {
		    statVec[ii]+=(ixIn[jj])*(X_T_times_Y[jj]);
		  } 
		statVec[ii]=statVec[ii]/denom;
		if(alternative==0) statVec[ii]*=statVec[ii]; 
		ixIn[ii]=0.0;
	      }     
	  }
	ixMax=max_element(statVec.begin(),statVec.end())-statVec.begin();
	statNewTemp=*max_element(statVec.begin(),statVec.end());
	if(statNewTemp-statOld>=Q)
	  {
	    ixIn[ixMax]=1.0;
	    ixOut[ixMax]=0.0;
	    statOld=statNewTemp;
	    cont=1;
	    noOut--;
	    for(ii=0;ii<statVec.size();ii++) 
	      {
		statVec[ii]=0.0;
	      }
	  }
	if(statNewTemp-statOld<Q)
	  cont=0;
      }
    return statOld;
  }
}
