#include <R.h>
#include "statist.h"

using namespace std;

#ifndef QUAD
#define QUAD(x) ((x)*(x))
#endif
#ifndef GAUSS
#define GAUSS(x,mu,sigma) (1/(sqrt(2*PI)*sigma)*exp(-QUAD((x-mu)/(1.0*sigma))/2.0))
#endif

double lts (double *const data,
	    int    n,
	    int    h,
	    int    *const drinnen)
{
    double minmean, mindist, tempmean, tempdist;
    int index[n];
    int i, j, start;

    for (i=0; i<n; i++)
	index[i]=i;
    
    indexsort(index,data,n,0,n-1,1);

    minmean = indexmean(index,data,n-h);
    mindist = 0;
    start  = 0;
    for (j=0;j<n-h;j++)
	mindist += QUAD(data[index[j]]-minmean);
	
    for (i=1;i<h+1;i++){
	tempmean = indexmean(index,data+i,n-h);
	tempdist = 0;
	for (j=0;j<n-h;j++)
	    tempdist += QUAD(data[index[i+j]]-tempmean);
	
	if (tempdist < mindist){
	    mindist=tempdist;
	    minmean=tempmean;
	    start=i;
	}
    }

    for (i=0;i<n;i++){
	if (i>=start && i < start+n-h)
	    drinnen[index[i]]=1;
	else
	    drinnen[index[i]]=0;
    }

    return minmean;
}

extern "C" {

    void c_cggm (double *const b,
		 int    *const nrow,
		 int    *const ncol,
		 double *const g,
		 double *const h,
		 double *const result)
    {
	int i,j,k,u,v,richtung,richtungalt,w,n;
	int z, s;
	double y,f0,f1,f2,fhilf,diff,halbschritt;

	n = ((*nrow)>(*ncol)?(*nrow):(*ncol));
	w = (int) floor((*h)*n);

	// Kopie des Bildes mit vergrößertem Rand
	double eb[(*nrow)+2*w][(*ncol)+2*w];

	z=(*nrow);
	s=(*ncol);
	// vergrößerten Rand berechnen
	for(i=0; i<(*nrow)+2*w; i++)
	    for(j=0; j<(*ncol)+2*w; j++)
		eb[i][j]=
		    b[(i<=w ? 0 : (i>=(*nrow)+w ? (*nrow)-1 : i-w))
		       +((*nrow)*
		     ((j<=w ? 0 : (j>=(*ncol)+w ? (*ncol)-1 : j-w))))];

	for (i=0; i<z; i++){
	    Rprintf("Row: %i / %i\n",i+1,z);
	    for (j=0; j<s; j++){
		y=b[i+(*nrow)*j];
		diff=1;
		richtung=0;
		for(k=0;k<=30 && fabs(diff)>0.01;k++){
		    f0=f1=f2=0;
		    for(u=0; u<=2*w; u++)
			for(v=0; v<=2*w; v++){
			    fhilf=-GAUSS(eb[i+u][j+v],y,(*g))*
				GAUSS(u,w,(*h)*n)*GAUSS(v,w,(*h)*n);
			    f0+=fhilf;
			    //erste Ableitung
			    f1+=(eb[i+u][j+v]-y)/(1.0*QUAD(*g))*fhilf;   
			    // zweite Ableitung
			    f2+=(-1.0/QUAD(*g)+QUAD(eb[i+u][j+v]-y)/
				 (1.0*pow((*g),4)))*fhilf;
			}
		    
		    if(f2<=0){ //concave case
			richtung=(f1>0?-1:1);
			diff=richtung*(*g);
			y+=diff;
		    }
		    else{  //convex case
			richtungalt=richtung;
			richtung= -f1/f2>0 ? 1 : -1;
			halbschritt= richtung*richtungalt>=0 ? 1 : 0.5;
			if(fabs(f1/f2)>(*g)){
			    diff=halbschritt*richtung*(*g);
			    y+=diff;
			}
			else{
			    diff=-halbschritt*f1/f2;
			    y+=diff;
			}
		    }
		}
		result[i+(*nrow)*j]=y;
	    }
	}
	
	return;
    }
    
    void c_cggm_lts (double *const b,
		     int    *const nrow,
		     int    *const ncol,
		     double *const g,
		     double *const h,
		     double *const anteil,
		     double *const result)
    {
	int i,j,k,u,v,richtung,richtungalt,ant;
	int z, s, w,n;
	double y,f0,f1,f2,fhilf,diff,halbschritt;

	n = ((*nrow)>(*ncol)?(*nrow):(*ncol));
	w = (int) floor((*h)*n);

	// Kopie des Bildes mit vergrößertem Rand
	double eb[(*nrow)+2*w][(*ncol)+2*w];

	double fenster[(2*w+1)*(2*w+1)];
	int drinnen[(2*w+1)*(2*w+1)];
	
	z=(*nrow);
	s=(*ncol);
	// vergrößerten Rand berechnen
	for(i=0; i<(*nrow)+2*w; i++)
	    for(j=0; j<(*ncol)+2*w; j++)
		eb[i][j]=
		    b[(i<=w ? 0 : (i>=(*nrow)+w ? (*nrow)-1 : i-w))
		       +((*nrow)*
		     ((j<=w ? 0 : (j>=(*ncol)+w ? (*ncol)-1 : j-w))))];
	

	ant=(int) floor((*anteil)*QUAD((2*w+1)));

	for (i=0; i<z; i++){
	    Rprintf("Row: %i / %i\n",i+1,z);
	    for (j=0; j<s; j++){
		
		for(u=0; u<2*w+1; u++)
		    for(v=0; v<2*w+1; v++)
			fenster[(2*w+1)*u+v]=eb[i+u][j+v];
		
		lts(fenster,QUAD((2*w+1)),ant,drinnen);

		y=b[i+(*nrow)*j];
		diff=1;
		richtung=0;
		for(k=0;k<=30 && fabs(diff)>0.01;k++){
		    f0=f1=f2=0;
    
		    for(u=0; u<2*w+1; u++)
			for(v=0; v<2*w+1; v++){
			    if(drinnen[(2*w+1)*u+v]==1){
				fhilf=-GAUSS(eb[i+u][j+v],y,(*g))*
				    GAUSS(u,w,(*h)*n)*GAUSS(v,w,(*h)*n);
				f0+=fhilf;
				//erste Ableitung
				f1+=(eb[i+u][j+v]-y)/(1.0*QUAD(*g))*fhilf;  
				// zweite Ableitung
				f2+=(-1.0/QUAD(*g)+QUAD(eb[i+u][j+v]-y)/
				     (1.0*pow(*g,4)))*fhilf;
			    }
			}
		    
		    if(f2<=0){ //concave case
			richtung=(f1>0?-1:1);
			diff=richtung*(*g);
			y+=diff;
		    }
		    else{  //convex case
			richtungalt=richtung;
			richtung= -f1/f2>0 ? 1 : -1;
			halbschritt= richtung*richtungalt>=0 ? 1 : 0.5;
			if(fabs(f1/f2)>(*g)){
			    diff=halbschritt*richtung*(*g);
			    y+=diff;
			}
			else{
			    diff=-halbschritt*f1/f2;
			    y+=diff;
			}
		    }
		}
		result[i+(*nrow)*j]=y;
	    }
	}
	return;
    }
}
