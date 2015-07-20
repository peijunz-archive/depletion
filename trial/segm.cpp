#include <stdio.h>
#include <math.h>
typedef struct{
    double r;
    double p[2];
    double q[2];
} segm;
typedef segm* link;
link initsegm(double r, double* p,double* q){
    link pt=(link) malloc(sizeof(segm));
    pt->p[0]=p[0];
    pt->p[1]=p[1];
    if(r==0){
        pt->q[0]=q[0];
        pt->q[1]=q[1];
    }
    else if(r>0){

    }
    else{

    }
    pt->r=abs(r)
}
