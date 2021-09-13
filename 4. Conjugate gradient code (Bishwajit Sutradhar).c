#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define nx 101
#define ny 101
int main()
{
    int i,j,n=nx*ny;
    FILE *f,*f1,*f2;
    f=fopen("CG_CONTOUR_DATA.txt","w");
    f1=fopen("CG_MID-VERTICAL.txt","w");
    f2=fopen("CG_MID-HORIZONTAL.txt","w");
    double AW[n],AS[n],AP[n],AN[n],AE[n];
    double T[n],T1[n],B[n],rk[n],dk[n];
    double Ad[n],rk1[n],dk1[n],Be[n];
    double Tfinal[(nx+2)*(ny+2)],Tanal[(nx+2)*(ny+2)];
    double sum,error=1.0e-6;
    double beta = 1.0,dx,x,y,alpha;
    dx=1.0/(nx+1);

    //Forming AW-AE
    for(i=0;i<n;i++)
    {
        AP[i]=-4;
        {
            if(i<n-nx)
                AE[i]=1.0;
            else
                AE[i]=0.0;
        }
        {
            if(i>=nx)
                AW[i]=1.0;
            else
                AW[i]=0.0;
        }
        {
            if(i%nx==nx-1)
                AN[i]=0.0;
            else
                AN[i]=1.0;
        }
        {
            if(i>0)
            {
                if(i%nx==0)
                    AS[i]=0.0;
                else
                    AS[i]=1.0;
            }
            else
                AS[i]=0;
        }
    }
    //forming b vector
    for(i=0;i<n;i++)
    {
        x=((i/nx)+1)*dx;
        y=((i%nx)+1)*dx;
        B[i]=dx*dx*(-8.0*pow(M_PI,2))*sin(2.0*M_PI*x)*sin(2.0*M_PI*y);
    }
    //initialize
    double Ax;
    for(i=0;i<n;i++)
        T[i]=0.0;
    for(i=0;i<n;i++)
    {
        Ax=0.0;
        if(i==0)
            Ax=AP[i]*T[i] + AN[i]*T[i+1] +AE[i]*T[i+nx];
        else if(i>0 && i<nx)
            Ax=AS[i]*T[i-1]+AP[i]*T[i] + AN[i]*T[i+1] +AE[i]*T[i+nx];
        else if(i>=nx && i<nx*ny-nx)
            Ax=AW[i]*T[i-nx]+AS[i]*T[i-1]+AP[i]*T[i] + AN[i]*T[i+1] +AE[i]*T[i+nx];
        else if(i>=nx*ny-nx && i<nx*ny-1)
            Ax=AW[i]*T[i-nx]+AS[i]*T[i-1]+AP[i]*T[i] + AN[i]*T[i+1];
        else
            Ax=AW[i]*T[i-nx]+AS[i]*T[i-1]+AP[i]*T[i];
        rk[i]=B[i]-Ax;
        dk[i]=rk[i];
    }

    //iteration
    do
    {
        double rtd=0, dtad=0;
        for(i=0;i<n;i++)
        {
            Ad[i]=0.0;
            rtd=rtd+rk[i]*dk[i];
            if(i==0)
                Ad[i]=AP[i]*dk[i] + AN[i]*dk[i+1] +AE[i]*dk[i+nx];
            else if(i>0 && i<nx)
                Ad[i]=AS[i]*dk[i-1]+AP[i]*dk[i] + AN[i]*dk[i+1] +AE[i]*dk[i+nx];
            else if(i>=nx && i<nx*ny-nx)
                Ad[i]=AW[i]*dk[i-nx]+AS[i]*dk[i-1]+AP[i]*dk[i] + AN[i]*dk[i+1] +AE[i]*dk[i+nx];
            else if(i>=nx*ny-nx && i<nx*ny-1)
                Ad[i]=AW[i]*dk[i-nx]+AS[i]*dk[i-1]+AP[i]*dk[i] + AN[i]*dk[i+1];
            else
                Ad[i]=AW[i]*dk[i-nx]+AS[i]*dk[i-1]+AP[i]*dk[i];
            dtad=dtad+dk[i]*Ad[i];
        }
        alpha=(rtd/dtad);

        double rtr=0,dtr=0;
        for(i=0;i<n;i++)
        {
            T1[i]=T[i]+alpha*dk[i];
            rk1[i]=rk[i]-alpha*Ad[i];
        }
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                rtr=rtr+rk1[j]*rk1[j];
                dtr=dtr+dk[j]*rk[j];
            }
            Be[i]=rtr/dtr;
        }
        for(i=0;i<n;i++)
            dk1[i]=rk1[i]+Be[i]*dk[i];
        for(i=0;i<n;i++)
        {
            rk[i]=rk1[i];
            dk[i]=dk1[i];
        }
        //CALCULATING ERROR
        sum=0.0;
        for(i=0;i<n;i++)
            sum=sum+fabs(T1[i]-T[i]);
        for(i=0;i<n;i++)
            T[i]=T1[i];
    }while(sum>error);
    //arranging temperature values at all grid points
    for(i=0; i<(nx+2);i++)
    {
        for(j=0;j<ny+2;j++)
        {
            if(i==0 || i==nx+1 ||j==0 ||j==ny+1)
                Tfinal[i*(nx+2)+j]=0.0;
            else
                Tfinal[i*(nx+2)+j]=T1[(i-1)*nx+j-1];
        }
    }
    fprintf(f,"x\ty\tT\tT(actual)\n");
    for(i=0;i<(nx+2)*(ny+2);i++)
    {
        x=(i/(nx+2))*dx;
        y=(i%(nx+2))*dx;
        Tanal[i]=sin(2.0*M_PI*x)*sin(2.0*M_PI*y);
        fprintf(f,"%f\t%f\t%f\t%f\n",x,y,Tfinal[i],Tanal[i]);
    }
    fprintf(f1,"y\tT\tT(actual)\n");
    fprintf(f2,"x\tT\tT(actual)\n");
    for(i=0;i<nx+2;i++)
    {
        x=i*dx;
        y=i*dx;
        fprintf(f1,"%f\t%f\t%f\n",y,Tfinal[51*(nx+2)+i],Tanal[51*(nx+2)+i]);
        fprintf(f2,"%f\t%f\t%f\n",x,Tfinal[i*(nx+2)+51],Tanal[i*(nx+2)+51]);
    }
}
