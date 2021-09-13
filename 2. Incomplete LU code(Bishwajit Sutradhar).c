#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define nx 101
#define ny 101
int main()
{
    int n=nx*ny;
    FILE *f,*f1,*f2;
    f=fopen("ILU_CONTOUR_DATA.txt","w");
    f1=fopen("ILU_MID-VERTICAL.txt","w");
    f2=fopen("ILU_MID-HORIZONTAL.txt","w");
    double AW[n],AS[n],AP[n],AN[n],AE[n];
    double LW[n],LS[n],LP[n],UN[n],UE[n],MNW[n],MSE[n],Y[n],T[n],B[n];
    double Tfinal[(nx+2)*(ny+2)],Tanal[(nx+2)*(ny+2)];
    double sum,error=0.0000001,temp[n],BN[n],x,y;
    double ana,t1,t2,t3,t4;
    int i,j,N;
    double dx;
    dx=1.0/(nx+1);

    //finding AW-AE
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
    //FINDING LW-LP AND UN,UE
    for(i=0;i<n;i++)
    {
        LW[i]=AW[i];
        LS[i]=AS[i];
        LP[i]=AP[i];
        if(i>0)
        {
            LP[i]=AP[i]-LS[i]*UN[i-1];
            if(i>=nx)
                LP[i]=AP[i]-LW[i]*UE[i-nx]-LS[i]*UN[i-1];
        }
        UN[i]=AN[i]/LP[i];
        UE[i]=AE[i]/LP[i];
    }
    //FINDING MNW MSE
    for(i=0;i<n;i++)
    {
        MSE[i]=0.0;
        MNW[i]=0.0;
        if(i>0)
        {
            MSE[i]=LS[i]*UE[i-1];
            if(i>=nx)
                MNW[i]=LW[i]*UN[i-nx];
        }
    }
    //INTIAL SOLUTION
    for(i=0;i<n;i++)
    {
        T[i]=0.0;
        temp[i]=0.0;
    }
    //iteration
    do
    {
        //boundary condition
        for(i=0;i<n;i++)
        {
            if(i%nx==nx-1)
                B[i]=-1.0;
            else
                B[i]=0.0;
        }
       //finding b+Nt
        for(i=0;i<nx-1;i++)
            BN[i]=MSE[i]*T[i+nx-1];
        for(i=nx-1;i<n-nx+1;i++)
            BN[i]=MSE[i]*T[i+nx-1] + MNW[i]*T[i-nx+1];
        for(i=n-nx+1;i<n;i++)
            BN[i]=MNW[i]*T[i-nx+1];
        for(i=0;i<n;i++)
            B[i]=B[i]+BN[i];
        //finding y using forward substitution
        for(i=0;i<n;i++)
        {
            if(i<1)
                Y[i]=B[i]/LP[i];
            else if(i>=1 && i<nx)
                Y[i]=(B[i]-LS[i]*Y[i-1])/LP[i];
            else
                Y[i]=(B[i]-LS[i]*Y[i-1]-LW[i]*Y[i-nx])/LP[i];
        }
        //finding T using backward substitution
        for(i=n-1;i>=0;i--)
        {
            if(i==n-1)
                T[i]=Y[i];
            else if(i<n-1 && i>n-1-nx)
                T[i]=Y[i]-UN[i]*T[i+1];
            else
                T[i]=Y[i]-UN[i]*T[i+1]-UE[i]*T[i+nx];
        }
        //calculating error
        sum=0.0;
        for(i=0;i<n;i++)
            sum=sum+fabs(T[i]-temp[i]);
        //driving iteration
        for(i=0;i<nx*ny;i++)
            temp[i]=T[i];
    }while(sum>error);
    //arranging temperature values at all grid points
    for(i=0; i<(nx+2);i++)
    {
        for(j=0;j<ny+2;j++)
        {
            if(j==ny+1)
                Tfinal[i*(nx+2)+j]=1.0;
            else
            {
                if(i==0 || i==nx+1 ||j==0)
                    Tfinal[i*(nx+2)+j]=0.0;
                else
                    Tfinal[i*(nx+2)+j]=T[(i-1)*nx+j-1];
            }
        }
    }
    //analytical solution.
    for(i=0;i<(nx+2)*(ny+2);i++)
    {
        x=(i/(nx+2))*dx;
        y=(i%(nx+2))*dx;
        ana=0.0;
        for(N=1;N<100;N++)
        {
            t1=((pow(-1,N+1)+1)/(double)N);
            t2=sin((double)N*M_PI*x);
            t3=sinh((double)N*M_PI*y);
            t4=sinh(M_PI*(double)N);
            ana=ana+ (t1*t2*t3)/t4;
        }
        Tanal[i]=ana*(2.0/M_PI);
    }
    //writing the data
    fprintf(f,"x\ty\tT\tT(actual)\n");
    for(i=0;i<(nx+2)*(ny+2);i++)
    {
        x=(i/(nx+2))*dx;
        y=(i%(nx+2))*dx;
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
