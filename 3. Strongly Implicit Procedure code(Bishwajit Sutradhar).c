#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define nx 100
#define ny 100
int main()
{
    int n=nx*ny;
    FILE *f,*f1,*f2;
    f=fopen("SIP_CONTOUR_DATA.txt","w");
    f1=fopen("SIP_MID-VERTICAL.txt","w");
    f2=fopen("SIP_MID-HORIZONTAL.txt","w");
    double AW[n],AS[n],AP[n],AN[n],AE[n];
    double LW[n],LS[n],LP[n],UN[n],UE[n],MNW[n],MSE[n],Y[n],T[n],B[n];
    double NP[n],NN[n],NS[n],NE[n],NW[n];
    double Tfinal[(nx+2)*(ny+2)],Tanal[(nx+2)*(ny+2)];
    double sum,error=0.0000001,temp[n],BN[n],x,y,alpha = 0.1;
    double t1,t2,t3;
    int i,j;
    double dx;
    dx=2.0/(nx+1);

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
        if(i>=nx)
            LW[i]=AW[i]/(1+alpha*UN[i-nx]);

        LS[i]=AS[i];
        if(i>=1)
            LS[i]=AS[i]/(1+alpha*UE[i-1]);

        LP[i]=AP[i];
        if(i>0)
        {
            LP[i]=AP[i]+alpha*LS[i]*UE[i-1] - LS[i]*UN[i-1];
            if(i>=nx)
                LP[i]=AP[i]+alpha*(LW[i]*UN[i-nx] + LS[i]*UE[i-1]) - LW[i]*UE[i-nx] - LS[i]*UN[i-1];;
        }

        UN[i]=AN[i]/LP[i];
        if(i>=nx)
            UN[i]=(AN[i]-alpha*LW[i]*UN[i-nx])/LP[i];

        UE[i]=AE[i]/LP[i];
        if(i>0)
            UE[i]=(AE[i]-alpha*LS[i]*UE[i-1])/LP[i];
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
    //FINDING N MATRIX VALUES
    for(i=0;i<n;i++)
    {
        NP[i]=alpha*(MNW[i]+MSE[i]);
        NN[i]=-alpha*MNW[i];
        NS[i]=-alpha*MSE[i];
        NE[i]=-alpha*MSE[i];
        NW[i]=-alpha*MNW[i];
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
            x=((i/nx)+1)*dx;
            if(i%nx==nx-1)
                B[i]=-2.0*sin((M_PI*x)/2.0);
            else
                B[i]=0.0;
        }
       //finding b+Nt
        BN[0]=NP[0]*T[0] + NN[0]*T[1] + MSE[0]*T[nx-1] +NE[0]*T[nx];
        for(i=1;i<nx-1;i++)
            BN[i]=NS[i]*T[i-1]+NP[i]*T[i] + NN[i]*T[i+1] + MSE[i]*T[i+nx-1] +NE[i]*T[i+nx];
        BN[nx-1]=MNW[nx-1]*T[0]+NS[nx-1]*T[nx-2]+NP[nx-1]*T[nx-1] + NN[nx-1]*T[nx] + MSE[i]*T[2*nx-2] +NE[nx-1]*T[2*nx-1];
        for(i=nx;i<n-nx;i++)
            BN[i]=NW[i]*T[i-nx]+MNW[i]*T[i-nx+1]+ NS[i]*T[i-1]+NP[i]*T[i] + NN[i]*T[i+1] + MSE[i]*T[i+nx-1] +NE[i]*T[i+nx];
        for(i=n-nx;i<n-nx+1;i++)
            BN[i]=NW[i]*T[i-nx]+MNW[i]*T[i-nx+1]+ NS[i]*T[i-1]+NP[i]*T[i] + NN[i]*T[i+1] + MSE[i]*T[i+nx-1];
        for(i=n-nx+1;i<n-1;i++)
            BN[i]=NW[i]*T[i-nx]+MNW[i]*T[i-nx+1]+ NS[i]*T[i-1]+NP[i]*T[i] + NN[i]*T[i+1];
        for(i=n-1;i<n;i++)
            BN[i]=NW[i]*T[i-nx]+MNW[i]*T[i-nx+1]+ NS[i]*T[i-1]+NP[i]*T[i];
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
            x=i*dx;
            if(j==ny+1)
                Tfinal[i*(nx+2)+j]=2.0*sin((M_PI*x)/2.0);
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
        {
            t1=sin((M_PI*x)/2.0);
            t2=sinh((M_PI*y)/2.0);
            t3=sinh(M_PI);
        }
        Tanal[i]= 2.0*(t2/t3)*t1;
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
