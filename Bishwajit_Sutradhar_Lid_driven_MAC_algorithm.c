#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define n 129

int main()
{
    FILE *f,*f1,*f2;
    int i,j,ite=0;
    double u[n][n+1]={0.},v[n+1][n]={0.},p[n+1][n+1]={0.};  //initializing
    double rhsu[n][n+1],rhsv[n+1][n];
    double uavg[n][n],vavg[n][n],pavg[n][n];
    double u21,u22,uv1,uv2,uv3,v21,v22,term1,term2;
    double psi[n][n],psinew[n][n];
    double unew,pnew,error=0.00001,re=400.0,dt=0.001;
    double dx=1.0/(n-1);
    double sumu,sump,sumpsi;
    //u BC (making avg value 1 at the boundary)
    for (i=0;i<=(n-1);i++)
    {
        u[i][n-1]=1.0;
        u[i][n]=1.0;
    }

    //Performing iteration
    printf("Performing iterations...\n\n");
    do
    {
        ite=ite+1;
        sumu=0.0;
        for(i=0;i<=(n-1);i++)
        {
            for(int j=0;j<=n;j++)
            {
                rhsu[i][j]=0.0;
            }
        }
        for(i=0;i<=n;i++)
        {
            for(j=0;j<=(n-1);j++)
            {
                rhsv[i][j]=0.0;
            }
        }
        //calculating rshu
        for(i=1;i<=(n-2);i++)
        {
            for(j=1;j<=(n-1);j++)
            {
                u21=(u[i][j]+u[i-1][j])*(u[i][j]+u[i-1][j])*0.25;
                u22=(u[i+1][j]+u[i][j])*(u[i+1][j]+u[i][j])*0.25;
                uv1=(u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])*0.25;
                uv2=(u[i][j]+u[i][j-1])*(v[i][j-1]+v[i+1][j-1])*0.25;
                term1=(dt/(dx*dx*re))*(u[i-1][j]-2.0*u[i][j]+u[i+1][j]);
                term2=(dt/(re*dx*dx))*(u[i][j-1]-2.0*u[i][j]+u[i][j+1]);
                rhsu[i][j]=u[i][j] - (dt/dx)*(u22-u21) - (dt/dx)*(uv1-uv2) + term1 + term2;
            }
        }
        //calculating rhsv
        for(i=1;i<=(n-1);i++)
        {
            for(j=1;j<=(n-2);j++)
            {
                uv1=(u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])*0.25;
                uv3=(u[i-1][j]+u[i-1][j+1])*(v[i][j]+v[i-1][j])*0.25;
                v21=(v[i][j]+v[i][j-1])*(v[i][j]+v[i][j-1])*0.25;
                v22=(v[i][j+1]+v[i][j])*(v[i][j+1]+v[i][j])*0.25;
                term1=(dt/(dx*dx*re))*(v[i-1][j]-2.0*v[i][j]+v[i+1][j]);
                term2=(dt/(re*dx*dx))*(v[i][j-1]-2.*v[i][j]+v[i][j+1]);
                rhsv[i][j]=v[i][j] - (dt/dx)*(uv1-uv3) - (dt/dx)*(v22-v21) + term1 + term2;
            }
        }
        //iteration for pressure
        do
        {
            sump=0.0;
            for(i=1;i<=(n-1);i++)
            {
                for(j=1;j<=(n-1);j++)
                {
                    pnew=0.25*((p[i+1][j]+p[i-1][j]+(p[i][j+1]+p[i][j-1])) - (dx*dx/dt)*((rhsu[i][j]-rhsu[i-1][j])/dx
                     + (rhsv[i][j]-rhsv[i][j-1])/dx));
                    sump=sump+fabs(pnew-p[i][j]);
                    p[i][j]=pnew;
                }
            }
            //BC for pressure
            for(i=1;i<=(n-1);i++)
            {
                p[i][0]=p[i][1];
                p[i][n]=p[i][n-1];
            }
            for(j=0;j<=n;j++)
            {
                p[n][j]=p[n-1][j];
                p[0][j]=p[1][j];
            }
        }while (sump>=0.001);
        // u velocity calculation
        for(i=1;i<=(n-2);i++)
        {
            for(j=1;j<=(n-1);j++)
            {
                unew=-(dt/dx)*(p[i+1][j]-p[i][j]) + rhsu[i][j];
                sumu=sumu+fabs(unew-u[i][j]);
                u[i][j]=unew;
            }
        }
        //BC for u velocity
        for(i=0;i<=(n-1);i++)
        {
            u[i][n]=2.0-u[i][n-1];
            u[i][0]=-u[i][1];
        }
        //v velocity calculation
        for(i=1;i<=(n-1);i++)
        {
            for(j=1;j<=(n-2);j++)
            {
                v[i][j]=-(dt/dx)*(p[i][j+1]-p[i][j]) + rhsv[i][j];
            }
        }
        //BC for v velocity
        for(j=0;j<=(n-1);j++)
        {
            v[n][j]=-v[n-1][j];
            v[0][j]=-v[1][j];
        }
        if(ite%1000==0)
            printf("iteration=%d\t error=%f\n",ite,sumu);
    } while (sumu>error);

    //finding avg(interpolating) at actual grid points
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            uavg[i][j]=0.5*(u[i][j]+u[i][j+1]);
            vavg[i][j]=0.5*(v[i][j]+v[i+1][j]);
            pavg[i][j]=0.25*(p[i][j]+p[i+1][j]+p[i+1][j+1]+p[i][j+1]);
        }
    }
    //stream function calculation
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
            psi[i][j]=0.0;
            psinew[i][j]=0.0;
    }
    do
    {
        sumpsi=0.0;
        for(i=0;i<n;i++)
        {
            for(j=2;j<n;j++)
                psinew[i][j]=2.0*dx*uavg[i][j-1]+psi[i][j-2];
        }
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                sumpsi=sumpsi+fabs(psinew[i][j]-psi[i][j]);
            }
        }
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
                psi[i][j]=psinew[i][j];
        }
    }while(sumpsi>error);
    f=fopen("contours_MAC.txt","w");
    fprintf(f,"x\ty\tpsi\tu\tv\tp\n");
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
            fprintf(f,"%f\t%f\t%f\t%f\t%f\t%f\n",i*dx,j*dx,psi[i][j],uavg[i][j],vavg[i][j],pavg[i][j]);
    }
    f1=fopen("u_centre_MAC.txt","w");
    fprintf(f1,"u\ty\n");
    for(i=0;i<n;i++)
        fprintf(f1,"%f\t%f\n",uavg[64][i],i*dx);
    f2=fopen("v_centre_MAC.txt","w");
    fprintf(f2,"v\ty\n");
    for(i=0;i<n;i++)
        fprintf(f2,"%f\t%f\n",vavg[i][64],i*dx);
    fclose(f);
    fclose(f1);
    fclose(f2);
    return 0;
}
