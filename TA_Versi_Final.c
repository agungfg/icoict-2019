#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"
#define Nx 6400

double minDeltaT(double *deltaT,double *T, double *V, int N,double dx){
    int i;
    double minimumdeltaT;
    #pragma omp parallel for
    for(i=1;i<=N;i++){
            deltaT[i]    = 0.5 * ( dx  /  ((pow(T[i],0.5))  +  V[i])  )   ;
    }

    minimumdeltaT=deltaT[1];
    for(i=1;i<=N;i++){
        if(deltaT[i]<minimumdeltaT){
             minimumdeltaT=deltaT[i];
        }
    }
    return minimumdeltaT;
}

double leftBoundary(double *P,double *T, double *V, int N){
    P[1]      = 1;
    V[1]      = 2*V[2] - V[3];
    T[1]      = 1;
}

double rightBoundary(double *P,double *T, double *V, int N){
    P[N]      = 2*P[N-1]-P[N-2];
    V[N]      = 2*V[N-1]-V[N-2];
    T[N]      = 2*T[N-1]-T[N-2];
}

double leftBoundaryCorPre(double *P,double *T, double *V){
    P[1]      = 1;
    V[1]      = 2*V[2] - V[3];
    T[1]      = 1;

}

void Plot(FILE* pipe,double *P,double *T, double *V,double *x, int N){
    int i;
    fprintf(pipe, "set multiplot\n");
    fprintf(pipe, "plot '-' with lines lc rgb 'blue'\n");
    for(i=1;i<=N;i++){
            fprintf(pipe, "%g %g\n",x[i],P[i]);
    }
    fprintf(pipe, "end\n");
    fprintf(pipe, "plot '-' with lines lc rgb 'green'\n");
    for(i=1;i<=N;i++){
            fprintf(pipe, "%g %g\n",x[i],T[i]);
    }
    fprintf(pipe, "end\n");
    fprintf(pipe, "plot '-' with lines lc rgb 'red'\n");
    for(i=1;i<=N;i++){
            fprintf(pipe, "%g %g\n",x[i],V[i]);
    }
    fprintf(pipe, "end\n");
    fprintf(pipe, "pause 0.01\n");
    fprintf(pipe, "unset multiplot\n");
    fflush(pipe);
}

void Dealokasi(double *x, double *P_New, double *T_New, double *V_New,  double *A,  double *P_Pre, double *V_Pre, double *T_Pre,double *P_Cor, double *V_Cor, double *T_Cor, double *P_FinalPre, double *V_FinalPre, double *T_FinalPre){
    free(x);
    free(P_New);
    free(T_New);
    free(V_New);
    free(A);
    free(P_Pre);
    free(V_Pre);
    free(T_Pre);
    free(P_Cor);
    free(V_Cor);
    free(T_Cor);
    free(P_FinalPre);
    free(V_FinalPre);
    free(T_FinalPre);
}

int main(){

    //Deklarasi dan Alokasi
    double *x       		= (double *)malloc((Nx+1)*sizeof(double));
    double *A       		= (double *)malloc((Nx+1)*sizeof(double));
    double *deltaT  		= (double *)malloc((Nx+1)*sizeof(double));
    double *P_Cor      		= (double *)malloc((Nx+1)*sizeof(double));
    double *T_Cor      		= (double *)malloc((Nx+1)*sizeof(double));
    double *V_Cor      		= (double *)malloc((Nx+1)*sizeof(double));
    double *P_Pre      		= (double *)malloc((Nx+1)*sizeof(double));
    double *T_Pre      		= (double *)malloc((Nx+1)*sizeof(double));
    double *V_Pre      		= (double *)malloc((Nx+1)*sizeof(double));
    double *P_FinalPre 		= (double *)malloc((Nx+1)*sizeof(double));
    double *T_FinalPre 		= (double *)malloc((Nx+1)*sizeof(double));
    double *V_FinalPre 		= (double *)malloc((Nx+1)*sizeof(double));
    double *P_Old      		= (double *)malloc((Nx+1)*sizeof(double));
    double *T_Old      		= (double *)malloc((Nx+1)*sizeof(double));
    double *V_Old      		= (double *)malloc((Nx+1)*sizeof(double));
    double *P_New      		= (double *)malloc((Nx+1)*sizeof(double));
    double *T_New      		= (double *)malloc((Nx+1)*sizeof(double));
    double *V_New      		= (double *)malloc((Nx+1)*sizeof(double));
    int i;

    FILE* pipe = popen("gnuplot", "w");
    fprintf(pipe, "set yrange[0:3.5]\n");
    fprintf(pipe, "set xrange[0:3]\n");
    double wtime;
    double t;
    int L                   = 3;
    int N ;
    double Tfinal = 28.9516;

    double y ;
    double minimumdeltaT;
    double dx;

    int jumlahThread=1;
    double eksekusiWaktuSerial[100];
    double eksekusiWaktu2Thread[100];
    double eksekusiWaktu4Thread[100];
    double eksekusiWaktu8Thread[100];
    double gridx[10];
    gridx[1]   =   100;
    gridx[2]   =   200;
    gridx[3]   =   400;
    gridx[4]   =   800;
    gridx[5]   =   1600;
    gridx[6]   =   3200;
    gridx[7]   =   6400;

    int no=1;
    y                  = 1.4;

    while(jumlahThread<=8){

        dx						   = L/(gridx[no]-1);
        N                   	   = gridx[no];

        for(i=1;i<=N;i++){
            x[i]				   = (i-1)*dx;
            A[i]                   = 1 + 2.2*pow(x[i]-1.5,2);
            P_Old[i]               = 1 - 0.3146*x[i];
            T_Old[i]               = 1 - 0.2314*x[i];
            V_Old[i]               = (0.1+1.09*x[i])*pow(T_Old[i],0.5);
        }

        omp_set_num_threads(jumlahThread);
        wtime   =  omp_get_wtime();

        t       = 0.0;

        while(t<Tfinal){
            //Predictor
            #pragma omp parallel for
            for(i=2;i<=N-1;i++){
                P_Pre[i]    		=  (( -P_Old[i]*( V_Old[i+1] - V_Old[i] ) ) / dx ) - ( P_Old[i]*V_Old[i]*(log(A[i+1])-log(A[i]))/dx) - (V_Old[i]*(P_Old[i+1]-P_Old[i])/dx);
                V_Pre[i]    		=  (( -V_Old[i]*( V_Old[i+1] - V_Old[i] ) ) / dx ) - ( (1/y) * ( ( ( T_Old[i+1]- T_Old[i]) / dx ) + (T_Old[i]/P_Old[i])*(P_Old[i+1]-P_Old[i])/dx));
                T_Pre[i]    		= -V_Old[i]*(T_Old[i+1]-T_Old[i])/dx  - (y-1)*T_Old[i]*( (V_Old[i+1]- V_Old[i])/dx + V_Old[i]*(log(A[i+1])-log(A[i]))/dx);
            }
            minimumdeltaT=minDeltaT(deltaT,T_Old,V_Old,N,dx);
            //Final Predicted
            #pragma omp parallel for
            for(i=2;i<=N-1;i++){
                P_FinalPre[i]       = P_Old[i] + P_Pre[i]*minimumdeltaT;
                V_FinalPre[i]       = V_Old[i] + V_Pre[i]*minimumdeltaT;
                T_FinalPre[i]       = T_Old[i] + T_Pre[i]*minimumdeltaT;
            }
            leftBoundaryCorPre(P_FinalPre,T_FinalPre,V_FinalPre);
            //Corrector
            #pragma omp parallel for
            for(i=2;i<=N-1;i++){
                P_Cor[i]    		=  (( -P_FinalPre[i]*( V_FinalPre[i] - V_FinalPre[i-1] ) ) / dx ) - ( P_FinalPre[i]*V_FinalPre[i]*(log(A[i])-log(A[i-1]))/dx) - (V_FinalPre[i]*(P_FinalPre[i]-P_FinalPre[i-1])/dx);
                V_Cor[i]    		=  (( -V_FinalPre[i]*( V_FinalPre[i] - V_FinalPre[i-1] ) ) / dx ) - ( (1/y) * ( ( ( T_FinalPre[i]- T_FinalPre[i-1]) / dx ) + (T_FinalPre[i]/P_FinalPre[i])*(P_FinalPre[i]-P_FinalPre[i-1])/dx));
                T_Cor[i]    		= -V_FinalPre[i]*(T_FinalPre[i]-T_FinalPre[i-1])/dx  - (y-1)*T_FinalPre[i]*( (V_FinalPre[i]- V_FinalPre[i-1])/dx + V_FinalPre[i]*(log(A[i])-log(A[i-1]))/dx);
            }
            //Final Corrected
            #pragma omp parallel for
            for(i=2;i<=N-1;i++){
                P_New[i]      		= P_Old[i] + 0.5*(P_Pre[i]+P_Cor[i])*minimumdeltaT;
                V_New[i]      		= V_Old[i] + 0.5*(V_Pre[i]+V_Cor[i])*minimumdeltaT;
                T_New[i]      		= T_Old[i] + 0.5*(T_Pre[i]+T_Cor[i])*minimumdeltaT;
            }
            leftBoundary(P_New,T_New,V_New,N);
            rightBoundary(P_New,T_New,V_New,N);
            //Plot(pipe,P_New,T_New,V_New,x,N);
            t=t+minimumdeltaT;
            P_Old=P_New;
            T_Old=T_New;
            V_Old=V_New;
        }

        if(jumlahThread==1){
            eksekusiWaktuSerial[no]=omp_get_wtime()-wtime;
        }else if(jumlahThread==2){
            eksekusiWaktu2Thread[no]=omp_get_wtime()-wtime;
        }else if(jumlahThread==4){
            eksekusiWaktu4Thread[no]=omp_get_wtime()-wtime;
        }else if(jumlahThread==8){
            eksekusiWaktu8Thread[no]=omp_get_wtime()-wtime;
        }
        if(no==7){
            jumlahThread=jumlahThread*2;
            no=0;
        }
        no=no+1;
    }
    //Plot(pipe,P_New,T_New,V_New,x,N);
    printf("|----|----------|-------------|-------------|-------------|-------------|\n");
    printf("| No | Delta X  | Time Serial | Time 2 Core | Time 4 Core | Time 8 Core |\n");
    printf("|----|----------|-------------|-------------|-------------|-------------|\n");
    for(i=1;i<=7;i++){
        printf("|%3.0i |%9.0f |%12.8f |%12.8f |%12.8f |%12.8f |\n",i,gridx[i],eksekusiWaktuSerial[i],eksekusiWaktu2Thread[i],eksekusiWaktu4Thread[i],eksekusiWaktu8Thread[i]);
    }
    printf("|----|----------|-------------|-------------|-------------|-------------|\n");
    getchar();
    Dealokasi(x,P_New,T_New,V_New,A,P_Pre,V_Pre,T_Pre,P_Cor,V_Cor,T_Cor,P_FinalPre,V_FinalPre,T_FinalPre);
}
