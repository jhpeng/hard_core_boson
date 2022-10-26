#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

// physics parameter
double Beta,Epsilon,Mu;

// system size
int Lx,Ly,Lz,Lt;
int Nsite;

// lattice parameter
int Nt;
int* Edge=NULL;

// dynamic variables
int* Node=NULL;
int* WorldLine=NULL;


void construct_cubic_lattice(int* edge, int lx, int ly, int lz, int lt) {
    int volume=lx*ly*lz;
    int area=lx*ly;
    for(int t=0;t<lt;t++) {
        for(int z=0;z<lz;z++) {
            for(int y=0;y<ly;y++) {
                for(int x=0;x<lx;x++) {
                    int i0 = t*volume+z*area+y*lx+x;

                    edge[i0+0] = ((t+1)%lt)*volume+z*area+y*lx+x;
                    edge[i0+1] = ((t+1)%lt)*volume+((z+1)%lz)*area+y*lx+x;
                    edge[i0+2] = ((t+1)%lt)*volume+z*area+((y+1)%ly)*lx+x;
                    edge[i0+3] = ((t+1)%lt)*volume+z*area+y*lx+((x+1)%lx);
                    edge[i0+4] = ((t+1)%lt)*volume+((z-1+lz)%lz)*area+y*lx+x;
                    edge[i0+5] = ((t+1)%lt)*volume+z*area+((y-1+ly)%ly)*lx+x;
                    edge[i0+6] = ((t+1)%lt)*volume+z*area+y*lx+((x-1+lx)%lx);

                    edge[i0+7] = ((t-1+lt)%lt)*volume+z*area+y*lx+x;
                    edge[i0+8] = ((t-1+lt)%lt)*volume+((z+1)%lz)*area+y*lx+x;
                    edge[i0+9] = ((t-1+lt)%lt)*volume+z*area+((y+1)%ly)*lx+x;
                    edge[i0+10] = ((t-1+lt)%lt)*volume+z*area+y*lx+((x+1)%lx);
                    edge[i0+11] = ((t-1+lt)%lt)*volume+((z-1+lz)%lz)*area+y*lx+x;
                    edge[i0+12] = ((t-1+lt)%lt)*volume+z*area+((y-1+ly)%ly)*lx+x;
                    edge[i0+13] = ((t-1+lt)%lt)*volume+z*area+y*lx+((x-1+lx)%lx);
                }
            }
        }
    }
}

//void worm_update(int* node, int* worldline, int* edge, int nsite, int nt, gsl_rng* rng) {
//    int start_point = (int)(gsl_rng_uniform_pos(rng)*nsite);
//}

int main(int argc, char** argv) {
    // optain arguments
    Lx=atoi(argv[1]);
    Ly=atoi(argv[2]);
    Lz=atoi(argv[3]);
    Beta=atof(argv[4]);
    Epsilon=atof(argv[5]);

    // setup useful parameters
    Lt=(int)(Beta/Epsilon);
    Nsite=Lx*Ly*Lz*Lt;
    Nt=6;

    // allocate memory
    Edge = (int*)malloc(sizeof(int)*2*(1+Nt)*Nsite);
    Node = (int*)malloc(sizeof(int)*Nsite);
    WorldLine = (int*)malloc(sizeof(int)*2*Nsite);

    // initialization
    construct_cubic_lattice(Edge,Lx,Ly,Lz,Lt);
    for(int i=0;i<Nsite;i++) {
        Node[i]=0;
        WorldLine[2*i+0] = -1;
        WorldLine[2*i+1] = -1;
    }



    // free memory
    free(Edge);
    free(Node);
    free(WorldLine);
}
