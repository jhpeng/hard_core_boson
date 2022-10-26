#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

                    edge[7*i0+0] = ((t+1)%lt)*volume+z*area+y*lx+x;
                    edge[7*i0+1] = ((t+1)%lt)*volume+((z+1)%lz)*area+y*lx+x;
                    edge[7*i0+2] = ((t+1)%lt)*volume+z*area+((y+1)%ly)*lx+x;
                    edge[7*i0+3] = ((t+1)%lt)*volume+z*area+y*lx+((x+1)%lx);
                    edge[7*i0+4] = ((t+1)%lt)*volume+((z-1+lz)%lz)*area+y*lx+x;
                    edge[7*i0+5] = ((t+1)%lt)*volume+z*area+((y-1+ly)%ly)*lx+x;
                    edge[7*i0+6] = ((t+1)%lt)*volume+z*area+y*lx+((x-1+lx)%lx);
                }
            }
        }
    }
}

void worm_update(int* node, int* worldline, int* edge, int nsite, int nt, double epsilon, double mu, gsl_rng* rng) {
    int tail = (int)(gsl_rng_uniform_pos(rng)*nsite);
    int head = tail;
    int next;
    int check=1;

    double p0=exp(mu*epsilon);
    double p1=exp(-1.0*mu*epsilon);
    double pt=nt*epsilon;

    // direction : 0
    //      going up and create particle
    // direction : 1
    //      going down and destroy particle
    int direction=1;
    while(check) {
        if(node[head]==0) {
            if(gsl_rng_uniform_pos(rng)<p0) {
                direction=0;
            } else {
                direction=1;
            }
        } else {
            if(gsl_rng_uniform_pos(rng)<p1) {
                direction=1;
            } else {
                direction=0;
            }
        }

        if(direction==0 && worldline[2*head]==-1) {
            if(gsl_rng_uniform_pos(rng)<pt) {
                int index = (int)(gsl_rng_uniform_pos(rng)*nt);
                next = edge[(nt+1)*head+1+index];
            } else {
                next = edge[(nt+1)*head];
            }

            if(worldline[2*next+1]==-1) {
                worldline[2*head+0]=next;
                worldline[2*next+1]=head;
                head=next;
                node[head]=1;
            } else {
                int temp=worldline[2*next+1];
                worldline[2*head+0]=next;
                worldline[2*next+1]=head;
                head=temp;
                worldline[2*head+0]=-1;
            }
        }
        if(direction==1 && worldline[2*head+1]!=-1) {
            node[head]=0;
            next = worldline[2*head+1];
            worldline[2*head+1]=-1;
            worldline[2*next+0]=-1;
            head=next;
        }

        printf("%d %d \n",head,tail);

        if(head==tail) check=0;
    }
}

int main(int argc, char** argv) {
    // optain arguments
    Lx=atoi(argv[1]);
    Ly=atoi(argv[2]);
    Lz=atoi(argv[3]);
    Beta=atof(argv[4]);
    Epsilon=atof(argv[5]);
    Mu=atof(argv[6]);
    int thermal=atoi(argv[7]);
    int block_size=atoi(argv[8]);
    int nblock=atoi(argv[9]);
    int seed=atoi(argv[10]);


    // setup useful parameters
    Lt=(int)(Beta/Epsilon);
    Nsite=Lx*Ly*Lz*Lt;
    Nt=6;

    // allocate memory
    Edge = (int*)malloc(sizeof(int)*(1+Nt)*Nsite);
    Node = (int*)malloc(sizeof(int)*Nsite);
    WorldLine = (int*)malloc(sizeof(int)*2*Nsite);

    // initialization
    construct_cubic_lattice(Edge,Lx,Ly,Lz,Lt);
    for(int i=0;i<Nsite;i++) {
        Node[i]=0;
        WorldLine[2*i+0] = -1;
        WorldLine[2*i+1] = -1;
    }

    // setup random number generator
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    // thermaliztion
    for(int i=0;i<thermal;i++) {
        worm_update(Node,WorldLine,Edge,Nsite,Nt,Epsilon,Mu,rng);
    }

    for(int j=0;j<nblock;j++) {
        for(int i=0;i<block_size;i++) {
            worm_update(Node,WorldLine,Edge,Nsite,Nt,Epsilon,Mu,rng);
        }
    }

    for(int i=0;i<Nsite;i++) {
        if(Node[i]) {
            printf("%d %d \n",WorldLine[2*i],WorldLine[2*i+1]);
        }
    }
    // free memory
    free(Edge);
    free(Node);
    free(WorldLine);
}
