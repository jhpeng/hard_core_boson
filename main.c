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
    int nt=6;
    for(int t=0;t<lt;t++) {
        for(int z=0;z<lz;z++) {
            for(int y=0;y<ly;y++) {
                for(int x=0;x<lx;x++) {
                    int i0 = t*volume+z*area+y*lx+x;

                    edge[(nt+1)*i0+0] = ((t+1)%lt)*volume+z*area+y*lx+x;
                    edge[(nt+1)*i0+1] = ((t+1)%lt)*volume+((z+1)%lz)*area+y*lx+x;
                    edge[(nt+1)*i0+2] = ((t+1)%lt)*volume+z*area+((y+1)%ly)*lx+x;
                    edge[(nt+1)*i0+3] = ((t+1)%lt)*volume+z*area+y*lx+((x+1)%lx);
                    edge[(nt+1)*i0+4] = ((t+1)%lt)*volume+((z-1+lz)%lz)*area+y*lx+x;
                    edge[(nt+1)*i0+5] = ((t+1)%lt)*volume+z*area+((y-1+ly)%ly)*lx+x;
                    edge[(nt+1)*i0+6] = ((t+1)%lt)*volume+z*area+y*lx+((x-1+lx)%lx);

                    if(0) {
                        printf("%d | %d : %d | ",t,i0%volume,(t+1)%lt);
                        for(int i=0;i<7;i++) {
                            printf("%d ",edge[(nt+1)*i0+i]%volume);
                        }
                        printf("\n");
                    }
                }
            }
        }
    }
}

void print_conf(int* node, int* worldline, int nsite) {
    printf("------------------------------------\n");
    for(int i=0;i<nsite;i++) {
        printf("i=%d \t %d | [%d,%d] \n",i,node[i],worldline[2*i],worldline[2*i+1]);
    }
}

static int number_of_moving;

#if(0)
void worm_update(int* node, int* worldline, int* edge, int nsite, int nt, gsl_rng* rng) {
    double epsilon = Epsilon;
    double mu = Mu;

    int tail = (int)(gsl_rng_uniform_pos(rng)*nsite);
    int head = tail;
    int next;
    int check=1;

    double p0=exp(mu*epsilon);
    double p1=exp(-1.0*mu*epsilon);
    double pt=nt*epsilon;

    number_of_moving=-1;

    // direction : 0
    //      moving forward in temperal direction
    // direction : 1
    //      moving backward in temperal direction
    int direction=0;
    if(node[tail]) direction=1;
    while(check) {
        if(direction) {
            if(gsl_rng_uniform_pos(rng)>p1) direction=0;
        } else {
            if(gsl_rng_uniform_pos(rng)>p0) direction=1;
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
            next = worldline[2*head+1];
            worldline[2*head+1]=-1;
            worldline[2*next+0]=-1;
            
            if(worldline[2*head+0]==-1) node[head]=0;
            head=next;
        }

        //printf("%d %d \n",head,tail);

        if(head==tail) check=0;

        direction = (int)(gsl_rng_uniform_pos(rng)*2);
        number_of_moving++;
    }

    printf("%d\n",number_of_moving);
}
#endif

void worm_update(int* node, int* worldline, int* edge, int nsite, int nt, gsl_rng* rng) {
    double epsilon = Epsilon;
    double mu = Mu;

    int tail = (int)(gsl_rng_uniform_pos(rng)*nsite);
    int head = tail;
    int next;
    int check=1;

    double p0=exp(mu*epsilon);
    double p1=exp(-1.0*mu*epsilon);
    double pt=nt*epsilon;

    number_of_moving=-1;

    // direction : 0
    //      moving forward in temperal direction
    // direction : 1
    //      moving backward in temperal direction
    int direction=0;
    if(node[tail]) direction=1;
    while(check) {
        if(direction) {
            if(gsl_rng_uniform_pos(rng)<p1) {
                if(worldline[2*head+1]!=-1) {
                    next = worldline[2*head+1];
                    worldline[2*head+1]=-1;
                    worldline[2*next+0]=-1;
            
                    if(worldline[2*head+0]==-1) node[head]=0;
                        head=next;
                }
            }
        } else {
            if(gsl_rng_uniform_pos(rng)<p0) {
                if(worldline[2*head]==-1) {
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
            } else if(worldline[2*head+1]!=-1) {
                next = worldline[2*head+1];
                worldline[2*head+1]=-1;
                worldline[2*next+0]=-1;
            
                if(worldline[2*head+0]==-1) node[head]=0;
                    head=next;
            }
        }


        //printf("%d %d \n",head,tail);

        if(head==tail) check=0;

        direction = (int)(gsl_rng_uniform_pos(rng)*2);
        number_of_moving++;
    }
}

static double nparticle_ave=0;
static double energy_ave=0;
static double energy_t_ave=0;
static double winding_square_x_ave=0;
static double winding_square_y_ave=0;
static double winding_square_z_ave=0;

static int measure_count=0;
void measure(int* node, int* worldline, int* edge, int nsite, int nt, int block_size) {
    double epsilon = Epsilon;
    double beta = Beta;
    double mu = Mu;

    int nhx=0;
    int nhy=0;
    int nhz=0;
    int nht=0;
    int nhs=0;

    for(int i=0;i<nsite;i++) {
        if(node[i]) {
            int k=0;
            for(int j=0;j<(nt+1);j++) {
                if(worldline[2*i]==edge[(nt+1)*i+j]) {
                    k=j;
                    break;
                }
            }

            if(k==0) {
                nht++;
            } else if(k==1) {
                nhz++;
                nhs++;
            } else if(k==2) {
                nhy++;
                nhs++;
            } else if(k==3) {
                nhx++;
                nhs++;
            } else if(k==4) {
                nhz--;
                nhs++;
            } else if(k==5) {
                nhy--;
                nhs++;
            } else if(k==6) {
                nhx--;
                nhs++;
            }
        }
    }

    int lt=Lt;
    int lx=Lx;
    int ly=Ly;
    int lz=Lz;

    int volume=lx*ly*lz;
    int nparticle=0;
    for(int i=0;i<volume;i++) {
        nparticle += node[i];
    }

    double energy_t = -(1.0/beta)*nhs+nt*nht/(lt*(1-nt*epsilon));
    double energy = energy_t-mu*nparticle;
    //double energy_t = -(1.0/beta)*nhs+nt*nparticle;
    //double energy = -(1.0/beta)*nhs+(nt-mu)*nparticle;
    double wx = nhx/lx;
    double wy = nhy/ly;
    double wz = nhz/lz;

    nparticle_ave += nparticle;
    energy_ave += energy;
    energy_t_ave += energy_t;
    winding_square_x_ave += wx*wx;
    winding_square_y_ave += wy*wy;
    winding_square_z_ave += wz*wz;

    measure_count++;

    if(measure_count==block_size) {
        nparticle_ave = nparticle_ave/block_size;
        energy_ave = energy_ave/block_size;
        energy_t_ave = energy_t_ave/block_size;
        winding_square_x_ave = winding_square_x_ave/block_size;
        winding_square_y_ave = winding_square_y_ave/block_size;
        winding_square_z_ave = winding_square_z_ave/block_size;

        // print to stdout
        if(1) {
            printf("----------------------------\n");
            printf("n   : %.12e \n",nparticle_ave);
            printf("e   : %.12e \n",energy_ave);
            printf("et  : %.12e \n",energy_t_ave);
            printf("wx2 : %.12e \n",winding_square_x_ave);
            printf("wy2 : %.12e \n",winding_square_y_ave);
            printf("wz2 : %.12e \n",winding_square_z_ave);
        }
        
        // print to file
        FILE* dfile = fopen("data.txt","a");
        fprintf(dfile,"%.12e %.12e %.12e %.12e %.12e %.12e \n",nparticle_ave,energy_ave,energy_t_ave,winding_square_x_ave,winding_square_y_ave,winding_square_z_ave);
        fclose(dfile);


        nparticle_ave = 0;
        energy_ave = 0;
        energy_t_ave = 0;
        winding_square_x_ave = 0;
        winding_square_y_ave = 0;
        winding_square_z_ave = 0;

        measure_count=0;
    }
}


int main(int argc, char** argv) {
    // optain arguments
    Lx=atoi(argv[1]);
    Ly=atoi(argv[2]);
    Lz=atoi(argv[3]);
    Mu=atof(argv[4]);
    Beta=atof(argv[5]);
    Epsilon=atof(argv[6]);
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
        worm_update(Node,WorldLine,Edge,Nsite,Nt,rng);
    }

    for(int j=0;j<nblock;j++) {
        for(int i=0;i<block_size;i++) {
            for(int k=0;k<10;k++)
                worm_update(Node,WorldLine,Edge,Nsite,Nt,rng);

            measure(Node,WorldLine,Edge,Nsite,Nt,block_size);
        }
    }

    if(0) {
        worm_update(Node,WorldLine,Edge,Nsite,Nt,rng);
        print_conf(Node,WorldLine,Nsite);
        worm_update(Node,WorldLine,Edge,Nsite,Nt,rng);
        print_conf(Node,WorldLine,Nsite);
        worm_update(Node,WorldLine,Edge,Nsite,Nt,rng);
        print_conf(Node,WorldLine,Nsite);
    }

    // free memory
    gsl_rng_free(rng);
    free(Edge);
    free(Node);
    free(WorldLine);
}
