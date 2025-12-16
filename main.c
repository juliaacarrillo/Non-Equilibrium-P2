#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define pi 3.14159265358979323846264
/* VARIABLES GENERADOR DE PARISI-RAPUANO */
#define Frec_Max 1000
#define NormRANu (2.3283063671E-10F)
#define N_data 5000
#define NMAX 50

/* PARAMETROS DE LA SIMULACION */
#define NUM_STEPS 5000        // Número total de pasos (ANTES 100000), como pide el enunciado (100 * dt).
#define L 100
#define N_particles 1000
#define delta_t 0.1
#define gamma 1.0
#define D 1.0 // Coeficiente de difusión
#define NSIM 20 // Número de simulaciones para promediar estadísticamente el MSD (antes era 25)

#define NF 4 
double forces[NF]={0.001,0.01,0.1,1.0};

/* FUNCTIONS */
unsigned int irr[256]; unsigned int ir1; unsigned char ind_ran,ig1,ig2,ig3;
float Random(void ); void ini_ran(int );
void BoxMuller (double* num1, double* num2);
void Normal_distribution(const char *filename);
double PBC(double pos);
void Brownian_N(double tau, int flag_ini, const char *filename );
void Brownian_MSD(double tau, int flag_ini, double *MSD_total);
void Brownian_MSD_force(double tau, int flag_ini, double f, double *eps, double *MSD_total, double *DX_total);

int main(){
    // Semilla del generador aleatorio
    ini_ran(time(NULL));
    double taus[] = {80.0,90.0,100.0};
    int ntaus= sizeof(taus)/sizeof(taus[0]);
    int j; 
    double forces_mu[] = {0.001, 1.0};       // fuerzas para μ
    double forces_D[] = { 0.001, 0.01, 0.1, 1.0}; // fuerzas para D si quieres probar
    int Nf_mu = sizeof(forces_mu)/sizeof(forces_mu[0]);
    int Nf_D  = sizeof(forces_D)/sizeof(forces_D[0]);
    int flag_ini=1;

    //int sim,i,fi; 
    char name_msd[NMAX], name_dx[NMAX];
    // Para acomular el valor promediado
    //double MSD_total[NUM_STEPS+1] = {0.0};
    char output_filename[NMAX]; // Un buffer para guardar el string

    /* SIMULATION */

    // BUCLE EN TAUS
    FILE *f;
    /*
    for(j=0;j<ntaus;j++){
        // Inicializar msd
        for(int i=0;i<=NUM_STEPS;i++){
            MSD_total[i]=0.0;
        }

        for (int i=0;i<NSIM;i++){
            Brownian_MSD(taus[j],flag_ini,MSD_total); 
        }
        // Guardamos el promedio en un archivo
        sprintf(output_filename, "MSD_Taus_%.1f.txt", taus[j]);
        f = fopen(output_filename,"w");

        for (int i=0;i<=NUM_STEPS;i++){
            double t = i*delta_t;
            MSD_total[i]/=NSIM; // Promediamos
            fprintf(f,"%lf\t%lf\n",t,MSD_total[i]);
        }

    }
    fclose(f);
    printf("¡Promedio de %d simulaciones completado!\n", NSIM); 
    for(j=0;j<ntaus;j++){
        sprintf(output_filename, "trajectories_tau_%.1f.txt", taus[j]);
        printf("Simulando con tau = %.2f -> %s\n", taus[j], output_filename);
        Brownian_N(taus[j],flag_ini,output_filename);
    }
    */
    /*
    for (fi=0;fi<NF;fi++){
        double f= forces[fi];
        printf("Simulando f = %g === \n",f);

        double MSD_tot[NUM_STEPS+1];
        double DX_tot[NUM_STEPS+1];
        for (i=0;i<=NUM_STEPS;i++) { MSD_tot[i]=0.0; DX_tot[i]=0.0;}

        for (sim = 0; sim<NSIM; sim++){

            double *eps= malloc (sizeof(double)*N_particles);
            for(int p=0; p<N_particles;p++){
                eps[p]= (Random()<0.5)? -1.0 :1.0; 
       
            }
            Brownian_MSD_force(10.0,1,f,eps,MSD_tot,DX_tot);
            free(eps);
        }
        for(i=0;i<=NUM_STEPS;i++){
            MSD_tot[i]/=(double)NSIM;
            DX_tot[i]/=(double)NSIM;
        }
        // Guardar resultados en archivos (dos archivos por valor de f)
        sprintf(name_msd, "MSD_x_f_%.3f.txt", f);
        sprintf(name_dx,  "DX_f_%.3f.txt",  f);

        FILE *fm = fopen(name_msd, "w");
        FILE *fd = fopen(name_dx,  "w");

        if (!fm || !fd) { perror("Error opening output files"); exit(1); }

        for (i = 0; i <= NUM_STEPS; i++) {
            double t = i * delta_t;
            fprintf(fm, "%g\t%g\n", t, MSD_tot[i]);
            fprintf(fd, "%g\t%g\n", t, DX_tot[i]);
        }
        fclose(fm); fclose(fd);
        printf("Guardado: %s , %s\n", name_msd, name_dx);
    }*/
   // FORCES AND TAUS 
   /* --- Guardar MSD y DX para cada tau y cada f --- */
    FILE *fout_D = fopen("D_vs_tau_new.txt","a");
    for(int j=0;j<ntaus;j++){
        double tau = taus[j];
        for(int fi=0; fi<Nf_D; fi++){
            double f = forces_D[fi];
            printf("Calculando D para tau=%.2f, f=%g...\n", tau, f); // <-- PRINT DE PROGRESO
            double *MSD_tot = malloc(sizeof(double)*(NUM_STEPS+1));
            for(int i=0;i<=NUM_STEPS;i++) MSD_tot[i]=0.0;

            // Promediar NSIM simulaciones
            for(int sim=0; sim<NSIM; sim++){
                double *eps = malloc(sizeof(double)*N_particles);
                for(int p=0;p<N_particles;p++) eps[p] = (Random()<0.5)? -1.0 : 1.0;
                double *DX_dummy = malloc(sizeof(double)*(NUM_STEPS+1));
                for(int i=0;i<=NUM_STEPS;i++) DX_dummy[i]=0.0;

                Brownian_MSD_force(tau, flag_ini, f, eps, MSD_tot, DX_dummy);

                free(eps); free(DX_dummy);
            }
            for(int i=0;i<=NUM_STEPS;i++) MSD_tot[i] /= (double)NSIM;

            double t_total = NUM_STEPS*delta_t;
            double D_val = MSD_tot[NUM_STEPS]/(2.0*t_total);
            fprintf(fout_D,"%g\t%g\t%g\n", tau, f, D_val);

            free(MSD_tot);
        }
    }
    fclose(fout_D);

    /*

    FILE *fout_mu = fopen("mu_vs_tau.txt","a");
    if (!fout_mu) {
        perror("No se pudo abrir mu_vs_tau.txt");
        return 1;
    }


    double *DX_tot = malloc(sizeof(double)*(NUM_STEPS+1));
    double *MSD_dummy = malloc(sizeof(double)*(NUM_STEPS+1));
    double *eps = malloc(sizeof(double)*N_particles);
    if (!DX_tot || !MSD_dummy || !eps) {
        fprintf(stderr,"Fallo malloc\n");
        return 1;
    }


    for (int j = 0; j < ntaus; j++) {
        double tau = taus[j];
        for (int fi = 0; fi < Nf_mu; fi++) {
            double f = forces_mu[fi];

            for (int i = 0; i <= NUM_STEPS; i++) {
                DX_tot[i] = 0.0;
            }

            for (int sim = 0; sim < NSIM; sim++) {

                for (int p = 0; p < N_particles; p++) {
                    eps[p] = (Random() < 0.5f) ? -1.0 : 1.0;
                }

                for (int i = 0; i <= NUM_STEPS; i++) MSD_dummy[i] = 0.0;

                Brownian_MSD_force(tau, flag_ini, f, eps, MSD_dummy, DX_tot);
            }


            for (int i = 0; i <= NUM_STEPS; i++) DX_tot[i] /= (double)NSIM;


            double t_total = NUM_STEPS * delta_t;
            double mu_val = DX_tot[NUM_STEPS] / (f * t_total);

            fprintf(fout_mu, "%g\t%g\t%g\n", tau, f, mu_val);
            printf("tau=%.2f f=%g mu=%g\n", tau, f, mu_val);
        }
    }

    free(DX_tot);
    free(MSD_dummy);
    free(eps);
    fclose(fout_mu);
    */
    return 0;
}


/* Generacion de numeros aleatorios */

float Random(void){ //genera un numero aleatorio con distribucion plana en [0,1)
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;
    //printf("r=%f\n",r);
    return r;
}

/* Semilla del generador  */

void ini_ran(int SEMILLA){
    int INI,FACTOR,SUM,i;
    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;
    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
}

/* Algoritmo de Box Muller para la generación de dos números aleatorios con distribución normal 
    a partir de dos distribuidos uniformemente */

void BoxMuller(double* num1, double* num2){
    double parte1, parte2;
    parte1=sqrt(-2*log(Random()));
    parte2=2*pi*Random();

    *num1=parte1*cos(parte2);
    *num2=parte1*sin(parte2);
}

/* Generación de la distribucion */ 
void Normal_distribution(const char *filename){
    FILE *f; 
    int i;
    f=fopen(filename,"w");

    /** RANDOM NUMBERS GENERATION **/
    for(i=0;i<N_data;i++){
        double n1,n2;
        BoxMuller(&n1,&n2);
        fprintf(f,"%lf\n%lf\n",n1,n2);
    }
    fclose(f);
    printf("guardado!"); 
}

/* Periodic boundary conditions */
double PBC(double pos){
    while (pos > L/ 2.0) {
        pos -= L;
    }
    while (pos < -L / 2.0) {
        pos += L;
    }
    return pos;
}

// Simulación de N partículas con movimiento browniano para plotear trayectorias 
void Brownian_N(double tau, int flag_ini, const char *filename) {
    char filename_local[NMAX];
    sprintf(filename_local, "xi_tau_%.1f.txt", tau); // Archivo para guardar xi
    FILE *g = fopen(filename, "w");
    FILE *c = fopen(filename_local, "a");
    if (!g) {
        printf("Error al abrir archivo %s\n", filename);
        exit(1);
    }

    double x[N_particles], xi[N_particles];
    double xi_sum;
    double dx, extra;
    double t = 0.0;

    // Inicialización
    for (int p = 0; p < N_particles; p++) {
        // NUEVA INICIALIZACIÓN
        BoxMuller(&dx, &extra); // dx ~ N(0,1)
        xi[p] = sqrt(D / tau) * dx;
        if (flag_ini == 0)
            x[p] = 0.0;
        else
            x[p] = (Random() - 0.5) * L;
    }

    // Guardar posiciones iniciales
    for (int p = 0; p < N_particles; p++) fprintf(g, "%lf\t", x[p]);
    fprintf(g, "\n");

    // Bucle temporal
    for (int i = 0; i < NUM_STEPS; i++) {
        t += delta_t;
        xi_sum=0.0; //reseteo 
        for (int p = 0; p < N_particles; p++) {
            BoxMuller(&dx, &extra);
            xi[p] = xi[p] * (1.0 - delta_t / tau)+ sqrt(2.0 * D * delta_t / (tau * tau)) * dx;
            x[p] += delta_t / gamma * xi[p];
            //PBD
            x[p] = PBC(x[p]);
            xi_sum+=xi[p];
            fprintf(g, "%lf\t", x[p]);
            fprintf(c,"%lf\t",xi[p]); 
        }
        fprintf(g, "\n");
        fprintf(c,"\n");
            }
    fclose(g);
    fclose(c);
}


/* Simulation of N particles with brownian motion */
void Brownian_MSD(double tau,  int flag_ini, double *MSD_total){

    double x[N_particles],xo[N_particles];

    double dx,extra,t,msd;
    double xi[N_particles]; // Variables para el ruido gaussiano
    int i,p; 
    t=0.0;

    for(p=0;p<N_particles;p++){
        xi[p]=0;
    }

    // Damos una posición inicial a CADA partícula
    for (p = 0; p < N_particles; p++) {
        if (flag_ini == 0) {
            x[p] = 0.0;
        } else {
            // Posición aleatoria dentro de la caja [-L/2, L/2]
            x[p] = (Random() - 0.5) * L;
        }
        // Guardar la posición inicial
        xo[p]=x[p];
    }

    // --- BUCLE DE TIEMPO (Principal) ---
    for (i = 0; i <= NUM_STEPS; i++) {
        msd = 0.0; //para cada t hay uno
        // --- Escritura de datos ---
        // Escribimos un "bloque" de datos para este paso de tiempo
        //fprintf(g, "# t = %lf\n", t);

        // CALCULO MSD(t)
        for (p = 0; p < N_particles; p++) {
            // Formato: ID_Particula   X   Y
            double dxr = x[p]-xo[p];

            // CONDICIONES DE CONTORNO PARA EL DESPLAZAMIENTO
            //dxr -=L*round(dxr/L);
            //dyr -=L*round(dyr/L);

            msd+= dxr*dxr;
        }
        msd/=N_particles;
        MSD_total[i] +=msd; // Acumulamos el valor del MSD para este tiempo

        // Salir del bucle DESPUÉS de escribir el último paso
        if (i == NUM_STEPS) break; 

        // MOVIMIENTO DE LAS PARTICULAS
        t += delta_t; // Incrementamos el tiempo para el *próximo* paso
        
        for (p = 0; p < N_particles; p++) {

            BoxMuller(&dx, &extra); // Generamos dos números normales (el segundo no se usa)

            // 2. Mover la partícula (Euler-Maruyama con Ornstein-Uhlenbeck process)
            xi[p] = xi[p]* (1.0 - delta_t / tau) + sqrt(2.0 * D * delta_t /(tau*tau)) * dx; // new xi
            x[p] += delta_t/gamma*xi[p];  
            
        }
    }

    //printf("¡Simulacion de N partículas completada! La trayectoria se ha guardado\n");
}

/* Simulación de N partículas (no interactuantes) con ruido OU y fuerza f (en x).
   Acumula en MSD_total[i] y DX_total[i] (i=0..NUM_STEPS).
   eps[] debe ser de tamaño N_particles y contener ±1 (double).
*/
void Brownian_MSD_force(double tau, int flag_ini, double f, double *eps, double *MSD_total, double *DX_total) {

    double x[N_particles], xo[N_particles];
    double xu[N_particles], xuo[N_particles]; //sin PBC
    double xi[N_particles];
    double dx, extra,dxp, delta;
    int p, i;
    double t = 0.0;

    /* Inicialización variables */
    for (p = 0; p < N_particles; p++) {
        BoxMuller(&dx, &extra); // dx ~ N(0,1)
        xi[p] = sqrt(D / tau) * dx;
        if (flag_ini == 0) x[p] = 0.0;
        else x[p] = (Random() - 0.5) * L;
        xo[p] = x[p]; // referencia inicial
        xu[p] = x[p];
        xuo[p] = x[p];
    }

    /* Paso t=0: computar y acumular */
    double msd0 = 0.0;
    double dxavg0 = 0.0;
    for (p = 0; p < N_particles; p++) {
        double dxp = x[p] - xo[p]; // =0
        msd0 += dxp * dxp;
        dxavg0 += eps[p] * dxp;
    }
    msd0 /= N_particles; dxavg0 /= N_particles;
    MSD_total[0] += msd0;
    DX_total[0] += dxavg0;

    /* Bucle temporal */
    for (i = 1; i <= NUM_STEPS; i++) {
        // avanzar tiempo
        t += delta_t;

        // Actualizar partículas
        for (p = 0; p < N_particles; p++) {
            BoxMuller(&dx, &extra); // dx ~ N(0,1)
            xi[p] = xi[p] * (1.0 - delta_t / tau) + sqrt(2.0 * D * delta_t / (tau * tau)) * dx;
            delta = delta_t/gamma *(xi[p]+eps[p]*f);
            // Añadimos la fuerza f en x, con signo eps[p]
            xu[p] += delta;

            // Para la visualizacion
            x[p] += delta; 
            x[p] = PBC(x[p]);

        }

        // Calcular MSD_x y DX para este paso y acumular
        double msd = 0.0;
        double dxavg = 0.0;
        for (p = 0; p < N_particles; p++) {
            dxp = xu[p] - xuo[p];
            // Ajuste PBC para el desplazamiento (opcional en no envuelto; PBC aquí es suficiente porque xo en [-L/2,L/2])
            // dxp = dxp - L*round(dxp/L);  // si usas round() actívalo con math.h
            msd += dxp * dxp;
            dxavg += eps[p] * dxp;
        }
        msd /= N_particles;
        dxavg /= N_particles;

        MSD_total[i] += msd;
        DX_total[i] += dxavg;
    }
}
