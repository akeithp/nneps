/**
 * Magister en Informatica
 * Universidad Austral de Chile
 * Computacion de Alto Rendimiento
 * Prof. Dr. Héctor Ferrada
 * 
 * Alan Keith Paz
 * Due Date: 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <float.h>
#include <omp.h>
#include <random>

using namespace std;

#define PRINT 0
#define TEST 0
//#define NTHREADS 4

uint NTHREADS;	// number of threads
ulong REPET = 10;

float xup1, xup2, xle2, xle3, xlo3, xlo4, xri4, xri1;
float yup1, yup2, yle2, yle3, ylo3, ylo4, yri4, yri1;
float xc1, xc2, xc3, xc4;
float yc1, yc2, yc3, yc4;

ulong ri1, up1; 			
ulong up2, le2; 			
ulong le3, lo3; 			
ulong lo4, ri4;
ulong c1, c2, c3, c4;


typedef struct{
	float *X, *Y;	// array of float points X and Y
	bool NORMAL;	// flag probability: 1 = NORMAL, 0 = UNIFORM
	float sigma;	// for normal distribution probability function
	float mean;		// for normal distribution probability function
	ulong n;
	float minx;		// for uniform distribution
	float maxx;		// for uniform distribution 
	float miny;		// for uniform distribution 
	float maxy;		// for uniform distribution 
	float r;		// for circunsference (radious)
	float p;		// for circunsference (probability)
	//point *Q1;	// cuadrante 1
	//point *Q2;	// cuadrante 2
	//point *Q3;	// cuadrante 3
	//point *Q4;	// cuadrante 4*/

} pointSet;

struct Compare {float val; ulong index; };
#pragma omp declare reduction(minim : struct Compare : omp_out = omp_in.val < omp_out.val ? omp_in : omp_out) 
#pragma omp declare reduction(maxim : struct Compare : omp_out = omp_in.val > omp_out.val ? omp_in : omp_out) 

void genArrays(pointSet *points);
void genArraysCir(pointSet *points);
void divideArrays(pointSet *points);
void runEPS(pointSet *points);
void runEPSPDAC(pointSet *points);
void runEPSPReduce(pointSet *points);
void filteredPercent(pointSet *points);
void testPoints(pointSet *points); //comprobar correctitud de los puntos encontrados (secuencial vs paralelo)

int main(int argc, char const *argv[]){
	double t1,t2;
	float avgTime;
	char aFile[400];
	char str[100];
	uint i;
	bool flagP;		// flag para ejecución paralela
	bool flagDR;	// flag para ejecución Divide and Conquer o Reduce
	pointSet *points = new pointSet();
	
	flagP = atoi(argv[1]);
	points->n = atoi(argv[2]);
	points->NORMAL = atoi(argv[3]);
	if(points->NORMAL){
		points->mean = atof(argv[4]);
		points->sigma = atof(argv[5]);
	}
	else{
		points->minx = atof(argv[4]);
		points->maxx = atof(argv[5]);
		points->miny = atof(argv[6]);
		points->maxy = atof(argv[7]);
	}
	/*
	 * Parallel Section
	 * @argv[6] number of threads
	 */
	if(flagP){
		if(points->NORMAL) 
			NTHREADS = atoi(argv[6]);
		else 
			NTHREADS = atoi(argv[8]);
		omp_set_num_threads(NTHREADS);
		flagDR = atoi(argv[9]);
	}
	
	// OJO: Corregir comprobación de parámetros de entrada y mensaje de error
	/*if(flagP && points->NORMAL && argc < 7){
		cout << "Execution Error! call: ./prog <PARALLEL> <n> <NORMAL flag> [<mean>] [<sigma>] <threads> <REPETS FOR TEST" << endl;
		exit(EXIT_FAILURE);
	}*/
	if(PRINT){
		cout << "Parameters..." << endl;
		cout << "n = " << points->n << endl;
		cout << "NORMAL flag = " << points->NORMAL << endl;
		cout << "Number of Threads = " << NTHREADS << endl;
		if (points->NORMAL){
			cout << "mean = " << points->mean << endl;
			cout << "sigma = " << points->sigma << endl;
		}
	}
	//float t; // = clock();
	
	genArrays(points);
	
	//t = (clock() - t)/CLOCKS_PER_SEC;	// seconds
	//cout << "Construction Time = " << t << " Seconds" << endl;
	
	/*
	 * Serial version
	 */
	 if(!flagP){
		//t = clock();
		t1=omp_get_wtime(); 
		for(i=0; i<REPET; i++){
			runEPS(points);
		}
		t2=omp_get_wtime();
		avgTime = (t2 - t1)/REPET;
		//t = (clock() - t)/CLOCKS_PER_SEC;	// seconds
		//cout << "Serial Time = " << t/REPET << " Seconds" << endl;
		if(PRINT){
			cout << "Average CPU time per execution, Sequential = " << avgTime*1000000.0 << " Microseconds" << endl; //revisar unidad tiempos
		}
		strcpy(aFile, "./RESULTS/");
		strcpy(str, "");
		sprintf(str, "EPS%i", points->NORMAL); //%ld, points->n
		strcat(aFile, str);
		cout << "Resume File: " << aFile << endl;
		FILE *fp = fopen(aFile, "a+" );
		if (points->NORMAL){
			// [n] [avg eps-time/exec] [esperanza] [varianza]
			fprintf(fp, "%ld %f %f %f\n", points->n, (avgTime*1000000.0), points->mean, points->sigma);
		}else{
			// [n] [nOcc] [avg bs-time/exec]
			fprintf(fp, "%ld %f\n", points->n, (avgTime*1000000.0));
		}
		fclose(fp);
	}
	
	/*
	 * Parallel version
	 */
	 if(flagP){
		 if(flagDR){
			t1=omp_get_wtime(); 
			for(i=0; i<REPET; i++){
				runEPSPDAC(points);
			}
			t2=omp_get_wtime();
			avgTime = (t2 - t1)/REPET;
			if(PRINT){
				cout << "Average CPU time per execution, Parallel DAC= " << avgTime*1000000.0 << " Microseconds" << endl;
			}
			strcpy(aFile, "./RESULTSP/");
			strcpy(str, "");
			sprintf(str, "EPSParallelD%i%d", points->NORMAL, NTHREADS); //%ld, points->n
			strcat(aFile, str);
			cout << "Resume File: " << aFile << endl;
			FILE *fp = fopen(aFile, "a+" );
			if (points->NORMAL){
				// [n] [avg epsP-time/exec] [esperanza] [varianza]
				fprintf(fp, "%ld %f %f %f\n", points->n, (avgTime*1000000.0), points->mean, points->sigma);
			}else{
				// [n] [avg bs-time/exec]
				fprintf(fp, "%ld %f\n", points->n, (avgTime*1000000.0));
			}
			fclose(fp);
		}
		else{
			t1=omp_get_wtime(); 
			for(i=0; i<REPET; i++){
				runEPSPReduce(points);
			}
			t2=omp_get_wtime();
			avgTime = (t2 - t1)/REPET;
			if(PRINT){
				cout << "Average CPU time per execution, Parallel Reduction= " << avgTime*1000000.0 << " Microseconds" << endl;
			}
			strcpy(aFile, "./RESULTSP/");
			strcpy(str, "");
			sprintf(str, "EPSParallelR%i%d", points->NORMAL, NTHREADS); //%ld, points->n
			strcat(aFile, str);
			cout << "Resume File: " << aFile << endl;
			FILE *fp = fopen(aFile, "a+" );
			if (points->NORMAL){
				// [n] [avg epsP-time/exec] [esperanza] [varianza]
				fprintf(fp, "%ld %f %f %f\n", points->n, (avgTime*1000000.0), points->mean, points->sigma);
			}else{
				// [n] [avg bs-time/exec]
				fprintf(fp, "%ld %f\n", points->n, (avgTime*1000000.0));
			}
			fclose(fp);
			if(TEST){
				testPoints(points);
			}
		}
		if(TEST){
			testPoints(points);
		}
	}
	filteredPercent(points);
	return 0;
}

// generar arrays
void genArrays(pointSet *points){
	ulong i;
	float numX, numY;

	if (points->NORMAL){ // DISTRIBUCION NORMAL
		points->X = new float[points->n]; // +1?
		points->Y = new float[points->n]; // +1?
		default_random_engine generator;
		normal_distribution<double> distribution(points->mean, points->sigma);	// (mean, stddev)

		for (i=0; i<points->n; i++){
			numX = distribution(generator);
	    	numY = distribution(generator);
	    	points->X[i] = numX;
	    	points->Y[i] = numY;
		}
	}else{ // DISTRIBUCION UNIFORME
		points->X = new float[points->n]; // +1?
		points->Y = new float[points->n]; // +1?
		float varx = points->maxx - points->minx;
		float vary = points->maxy - points->miny;
		for (i=0; i<points->n; i++){
			points->X[i] = (float)rand()/RAND_MAX*varx+points->minx;
			points->Y[i] = (float)rand()/RAND_MAX*vary+points->miny;
		}
	}

}

/**
 * Generar Arreglo para el caso de circunsferencia
 */
void genArraysCir(pointSet *points){
	
	// PENDIENTE
}

/**
 * Divide los arrays por cuadrante
 */
void divideArrays(pointSet *points){
	//ulong i, q1, q2, q3, q4;
	//q1 = q2 = q3 = q4 = 0;
	/*for(i=0; i<points->n;i++){
		if( (points->X[i].x >= 0) && (points->X[i].y >= 0) ){
			points->Q1[q1] = points[i];
			q1++;
		}
		if( (points->X[i].x >= 0) && (points->X[i].y < 0) ){
			points->Q2[q2] = points[i];
			q2++;
		}
		if( (points->X[i].x < 0) && (points->X[i].y < 0) ){
			points->Q3[q3] = points[i];
			q3++;
		}
		if( (points->X[i].x < 0) && (points->X[i].y >= 0) ){
			points->Q4[q4] = points[i];
			q4++;
		}			
	}*/
}

/**
 *  Busca los puntos extremos de un array (Secuencial)
 */
void runEPS(pointSet *points){
	float xi,yi;
	bool seei=true;

	c1 = c2 = c3 = c4 = points->n+1;
	xc1 = xc4 = yc1 = yc2 = -1*FLT_MAX;
	xc2 = xc3 = yc3 = yc4 = FLT_MAX;
	xri1 = xup1 = xup2 = xle2 = xle3 = xlo3 = xlo4 = xri4 = points->X[0];
	yri1 = yup1 = yup2 = yle2 = yle3 = ylo3 = ylo4 = yri4 = points->Y[0];
	ri1 = up1 = up2 = le2 = le3 = lo3 = lo4 = ri4 = 0; 
	for(ulong i=1; i<points->n; i++, seei=true){
		xi=points->X[i];yi=points->Y[i];
		if(xi < xle2){
			if (xc2-yc2 > xle2-yle2){
				c2 = le2;
				xc2=xle2;
				yc2=yle2;
			}
			if (xc3+yc3 > xle3+yle3){
				c3 = le3;
				xc3=xle3;
				yc3=yle3;
			}
			xle2 = xle3 = xi;
			yle2 = yle3 = yi;
			le2 = le3 = i;
			seei=false;
		}else{
			if(xi == xle2){
				if(yi > yle2){
					le2 = i;
					xle2 = xi;
					yle2 = yi;
					seei=false;
				}else{
					if(yi < yle3){
						le3 = i;
						xle3 = xi;
						yle3 = yi;
						seei=false;
					}
				}
			}else{
				if(xi > xri1){
					if (xc1+yc1 < xri1+yri1){
						c1 = ri1;
						xc1=xri1;
						yc1=yri1;
					}
					if (yc4-xc4 > yri4-xri4){
						c4 = ri4;
						xc4=xri4;
						yc4=yri4;
					}
					ri1 = ri4 = i;
					xri1 = xri4 = xi;
					yri1 = yri4 = yi;
					seei=false;
				}else{
					if(xi == xri1){
						if(yi > yri1){
							ri1 = i;
							xri1 = xi;
							yri1 = yi;
							seei=false;
						}else{
							if(yi < yri4){
								ri4 = i;
								xri4 = xi;
								yri4 = yi;
								seei=false;
							}
						}
					}
				}
			}
		}


		if(yi < ylo3){
			if (xc3+yc3 > xlo3+ylo3){
				xc3=xlo3;
				yc3=ylo3;
				c3 = lo3;
			}
			if (yc4-xc4 > ylo4-xlo4){
				xc4=xlo4;
				yc4=ylo4;
				c4 = lo4;
			}

			xlo3 = xlo4 = xi;
			ylo3 = ylo4 = yi;
			lo3 = lo4 = i;
			seei=false;
		}else{
			if(yi == ylo3){
				if(xi < xlo3){
					xlo3 = xi;
					ylo3 = yi;
					lo3 = i;
					seei=false;
				}else{
					if(xi > xlo4){
						xlo4 = xi;
						ylo4 = yi;
						lo4 = i;
						seei=false;
					}
				}
			}else{
				if(yi > yup2){
					if (xc1+yc1 < xup1+yup1){
						xc1=xup1;
						yc1=yup1;
						c1 = up1;
					}
					if (xc2-yc2 > xup2-yup2){
						xc2=xup2;
						yc2=yup2;
						c2 = up2;
					}

					xup2 = xup1 = xi;
					yup2 = yup1 = yi;
					up1 = up2 = i;
					seei=false;
				}else{
					if(yi == yup2){
						if(xi < xup2){
							xup2 = xi;
							yup2 = yi;
							up2 = i;
							seei=false;
						}else{
							if(xi > xup1){
								xup1 = xi;
								yup1 = yi;
								up1 = i;
								seei=false;
							}
						}
					}
				}
			}

		}

		if (seei){
			if (xc1+yc1 < xi+yi){
				c1=i;
				xc1=xi;
				yc1=yi;
			}else{
				if (xc2-yc2 > xi-yi){
					c2=i;
					xc2=xi;
					yc2=yi;
				}else{
					if (xc3+yc3 > xi+yi){
						c3=i;
						xc3=xi;
						yc3=yi;
					}else{
						if (yc4-xc4 > yi-xi){
							c4=i;
							xc4=xi;
							yc4=yi;
						}
					}
				}
			}
		}
	}
}

/** 
 * Búsqueda de puntos extremos en paralelo Divide and Conquer
 */
void runEPSPDAC(pointSet *points){
	float xup1p[NTHREADS], xup2p[NTHREADS], xle2p[NTHREADS], xle3p[NTHREADS], xlo3p[NTHREADS], xlo4p[NTHREADS], xri4p[NTHREADS], xri1p[NTHREADS];
	float yup1p[NTHREADS], yup2p[NTHREADS], yle2p[NTHREADS], yle3p[NTHREADS], ylo3p[NTHREADS], ylo4p[NTHREADS], yri4p[NTHREADS], yri1p[NTHREADS];
	float xc1p[NTHREADS], xc2p[NTHREADS], xc3p[NTHREADS], xc4p[NTHREADS];
	float yc1p[NTHREADS], yc2p[NTHREADS], yc3p[NTHREADS], yc4p[NTHREADS];
	ulong ri1p[NTHREADS], up1p[NTHREADS];
	ulong up2p[NTHREADS], le2p[NTHREADS];
	ulong le3p[NTHREADS], lo3p[NTHREADS];
	ulong lo4p[NTHREADS], ri4p[NTHREADS];
	ulong c1p[NTHREADS], c2p[NTHREADS], c3p[NTHREADS], c4p[NTHREADS];

	#pragma omp parallel
	{
		uint id;
		uint step, start, stop;
		ulong i, j;
		float xi,yi;
		step = (int)points->n/omp_get_num_threads();
		id = omp_get_thread_num();
		start = id * step;
		if (id != (NTHREADS-1))
			stop = start + step;
		else
			stop = points->n;

		xri1p[id] = xup1p[id] = xup2p[id] = xle2p[id] = xle3p[id] = xlo3p[id] = xlo4p[id] = xri4p[id] = points->X[start];
		yri1p[id] = yup1p[id] = yup2p[id] = yle2p[id] = yle3p[id] = ylo3p[id] = ylo4p[id] = yri4p[id] = points->Y[start];
		ri1p[id] = up1p[id] = up2p[id] = le2p[id] = le3p[id] = lo3p[id] = lo4p[id] = ri4p[id] = 0;
		#pragma omp parallel for private(i, xi, yi) shared(xri1p,xup1p,xup2p,xle2p,xle3p,xlo3p,xlo4p,xri4p,yri1p,yup1p,yup2p,yle2p,yle3p,ylo3p,ylo4p,yri4p,ri1p,up1p,up2p,le2p,le3p,lo3p,lo4p,ri4p)
		for(i=start+1; i<stop; i++){
			xi=points->X[i];yi=points->Y[i];
			if(xi < xle2p[id]){
				xle2p[id] = xle3p[id] = xi;
				yle2p[id] = yle3p[id] = yi;
				le2p[id] = le3p[id] = i;
			}else{
				if(xi == xle2p[id]){
					if(yi > yle2p[id]){
						le2p[id] = i;
						xle2p[id] = xi;
						yle2p[id] = yi;
					}else{
						if(yi < yle3p[id]){
							le3p[id] = i;
							xle3p[id] = xi;
							yle3p[id] = yi;
						}
					}
				}else{
					if(xi > xri1p[id]){
						ri1p[id] = ri4p[id] = i;
						xri1p[id] = xri4p[id] = xi;
						yri1p[id] = yri4p[id] = yi;
					}else{
						if(xi == xri1p[id]){
							if(yi > yri1p[id]){
								ri1p[id] = i;
								xri1p[id] = xi;
								yri1p[id] = yi;
							}else{
								if(yi < yri4p[id]){
									ri4p[id] = i;
									xri4p[id] = xi;
									yri4p[id] = yi;
								}
							}
						}
					}
				}
			}


			if(yi < ylo3p[id]){
				xlo3p[id] = xlo4p[id] = xi;
				ylo3p[id] = ylo4p[id] = yi;
				lo3p[id] = lo4p[id] = i;
			}else{
				if(yi == ylo3p[id]){
					if(xi < xlo3p[id]){
						xlo3p[id] = xi;
						ylo3p[id] = yi;
						lo3p[id] = i;
					}else{
						if(xi > xlo4p[id]){
							xlo4p[id] = xi;
							ylo4p[id] = yi;
							lo4p[id] = i;
						}
					}
				}else{
					if(yi > yup2p[id]){
						xup2p[id] = xup1p[id] = xi;
						yup2p[id] = yup1p[id] = yi;
						up1p[id] = up2p[id] = i;
					}else{
						if(yi == yup2p[id]){
							if(xi < xup2p[id]){
								xup2p[id] = xi;
								yup2p[id] = yi;
								up2p[id] = i;
							}else{
								if(xi > xup1p[id]){
									xup1p[id] = xi;
									yup1p[id] = yi;
									up1p[id] = i;
								}
							}
						}
					}
				}

			}
		}
		#pragma omp barrier
		
		xri1=xri1p[0]; xup1=xup1p[0]; xup2=xup2p[0]; xle2=xle2p[0]; xle3=xle3p[0]; xlo3=xlo3p[0]; xlo4=xlo4p[0]; xri4=xri4p[0];
		yri1=yri1p[0]; yup1=yup1p[0]; yup2=yup2p[0]; yle2=yle2p[0]; yle3=yle3p[0]; ylo3=ylo3p[0]; ylo4=ylo4p[0]; yri4=yri4p[0];
		ri1=ri1p[0]; up1=up1p[0]; up2=up2p[0]; le2=le2p[0]; le3=le3p[0]; lo3=lo3p[0]; lo4=lo4p[0]; ri4=ri4p[0];
		#pragma omp critical
		{
			for(uint i=1; i<NTHREADS; i++)
			{
				if(xle2p[i] < xle2){
					xle2 = xle2p[i];
					yle2 = yle2p[i];
					le2 = le2p[i];
				}
				if(xri1p[i] > xri1){
					xri1 = xri1p[i];
					yri1 = yri1p[i];
					ri1 = ri1p[i];
				}
				if(yup1p[i] > yup1){
					yup1 = yup1p[i];
					xup1 = xup1p[i];
					up1 = up1p[i];
				}
				if(yup2p[i] > yup2){
					yup2 = yup2p[i];
					xup2 = xup2p[i];
					up2 = up2p[i];
				}
			
				if(xle3p[i] < xle3){
					xle3 = xle3p[i];
					yle3 = yle3p[i];
					le3 = le3p[i];
				}
				if(ylo3p[i] < ylo3){
					ylo3 = ylo3p[i];
					xlo3 = xlo3p[i];
					lo3 = lo3p[i];
				}
				if(ylo4p[i] < ylo4){
					ylo4 = ylo4p[i];
					xlo4 = xlo4p[i];
					lo4 = lo4p[i];
				}
				if(xri4p[i] > xri4){
					xri4 = xri4p[i];
					yri4 = yri4p[i];
					ri4 = ri4p[i]; 
				}
			}
		}
		#pragma omp barrier

		xc1p[id] = xc4p[id] = yc1p[id] = yc2p[id] = -1*FLT_MAX;
		xc2p[id] = xc3p[id] = yc3p[id] = yc4p[id] = FLT_MAX;
		#pragma omp parallel for private(j, xi, yi) shared(c1p, c2p, c3p, c4p)
		for(j=start+1; j<stop; j++){
			xi=points->X[j];yi=points->Y[j];
			if(xi>xup1 && yi>yri1){
				if (xc1p[id]+yc1p[id] < xi+yi){
					c1p[id]=j;
					xc1p[id]=xi;
					yc1p[id]=yi;
				}
			}
			if(xi<xup2 && yi>yle2){
				if (xc2p[id]-yc2p[id] > xi-yi){
					c2p[id]=j;
					xc2p[id]=xi;
					yc2p[id]=yi;
				}
			}
			if(xi<xlo3 && yi<yle3){
				if (xc3p[id]+yc3p[id] > xi+yi){
					c3p[id]=j;
					xc3p[id]=xi;
					yc3p[id]=yi;
				}
			}
			if(xi>xlo4 && yi <yri4){
				if (yc4p[id]-xc4p[id] > yi-xi){
					c4p[id]=j;
					xc4p[id]=xi;
					yc4p[id]=yi;
				}
			}
		}

		#pragma omp critical
		{
			xc1=xc1p[0]; xc2=xc2p[0]; xc3=xc3p[0]; xc4=xc4p[0];
			yc1=yc1p[0]; yc2=yc2p[0]; yc3=yc3p[0]; yc4=yc4p[0];
			for(uint k=1; k<NTHREADS; k++){
				if(xc1+yc1 < xc1p[k]+yc1p[k]){
					c1=c1p[k];
					xc1=xc1p[k];
					yc1=yc1p[k];
				}
				if(xc2-yc2 > xc2p[k]-yc2p[k]){
					c2=c2p[k];
					xc2=xc2p[k];
					yc2=yc2p[k];
				}
				if(xc3+yc3 > xc3p[k]+yc3p[k]){
					c3=c3p[k];
					xc3=xc3p[k];
					yc3=yc3p[k];
				}
				if(yc4-xc4 > yc4p[k]-xc4p[k]){
					c4=c4p[k];
					xc4=xc4p[k];
					yc4=yc4p[k];
				}
			}
		}
	}

}

/** 
 * Búsqueda de puntos extremos en paralelo Divide and Conquer
 */
void runEPSPDAC2(pointSet *points){

	#pragma omp parallel
	{
		uint id;
		uint step, start, stop;
		ulong i, j;
		float xi,yi;
		bool seei;
		float xup1p[NTHREADS], xup2p[NTHREADS], xle2p[NTHREADS], xle3p[NTHREADS], xlo3p[NTHREADS], xlo4p[NTHREADS], xri4p[NTHREADS], xri1p[NTHREADS];
		float yup1p[NTHREADS], yup2p[NTHREADS], yle2p[NTHREADS], yle3p[NTHREADS], ylo3p[NTHREADS], ylo4p[NTHREADS], yri4p[NTHREADS], yri1p[NTHREADS];
		float xc1p[NTHREADS], xc2p[NTHREADS], xc3p[NTHREADS], xc4p[NTHREADS];
		float yc1p[NTHREADS], yc2p[NTHREADS], yc3p[NTHREADS], yc4p[NTHREADS];
		ulong ri1p[NTHREADS], up1p[NTHREADS];
		ulong up2p[NTHREADS], le2p[NTHREADS];
		ulong le3p[NTHREADS], lo3p[NTHREADS];
		ulong lo4p[NTHREADS], ri4p[NTHREADS];
		ulong c1p[NTHREADS], c2p[NTHREADS], c3p[NTHREADS], c4p[NTHREADS];
		
		step = (int)points->n/omp_get_num_threads();
		id = omp_get_thread_num();
		start = id * step;
		if (id != (NTHREADS-1))
			stop = start + step;
		else
			stop = points->n;

		xri1p[id] = xup1p[id] = xup2p[id] = xle2p[id] = xle3p[id] = xlo3p[id] = xlo4p[id] = xri4p[id] = points->X[start];
		yri1p[id] = yup1p[id] = yup2p[id] = yle2p[id] = yle3p[id] = ylo3p[id] = ylo4p[id] = yri4p[id] = points->Y[start];
		ri1p[id] = up1p[id] = up2p[id] = le2p[id] = le3p[id] = lo3p[id] = lo4p[id] = ri4p[id] = 0;
		#pragma omp parallel for private(i, xi, yi, seei) shared(xri1p,xup1p,xup2p,xle2p,xle3p,xlo3p,xlo4p,xri4p,yri1p,yup1p,yup2p,yle2p,yle3p,ylo3p,ylo4p,yri4p,ri1p,up1p,up2p,le2p,le3p,lo3p,lo4p,ri4p)
		for(i=start+1; i<stop; i++){
			seei=true;
			xi=points->X[i];yi=points->Y[i];
			if(xi < xle2p[id]){
				if (xc2p[id]-yc2p[id] > xle2p[id]-yle2p[id]){
					c2p[id] = le2p[id];
					xc2p[id]=xle2p[id];
					yc2p[id]=yle2p[id];
				}
				if (xc3p[id]+yc3p[id] > xle3p[id]+yle3p[id]){
					c3p[id] = le3p[id];
					xc3p[id]=xle3p[id];
					yc3p[id]=yle3p[id];
				}
				xle2p[id] = xle3p[id] = xi;
				yle2p[id] = yle3p[id] = yi;
				le2p[id] = le3p[id] = i;
				seei=false;
			}else{
				if(xi == xle2p[id]){ 
					if(yi > yle2p[id]){
						le2p[id] = i;
						xle2p[id] = xi;
						yle2p[id] = yi;
						seei=false;
					}else{
						if(yi < yle3p[id]){
							le3p[id] = i;
							xle3p[id] = xi;
							yle3p[id] = yi;
							seei=false;
						}
					}
				}else{
					if(xi > xri1p[id]){
						if (xc1p[id]+yc1p[id] < xri1p[id]+yri1p[id]){
							c1p[id] = ri1p[id];
							xc1p[id]=xri1p[id];
							yc1p[id]=yri1p[id];
						}
						if (yc4p[id]-xc4p[id] > yri4p[id]-xri4p[id]){
							c4p[id] = ri4p[id];
							xc4p[id]=xri4p[id];
							yc4p[id]=yri4p[id];
						}
						ri1p[id] = ri4p[id] = i;
						xri1p[id] = xri4p[id] = xi;
						yri1p[id] = yri4p[id] = yi;
						seei=false;
					}else{
						if(xi == xri1p[id]){
							if(yi > yri1p[id]){
								ri1p[id] = i;
								xri1p[id] = xi;
								yri1p[id] = yi;
								seei=false;
							}else{
								if(yi < yri4p[id]){
									ri4p[id] = i;
									xri4p[id] = xi;
									yri4p[id] = yi;
									seei=false;
								}
							}
						}
					}
				}
			}


			if(yi < ylo3p[id]){
				if (xc3p[id]+yc3p[id] > xlo3p[id]+ylo3p[id]){
					xc3p[id]=xlo3p[id];
					yc3p[id]=ylo3p[id];
					c3p[id] = lo3p[id];
				}
				if (yc4p[id]-xc4p[id] > ylo4p[id]-xlo4p[id]){
					xc4p[id]=xlo4p[id];
					yc4p[id]=ylo4p[id];
					c4p[id] = lo4p[id];
				}

				xlo3p[id] = xlo4p[id] = xi;
				ylo3p[id] = ylo4p[id] = yi;
				lo3p[id] = lo4p[id] = i;
				seei=false;
			}else{
				if(yi == ylo3p[id]){
					if(xi < xlo3p[id]){
						xlo3p[id] = xi;
						ylo3p[id] = yi;
						lo3p[id] = i;
						seei=false;
					}else{
						if(xi > xlo4p[id]){
							xlo4p[id] = xi;
							ylo4p[id] = yi;
							lo4p[id] = i;
							seei=false;
						}
					}
				}else{
					if(yi > yup2p[id]){
						if (xc1p[id]+yc1p[id] < xup1p[id]+yup1p[id]){
							xc1p[id]=xup1p[id];
							yc1p[id]=yup1p[id];
							c1p[id] = up1p[id];
						}
						if (xc2p[id]-yc2p[id] > xup2p[id]-yup2p[id]){
							xc2p[id]=xup2p[id];
							yc2p[id]=yup2p[id];
							c2p[id] = up2p[id];
						}

						xup2p[id] = xup1p[id] = xi;
						yup2p[id] = yup1p[id] = yi;
						up1p[id] = up2p[id] = i;
						seei=false;
					}else{
						if(yi == yup2p[id]){
							if(xi < xup2p[id]){
								xup2p[id] = xi;
								yup2p[id] = yi;
								up2p[id] = i;
								seei=false;
							}else{
								if(xi > xup1p[id]){
									xup1p[id] = xi;
									yup1p[id] = yi;
									up1p[id] = i;
									seei=false;
								}
							}
						}
					}
				}

			}
		}
		#pragma omp barrier
		
		xri1=xri1p[0]; xup1=xup1p[0]; xup2=xup2p[0]; xle2=xle2p[0]; xle3=xle3p[0]; xlo3=xlo3p[0]; xlo4=xlo4p[0]; xri4=xri4p[0];
		yri1=yri1p[0]; yup1=yup1p[0]; yup2=yup2p[0]; yle2=yle2p[0]; yle3=yle3p[0]; ylo3=ylo3p[0]; ylo4=ylo4p[0]; yri4=yri4p[0];
		ri1=ri1p[0]; up1=up1p[0]; up2=up2p[0]; le2=le2p[0]; le3=le3p[0]; lo3=lo3p[0]; lo4=lo4p[0]; ri4=ri4p[0];
		#pragma omp critical
		{
			for(uint i=1; i<NTHREADS; i++)
			{
				if(xle2p[i] < xle2){
					xle2 = xle2p[i];
					yle2 = yle2p[i];
					le2 = le2p[i];
				}
				if(xri1p[i] > xri1){
					xri1 = xri1p[i];
					yri1 = yri1p[i];
					ri1 = ri1p[i];
				}
				if(yup1p[i] > yup1){
					yup1 = yup1p[i];
					xup1 = xup1p[i];
					up1 = up1p[i];
				}
				if(yup2p[i] > yup2){
					yup2 = yup2p[i];
					xup2 = xup2p[i];
					up2 = up2p[i];
				}
			
				if(xle3p[i] < xle3){
					xle3 = xle3p[i];
					yle3 = yle3p[i];
					le3 = le3p[i];
				}
				if(ylo3p[i] < ylo3){
					ylo3 = ylo3p[i];
					xlo3 = xlo3p[i];
					lo3 = lo3p[i];
				}
				if(ylo4p[i] < ylo4){
					ylo4 = ylo4p[i];
					xlo4 = xlo4p[i];
					lo4 = lo4p[i];
				}
				if(xri4p[i] > xri4){
					xri4 = xri4p[i];
					yri4 = yri4p[i];
					ri4 = ri4p[i]; 
				}
			}
		}
		#pragma omp barrier

		#pragma omp critical
		{
			xc1=xc1p[0]; xc2=xc2p[0]; xc3=xc3p[0]; xc4=xc4p[0];
			yc1=yc1p[0]; yc2=yc2p[0]; yc3=yc3p[0]; yc4=yc4p[0];
			for(uint k=1; k<NTHREADS; k++){
				if(xc1+yc1 < xc1p[k]+yc1p[k]){
					c1=c1p[k];
					xc1=xc1p[k];
					yc1=yc1p[k];
				}
				if(xc2-yc2 > xc2p[k]-yc2p[k]){
					c2=c2p[k];
					xc2=xc2p[k];
					yc2=yc2p[k];
				}
				if(xc3+yc3 > xc3p[k]+yc3p[k]){
					c3=c3p[k];
					xc3=xc3p[k];
					yc3=yc3p[k];
				}
				if(yc4-xc4 > yc4p[k]-xc4p[k]){
					c4=c4p[k];
					xc4=xc4p[k];
					yc4=yc4p[k];
				}
			}
		}
	}

}

/**
 * Cálculo de Puntos extremos en paralelo Reduction
 */
void runEPSPReduce(pointSet *points){

	#pragma omp parallel 
	{
		ulong i;
		xri1 = xup1 = xup2 = xle2 = xle3 = xlo3 = xlo4 = xri4 = points->X[0];
		yri1 = yup1 = yup2 = yle2 = yle3 = ylo3 = ylo4 = yri4 = points->Y[0];
		xc1 = xc4 = yc1 = yc2 = -1*FLT_MAX;
		xc2 = xc3 = yc3 = yc4 = FLT_MAX;
		ri1 = up1 = up2 = le2 = le3 = lo3 = lo4 = ri4 = 0; 
		float dm1, dm2, dm3, dm4;
		dm1 = -1*FLT_MAX;
		dm2 = dm3 = dm4 = FLT_MAX;
		float xi,yi;
		bool seei=true;
		#pragma omp parallel for private(i, seei) reduction(max:xri1,yup1,yup2,xri4,dm1) reduction(min:xle2,xle3,ylo3,ylo4,dm2,dm3,dm4) //shared(c1,c2,c3,c4)
		for(i=1; i<points->n; i++){
		seei = true;
		xi=points->X[i];yi=points->Y[i];
		if(xi < xle2){
			if (xc2-yc2 > xle2-yle2){
				dm2 = xle2-yle2;
				c2 = le2;
				xc2=xle2;
				yc2=yle2;
			}
			if (xc3+yc3 > xle3+yle3){
				dm3=xle2+yle3;
				c3 = le3;
				xc3=xle3;
				yc3=yle3;
			}
			xle2 = xle3 = xi;
			yle2 = yle3 = yi;
			le2 = le3 = i;
			seei=false;
		}else{
			if(xi == xle2){
				if(yi > yle2){
					le2 = i;
					xle2 = xi;
					yle2 = yi;
					seei=false;
				}else{
					if(yi < yle3){
						le3 = i;
						xle3 = xi;
						yle3 = yi;
						seei=false;
					}
				}
			}else{
				if(xi > xri1){
					if (xc1+yc1 < xri1+yri1){
						dm1=xri1+yri1;
						c1 = ri1;
						xc1=xri1;
						yc1=yri1;
					}
					if (yc4-xc4 > yri4-xri4){
						dm4=yri4-xri4;
						c4 = ri4;
						xc4=xri4;
						yc4=yri4;
					}
					ri1 = ri4 = i;
					xri1 = xri4 = xi;
					yri1 = yri4 = yi;
					seei=false;
				}else{
					if(xi == xri1){
						if(yi > yri1){
							ri1 = i;
							xri1 = xi;
							yri1 = yi;
							seei=false;
						}else{
							if(yi < yri4){
								ri4 = i;
								xri4 = xi;
								yri4 = yi;
								seei=false;
							}
						}
					}
				}
			}
		}


		if(yi < ylo3){
			if (xc3+yc3 > xlo3+ylo3){
				dm3=xlo3+ylo3;
				xc3=xlo3;
				yc3=ylo3;
				c3 = lo3;
			}
			if (yc4-xc4 > ylo4-xlo4){
				dm4=ylo4-xlo4;
				xc4=xlo4;
				yc4=ylo4;
				c4 = lo4;
			}

			xlo3 = xlo4 = xi;
			ylo3 = ylo4 = yi;
			lo3 = lo4 = i;
			seei=false;
		}else{
			if(yi == ylo3){
				if(xi < xlo3){
					xlo3 = xi;
					ylo3 = yi;
					lo3 = i;
					seei=false;
				}else{
					if(xi > xlo4){
						xlo4 = xi;
						ylo4 = yi;
						lo4 = i;
						seei=false;
					}
				}
			}else{
				if(yi > yup2){
					if (xc1+yc1 < xup1+yup1){
						dm1=xup1+yup1;
						xc1=xup1;
						yc1=yup1;
						c1 = up1;
					}
					if (xc2-yc2 > xup2-yup2){
						dm2=xup2-yup2;
						xc2=xup2;
						yc2=yup2;
						c2 = up2;
					}

					xup2 = xup1 = xi;
					yup2 = yup1 = yi;
					up1 = up2 = i;
					seei=false;
				}else{
					if(yi == yup2){
						if(xi < xup2){
							xup2 = xi;
							yup2 = yi;
							up2 = i;
							seei=false;
						}else{
							if(xi > xup1){
								xup1 = xi;
								yup1 = yi;
								up1 = i;
								seei=false;
							}
						}
					}
				}
			}

		}

		if (seei){
			if (xc1+yc1 <= xi+yi){
				dm1 = xi+yi;
				c1=i;
				xc1=xi;
				yc1=yi;
			}else{
				if (xc2-yc2 >= xi-yi){
					dm2=xi-yi;
					c2=i;
					xc2=xi;
					yc2=yi;
				}else{
					if (xc3+yc3 >= xi+yi){
						dm3=xi+yi;
						c3=i;
						xc3=xi;
						yc3=yi;
					}else{
						if (yc4-xc4 >= yi-xi){
							dm4=yi-xi;
							c4=i;
							xc4=xi;
							yc4=yi;
						}
					}
				}
			}
		}
	}

	}
}


/** 
 * Cálculo de Porcentaje de Puntos Filtrados
 * Corregir
 */
 void filteredPercent(pointSet *points){
	/*ulong i;
	ulong count = 0;
	float percent = 0;
	char aFile[400];
	char str[100];
	float a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6,b7, b8;
	a1 = (yup-yup_le)/(xup-xup_le);
	b1 = -(yup-yup_le)/(xup-xup_le) + yup_le;
	a2 = (yleft-yup_le)/(xleft-xup_le);
	b2 = -(yleft-yup_le)/(xleft-xup_le) + yup_le;
	a3 = (yleft-yle_do)/(xleft-xle_do);
	b3 = -(yleft-yle_do)/(xleft-xle_do) + yle_do;
	a4 = (ydown-yle_do)/(xdown-xle_do);
	b4 = -(ydown-yle_do)/(xdown-xle_do) + yle_do;
	a5 = (ydown-ydo_ri)/(xdown-xdo_ri);
	b5 = -(ydown-ydo_ri)/(xdown-xdo_ri) + ydo_ri;
	a6 = (yright-ydo_ri)/(xright-xdo_ri);
	b6 = -(yright-ydo_ri)/(xright-xdo_ri) + ydo_ri;
	a7 = (yright-yri_up)/(xright-xri_up);
	b7 = -(yright-yri_up)/(xright-xri_up) + yri_up;
	a8 = (yup-yri_up)/(xup-xri_up);
	b8 = -(yup-yri_up)/(xup-xri_up) + yri_up;
	for(i=0; i<points->n; i++){
		if(points->X[i] < ((points->Y[i]-b1)/a1) && points->Y[i] > (a1*points->X[i]+b1)){
			count++;
		}
		else if(points->X[i] < ((points->Y[i]-b2)/a2) && points->Y[i] > (a2*points->X[i]+b2)){
			count++;
		}
		else if(points->X[i] < ((points->Y[i]-b3)/a3) && points->Y[i] < (a3*points->X[i]+b3)){
			count++;
		}
		else if(points->X[i] < ((points->Y[i]-b4)/a4) && points->Y[i] < (a4*points->X[i]+b4)){
			count++;
		}
		else if(points->X[i] > ((points->Y[i]-b5)/a5) && points->Y[i] < (a5*points->X[i]+b5)){
			count++;
		}
		else if(points->X[i] > ((points->Y[i]-b6)/a6) && points->Y[i] < (a6*points->X[i]+b6)){
			count++;
		}
		else if(points->X[i] > ((points->Y[i]-b7)/a7) && points->Y[i] > (a7*points->X[i]+b7)){
			count++;
		}
		else if(points->X[i] > ((points->Y[i]-b8)/a8) && points->Y[i] > (a8*points->X[i]+b8)){
			count++;
		}
	}
	//cout << count << endl;
	long resto = points->n - count;
	//cout << resto << endl;
	percent = (float)resto/(float)points->n;
	//cout << percent << endl;
	strcpy(aFile, "./RESULTS/");
	strcpy(str, "");
	sprintf(str, "EPSFilteredPerc"); //%ld, points->n
	strcat(aFile, str);
	cout << "Resume File: " << aFile << endl;
	FILE *fp = fopen(aFile, "a+" );
	if(NTHREADS > 0){
		if (points->NORMAL){
			// [n] [filtered percent] [esperanza] [varianza]
			fprintf(fp, "%ld %f %d %d\n", points->n, percent, points->NORMAL, NTHREADS);
		}else{
			// [n] [filtered percent] [avg bs-time/exec]
			fprintf(fp, "%ld %f %d %d\n", points->n, percent, points->NORMAL, NTHREADS);
		}
		fclose(fp);
	}
	else{
		if (points->NORMAL){
			// [n] [filtered percent] [esperanza] [varianza]
			fprintf(fp, "%ld %f %d %d\n", points->n, percent, points->NORMAL, 0);
		}else{
			// [n] [filtered percent] [avg bs-time/exec]
			fprintf(fp, "%ld %f %d %d\n", points->n, percent, points->NORMAL, 0);
		}
		fclose(fp);
	}*/
}

/**
 * Revisar la correctitud del algoritmo paralelo
 */
 void testPoints(pointSet *points){
	double t1, t2;
	float avgTime;
	ulong i;
	float xup1s, xup2s, xle2s, xle3s, xlo3s, xlo4s, xri4s, xri1s;
	float yup1s, yup2s, yle2s, yle3s, ylo3s, ylo4s, yri4s, yri1s;
	float xc1s, xc2s, xc3s, xc4s;
	float yc1s, yc2s, yc3s, yc4s;
	
	ulong ri1s, up1s;
	ulong up2s, le2s;
	ulong le3s, lo3s;
	ulong lo4s, ri4s;
	ulong c1s, c2s, c3s, c4s;
	
	t1=omp_get_wtime(); 
	for(i=0; i<REPET; i++){
		runEPS(points);
	}
	t2=omp_get_wtime();
	avgTime = (t2 - t1)/REPET;
	cout << "Average CPU time per execution, Sequential = " << avgTime*1000000.0 << " Microseconds" << endl;
	
	xup1s=xup1; xup2s=xup2; xle2s=xle2; xle3s=xle3; xlo3s=xlo3; xlo4s=xlo4; xri4s=xri4; xri1s=xri1;
	yup1s=yup1; yup2s=yup2; yle2s=yle2; yle3s=yle3; ylo3s=ylo3; ylo4s=ylo4; yri4s=yri4; yri1s=yri1;
	xc1s=xc1; xc2s=xc2; xc3s=xc3; xc4s=xc4;
	yc1s=yc1; yc2s=yc2; yc3s=yc3; yc4s=yc4;
	ri1s=ri1; up1s=up1;
	up2s=up2; le2s=le2;
	le3s=le3; lo3s=lo3;
	lo4s=lo4; ri4s=ri4;
	c1s=c1; c2s=c2; c3s=c3; c4s=c4;
	
	t1=omp_get_wtime(); 
	for(i=0; i<REPET; i++){
		runEPSPDAC(points);
		//~ runEPSPReduce(points);
	}
	t2=omp_get_wtime();
	avgTime = (t2 - t1)/REPET;
	cout << "Average CPU time per execution, Parallel DAC = " << avgTime*1000000.0 << " Microseconds" << endl;
	
	t1=omp_get_wtime(); 
	for(i=0; i<REPET; i++){
		//~ runEPSPDAC(points);
		runEPSPReduce(points);
	}
	t2=omp_get_wtime();
	avgTime = (t2 - t1)/REPET;
	cout << "Average CPU time per execution, Parallel Reduction = " << avgTime*1000000.0 << " Microseconds" << endl;
	
	runEPSPDAC(points);
	
	if(up1 != up1s) cout << "ERROR, up1: "<< up1s << " != " << up1 <<" The Indexes aren't equal" << endl;
	if(up2 != up2s) cout << "ERROR, up2: "<< up2s << " != " << up2 <<" The Indexes aren't equal" << endl;
	if(le2 != le2s) cout << "ERROR, le2: "<< le2s << " != " << le2 <<" The Indexes aren't equal" << endl;
	if(le3 != le3s) cout << "ERROR, le3: "<< le3s << " != " << le3 <<" The Indexes aren't equal" << endl;
	if(lo3 != lo3s) cout << "ERROR, lo3: "<< lo3s << " != " << lo3 <<" The Indexes aren't equal" << endl;
	if(lo4 != lo4s) cout << "ERROR, lo4: "<< lo4s << " != " << lo4 <<" The Indexes aren't equal" << endl;
	if(ri4 != ri4s) cout << "ERROR, ri4: "<< ri4s << " != " << ri4 <<" The Indexes aren't equal" << endl;
	if(ri1 != ri1s) cout << "ERROR, ri1: "<< ri1s << " != " << ri1 <<" The Indexes aren't equal" << endl;
	if(c1 != c1s) cout << "ERROR, c1: "<< c1s << " != " << c1 <<" The Indexes aren't equal" << endl;
	if(c2 != c2s) cout << "ERROR, c2: "<< c2s << " != " << c2 <<" The Indexes aren't equal" << endl;
	if(c3 != c3s) cout << "ERROR, c3: "<< c3s << " != " << c3 <<" The Indexes aren't equal" << endl;
	if(c4 != c4s) cout << "ERROR, c4: "<< c4s << " != " << c4 <<" The Indexes aren't equal" << endl;
	if(xup1 != xup1s) cout << "ERROR, xup1: "<< xup1s << " != " << xup1 <<" The Coordinates aren't equal" << endl;
	if(yup1 != yup1s) cout << "ERROR, yup1: "<< yup1s << " != " << yup1 <<" The Coordinates aren't equal" << endl;
	if(xup2 != xup2s) cout << "ERROR, xup2: "<< xup2s << " != " << xup2 <<" The Coordinates aren't equal" << endl;
	if(yup2 != yup2s) cout << "ERROR, yup2: "<< yup2s << " != " << yup2 <<" The Coordinates aren't equal" << endl;
	if(xle2 != xle2s) cout << "ERROR, xle2: "<< xle2s << " != " << xle2 <<" The Coordinates aren't equal" << endl;
	if(yle2 != yle2s) cout << "ERROR, yle2: "<< yle2s << " != " << yle2 <<" The Coordinates aren't equal" << endl;
	if(xle3 != xle3s) cout << "ERROR, xle3: "<< xle3s << " != " << xle3 <<" The Coordinates aren't equal" << endl;
	if(yle3 != yle3s) cout << "ERROR, yle3: "<< yle3s << " != " << yle3 <<" The Coordinates aren't equal" << endl;
	if(xlo3 != xlo3s) cout << "ERROR, xlo3: "<< xlo3s << " != " << xlo3 <<" The Coordinates aren't equal" << endl;
	if(ylo3 != ylo3s) cout << "ERROR, ylo3: "<< ylo3s << " != " << ylo3 <<" The Coordinates aren't equal" << endl;
	if(xlo4 != xlo4s) cout << "ERROR, xlo4: "<< xlo4s << " != " << xlo4 <<" The Coordinates aren't equal" << endl;
	if(ylo4 != ylo4s) cout << "ERROR, ylo4: "<< ylo4s << " != " << ylo4 <<" The Coordinates aren't equal" << endl;
	if(xri4 != xri4s) cout << "ERROR, xri4: "<< xri4s << " != " << xri4 <<" The Coordinates aren't equal" << endl;
	if(yri4 != yri4s) cout << "ERROR, yri4: "<< yri4s << " != " << yri4 <<" The Coordinates aren't equal" << endl;
	if(xri1 != xri1s) cout << "ERROR, xri1: "<< xri1s << " != " << xri1 <<" The Coordinates aren't equal" << endl;
	if(yri1 != yri1s) cout << "ERROR, yri1: "<< yri1s << " != " << yri1 <<" The Coordinates aren't equal" << endl;
	if(xc1 != xc1s) cout << "ERROR, xc1: "<< xc1s << " != " << xc1 <<" The Coordinates aren't equal" << endl;
	if(yc1 != yc1s) cout << "ERROR, yc1: "<< yc1s << " != " << yc1 <<" The Coordinates aren't equal" << endl;
	if(xc2 != xc2s) cout << "ERROR, xc2: "<< xc2s << " != " << xc2 <<" The Coordinates aren't equal" << endl;
	if(yc2 != yc2s) cout << "ERROR, yc2: "<< yc2s << " != " << yc2 <<" The Coordinates aren't equal" << endl;
	if(xc3 != xc3s) cout << "ERROR, xc3: "<< xc3s << " != " << xc3 <<" The Coordinates aren't equal" << endl;
	if(yc3 != yc3s) cout << "ERROR, yc3: "<< yc3s << " != " << yc3 <<" The Coordinates aren't equal" << endl;
	if(xc4 != xc4s) cout << "ERROR, xc4: "<< xc4s << " != " << xc4 <<" The Coordinates aren't equal" << endl;
	if(yc4 != yc4s) cout << "ERROR, yc4: "<< yc4s << " != " << yc4 <<" The Coordinates aren't equal" << endl;
 }
