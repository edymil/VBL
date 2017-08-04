/*
 *  CellsSystem-E.cpp
 *  Sim3D
 *
 *  Created by Edoardo Milotti on 22/04/10.
 *  Copyright 2010 I.N.F.N.-Sezione di Trieste. All rights reserved.
 * 
 *  This file takes care of mechanical forces
 *
 */

#include "sim.h"

#include "InputFromFile.h"
#include "CellType.h"
#include "Environment.h"
#include "EnvironmentalSignals.h"
#include "geom-2.h"
#include "BloodVessel.h"
#include "CellsSystem.h"
#include "Utilities.h"


// ***************************************************************
// Cluster mechanics
// ***************************************************************
//
// Cell dynamics (this method calculates the shift in a timestep)
//
// The algorithm used by this method corresponds to a semi-implicit Euler integration (explicit on forces, implied
// on the dissipative part; Should approximate approximately to an implicit verlet method)
// 
// Note that in the current version, Dynamics is not an interface for the GetForces and NewPositionsAndVelocities methods
//
void CellsSystem::Dynamics( )
{


	// prima parte, inizializzazione e calcolo della forza 
	// il calcolo della forza si fa una volta per tutte, corrisponde alla parte esplicita dell'integrazione

	GetForces( );
	
	// calcolo di nuove velocità e posizioni
	NewPositionsAndVelocities( );


		
}

// Calculating force between cells in the current configuration
void CellsSystem::GetForces()
{

	// ridimensionamento dei vettori
	fx.resize(ncells);
	fy.resize(ncells);
	fz.resize(ncells);
	
	// inizializzazione dei vettori
	fx.assign(ncells,0);
	fy.assign(ncells,0);
	fz.assign(ncells,0);
		
	// la forza cellula-cellula si calcola solo se ci sono almeno due cellule
	if(ncells>1)
		{
#pragma omp parallel for
		for(unsigned long n=0; n<ncells; n++)
			{
			
			// parametri relativi alla cellula n-esima
			double rn = (type[n]->Get_packing_factor())*r[n];					// raggio della cellula pesato con il packing factor
			double yn = (type[n]->Get_YoungMod());								// modulo di Young della cellula
			double pn = (type[n]->Get_PoissonRatio());							// PoissonRatio della cellula
			
			double arn = (type[n]->Get_adhesion_range());	// range di adesione del tipo della cellula n-esima
			double adn = (type[n]->Get_adhesion_decay());	// decay rate dell'adesione della cellula k-esima

			int nneigh = neigh[n];				// numero di cellule adiacenti alla cellula n-esima
			
			// loop sulle cellule adiacenti alla cellula n-esima (si calcola solo se c'e' almeno una cellula vicina)
			if(nneigh > 0)
				{
				for( int k=0; k< nneigh; k++)					
					{
					int neighbor = vneigh[n][k];	// nome della k-esima cellula adiacente
					double dd = Distance(n,neighbor);	// distanza tra cellula n-esima e la k-esima cellula vicina

					// parametri relativi alla vicina k-esima
					double rk = (type[neighbor]->Get_packing_factor())*r[neighbor];	// raggio della vicina k-esima pesato
					double yk = (type[neighbor]->Get_YoungMod());			// modulo di Young della cellula
					double pk = (type[neighbor]->Get_PoissonRatio());		// PoissonRatio della cellula

					double ark = (type[neighbor]->Get_adhesion_range());	// range di adesione del tipo della vicina k-esima
					double adk = (type[neighbor]->Get_adhesion_decay());	// decay rate dell'adesione della vicina k-esima

					double kC = sqrt(rn*rk)*(rn+rk)/(0.75*((1-pn*pn)/yn+(1-pk*pk)/yk));	// calcolo della "costante elastica" derivata dal modello di Hertz
					
					// qui si calcola la forza, che pero' quando e' repulsiva satura ad un valore max. Qui si sceglie di definire il valore max sulla base 
					// delle caratteristiche cellulari. Una scelta alternativa potrebbe essere quella di definire a priori il valore max per tutte le cellule. 
					// Qual e' la scelta giusta? Qui prendo la prima, arbitrariamente. 
					
					const double saturation_distance = 0.25;                               // distanza di saturazione della forza repulsiva (micron)
					// const double saturation_distance = 0.03;							// distanza di saturazione della forza repulsiva
					double fm = kC * pow(fabs((rn+rk-dd)/(rn+rk)),(double)1.5);		// modulo della forza
					if(rn+rk-dd > saturation_distance) 
						fm = kC * pow( saturation_distance/(rn+rk), (double)1.5);			// la forza repulsiva satura ad un valore abbastanza piccolo ... 
					
					if( dd > (rn+rk) )													// nel caso in cui la distanza sia maggiore della somma dei raggi la forza e' attrattiva
						fm *= -0.5*( 1. - tanh( 0.5*(adn+adk)*( dd/(rn+rk)-(1.+0.5*(arn+ark)) ) ) );	// la forza cambia segno e viene modulata
						
					fx[n] += fm*(x[n]-x[neighbor])/dd;								// si somma la proiezione di fm a ciascuna delle componenti della forza totale
					fy[n] += fm*(y[n]-y[neighbor])/dd;
					fz[n] += fm*(z[n]-z[neighbor])/dd;

		// sezione per il dump diagnostico, scommentare se serve
		/*
					cout << "\n*************************\n" << endl;
					cout << "forza di interazione tra le cellule: " << n << " e " << neighbor << endl;
					cout << "in posizione (" << scientific << x[n] << ", " << y[n] << ", " << z[n];
					cout << ");    (" << x[neighbor] << ", " << y[neighbor] << ", " << z[neighbor] << ");    " << endl;
					cout << "raggi cellulari: " << scientific << rn << "; " << rk << endl;
					cout << "range di adesione: " << scientific << arn << "; " << ark << endl;
					cout << "decay rate adesione: " << scientific << arn << "; " << ark << endl;
					cout << "distanza tra cellule: " << scientific << dd << endl;
					cout << "kC: " << scientific << kC << endl;
					cout << "modulo della forza: " << scientific << fm << " (" << sqrt( SQR(fx[n]) + SQR(fy[n]) + SQR(fz[n]) ) << ") " << "\n" << endl;
		*/
					
					}
				}
			
			}
		}


}


// Calculating the position and velocity of each cell (it is calculated only if there are at least two cells)
void CellsSystem::NewPositionsAndVelocities( )
{
	
	// inizializzazione dei vettori
	vxnew = vx;
	vynew = vy;
	vznew = vz;
	
	
	loop_count = 0;	// conteggio dei passaggi nel loop
	
	if( ncells > 1 )
		{
				
		while( true )	// loop infinito, viene interrotto da un break quando si raggiunge la condizione di stop
			{
			
			double dvmax = 0;	// max modulo della diff di velocita' tra due iterazioni

#pragma omp parallel for ordered schedule(dynamic)
			for(unsigned long n=0; n<ncells; n++)	// loop sulle cellule
				{
				
				int nneigh = neigh[n];					// numero di cellule adiacenti

				double rn = r[n];					// raggio della cellula (ora serve al calcolo del coeff di attrito)
				
				// inizializzazioni
				
				// coefficienti di attrito
				// ATTENZIONE ... questo statement dovrebbe entrare nel fenotipo
				double gamma = 6.*PI*VISCOSITY_ENV*rn;
				double gamma_int = 0.1e15;	// viscosita' interna in pg/(micron s)
				
				// elementi di matrice: questi elementi vengono calcolati al volo nel loop, ma potrebbero venire valutati prima del loop una volta per tutte
				// in questo modo si utilizzerebbe memoria, ma il programma diventerebbe piu' veloce ... potrebbe essere un'ottimizzazione da fare in futuro
				double axx = 0;			
				double axy = 0.;
				double axz = 0.;
				double ayy = 0;
				double ayz = 0.;
				double azz = 0;
				
				// vettore "costante"
				double bx = 0.;
				double by = 0.;
				double bz = 0.;
				
				double mn = mass[n];				// massa della cellula
				
				for( int k=0; k< nneigh; k++)							// loop sulle cellule adiacenti
					{
					int neighbor = vneigh[n][k];		// nome della k-esima cellula adiacente
					
					double drx = (x[n]-x[neighbor]);				// componenti del raggio vettore tra le cellule n e neighbor
					double dry = (y[n]-y[neighbor]);
					double drz = (z[n]-z[neighbor]);
					double dd2 = drx*drx+dry*dry+drz*drz;			// distanza al quadrato tra cellula n-esima e la k-esima cellula vicina
										
					double sp = vx[neighbor]*drx + vy[neighbor]*dry + vz[neighbor]*drz;
					
					bx += sp*drx/dd2;
					by += sp*dry/dd2;
					bz += sp*drz/dd2;
					
					axx += drx*drx/dd2;
					axy += drx*dry/dd2;
					axz += drx*drz/dd2;
					ayy += dry*dry/dd2;
					ayz += dry*drz/dd2;
					azz += drz*drz/dd2;
					
					}
				
				bx = vx[n] + (gamma_int*bx + fx[n])*dt/mn;
				by = vy[n] + (gamma_int*by + fy[n])*dt/mn;
				bz = vz[n] + (gamma_int*bz + fz[n])*dt/mn;
				
				axx = 1. + (gamma + gamma_int*axx)*dt/mn;
				axy = gamma_int*axy*dt/mn;
				axz = gamma_int*axz*dt/mn;
				ayy = 1. + (gamma + gamma_int*ayy)*dt/mn;
				ayz = gamma_int*ayz*dt/mn;
				azz = 1. + (gamma + gamma_int*azz)*dt/mn;
				
				double det = axx*ayy*azz + 2.*axy*axz*ayz - axx*ayz*ayz - ayy*axz*axz - azz*axy*axy;
				
				vxnew[n] = ( (ayy*azz-ayz*ayz)*bx + (axz*ayz-axy*azz)*by + (axy*ayz-axz*ayy)*bz )/det;
				vynew[n] = ( (axz*ayz-axy*azz)*bx + (axx*azz-axz*axz)*by + (axy*axz-ayz*axx)*bz )/det;
				vznew[n] = ( (axy*ayz-axz*ayy)*bx + (axy*axz-ayz*axx)*by + (axx*ayy-axy*axy)*bz )/det;

				double dv = sqrt((vx[n]-vxnew[n])*(vx[n]-vxnew[n])+(vy[n]-vynew[n])*(vy[n]-vynew[n])+(vz[n]-vznew[n])*(vz[n]-vznew[n]));
				if(dv > dvmax) dvmax = dv;

// statement per il debugging, scommentare se necessario
/*				
				cout << "\n------------------------" << endl;
				cout << "cellula: " << n << endl;
				cout << "massa: " << scientific << mn << endl;
				cout << "coefficiente di attrito con il mezzo: " << scientific << 6.*PI*VISCOSITY_ENV*rn << endl;
				cout << "matrice M:" << endl;
				cout << axx << "\t " << axy << "\t " << axz << endl;
				cout << axy << "\t " << ayy << "\t " << ayz << endl;
				cout << axz << "\t " << ayz << "\t " << azz << endl;
				cout << "\nmatrice inversa: " << endl;
				cout << (ayy*azz-ayz*ayz)/det << "\t " << (axz*ayz-axy*azz)/det << "\t " << (axy*ayz-axz*ayy)/det << endl;
				cout << (axz*ayz-axy*azz)/det << "\t " << (axx*azz-axz*axz)/det << "\t " << (axy*axz-ayz*axx)/det << endl;
				cout << (axy*ayz-axz*ayy)/det << "\t " << (axy*axz-ayz*axx)/det << "\t " << (axx*ayy-axy*axy)/det << endl;
				cout << "\nvettore b" << endl;
				cout << bx << "\t " << by << "\t " << bz << "\n " << endl;
				cout << "v:    (" << scientific << vx[n] << ", " << vy[n] << ", " << vz[n] << ") " << endl;
				cout << "vnew: (" << scientific << vxnew[n] << ", " << vynew[n] << ", " << vznew[n] << ") " << endl;
*/				
				}
			
			loop_count++; 
			
			//if(loop_count > 100) exit(0);
			// if( dvmax < (1.e-11) ) break;	// la condizione di interruzione del loop e' che l'imprecisione nella determinazione della velocità sia minore di 0.01 nm/s
			if( dvmax < delta_vmax ) break;	// la condizione di interruzione del loop e' che l'imprecisione nella determinazione della velocità sia minore di 0.1 nm/s

			for(unsigned long n=0; n<ncells; n++)	// qui si copiano in v i nuovi valori vnew
				{
				vx[n] = vxnew[n];
				vy[n] = vynew[n];
				vz[n] = vznew[n];
				}
			
			}
		}
	else 
	// soluzione nel caso ci sia un'unica cellula: l'unico coeff di attrito e' quello con l'ambiente e l'unico indice e' n=0
		{
		double rn = r[0];					// raggio dell'unica cellula (ora serve al calcolo del coeff di attrito)
        double mn = mass[0];				// massa della cellula
        double gamma = 6.*PI*VISCOSITY_ENV*rn;
		
		vx[0] = vx[0]/(1.+gamma*dt/mn);
		vy[0] = vy[0]/(1.+gamma*dt/mn);
		vz[0] = vz[0]/(1.+gamma*dt/mn);
		
		}

	
		
	// calcolo delle nuove posizioni
	
	maxdr = 0;	// inizializzazione dello spostamento massimo nel sistema di cellule

	for(unsigned long n=0; n<ncells; n++)				// loop sulle cellule
		{

		double rn = r[n];			// raggio della cellula

		double dx = vx[n]*dt;						// variazioni delle coordinate
		double dy = vy[n]*dt;
		double dz = vz[n]*dt;
			
		double dr = sqrt( dx*dx + dy*dy + dz*dz);	// spostamento totale
		if(dr>maxdr) maxdr = dr;						// update dello spostamento massimo
			
		if( dr > rn )		// gestione di possibili errori nel calcolo della posizione
			{
			cout << "Errore in Cells::NewPositionsAndVelocities: " << endl;
			cout << "Al passo " << nstep;
			cout << " si e' verificato un problema nel calcolo della posizione della cellula " << n << "-esima" << endl;
			cout << "dx = " << scientific << dx << endl;
			cout << "dy = " << scientific << dy << endl;
			cout << "dz = " << scientific << dz << endl;
			cout << "stepsize: " << scientific << dr << endl;
			cout << "cell radius: " << scientific << r[n] << endl;
			cout << "packing factor: " << (type[n]->Get_packing_factor()) << endl;
			cout << "cellForcex = " << scientific << fx[n] << endl;
			cout << "cellForcey = " << scientific << fy[n] << endl;
			cout << "cellForcez = " << scientific << fz[n] << endl;
			exit(0);
			}
			
		x[n] += dx;				// calcolo delle nuove coordinate
		y[n] += dy;
		z[n] += dz;

		
		}
		
}


// dinamica dummy
void CellsSystem::DummyDynamics( )
{
	maxdr = 0;	// spostamento massimo
}


