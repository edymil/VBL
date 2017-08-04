/*
 *  CellSystem-B.cpp
 *  VBL
 *
 *  Created by Edoardo Milotti on 20/04/09.
 *  Copyright 2009 I.N.F.N.-Sezione di Trieste. All rights reserved.
 * 
 *  This file is to manage spread in the cell network.
 *  ONLY ONE FUNCTION: void CellSystem::Diff()
 *
 */

#include "sim.h"

#include "InputFromFile.h"
#include "CellType.h"
#include "Environment.h"
#include "EnvironmentalSignals.h"
//#include "Pt.h"
#include "Cell.h"
#include "CellSystem.h"
#include "Utilities.h"


// se questo parametro e' true allora vengono stampati anche messaggi d'errore aggiuntivi
#define	EXTENDED_ERRORLOG	false

// scommentare il seguente statement per includere anche i messaggi d'errore per il metodo della secante
// #define SECANT_IN_ERRORLOG


//
// Method that manages the spread in the cell network
//
// New values ​​are calculated using nonlinear equations that come from a general solution with the BDF method
// (implicitly implied Euler); Nonlinear equations are in turn resolved by a Newton method
// 
void CellSystem::Diff()
{

	// allocazione della memoria per i vettori temporanei necessari al calcolo della diffusione
	
	long double* volumeOld = new long double[ncell];				// vettore dei volumi cellulari (valore vecchio)
	long double* volumeNew = new long double[ncell];				// vettore dei volumi cellulari (valore nuovo)
	long double* volume_extraOld = new long double[ncell];			// volume dei volume extracellulari (valore vecchio)
	long double* volume_extraNew = new long double[ncell];			// volume dei volume extracellulari (valore nuovo)
																	// anche se il volume non e' una variabile dinamica, va immagazzinato in vettori 
																	// per poter effettuare le somme necessarie al calcolo della diffusione
																
	long double* MitOld = new long double[ncell];					// vettore dei mitocondri (valore vecchio)
	long double* MitNew = new long double[ncell];					// vettore dei mitocondri (valore nuovo)

	long double* pHiOld = new long double[ncell];					// vettore dei pH cellulari (valore vecchio)
	long double* pHiNew = new long double[ncell];					// vettore dei pH cellulari (valore nuovo)
	long double* pHOld = new long double[ncell];					// vettore dei pH extracellulari (valore vecchio)
	long double* pHNew = new long double[ncell];					// vettore dei pH extracellulari (valore nuovo)
																	// anche se il pH non e' una variabile dinamica, va immagazzinato in vettori 
																	// per poter effettuare le somme necessarie al calcolo della diffusione
	
	long double* mGinOld = new long double[ncell];					// vettore della massa di glucosio dentro le cellule (valore vecchio)
	long double* mGinNew = new long double[ncell];					// vettore della massa di glucosio dentro le cellule (valore nuovo)
	long double* mGextOld = new long double[ncell];					// vettore della massa di glucosio negli spazi extracellulari (valore vecchio)
	long double* mGextNew = new long double[ncell];					// vettore della massa di glucosio negli spazi extracellulari (valore nuovo)

	long double* mG6POld = new long double[ncell];					// vettore della massa di G6P nelle cellule (valore vecchio)
	long double* mG6PNew = new long double[ncell];					// vettore della massa di G6P nelle cellule (valore nuovo)

	long double* mO2Old = new long double[ncell];					// vettore della massa di O2 nelle cellule (valore vecchio)
	long double* mO2New = new long double[ncell];					// vettore della massa di O2 nelle cellule (valore nuovo)

	long double* StoreOld = new long double[ncell];					// vettore dello Store nelle cellule (valore vecchio)
	long double* StoreNew = new long double[ncell];					// vettore dello Store nelle cellule (valore nuovo)

	long double* mAinOld = new long double[ncell];					// vettore della massa di glutammina dentro le cellule (valore vecchio)
	long double* mAinNew = new long double[ncell];					// vettore della massa di glutammina dentro le cellule (valore nuovo)
	long double* mAextOld = new long double[ncell];					// vettore della massa di glutammina negli spazi extracellulari (valore vecchio)
	long double* mAextNew = new long double[ncell];					// vettore della massa di glutammina negli spazi extracellulari (valore nuovo)

	long double* mAcLinOld = new long double[ncell];				// vettore della massa di acido lattico dentro le cellule (valore vecchio)
	long double* mAcLinNew = new long double[ncell];				// vettore della massa di acido lattico dentro le cellule (valore nuovo)
	long double* mAcLextOld = new long double[ncell];				// vettore della massa di acido lattico negli spazi extracellulari (valore vecchio)
	long double* mAcLextNew = new long double[ncell];				// vettore della massa di acido lattico negli spazi extracellulari (valore nuovo)

	long double* ATPpOld = new long double[ncell];					// vettore dell'ATPp nelle cellule (valore vecchio)
	long double* ATPpNew = new long double[ncell];					// vettore dell'ATPp nelle cellule (valore nuovo)

	// altre allocazioni di memoria
	
	long double* protein = new long double[ncell];					// vettore della massa di proteina nelle cellule
	long double* pRb = new long double[ncell];						// vettore della massa di pRb nelle cellule
	long double* delta_protein = new long double[ncell];			// vettore della variazione di massa totale di proteine nelle cellule
	long double* ConcS = new long double[ncell];					// concentrazione del substrato S per il calcolo delle soglie dei checkpoints

	long double* DNA = new long double[ncell];						// vettore che contiene la lunghezza relativa della molecola di DNA (1 = intera molecola)
	long double* DNA_rate = new long double[ncell];					// vettore che contiene il rate di sintesi del DNA
	
	double* GAbsRate = new double[ncell];							// vettore che memorizza l'assorbimento di glucosio
	double* GConsRate = new double[ncell];							// vettore che memorizza il consumo di glucosio

	double* AAbsRate = new double[ncell];							// vettore che memorizza l'assorbimento di glutammina
	double* AConsRate = new double[ncell];							// vettore che memorizza il consumo di glutammina
	
	double* StoreFillRate = new double[ncell];						// vettore che memorizza il rate di riempiemento dello store
	double* StoreConsRate = new double[ncell];						// vettore che memorizza il rate di consumo dello store
	
	double* AcLRate = new double[ncell];							// vettore che memorizza il rate di produzione dell'acido lattico
	double* AcLOutRate = new double[ncell];							// vettore che memorizza il rate di espulsione dell'acido lattico
	
	long double* ATP_Ox = new long double[ncell];					// vettore che memorizza il rate ATP_Ox
	long double* ATP_NOx = new long double[ncell];					// vettore che memorizza il rate ATP_NOx
	long double* ATP2 = new long double[ncell];						// vettore che memorizza il rate ATP2
	long double* ATP3 = new long double[ncell];						// vettore che memorizza il rate ATP3
	long double* ConsATP = new long double[ncell];					// vettore che memorizza il rate ConsATP
	long double* ConsATP_1 = new long double[ncell];				// vettore che memorizza il rate ConsATP_1
	long double* ConsATP_2 = new long double[ncell];				// vettore che memorizza il rate ConsATP_2
	long double* ConsATP_3 = new long double[ncell];				// vettore che memorizza il rate ConsATP_3
	long double* ConsATP_5 = new long double[ncell];				// vettore che memorizza il rate ConsATP_5
	
	long double ATPprod;
	long double ATPcons;

	// assegnazioni preliminari fondamentali per il loop del calcolo di diffusione e metabolismo (vengono fatte comunque, anche per le cellule morte)

#pragma omp parallel for
	for(unsigned long n=0; n<ncell; n++)						
		{
		
		// controllo preliminare di consistenza
		if( CellVector[n].Get_phase() != dead ) 
			{
			int code = CellVector[n].CheckMVA();
			if(code < 0) errorlog_file << "Errore " << code << " all'inizio di CellSystem::Diff nel controllo di consistenza per la cellula " << n << "\n" << endl;
			}

		// volume
		volumeOld[n] = volumeNew[n] = CellVector[n].Get_volume();
		volume_extraOld[n] = volume_extraNew[n] = CellVector[n].Get_volume_extra();
		
		// mitocondri
		MitOld[n] = MitNew[n] = CellVector[n].Get_M();
		
		// pH
		pHiOld[n] = pHiNew[n] = CellVector[n].Get_pHi();
		pHOld[n] = pHNew[n] = CellVector[n].Get_pH();
				
		// glucosio
		mGinOld[n] = mGinNew[n] = CellVector[n].Get_G();
		mGextOld[n] = mGextNew[n] = CellVector[n].Get_G_extra();
		
		// G6P
		mG6POld[n] = mG6PNew[n] = CellVector[n].Get_G6P();
		
		// O2
		mO2Old[n] = mO2New[n] = CellVector[n].Get_O2();
		
		// Store
		StoreOld[n] = StoreNew[n] = CellVector[n].Get_store();
		
		// glutammina
		mAinOld[n] = mAinNew[n] = CellVector[n].Get_A();
		mAextOld[n] = mAextNew[n] = CellVector[n].Get_A_extra();
		
		// AcL
		mAcLinOld[n] = mAcLinNew[n] = CellVector[n].Get_AcL();
		mAcLextOld[n] = mAcLextNew[n] = CellVector[n].Get_AcL_extra();
		
		// ATP
		ATPpOld[n] = ATPpNew[n] = CellVector[n].Get_ATPp();

		}

	// altre assegnazioni
	
#pragma omp parallel for
	for(unsigned long n=0; n<ncell; n++)						
		{
		
		protein[n] = CellVector[n].Get_protein();
		DNA[n] = CellVector[n].Get_DNA();
		
		}
	
	// valori ambientali
	
	long double volume_envOld = Env.GetEnvironmentvolume();		// volume libero dell'ambiente
	long double volume_envNew = volume_envOld;
	
	long double mG_envOld = Env.GetEnvironmentG();				// massa di glucosio nell'ambiente
	long double mG_envNew = mG_envOld;
	
	long double mO2_envOld = Env.GetEnvironmentO2();			// massa di O2 nell'ambiente
	long double mO2_envNew = mO2_envOld;
	
	long double mA_envOld = Env.GetEnvironmentA();				// massa di glutammina nell'ambiente
	long double mA_envNew = mA_envOld;
	
	long double mAcL_envOld = Env.GetEnvironmentAcL();			// massa di acido lattico nell'ambiente
	long double mAcL_envNew = mAcL_envOld;
	

	// * loop per la soluzione delle equazioni nonlineari *
	
	bool isOK;													// variabile che ferma il loop
	ncalls++;													// update del numero di volte che si inizializza il loop
	nrepeats=0;													// inizializzazione del numero di ripetizioni del loop
	
	long double prec=eps;										// precisione ottenuta

	do {														// inizio del loop
	
	isOK = true;
	
	// prima si calcolano i nuovi valori ambientali 
	
	long double s0_env = 0.;									// somme che servono per il calcolo della diffusione nell'ambiente
	long double sG_env = 0.;
	// long double sO2_env = 0.;								// la somma per l'ossigeno e' attualmente inattiva perche' l'ossigeno viene fissato al valore atmosferico
	long double sA_env = 0.;
	long double sAcL_env = 0.;
	long double cell_volume = 0.;								// volume totale delle cellule 
	
#pragma omp parallel for ordered schedule(dynamic)
	for(unsigned long n=0; n<ncell; n++)						// loop sulle cellule per il calcolo delle somme
		{
		// volume totale occupato dalle cellule
		cell_volume += volumeOld[n];
		
		// le somme per il calcolo della diffusione nell'ambiente sono non-zero solo per cellule che appartengono all'alpha shape
		if( CellVector[n].Get_isonAS() )
			{
			long double g_env = CellVector[n].Get_g_env();
			s0_env += g_env;
			sG_env += g_env*(mGextOld[n]/volume_extraOld[n]);
			// sO2_env += g_env*(mO2Old[n]/volumeOld[n]);
			sA_env += g_env*(mAextOld[n]/volume_extraOld[n]);
			sAcL_env += g_env*(mAcLextOld[n]/volume_extraOld[n]);
			}
			
		}

	
	volume_envNew = Env_0.GetEnvironmentvolume() - cell_volume;

	if( flowON && ready2start )				// se si lavora con il bioreattore allora le equazioni per l'ambiente devono tenere conto di termini in piu' ... 
		{
		mG_envNew = ( Env.GetEnvironmentG() + dt*( Diff_W_G*sG_env + Env_0.GetEnvironmentG()/Env.GetEnvironmentvolume() * flowSignal.SignalValue(treal) ) )\
			/( 1. + dt*(Diff_W_G*s0_env)/Env.GetEnvironmentvolume() + dt * flowSignal.SignalValue(treal)/volume_envOld );
		
		// la concentrazione di ossigeno nell'ambiente e' fissa
		mO2_envNew = (Env_0.GetEnvironmentO2()/Env_0.GetEnvironmentvolume()) * volume_envOld;									
		
		mA_envNew = ( Env.GetEnvironmentA() + dt*( Diff_W_A*sA_env + Env_0.GetEnvironmentA()/Env.GetEnvironmentvolume() * flowSignal.SignalValue(treal) ) )\
			/( 1. + dt*(Diff_W_A*s0_env)/Env.GetEnvironmentvolume() + dt * flowSignal.SignalValue(treal)/volume_envOld );
		
		// l'acido lattico nella soluzione nutritiva pulita e' sempre zero
		mAcL_envNew = ( Env.GetEnvironmentAcL() + dt*( Diff_W_AcL*sAcL_env ) )\
			/( 1. + dt*(Diff_W_AcL*s0_env)/Env.GetEnvironmentvolume() + dt * flowSignal.SignalValue(treal)/volume_envOld );
		}
	else if( !flowON && ready2start )		// ... altrimenti si utilizzano le equazioni senza i termini aggiuntivi (questo e' vero anche per l'inizializzazione)
		{
		mG_envNew = ( Env.GetEnvironmentG() + dt*( Diff_W_G*sG_env ) )\
			/( 1. + dt*(Diff_W_G*s0_env)/volume_envOld );
		
		// la concentrazione di ossigeno nell'ambiente e' fissa
		mO2_envNew = (Env_0.GetEnvironmentO2()/Env_0.GetEnvironmentvolume()) * volume_envOld;									
		
		mA_envNew = ( Env.GetEnvironmentA() + dt*( Diff_W_A*sA_env ) )\
			/( 1. + dt*(Diff_W_A*s0_env)/volume_envOld );
		
		// l'acido lattico nella soluzione nutritiva pulita e' sempre zero
		mAcL_envNew = ( Env.GetEnvironmentAcL() + dt*( Diff_W_AcL*sAcL_env ) )\
			/( 1. + dt*(Diff_W_AcL*s0_env)/volume_envOld );
		}
	else if ( !ready2start )				// infine, durante la fase di inizializzazione si tengono fisse le concentrazioni ambientali
		{
		mG_envNew = (Env_0.GetEnvironmentG()/Env_0.GetEnvironmentvolume()) * volume_envOld;
		
		mO2_envNew = (Env_0.GetEnvironmentO2()/Env_0.GetEnvironmentvolume()) * volume_envOld;
		
		mA_envNew = (Env_0.GetEnvironmentA()/Env_0.GetEnvironmentvolume()) * volume_envOld;
		
		mAcL_envNew = (Env_0.GetEnvironmentAcL()/Env_0.GetEnvironmentvolume()) * volume_envOld;
		}

	// P.S. ... si noti che nell'implementazione attuale il bioreattore funziona solo introducendo liquido di coltura pulito, uguale a quello iniziale
	// se si vogliono introdurre dosi diverse (ad esempio valori di glutammina piu' alti rispetto al valore iniziale), allora si deve modificare questa parte 
	// si noti anche che il tempo di partenza (fase 0) degli impulsi variabili nel tempo è riferito al tempo di inizializzazione


	// qui si passa al calcolo dei nuovi valori per le cellule e spazi extracellulari

	O2Flow = 0.;	// inizializzazione della somma per il calcolo del rate di assorbimento dell'ossigeno

#pragma omp parallel for ordered schedule(dynamic)
	for(unsigned long n=0; n<ncell; n++)						// loop sulle cellule
		{
		
		
		
		// coefficienti
		long double vmax2 = CellVector[n].Get_type()->Get_VMAX_2();
		long double vmax22 = CellVector[n].Get_type()->Get_VMAX_22();
		long double vmaxP = CellVector[n].Get_type()->Get_VMAX_P();
		long double vmaxP_A = CellVector[n].Get_type()->Get_VMAX_P_A();
		long double vmaxP_ATP = CellVector[n].Get_type()->Get_VMAX_P_ATP();
		long double vmaxDNA = CellVector[n].Get_type()->Get_VMAX_DNA();
		long double v_WORK = CellVector[n].Get_type()->Get_v_WORK();
		long double C1 = CellVector[n].Get_type()->Get_C1();		
		long double C2 = CellVector[n].Get_type()->Get_C2();		
		long double Vmin = CellVector[n].Get_type()->Get_Vmin();		
		long double vmaxDNA_A = CellVector[n].Get_type()->Get_VMAX_DNA_A();
		long double vmaxDNA_ATP = CellVector[n].Get_type()->Get_VMAX_DNA_ATP();
		long double vmaxM = CellVector[n].Get_type()->Get_VMAX_M();
		long double vmaxM_A = CellVector[n].Get_type()->Get_VMAX_M_A();
		long double vmaxM_ATP = CellVector[n].Get_type()->Get_VMAX_M_ATP();
		long double Km1 = CellVector[n].Get_type()->Get_Km1();
		long double Km2 = CellVector[n].Get_type()->Get_Km2();
		long double Km22 = CellVector[n].Get_type()->Get_Km22();
		long double Kmc = CellVector[n].Get_type()->Get_Kmc();
		long double Kmd = CellVector[n].Get_type()->Get_Kmd();
		long double Ka = CellVector[n].Get_type()->Get_Ka();
		long double KmA = CellVector[n].Get_type()->Get_KmA();
		long double KmAL = CellVector[n].Get_type()->Get_KmAL();
		long double Kmp = CellVector[n].Get_type()->Get_Kmp();
		long double KmDNA = CellVector[n].Get_type()->Get_KmDNA();
		long double KmM = CellVector[n].Get_type()->Get_KmM();
		long double ATPSt = CellVector[n].Get_type()->Get_ATPSt();

		long double coeffg1 = CellVector[n].Get_type()->Get_coeffg1();
		long double coeffg2 = CellVector[n].Get_type()->Get_coeffg2();
		long double coeffg3 = CellVector[n].Get_type()->Get_coeffg3();
		long double coeffr1 = CellVector[n].Get_type()->Get_coeffr1();
		
		
		// quantita' derivate (attenzione, l'ordine e' importante)
		long double tpH = 0.5 * ( 1. + tanh( CellVector[n].Get_type()->Get_tph_slope()*(pHiOld[n] - CellVector[n].Get_type()->Get_tph_thr()) ) );
		long double tp11 = 0.5 * ( 1. + tanh( CellVector[n].Get_type()->Get_tp11_slope()*(pHiOld[n] - CellVector[n].Get_type()->Get_tp11_thr()) ) );
		long double a2c = 0.5 * ( 1. + tanh( CellVector[n].Get_type()->Get_a2c_slope()*(pHOld[n] - CellVector[n].Get_type()->Get_a2c_thr()) ) );
		long double c2a = 0.5 * ( 1. + tanh( CellVector[n].Get_type()->Get_c2a_slope()*(pHiOld[n] - CellVector[n].Get_type()->Get_c2a_thr()) ) );
		long double a2cA = 0.5 * ( 1. + tanh( CellVector[n].Get_type()->Get_a2cA_slope()*(pHOld[n] - CellVector[n].Get_type()->Get_a2cA_thr()) ) );
		long double c2aA = 0.5 * ( 1. + tanh( CellVector[n].Get_type()->Get_c2aA_slope()*(pHiOld[n] - CellVector[n].Get_type()->Get_c2aA_thr()) ) );
		long double a2cAcL = 2. - tanh( CellVector[n].Get_type()->Get_a2cAcL_slope()*( pHOld[n]-CellVector[n].Get_type()->Get_a2cAcL_thr() ) );
		long double c2aAcL = 2. - tanh( CellVector[n].Get_type()->Get_c2aAcL_slope()*( pHiOld[n]-CellVector[n].Get_type()->Get_c2aAcL_thr() ) );
		
		// derivata di a2cAcL rispetto la massa di acido lattico nel volume extracellulare: e' l'unica derivata che si deve calcolare al momento
		// ma questo dovrebbe cambiare nel momento in cui si inserira' il nuovo modello di acidita' cellulare
		long double da2cAcL = (1./(volume_extraOld[n]*BufCapEnv))/pow(cosh( CellVector[n].Get_type()->Get_a2cAcL_slope()*( pHOld[n]-CellVector[n].Get_type()->Get_a2cAcL_thr() ) ),2);
		
		
		if( CellVector[n].Get_phase() != dead )
			pHiNew[n] = pHi_STANDARD;	// ***** al momento il pH interno e' fisso *****
										// ***** questo significa anche che le derivate rispetto pHi sono nulle !!! *****
		else
			pHiNew[n] = pHOld[n];		// se la cellula è morta ha un pH interno uguale a quello esterno ... 
		
		pHNew[n] = 7.5443-( mAcLextOld[n]/volume_extraOld[n])/BufCapEnv;
		
		long double r = pow(3.*volumeOld[n]/(4.*PI), (long double)1./3.);
		
		// long double surf_M_corr;
		// if( CellVector[n].Get_phase() == M_phase )
		//	surf_M_corr = 1.+0.129960524947437*( tanh( 5.*( CellVector[n].Get_phase_age()/CellVector[n].Get_M_T() - 0.5 ) ) + 1. );	// fattore correttivo per la mitosi
		// else 
		// 	surf_M_corr = 1.;
		//	
		// long double surface = 4.*PI*r*r*surf_M_corr;
		
		long double surface = 4.*PI*r*r;

		
		if( CellVector[n].Get_phase() != dead )
			{
			// se la cellula e' viva il volume cambia
			volumeNew[n] = Vmin * (1+DNA[n]) + C2*MitOld[n] + C1*ATPpOld[n];		
			volume_extraNew[n] = surface*CellVector[n].Get_type()->Get_extvolume_thickness()*CellVector[n].Get_type()->Get_extvolume_fraction();
			}
		else
			{
			// se la cellula e' morta il volume resta fisso (viene ridotto solo alla fine di questo metodo)
			volumeNew[n] = volumeOld[n];
			volume_extraNew[n] = volume_extraOld[n];
			}
		
		long double ATP_St = ATPSt;
		long double SensO2 = 0.;
		long double dSensO2_O2 = 0.;
		long double SensATP = 0.;
		
		if( CellVector[n].Get_phase() != dead )
			{
			SensO2 = mO2Old[n]/(volumeOld[n]*CellVector[n].Get_type()->Get_KmO2() + mO2Old[n]);
			dSensO2_O2 = volumeOld[n]*CellVector[n].Get_type()->Get_KmO2()/pow(volumeOld[n]*CellVector[n].Get_type()->Get_KmO2() + mO2Old[n],2);
			ATP_Ox[n] = 30.*(PM_ATP/PM_G)*SensO2 * ( coeffg2 + coeffg3*StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n]) ) * mG6POld[n];
			SensATP = 0.5*( 1 - tanh( 100.*(ATP_Ox[n]/ATP_St-1.) ) );
			ATP_NOx[n] = 2.*(PM_ATP/PM_G)*tpH*( coeffg1*mG6POld[n] + coeffr1*StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n]) );
			ATP2[n] = SensO2*SensATP*(ATP_St - ATP_Ox[n])*StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n]) ;
			ATP3[n] = (14./5.)*SensO2*SensATP*(ATP_St - ATP_Ox[n]) * mAinOld[n]/(volumeOld[n]*Kmd+mAinOld[n]);
			ConsATP[n] = 2.*(PM_ATP/PM_G) * tp11*tpH * coeffr1 * (1.-SensATP) * mAinOld[n]/(volumeOld[n]*Kmd+mAinOld[n]);
			}
		else
			ATP_Ox[n] = ATP_NOx[n] = ATP2[n] = ConsATP[n] = 0.;
		
		long double vp, dvp_A, dvp_ATPp;					// produzione di proteine senza il termine vmax 
		// la produzione di proteine e' bloccata in fase G0
		if(CellVector[n].Get_phase() != G0_phase && CellVector[n].Get_phase() != dead)	
			{
			vp = ATPpOld[n]*mAinOld[n]/(pow(volumeOld[n],2)*Kmp+ATPpOld[n]*mAinOld[n]);
			dvp_A = ATPpOld[n]*pow(volumeOld[n],2)*Kmp/pow(pow(volumeOld[n],2)*Kmp+ATPpOld[n]*mAinOld[n],2);
			dvp_ATPp = pow(volumeOld[n],2)*Kmp*mAinOld[n]/pow(pow(volumeOld[n],2)*Kmp+ATPpOld[n]*mAinOld[n],2);
			}
		else
			vp= dvp_A = dvp_ATPp = 0.;
		
		long double vDNA, dvDNA_A, dvDNA_ATPp;			// produzione di DNA senza il termine vmax
		if(CellVector[n].Get_phase() == S_phase)	// la produzione di DNA procede solo in fase S
			{
			vDNA = ATPpOld[n]*mAinOld[n]/(pow(volumeOld[n],2)*KmDNA+ATPpOld[n]*mAinOld[n]);
			dvDNA_A = ATPpOld[n]*pow(volumeOld[n],2)*KmDNA/pow(pow(volumeOld[n],2)*KmDNA+ATPpOld[n]*mAinOld[n],2);
			dvDNA_ATPp = pow(volumeOld[n],2)*KmDNA*mAinOld[n]/pow(pow(volumeOld[n],2)*KmDNA+ATPpOld[n]*mAinOld[n],2);
			}
		else
			vDNA = dvDNA_A = dvDNA_ATPp = 0.;
		DNA_rate[n] = vmaxDNA*vDNA;
		
		// produzione di mtDNA: termine comune a mitocondri, glutammina e ATPp senza le vmax, e sue derivate rispetto A e ATPp
		long double vM, dvM_A, dvM_ATPp;
		if(CellVector[n].Get_phase() != dead)
			{
			vM = ATPpOld[n]*mAinOld[n]/(pow(volumeOld[n],2)*KmM+ATPpOld[n]*mAinOld[n]);
			dvM_A = ATPpOld[n]*pow(volumeOld[n],2)*KmM/pow(pow(volumeOld[n],2)*KmM+ATPpOld[n]*mAinOld[n],2);
			dvM_ATPp = pow(volumeOld[n],2)*KmM*mAinOld[n]/pow(pow(volumeOld[n],2)*KmM+ATPpOld[n]*mAinOld[n],2);
			}
		else
			vM = dvM_A = dvM_ATPp = 0.;



		long double h = 0.5*(1.3*(1.-mO2Old[n]/(volumeOld[n]*O2st))+1.)*(1.+tanh(100.*(1-mO2Old[n]/(volumeOld[n]*O2st)))) + 0.5*(1-tanh(100.*(1.-mO2Old[n]/(volumeOld[n]*O2st))));
		
		long double v1max = CellVector[n].Get_type()->Get_VMAX_1()*surface*h;
		long double vmaxA = CellVector[n].Get_type()->Get_VMAX_A()*surface;
		long double vmaxAcL = CellVector[n].Get_type()->Get_VmaxAL0() * surface;
		
		// fattore geometrico verso l'ambiente
		long double g_env = CellVector[n].Get_g_env();

		// somme 
		long double s0 = 0.;
		long double sG = 0.;
		long double sO2 = 0.;
		long double sA = 0.;
		long double sAcL = 0.;
		
		// calcolo delle somme per la diffusione nel cluster di cellule
		for(int nv=0; nv<CellVector[n].Get_neigh(); nv++)		
			{
			int name = (CellVector[n].Get_vneigh())[nv];
			long double gnk = CellVector[n].Get_gnk()[nv];
			s0 += gnk;
			sG += gnk*mGextOld[name]/volume_extraOld[name];			
			sO2 += gnk*mO2Old[name]/volumeOld[name];			
			sA += gnk*mAextOld[name]/volume_extraOld[name];		
			sAcL += gnk*mAcLextOld[name]/volume_extraOld[name];			
			}
		


		// *** glucosio

		long double M_Gnow, DM_Gnow, T_Gnow, DT_Gnow_in, DT_Gnow_ext;	// valori che servono al calcolo con il metodo di Newton
		long double M_Gnow_0, T_Gnow_0in, T_Gnow_0ext;					// valori aggiuntivi per il calcolo con il metodo della secante
	
		if( CellVector[n].Get_phase() != dead )
			{
			// funzione metabolica	
			M_Gnow = -vmax2*mGinOld[n]*mGinOld[n]/((volumeOld[n]*Km2+mGinOld[n])*(volumeOld[n]*Ka+mGinOld[n])) \
					- vmax22*mGinOld[n]*mGinOld[n]/((volumeOld[n]*Km22+mGinOld[n])*(volumeOld[n]*Ka+mGinOld[n]));
			// derivata funzione metabolica rispetto mGin calcolata in mGinOld
			DM_Gnow = -2*vmax2*mGinOld[n]/((volumeOld[n]*Km2+mGinOld[n])*(volumeOld[n]*Ka+mGinOld[n])) \
					+ vmax2*mGinOld[n]*mGinOld[n]/(pow(volumeOld[n]*Km2+mGinOld[n],2)*(volumeOld[n]*Ka+mGinOld[n])) \
					+ vmax2*mGinOld[n]*mGinOld[n]/((volumeOld[n]*Km2+mGinOld[n])*pow(volumeOld[n]*Ka+mGinOld[n],2)) \
					- 2*vmax22*mGinOld[n]/((volumeOld[n]*Km22+mGinOld[n])*(volumeOld[n]*Ka+mGinOld[n])) \
					+ vmax22*mGinOld[n]*mGinOld[n]/(pow(volumeOld[n]*Km22+mGinOld[n],2)*(volumeOld[n]*Ka+mGinOld[n])) \
					+ vmax22*mGinOld[n]*mGinOld[n]/((volumeOld[n]*Km22+mGinOld[n])*pow(volumeOld[n]*Ka+mGinOld[n],2));
			
			M_Gnow_0 = 0.;	// valore della funzione metabolica in mGin = 0
					
			// funzione di trasporto
			T_Gnow = a2c*v1max*mGextOld[n]/(volume_extraOld[n]*Km1+mGextOld[n]) - c2a*v1max*mGinOld[n]/(volumeOld[n]*Km1+mGinOld[n]);
			// derivata della funzione di trasporto rispetto mGin
			DT_Gnow_in = -c2a*v1max*volumeOld[n]*Km1/pow(volumeOld[n]*Km1+mGinOld[n],2);
			// derivata della funzione di trasporto rispetto mGext
			DT_Gnow_ext = a2c*v1max*volume_extraOld[n]*Km1/pow(volume_extraOld[n]*Km1+mGextOld[n],2);
			
			// valori speciali della funzione di trasporto
			T_Gnow_0in = a2c*v1max*mGextOld[n]/(volume_extraOld[n]*Km1+mGextOld[n]);	// mGin = 0
			T_Gnow_0ext = - c2a*v1max*mGinOld[n]/(volumeOld[n]*Km1+mGinOld[n]);			// mGext = 0
			
			}
		else
			{
			M_Gnow = DM_Gnow = T_Gnow = DT_Gnow_in = DT_Gnow_ext = 0.;	// se la cellula e' morta non c'e' alcuna attivita' metabolica o di trasporto
			M_Gnow_0 = T_Gnow_0in = T_Gnow_0ext = 0; 
			}
		
		// 1. calcolo del nuovo valore di mGin
		if( CellVector[n].Get_phase() != dead )
			mGinNew[n] = mGinOld[n] - ( mGinOld[n] - ( CellVector[n].Get_G() + dt*( M_Gnow + T_Gnow ) ) )/( 1. - dt*( DM_Gnow + DT_Gnow_in ) );
		else
			mGinNew[n] = mGinOld[n];

		if( mGinNew[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
			{
			long double fhi = mGinOld[n] - ( CellVector[n].Get_G() + dt*( M_Gnow + T_Gnow ) );	// f valutata in mGinOld[n]
			long double flo = - ( CellVector[n].Get_G() + dt*( M_Gnow_0 + T_Gnow_0in ) );		// f valutata in mGin = 0
			long double xsol = mGinOld[n]*flo/(flo-fhi);
			mGinNew[n] = xsol;
#ifdef SECANT_IN_ERRORLOG
			errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
			errorlog_file << " e' necessario passare al metodo della secante" << endl;
			errorlog_file << scientific << "mGinOld["<< n << "] = " << mGinOld[n] << "; mGinNew["<< n << "] = " << mGinNew[n] << "\n" <<endl;
#endif
			}
		if( mGinNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
			{
			errorlog_file << "ATTENZIONE: al passo " << nstep << " mGinNew["<< n << "] = " << scientific << mGinNew[n] << " < 0" << endl;
			errorlog_file << "in questa iterazione si mantiene mGinNew[n] = mGinOld[n] = " << mGinOld[n] << "\n" << endl;
			mGinNew[n] = mGinOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
			}

		// 2. calcolo del nuovo valore di mGext
		long double cG = CellVector[n].Get_G_extra() + dt*Diff_ES_G*sG + dt*Diff_W_G*(mG_envOld/volume_envOld)*g_env;	// calcolo della parte fissa
		
		mGextNew[n] = mGextOld[n] - ( mGextOld[n]*(1+Diff_ES_G*dt*s0/volume_extraOld[n]+Diff_W_G*dt*g_env/volume_extraOld[n]) + dt*T_Gnow - cG )\
				/( 1. + Diff_ES_G*dt*s0/volume_extraOld[n] + Diff_W_G*dt*g_env/volume_extraOld[n] + dt*DT_Gnow_ext );
		if(mGextNew[n] < 0)	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
			{
			long double fhi = ( mGextOld[n]*(1+Diff_ES_G*dt*s0/volume_extraOld[n]+Diff_W_G*dt*g_env/volume_extraOld[n]) + dt*T_Gnow - cG );
			long double flo = ( dt*T_Gnow_0ext - cG );
			long double xsol = mGextOld[n]*flo/(flo-fhi);
			mGextNew[n] = xsol;
#ifdef SECANT_IN_ERRORLOG
			errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
			errorlog_file << " e' necessario passare al metodo della secante" << endl;
			errorlog_file << scientific << "mGextOld["<< n << "] = " << mGextOld[n] << "; mGextNew["<< n << "] = " << mGextNew[n] << "\n" <<endl;
#endif
			}
		if( mGextNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
			{
			errorlog_file << "ATTENZIONE: al passo " << nstep << " mGextNew["<< n << "] = " << scientific << mGextNew[n] << " < 0" << endl;
			errorlog_file << "in questa iterazione si mantiene mGextNew[n] = mGextOld[n] = " << mGextOld[n] << "\n" << endl;
			mGextNew[n] = mGextOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
			}
		
		// memorizzazione di consumo e trasporto
		GConsRate[n] = (double)M_Gnow;
		GAbsRate[n] = (double)T_Gnow;
	

		
		// *** G6P (soluzione esatta dell'equazione iterativa derivata dal metodo di Eulero)
		if( CellVector[n].Get_phase() != dead )
			mG6PNew[n] = (mG6POld[n] - dt*M_Gnow)/( 1. + dt*( tpH*coeffg1 + coeffg2*SensO2 + coeffg3) );
		else
			mG6PNew[n] = mG6POld[n];



		// *** O2
				
		long double M_O2now, DM_O2now;								// valori che servono al calcolo con il metodo di Newton
		long double M_O2now_0;										// valore aggiuntivo per il calcolo con il metodo della secante
		
		if( CellVector[n].Get_phase() != dead )
			{
			M_O2now = -6.*(PM_O2/PM_G) * SensO2 * ( coeffg2*mG6POld[n] \
				+ coeffg3*mG6POld[n]*(StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n])) \
				+ SensATP*0.033333*(PM_G/PM_ATP)*( ATP_St - ATP_Ox[n] ) * (StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n])) \
				+ 0.1*SensATP*(PM_G/PM_ATP)*( ATP_St - ATP_Ox[n] ) * (mAinOld[n]/(volumeOld[n]*Kmd+mAinOld[n]))  );
			DM_O2now = -6.*(PM_O2/PM_G) * dSensO2_O2 * ( coeffg2*mG6POld[n] \
				+ coeffg3*mG6POld[n]*(StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n])) \
				+ SensATP*0.033333*(PM_G/PM_ATP)*( ATP_St - ATP_Ox[n] ) * (StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n])) \
				+ 0.1*SensATP*(PM_G/PM_ATP)*( ATP_St - ATP_Ox[n] ) * (mAinOld[n]/(volumeOld[n]*Kmd+mAinOld[n]))  );
				
			M_O2now_0 = 0.;	// valore della funzione metabolica in mO2 = 0
			
			}
		else
			{
			M_O2now = DM_O2now = 0.;	// se la cellula e' morta non c'e' alcun consumo metabolico
			M_O2now_0 = 0.;
			}
			
		O2Flow += M_O2now;
			
		// calcolo del nuovo valore
		long double cO2 = CellVector[n].Get_O2() + dt*Diff_ES_O2*sO2 + dt*Diff_W_O2*(mO2_envOld/volume_envOld)*g_env;	// calcolo della parte fissa
		
		mO2New[n] = mO2Old[n] - ( mO2Old[n] * ( 1+Diff_ES_O2*dt*s0/volumeOld[n]+Diff_W_O2*dt*g_env/volumeOld[n]) - dt*M_O2now - cO2 )\
				/( 1 + Diff_ES_O2*dt*s0/volumeOld[n] + Diff_W_O2*dt*g_env/volumeOld[n] - dt*DM_O2now);
				
		if( mO2New[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
			{
			long double fhi = mO2Old[n] * ( 1+Diff_ES_O2*dt*s0/volumeOld[n]+Diff_W_O2*dt*g_env/volumeOld[n]) - dt*M_O2now - cO2;
			long double flo = - dt*M_O2now_0 - cO2;
			long double xsol = mO2Old[n]*flo/(flo-fhi);
			mO2New[n] = xsol;
#ifdef SECANT_IN_ERRORLOG
			errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
			errorlog_file << " e' necessario passare al metodo della secante" << endl;
			errorlog_file << scientific << "mO2Old["<< n << "] = " << mO2Old[n] << "; mO2New["<< n << "] = " << mO2New[n] << "\n" <<endl;
#endif
			}
		if( mO2New[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
			{
			errorlog_file << "ATTENZIONE: al passo " << nstep << " mO2New["<< n << "] = " << scientific << mO2New[n] << " < 0" << endl;
			errorlog_file << "in questa iterazione si mantiene mO2New[n] = mO2Old[n] = " << mO2Old[n] << "\n" << endl;
			mO2New[n] = mO2Old[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
			}


		// *** glutammina
		
		long double M_Anow, DM_Anow, T_Anow, DT_Anow_in, DT_Anow_ext;	// valori che servono al calcolo con il metodo di Newton
		long double M_Anow_0, T_Anow_0in, T_Anow_0ext;					// valori aggiuntivi per il calcolo con il metodo della secante
		
		if( CellVector[n].Get_phase() != dead )
			{
			// funzione metabolica
			M_Anow = -( tp11*tpH*coeffr1*(1-SensATP) + 0.1*SensO2*SensATP*(PM_G/PM_ATP)*( ATP_St - ATP_Ox[n] ) ) * (mAinOld[n]/(volumeOld[n]*Kmd+mAinOld[n])) \
				- vmaxP_A*vp - vmaxDNA_A*vDNA - vmaxM_A*vM;
			// derivata funzione metabolica rispetto mAin calcolata in mAinOld
			DM_Anow = -( tp11*tpH*coeffr1*(1-SensATP) + 0.1*SensO2*SensATP*(PM_G/PM_ATP)*( ATP_St - ATP_Ox[n] ) ) * (volumeOld[n]*Kmd/pow(volumeOld[n]*Kmd+mAinOld[n],2)) \
				- vmaxP_A*dvp_A - vmaxDNA_A*dvDNA_A - vmaxM_A*dvM_A;
				
			M_Anow_0 = 0;	// valore della funzione metabolica in mAin = 0
			
			// funzione di trasporto
			T_Anow = a2cA*vmaxA * mAextOld[n]/(volume_extraOld[n]*KmA+mAextOld[n]) - c2aA*vmaxA * mAinOld[n]/(volumeOld[n]*KmA+mAinOld[n]);
			// derivata della funzione di trasporto rispetto mAin
			DT_Anow_in = - c2aA*vmaxA * volumeOld[n]*KmA/pow(volumeOld[n]*KmA+mAinOld[n],2);
			// derivata della funzione di trasporto rispetto mAext
			DT_Anow_ext = a2cA*vmaxA * volume_extraOld[n]*KmA/pow(volume_extraOld[n]*KmA+mAextOld[n],2);

			// valori speciali della funzione di trasporto
			T_Anow_0in = a2cA*vmaxA * mAextOld[n]/(volume_extraOld[n]*KmA+mAextOld[n]);	// mAin = 0
			T_Anow_0ext = - c2aA*vmaxA * mAinOld[n]/(volumeOld[n]*KmA+mAinOld[n]);		// mAext = 0
			
			}
		else
			{
			M_Anow = DM_Anow = T_Anow = DT_Anow_in = DT_Anow_ext = 0.;
			M_Anow_0 = T_Anow_0in = T_Anow_0ext = 0.;
			}
		
		// 1. calcolo del nuovo valore di mAin
		if( CellVector[n].Get_phase() != dead )
			mAinNew[n] = mAinOld[n] - ( mAinOld[n] - ( CellVector[n].Get_A() + dt*( M_Anow + T_Anow ) ) )/( 1. - dt*( DM_Anow + DT_Anow_in ) );
		else
			mAinNew[n] = mAinOld[n];
			
		if( mAinNew[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
			{
			long double fhi = mAinOld[n] - ( CellVector[n].Get_A() + dt*( M_Anow + T_Anow ) );
			long double flo = - ( CellVector[n].Get_A() + dt*( M_Anow_0 + T_Anow_0in ) );
			long double xsol = mAinOld[n]*flo/(flo-fhi);
			mAinNew[n] = xsol;
#ifdef SECANT_IN_ERRORLOG
			errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
			errorlog_file << " e' necessario passare al metodo della secante" << endl;
			errorlog_file << scientific << "mAinOld["<< n << "] = " << mAinOld[n] << "; mAinNew["<< n << "] = " << mAinNew[n] << "\n" <<endl;
#endif
			}
		if( mAinNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
			{
			errorlog_file << "ATTENZIONE: al passo " << nstep << " mAinNew["<< n << "] = " << scientific << mAinNew[n] << " < 0" << endl;
			errorlog_file << "in questa iterazione si mantiene mAinNew[n] = mAinOld[n] = " << mAinOld[n] << "\n" << endl;
			mAinNew[n] = mAinOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
			}
		
		// 2. calcolo del nuovo valore di mAext
		long double cA = CellVector[n].Get_A_extra() + dt*Diff_ES_A*sA + dt*Diff_W_A*(mA_envOld/volume_envOld)*g_env;	// calcolo della parte fissa
		
		mAextNew[n] = mAextOld[n] - ( mAextOld[n]*(1+Diff_ES_A*dt*s0/volume_extraOld[n]+Diff_W_A*dt*g_env/volume_extraOld[n]) + dt*T_Anow - cA )\
				/( 1. + Diff_ES_A*dt*s0/volume_extraOld[n] + Diff_W_A*dt*g_env/volume_extraOld[n] + dt*DT_Anow_ext );
		if( mAextNew[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
			{
			long double fhi = ( mAextOld[n]*(1+Diff_ES_A*dt*s0/volume_extraOld[n]+Diff_W_A*dt*g_env/volume_extraOld[n]) + dt*T_Anow - cA );
			long double flo = ( dt*T_Anow_0ext - cA );
			long double xsol = mAextOld[n]*flo/(flo-fhi);
			mAextNew[n] = xsol;
#ifdef SECANT_IN_ERRORLOG
			errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
			errorlog_file << " e' necessario passare al metodo della secante" << endl;
			errorlog_file << scientific << "mAextOld["<< n << "] = " << mAextOld[n] << "; mAextNew["<< n << "] = " << mAextNew[n] << "\n" <<endl;
#endif
			}
		if( mAextNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
			{
			errorlog_file << "ATTENZIONE: al passo " << nstep << " mAextNew["<< n << "] = " << scientific << mAextNew[n] << " < 0" << endl;
			errorlog_file << "in questa iterazione si mantiene mAextNew[n] = mAextOld[n] = " << mAextOld[n] << "\n" << endl;
			mAextNew[n] = mAextOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
			}

		// memorizzazione di consumo e trasporto
		AConsRate[n] = (double)M_Anow;
		AAbsRate[n] = (double)T_Anow;


		
		
		// *** store (soluzione esatta dell'equazione iterativa derivata dal metodo di Eulero)

		if( CellVector[n].Get_phase() != dead )
			{
			long double Ast = CellVector[n].Get_store() + dt * ( coeffg3*mG6POld[n] + tp11*tpH*coeffr1*(mAinOld[n]/(volumeOld[n]*Kmd+mAinOld[n]))*(1-SensATP) );
			long double Bst = -dt * ( tpH*coeffr1 + coeffg3*mG6POld[n]*SensO2 + SensO2*SensATP*0.033333*(PM_G/PM_ATP)*( ATP_St - ATP_Ox[n] ) );
			long double Kst = volumeOld[n]*Kmc;

			StoreNew[n] = 0.5*( (Ast+Bst-Kst) + sqrt( pow(Ast+Bst-Kst,2) + 4*Ast*Kst ) );
			
			// memorizzazione di riempimento e consumo
			StoreFillRate[n] = (double)( coeffg3*mG6POld[n] + tp11*tpH*coeffr1*(mAinOld[n]/(volumeOld[n]*Kmd+mAinOld[n]))*(1-SensATP) );
			StoreConsRate[n] = (double)(-( tpH*coeffr1 + coeffg3*mG6POld[n]*SensO2 + SensO2*SensATP*0.033333*(PM_G/PM_ATP)*( ATP_St - ATP_Ox[n] ) ) * (StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n])));
			}
		else
			{
			StoreNew[n] = StoreOld[n];
			StoreFillRate[n] = StoreConsRate[n] = 0.;
			}

		// *** AcL 
		
		long double M_AcLnow, DM_AcLnow, T_AcLnow, DT_AcLnow_in, DT_AcLnow_ext;
		long double M_AcLnow_0, T_AcLnow_0in, T_AcLnow_0ext;
		
		if( CellVector[n].Get_phase() != dead )
			{
			// funzione metabolica
			M_AcLnow = 2.*tpH*( coeffg1*mG6POld[n] + coeffr1*(StoreOld[n]/(volumeOld[n]*Kmc+StoreOld[n])) );
			
			// PRESTARE GRANDE ATTENZIONE ALLE DERIVATE RISPETTO AcL DEI COEFFICIENTI CHE REGOLANO IL TRASPORTO ... 
			// qui sotto SI ASSUME CHE LE DERIVATE RISPETTO pHi SIANO NULLE perche' al momento pHi e' tenuto costante e non cambia in funzione di mAcLin.
			
			// derivata funzione metabolica rispetto mAcLin calcolata in mAcLinOld			
			DM_AcLnow = 0.;
			
			M_AcLnow_0 = M_AcLnow;	// valore della funzione metabolica in mAcLin = 0: nella formulazione attuale la produzione di AcL non dipende dal contenuto di AcL della cellula
			
			// funzione di trasporto
			T_AcLnow = a2cAcL*vmaxAcL * mAcLextOld[n]/(volume_extraOld[n]*KmAL+mAcLextOld[n]) - c2aAcL*vmaxAcL * mAcLinOld[n]/(volumeOld[n]*KmAL+mAcLinOld[n]);
			// derivata della funzione di trasporto rispetto mAcLin
			DT_AcLnow_in = - c2aAcL*vmaxAcL * volumeOld[n]*KmAL/pow(volumeOld[n]*KmAL+mAcLinOld[n],2);
			// derivata della funzione di trasporto rispetto mAcLext			
			DT_AcLnow_ext = da2cAcL*vmaxAcL * mAcLextOld[n]/(volume_extraOld[n]*KmAL+mAcLextOld[n]) + \
							a2cAcL*vmaxAcL * volume_extraOld[n]*KmAL/pow(volume_extraOld[n]*KmAL+mAcLextOld[n],2);
			
			// valori speciali della funzione di trasporto
			T_AcLnow_0in = a2cAcL*vmaxAcL * mAcLextOld[n]/(volume_extraOld[n]*KmAL+mAcLextOld[n]);	// mAcLin = 0
			T_AcLnow_0ext = - c2aAcL*vmaxAcL * mAcLinOld[n]/(volumeOld[n]*KmAL+mAcLinOld[n]);		// mAcLext = 0
			
			}
		else
			{
			M_AcLnow = DM_AcLnow = T_AcLnow = DT_AcLnow_in = DT_AcLnow_ext = 0.;
			M_AcLnow_0 = T_AcLnow_0in = T_AcLnow_0ext = 0;
			}

		// 1. calcolo del nuovo valore di mAcLin
		if( CellVector[n].Get_phase() != dead )
			mAcLinNew[n] = mAcLinOld[n] - ( mAcLinOld[n] - ( CellVector[n].Get_AcL() + dt*( M_AcLnow + T_AcLnow ) ) )/( 1 - dt*( DM_AcLnow + DT_AcLnow_in ) );
		else
			mAcLinNew[n] = mAcLinOld[n];
			
		if( mAcLinNew[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
			{
			long double fhi = mAcLinOld[n] - ( CellVector[n].Get_AcL() + dt*( M_AcLnow + T_AcLnow ) );	// f valutata in mAcLinOld[n]
			long double flo = - ( CellVector[n].Get_AcL() + dt*( M_AcLnow_0 + T_AcLnow_0in ) );			// f valutata in mAcLin = 0
			long double xsol = mAcLinOld[n]*flo/(flo-fhi);
			mAcLinNew[n] = xsol;
#ifdef SECANT_IN_ERRORLOG
			errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
			errorlog_file << " e' necessario passare al metodo della secante" << endl;
			errorlog_file << scientific << "mAcLinOld["<< n << "] = " << mAcLinOld[n] << "; mAcLinNew["<< n << "] = " << mAcLinNew[n] << "\n" <<endl;
#endif
			}
		if( mAcLinNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
			{
			errorlog_file << "ATTENZIONE: al passo " << nstep << " mAcLinNew["<< n << "] = " << scientific << mAcLinNew[n] << " < 0" << endl;
			errorlog_file << "in questa iterazione si mantiene mAcLinNew[n] = mAcLinOld[n] = " << mAcLinOld[n] << "\n" << endl;
			mAcLinNew[n] = mAcLinOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
			}
		
		// 2. calcolo del nuovo valore di mAcLext
		long double cAcL = CellVector[n].Get_AcL_extra() + dt*Diff_ES_AcL*sAcL + dt*Diff_W_AcL*(mAcL_envOld/volume_envOld)*g_env;	// calcolo della parte fissa
		
		mAcLextNew[n] = mAcLextOld[n] - ( mAcLextOld[n]*(1+Diff_ES_AcL*dt*s0/volume_extraOld[n]+Diff_W_AcL*dt*g_env/volume_extraOld[n]) + dt*T_AcLnow - cAcL )\
				/( 1. + Diff_ES_AcL*dt*s0/volume_extraOld[n] + Diff_W_AcL*dt*g_env/volume_extraOld[n] + dt*DT_AcLnow_ext );
		if( mAcLextNew[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
			{
			long double fhi = mAcLextOld[n]*(1+Diff_ES_AcL*dt*s0/volume_extraOld[n]+Diff_W_AcL*dt*g_env/volume_extraOld[n]) + dt*T_AcLnow - cAcL;
			long double flo = dt*T_AcLnow_0ext - cAcL;
			long double xsol = mAcLextOld[n]*flo/(flo-fhi);
			mAcLextNew[n] = xsol;
#ifdef SECANT_IN_ERRORLOG
			errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
			errorlog_file << " e' necessario passare al metodo della secante" << endl;
			errorlog_file << scientific << "mAcLextOld["<< n << "] = " << mAcLextOld[n] << "; mAcLextNew["<< n << "] = " << mAcLextNew[n] << "\n" <<endl;
#endif
			}
		if( mAcLextNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
			{
			errorlog_file << "ATTENZIONE: al passo " << nstep << " mAcLextNew["<< n << "] = " << scientific << mAcLextNew[n] << " < 0" << endl;
			errorlog_file << "in questa iterazione si mantiene mAcLextNew[n] = mAcLextOld[n] = " << mAcLextOld[n] << "\n" << endl;
			mAcLextNew[n] = mAcLextOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
			}

		// memorizzazione di consumo e trasporto
		AcLRate[n] = (double)M_AcLnow;
		AcLOutRate[n] = (double)T_AcLnow;



		// *** ATP
		
		if( CellVector[n].Get_phase() != dead )
			{
			
			long double M_ATPnow, DM_ATPnow;								// valori che servono al calcolo con il metodo di Newton
			
			// memorizzazione di produzione e consumo
			ConsATP_1[n] = v_WORK*C1*ATPpOld[n];
			ConsATP_2[n] = vmaxP_ATP*vp;
			ConsATP_3[n] = vmaxDNA_ATP*vDNA;
			ConsATP_5[n] = vmaxM_ATP*vM;
			
			ATPprod = ATP_Ox[n] + ATP_NOx[n] + ATP2[n] + ATP3[n];
			ATPcons = ConsATP[n] + ConsATP_1[n] + ConsATP_2[n] + ConsATP_3[n] + ConsATP_5[n];

			// funzione metabolica
			M_ATPnow = ATPprod - ATPcons;
			// derivata funzione metabolica rispetto ATPp calcolata in ATPpOld
			DM_ATPnow = - v_WORK*C1 - vmaxP_ATP*dvp_ATPp - vmaxDNA_ATP*dvDNA_ATPp - vmaxM_ATP*dvM_ATPp; 

			// calcolo del nuovo valore di ATPp
			ATPpNew[n] = ATPpOld[n] - ( ATPpOld[n] - CellVector[n].Get_ATPp() - dt*M_ATPnow )/(1.-dt*DM_ATPnow);
			if( ATPpNew[n] < 0 )	// se si arriva a questo punto vuol dire che il metodo di Newton ha dei problemi ... 
				{
				long double fhi = ATPpOld[n] - CellVector[n].Get_ATPp() - dt*M_ATPnow;
				long double flo = - ( CellVector[n].Get_ATPp() + dt*(ATP_NOx[n]) );
				long double xsol = ATPpOld[n]*flo/(flo-fhi);
				ATPpNew[n] = xsol;
#ifdef SECANT_IN_ERRORLOG
				errorlog_file << "ATTENZIONE: al passo " << nstep << ", iterazione " << nrepeats << ", cellula " << n;
				errorlog_file << " e' necessario passare al metodo della secante" << endl;
				errorlog_file << scientific << "ATPpOld["<< n << "] = " << ATPpOld[n] << "; ATPpNew["<< n << "] = " << ATPpNew[n] << "\n" <<endl;
#endif
				}
			if( ATPpNew[n] < 0 )	// se si arriva a questo punto vuol dire che anche il metodo della secante ha dei problemi ... 
				{
				errorlog_file << "ATTENZIONE: al passo " << nstep << " ATPpNew["<< n << "] = " << scientific << ATPpNew[n] << " < 0" << endl;
				errorlog_file << "in questa iterazione si mantiene ATPpNew[n] = ATPpOld[n] = " << ATPpOld[n] << "\n" << endl;
				errorlog_file << "Altre informazioni sullo stato cellulare attuale: " << endl;
				errorlog_file << "eta' cellulare " << CellVector[n].Get_age() << endl;
				errorlog_file << "eta' di fase cellulare " << CellVector[n].Get_phase_age() << endl;
				errorlog_file << "ATP_St " << CellVector[n].Get_ATP_St() << endl;
				errorlog_file << "ATP_Ox " << CellVector[n].Get_ATP_Ox() << endl;
				errorlog_file << "ATP_NOx " << CellVector[n].Get_ATP_NOx() << endl;
				errorlog_file << "ATP2 " << CellVector[n].Get_ATP2() << endl;
				errorlog_file << "ATP3 " << CellVector[n].Get_ATP3() << endl;
				errorlog_file << "ConsATP " << CellVector[n].Get_ConsATP() << endl;
				errorlog_file << "ConsATP_1 " << CellVector[n].Get_ConsATP_1() << endl;
				errorlog_file << "ConsATP_2 " << CellVector[n].Get_ConsATP_2() << endl;
				errorlog_file << "ConsATP_3 " << CellVector[n].Get_ConsATP_3() << endl;
				errorlog_file << "ConsATP_4 " << CellVector[n].Get_ConsATP_4() << endl;
				errorlog_file << "ConsATP_5 " << CellVector[n].Get_ConsATP_5() << endl;
				errorlog_file << "ATPtot " << CellVector[n].Get_ATPtot() << endl;
				errorlog_file << "ATPp " << CellVector[n].Get_ATPp() << endl;
				errorlog_file << "ATPmin " << CellVector[n].Get_ATPmin() << endl;
				errorlog_file << endl;
				ATPpNew[n] = ATPpOld[n];	// allora si mantiene temporanemente il vecchio valore e si confida nella prossima iterazione ... 
				}
			// long double deltaATPp = ATPpNew[n]-CellVector[n].Get_ATPp();
			
			// proliferazione dei mitocondri (la proliferazione dei mitocondri e' bloccata in fase G0)
			if(CellVector[n].Get_phase() != G0_phase)
				MitNew[n] = CellVector[n].Get_M() + vmaxM*vM*dt; 
			else
				MitNew[n] = CellVector[n].Get_M();
			}
		else
			{
			ConsATP_1[n] = 0.;
			ConsATP_2[n] = 0.;
			ConsATP_3[n] = 0.;
			ConsATP_5[n] = 0.;
			ATPpNew[n] = ATPpOld[n];
			MitNew[n] = CellVector[n].Get_M();
			}
			
		
		// *** proteine
		
		if( CellVector[n].Get_phase() != dead )
			{
			delta_protein[n] = vmaxP*vp*dt;
			protein[n] = CellVector[n].Get_protein() + delta_protein[n];
			
			if( CellVector[n].Get_phase() == G2_phase || CellVector[n].Get_phase() == M_phase)
				pRb[n] = CellVector[n].Get_pRb() + CellVector[n].Get_type()->Get_pRb_fraction()*delta_protein[n];
			else
				pRb[n] = CellVector[n].Get_pRb();
				
				
			// calcolo del numero di molecole di pRb attive
			// calcolo approssimato di p, di Pk_pRb e del numero di molecole attive NpRbk
			// 
			// Qui si assume un processo di fosforilazione di questo tipo: 
			// 1. la ciclina D si lega all'ATP 
			// 2. il complesso si lega alla kinasi
			// 3. tutto questo complessone va a fosforilare la pRb sui singoli siti
			// 
			// diamo per scontato che i passi 1 e 2 possano procedere (in altre parole assumiamo che l'ATP necessario sia 
			// molto superiore a quello della soglia di morte)
			// Quello che chiamiamo ConcCyclinD e' in realta' la concentrazione del complesso fosforilante

			long double ConcpRb = pRb[n]/(PM_pRb*1000. * CellVector[n].Get_volume());			// concentrazione MOLARE (!!!!) della pRb nella cellula
			long double ConcCyclinD = CellVector[n].Get_cyclinD()/( PM_cyclinD*1000. * CellVector[n].Get_volume() ); // concentrazione MOLARE (!!!!) della ciclina D nella cellula
			long double ConcCyclinE = CellVector[n].Get_cyclinE()/( PM_cyclinE*1000. * CellVector[n].Get_volume() ); // concentrazione MOLARE (!!!!) della ciclina E nella cellula

			int N_pRb = CellVector[n].Get_type()->Get_N_pRb();		// numero di siti di fosforilazione della pRb
			int k_pRb = CellVector[n].Get_type()->Get_k_pRb();		// numero di siti di soglia
			long double pRb_ONOFFratio = CellVector[n].Get_type()->Get_pRb_ONOFFratio();

			long double p_pRb = ( (N_pRb * ConcpRb + (ConcCyclinD+ConcCyclinE) + pRb_ONOFFratio) - \
					sqrt( pow(N_pRb * ConcpRb + (ConcCyclinD+ConcCyclinE) + pRb_ONOFFratio,2) - 4.*N_pRb*ConcpRb*(ConcCyclinD+ConcCyclinE) ) )/(2. * N_pRb * ConcpRb);
			
			long double Pk_pRb = 0.;
			for(int l=k_pRb; l<=N_pRb; l++)
				{
				Pk_pRb += bico(N_pRb,k_pRb) * pow(p_pRb,l)*pow(1.-p_pRb,N_pRb-l);
				}
			
			long double ConcE = ConcpRb*Pk_pRb;	// concentrazione MOLARE dell'enzima

			// reazione downstream (MM)
			ConcS[n] = 0.5*( (CellVector[n].Get_ConcS() - dt*ConcE*CellVector[n].Get_type()->Get_k3MM() - CellVector[n].Get_type()->Get_KmMM()) \
				+ sqrt( pow(CellVector[n].Get_ConcS() - dt*ConcE*CellVector[n].Get_type()->Get_k3MM() - CellVector[n].Get_type()->Get_KmMM(),2) \
						+ 4.*CellVector[n].Get_ConcS()*CellVector[n].Get_type()->Get_KmMM() ) );
			}
		else
			{
			delta_protein[n] = 0.;
			protein[n] = CellVector[n].Get_protein();
			pRb[n] = CellVector[n].Get_pRb();
			ConcS[n] = CellVector[n].Get_ConcS();
			}
		
		
		// *** DNA
		
		if(CellVector[n].Get_phase() == S_phase)
			DNA[n] = CellVector[n].Get_DNA() + DNA_rate[n]*dt;
		else
			DNA[n] = CellVector[n].Get_DNA();

															

		}	// fine del loop sulle cellule

	// controllo di convergenza per le variabili ambientali
	if( fabs(mG_envOld - mG_envNew) > eps*0.5*(fabs(mG_envOld)+fabs(mG_envNew))+TOL ) 
		{
		isOK = false;
		convergence_fail[0]++;
		long double newprec = 2.*fabs(mG_envOld - mG_envNew)/(fabs(mG_envOld)+fabs(mG_envNew));
		prec = newprec > prec ? newprec : prec;
		}
//	if( fabs(mO2_envOld - mO2_envNew) > eps*0.5*(fabs(mO2_envOld)+fabs(mO2_envNew))+TOL ) 
//		{
//		isOK = false;
//		convergence_fail[1]++;
//		long double newprec = 2.*fabs(mO2_envOld - mO2_envNew)/(fabs(mO2_envOld)+fabs(mO2_envNew));
//		prec = newprec > prec ? newprec : prec;
//		}
	if( fabs(mA_envOld - mA_envNew) > eps*0.5*(fabs(mA_envOld)+fabs(mA_envNew))+TOL ) 
		{
		isOK = false;
		convergence_fail[2]++;
		long double newprec = 2.*fabs(mA_envOld - mA_envNew)/(fabs(mA_envOld)+fabs(mA_envNew));
		prec = newprec > prec ? newprec : prec;
		}

	if( fabs(mAcL_envOld - mAcL_envNew) > eps*0.5*(fabs(mAcL_envOld)+fabs(mAcL_envNew))+TOL )
		{
		isOK = false;
		convergence_fail[3]++;
		long double newprec = 2*fabs(mAcL_envOld - mAcL_envNew)/(fabs(mAcL_envOld)+fabs(mAcL_envNew));
		prec = newprec > prec ? newprec : prec;
		}

	// assegnazione delle variabili ambientali
	volume_envOld = volume_envNew;
	mG_envOld = mG_envNew;
	mO2_envOld = mO2_envNew;
	mA_envOld = mA_envNew;
	mAcL_envOld = mAcL_envNew;

	int ncell_fails = 0;										// variabile che conta il numero di cellule per cui fallisce la convergenza
	
#pragma omp parallel for ordered schedule(dynamic)
	for(unsigned long n=0; n<ncell; n++)						// assegnazione dei nuovi valori in preparazione alla prossima iterazione
		{
		
		bool cellisOK = true;
		
		// if( CellVector[n].Get_phase() != dead )					// le istruzioni di controllo e assegnazione si eseguono solo se la cellula e' viva
			{
			// prima si controlla se l'algoritmo e' arrivato a convergenza ...
			if( fabs(mGinOld[n] - mGinNew[n]) > eps*0.5*(fabs(mGinOld[n])+fabs(mGinNew[n]))+TOL ) 
				{
				isOK = false;
				cellisOK = false;
				convergence_fail[4]++;
				long double newprec = 2.*fabs(mGinOld[n] - mGinNew[n])/(fabs(mGinOld[n])+fabs(mGinNew[n]));
				prec = newprec > prec ? newprec : prec;
				if ( fabs(mGinOld[n] - mGinNew[n]) > 0.5*(fabs(mGinOld[n])+fabs(mGinNew[n])) && EXTENDED_ERRORLOG ) 
					{
					errorlog_file << "CellSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << endl;
					errorlog_file << "mGin: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
					errorlog_file << "\tdifferenza : " << 100.*(mGinNew[n]-mGinOld[n])/mGinOld[n] << "%\n";
					errorlog_file << "\tfase : " << CellVector[n].Get_phase() << endl;
					errorlog_file << "\teta' di fase: " << CellVector[n].Get_phase_age() << " s\n" << endl;
					}
				}
			if( fabs(mGextOld[n] - mGextNew[n]) > eps*0.5*(fabs(mGextOld[n])+fabs(mGextNew[n]))+TOL )
				{
				isOK = false;
				cellisOK = false;
				convergence_fail[5]++;
				long double newprec = 2.*fabs(mGextOld[n] - mGextNew[n])/(fabs(mGextOld[n])+fabs(mGextNew[n]));
				prec = newprec > prec ? newprec : prec;
				if ( fabs(mGextOld[n] - mGextNew[n]) > 0.5*(fabs(mGextOld[n])+fabs(mGextNew[n])) && EXTENDED_ERRORLOG ) 
					{
					errorlog_file << "CellSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << endl;
					errorlog_file << "mGext: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
					errorlog_file << "\tdifferenza : " << 100.*(mGextNew[n]-mGextOld[n])/mGextOld[n] << "%\n";
					errorlog_file << "\tfase : " << CellVector[n].Get_phase() << endl;
					errorlog_file << "\teta' di fase: " << CellVector[n].Get_phase_age() << " s\n" << endl;
					}
				}
			if( fabs(mO2Old[n] - mO2New[n]) > eps*0.5*(fabs(mO2Old[n])+fabs(mO2New[n]))+TOL )
				{
				isOK = false;
				cellisOK = false;
				convergence_fail[6]++;
				long double newprec = 2.*fabs(mO2Old[n] - mO2New[n])/(fabs(mO2Old[n])+fabs(mO2New[n]));
				prec = newprec > prec ? newprec : prec;
				if ( fabs(mO2Old[n] - mO2New[n]) > 0.5*(fabs(mO2Old[n])+fabs(mO2New[n])) && EXTENDED_ERRORLOG ) 
					{
					errorlog_file << "CellSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << endl;
					errorlog_file << "mO2: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
					errorlog_file << "\tdifferenza : " << 100.*(mO2New[n]-mO2Old[n])/mO2Old[n] << "%\n";
					errorlog_file << "\tfase : " << CellVector[n].Get_phase() << endl;
					errorlog_file << "\teta' di fase: " << CellVector[n].Get_phase_age() << " s\n" << endl;
					}
				}
			if( fabs(mAinOld[n] - mAinNew[n]) > eps*0.5*(fabs(mAinOld[n])+fabs(mAinNew[n]))+TOL )
				{
				isOK = false;
				cellisOK = false;
				convergence_fail[7]++;
				long double newprec = 2.*fabs(mAinOld[n] - mAinNew[n])/(fabs(mAinOld[n])+fabs(mAinNew[n]));
				prec = newprec > prec ? newprec : prec;
				if ( fabs(mAinOld[n] - mAinNew[n]) > 0.5*(fabs(mAinOld[n])+fabs(mAinNew[n])) && EXTENDED_ERRORLOG ) 
					{
					errorlog_file << "CellSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << endl;
					errorlog_file << "mAin: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
					errorlog_file << "\tdifferenza : " << 100.*(mAinNew[n]-mAinOld[n])/mAinOld[n] << "%\n";
					errorlog_file << "\tfase : " << CellVector[n].Get_phase() << endl;
					errorlog_file << "\teta' di fase: " << CellVector[n].Get_phase_age() << " s\n" << endl;
					}
				}
			if( fabs(mAextOld[n] - mAextNew[n]) > eps*0.5*(fabs(mAextOld[n])+fabs(mAextNew[n]))+TOL )
				{
				isOK = false;
				cellisOK = false;
				convergence_fail[8]++;
				long double newprec = 2.*fabs(mAextOld[n] - mAextNew[n])/(fabs(mAextOld[n])+fabs(mAextNew[n]));
				prec = newprec > prec ? newprec : prec;
				if ( fabs(mAextOld[n] - mAextNew[n]) > 0.5*(fabs(mAextOld[n])+fabs(mAextNew[n])) && EXTENDED_ERRORLOG ) 
					{
					errorlog_file << "CellSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << endl;
					errorlog_file << "mAext: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
					errorlog_file << "\tdifferenza : " << 100.*(mAextNew[n]-mAextOld[n])/mAextOld[n] << "%\n";
					errorlog_file << "\tfase : " << CellVector[n].Get_phase() << endl;
					errorlog_file << "\teta' di fase: " << CellVector[n].Get_phase_age() << " s\n" << endl;
					}
				}
			if( fabs(mAcLinOld[n] - mAcLinNew[n]) > eps*0.5*(fabs(mAcLinOld[n])+fabs(mAcLinNew[n]))+TOL )
				{
				isOK = false;
				cellisOK = false;
				convergence_fail[9]++;
				long double newprec = 2.*fabs(mAcLinOld[n] - mAcLinNew[n])/(fabs(mAcLinOld[n])+fabs(mAcLinNew[n]));
				prec = newprec > prec ? newprec : prec;
				if ( fabs(mAcLinOld[n] - mAcLinNew[n]) > 0.5*fabs(fabs(mAcLinOld[n])+fabs(mAcLinNew[n])) && EXTENDED_ERRORLOG ) 
					{
					errorlog_file << "CellSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << endl;
					errorlog_file << "mAcLin: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
					errorlog_file << "\tdifferenza : " << 100.*(mAcLinNew[n]-mAcLinOld[n])/mAcLinOld[n] << "%\n";
					errorlog_file << "\tfase : " << CellVector[n].Get_phase() << endl;
					errorlog_file << "\teta' di fase: " << CellVector[n].Get_phase_age() << " s\n" << endl;
					}
				}
			if( fabs(mAcLextOld[n] - mAcLextNew[n]) > eps*0.5*(fabs(mAcLextOld[n])+fabs(mAcLextNew[n]))+TOL )
				{
				isOK = false;
				cellisOK = false;
				convergence_fail[10]++;
				long double newprec = 2.*fabs(mAcLextOld[n] - mAcLextNew[n])/(fabs(mAcLextOld[n])+fabs(mAcLextNew[n]));
				prec = newprec > prec ? newprec : prec;
				if ( fabs(mAcLextOld[n] - mAcLextNew[n]) > 0.5*fabs(fabs(mAcLextOld[n])+fabs(mAcLextNew[n])) && EXTENDED_ERRORLOG ) 
					{
					errorlog_file << "CellSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << endl;
					errorlog_file << "mAcLext: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
					errorlog_file << "\tdifferenza : " << 100.*(mAcLextNew[n]-mAcLextOld[n])/mAcLextOld[n] << "%\n";
					errorlog_file << "\tfase : " << CellVector[n].Get_phase() << endl;
					errorlog_file << "\teta' di fase: " << CellVector[n].Get_phase_age() << " s\n" << endl;
					}
				}
			if( fabs(ATPpOld[n] - ATPpNew[n]) > eps*0.5*(fabs(ATPpOld[n])+fabs(ATPpNew[n]))+TOL )
				{
				isOK = false;
				cellisOK = false;
				convergence_fail[11]++;
				long double newprec = 2.*fabs(ATPpOld[n] - ATPpNew[n])/(fabs(ATPpOld[n])+fabs(ATPpNew[n]));
				prec = newprec > prec ? newprec : prec;
				if ( fabs(ATPpOld[n] - ATPpNew[n]) > 0.5*(fabs(ATPpOld[n])+fabs(ATPpNew[n])) && EXTENDED_ERRORLOG ) 
					{
					errorlog_file << "CellSystem::Diff() - passo " << Get_nstep() << ", iterazione " << nrepeats << endl;
					errorlog_file << "ATPp: Differenza maggiore del 50% nel calcolo iterativo per la cellula " << n <<"-esima\n";
					errorlog_file << "\tdifferenza : " << 100.*(ATPpNew[n]-ATPpOld[n])/ATPpOld[n] << "%\n";
					errorlog_file << "\tfase : " << CellVector[n].Get_phase() << endl;
					errorlog_file << "\teta' di fase: " << CellVector[n].Get_phase_age() << " s\n" << endl;
					}
				}
			
			if( !cellisOK ) ncell_fails++;
			
			// qui si rimpiazzano le vecchie approssimazioni con le nuove
			volumeOld[n] = volumeNew[n];
			volume_extraOld[n] = volume_extraNew[n];
			MitOld[n] = MitNew[n];
			pHiOld[n] = pHiNew[n];
			pHOld[n] = pHNew[n];
			mGinOld[n] = mGinNew[n];
			mGextOld[n] = mGextNew[n];
			mG6POld[n] = mG6PNew[n];
			mO2Old[n] = mO2New[n];
			StoreOld[n] = StoreNew[n];
			mAinOld[n] = mAinNew[n];
			mAextOld[n] = mAextNew[n];
			mAcLinOld[n] = mAcLinNew[n];
			mAcLextOld[n] = mAcLextNew[n];
			ATPpOld[n] = ATPpNew[n];
			}
		
		}
		

	nrepeats++;
	if( nrepeats >= MAXREPEATS )
		{
		errorlog_file << "\n*** ATTENZIONE: al passo " << nstep << " potenziale problema di convergenza dell'algoritmo in CellSystem::Diff ***" << endl;
		errorlog_file << "precisione limitata a " << scientific << prec << " per " << ncell_fails << " cellule" << "\n" << endl;
		}
	
	} while( !isOK && nrepeats < MAXREPEATS );
	// fine del loop per la soluzione delle equazioni



	// calcolo del flusso di AcL nell'ambiente
	AcLFlow = ( mAcL_envOld - Env.GetEnvironmentAcL() )/dt;

	// assegnazione delle variabili ambientali
	Env.SetEnvironmentvolume(volume_envOld);
	Env.SetEnvironmentG(mG_envOld);
	Env.SetEnvironmentO2(mO2_envOld);
	Env.SetEnvironmentA(mA_envOld);
	Env.SetEnvironmentAcL(mAcL_envOld);
	Env.SetEnvironmentpH();
	


	for(unsigned long n=0; n<ncell; n++)						// assegnazioni finali
		{
		
		// assegnazioni comuni a cellule vive e a cellule morte: si tratta delle variabili relative alla diffusione che continua anche 
		// nel caso delle cellule morte
		CellVector[n].Set_pH(pHOld[n]);
		CellVector[n].Set_G_extra(mGextOld[n]);
		CellVector[n].Set_O2(mO2Old[n]);
		CellVector[n].Set_A_extra(mAextOld[n]);
		CellVector[n].Set_AcL_extra(mAcLextOld[n]);
		
		// assegnazioni che si fanno solo nel caso in cui la cellula e' viva
		if( CellVector[n].Get_phase() != dead )					
			{

			CellVector[n].Set_pHi(pHiOld[n]);
			CellVector[n].Set_G(mGinOld[n]);
			CellVector[n].Set_G6P(mG6POld[n]);
			CellVector[n].Set_store(StoreOld[n]);
			CellVector[n].Set_A(mAinOld[n]);
			CellVector[n].Set_AcL(mAcLinOld[n]);

			
			// ******* ATTENZIONE !!! ******* 
			// questa istruzione definisce e memorizza anche raggio, superficie, volume,vol_extra, massa e ATPmin
			// nel caso delle cellule morte si usa Set_volume (v. sotto)
			// ATTENZIONE: l'istruzione di setting del DNA va messa PRIMA di quella del volume, etc. perche' il valore del 
			// DNA calcolato viene utilizzato per il setting del volume, altrimenti non c'e' consistenza ... 
			CellVector[n].Set_DNA(DNA[n]);
			CellVector[n].Set_M_and_ATPp(MitOld[n], ATPpOld[n]);

			
			// altre variabili e quantita' collegate
			CellVector[n].Set_protein(protein[n]);
			CellVector[n].Set_prot_rate(delta_protein[n]/dt);
			CellVector[n].Set_DNA_rate(DNA_rate[n]);
						
			CellVector[n].Set_pRb(pRb[n]);
			
			if( CellVector[n].Get_phase() == G1m_phase)
				{
				long double cyclinD = CellVector[n].Get_cyclinD() + CellVector[n].Get_type()->Get_cyclinD_fraction()*delta_protein[n];
				CellVector[n].Set_cyclinD( cyclinD );
				}
			if( CellVector[n].Get_phase() == G1p_phase)
				{
				long double cyclinE = CellVector[n].Get_cyclinE() + CellVector[n].Get_type()->Get_cyclinE_fraction()*delta_protein[n];
				CellVector[n].Set_cyclinE( cyclinE );
				}
			if( CellVector[n].Get_phase() == G2_phase)
				{
				long double cyclinX = CellVector[n].Get_cyclinX() + CellVector[n].Get_type()->Get_cyclinX_fraction()*delta_protein[n];
				CellVector[n].Set_cyclinX( cyclinX);
				}
			
			CellVector[n].Set_concS(ConcS[n]);
						
			// variabili ausiliarie 
			CellVector[n].Set_GAbsRate(GAbsRate[n]);
			CellVector[n].Set_GConsRate(GConsRate[n]);
			CellVector[n].Set_AAbsRate(AAbsRate[n]);
			CellVector[n].Set_AConsRate(AConsRate[n]);
			long double newATPprod = CellVector[n].Get_ATPprod() + (ATP_Ox[n] + ATP_NOx[n] + ATP2[n] + ATP3[n])*dt;
			CellVector[n].Set_ATPprod(newATPprod);
			long double newATPcons = CellVector[n].Get_ATPcons() + (ConsATP[n] + ConsATP_1[n] + ConsATP_2[n] + ConsATP_3[n] + ConsATP_5[n])*dt;
			CellVector[n].Set_ATPcons(newATPcons);
			
			CellVector[n].Set_ATPtot( (ATP_Ox[n] + ATP_NOx[n] + ATP2[n] + ATP3[n]) - (ConsATP[n] + ConsATP_1[n] + ConsATP_2[n] + ConsATP_3[n] + ConsATP_5[n]) );

			CellVector[n].Set_ATP_Ox(ATP_Ox[n]);
			CellVector[n].Set_ATP_NOx(ATP_NOx[n]);
			CellVector[n].Set_ATP2(ATP2[n]);
			CellVector[n].Set_ATP3(ATP3[n]);
			CellVector[n].Set_ConsATP(ConsATP[n]);
			CellVector[n].Set_ConsATP_1(ConsATP_1[n]);
			CellVector[n].Set_ConsATP_2(ConsATP_2[n]);
			CellVector[n].Set_ConsATP_3(ConsATP_3[n]);
			CellVector[n].Set_ConsATP_5(ConsATP_5[n]);


			CellVector[n].Set_StoreFillRate(StoreFillRate[n]);
			CellVector[n].Set_StoreConsRate(StoreConsRate[n]);
			CellVector[n].Set_AcLRate(AcLRate[n]);
			CellVector[n].Set_AcLOutRate(AcLOutRate[n]);

			}
		else	
		// la cellula e' morta e ci sono solo processi residuali 												
			{
			
			// il volume cambia solo se la cellula non si è ancora ridotta al solo nucleo
			if( volumeOld[n] > CellVector[n].Get_type()->Get_Vmin() )
				{
				// questo passo corrisponde alla soluzione esatta dell'equazione per il volume
				long double new_volume = volumeOld[n]*exp(-dt*CellVector[n].Get_type()->Get_DVap());
				CellVector[n].Set_volume(new_volume);		// questa istruzione aggiorna anche raggio, superficie, vol_extra e massa
				}
			
			}
			
		// controllo finale di consistenza
		if( CellVector[n].Get_phase() != dead ) 
			{
			int code = CellVector[n].CheckMVA();
			if(code < 0) errorlog_file << "Errore " << code << " alla fine di CellSystem::Diff nel controllo di consistenza per la cellula " << n << "\n" << endl;
			}
			
		
		}


// cleanup finale

	delete[] volumeOld;
	delete[] volumeNew;
	delete[] volume_extraOld;
	delete[] volume_extraNew;

	delete[] MitOld;
	delete[] MitNew;
	
	delete[] pHiOld;
	delete[] pHiNew;
	delete[] pHOld;
	delete[] pHNew;
	
	delete[] mGinOld;
	delete[] mGinNew;
	delete[] mGextOld;
	delete[] mGextNew;
	
	delete[] mG6POld;
	delete[] mG6PNew;
	
	delete[] mO2Old;
	delete[] mO2New;
	
	delete[] StoreOld;
	delete[] StoreNew;
	
	delete[] mAinOld;
	delete[] mAinNew;
	delete[] mAextOld;
	delete[] mAextNew;
	
	delete[] mAcLinOld;
	delete[] mAcLinNew;
	delete[] mAcLextOld;
	delete[] mAcLextNew;

	delete[] ATPpOld;
	delete[] ATPpNew;
	
	delete[] protein;
	delete[] pRb;
	delete[] delta_protein;
	delete[] ConcS;
	
	delete[] DNA;
	delete[] DNA_rate;
	
	delete[] GAbsRate;
	delete[] GConsRate;
	delete[] AAbsRate;
	delete[] AConsRate;

	delete[] StoreFillRate;
	delete[] StoreConsRate;
	delete[] AcLRate;
	delete[] AcLOutRate;
	
	delete[] ATP_Ox;
	delete[] ATP_NOx;
	delete[] ATP2;
	delete[] ATP3;
	delete[] ConsATP;
	delete[] ConsATP_1;
	delete[] ConsATP_2;
	delete[] ConsATP_3;
	delete[] ConsATP_5;
	
	
}


