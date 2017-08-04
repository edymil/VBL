/*
 *  CellsSystem-A.cpp
 *  Sim3D
 *
 *  Created by Edoardo Milotti on 20/04/10.
 *  Copyright 2010 I.N.F.N.-Sezione di Trieste. All rights reserved.
 * 
 *  This file is for initialization of stuff.
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

#include <sys/utsname.h>	// header per i metodi di identificazione della macchina

extern double EventTime;
extern bool eventON;
extern double pAlt;

//  ******************** InitializeCellsSystem ********************

// ***************************************************************
// Run initialization from terminal
//
// NB: this method is NOT called in case of continuation of run ...
// ***************************************************************
//
void CellsSystem::InitializeCellsSystem( bool terminal )
{

	unsigned long n;
	int ft;
	
	ifstream commands;

// se la lettura viene fatta da file si apre il file di comandi
	if( !terminal )
		commands.open( commandFile.c_str() );
		
// inizializzazione del seme del generatore di numeri casuali 

	cout << "The seed of the generator and' " << idum << endl;

	
// queste dichiarazioni servono alla definizione di flusso e dose di radiazione: si tratta di una struttura un po' ridondante che 
// tengo solo per chiarezza	
	EnvironmentalSignalType flow_type;
	double flow_min;
	double flow_max;
	double flow_period;
	double flow_start;
	double flow_stop;
	double flow_tON;
	double flow_tOFF;

	EnvironmentalSignalType dose_rate_type;
	double dose_rate_min;
	double dose_rate_max;
	double dose_rate_period;
	double dose_rate_start;
	double dose_rate_stop;
	double dose_rate_tON;
	double dose_rate_tOFF;
// fine delle dichiarazioni per flusso e dose


// la prima stampa dopo l'inizializzazione produce anche un header
	// first_print = true;
	// first_print2file = true;
	
// setup della simulazione
	if( terminal )
		{
		cout << endl;

		// cout << "Simulazione con cellule disperse (0), full 3D (1), oppure con configurazione fissa (2) ? ";
		cout << "Simulation with dispersed cells (0) or full 3D (1) ? ";
		cin >> sim_type;
		
		}
	else 
		{
		commands >> sim_type;
		}

	switch(sim_type)
		{
		case Disperse:
		cout << "E' stata selezionata la simulazione con cellule disperse\n" << endl;
		break;
		
		case Full3D:
		cout << "E' stata selezionata la simulazione full 3D\n" << endl;
		maxdr = 0.;
		break;
		
		case FixedConfig:
		//cout << "E' stata selezionata la simulazione con configurazione fissa\n" << endl;
		break;
						
		default:
		break; 		 
		}

	if( sim_type == Disperse || sim_type == Full3D )
		{
		if( terminal )
			{
			cout << "Numero di cellule iniziali: ";
			cin >> nstart;
			}
		else 
			{
			commands >> nstart;
			cout << "Numero di cellule iniziali: " << nstart << endl;
			}
		}

// scelta del tipo di distribuzione spaziale

// se le cellule sono disperse, vengono sparse nell'intero volume	
	if( sim_type == Disperse )
		{
		initial_cell_dist = Dist3D;
		}

// se la simulazione e' di tipo 3D allora si distinguono diversi tipi di simulazione
	if( sim_type == Full3D )
		{
		// se c'e' solo una cellula, allora non ci sono scelte possibili ... 
		if( nstart == 1 )
			{
			initial_cell_dist = OneCell;
			}
		// se ci sono piu' cellule allora si prende la distribuzione 2D oppure quella 3D
		else 
			{
			if( terminal )
				{
				cout << "Tipo di distribuzione delle cellule (2 = distribuzione 2D, 3 = distribuzione 3D): ";
				cin >> initial_cell_dist;
				}
			else 
				{
				commands >> initial_cell_dist;
				if(initial_cell_dist == 2) 
					cout << "E' stata selezionata la distribuzione di cellule 2D" << endl;
				else if (initial_cell_dist == 3)
					cout << "E' stata selezionata la distribuzione di cellule 3D" << endl;
				else 
					cout << "ATTENZIONE: distribuzione di cellule indefinita!" << endl;
					}
			}
		}


	if( terminal )
		{
		cout << "Stepsize (s): ";
		cin >> dt;
		cout << "Precisione per il metodo di diffusione trasporto e metabolismo (-1 = default = 1e-5): ";
		cin >> eps;
		if (eps < 0) eps = 1e-5; // precisione relativa
		cout << "Precisione nella determinazione della velocita' (micron/s) (-1 = default = 1e-5 micron/s): ";
		cin >> delta_vmax;
		if( delta_vmax < 0 ) delta_vmax = 1.e-5; // micron/s
		cout << "Durata dell'inizializzazione (s) (-1 = default = 1e6 s): ";
		cin >> t_ini;
		if(t_ini < 0) t_ini = 1000000.; // s
		cout << "Durata della simulazione (s): ";
		cin >> tmax;
		cout << "Maximum CPU time for a fraction of run (s) (Select if you break the run in multiple parts, if <= 0 = idle): ";
		cin >> t_CPU_max;
        cout << "Stepsize per slow motion (s): (0 = no slow motion) ";
        cin >> dt_sm;
        if(dt_sm > 0)
            {
            cout << "Slow motion start time (s): ";
            cin >> tsm_start;
            cout << "Slow motion stop time (s): ";
            cin >> tsm_stop;
            }
        else
            {
            tsm_start = 0;
            tsm_stop = 0;
            }
		}
	else 
		{
		commands >> dt;
		cout << "Stepsize (s): " << dt << endl;
		commands >> eps;
		if (eps < 0) eps = 1e-2; // precisione relativa
		cout << "Precisione per il metodo di diffusione trasporto e metabolismo: " << eps << endl;
		commands >> delta_vmax;
		if( delta_vmax < 0 ) delta_vmax = 1.e-4; // micron/s
		cout << "Precisione nella determinazione della velocita' (micron/s): " << delta_vmax << endl;
		commands >> t_ini;
		if(t_ini < 0) t_ini = 1000000.; // s
		cout << "Durata dell'inizializzazione (s): " << t_ini << endl;
		commands >> tmax;
		cout << "Durata della simulazione (s): " << tmax << endl;
		commands >> t_CPU_max;
		cout << "Massimo CPU time per una frazione di run (s): " << t_CPU_max << endl;
        commands >> dt_sm;
        cout << "Stepsize per slow motion (s): " << dt_sm << endl;
        if(dt_sm > 0)
            {
            commands >> tsm_start;
            cout << "Slow motion start time (s): " << tsm_start << endl;
            commands >> tsm_stop;
            cout << "Slow motion stop time (s): " << tsm_stop << endl;
            }
        else
            {
            tsm_start = 0;
            tsm_stop = 0;
            }
		}

	nmax = (unsigned long)floor((tmax+t_ini)/dt);

/*		
	cout << "Selezione dell'ambiente (0 = Standard, 1 = Modificato): ";
	
	int EnvType;
	cin >> EnvType;
	
	Environment input_env((EnvironmentTypeSelector)EnvType);		// qui si definisce l'ambiente
*/
	// if( terminal )
		{
		cout << "Ambiente letto dal file " << Get_EnvironmentFile() << "\n" << endl;
		}
		
	Environment input_env( EnvironmentFile );		// qui si definisce l'ambiente
	
	Env_0 = input_env;						// l'ambiente definito localmente viene copiato nella variabile di CellsSystem
	Env = Env_0;							// l'ambiente attuale e' uguale all'ambiente iniziale


// il codice seguente serve a defire le caratteristiche del flusso nell'ambiente

	if( terminal )
		{
		cout << "\nFlusso nullo (0) o non nullo (1)? ";
		cin >> flowON;
		}
	else 
		{
		commands >> flowON;
		}
	if(!flowON) cout << "In questo run il flusso e' nullo" << endl;

	if(!flowON)
		{
		flow_type = NullSignal;
		flowSignal = EnvironmentalSignal( );
    oxygenflowON = 0;
		}
	else
		{
		if( terminal )
			{
			cout << "Tipo di segnale (0 = costante, 1 = sinusoidale, 2 = onda quadra, 3 = impulso singolo): ";
			cin >> ft;
      cout << "Si modula anche l'ossigeno? (0 = NO, 1 = SI) ";
      cin >> oxygenflowON;
			}
		else 
			{
			commands >> ft;
			switch (ft) 
				{
				case 0:
					cout << "Segnale costante" << endl;
					break;
				case 1:
					cout << "Segnale sinusoidale" << endl;
					break;				
				case 2:
					cout << "Onda quadra" << endl;
					break;
				case 3:
					cout << "Impulso singolo" << endl;
					break;
				default:
					break;
				}
            commands >> oxygenflowON;
			}

		flow_type = (EnvironmentalSignalType)ft;
		
		switch(flow_type)
			{
			case ConstantSignal:
			if( terminal )
				{
				cout << "ampiezza del flusso costante (microlitri/min): ";
				cin >> flow_min;
				}
			else 
				{
				commands >> flow_min;
				cout << "ampiezza del flusso costante (microlitri/min): " << flow_min <<endl;
				}

			flow_min *= (1.e9)/60.;			// conversione da microlitri/min a micron^3/s
			flowSignal = EnvironmentalSignal( flow_min );
			break;
			
			case SineSignal:
			if( terminal )
				{
				cout << "flusso minimo (microlitri/min): ";
				cin >> flow_min;
				cout << "flusso massimo (microlitri/min): ";
				cin >> flow_max;
				cout << "Istante di inizio della modulazione (s): ";
				cin >> flow_start;
				cout << "Periodo di modulazione (s): ";
				cin >> flow_period;
				}
			else 
				{
				commands >> flow_min;
				cout << "flusso minimo (microlitri/min): " << flow_min << endl;
				commands >> flow_max;
				cout << "flusso massimo (microlitri/min): " << flow_max << endl;
				commands >> flow_start;
				cout << "Istante di inizio della modulazione (s): " << flow_start << endl;
				commands >> flow_period;
				cout << "Periodo di modulazione (s): " << flow_period << endl;
				}

			flow_min *= (1.e9)/60.;			// conversione da microlitri/min a micron^3/s
			flow_max *= (1.e9)/60.;			// conversione da microlitri/min a micron^3/s
			flowSignal = EnvironmentalSignal( flow_type, flow_min, flow_max, flow_start, flow_period );
			break;
					
			case SquareSignal:
			if( terminal )
				{
				cout << "flusso minimo (microlitri/min): ";
				cin >> flow_min;
				cout << "flusso massimo (microlitri/min): ";
				cin >> flow_max;
				cout << "Istante di inizio della modulazione (s): ";
				cin >> flow_start;
				cout << "Durata della fase ON (s): ";
				cin >> flow_tON;
				cout << "Durata della fase OFF (s): ";
				cin >> flow_tOFF;
				}
			else 
				{
				commands >> flow_min;
				cout << "flusso minimo (microlitri/min): " << flow_min << endl;
				commands >> flow_max;
				cout << "flusso massimo (microlitri/min): " << flow_max << endl;
				commands >> flow_start;
				cout << "Istante di inizio della modulazione (s): " << flow_start << endl;
				commands >> flow_tON;
				cout << "Durata della fase ON (s): " << flow_tON << endl;
				commands >> flow_tOFF;
				cout << "Durata della fase OFF (s): " << flow_tOFF << endl;
				}

			flow_min *= (1.e9)/60.;			// conversione da microlitri/min a micron^3/s
			flow_max *= (1.e9)/60.;			// conversione da microlitri/min a micron^3/s
			flowSignal = EnvironmentalSignal( flow_type, flow_min, flow_max, flow_start, flow_tON, flow_tOFF );
			break;
			
			case Pulse:
			if( terminal )
				{
				cout << "flusso minimo (microlitri/min): ";
				cin >> flow_min;
				cout << "flusso massimo (microlitri/min): ";
				cin >> flow_max;
				cout << "Istante di inizio dell'impulso (s): ";
				cin >> flow_start;
				cout << "Istante di fine dell'impulso (s): ";
				cin >> flow_stop;
				}
			else 
				{
				commands >> flow_min;
				cout << "flusso minimo (microlitri/min): " << flow_min << endl;
				commands >> flow_max;
				cout << "flusso massimo (microlitri/min): " << flow_max << endl;
				commands >> flow_start;
				cout << "Istante di inizio dell'impulso (s): " << flow_start << endl;
				commands >> flow_stop;
				cout << "Istante di fine dell'impulso (s): " << flow_stop << endl;
				}

			flow_min *= (1.e9)/60.;			// conversione da microlitri/min a micron^3/s
			flow_max *= (1.e9)/60.;			// conversione da microlitri/min a micron^3/s
			flowSignal = EnvironmentalSignal( flow_type, flow_min, flow_max, flow_start, flow_stop );
			break;
			
			case UserDefined:
			break;
			
			default:
			break; 
			 
			}
		}

// la sezione che segue è del tutto analoga alla precedente e serve a definire la dose di radiazione somministrata alle cellule	
	if( terminal )
		{
		cout << "\nDose di radiazione nulla (0) o non nulla (1)? ";
		cin >> doseON;
		}
	else 
		{
		commands >> doseON;
		}
	if(!doseON) cout << "In questo run la dose e' nulla" << endl;

	if(!doseON)
		{
		dose_rate_type = NullSignal;
		dose_rateSignal = EnvironmentalSignal( );
		}
	else
		{
		if( terminal )
			{
			cout << "Tipo di segnale (0 = costante, 1 = sinusoidale, 2 = onda quadra, 3 = impulso singolo): ";
			cin >> ft;
			}
		else 
			{
			commands >> ft;
			switch (ft) 
				{
				case 0:
					cout << "Segnale costante" << endl;
					break;
				case 1:
					cout << "Segnale sinusoidale" << endl;
					break;				
				case 2:
					cout << "Onda quadra" << endl;
					break;
				case 3:
					cout << "Impulso singolo" << endl;
					break;
				default:
					break;
				}
			}

		dose_rate_type = (EnvironmentalSignalType)ft;
		
		switch(dose_rate_type)
			{
			case ConstantSignal:
			if( terminal )
				{
				cout << "ampiezza della dose-rate costante (Gy/s): ";
				cin >> dose_rate_min;
				}
			else 
				{
				commands >> dose_rate_min;
				cout << "ampiezza della dose-rate costante (Gy/s): " << dose_rate_min << endl;
				}

			dose_rateSignal = EnvironmentalSignal( dose_rate_min );
			break;
			
			case SineSignal:
			if( terminal )
				{
				cout << "dose-rate minima (Gy/s): ";
				cin >> dose_rate_min;
				cout << "dose-rate massima (Gy/s): ";
				cin >> dose_rate_max;
				cout << "Istante di inizio della modulazione (s): ";
				cin >> dose_rate_start;
				cout << "Periodo di modulazione (s): ";
				cin >> dose_rate_period;
				}
			else 
				{
				commands >> dose_rate_min;
				cout << "dose-rate minima (Gy/s): " << dose_rate_min << endl;
				commands >> dose_rate_max;
				cout << "dose-rate massima (Gy/s): " << dose_rate_max << endl;
				commands >> dose_rate_start;
				cout << "Istante di inizio della modulazione (s): " << dose_rate_start << endl;
				commands >> dose_rate_period;
				cout << "Periodo di modulazione (s): " << dose_rate_period << endl;
				}

			dose_rateSignal = EnvironmentalSignal( dose_rate_type, dose_rate_min, dose_rate_max, dose_rate_start, dose_rate_period );
			break;
					
			case SquareSignal:
			if( terminal )
				{
				cout << "dose-rate minima (Gy/s): ";
				cin >> dose_rate_min;
				cout << "dose-rate massima (Gy/s): ";
				cin >> dose_rate_max;
				cout << "Istante di inizio della modulazione (s): ";
				cin >> dose_rate_start;
				cout << "Durata della fase ON (s): ";
				cin >> dose_rate_tON;
				cout << "Durata della fase OFF (s): ";
				cin >> dose_rate_tOFF;
				}
			else 
				{
				commands >> dose_rate_min;
				cout << "dose-rate minima (Gy/s): " << dose_rate_min  << endl;
				commands >> dose_rate_max;
				cout << "dose-rate massima (Gy/s): " << dose_rate_max  << endl;
				commands >> dose_rate_start;
				cout << "Istante di inizio della modulazione (s): " << dose_rate_start  << endl;
				commands >> dose_rate_tON;
				cout << "Durata della fase ON (s): " << dose_rate_tON  << endl;
				commands >> dose_rate_tOFF;
				cout << "Durata della fase OFF (s): " << dose_rate_tOFF  << endl;
				}

			dose_rateSignal = EnvironmentalSignal( dose_rate_type, dose_rate_min, dose_rate_max, dose_rate_start, dose_rate_tON, dose_rate_tOFF );
			break;
			
			case Pulse:
			if( terminal )
				{
				cout << "dose-rate minima (Gy/s): ";
				cin >> dose_rate_min;
				cout << "dose-rate massima (Gy/s): ";
				cin >> dose_rate_max;
				cout << "Istante di inizio dell'impulso (s): ";
				cin >> dose_rate_start;
				cout << "Istante di fine dell'impulso (s): ";
				cin >> dose_rate_stop;
				}
			else 
				{
				commands >> dose_rate_min;
				cout << "dose-rate minima (Gy/s): " << dose_rate_min << endl;
				commands >> dose_rate_max;
				cout << "dose-rate massima (Gy/s): " << dose_rate_max << endl;
				commands >> dose_rate_start;
				cout << "Istante di inizio dell'impulso (s): " << dose_rate_start << endl;
				commands >> dose_rate_stop;
				cout << "Istante di fine dell'impulso (s): " << dose_rate_stop << endl;
				}

			dose_rateSignal = EnvironmentalSignal( dose_rate_type, dose_rate_min, dose_rate_max, dose_rate_start, dose_rate_stop );
			break;
			
			case UserDefined:
			break;
			
			default:
			break; 
			 
			}
		}

	
	
	//CellType type;							// Here is defined the standard phenotype
	CellType initial_type( CellTypeFile );		// Here is defined the file phenotype
  initial_type.Set_name( 0 );                 // Here is the name of the standard type
  CellType alt_type( CellTypeFileAlt );		// Here is the alternative phenotype of file
  alt_type.Set_name( 1 );                     // Here is the name of the alternative type
    
	
	// cout << "\nTipo iniziale: \n" << endl;	// stampa del tipo iniziale
	// cout << initial_type;
	
	CellTypeVector.push_back( initial_type );	// si immagazzina il tipo appena definito nel vettore dei tipi cellulari
                                                // per assicurarne la permanenza
												// si noti che questo primo tipo cellulare sta in posizione 0 nel vettore
    ntypes = 1;                                 // qui si inizializza ntypes (conteggio dei tipi già definiti)
	
	CellTypeVector.push_back( alt_type );       // si immagazzina il tipo alternativo in fondo al vettore, in posizione 1
    ntypes++;                                   // si incrementa ntypes
    
    lastname = 0;
	AddInitializedCell(idum, &CellTypeVector[0], &Env); 
	
	// la cellula viene immediatamente replicata (nstart-1 copie, in modo da avere un totale di nstart cellule)
	if(nstart > 1)
		ReplicateCell( 0 , nstart-1 );
	
	
	// if( terminal )
		{
		cout << "\nAmbiente iniziale: " << endl;
		cout << Env << endl;						// printout delle condizioni ambientali
		}
	
	double xmin_s = Env.Get_xmin();			// qui si recupera la forma del volume
	double xmax_s = Env.Get_xmax();	
	double ymin_s = Env.Get_ymin();
	double ymax_s = Env.Get_ymax();	
	double zmin_s = Env.Get_zmin();
	double zmax_s = Env.Get_zmax();	

	if( sim_type == Disperse || sim_type == Full3D )
		{
		
		double delta_x = 0, delta_y = 0, delta_z = 0.;	// semiintervalli per ciascuna coordinata
		double x_0, y_0, z_0;							// centro del parallelepipedo
		
		x_0 = 0.5*(xmax_s+xmin_s);
		y_0 = 0.5*(ymax_s+ymin_s);
		z_0 = 0.5*(ymax_s+ymin_s);
		
		switch (initial_cell_dist) 
			{
			case OneCell:
				delta_x = delta_y = delta_z = 0.;
				break;
			case Dist2D:
				delta_x = 0.5*(xmax_s-xmin_s);
				delta_y = 0.5*(ymax_s-ymin_s);
				delta_z = 0.01*(zmax_s-zmin_s);
				break;
			case Dist3D:
				delta_x = 0.5*(xmax_s-xmin_s);
				delta_y = 0.5*(ymax_s-ymin_s);
				delta_z = 0.5*(zmax_s-zmin_s);
				break;
			default:
				break;
			}

		for(n=0; n<nstart; n++)
			{
			x[n] = x_0 + 0.1*delta_x*(2.*ran2(idum)-1.);	// inizializzazione dei centri delle cellule
			y[n] = y_0 + 0.1*delta_y*(2.*ran2(idum)-1.);
			z[n] = z_0 + 0.1*delta_z*(2.*ran2(idum)-1.);
			vx[n] =  0.;									// inizializzazione delle velocita'
			vy[n] =  0.;
			vz[n] =  0.;
			neigh[n] =  0;									// inizialmente non ci sono vicini definiti		
			}
			
		}


// in questa posizione c'era il codice per la configurazione fissa; il codice non e' piu' mantenuto ed e' ora in un frammento a parte

	if( terminal )
		{
		cout << "\nNumer of dt: " << dt << "steps Seconds between two file records: ";
		cin >> nprint;
		}
	else 
		{
		commands >> nprint;
		cout << "\nNumero di passi da " << dt << " secondi tra due record su file: " << nprint << endl;
		}

	if( terminal )
		{
		cout << "\nNumero di passi da " << dt << " Number of 1-second steps between two on-screen recordings: ";
		cin >> nscreen;
		}
	else 
		{
		commands >> nscreen;
		cout << "\nNumero di passi da " << dt << " secondi tra due record su schermo: " << nscreen << endl;
		}
    
// gestione eventi speciali
    if( terminal )
        {
        cout << "\nEvento speciale? (NO = 0, SI = 1) ";
        cin >> eventON;
        }
    else
        {
        commands >> eventON;
        }
    
    if( eventON )
        {
        cout << "In questo run c'e' un evento speciale" << endl;
            
        if( terminal )
            {
            cout << "Epoca dell'evento speciale: ";
            cin >> EventTime;
            }
        else
            {
            commands >> EventTime;
            cout << "\nEpoca dell'evento speciale: " << EventTime << endl;
            }

        if( terminal )
            {
            cout << "probabilita' dell'evento speciale (2 = si verifica un sono evento): ";
            cin >> pAlt;
            }
        else
            {
            commands >> pAlt;
            cout << "\nprobabilita' dell'evento speciale (2 = si verifica un sono evento): " << pAlt << endl;
            }
        }
    else
        {
        EventTime = 0.;
        pAlt = 0;
        }

		
	ncells = nstart;						// numero di cellule alla partenza
	//ntypes = 1;								// un solo tipo cellulare, per ora
	
	t = treal = 0;							// inizializzazione dei timers del tempo di simulazione
	nstep = 0;								// inizializzazione del passo
    
    slow_motion = false;                    // inizialmente lo slow motion non c'e'

	ready2start = false;					// il programma non e' ancora pronto a partire
  faketumAtCurrentTime = false;
}


//  ******************** RunDefinition ********************

// questo metodo definisce i nomi dei files 
// 
void CellsSystem::RunDefinition(  ) 
{

	
	
/*

// lettura del nome del directory di lettura-scrittura
	ifstream filein( "./prefix.txt" );
	filein >> dir; 
	filein.close();
	
*/

// qui si fa l'update del numero del run	
	fstream runfile( "./runs/run_number.txt" );			// si cerca di aprire in lettura il file run_number.txt
	if( !runfile )	// se non riesce ... 
		{
		
		ofstream fileout( "./runs/run_number.txt" );	// ... allora si crea il file e ci si mette il numero 1
		run = 1;
		fileout << "name" << endl;
		fileout << run << endl;
		fileout.close();
		
		}
	else											// se l'apertura del file riesce ... 
		{
		
		runfile >> machine;							// si legge il nome della macchina
		runfile >> run;								// si legge l'ultimo numero di run
		runfile.close();							// si chiude il file in lettura ...
		run++;										// ... si incrementa il numero di run
		runfile.open( "./runs/run_number.txt");		// ... quindi si sovrascrive il vecchio file 
		runfile << machine << endl;					// con il nome della macchina
		runfile << run << endl;						// e con il nuovo numero di run
		runfile.close();
		}

// qui si definisce il nome del directory di output
	char str[100];
	sprintf(str,"./runs/run_%d-%s/",run,machine.c_str());
	//dir = dir+str;
	dir = str;
	cout << endl << "I/O directory: " << dir << "\n\n";

// qui si crea il directory di output
	string command;
	command = "mkdir ";
	command.append(dir);
	system( command.c_str() );
	
// qui si crea il directory per i files con i records dettagliati della configurazione
	command = "mkdir "+dir+"configuration/";
	system( command.c_str() ); 
	
// qui si crea il directory per i files con i records dettagliati della configurazione (binary write)
	command = "mkdir "+dir+"configuration_b/";
	system( command.c_str() ); 
	
// qui si crea il directory per i files con i records dettagliati del flusso (binary write)
	command = "mkdir "+dir+"flows_b/";
	system( command.c_str() ); 
	
	
	
// apertura del file stream per l'output su file delle impostazioni del run e quindi output
	string runfilename = dir+"run_data.txt";
	ofstream run_data( runfilename.c_str() );
	run_data << "Parametri del run " << run << endl << endl;
	
// printout dei parametri ottenuti da terminale o command file
	run_data << "Timestep (s): " << dt << endl;
	run_data << "Precisione per CellsSystem::Diff: " << eps << endl;
	run_data << "Precisione nella determinazione delle velocita (micron/s): " << delta_vmax << endl;
	run_data << "Durata inizializzazione (s): " << t_ini << endl;
	run_data << "Durata max della simulazione (s) " << tmax << endl;
	run_data << endl;


// printout dei dati ambientali
	run_data << "Ambiente iniziale: \n";
	run_data << Env_0 << endl;
	run_data << endl;
/*
	run_data << "Ambiente attuale: \n";
	run_data << CellsSystem->Env << endl;
	run_data << endl;
*/	

// printout di flusso e dose-rate
	if( flowON )
		{
		run_data << "Flusso -- non nullo -- con le seguenti caratteristiche (ampiezze in m^3/s):  " << endl;
		run_data << flowSignal << endl;
		}
	else
		run_data << "Flusso -- nullo -- " << endl;

	if( doseON )
		{
		run_data << "Dose di radiazione -- non nulla -- con le seguenti caratteristiche (ampiezze in Gy/s):  " << endl;
		run_data << dose_rateSignal << endl;
		}
	else
		run_data << "Dose di radiazione -- nulla --" << endl;

// printout dell'evento speciale 
       if( eventON )
      		{
        	run_data<< "In questo run c'e' un evento speciale" << endl;
                run_data<< "Epoca dell'evento speciale: " << EventTime << endl;
		run_data<< "\nprobabilita' dell'evento speciale: " << pAlt << endl;
                }
      else
                run_data<< "In questo run non c'e' alcun evento speciale" << endl;
          
	
// printout del tipo cellulare	
	run_data << "\nFenotipo cellulare iniziale (cellula 0): \n" << endl;
	run_data << *type[0] << endl;
	
// printout della prima cellula
	run_data << "\nCellula iniziale (cellula 0):\n" << endl;
	PrintCell( run_data, 0 );

// printout del numero di cellule iniziali		
	run_data << "\nNumero di cellule iniziali: " << nstart << endl;

// chiusura del file di logging del run	
	run_data.close();

// infine si restituiscono i nomi dei files su cui fare l'output dei dati
	output_filename = dir+"dump.txt";
	log_filename = dir+"logfile.txt";
	configuration_filename = dir+"configuration/";
	configuration_b_filename = dir+"configuration_b/";
	flow_b_filename = dir+"flows_b/";
	screen_dump_filename = dir+"screen_dump.txt";
	errorlog_filename = dir+"error_log.txt";
	cell_filename = dir+"cell_log.txt";
	env_filename = dir+"env_log.txt";
	convlog_filename = dir+"convergence_log.txt";
	
	cellsys_out_filename = dir+"cellsys_0.bin";

// apertura dei files di output
	output_file.open( output_filename.c_str() ); 
	log_file.open( log_filename.c_str() ); 
	screen_dump_file.open( screen_dump_filename.c_str() ); 
	errorlog_file.open( errorlog_filename.c_str() ); 
	cell_file.open( cell_filename.c_str() ); 
	env_file.open( env_filename.c_str() ); 
	convlog_file.open( convlog_filename.c_str() ); 
	
// inizializzazione del numero di configurazioni scritte
	
	part = 0;									// la parte attuale e' part = 0

	string partfilename = dir+"part.txt";
	ofstream partfile( partfilename.c_str() );	// qui si crea il file che contiene il numero di parte e si scrive la parte attuale (0)
	partfile << part << endl;
	partfile.close();

// stampa iniziale su error log
	
	errorlog_file << "Errors and warnings for run " << run << "\n" << endl;
	

}

// 
// versione overloaded che definisce i nomi dei files nel caso di continuazione di un run
// 
void CellsSystem::RunDefinition( string run_name ) 
{

/*

// lettura del nome del directory di lettura-scrittura
	ifstream filein( "./prefix.txt" );
	filein >> dir; 
	filein.close();
	
*/


// qui si definisce il nome del directory di output
	dir = "./runs/"+run_name+"/";
	cout << endl << "I/O directory: " << dir << "\n\n";
	
	
// qui si definisce la parte attuale
	string partfilename = dir+"part.txt";
	fstream partfile( partfilename.c_str() );	// qui si legge il file che contiene il numero di parte e dell'ultima configurazione scritta
	partfile >> part;
	partfile.close();
	
	char partname[100];						// partname e' la stringa corrispondente alla parte
	sprintf(partname,"%d",part);
	
	cellsys_in_filename = dir+"cellsys_"+partname+".bin";
	
	part++;										// qui si incrementa il numero di parte e si riscrive il file
	
	partfile.open( partfilename.c_str() );
	partfile << part << endl;
	partfile.close();
	
	sprintf(partname,"%d",part);				// stringa con il nuovo numero di parte


// infine si restituiscono i nomi dei files su cui fare l'output dei dati
	output_filename = dir+"dump_"+partname+".txt";
	log_filename = dir+"logfile_"+partname+".txt";
	configuration_filename = dir+"configuration/";
	configuration_b_filename = dir+"configuration_b/";
	flow_b_filename = dir+"flows_b/";
	screen_dump_filename = dir+"screen_dump_"+partname+".txt";
	errorlog_filename = dir+"error_log_"+partname+".txt";
	cell_filename = dir+"cell_log_"+partname+".txt";
	env_filename = dir+"env_log_"+partname+".txt";
	convlog_filename = dir+"convergence_log_"+partname+".txt";
	
	cellsys_out_filename = dir+"cellsys_"+partname+".bin";
	
	
// apertura dei files di output
	output_file.open( output_filename.c_str() ); 
	log_file.open( log_filename.c_str() ); 
	screen_dump_file.open( screen_dump_filename.c_str() ); 
	errorlog_file.open( errorlog_filename.c_str() ); 
	cell_file.open( cell_filename.c_str() ); 
	env_file.open( env_filename.c_str() ); 
	convlog_file.open( convlog_filename.c_str() ); 
	

// stampa iniziale su error log
	
	errorlog_file << "Errors and warnings for run " << run << "\n" << endl;
	

}


//  ******************** WriteCellsSystem ********************

// questo metodo scrive il CellsSystem per poterlo ricaricare e riutilizzare
// 
void CellsSystem::WriteCellsSystem( )
{

	ofstream stream;
	
 	stream.open( cellsys_out_filename.c_str(), ios::binary );

// parametri per la definizione del run	
	stream.write( (char*)(&sim_type), sizeof(int) );
	stream.write( (char*)(&run), sizeof(int) );
	stream.write( (char*)(&dt), sizeof(double) );
	stream.write( (char*)(&dt_sm), sizeof(double) );
	stream.write( (char*)(&t), sizeof(double) );
	stream.write( (char*)(&t_ini), sizeof(double) );
	stream.write( (char*)(&treal), sizeof(double) );
	stream.write( (char*)(&tmax), sizeof(double) );
	stream.write( (char*)(&tsm_start), sizeof(double) );
	stream.write( (char*)(&tsm_stop), sizeof(double) );
	stream.write( (char*)(&slow_motion), sizeof(bool) );
	
	stream.write( (char*)(&t_CPU_max), sizeof(double) );

	stream.write( (char*)(&nstep), sizeof(unsigned long) );
	stream.write( (char*)(&nstep_start), sizeof(unsigned long) );
	stream.write( (char*)(&nmax), sizeof(unsigned long) );
	stream.write( (char*)(&idum), sizeof(int) );

	stream.write( (char*)(&nprint), sizeof(unsigned long) );
	stream.write( (char*)(&nscreen), sizeof(unsigned long) );
	stream.write( (char*)(&nconfiguration), sizeof(unsigned long) );

	stream.write( (char*)(&eps), sizeof(double) );
	stream.write( (char*)(&delta_vmax), sizeof(double) );

	stream.write( (char*)(&ncells), sizeof(unsigned long) );
	stream.write( (char*)(&alive), sizeof(unsigned long) );
	stream.write( (char*)(&ntypes), sizeof(unsigned long) );

// parametri ambientali
	// ambiente iniziale Env_0
	Env_0.WriteEnvironment( stream );
	
	// ambiente attuale Env
	Env.WriteEnvironment( stream );

	// flussi
	stream.write( (char*)(&flowON), sizeof(bool) );
	stream.write( (char*)(&doseON), sizeof(bool) );
	
	// segnali
	flowSignal.WriteEnvironmentalSignal( stream );
	dose_rateSignal.WriteEnvironmentalSignal( stream );
	
	
// vettore dei tipi cellulari
	for(unsigned long int k=0; k<ntypes; k++)
		{
		CellTypeVector[k].WriteCellType( stream );
		}
	
	// maxdr e i dati sul flusso globale di O2 e AcL vengono ricalcolati 
	
	stream.write( (char*)(&name[0]) , ncells*sizeof( unsigned long ) );
	stream.write( (char*)(&mark[0]) , ncells*sizeof( int ) );
	for(unsigned long int n = 0; n<ncells; n++)
		{
		unsigned long int k = (unsigned long int)(type[n]-&CellTypeVector[0]);
		stream.write( (char*)(&k) , sizeof( unsigned long int ) );
		// cout << k << endl;
		}
	
	stream.write( (char*)(&Temperature[0]) , ncells*sizeof( double ) );
	
	stream.write( (char*)(&phase[0]) , ncells*sizeof( CellPhase ) );
	stream.write( (char*)(&death_condition[0]) , ncells*sizeof( int ) );
	stream.write( (char*)(&age[0]) , ncells*sizeof( float ) );
	stream.write( (char*)(&phase_age[0]) , ncells*sizeof( float ) );
	stream.write( (char*)(&age_mother[0]) , ncells*sizeof( float ) );
	stream.write( (char*)(&n_mitosis[0]) , ncells*sizeof( int ) );
	
	stream.write( (char*)(&x[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&y[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&z[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&vx[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&vy[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&vz[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&r[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&surface[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&volume[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&mass[0]) , ncells*sizeof( double ) );
	
	// le variabili dinamiche vxnew, fx, etc non vengono scritte ma ricalcolate dal metodo per la dinamica
	
	stream.write( (char*)(&volume_extra[0]) , ncells*sizeof( double ) );
	
	// le variabili geometriche che si calcolano con CGAL non vengono scritte ma ricalcolate al volo dopo la lettura

	stream.write( (char*)(&M[0]) , ncells*sizeof( double ) );
	
	stream.write( (char*)(&G[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&G6P[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&O2[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&store[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&A[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&AcL[0]) , ncells*sizeof( double ) );

	stream.write( (char*)(&h[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&pHi[0]) , ncells*sizeof( double ) );

	stream.write( (char*)(&protein[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&prot_rate[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&DNA[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&DNA_rate[0]) , ncells*sizeof( double ) );

	// i rates non vengono scritti ma ricalcolati dal metodo del metabolismo

	stream.write( (char*)(&ATP_St[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ATP_Ox[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ATP_NOx[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ATP2[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ATP3[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ConsATP[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ConsATP_1[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ConsATP_2[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ConsATP_3[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ConsATP_4[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ConsATP_5[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ATPtot[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ATPp[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ATPmin[0]) , ncells*sizeof( double ) );

	stream.write( (char*)(&ATPstart[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ATPprod[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ATPcons[0]) , ncells*sizeof( double ) );

	stream.write( (char*)(&G_extra[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&A_extra[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&AcL_extra[0]) , ncells*sizeof( double ) );

	stream.write( (char*)(&SensO2[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&ConsO[0]) , ncells*sizeof( double ) );

	stream.write( (char*)(&DNA_spread[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&M_T[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&pRb[0]) , ncells*sizeof( double ) );
	
	stream.write( (char*)(&ConcS[0]) , ncells*sizeof( double ) );
	
	stream.write( (char*)(&cyclinD[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&cyclinE[0]) , ncells*sizeof( double ) );
	stream.write( (char*)(&cyclinX[0]) , ncells*sizeof( double ) );
	
	stream.write( (char*)(&NpRbk[0]) , ncells*sizeof( double ) );
	
	stream.write( (char*)(&DNA_spread[0]) , ncells*sizeof( double ) );
	
	// le copie per i calcoli in Diff non vengono scritte su file
	
	stream.close();
		
}

//  ******************** ReadCellsSystem ********************

// questo metodo legge il CellsSystem per poterlo ricaricare e riutilizzare
// 
void CellsSystem::ReadCellsSystem( )
{ 

	ifstream stream;
	
	stream.open( cellsys_in_filename.c_str(), ios::binary );
	if( !stream.is_open() )
		{
		cout << "Impossibile aprire il file " << cellsys_in_filename << ", exiting ... " << endl;
		exit(-1);
		}

// parametri per la definizione del run	
	stream.read( (char*)(&sim_type), sizeof(int) );
	cout << "Tipo di simulazione: " << sim_type << endl;
	stream.read( (char*)(&run), sizeof(int) );
	cout << "Run: " << run << endl;
	stream.read( (char*)(&dt), sizeof(double) );
	cout << "dt (s): " << dt << endl;
	stream.read( (char*)(&dt_sm), sizeof(double) );
	cout << "dt_sm (s): " << dt_sm << endl;
	stream.read( (char*)(&t), sizeof(double) );
	stream.read( (char*)(&t_ini), sizeof(double) );
	stream.read( (char*)(&treal), sizeof(double) );
	cout << "tempo simulato raggiunto (s): " << treal << endl;
	stream.read( (char*)(&tmax), sizeof(double) );
	stream.read( (char*)(&tsm_start), sizeof(double) );
	stream.read( (char*)(&tsm_stop), sizeof(double) );
	stream.read( (char*)(&slow_motion), sizeof(bool) );

	stream.read( (char*)(&t_CPU_max), sizeof(double) );
	cout << "tempo di CPU max per una frazione di run (s): " << t_CPU_max << endl;

	stream.read( (char*)(&nstep), sizeof(unsigned long) );
	stream.read( (char*)(&nstep_start), sizeof(unsigned long) );
	stream.read( (char*)(&nmax), sizeof(unsigned long) );
	stream.read( (char*)(&idum), sizeof(int) );

	stream.read( (char*)(&nprint), sizeof(unsigned long) );
	stream.read( (char*)(&nscreen), sizeof(unsigned long) );
	stream.read( (char*)(&nconfiguration), sizeof(unsigned long) );

	stream.read( (char*)(&eps), sizeof(double) );
	stream.read( (char*)(&delta_vmax), sizeof(double) );

	ncells = 0;
	unsigned long old_ncells;
	stream.read( (char*)(&old_ncells), sizeof(unsigned long) );
	cout << "numero di cellule: " << old_ncells << endl;
	
	stream.read( (char*)(&alive), sizeof(unsigned long) );
	cout << "di cui vive: " << alive << endl;
	stream.read( (char*)(&ntypes), sizeof(unsigned long) );
	cout << "numero di fenotipi: " << ntypes << endl;

// parametri ambientali
	// ambiente iniziale Env_0
	Env_0.ReadEnvironment( stream );
	cout << "volume ambiente iniziale (microlitri) " << 1e-9*Env_0.volume << endl;
	
	// ambiente attuale Env
	Env.ReadEnvironment( stream );
	cout << "volume ambiente attuale (microlitri) " << 1e-9*Env.volume << endl;

	// flussi
	stream.read( (char*)(&flowON), sizeof(bool) );
	stream.read( (char*)(&doseON), sizeof(bool) );
	
	// segnali
	flowSignal.ReadEnvironmentalSignal( stream );
	dose_rateSignal.ReadEnvironmentalSignal( stream );
	
	cout << "flussi e segnali letti ... " << endl;
	
	
// vettore dei tipi cellulari

	CellTypeVector.resize( ntypes );
	
	for(unsigned long int k=0; k<ntypes; k++)
		{
		CellTypeVector[k].ReadCellType( stream );
		}
	
	cout << "lettura dei tipi cellulari completata " << endl;
	
	// maxdr e i dati sul flusso globale di O2 e AcL vengono ricalcolati 
	
	AddCells( old_ncells );		// qui si alloca spazio ai vettori che descrivono le cellule (da qui in poi e' definito ncells );
	
	stream.read( (char*)(&name[0]) , ncells*sizeof( unsigned long ) );
	stream.read( (char*)(&mark[0]) , ncells*sizeof( int ) );
	
	unsigned long int ktype;
	for(unsigned long int n = 0; n<ncells; n++)
		{
		stream.read( (char*)(&ktype) , sizeof( unsigned long int ) );
		type[n] = &(CellTypeVector[ktype]);
		}
	
	stream.read( (char*)(&Temperature[0]) , ncells*sizeof( double ) );
	
	stream.read( (char*)(&phase[0]) , ncells*sizeof( CellPhase ) );
	stream.read( (char*)(&death_condition[0]) , ncells*sizeof( int ) );
	stream.read( (char*)(&age[0]) , ncells*sizeof( float ) );
	stream.read( (char*)(&phase_age[0]) , ncells*sizeof( float ) );
	stream.read( (char*)(&age_mother[0]) , ncells*sizeof( float ) );
	stream.read( (char*)(&n_mitosis[0]) , ncells*sizeof( int ) );
	
	stream.read( (char*)(&x[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&y[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&z[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&vx[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&vy[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&vz[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&r[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&surface[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&volume[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&mass[0]) , ncells*sizeof( double ) );
	
	// le variabili dinamiche vxnew, fx, etc non vengono scritte ma ricalcolate dal metodo per la dinamica
	
	stream.read( (char*)(&volume_extra[0]) , ncells*sizeof( double ) );
	
	// le variabili geometriche che si calcolano con CGAL non vengono scritte ma ricalcolate al volo dopo la lettura

	stream.read( (char*)(&M[0]) , ncells*sizeof( double ) );
	
	stream.read( (char*)(&G[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&G6P[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&O2[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&store[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&A[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&AcL[0]) , ncells*sizeof( double ) );

	stream.read( (char*)(&h[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&pHi[0]) , ncells*sizeof( double ) );

	stream.read( (char*)(&protein[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&prot_rate[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&DNA[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&DNA_rate[0]) , ncells*sizeof( double ) );

	// i rates non vengono scritti ma ricalcolati dal metodo del metabolismo

	stream.read( (char*)(&ATP_St[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ATP_Ox[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ATP_NOx[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ATP2[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ATP3[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ConsATP[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ConsATP_1[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ConsATP_2[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ConsATP_3[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ConsATP_4[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ConsATP_5[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ATPtot[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ATPp[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ATPmin[0]) , ncells*sizeof( double ) );

	stream.read( (char*)(&ATPstart[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ATPprod[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ATPcons[0]) , ncells*sizeof( double ) );

	stream.read( (char*)(&G_extra[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&A_extra[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&AcL_extra[0]) , ncells*sizeof( double ) );

	stream.read( (char*)(&SensO2[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&ConsO[0]) , ncells*sizeof( double ) );

	stream.read( (char*)(&DNA_spread[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&M_T[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&pRb[0]) , ncells*sizeof( double ) );
	
	stream.read( (char*)(&ConcS[0]) , ncells*sizeof( double ) );
	
	stream.read( (char*)(&cyclinD[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&cyclinE[0]) , ncells*sizeof( double ) );
	stream.read( (char*)(&cyclinX[0]) , ncells*sizeof( double ) );
	
	stream.read( (char*)(&NpRbk[0]) , ncells*sizeof( double ) );
	
	stream.read( (char*)(&DNA_spread[0]) , ncells*sizeof( double ) );
	
	// le copie per i calcoli in Diff non vengono scritte su file
	
	stream.close();

}


//  ******************** TimersAdvance ********************

// funzione per l'avanzamento dei contatori del run 
// 
bool CellsSystem::TimersAdvance( )
{

    double timestep; 
    
    // selezione del timestep
    if( treal < tsm_start || treal >= tsm_stop)      // in questo intervallo ... slow motion
        {
        timestep = dt;
        slow_motion = false;
        }
    else
        {
        timestep = dt_sm;
        slow_motion = true;
        }  
    
	t += timestep;							// update del tempo di simulazione
	time_from_CGAL += timestep;				// update del tempo di simulazione dall'ultima chiamata a CGAL
	if(ready2start) treal += timestep;		// update del tempo che e' passato dalla partenza vera e propria della simulazione
	
	nstep++;
	
	if(!ready2start && t > t_ini)			// condizione di fine dell'inizializzazione dello stato cellulare
		{
		cout << "\nFine dell'inizializzazione\n" << endl;
		ready2start = true;
		nstep_start = nstep;				// qui si memorizza il numero del passo che corrisponde alla fine dell'inizializzazione
		}	

	if(treal < tmax && ( t_CPU < t_CPU_max || t_CPU_max <= 0 ) )
		return true;
	else 
		return false;
}

bool CellsSystem::TimersAdvanceUntil( double endtime)
{

    double timestep; 
    
    // selezione del timestep
    if( treal < tsm_start || treal >= tsm_stop)      // in questo intervallo ... slow motion
        {
        timestep = dt;
        slow_motion = false;
        }
    else
        {
        timestep = dt_sm;
        slow_motion = true;
        }  
    
	t += timestep;							// update del tempo di simulazione
	time_from_CGAL += timestep;				// update del tempo di simulazione dall'ultima chiamata a CGAL
	if(ready2start) treal += timestep;		// update del tempo che e' passato dalla partenza vera e propria della simulazione
	
	nstep++;
	
	if(!ready2start && t > t_ini)			// condizione di fine dell'inizializzazione dello stato cellulare
		{
		cout << "\nFine dell'inizializzazione\n" << endl;
		ready2start = true;
		nstep_start = nstep;				// qui si memorizza il numero del passo che corrisponde alla fine dell'inizializzazione
		}	

	if(treal < endtime && ( t_CPU < t_CPU_max || t_CPU_max <= 0 ) )
		return true;
	else 
		return false;
}

//  ******************** CPU_timer ********************
//
// funzione che restituisce il tempo di CPU in secondi
//
// 
void CellsSystem::CPU_timer( timer_button button )
{
#pragma omp critical
	{
	switch (button) 
		{
		
		case Start_timer:
#ifdef _OPENMP
			CPU_time_0 = CPU_time_1 = omp_get_wtime();
#else
			CPU_time_0 = CPU_time_1 = clock();
#endif
			t_CPU = 0.;
			t_CPU_int = 0.;
			t_CPU_0 = 0.;
			delta_t_CPU = 0.;
			break;
			
		case Start_intertime:
#ifdef _OPENMP
			CPU_time_1 = omp_get_wtime();
			t_CPU += CPU_time_1-CPU_time_0;
#else
			CPU_time_1 = clock();
			t_CPU += ( (double)(CPU_time_1-CPU_time_0) )/CLOCKS_PER_SEC;
#endif
			CPU_time_0 = CPU_time_1;
			break;
			
		case Stop_intertime:
#ifdef _OPENMP
			CPU_time_1 = omp_get_wtime();
			t_CPU_int += (CPU_time_1-CPU_time_0);
			t_CPU += (CPU_time_1-CPU_time_0);
#else
			CPU_time_1 = clock();
			t_CPU_int += ( (double)(CPU_time_1-CPU_time_0) )/CLOCKS_PER_SEC;
			t_CPU += ( (double)(CPU_time_1-CPU_time_0) )/CLOCKS_PER_SEC;
#endif
			CPU_time_0 = CPU_time_1;
			break;
			
		case Restart_timer:
#ifdef _OPENMP
			CPU_time_1 = omp_get_wtime();
			t_CPU += (CPU_time_1-CPU_time_0);
#else
			CPU_time_1 = clock();
			t_CPU += ( (double)(CPU_time_1-CPU_time_0) )/CLOCKS_PER_SEC;
#endif
			delta_t_CPU = t_CPU-t_CPU_0;
			t_CPU_0 = t_CPU;
			CPU_time_0 = CPU_time_1;
			break;
			
		case Clear_intertime:
#ifdef _OPENMP
			CPU_time_1 = omp_get_wtime();
			t_CPU += (CPU_time_1-CPU_time_0);
#else
			CPU_time_1 = clock();
			t_CPU += ( (double)(CPU_time_1-CPU_time_0) )/CLOCKS_PER_SEC;
#endif
			t_CPU_int = 0.;
			CPU_time_0 = CPU_time_1;
			break;
			
		default:
			break;
			
		}
	}
		
}



//  ******************** Timing ********************
//
// funzione che restituisce il tempo trascorso in secondi
//
// 
double CellsSystem::Timing( bool reset )
{

	if(reset) 
		{
		time(&time_old);
		timing = 0.;
		}
	else
		{
		time(&time_now);
		timing = difftime(time_now, time_old);
		time_old = time_now;
		}
		
	return timing;

}

//  ******************** Printout ********************
//
//		summary printout
//
// ***************************************************************
//
void CellsSystem::Printout()
{

	static bool first_print=true;

// stampa degli header nel caso che questa sia la prima volta che si stampa
	if(first_print)
		{
		first_print = false;
		
		cout << "\n   nstep     CPU time  (%diff)       time |         t       treal |lastname ncells alive|   %G0  %G1m  %G1p    %S   %G2   %M  %dead  |    volume" \
			<< "    sph.r       %V |    av_vol     av_r   av_M   av_ATP |   env_pH" << endl;
		
		screen_dump_file << "nstep\ttotal CPU time\tdelta CPU time\t(%diff)" \
			<< "\ttime\tt\ttreal\tnumtot\tncells\talive\t%G0\t%G1m\t%G1p\t%S\t%G2\t%M\t%dead\tvolume" \
			<< "\tspheroid radius\t%V\tav_vol\tav_r\tav_M\tav_ATP\tenv_pH" \
			<< "\tDC0\tDC1\tDC2\tDC3\tDC4\tDC5\tDC6" \
			<< "\tidum\tnrepeats_average\tnrepeats_max\tnrepeats_min" \
			<< "\tloop_count_average\tloop_count_max\tloop_count_min" \
			<< "\tn_mitoses_average\tn_mitoses_max\tn_mitoses_min" \
			<< "\tmin_Gin\tmax_Gin\tmin_Gext\tmax_Gext\tmin_O2\tmax_O2\tmin_Ain\tmax_Ain\tmin_Aext\tmax_Aext" \
			<< "\tmin_AcLin\tmax_AcLin\tmin_AcLext\tmax_AcLext\tmin_extvolume\tmax_extvolume" << endl;
		
		convlog_file << "nstep\tG_env\tO2env\tA_env\tAcL_env\tGin\tGext\tO2\tAin\tAext\tAcLin\tAcLext\tATPp" << endl;
		}

// variabili locali per le statistiche	
	CellPhase fase;
	int phase_counter[Nphase];
	for(int k=0; k<Nphase; k++) phase_counter[k] = 0;
	
	int death_cond[7];
	for(int k=0; k<7; k++) death_cond[k] = 0;
	
	double av_mitocondri = 0;
	double volume_totale = 0;
	double volume_totale_vive = 0;
	double ATP_totale = 0;
	double ATP_prod_tot = 0;
	double ATP_cons_tot = 0;
	
// loop sulle cellule per le statistiche	
	for(unsigned long n=0; n<ncells; n++)
		{
		
		volume_totale += volume[n];
		
		fase = phase[n];
		phase_counter[fase]++;
		
		if(fase != dead)
			{
			av_mitocondri += M[n];	
			volume_totale_vive += volume[n];	
			ATP_totale += ATPp[n];
			if(ready2start)
				{
				ATP_prod_tot += ATPprod[n];
				ATP_cons_tot += ATPcons[n];
				}
			
			}
		else 
			{
			int mask = 1;
			for(int k=0; k<7; k++)
				{
				int index = death_condition[n] & (mask << k);
				if(index) death_cond[k]++;
				}
			}

		}

	av_mitocondri /= alive;			// numero medio di mitocondri per cellula
	ATP_totale /= alive;			// ATPp medio per cellula
	
// output su schermo

	cout << setw(8) << right << nstep;																					// numero di passi
	cout << " " << setw(12) << setprecision(3) << fixed << right << showpoint << delta_t_CPU;							// tempo di CPU
	cout << " " << setw(8) << setprecision(1) << fixed << right << showpoint << (delta_t_CPU>0 ? 100.*t_CPU_int/delta_t_CPU : 0);	// frazione di tempo di CPU passata in diff
	cout << " " << setw(10) << setprecision(0) << fixed << right << showpoint << timing;								// tempo reale
	cout << " | " << setw(10) << setprecision(3) << fixed << right << showpoint << t/3600.;								// tempo di simulazione
    if(!slow_motion)
        {
        cout << " " << setw(10) << setprecision(3) << fixed << right << showpoint << treal/3600.;							// tempo dalla partenza
        cout << " ";
        }
    else
        {
        cout << " " << setw(10) << setprecision(6) << fixed << right << showpoint << treal-tsm_start;							// tempo dalla partenza
        cout << "s";
        }
	cout << "| " << setw(6) << right << lastname;																		// numero totale di cellule 
	cout << setw(6) << right << ncells;																		// numero totale di cellule
	cout << " " << setw(6) << right << alive << " |";																	// numero di cellule vive
	for(int k=0; k<Nphase-1; k++) cout << " " << setw(5) << setprecision(1) << fixed << right << 100.*((double)phase_counter[k])/alive; // % di cellule vive in ciascuna fase cellulare
	cout << " " << setw(5) << setprecision(1) << fixed << right << 100*((double)phase_counter[Nphase-1])/ncells;	// % di cellule morte rispetto al totale
	cout << "  | " << setw(9) << setprecision(3) << scientific << right << volume_totale;								// volume totale occupato dalle cellule
	
	cout << " " << setw(8) << setprecision(2) << fixed << right << pow( 3.*volume_totale/(4.*PI), (double)0.333333333);		// raggio dello sferoide in micron
	cout << " " << setw(8) << setprecision(5) << fixed << right << 100.*volume_totale/Env.GetEnvironmentvolume0();		// % volume occupato dalle cellule
	cout << " | " << setw(9) << setprecision(3) << scientific << right << volume_totale_vive/alive;						// volume medio di una cellula viva
	cout << " " << setw(8) << setprecision(2) << fixed << right << pow(3*volume_totale_vive/(alive*4*PI),(double)1./3.);		// raggio medio di una cellula viva
	cout << " " << setw(6) << setprecision(1) << fixed << right << av_mitocondri;										// numero medio di mitocondri
	cout << " " << setw(8) << setprecision(2) << scientific << right << ATP_totale;										// ATPp medio per cellula

	cout << " | " << setw(8) << setprecision(6) << fixed << right << Env.GetEnvironmentpH();							// pH ambientale
	cout << " | ( ";
	for(int k=0; k<7; k++) cout << death_cond[k] << " ";																// death condition
	cout << ")";
	cout << endl;
	
	cout.flush();
	
	// screen dump su file
	
	screen_dump_file << nstep;													// numero di passi
	screen_dump_file << "\t" << t_CPU;											// tempo di CPU integrato
	screen_dump_file << "\t" << delta_t_CPU;									// tempo di CPU
	screen_dump_file << "\t" << (delta_t_CPU>0 ? 100.*t_CPU_int/delta_t_CPU : 0);	// frazione di tempo di CPU passata in diff
	screen_dump_file << "\t" << timing;											// tempo reale
	screen_dump_file << "\t" << fixed << t/3600.;								// tempo di simulazione
	screen_dump_file << "\t" << fixed << treal/3600.;							// tempo dalla partenza
	screen_dump_file << "\t" << lastname;										// numero totale di cellule
	screen_dump_file << "\t" << ncells;											// numero totale di cellule
	screen_dump_file << "\t" << alive;											// numero di cellule vive
	for(int k=0; k<Nphase-1; k++) screen_dump_file << "\t" << 100*((double)phase_counter[k])/alive;	// conteggio in ciascuna fase cellulare
	screen_dump_file << "\t" << 100*((double)phase_counter[Nphase-1])/ncells;	// % di cellule morte rispetto al totale
	screen_dump_file << "\t" << scientific << volume_totale;					// volume totale occupato dalle cellule
	
	screen_dump_file << "\t" << pow( 3.*volume_totale/(4.*PI), (double)0.333333333);// raggio dello sferoide in micron
	screen_dump_file << "\t" << 100.*volume_totale/Env.GetEnvironmentvolume0();	// % volume occupato dalle cellule
	screen_dump_file << "\t" << scientific << volume_totale_vive/alive;			// volume medio di una cellula viva
	screen_dump_file << "\t" << pow(3*volume_totale_vive/(alive*4*PI), (double)1./3.);// raggio medio di una cellula viva
	screen_dump_file << "\t" << av_mitocondri;									// numero medio di mitocondri
	screen_dump_file << "\t" << scientific << ATP_totale;						// ATPp medio per cellula
	screen_dump_file << "\t" << Env.GetEnvironmentpH();							// pH ambientale
	for(int k=0; k<7; k++) screen_dump_file << "\t" << death_cond[k];			// death condition
	screen_dump_file << "\t" << Get_idum();										// valore attuale del seme del generatore di numeri casuali
	screen_dump_file << "\t" << nrepeats_average;								// numero medio di iterazioni del metodo di diffusione
	screen_dump_file << "\t" << nrepeats_max;									// numero max di iterazioni del metodo di diffusione
	screen_dump_file << "\t" << nrepeats_min;									// numero min di iterazioni del metodo di diffusione
	screen_dump_file << "\t" << loop_count_average;								// numero medio di iterazioni del metodo di calcolo della posizione
	screen_dump_file << "\t" << loop_count_max;									// numero max di iterazioni del metodo di calcolo della posizione
	screen_dump_file << "\t" << loop_count_min;									// numero min di iterazioni del metodo di calcolo della posizione

	screen_dump_file << "\t" << n_mitoses_average;								// numero medio di mitosi passo-passo
	screen_dump_file << "\t" << n_mitoses_max;									// numero max di mitosi passo-passo
	screen_dump_file << "\t" << n_mitoses_min;									// numero min di mitosi passo-passo

	screen_dump_file << "\t" << scientific << min_Gin;							// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << max_Gin;							// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << min_Gext;							// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << max_Gext;							// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << min_O2;							// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << max_O2;							// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << min_Ain;							// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << max_Ain;							// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << min_Aext;							// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << max_Aext;							// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << min_AcLin;						// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << max_AcLin;						// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << min_AcLext;						// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << max_AcLext;						// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << min_extvolume;					// altre statistiche estreme passo-passo (da StepStat)
	screen_dump_file << "\t" << scientific << max_extvolume;					// altre statistiche estreme passo-passo (da StepStat)

	screen_dump_file << endl;
	
	screen_dump_file.flush();
	
	
	// printout dei dati relativi alla convergenza
	convlog_file << nstep << "\t";
	for(int k=0; k<NCONV_TEST-1; k++)
		convlog_file << convergence_fail[k] << "\t";
	convlog_file << convergence_fail[NCONV_TEST-1] << endl;
	
	convlog_file.flush();


}

//  ******************** Print2file ********************
//
//		summary printout to file
//
// ***************************************************************
//
void CellsSystem::Print2file()
{

	unsigned long n;
	int k;
	
	static bool first_print2file=true;

// stampa dell'header nel caso che questa sia la prima volta che si stampa
	if(first_print2file)
		{
		first_print2file = false;
		output_file << "nstep\t t\t treal\t flow\t dose_rate\t ncells\t alive\t G0\t G1p\t G1m\t S\t G2\t M\t dead\t %S\t" \
			<< "%vol\t %vol_alive\t av_vol\t av_r\t av_M\t dev.st._M\t min_M\t max_M\t av_DNA\t dev.st._DNA\t min_DNA\t max_DNA" \
			<< "\t av_ATPp\t ATP_prod\t ATP_cons\t env_vol" \
			<< "\t env_Gconc\t env_Aconc" \
			<< "\t env_AcLconc\t env_pH\t av_mother_age\t av_GAbsRate\t av_GConsRate\t av_AAbsRate\t av_AConsRate\t av_StoreFillRate "\
			<< "\t av_StoreConsRate\t av_AcLRate\t av_AcLOutRate\t av_ATP_Ox\t av_ATP_NOx\t av_ATP2\t av_ATP3"\
			<< "\t av_ConsATP\t av_ConsATP_1\t av_ConsATP_2\t av_ConsATP_3\t av_ConsATP_5\t av_ATPtot\t av_G\t av_A\t av_O2\t av_AcL"\
			<< "\t AcLFlow\t O2flow\t av_neigh\t dev.st._neigh\t min_neigh\t max_neigh" << endl;
		}

// variabili locali per le statistiche	
	int phase_counter[Nphase];
	for(k=0; k<Nphase; k++) phase_counter[k] = 0;
	
	double av_mitocondri = 0;
	double av2_mitocondri = 0;
	double max_mitocondri = 0;
	double min_mitocondri = 10000;
	double volume_totale = 0;
	double volume_totale_vive = 0;
	double ATP_totale = 0;
	// double pRb_totale = 0;
	double ATP_prod_tot = 0;
	double ATP_cons_tot = 0;
	
	double av_age_mother = 0.;
	
	// rates medi
	double av_GAbsRate = 0.;		// rate di assorbimento del glucosio
	double av_GConsRate = 0.;		// rate di consumo del glucosio
	double av_AAbsRate = 0.;		// rate di assorbimento della glutammina
	double av_AConsRate = 0.;		// rate di consumo del glutammmina
	double av_StoreFillRate = 0.;	// rate di riempiemento dello store
	double av_StoreConsRate = 0.;	// rate di consumo dello store
	double av_AcLRate = 0.;		// rate di produzione dell'acido lattico
	double av_AcLOutRate = 0.;		// rate di espulsione dell'acido lattico
	
	double av_ATP_Ox = 0.;			// rate di produzione di ATP tramite fosforilazione ossidativa
	double av_ATP_NOx = 0.;		// rate di produzione di ATP tramite glicolisi anaerobica
	double av_ATP2 = 0.;			// rate di produzione di ATP dallo store
	double av_ATP3 = 0.;			// rate di produzione di ATP dagli altri nutrienti
	double av_ConsATP = 0.;		// rate di consumo di ATP associato al metabolismo
	double av_ConsATP_1 = 0.;		// rate di consumo di ATP associato ai mitocondri
	double av_ConsATP_2 = 0.;		// rate di consumo di ATP associato alla produzione di proteine
	double av_ConsATP_3 = 0.;		// rate di consumo di ATP associato alla produzione di DNA
	//double av_ConsATP_4 = 0.;	// rate di consumo di ATP associato alla pompa protonica
	double av_ConsATP_5 = 0.;		// rate di consumo di ATP associato alla produzione di mtDNA
	double av_ATPtot = 0.;			// rate totale di variazione dell'ATP = somma dei rates
	
	double av_Gin = 0.;			// valore medio del gluocosio interno
	double av_Ain = 0.;			// valore medio della glutammina interna
	double av_O2 = 0.;				// valore medio dell'ossigeno
	double av_AcLin = 0.;			// valore medio dell'AcL
	
	double av_DNA = 0;				// valore medio del DNA
	double av2_DNA = 0;
	double min_DNA = 10;
	double max_DNA = -10;
	
	// geometria
	double av_neigh = 0;			// numero medio di vicini
	double av2_neigh = 0;
	int min_neigh = 10000;
	int max_neigh = 0;
	
// loop sulle cellule per le statistiche	
	for(n=0; n<ncells; n++)
		{
		
		volume_totale += volume[n];
		
		CellPhase fase = phase[n];
		phase_counter[fase]++;
		
		if(fase != dead)
			{
			av_mitocondri += M[n];	
			av2_mitocondri += SQR( M[n] );	
			if( M[n] > max_mitocondri ) max_mitocondri = M[n];
			if( M[n] < min_mitocondri ) min_mitocondri = M[n];
			volume_totale_vive += volume[n];	
			ATP_totale += ATPp[n];
			// pRb_totale += cellula.Get_pRb();
			if(ready2start)	// quando la simulazione e' pronta allora si somma produzione e consumo di ATP di tutte le cellule
				{
				ATP_prod_tot += ATPprod[n];
				ATP_cons_tot += ATPcons[n];
				}
				
			av_age_mother += age_mother[n];
			
			av_GAbsRate += GAbsRate[n];
			av_GConsRate += GConsRate[n];
			av_AAbsRate += AAbsRate[n];
			av_AConsRate += AConsRate[n];
			av_StoreFillRate += StoreFillRate[n];
			av_StoreConsRate += StoreConsRate[n];
			av_AcLRate += AcLRate[n];
			av_AcLOutRate += AcLOutRate[n];
			
			av_ATP_Ox += ATP_Ox[n];
			av_ATP_NOx += ATP_NOx[n];
			av_ATP2 += ATP2[n];
			av_ATP3 += ATP3[n];
			av_ConsATP += ConsATP[n];
			av_ConsATP_1 += ConsATP_1[n];
			av_ConsATP_2 += ConsATP_2[n];
			av_ConsATP_3 += ConsATP_3[n];
			av_ConsATP_5 += ConsATP_5[n];
			av_ATPtot += ATPtot[n];
			
			av_Gin += G[n];
			av_Ain += A[n];
			av_O2 += O2[n];
			av_AcLin += AcL[n];
						
			av_DNA += DNA[n];	
			av2_DNA += SQR( DNA[n] );	
			if( DNA[n] > max_DNA ) max_DNA = DNA[n];
			if( DNA[n] < min_DNA ) min_DNA = DNA[n];

			}
		
		// l'informazione geometrica viene calcolata per tutte le cellule, anche quelle morte
		
		av_neigh += neigh[n];	
		av2_neigh += SQR( neigh[n] );	
		if( neigh[n] > max_neigh ) max_neigh = neigh[n];
		if( neigh[n] < min_neigh ) min_neigh = neigh[n];

		}

	av_mitocondri /= alive;			// numero medio di mitocondri per cellula
	av2_mitocondri /= alive;		// secondo momento del numero di mitocondri per cellula
	av_DNA /= alive;				// quantita' media di DNA per cellula
	av2_DNA /= alive;				// secondo momento del DNA per cellula
	ATP_totale /= alive;			// ATPp medio per cellula
	// pRb_totale /= alive;			// pRb media per cellula
	av_age_mother /= alive;			// eta' media della madre
	av_age_mother /= 86400.;		// l'eta' media della madre viene convertita in giorni
	
	av_GAbsRate /= alive;			// valore medio del rate di assorbimento del glucosio
	av_GConsRate /= alive;			// valore medio del rate di consumo del glucosio
	av_AAbsRate /= alive;			// valore medio del rate di assorbimento della glutammina
	av_AConsRate /= alive;			// valore medio del rate di consumo del glutammmina
	av_StoreFillRate /= alive;		// valore medio del rate di riempiemento dello store
	av_StoreConsRate /= alive;		// valore medio del rate di consumo dello store
	av_AcLRate /= alive;			// valore medio del rate di produzione dell'acido lattico
	av_AcLOutRate /= alive;			// valore medio del rate di espulsione dell'acido lattico

	av_ATP_Ox /= alive;				// valore medio del rate di produzione di ATP tramite fosforilazione ossidativa
	av_ATP_NOx /= alive;			// valore medio del rate di produzione di ATP tramite glicolisi anaerobica
	av_ATP2 /= alive;				// valore medio del rate di produzione di ATP dallo store
	av_ATP3 /= alive;				// valore medio del rate di produzione di ATP dagli altri nutrienti
	av_ConsATP /= alive;			// valore medio del rate di consumo di ATP associato al metabolismo
	av_ConsATP_1 /= alive;			// valore medio del rate di consumo di ATP associato ai mitocondri
	av_ConsATP_2 /= alive;			// valore medio del rate di consumo di ATP associato alla produzione di proteine
	av_ConsATP_3 /= alive;			// valore medio del rate di consumo di ATP associato alla produzione di DNA
	av_ConsATP_5 /= alive;			// valore medio del rate di consumo di ATP associato alla produzione di mtDNA
	av_ATPtot /= alive;				// valore medio del rate totale di variazione dell'ATP = somma dei rates
	
	av_Gin /= alive;				// valore medio del glucosio interno
	av_Ain /= alive;				// valore medio della glutammina interna
	av_O2 /= alive;					// valore medio dell'ossigeno
	av_AcLin /= alive;				// valore medio dell'acido lattico interno

	av_neigh /= alive;				// numero medio di vicini per cellula
	av2_neigh /= alive;				// secondo momento del numero di vicini per cellula

// output	

	output_file << nstep;												// numero di passi
	output_file << "\t" << t/3600.;										// tempo di simulazione
	output_file << "\t" << treal/3600.;									// tempo dalla partenza
	if(ready2start)
		{
		output_file << "\t" << flowSignal.SignalValue(treal)/((1.e-9)/60.);	// flusso attuale (in microlitri/min)
		output_file << "\t" << Env.GetEnvironmentDoseRate();			// dose (Gy/s)
		// output_file << "\t" << setprecision(6) << right << showpoint << dose_rateSignal.SignalValue(t);
		}
	else
		{
		output_file << "\t" << 0.;										// flusso attuale (in microlitri/min)
		output_file << "\t" << 0.;										// flusso attuale (in microlitri/min)
		}
		
	output_file << "\t" << ncells;													// numero totale di cellule
	output_file << "\t" << alive;													// numero di cellule vive
	for(k=0; k<Nphase; k++) output_file << "\t" << phase_counter[k];				// conteggio in ciascuna fase cellulare
	output_file << "\t" << ((double)phase_counter[S_phase])/alive;				// % di cellule in fase S
	
	output_file << "\t" << 100.*volume_totale/Env.GetEnvironmentvolume0();			// % volume occupato dalle cellule
	output_file << "\t" << 100.*volume_totale_vive/Env.GetEnvironmentvolume0();		// % volume occupato dalle cellule vive
	output_file << "\t" << scientific << volume_totale_vive/alive;					// volume medio di una cellula viva
	output_file << "\t" << scientific << pow(3*volume_totale_vive/(alive*4*PI),(double)1./3.);	// raggio medio di una cellula viva
	output_file << "\t" << av_mitocondri;											// numero medio di mitocondri
	output_file << "\t" << sqrt(fabs(av2_mitocondri-SQR(av_mitocondri)));			// dev. st. del numero di mitocondri
	output_file << "\t" << min_mitocondri;											// numero min di mitocondri
	output_file << "\t" << max_mitocondri;											// numero max di mitocondri
	output_file << "\t" << scientific << av_DNA;									// DNA medio
	output_file << "\t" << scientific << sqrt(fabs(av2_DNA-SQR(av_DNA)));			// dev. st. del DNA
	output_file << "\t" << scientific << min_DNA;									// DNA minimo
	output_file << "\t" << scientific << max_DNA;									// DNA massimo
	output_file << "\t" << scientific << ATP_totale;								// ATPp medio per cellula
	output_file << "\t" << scientific << ATP_prod_tot;								// ATP consumato in totale
	output_file << "\t" << scientific << ATP_cons_tot;								// ATP prodotto in totale
	// output_file << "\t" << setw(6) << right << pRb_totale;						// pRb media per cellula
	
	output_file << "\t" << scientific << Env.GetEnvironmentvolume();				// volume libero dell'ambiente
	output_file << "\t" << Env.GetEnvironmentG()/Env.GetEnvironmentvolume();		// concentrazione di glucosio nell'ambiente
	output_file << "\t" << Env.GetEnvironmentA()/Env.GetEnvironmentvolume();		// concentrazione di altri nutrienti nell'ambiente
	output_file << "\t" << Env.GetEnvironmentAcL()/Env.GetEnvironmentvolume();		// concentrazione di AcL nell'ambiente
	output_file << "\t" << Env.GetEnvironmentpH();									// pH ambientale
	output_file << "\t" << av_age_mother;											// eta' media della madre in giorni
	
	// valori medi di alcuni rates
	output_file << "\t" << scientific << av_GAbsRate;								// rate di assorbimento del glucosio
	output_file << "\t" << scientific << av_GConsRate;								// rate di consumo del glucosio
	output_file << "\t" << scientific << av_AAbsRate;								// rate di assorbimento della glutammina
	output_file << "\t" << scientific << av_AConsRate;								// rate di consumo del glutammmina
	output_file << "\t" << scientific << av_StoreFillRate;							// rate di riempiemento dello store
	output_file << "\t" << scientific << av_StoreConsRate;							// rate di consumo dello store
	output_file << "\t" << scientific << av_AcLRate;								// rate di produzione dell'acido lattico
	output_file << "\t" << scientific << av_AcLOutRate;								// rate di espulsione dell'acido lattico

	output_file << "\t" << scientific << av_ATP_Ox;									// rate di produzione di ATP tramite fosforilazione ossidativa
	output_file << "\t" << scientific << av_ATP_NOx;								// rate di produzione di ATP tramite glicolisi anaerobica
	output_file << "\t" << scientific << av_ATP2;									// rate di produzione di ATP dallo store
	output_file << "\t" << scientific << av_ATP3;									// rate di produzione di ATP dagli altri nutrienti
	output_file << "\t" << scientific << av_ConsATP;								// rate di consumo di ATP associato al metabolismo
	output_file << "\t" << scientific << av_ConsATP_1;								// rate di consumo di ATP associato ai mitocondri
	output_file << "\t" << scientific << av_ConsATP_2;								// rate di consumo di ATP associato alla produzione di proteine
	output_file << "\t" << scientific << av_ConsATP_3;								// rate di consumo di ATP associato alla produzione di DNA
	output_file << "\t" << scientific << av_ConsATP_5;								// rate di consumo di ATP associato alla produzione di mtDNA
	output_file << "\t" << scientific << av_ATPtot;									// rate totale di variazione dell'ATP = somma dei rates

	output_file << "\t" << scientific << av_Gin;									// valore medio del glucosio interno
	output_file << "\t" << scientific << av_Ain;									// valore medio della glutammina interna
	output_file << "\t" << scientific << av_O2;										// valore medio dell'ossigeno 
	output_file << "\t" << scientific << av_AcLin;									// valore medio dell'AcL interno

	output_file << "\t" << scientific << AcLFlow;									// flusso di acido lattico verso l'ambiente
	output_file << "\t" << scientific << O2Flow;									// assorbimento totale di ossigeno
	
	output_file << "\t" << av_neigh;												// numero medio di vicini
	output_file << "\t" << sqrt(fabs(av2_neigh-SQR(av_neigh)));						// dev. st. del numero di vicini
	output_file << "\t" << min_neigh;												// numero min di vicini
	output_file << "\t" << max_neigh;												// numero max di vicini

	output_file << endl;
	
	output_file.flush();
	
}

//  ******************** Print2logfile ********************
//
// questo metodo stampa un record su logfile per scopi di debugging
//
// ***************************************************************
//
void CellsSystem::Print2logfile(string str)
{
	log_file << "\n" << str << "\n" << endl;
	log_file << "\n\n***** Step " << nstep << " *****\n\n" << endl;
	log_file << "Environment \n" << endl;
	log_file << Env << "\n\n" << endl;
	log_file << "*** ncells = " << ncells << " ***\n" << endl;
	log_file << "First cell (0)\n" << endl;
	PrintCell( log_file, 0 );
	log_file << "\n\n" << endl;
	log_file << "Cell type (0)\n" << endl;
	log_file << *(type[0]) << endl;

	if(ncells > 1) 
		{
		log_file << "Second cell" << endl;
		PrintCell( log_file, 1 );
		log_file << "\n\n" << endl;
		}
	if(ncells > 2) 
		{
		log_file << "Last cell (" << ncells-1 << ")\n" << endl;
		PrintCell( log_file, ncells-1 );
		log_file << "\n\n" << endl;
		}

}

//  ******************** PrintAll2logfile ********************
//
// questo metodo stampa TUTTE LE CELLULE su logfile per scopi di debugging (ATTENZIONE!!!)
//
// ***************************************************************
//
void CellsSystem::PrintAll2logfile(string str)
{
	log_file << "\n" << str << "\n" << endl;
	log_file << "\n\n***** Step " << nstep << " *****\n\n" << endl;
	log_file << "Environment \n" << endl;
	log_file << Env << "\n\n" << endl;
	log_file << "*** ncells = " << ncells << " ***\n" << endl;
	
	for(unsigned int n=0; n<ncells; n++)
		{
		log_file << "Cell " << n << "\n" << endl;
		PrintCell( log_file, 0 );
		log_file << "\n\n" << endl;
		}

	log_file << "\n\n" << endl;

}


//  ******************** PrintConfiguration ********************
//
// questo metodo scrive una configurazione su file
//
// ***************************************************************
//
void CellsSystem::PrintConfiguration(bool isBinary)
{

	char str[100];
	if(isBinary)
		{
		if(!slow_motion) 
            sprintf(str,"%lu.bin",nconfiguration);
        else
            sprintf(str,"%lu_sm.bin",nconfiguration);

		string configuration_record = configuration_b_filename+str;
		configuration_b_file.open( configuration_record.c_str(), ios::binary );
		
		configuration_b_file.write( (char*)(&sim_type), sizeof(int) );
		}
	else
		{
		if(!slow_motion) 
            sprintf(str,"%lu.txt",nconfiguration);
        else
            sprintf(str,"%lu_sm.txt",nconfiguration);

		string configuration_record = configuration_filename+str;
		configuration_file.open( configuration_record.c_str() );
		
		configuration_file << sim_type << endl;
		}
	
	PrintHeader(isBinary);		// header
	PrintEnv(isBinary);			// concentrazioni ambientali di G, O2, A, AcL
	PrintPoints(isBinary);		// posizioni dei centri
	PrintLinks(isBinary);		// links (triangolazione di Delaunay)
	PrintCHFlag(isBinary);		// flag che indica se la cellula e' sul convex hull
	PrintASFlag(isBinary);		// flag che indica se la cellula e' sull'alpha shape
	PrintR(isBinary);			// raggi cellulari
	PrintPhase(isBinary);		// fase
    PrintType(isBinary);    // cell type
	PrintVar(isBinary);			// stampa di un sottoinsieme di variabili cellulari
								// (concentrazione di O2, G, Gextra, A, Aextra, AcL, AcLextra, ATPp, AcLrate, pHextra)
	
	if(isBinary)
		configuration_b_file.close();
	else
		configuration_file.close();
	
}

//  ******************** PrintHeader ********************

// printout dell'header di una singola configurazione
//
void CellsSystem::PrintHeader(bool isBinary)
{
	if( isBinary )
		{
		configuration_b_file.write( (char*)(&nstep), sizeof(nstep) );
		double treal_d = (double )treal;
		configuration_b_file.write( (char*)(&treal_d), sizeof(treal_d) );
		configuration_b_file.write( (char*)(&ncells), sizeof(ncells) );
		}
	else
		configuration_file << nstep << " " << treal << " " << ncells << endl;
}


//  ******************** PrintPoints ********************

// printout dei punti su file 
// i cast espliciti servono ad evitare warnings del compilatore (e comunque rendono chiara la scelta fatta)
//
void CellsSystem::PrintPoints(bool isBinary)
{
	for(unsigned long k=0; k<ncells; k++)
		{
		double xpr = (double)x[k];
		double ypr = (double)y[k];
		double zpr = (double)z[k];
		double vxpr = (double)vx[k];
		double vypr = (double)vy[k];
		double vzpr = (double)vz[k];
	
		if( isBinary )
			{
			configuration_b_file.write( (char*)(&xpr), sizeof( double) );
			configuration_b_file.write( (char*)(&ypr), sizeof( double) );
			configuration_b_file.write( (char*)(&zpr), sizeof( double) );
			configuration_b_file.write( (char*)(&vxpr), sizeof( double) );
			configuration_b_file.write( (char*)(&vypr), sizeof( double) );
			configuration_b_file.write( (char*)(&vzpr), sizeof( double) );
			}
		else
			{
			configuration_file << xpr << " " << ypr << " " << zpr << endl;
			configuration_file << vxpr << " " << vypr << " " << vzpr << endl;
			}
		}
}


//  ******************** PrintLinks ********************

// questo metodo stampa solo i links su file
void CellsSystem::PrintLinks(bool isBinary)
{
		
	for(unsigned long k=0; k<ncells; k++)
		{
		
		if( isBinary )
			{
			configuration_b_file.write( (char*)(&neigh[k]), sizeof(int) );
			for( int kk=0; kk< neigh[k]; kk++)
				configuration_b_file.write( (char*)(&vneigh[k][kk]), sizeof(int) );
			}
		else
			{
			configuration_file << neigh[k];
			for( int kk=0; kk< neigh[k]; kk++)
				configuration_file << " " << vneigh[k][kk];
			configuration_file << endl;
			}
			
		}

}


//  ******************** PrintCHFlag ********************

// questo metodo stampa solo la flag di appartenenza al CH su file
void CellsSystem::PrintCHFlag(bool isBinary)
{
	for(unsigned long k=0; k<ncells; k++)
		{
		
		bool IsOnCH = isonCH[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&IsOnCH), sizeof(bool) );
		else
			configuration_file << " " << IsOnCH;
			
		}
		
	if( !isBinary ) configuration_file << endl;

}

//  ******************** PrintASFlag ********************

// questo metodo stampa solo la flag di appartenenza all'AS su file
void CellsSystem::PrintASFlag(bool isBinary)
{
	for(unsigned long k=0; k<ncells; k++)
		{
		
		bool IsOnAS = isonAS[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&IsOnAS), sizeof(bool) );
		else
			configuration_file << " " << IsOnAS;
			
		}
		
	if( !isBinary ) configuration_file << endl;

}

//  ******************** PrintBVFlag ********************

// this method prints on file the isonBV flag
void CellsSystem::PrintBVFlag(bool isBinary)
{
	for(unsigned long k=0; k<ncells; k++)
		{
		
		int IsOnBV = isonBV[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&IsOnBV), sizeof(int) );
		else
			configuration_file << " " << IsOnBV;
			
		}
		
	if( !isBinary ) configuration_file << endl;

}

//  ******************** PrintR ********************

// questo metodo stampa solo i raggi cellulari su file
void CellsSystem::PrintR(bool isBinary)
{
	for(unsigned long k=0; k<ncells; k++)
		{
		double rpr = (double)r[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&rpr), sizeof(double) );
		else
			configuration_file << " " << rpr;
		}
	if( !isBinary ) configuration_file << endl;

}


//  ******************** PrintPhase ********************

// questo metodo stampa il codice della fase cellulare
void CellsSystem::PrintPhase(bool isBinary)
{
	for(unsigned long k=0; k<ncells; k++)
		{
		if( isBinary )
			configuration_b_file.write( (char*)(&phase[k]), sizeof(int) );
		else		
			configuration_file << " " << phase[k];
		}
	if( !isBinary ) configuration_file << endl;

}

//  ******************** PrintType ********************

// questo metodo stampa il nome del CellType
void CellsSystem::PrintType(bool isBinary)
{
    for(unsigned long k=0; k<ncells; k++)
    {
        int name = type[k]->Get_name();
        
        if( isBinary )
            configuration_b_file.write( (char*)(&name), sizeof(int) );
        else
            configuration_file << " " << name;
    }
    if( !isBinary ) configuration_file << endl;
    
}


//  ******************** PrintVar ********************

// questo metodo stampa variabili importanti dentro le cellule
void CellsSystem::PrintVar(bool isBinary)
{
	// G
	for(unsigned long k=0; k<ncells; k++)
		{
		double Gpr = (double)G[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&Gpr), sizeof( double) );
		else
			configuration_file << " " << Gpr;
		}
	if( !isBinary ) configuration_file << endl;
	
	// G extracellulare
	for(unsigned long k=0; k<ncells; k++)
		{
		double concGextra = (double)(G_extra[k]/volume_extra[k]);
		if( isBinary )
			configuration_b_file.write( (char*)(&concGextra), sizeof( double) );
		else
			configuration_file << " " << concGextra;
		}
	if( !isBinary ) configuration_file << endl;
	
	// rate di assorbimento del glucosio
	for(unsigned long k=0; k<ncells; k++)
		{
		double GAbsRatepr = (double)GAbsRate[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&GAbsRatepr), sizeof( double) );
		else
			configuration_file << " " << GAbsRatepr;
		}
	if( !isBinary ) configuration_file << endl;
	
	// rate di consumo del glucosio
	for(unsigned long k=0; k<ncells; k++)
		{
		double GConsRatepr = (double)GConsRate[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&GConsRatepr), sizeof( double) );
		else
			configuration_file << " " << GConsRatepr;
		}
	if( !isBinary ) configuration_file << endl;
	
	// G6P
	for(unsigned long k=0; k<ncells; k++)
		{
		double G6Ppr = (double)G6P[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&G6Ppr), sizeof( double) );
		else
			configuration_file << " " << G6Ppr;
		}
	if( !isBinary ) configuration_file << endl;
	
	// O2
	for(unsigned long k=0; k<ncells; k++)
		{
		double concO2 = (double)(O2[k]/volume[k]);
		if( isBinary )
			configuration_b_file.write( (char*)(&concO2), sizeof( double) );
		else
			configuration_file << " " << concO2;
		}
	if( !isBinary ) configuration_file << endl;

	// Store
	for(unsigned long k=0; k<ncells; k++)
		{
		double storepr = (double)store[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&storepr), sizeof( double) );
		else
			configuration_file << " " << storepr;
		}
	if( !isBinary ) configuration_file << endl;

	// A
	for(unsigned long k=0; k<ncells; k++)
		{
		double Apr = (double)A[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&Apr), sizeof( double) );
		else
			configuration_file << " " << Apr;
		}
	if( !isBinary ) configuration_file << endl;
	
	// A extracellulare
	for(unsigned long k=0; k<ncells; k++)
		{
		double concAextra = (double)(A_extra[k]/volume_extra[k]);
		if( isBinary )
			configuration_b_file.write( (char*)(&concAextra), sizeof( double) );
		else
			configuration_file << " " << concAextra;
		}
	if( !isBinary ) configuration_file << endl;
	
	// rate di assorbimento della glutammina
	for(unsigned long k=0; k<ncells; k++)
		{
		double AAbsRatepr = (double)AAbsRate[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&AAbsRatepr), sizeof( double) );
		else
			configuration_file << " " << AAbsRatepr;
		}
	if( !isBinary ) configuration_file << endl;
	
	// rate di consumo della glutammina
	for(unsigned long k=0; k<ncells; k++)
		{
		double AConsRatepr = (double)AConsRate[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&AConsRatepr), sizeof( double) );
		else
			configuration_file << " " << AConsRatepr;
		}
	if( !isBinary ) configuration_file << endl;
	
	// AcL
	for(unsigned long k=0; k<ncells; k++)
		{
		double AcLpr = (double)AcL[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&AcLpr), sizeof( double) );
		else
			configuration_file << " " << AcLpr;
		}
	if( !isBinary ) configuration_file << endl;

	// AcL extracellulare
	for(unsigned long k=0; k<ncells; k++)
		{
		double concAcL = (double)(AcL_extra[k]/volume_extra[k]);
		if( isBinary )
			configuration_b_file.write( (char*)(&concAcL), sizeof( double) );
		else
			configuration_file << " " << concAcL;
		}
	if( !isBinary ) configuration_file << endl;

	// ATPp
	for(unsigned long k=0; k<ncells; k++)
		{
		double ATPppr = (double)ATPp[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&ATPppr), sizeof(double) );
		else
			configuration_file << " " << ATPppr;
		}
	if( !isBinary ) configuration_file << endl;

	// Mitocondri
	for(unsigned long k=0; k<ncells; k++)
		{
		double Mpr = (double)M[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&Mpr), sizeof( double) );
		else
			configuration_file << " " << Mpr;
		}
	if( !isBinary ) configuration_file << endl;

	// pH extracellulare
	for(unsigned long k=0; k<ncells; k++)
		{
		double pHpr = (double)pH[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&pHpr), sizeof( double) );
		else
			configuration_file << " " << pHpr;
		}
	if( !isBinary ) configuration_file << endl;
	
	// protein
	for(unsigned long k=0; k<ncells; k++)
		{
		double proteinpr = (double)protein[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&proteinpr), sizeof( double) );
		else
			configuration_file << " " << proteinpr;
		}
	if( !isBinary ) configuration_file << endl;

	// prot_rate
	for(unsigned long k=0; k<ncells; k++)
		{
		double prot_ratepr = (double)prot_rate[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&prot_ratepr), sizeof( double) );
		else
			configuration_file << " " << prot_ratepr;
		}
	if( !isBinary ) configuration_file << endl;

	// pRb
	for(unsigned long k=0; k<ncells; k++)
		{
		double pRbpr = (double)pRb[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&pRbpr), sizeof( double) );
		else
			configuration_file << " " << pRbpr;
		}
	if( !isBinary ) configuration_file << endl;

	// ConcS
	for(unsigned long k=0; k<ncells; k++)
		{
		double ConcSpr = (double)ConcS[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&ConcSpr), sizeof( double) );
		else
			configuration_file << " " << ConcSpr;
		}
	if( !isBinary ) configuration_file << endl;

	// DNA_rate
	for(unsigned long k=0; k<ncells; k++)
		{
		double DNA_ratepr = (double)DNA_rate[k];
		if( isBinary )
			configuration_b_file.write( (char*)(&DNA_ratepr), sizeof( double) );
		else
			configuration_file << " " << DNA_ratepr;
		}
	if( !isBinary ) configuration_file << endl;

	// age
	for(unsigned long k=0; k<ncells; k++)
		{
		if( isBinary )
			configuration_b_file.write( (char*)(&age[k]), sizeof(float) );
		else
			configuration_file << " " << age[k];
		}
	if( !isBinary ) configuration_file << endl;

	// phase age
	for(unsigned long k=0; k<ncells; k++)
		{
		if( isBinary )
			configuration_b_file.write( (char*)(&phase_age[k]), sizeof(float) );
		else
			configuration_file << " " << phase_age[k];
		}
	if( !isBinary ) configuration_file << endl;
	

/*
	// CO2
	for(unsigned long k=0; k<ncell; k++)
		{
		configuration_file << " " << CellVector[k].Get_CO2()/CellVector[k].Get_volume();
		}
	configuration_file << endl;

	// AcLRate
	for(unsigned long k=0; k<ncell; k++)
		{
		configuration_file << " " << CellVector[k].Get_AcLRate();
		}
	configuration_file << endl;

	// volume extracellulare
	for(unsigned long k=0; k<ncell; k++)
		{
		configuration_file << " " << CellVector[k].Get_volume_extra();
		}
	configuration_file << endl;
	
*/

}

//  ******************** PrintEnv ********************

// questo metodo stampa i dati ambientali essenziali
void CellsSystem::PrintEnv(bool isBinary)
{
	double EnvVol = (double)Env.GetEnvironmentvolume(); 
	double concG = (double)(Env.GetEnvironmentG()/Env.GetEnvironmentvolume());
	double concO2 = (double)(Env.GetEnvironmentO2()/Env.GetEnvironmentvolume());
	double concA = (double)(Env.GetEnvironmentA()/Env.GetEnvironmentvolume());
	double concAcL = (double)(Env.GetEnvironmentAcL()/Env.GetEnvironmentvolume());
	
	if( isBinary )
		{
		configuration_b_file.write((char*)(&EnvVol), sizeof( double) );
		configuration_b_file.write((char*)(&concG), sizeof( double) );
		configuration_b_file.write((char*)(&concO2), sizeof( double) );
		configuration_b_file.write((char*)(&concA), sizeof( double) );
		configuration_b_file.write((char*)(&concAcL), sizeof( double) );
		}
	else
		{
		configuration_file << EnvVol << " " << concG << " " << concO2 << " " << concA << " " << concAcL << endl;
		}
	
}

//  ******************** StepStat ********************
//
// funzione che calcola le statistiche ad ogni passo
//
// 
void CellsSystem::StepStat( bool reset_stat )
{
	
	if(reset_stat)
		{
		ncalls = 0.;
		nrepeats_average = 0.;
		nrepeats_max = 0;
		nrepeats_min = MAXREPEATS;
		
		for(int k=0; k<NCONV_TEST; k++)
			convergence_fail[k] = 0;
		
		loop_count_average = 0;
		loop_count_max = 0;
		loop_count_min = 0;
		
		n_mitoses_average = 0.;
		n_mitoses_max = 0;
		n_mitoses_min = ncells;

		min_Gin = 1.;
		max_Gin = -1;
		min_Gext = 1.;
		max_Gext = -1;
		min_O2 = 1.;
		max_O2 = -1;
		min_Ain = 1.;
		max_Ain = -1;
		min_Aext = 1.;
		max_Aext = -1;
		min_AcLin = 1.;
		max_AcLin = -1;
		min_AcLext = 1.;
		max_AcLext = -1;
		min_extvolume = 1.;
		max_extvolume = -1;

		}
	else 
		{
		nrepeats_average = ( (ncalls-1)*nrepeats_average + nrepeats )/ncalls;
		if( nrepeats > nrepeats_max ) nrepeats_max = nrepeats;
		if( nrepeats < nrepeats_min ) nrepeats_min = nrepeats;
		
		loop_count_average = ( (ncalls-1)*loop_count_average + loop_count )/ncalls;
		if( loop_count > loop_count_max ) loop_count_max = loop_count;
		if( loop_count > 0 && loop_count < loop_count_min ) 
			loop_count_min = loop_count;	
		else if( loop_count == 0 )
			loop_count_min = loop_count;
			
		n_mitoses_average = ( (ncalls-1)*n_mitoses_average + n_mitoses )/ncalls;
		if( n_mitoses > n_mitoses_max ) n_mitoses_max = n_mitoses;
		if( n_mitoses < n_mitoses_min ) n_mitoses_min = n_mitoses;
			
		for(unsigned long n=0; n<ncells; n++)
			{
			if( G[n] < min_Gin ) min_Gin = G[n];
			if( G[n] > max_Gin ) max_Gin = G[n];

			if( G_extra[n] < min_Gext ) min_Gext = G_extra[n];
			if( G_extra[n] > max_Gext ) max_Gext = G_extra[n];

			if( O2[n] < min_O2 ) min_O2 = O2[n];
			if( O2[n] > max_O2 ) max_O2 = O2[n];

			if( A[n] < min_Ain ) min_Ain = A[n];
			if( A[n] > max_Ain ) max_Ain = A[n];

			if( A_extra[n] < min_Aext ) min_Aext = A_extra[n];
			if( A_extra[n] > max_Aext ) max_Aext = A_extra[n];

			if( AcL[n] < min_AcLin ) min_AcLin = AcL[n];
			if( AcL[n] > max_AcLin ) max_AcLin = AcL[n];

			if( AcL_extra[n] < min_AcLext ) min_AcLext = AcL_extra[n];
			if( AcL_extra[n] > max_AcLext ) max_AcLext = AcL_extra[n];

			if( volume_extra[n] < min_extvolume ) min_extvolume = volume_extra[n];
			if( volume_extra[n] > max_extvolume ) max_extvolume = volume_extra[n];
			}
		
		}

}


