#ifndef BLOODVESSEL_H
#define BLOODVESSEL_H  // header guard
//
//  BloodVessel.h
//  Sim3D-v3-S
//
//  Created by Edoardo Milotti on 25/7/2017.
//
//


class BloodVessel
{
    // friends
    
    // friend class CellsSystem;
    
    // printout function
    friend ostream& operator<<(ostream& s, BloodVessel& cBloodVessel);
    
private:
    
    vector<double> a;           // starting position (start)
    vector<double> b;           // end position (end)
    vector<double> va;          // velocità della posizione di inizio
    vector<double> vb;          // velocità della posizione di fine
    double R;                   // raggio
    double vR;                  // velocità del raggio
    double O2start;				// O2 concentration at start
    double O2end;				// O2 concentration at end
    double CO2start;			// CO2 concentration at start
    double CO2end;				// CO2 concentration at end
    double G;					// glucose concentration
    double A;					// other nutrients concentration 
    double AcL;					// lactate concentration 
    
public:
    
    // costruttori
    // costruttore di default (nessun attributo specificato)
    BloodVessel();
    // costruttore con valori specificati
    BloodVessel(const vector<double> ca, const vector<double> cb, const vector<double> cva, const vector<double> cvb, const double cR, const double cvR, const double cO2start, const double cO2end, const double cCO2start, const double cCO2end, const double cG, const double cA, const double cAcL);
    // costruttore copia
    BloodVessel(const BloodVessel& bv);
    
    // nessun distruttore, utilizzo il default del compilatore
    
    // overloaded =
    BloodVessel& operator=(const BloodVessel& bv);
    
    // setters
    void SetBloodVesselR(const double newR) { R = newR; };
    void SetBloodVesselvR(const double newvR) { vR = newvR; };
    void SetBloodVessela(const vector<double> newa) { a = newa; };
    void SetBloodVesselb(const vector<double> newb) { b = newb; };
    void SetBloodVesselva(const vector<double> newva) { va = newva; };
    void SetBloodVesselvb(const vector<double> newvb) { vb = newvb; };
    void SetBloodVesselak( const double newak, const int k ) { a[k] = newak; };
    void SetBloodVesselbk( const double newbk, const int k ) { b[k] = newbk; };
    void SetBloodVesselvak( const double newvak, const int k ) { va[k] = newvak; };
    void SetBloodVesselvbk( const double newvbk, const int k ) { vb[k] = newvbk; };
    void SetBloodVesselO2start(const double newO2start) { O2start = newO2start; };
    void SetBloodVesselO2end(const double newO2end) { O2end = newO2end; };
    void SetBloodVesselCO2start(const double newCO2start) { CO2start = newCO2start; };
    void SetBloodVesselCO2end(const double newCO2end) { CO2end = newCO2end; };
    void SetBloodVesselG(const double newG) { G = newG; };
    void SetBloodVesselA(const double newA) { A = newA; };
    void SetBloodVesselAcL(const double newAcL) { AcL = newAcL; };
    
    // getters
    vector<double> GetBloodVessela() { return a; };
    vector<double> GetBloodVesselb() { return b; };
    vector<double> GetBloodVesselva() { return va; };
    vector<double> GetBloodVesselvb() { return vb; };
    double GetBloodVesselR() { return R; };
    double GetBloodVesselvR() { return vR; };
    double GetBloodVesselO2start() { return O2start; };
    double GetBloodVesselO2end() { return O2end; };
    double GetBloodVesselCO2start() { return CO2start; };
    double GetBloodVesselCO2end() { return CO2end; };
    double GetBloodVesselG() { return G; };
    double GetBloodVesselA() { return A; };
    double GetBloodVesselAcL() { return AcL; };
    
    
    // altri metodi
    double GetBloodVesselLength() { return sqrt( (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]) ); };
    double GetBloodVesselVolume() { return M_PI * R * R * sqrt( (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]) ); };
    double DistanceFromVessel( const vector<double> y, double* x0 );
    
    // read and write binario
        
    
};

#endif //#ifndef BLOODVESSEL_H


