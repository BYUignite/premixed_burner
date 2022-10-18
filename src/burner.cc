
#include "cantera/onedim.h"
#include "cantera/base/stringUtils.h"
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include "sootModel.h"

using namespace Cantera;
using namespace soot;

using fmt::print;
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using std::setw;
using std::scientific;
using std::setprecision;
using std::cout;
using std::endl;
using std::to_string;

//////////////////////////////////////////////////////////////////////////////////////////

int main() {

    bool dosoot = true;                                // change just this to turn soot on and off

    nucleationMech  n = nucleationMech::LIN;            // Nucleation: NONE, LL, LIN, PAH
    growthMech      g = growthMech::LIN;                // Surface growth: NONE, LL, LIN, HACA
    oxidationMech   x = oxidationMech::LL;              // Oxidation: NONE, LL, LEE_NEOH, NSC_NEOH, HACA
    coagulationMech c = coagulationMech::FM;            // Coagulation: NONE, LL, FUCHS, FRENK
    psdMech         PSD = psdMech::MONO;                // PSD mechanisms: MONO, LOGN, QMOM, MOMIC
    int             nsoot = 2;                          // number of soot moments
    sootModel       SM = sootModel(PSD, nsoot, n, g, x, c);
    state           SS = state(nsoot);
    vector<double>  sootScales{1E16, 0.01};
    SS.setSootScales(sootScales);

    ///////////// setup

    int loglevel = 1;
    bool refine_grid = true;

    double Tin     = 300.0;                       // K
    double P       = 101325.0;                    // Pa
    string Xstr    = "C2H4:2.34, O2:3, N2:11.28";
    double uin     = 0.0673;                      // m/sec
    double Ldomain = 0.03;                        // domain size (m)

    auto sol = newSolution("gri30.yaml", "gri30", "None");
    auto gas = sol->thermo();

    int nsp  = gas->nSpecies();

    gas->setState_TPX(Tin, P, Xstr);
    vector<double> X(nsp);
    gas->getMoleFractions(&X[0]);
    double mdot   = uin * gas->density();

    vector<double> yin(nsp);
    gas->getMassFractions(&yin[0]);

    gas->equilibrate("HP");
    vector<double> yout(nsp);
    gas->getMassFractions(&yout[0]);
    double rho_out = gas->density();
    double Tad = gas->temperature();


    ///////////// build each domain

    //----------- create the inlet

    Inlet1D inlet;

    inlet.setMoleFractions(&X[0]);
    inlet.setMdot(mdot);
    inlet.setTemperature(Tin);
    //inlet.setSpreadRate(10.);

    //----------- create the flow

    int nscant = dosoot ? nsoot : 0;
    StFlow flow(gas, 1, 1, nscant, SM, SS);   // soot
    flow.setAxisymmetricFlow();

    //----------- initial grid

    vector<double> z{0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0};
    for(int i=0; i<z.size(); i++)
        z[i] *= Ldomain;
    flow.setupGrid(z.size(), &z[0]);

    //----------- temperature profile

    ifstream ifile("zT.dat");
    vector_fp zTz, zTT;
    double a, b;
    for(;;){
        ifile >> a >> b;          // z (m), T (K)
        if(ifile.eof()) break;
        zTz.push_back(a);
        zTT.push_back(b);
    }
    ifile.close();
    for(int i=0; i<zTz.size(); i++)
        zTz[i] /= zTz.back();

    flow.setFixedTempProfile(zTz, zTT);

    //----------- set transport and kinetics objects

    sol->setTransport("mixture-averaged");
    flow.setTransport(*sol->transport());
    flow.setKinetics(*sol->kinetics());
    flow.setPressure(P);

    //----------- create the outlet

    Outlet1D outlet;

    ///////////// create the container and insert the domains

    std::vector<Domain1D*> domains { &inlet, &flow, &outlet };

    Sim1D flame(domains);

    //----------- Supply initial guess, flame parameters

    vector_fp locs{0.0, 0.3, 0.7, 1.0};
    vector_fp value;

    double uout = inlet.mdot()/rho_out;
    value = {uin, uin, uout, uout};
    flame.setInitialGuess("velocity",locs,value);
    value = {Tin, Tin, Tad, Tad};
    flame.setInitialGuess("T",locs,value);

    for (size_t i=0; i<nsp; i++) {
        value = {yin[i], yin[i], yout[i], yout[i]};
        flame.setInitialGuess(gas->speciesName(i),locs,value);
    }

    double ndens0 = 1;                    // scaled values (code solves and outputs scaled Mhat0, Mhat1)
    double rhoYs0 = 1;                    // the scales are hard coded in Cantera StFlow.cpp: Mhat0scale=1E16, Mhat1scale = 0.01. 
    value = {0.0, 0.0, ndens0, ndens0};   // so, to get fv: Mhat1 * 0.01 * rho / rhosoot
    flame.setInitialGuess("soot_var_0",locs,value);
    value = {0.0, 0.0, rhoYs0, rhoYs0};
    flame.setInitialGuess("soot_var_1",locs,value);

    inlet.setMoleFractions(&X[0]);
    inlet.setMdot(mdot);
    inlet.setTemperature(Tin);

    flame.showSolution();

    int flowdomain = 1;
    double ratio = 3.0; // 10.0;
    double slope = 0.3; // 0.08;
    double curve = 1.0; // 0.1;

    flame.setRefineCriteria(flowdomain,ratio,slope,curve);

    ///////////// solve the problem

    flame.solve(loglevel,refine_grid);

    ///////////// recover soot from flame 

    vector_fp rhovec;
    vector<vector_fp> soot(nsoot);

    for (size_t n = 0; n < flow.nPoints(); n++) {
        double T = flame.value(flowdomain,flow.componentIndex("T"),n);
        vector<double> y(nsp);
        for(size_t k=0; k<nsp; k++)
            y[k] = flame.value(flowdomain, flow.componentIndex(gas->speciesName(k)), n);
        gas->setState_TPX(T, 101325, &y[0]);
        rhovec.push_back(gas->density());
        if(dosoot)
            for(size_t k=0; k<nsoot; k++)
                soot[k].push_back(flame.value(flowdomain, flow.componentIndex("soot_var_"+to_string(k)), n)*gas->density()*SS.sootScales[k]);
    }

    ///////////// output


    ofstream ofile("burner.out");

    ofile << "# ";
    ofile << setw(18) << "z_(m) ";
    ofile << setw(19) << "V_(m/s) ";
    ofile << setw(19) << "T_(K) ";
    ofile << setw(19) << "rho_(kg/m3) ";
    ofile << setw(19) << "Y_CO ";
    ofile << setw(19) << "Y_CO2 ";
    ofile << setw(19) << "Y_C2H2 ";
    if(dosoot)
        for(size_t k=0; k<nsoot; k++)
            ofile << setw(19) << "S_" + to_string(k) + " ";

    ofile << scientific;
    ofile << setprecision(10);
    for(size_t n = 0; n<flow.nPoints(); n++) {
        ofile << endl 
              << setw(19) << flow.grid(n) 
              << setw(19) << flame.value(flowdomain, flow.componentIndex("velocity"),n)
              << setw(19) << flame.value(flowdomain, flow.componentIndex("T"),n)
              << setw(19) << rhovec[n]
              << setw(19) << flame.value(flowdomain, flow.componentIndex("CO"),n)
              << setw(19) << flame.value(flowdomain, flow.componentIndex("CO2"),n)
              << setw(19) << flame.value(flowdomain, flow.componentIndex("C2H2"),n);
        if(dosoot)
            for(size_t k=0; k<nsoot; k++)
                ofile << setw(19) << soot[k][n];
    }
    ofile.close();

    flame.showSolution();

    return 0;
}

