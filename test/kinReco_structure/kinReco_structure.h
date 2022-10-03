
#ifndef TTBAR_KINRECO_H
#define TTBAR_KINRECO_H

// C++ library or ROOT header files
#include <TMath.h>
#include <Math/Polynomial.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <vector>
#include <complex>

// debugging level (0 for silence, > 0 for some messages)
int gDebug = 0;

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>> ZSolutionKinRecoDilepton >>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// structure to store one solution of the kinematic reconstruction
struct ZSolutionKinRecoDilepton
{
  // constructor
  // (set weight to -1 by default)
  ZSolutionKinRecoDilepton(): zWeight(-1.0) {;}
  TLorentzVector zT, zTbar;
  int zBTag;
  double zWeight;
};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ZSolutionKinRecoDilepton* SolveKinRecoDilepton( const TLorentzVector& lm, 
                                                const TLorentzVector& lp, 
                                                const TLorentzVector& b, 
                                                const TLorentzVector& bbar, 
                                                const double metX, 
                                                const double metY, 
                                                TH1D* hInacc = NULL, 
                                                int* ambiguity = NULL)
{
    //...

    // store and return best solution as ZSolutionKinRecoDilepton instance
    ZSolutionKinRecoDilepton* solution = new ZSolutionKinRecoDilepton;

    solution->zT = (nuBest + lp + b);
    solution->zTbar = (nubarBest + lm + bbar);
    solution->zWeight = weightBest;
    return solution;
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int KinRecoDilepton(const TLorentzVector& lm, 
                    const TLorentzVector& mp, 
                    const std::vector<TLorentzVector>& jets, 
                    const double metX, 
                    const double metY, 
                    TLorentzVector& t, 
                    TLorentzVector& tbar, 
                    TH1D* hInacc = NULL, 
                    TH1D* hAmbig = NULL)
{
    int solved = 0;   
    //...

    // container with solutions
    std::vector<ZSolutionKinRecoDilepton*> vSolutions;

    // loop jets
        ZSolutionKinRecoDilepton* solution = SolveKinRecoDilepton(lm, mp, jetB, jetBbar, metX, metY, hInacc, ambiguity);
        vSolutions.push_back(solution);

    for(std::vector<ZSolutionKinRecoDilepton*>::iterator it = vSolutions.begin(); it != vSolutions.end(); it++)
    {
        //...
        solved = 1;
    }
    
    return solved
}


#endif