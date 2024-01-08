#include "dependentVariables.hpp"
#include <cmath>

using namespace std;

void dependentVariables(solverConfig &cfg , mesh &msh , variables &v , matrix& mat_ns)
{
    vector<flow_float>& Ux = v.c["Ux"];
    vector<flow_float>& Uy = v.c["Uy"];
    vector<flow_float>& Uz = v.c["Uz"];

    vector<flow_float>& P = v.c["P"];
    vector<flow_float>& T = v.c["T"];

    vector<flow_float>& ro = v.c["ro"];
    vector<flow_float>& roUx = v.c["roUx"];
    vector<flow_float>& roUy = v.c["roUy"];
    vector<flow_float>& roUz = v.c["roUz"];
    vector<flow_float>& roe  = v.c["roe"];
    vector<flow_float>& sonic= v.c["sonic"];
    vector<flow_float>& Ht= v.c["Ht"];

    vector<flow_float>& ros= v.c["ros"];

    flow_float ds;
    flow_float ga = cfg.gamma;

    flow_float ros_sum = 0.0;
    flow_float rok_sum = 0.0;
    flow_float ek;

    flow_float intE;

    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        Ux[ic] = roUx[ic]/ro[ic];
        Uy[ic] = roUy[ic]/ro[ic];
        Uz[ic] = roUz[ic]/ro[ic];

        ek = 0.5*(Ux[ic]*Ux[ic] +Uy[ic]*Uy[ic] +Uz[ic]*Uz[ic]);
        intE =(roe[ic]/ro[ic] -ek);
        T[ic] =intE/(cfg.cp/cfg.gamma);

        P[ic] =(cfg.gamma-1.0)*(roe[ic]-ro[ic]*ek);
        Ht[ic] = roe[ic]/ro[ic] + P[ic]/ro[ic];

        sonic[ic] = sqrt(cfg.gamma*P[ic]/ro[ic]);

        ds = (cfg.cp/cfg.gamma)*log(T[ic]/cfg.Tref) - ((ga-1.0)/ga)*log(P[ic]/cfg.Pref);
        ros[ic] = ro[ic]*ds;

        ros_sum += ros[ic]*msh.cells[ic].volume;
        rok_sum += ro[ic]*ek*msh.cells[ic].volume;
    }

    std::cout << "ros total : " << ros_sum << std::endl;
    std::cout << "rok total : " << rok_sum << std::endl;
    
}
