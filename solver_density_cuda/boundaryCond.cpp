#include "boundaryCond.hpp"

#include "cuda_nagare/cudaWrapper.cuh"
#include "cuda_nagare/boundaryCond_d.cuh"

using namespace std;

bcondConfFormat::bcondConfFormat(){};

void readBcondConfig(solverConfig& cfg , vector<bcond>& bconds)
{
    map<int,bcondConfFormat> bcondConfMap;

    std::string bcondConfigFileName  = "bcondConfig.yaml";

    try {
        YAML::Node config = YAML::LoadFile(bcondConfigFileName);

        for(auto bcInYaml : config) {
            bcondConfFormat bcf;

            std::string bname = bcInYaml.first.as<std::string>();
            cout << "bname = " << bname << endl;

            // akirameta
            //if (bcInYaml.find("physID") == bcInYaml.end()) {
            //    // not found
            //    cerr << "Error: couldn't find physID of " << bname << " in bcondConfig.yaml\n";
            //    exit();
            //}

            int physID = bcInYaml.second["physID"].as<int>();

            std::map<std::string, int> inputInts_temp;
            for (auto bcInYaml_int : bcInYaml.second["ints"])
            {
                std::string key = bcInYaml_int.first.as<std::string>();
                int val = bcInYaml_int.second.as<int>();
                inputInts_temp[key] = val;
            }

            std::map<std::string, flow_float> inputFloats_temp;
            for (auto bcInYaml_float : bcInYaml.second["floats"])
            {
                std::string key = bcInYaml_float.first.as<std::string>();
                flow_float val = bcInYaml_float.second.as<flow_float>();
                inputFloats_temp[key] = val;
            }

            std::string kind = bcInYaml.second["kind"].as<std::string>();


            bcf.physID = physID;
            bcf.kind = kind;
            bcf.inputInts = inputInts_temp;
            bcf.inputFloats = inputFloats_temp;

            bcondConfMap[bcf.physID] = bcf;

            cout << "physID=" <<  bcf.physID << endl;
            cout << "kind=" << bcf.kind<< endl;
        }

    } catch(const YAML::BadFile& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);

    } catch(const YAML::ParserException& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);
    }

    // set bconds from cofig file (YAML file)
    for (bcond& bc : bconds)
    {
        cout << "bc.physID=" << bc.physID << endl;

        bcondConfFormat bcf;
        if (bcondConfMap.find(bc.physID) == bcondConfMap.end()) {
            // not found
            cerr << "Error: couldn't find physID " << bc.physID << " in bcondConfig.yaml\n";
            exit;
        } else {
            bcf = bcondConfMap[bc.physID];
        }

        cout << "bcf.kind=" << bcf.kind << endl;

        bc.bcondKind = bcf.kind;
        bc.inputInts = bcf.inputInts;
        bc.inputFloats= bcf.inputFloats;

        // copy value types to bcond
        map<string,int> valueTypes = bcf.valueTypesOfBC[bcf.kind];
        for (auto& vt : valueTypes)
        {
            bc.valueTypes[vt.first] = vt.second;
        }

        bc.bcondInitVariables(cfg.gpu); // allocate and set boundary variables
    }
};

void wall_isothermal(solverConfig& cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p )
{
//    //IntName   : none
//    //FloatName : T 
//    //geom_int ib_loc = 0;
//    //for (auto ib : bc.iBPlanes)
//    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
//    {
//        geom_int ip = bc.iPlanes[ib];
//        var.p["T"][ip] = bc.bvar["T"][ib];
//    }
//
//    vector<geom_float> pcent(3);
//    vector<geom_float> c1cent(3);
//    geom_int ic0;
//    geom_float dn;
//
//    geom_float ss;
//    vector<geom_float> sv(3);
//
//    flow_float temp;
//
//
//    // calculate diffusion term
//    for (geom_int& ip : bc.iPlanes)
//    {
//        ic0     = msh.planes[ip].iCells[0];
//        sv      = msh.planes[ip].surfVect;
//        ss      = msh.planes[ip].surfArea;
//        pcent   = msh.planes[ip].centCoords;
//    
//        c1cent  = msh.cells[ic0].centCoords;
//    
//        dn = ( (pcent[0] - c1cent[0])*sv[0]
//              +(pcent[1] - c1cent[1])*sv[1]
//              +(pcent[2] - c1cent[2])*sv[2] )/ss;
//    
//        var.c["diffT"][ic0] += (var.p["T"][ip] - var.c["T"][ic0])/dn*ss;
//
//        mat_p.lhs[ic0][0] += ss/dn;
//        mat_p.rhs[ic0] += var.p["T"][ip]/dn*ss; 
//    }
//
};

void wall(solverConfig& cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p )
{

    flow_float roc;
    flow_float roec;

    flow_float Uxb;
    flow_float Uyb;
    flow_float Uzb;

    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
    {
        geom_int ip = bc.iPlanes[ib];
        geom_int ic0 = bc.iCells[ib];

        Uxb = bc.bvar["Ux"][ib];
        Uyb = bc.bvar["Uy"][ib];
        Uzb = bc.bvar["Uz"][ib];

        roc = var.c["ro"][ic0];
        roec = var.c["roe"][ic0];

        bc.bvar["ro"][ib] = roc;
        bc.bvar["Ux"][ib] = Uxb;
        bc.bvar["Uy"][ib] = Uyb;
        bc.bvar["Uz"][ib] = Uzb;
        bc.bvar["roUx"][ib] = roc*Uxb;
        bc.bvar["roUy"][ib] = roc*Uyb;
        bc.bvar["roUz"][ib] = roc*Uzb;
        bc.bvar["roe"][ib] = roec;
        bc.bvar["Ps"][ib]  = (cfg.gamma-1.0)*(roec-0.5*roc*(Uxb*Uxb+Uyb*Uyb+Uzb*Uzb));
    }
};

void inlet_uniformVelocity(solverConfig& cfg , bcond& bc , mesh& msh , variables& v , matrix& mat_p )
{
    flow_float rob;
    flow_float Uxb;
    flow_float Uyb;
    flow_float Uzb;
    flow_float Psb;

    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
    {
        geom_int ip = bc.iPlanes[ib];
        geom_int ic = bc.iCells[ib];
        geom_int ig = bc.iCells_ghst[ib];

        rob = bc.bvar["ro"][ib]; // ro is given
        Uxb = bc.bvar["Ux"][ib]; // U is given
        Uyb = bc.bvar["Uy"][ib];
        Uzb = bc.bvar["Uz"][ib];
        Psb = bc.bvar["Ps"][ib];

        v.c["ro"][ig]    = rob;
        v.c["roUx"][ig]  = rob*Uxb;
        v.c["roUy"][ig]  = rob*Uyb;
        v.c["roUz"][ig]  = rob*Uzb;
        v.c["roe"][ig]   = Psb/(cfg.gamma-1.0) + 0.5*rob*(Uxb*Uxb +Uyb*Uyb +Uzb*Uzb);
        v.c["P"][ig]     = Psb;
        v.c["Ux"][ig]    = Uxb;
        v.c["Uy"][ig]    = Uyb;
        v.c["Uz"][ig]    = Uzb;
        v.c["Ht"][ig]    = v.c["roe"][ig]/rob + Psb/rob;
        v.c["sonic"][ig] = sqrt(cfg.gamma*Psb/rob);

    }
}

void inlet_Pressure(solverConfig& cfg , bcond& bc , mesh& msh , variables& v , matrix& mat_p )
{
    geom_int ip;
    geom_int ic;

    flow_float Ux_c;
    flow_float Uy_c;
    flow_float Uz_c;
    flow_float Umag_c;

    flow_float mach_c;
    flow_float mach_b;
    flow_float mach_new;

    flow_float T_c;
    flow_float P_c;

    flow_float Pt_b;
    flow_float Tt_b;

    flow_float Ps_new;
    flow_float Ts_new;
    flow_float ro_new;
    flow_float sonic_new;

    flow_float sx;
    flow_float sy;
    flow_float sz;

    flow_float rf=0.5;

    geom_int ig;

    // specify total pressure and total temperature
    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
    {
        ip = bc.iPlanes[ib];
        ic = bc.iCells [ib];
        ig = bc.iCells_ghst[ib];

        vector<flow_float> sv = msh.planes[ip].surfVect;
        flow_float ss = msh.planes[ip].surfArea;

        Pt_b = bc.bvar["Pt"][ib];
        Tt_b = bc.bvar["Tt"][ib];

        Ux_c = v.c["Ux"][ic];
        Uy_c = v.c["Uy"][ic];
        Uz_c = v.c["Uz"][ic];

        Umag_c = sqrt(Ux_c*Ux_c + Uy_c*Uy_c + Uz_c*Uz_c);
        mach_c = Umag_c/v.c["sonic"][ic];

        T_c  = v.c["T"][ic]; //static
        P_c  = v.c["P"][ic]; //static

        flow_float ga = cfg.gamma;

        mach_b = sqrt((pow((P_c/Pt_b),-(ga-1.0)/ga) -1.0)*2.0/(ga-1.0));

        mach_new = rf*mach_b + (1.0-rf)*mach_c;
        //mach_new = mach_b ;

        Ps_new = P_c;
        Ts_new = Tt_b/(1.0+0.5*(cfg.gamma-1.0)*mach_new*mach_new);
        sonic_new = sqrt((cfg.gamma-1.0)*cfg.cp*Ts_new);
        ro_new = cfg.gamma*Ps_new/((cfg.gamma-1.0)*cfg.cp*Ts_new);

        //flow_float Uxnew = (Ux_c/Umag_c)*mach_new*sonic_new;
        //flow_float Uynew = (Uy_c/Umag_c)*mach_new*sonic_new;
        //flow_float Uznew = (Uz_c/Umag_c)*mach_new*sonic_new;

        flow_float Uxnew = -mach_new*sonic_new*sv[0]/ss;
        flow_float Uynew = -mach_new*sonic_new*sv[1]/ss;
        flow_float Uznew = -mach_new*sonic_new*sv[2]/ss;

        v.c["ro"][ig]   = ro_new;
        v.c["P"][ig]    = Ps_new;
        v.c["Ux"][ig]   = Uxnew;
        v.c["Uy"][ig]   = Uynew;
        v.c["Uz"][ig]   = Uznew;
        v.c["roUx"][ig] = ro_new*Uxnew;
        v.c["roUy"][ig] = ro_new*Uynew;
        v.c["roUz"][ig] = ro_new*Uznew;
        v.c["roe"][ig]  = Ps_new/(cfg.gamma-1.0) + 0.5*ro_new*(Uxnew*Uxnew+Uynew*Uynew+Uznew*Uznew);
        v.c["Ht"][ig]    = v.c["roe"][ig]/ro_new + Ps_new/ro_new;
        v.c["sonic"][ig] = sqrt(cfg.gamma*Ps_new/ro_new);
    }
}

void outlet_statPress(solverConfig& cfg , bcond& bc , mesh& msh , variables& v , matrix& mat_p )
{
    vector<geom_float> sv;
    geom_float ss;

    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
    {
        geom_int ip = bc.iPlanes[ib];
        geom_int ic = bc.iCells [ib];
        geom_int ig = bc.iCells_ghst[ib];

        sv = msh.planes[ip].surfVect;
        ss = msh.planes[ip].surfArea;

        flow_float Pnew = bc.bvar["Ps"][ib];
        flow_float roc = v.c["ro"][ic];
        flow_float Uxc = v.c["Ux"][ic];
        flow_float Uyc = v.c["Uy"][ic];
        flow_float Uzc = v.c["Uz"][ic];
        flow_float Umagc = sqrt(Uxc*Uxc+Uyc*Uyc+Uzc*Uzc);
        flow_float roec = v.c["roe"][ic];
        flow_float Umag_new;

        if (Uxc*sv[0]+Uyc*sv[1]+Uzc*sv[2] > 0) {
            Umag_new = sqrt((roec-Pnew/(cfg.gamma-1.0))*2.0/roc);
        } else {
            Umag_new = 0.0;
        }

        //cout << "bvar ps=" << Pnew << endl;
        //cout << "bvar roc=" << roc << endl;
        //cout << "bvar Umag_new=" << Umag_new << endl;
        //cout << "bvar Umag_c=" << Umagc << endl;
        //cout << "bvar roec=" << roec << endl;

        v.c["ro"][ig]   = roc;
        v.c["P"][ig]    = Pnew; // already set
        v.c["Ux"][ig]   = Uxc*Umag_new/Umagc;
        v.c["Uy"][ig]   = Uyc*Umag_new/Umagc;
        v.c["Uz"][ig]   = Uzc*Umag_new/Umagc;
        v.c["roUx"][ig] = roc*Uxc*Umag_new/Umagc;
        v.c["roUy"][ig] = roc*Uyc*Umag_new/Umagc;
        v.c["roUz"][ig] = roc*Uzc*Umag_new/Umagc;
        v.c["roe"][ig]  = roec;
        v.c["Ht"][ig]    = v.c["roe"][ig]/roc + Pnew/roc;
        v.c["sonic"][ig] = sqrt(cfg.gamma*Pnew/roc);
    }
}

void outflow(solverConfig& cfg , bcond& bc , mesh& msh , variables& v , matrix& mat_p )
{
    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
    {
        geom_int ip = bc.iPlanes[ib];
        geom_int ic = bc.iCells[ib];
        geom_int ig = bc.iCells_ghst[ib];

        //bc.bvar["ro"][ib]   = v.c["ro"][ic];
        //bc.bvar["Ps"][ib]   = v.c["P"][ic];
        //bc.bvar["Ux"][ib]   = v.c["Ux"][ic];
        //bc.bvar["Uy"][ib]   = v.c["Uy"][ic];
        //bc.bvar["Uz"][ib]   = v.c["Uz"][ic];
        //bc.bvar["roUx"][ib] = v.c["roUx"][ic];
        //bc.bvar["roUy"][ib] = v.c["roUy"][ic];
        //bc.bvar["roUz"][ib] = v.c["roUz"][ic];
        //bc.bvar["roe"][ib]  = v.c["roe"][ic];
        v.c["ro"][ig]   = v.c["ro"][ic];
        v.c["P"][ig]    = v.c["P"][ic];
        v.c["Ux"][ig]   = v.c["Ux"][ic];
        v.c["Uy"][ig]   = v.c["Uy"][ic];
        v.c["Uz"][ig]   = v.c["Uz"][ic];
        v.c["roUx"][ig] = v.c["roUx"][ic];
        v.c["roUy"][ig] = v.c["roUy"][ic];
        v.c["roUz"][ig] = v.c["roUz"][ic];
        v.c["roe"][ig]  = v.c["roe"][ic];
        v.c["Ht"][ig]   = (v.c["roe"][ig] + v.c["P"][ig])/v.c["ro"][ig];
        v.c["sonic"][ig]= sqrt(cfg.gamma*v.c["P"][ig]/v.c["ro"][ig]);

    }
}

void slip(solverConfig& cfg , bcond& bc , mesh& msh , variables& v , matrix& mat_p )
{
    vector<flow_float>& Ux  = v.c["Ux"];
    vector<flow_float>& Uy  = v.c["Uy"];
    vector<flow_float>& Uz  = v.c["Uz"];
    vector<flow_float>& ro  = v.c["ro"];
    vector<flow_float>& P   = v.c["P"];
    vector<flow_float>& roe = v.c["roe"];

    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
    {
        geom_int ip = bc.iPlanes[ib];
        geom_int ic = bc.iCells[ib];
        geom_int ig = bc.iCells_ghst[ib];

        flow_float Un;
        vector<flow_float> sv = msh.planes[ip].surfVect;
        flow_float ss = msh.planes[ip].surfArea;

        Un =  (sv[0]*Ux[ic] + sv[1]*Uy[ic] + sv[2]*Uz[ic])/ss;

        v.c["ro"][ig]   = ro[ic];
        v.c["P"][ig]    = P[ic];
        v.c["Ux"][ig]   = Ux[ic]- 2*Un*sv[0]/ss;
        v.c["Uy"][ig]   = Uy[ic]- 2*Un*sv[1]/ss;
        v.c["Uz"][ig]   = Uz[ic]- 2*Un*sv[2]/ss;
        v.c["roUx"][ig] = ro[ic]*(Ux[ic]- 2*Un*sv[0]/ss);
        v.c["roUy"][ig] = ro[ic]*(Uy[ic]- 2*Un*sv[1]/ss);
        v.c["roUz"][ig] = ro[ic]*(Uz[ic]- 2*Un*sv[2]/ss);
        v.c["roe"][ig]  = P[ic]/(cfg.gamma-1.0) + 0.5*ro[ic]*(Ux[ic]*Ux[ic]+Uy[ic]*Uy[ic]+Uz[ic]*Uz[ic]);

        v.c["Ht"][ig]   = (v.c["roe"][ig] + v.c["P"][ig])/v.c["ro"][ig];
        v.c["sonic"][ig]= sqrt(cfg.gamma*v.c["P"][ig]/v.c["ro"][ig]);

        //cout << "in slip " << ic << " " << ig << " " << endl; 
        //cout << Un << " " << Ux[ic] << " " << Uy[ic] << " " << Uz[ic]<< endl;
        //cout << sv[0] << " " << sv[1] << " " << sv[2]<< endl;
        //cout << "Ux ig=" << v.c["Ux"][ig]  << endl; 
    }
}

void periodic(solverConfig& cfg , bcond& bc , mesh& msh , variables& v , matrix& mat_p )
{

    vector<flow_float>& ro    = v.c["ro"];
    vector<flow_float>& roUx  = v.c["roUx"];
    vector<flow_float>& roUy  = v.c["roUy"];
    vector<flow_float>& roUz  = v.c["roUz"];
    vector<flow_float>& roe   = v.c["roe"];

    vector<flow_float>& Ux    = v.c["Ux"];
    vector<flow_float>& Uy    = v.c["Uy"];
    vector<flow_float>& Uz    = v.c["Uz"];
    vector<flow_float>& P     = v.c["P"];

    vector<flow_float>& Ht    = v.c["Ht"];
    vector<flow_float>& sonic = v.c["sonic"];

    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
    {
        geom_int ip         = bc.iPlanes[ib];
        geom_int ip_partner = bc.bint["partnerPlnID"][ib];
        geom_int ic = bc.iCells[ib];
        geom_int ic_partner = msh.planes[ip_partner].iCells[0];
        geom_int ig = bc.iCells_ghst[ib];

        //cout << "ic = "<<  ic << " ic_p=" << ic_partner  << endl;

        ro[ig]   = ro[ic_partner];
        P[ig]    = P[ic_partner];
        Ux[ig]   = Ux[ic_partner];
        Uy[ig]   = Uy[ic_partner];
        Uz[ig]   = Uz[ic_partner];
        roUx[ig] = roUx[ic_partner];
        roUy[ig] = roUy[ic_partner];
        roUz[ig] = roUz[ic_partner];
        roe[ig]  = roe[ic_partner];

        Ht[ig]   = Ht[ic_partner];
        sonic[ig]= sonic[ic_partner];
    }
}

void applyBconds(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_p)
{
    if (cfg.gpu == 0) { // cpu
        for (auto& bc : msh.bconds)
        {
            if      (bc.bcondKind == "wall_isothermal") { wall_isothermal(cfg , bc , msh , var , mat_p); } 
            else if (bc.bcondKind == "slip") { slip(cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "wall") { wall(cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "inlet_uniformVelocity") { inlet_uniformVelocity(cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "inlet_Pressure") { inlet_Pressure(cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "outlet_statPress") { outlet_statPress(cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "outflow") { outflow(cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "periodic") { periodic(cfg , bc , msh , var , mat_p); }
            else {
                cerr << "Error: unknown bcondKind " << bc.bcondKind << endl;
                exit(EXIT_FAILURE);
            };
        }
    } else if (cfg.gpu ==1) { // gpu
        for (auto& bc : msh.bconds)
        {
//temp            //if      (bc.bcondKind == "wall_isothermal") { wall_isothermal_d(cfg , bc , msh , var , mat_p); } 
//temp            if      (bc.bcondKind == "slip") slip_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); 
//temp            else if (bc.bcondKind == "wall") wall_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); 
//temp            else if (bc.bcondKind == "inlet_uniformVelocity") { inlet_uniformVelocity_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
//temp            else if (bc.bcondKind == "outlet_statPress") { outlet_statPress_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
//temp            //else {
            //    cerr << "Error: unknown bcondKind " << bc.bcondKind << endl;
            //    exit(EXIT_FAILURE);
            //};
            cudaThreadSynchronize();
        }
//    } else if (cfg.gpu ==-1) { // gpu check
//        for (auto& bc : msh.bconds)
//        {
////            //if      (bc.bcondKind == "wall_isothermal") { wall_isothermal_d(cfg , bc , msh , var , mat_p); } 
////            if      (bc.bcondKind == "slip") slip_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); 
////            else if (bc.bcondKind == "wall") wall_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); 
////            else if (bc.bcondKind == "inlet_uniformVelocity") { inlet_uniformVelocity_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
////            else if (bc.bcondKind == "outlet_statPress") { outlet_statPress_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
////            //else {
////            //    cerr << "Error: unknown bcondKind " << bc.bcondKind << endl;
////            //    exit(EXIT_FAILURE);
////            //};
////            cudaThreadSynchronize();
//            if      (bc.bcondKind == "wall_isothermal") { wall_isothermal(cfg , bc , msh , var , mat_p); } 
//            else if (bc.bcondKind == "slip") { slip(cfg , bc , msh , var , mat_p); }
//            else if (bc.bcondKind == "wall") { wall(cfg , bc , msh , var , mat_p); }
//            else if (bc.bcondKind == "inlet_uniformVelocity") { inlet_uniformVelocity(cfg , bc , msh , var , mat_p); }
//            else if (bc.bcondKind == "outlet_statPress") { outlet_statPress(cfg , bc , msh , var , mat_p); }
//            else {
//                cerr << "Error: unknown bcondKind " << bc.bcondKind << endl;
//                exit(EXIT_FAILURE);
//            };
//
//        }
    }
}

