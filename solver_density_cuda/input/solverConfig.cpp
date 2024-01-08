#include "input/solverConfig.hpp"


solverConfig::solverConfig(){};

void solverConfig::read(std::string fname)
{
    try {
        this->solConfigFileName = fname;
        YAML::Node config = YAML::LoadFile(this->solConfigFileName);

        std::string meshFormat   = config["mesh"]["meshFormat"].as<std::string>();
        std::string meshFileName = config["mesh"]["meshFileName"].as<std::string>();
        std::string solver = config["solver"].as<std::string>();

        this->gpu = config["gpu"].as<int>();

        auto last = config["time"]["last"];
            int endTimeControl = last["control"].as<int>();
            int nStep = last["nStep"].as<int>();
            int endTime = last["time"].as<flow_float>();

        auto deltaT = config["time"]["deltaT"];
            int dtControl = deltaT["control"].as<int>();
            flow_float dt = deltaT["dt"].as<flow_float>();
            flow_float cfl = deltaT["cfl"].as<flow_float>();
            flow_float dt_max= deltaT["dt_max"].as<flow_float>();
            flow_float dt_min= deltaT["dt_min"].as<flow_float>();

        int outStepInterval = config["time"]["outStepInterval"].as<int>();
        int timeIntegration = config["time"]["timeIntegration"].as<int>();

        auto implicit = config["time"]["implicit"];
            int nLoop      = implicit["nLoop"].as<int>();
            flow_float cfl_pseudo = implicit["cfl_pseudo"].as<flow_float>();

        auto space = config["space"];
            int convMethod = space["convMethod"].as<int>();


        auto physProp = config["physProp"];
            int isCompressible = physProp["isCompressible"].as<int>();
            flow_float ro = physProp["ro"].as<flow_float>();
            flow_float visc = physProp["visc"].as<flow_float>();
            flow_float thermCond = physProp["thermCond"].as<flow_float>();
            flow_float cp = physProp["cp"].as<flow_float>();
            flow_float gamma = physProp["gamma"].as<flow_float>();


        this->initial = config["initial"].as<std::string>();

        std::cout << "Mesh Name : " << meshFileName << '\n';
        std::cout << "Mesh Name size : " << meshFileName.size() << '\n';

        std::cout << "Solver Name : " << solver << '\n';
        if        (endTimeControl == 0) { std::cout << "End Step : " << nStep << '\n';
        } else if (endTimeControl == 1) { std::cout << "End Time : " << endTime << '\n';
        } else {
            std::cerr << "something wrong in solver yaml" << std::endl;
            exit(EXIT_FAILURE);
        }

        if        (dtControl == 0) { std::cout << "delta Time : " << dt << '\n'; 
        } else if (dtControl == 1) { std::cout << "CFL : " << cfl << '\n';
        } else {
            std::cerr << "something wrong in solver yaml" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::cout << "convection method : " << convMethod<< '\n';

        this->meshFormat   = meshFormat;
        this->meshFileName = meshFileName;
        this->solver = solver;

        this->endTimeControl = endTimeControl;
        this->nStep = nStep;
        this->outStepInterval = outStepInterval;

        this->dtControl = dtControl;
        this->dt = dt;
        this->cfl = cfl;
        this->cfl_pseudo = cfl_pseudo;
        this->dt_max = dt_max;
        this->dt_min = dt_min;
        this->timeIntegration = timeIntegration;
        this->nLoop = nLoop;
        this->cfl_pseudo = cfl_pseudo;

        initTimeIntegrationScheme(timeIntegration);

        this->convMethod = convMethod;

        this->isCompressible = isCompressible ;
        this->ro = ro;
        this->visc = visc;
        this->thermCond = thermCond;
        this->cp = cp;
        this->gamma = gamma;

        this->initial = initial;

    } catch (const YAML::BadFile& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);

    } catch (const YAML::ParserException& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);
    }
}

void solverConfig::initTimeIntegrationScheme(int timeIntegration){

    int ns;
    switch (timeIntegration) {
        case 1: // 1st Euler Explicit
            ns = 1;
            this->isImplicit = 0;
            this->nLoop= ns;
            this->coef_N.resize(ns);
            this->coef_M.resize(ns);
            this->coef_Res.resize(ns);

            coef_N[0]   = 1.0;
            coef_M[0]   = 0.0;
            coef_Res[0] = 1.0;
            break;
        case 3: // 3rd Runge-Kutta Explicit
            ns = 3;
            this->isImplicit = 0;
            this->nLoop = ns;
            this->coef_N.resize(ns);
            this->coef_M.resize(ns);
            this->coef_Res.resize(ns);

            coef_N[0]   = 1.0 ; coef_N[1]   = 0.75 ; coef_N[2]   = 1.0/3.0;
            coef_M[0]   = 0.0 ; coef_M[1]   = 0.25 ; coef_M[2]   = 2.0/3.0;
            coef_Res[0] = 1.0 ; coef_Res[1] = 0.25 ; coef_Res[2] = 2.0/3.0;
            break;

        case 10: // dual time stepping by explicit scheme
            ns = this->nLoop;
            this->isImplicit = 1;
            nLoop = this->nLoop;

            this->coef_N.resize(nLoop);
            this->coef_M.resize(nLoop);
            this->coef_Res.resize(nLoop);

            for (int i=0; i<nLoop; i++) 
            {
                coef_N[i] = 1.0;
                coef_M[i] = 0.0;
                coef_Res[i] = 1.0/(nLoop-i);
            }

            break;

    }

}
