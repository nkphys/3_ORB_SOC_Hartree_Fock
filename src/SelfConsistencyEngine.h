#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "tensor_type.h"
#include <iomanip>

#ifndef SelfConsistencyEngine_H
#define SelfConsistencyEngine_H


class SelfConsistencyEngine{
public:
    SelfConsistencyEngine(Parameters& Parameters__, Coordinates& Coordinates__,
                          MFParams& MFParams__, Hamiltonian& Hamiltonian__,
                          Observables& Observables__)
        : Parameters_(Parameters__),Coordinates_(Coordinates__),
          MFParams_(MFParams__), Hamiltonian_(Hamiltonian__),
          Observables_(Observables__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns)
    {

    }

    void RUN_SelfConsistencyEngine();
    double Prob (double muu, double mu_new);
    double ProbCluster (double muu, double mu_new);
    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    Hamiltonian &Hamiltonian_;
    Observables &Observables_;
    const int lx_, ly_, ns_;

};

/*
 * ***********
 *  Functions in Class SelfConsistencyEngine ------
 *  ***********
*/

void SelfConsistencyEngine::RUN_SelfConsistencyEngine(){

    complex<double> zero(0.0,0.0);

    double muu_prev;
    double Curr_QuantE;
    double Prev_QuantE;
    double Curr_ClassicalE;
    double Prev_ClassicalE;

    cout << "Temperature = " << Parameters_.Temperature<<" is being done"<<endl;


    string File_Out_progress;

    double initial_mu_guess;
    int n_states_occupied_zeroT;

    //starting with a random guess
    File_Out_progress = "output.txt";
    ofstream file_out_progress(File_Out_progress.c_str());


    file_out_progress<< "Maximum no of self consistency iterations = "<<Parameters_.IterMax<<"."<<endl;
    file_out_progress<<"Convergence error targetted = "<<Parameters_.Convergence_Error<<endl;

    Parameters_.Dflag='V'; // flag to calculate only Eigenvalue


    file_out_progress<<"Iter"<<setw(15)<<
                       "Error_OP"<<setw(15)<<
                       "mu"<<setw(17)<<
                       "E_CL"<<setw(17)<<"E_Quantum"<<endl;



    /*
    Prev_ClassicalE = Hamiltonian_.GetCLEnergy();
    Hamiltonian_.InteractionsCreate();
    Hamiltonian_.Diagonalize(Parameters_.Dflag);
    n_states_occupied_zeroT=Parameters_.ns*Parameters_.Fill*2.0;
    initial_mu_guess=0.5*(Hamiltonian_.eigs_[n_states_occupied_zeroT-1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);
    Parameters_.mus=Hamiltonian_.chemicalpotential(initial_mu_guess,Parameters_.Fill);
    Prev_QuantE = Hamiltonian_.E_QM();
    Hamiltonian_.copy_eigs(1);
    cout<<"Initial Classical Energy = "<<Prev_ClassicalE<<endl;
    cout<<"Initial Quantum Energy = "<<Prev_QuantE<<endl;
    cout<<"Initial mu="<<Parameters_.mus<<endl;
    */


    double Error_OP=100;
    Parameters_.mu_old=Parameters_.mus;

    for(int count=0;count<Parameters_.IterMax;count++){
        if(
                Error_OP>Parameters_.Convergence_Error
                ){



            Hamiltonian_.InteractionsCreate();
           // Hamiltonian_.Ham_.print();
            Hamiltonian_.Check_Hermiticity();
            Hamiltonian_.Diagonalize(Parameters_.Dflag);
            if(count==0){
                n_states_occupied_zeroT=Parameters_.Total_Particles;
                initial_mu_guess=0.5*(Hamiltonian_.eigs_[n_states_occupied_zeroT-1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);
                Parameters_.mu_old=Hamiltonian_.chemicalpotential(initial_mu_guess,Parameters_.Total_Particles);
            }


            //Getting order params + mu
            Observables_.Calculate_Order_Params();
            n_states_occupied_zeroT=Parameters_.Total_Particles;
            initial_mu_guess=0.5*(Hamiltonian_.eigs_[n_states_occupied_zeroT-1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);

            Parameters_.mus=Hamiltonian_.chemicalpotential(initial_mu_guess,Parameters_.Total_Particles);


            Curr_QuantE = Hamiltonian_.E_QM();
            Curr_ClassicalE = Hamiltonian_.GetCLEnergy();


            Observables_.Get_OrderParameters_diffs();

            Error_OP=Observables_.Error_OP_;



            file_out_progress<<setprecision(15)<<count<<setw(30)<<
                               Error_OP<<setw(30)<<
                               Parameters_.mus<<setw(32)<<
                               Curr_ClassicalE<<setw(32)<<Curr_QuantE<<endl;

            Observables_.Update_OrderParameters(count);


        }

    }

    string File_Out_Local_OP = Parameters_.File_OPs_out;
    ofstream file_out_Local_OP(File_Out_Local_OP.c_str());
    file_out_Local_OP<<"#site     state_i   state_j   OP.real    OP.imag"<<endl;

    for(int site=0;site<ns_;site++){
        for(int state_i=0;state_i<6;state_i++){
            for(int state_j=0;state_j<6;state_j++){
            file_out_Local_OP<<site<<"       "<<state_i<<"          "<<state_j<<
                               "        "<<Observables_.OParams_[site][state_i][state_j].real()<<
                               "        "<<Observables_.OParams_[site][state_i][state_j].imag()<<endl;
        }
    }
    }


} // ---------



#endif // SelfConsistencyEngine_H
