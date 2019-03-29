#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class
#define Hamiltonian_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);


class Hamiltonian {
public:

    Hamiltonian(Parameters& Parameters__, Coordinates&  Coordinates__, MFParams& MFParams__ )
        :Parameters_(Parameters__),Coordinates_(Coordinates__),MFParams_(MFParams__)

    {
        Initialize();
        Hoppings();
        HTBCreate();
    }


    void Initialize();    //::DONE
    void Hoppings();        //::DONE
    double GetCLEnergy();    //::DONE
    void InteractionsCreate();   //::DONE
    void Check_Hermiticity();  //::DONE
    void HTBCreate();   //::DONE
    double chemicalpotential(double muin,double Particles);    //::DONE

    double TotalDensity();   //::DONE
    double E_QM();   //::DONE

    void Diagonalize(char option);   //::DONE
    void copy_eigs(int i);  //::DONE

    int convert_jm_to_int(string jm_val);

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    int lx_, ly_, ns_, orbs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx,Ty,Tpxpy,Tpxmy;
    vector<double> eigs_,eigs_saved_,sx_,sy_,sz_;

    double HS_factor;

};



double Hamiltonian::chemicalpotential(double muin,double Particles){


    double mu_out;
    double n1,N;
    double dMubydN;
    double nstate = eigs_.size();
    dMubydN = 0.0000005*(eigs_[nstate-1] - eigs_[0])/nstate;
    N=Particles;
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged=false;


    if(1==1){
        for(int i=0;i<50000;i++){
            n1=0.0;
            for(int j=0;j<nstate;j++){
                n1+=double(1.0/( exp( (eigs_[j]-mu_out)*Parameters_.beta ) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if(abs(N-n1)<double(0.0001)){
                //cout<<abs(N-n1)<<endl;
                converged=true;
                break;
            }
            else {
                mu_out += (N-n1)*dMubydN;
                //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;

            }
        }

        if(!converged){
            cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else{
            cout<<"mu converged, N = "<<n1<<endl;
        }

    }


    double mu1, mu2;
    double mu_temp = muin;
    //cout<<"mu_input = "<<mu_temp<<endl;
    if(1==2){
        mu1=eigs_[0];
        mu2=eigs_[nstate-1];
        for(int i=0;i<40000;i++){
            n1=0.0;
            for(int j=0;j<nstate;j++){
                n1+=double(1.0/( exp( (eigs_[j]-mu_temp)*Parameters_.beta ) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if(abs(N-n1)<double(0.0001)){
                //cout<<abs(N-n1)<<endl;
                converged=true;
                break;
            }
            else {
                if(n1 >N){
                    mu2=mu_temp;
                    mu_temp=0.5*(mu1 + mu_temp);
                }
                else{
                    mu1=mu_temp;
                    mu_temp=0.5*(mu2 + mu_temp);
                }

            }
            //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
        }

        if(!converged){
            cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else{
            cout<<"mu converged, N = "<<n1<<endl;
        }

        mu_out = mu_temp;
    }

    return mu_out;
} // ----------


void Hamiltonian::Initialize(){


    ly_=Parameters_.ly;
    lx_=Parameters_.lx;
    ns_=Parameters_.ns;


    int space=Coordinates_.no_dof_ ;

    HTB_.resize(space,space);
    Ham_.resize(space,space);
    eigs_.resize(space);
    eigs_saved_.resize(space);

} // ----------

double Hamiltonian::TotalDensity(){

    double n1=0.0;
    /*
    for(int j=0;j<eigs_.size();j++){
        n1 +=  1.0f/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    */
    return n1;

} // ----------



double Hamiltonian::E_QM(){

    double E=0.0;
    for(int j=0;j<eigs_.size();j++){
        //E +=  (eigs_[j]-Parameters_.mus)/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
        E +=  (eigs_[j])/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    return E;

} // ----------



double Hamiltonian::GetCLEnergy(){

    double EClassical=0.0;
    complex<double> Classical_Energy;
    Classical_Energy=zero_complex;

    int n_i, n_j;
    int op_i, op_j;
    int op_state_i, op_state_j;
    int state_i, state_j;

    Mat_2_string STATE_i_, STATE_j_, OP_STATE_i_, OP_STATE_j_;
    Mat_2_doub Signs_;
    Mat_1_Complex_doub VALUES_;

    Mat_1_string TEMP_STATE_i_, TEMP_STATE_j_, TEMP_OP_STATE_i_, TEMP_OP_STATE_j_;
    Mat_1_doub TEMP_Signs_;
    complex<double> temp_value_;

    for(int site=0;site<ns_;site++) {



        //-----H_{1/2}--------//
        for(state_i=4;state_i<6;state_i++){
            n_i=Coordinates_.Nc_dof(site,state_i);
            op_state_i=4 + abs((state_i - 4) - 1);
            op_i=Coordinates_.Nc_dof(site,op_state_i);

            for(state_j=4;state_j<6;state_j++){
                n_j=Coordinates_.Nc_dof(site,state_j);
                op_state_j=4 + abs((state_j - 4) - 1);
                op_j=Coordinates_.Nc_dof(site,op_state_j);

                //Value in front of (a^{dagger}_{ni}a_{nj})=OP[op_i][op_j]
                temp_value_ = one_complex*((pow(-1.0,(state_i-state_j)*1.0))*(Parameters_.U_onsite - ((4*Parameters_.J_Hund)/3.0)));
                Classical_Energy += temp_value_*MFParams_.OParams[site][op_state_i][op_state_j]*MFParams_.OParams[site][state_i][state_j];

            }
        }

        //--------------------//



        //-----H_{3/2}--------//
        STATE_i_.clear();STATE_j_.clear();
        OP_STATE_i_.clear();OP_STATE_j_.clear();
        Signs_.clear();
        VALUES_.clear();


        //Term for (U-JH)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2"};
        TEMP_STATE_j_ =    {"3by2_m1by2","3by2_1by2","3by2_1by2","3by2_m1by2","3by2_m3by2","3by2_3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_j_ = {"3by2_1by2","3by2_m1by2","3by2_m1by2","3by2_1by2","3by2_3by2","3by2_m3by2","3by2_m3by2","3by2_3by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - (Parameters_.J_Hund));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Term for U- (7JH/3) [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_3by2","3by2_1by2","3by2_3by2","3by2_m1by2","3by2_3by2","3by2_m1by2","3by2_3by2"};
        TEMP_STATE_j_ =    {"3by2_1by2","3by2_3by2","3by2_3by2","3by2_1by2","3by2_m1by2","3by2_3by2","3by2_3by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_3by2","3by2_1by2","3by2_3by2","3by2_1by2","3by2_3by2","3by2_m1by2","3by2_3by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_3by2","3by2_1by2","3by2_1by2","3by2_3by2","3by2_3by2","3by2_m1by2","3by2_m1by2","3by2_3by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - (7.0*Parameters_.J_Hund/3.0));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Term for U- (7JH/3) [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_m3by2","3by2_1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2"};
        TEMP_STATE_j_ =    {"3by2_1by2","3by2_m3by2","3by2_m3by2","3by2_1by2","3by2_m1by2","3by2_m3by2","3by2_m3by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_1by2","3by2_m3by2","3by2_1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m3by2","3by2_1by2","3by2_1by2","3by2_m3by2","3by2_m3by2","3by2_m1by2","3by2_m1by2","3by2_m3by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - (7.0*Parameters_.J_Hund/3.0));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);




        //Terms (-4JH/3)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_1by2","3by2_m1by2","3by2_m1by2","3by2_1by2"};
        TEMP_STATE_j_ =    {"3by2_1by2","3by2_m1by2","3by2_m1by2","3by2_1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m1by2","3by2_1by2","3by2_1by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m1by2","3by2_1by2","3by2_1by2","3by2_m1by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2"};
        TEMP_Signs_ = {-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real(-4.0*Parameters_.J_Hund/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //--------------------//



        //-----H_{mix}--------//

        //Terms (U-2JH)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - (2.0*Parameters_.J_Hund));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (U-8JH/3)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","3by2_m3by2","3by2_3by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_3by2","1by2_m1by2","1by2_1by2","3by2_m3by2","3by2_3by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m3by2","3by2_3by2","1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2"};
        TEMP_Signs_ = {1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - ((8.0*Parameters_.J_Hund)/3.0));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (U-7JH/3)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","3by2_m1by2","3by2_1by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m1by2","3by2_1by2","1by2_m1by2","1by2_1by2","3by2_m1by2","3by2_1by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m1by2","3by2_1by2","1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - ((7.0*Parameters_.J_Hund)/3.0));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (U-5JH/3)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2"};
        TEMP_Signs_ = {1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - ((5.0*Parameters_.J_Hund)/3.0));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (2JH/3)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_m1by2","3by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2"};
        TEMP_STATE_j_ =    {"3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2","1by2_m1by2","3by2_1by2","1by2_m1by2","3by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","1by2_m1by2","3by2_1by2","1by2_m1by2","1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2","3by2_1by2","1by2_m1by2","3by2_1by2","1by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((2.0*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (JH/sqrt(3)) [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","3by2_1by2","1by2_1by2","3by2_1by2","3by2_3by2","1by2_m1by2","1by2_m1by2","3by2_3by2"};
        TEMP_STATE_j_ =    {"3by2_3by2","1by2_m1by2","1by2_m1by2","3by2_3by2","1by2_1by2","3by2_1by2","1by2_1by2","3by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","1by2_1by2","3by2_1by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_3by2","1by2_m1by2"};
        TEMP_OP_STATE_j_ = {"1by2_m1by2","3by2_3by2","3by2_3by2","1by2_m1by2","3by2_1by2","1by2_1by2","3by2_1by2","1by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real(0.0);temp_value_.imag((Parameters_.J_Hund)/sqrt(3.0));
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (JH/sqrt(3)) [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_m1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2","3by2_m3by2","1by2_1by2","1by2_1by2","3by2_m3by2"};
        TEMP_STATE_j_ =    {"3by2_m3by2","1by2_1by2","1by2_1by2","3by2_m3by2","1by2_m1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m1by2","1by2_m1by2","3by2_m1by2","1by2_m1by2","1by2_1by2","3by2_m3by2","3by2_m3by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"1by2_1by2","3by2_m3by2","3by2_m3by2","1by2_1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2","1by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real(0.0);temp_value_.imag((Parameters_.J_Hund)/sqrt(3.0));
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (-5JH/3)[m=1by2]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","1by2_1by2","1by2_m1by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"1by2_1by2","1by2_m1by2","1by2_m1by2","1by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((-5.0*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (-5JH/3)[m=3by2]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2","1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2","1by2_1by2","1by2_m1by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"1by2_1by2","1by2_m1by2","1by2_m1by2","1by2_1by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((-5.0*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (JH*sqrt(2)/3) [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_m1by2","3by2_3by2","3by2_m1by2","3by2_3by2","1by2_m1by2","3by2_3by2","3by2_3by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","3by2_3by2","3by2_3by2","1by2_m1by2","3by2_m1by2","3by2_3by2","3by2_m1by2","3by2_3by2"};
        TEMP_OP_STATE_i_ = {"3by2_3by2","3by2_m1by2","3by2_3by2","3by2_m1by2","3by2_3by2","1by2_m1by2","1by2_m1by2","3by2_3by2"};
        TEMP_OP_STATE_j_ = {"3by2_3by2","1by2_m1by2","1by2_m1by2","3by2_3by2","3by2_3by2","3by2_m1by2","3by2_3by2","3by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Terms (JH*sqrt(2)/3) [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_m3by2","3by2_1by2","3by2_m3by2","1by2_1by2","3by2_m3by2","3by2_m3by2","1by2_1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","3by2_m3by2","3by2_m3by2","1by2_1by2","3by2_1by2","3by2_m3by2","3by2_1by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_1by2","3by2_m3by2","3by2_1by2","3by2_m3by2","1by2_1by2","1by2_1by2","3by2_m3by2"};
        TEMP_OP_STATE_j_ = {"3by2_m3by2","1by2_1by2","1by2_1by2","3by2_m3by2","3by2_m3by2","3by2_1by2","3by2_m3by2","3by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (-4*JH/3*sqrt(2)) [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_3by2","3by2_1by2","3by2_3by2","1by2_1by2","3by2_3by2","3by2_3by2","1by2_1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","3by2_3by2","3by2_3by2","1by2_1by2","3by2_1by2","3by2_3by2","3by2_1by2","3by2_3by2"};
        TEMP_OP_STATE_i_ = {"3by2_3by2","3by2_1by2","3by2_3by2","3by2_1by2","3by2_3by2","1by2_1by2","1by2_1by2","3by2_3by2"};
        TEMP_OP_STATE_j_ = {"3by2_3by2","1by2_1by2","1by2_1by2","3by2_3by2","3by2_3by2","3by2_1by2","3by2_3by2","3by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((-4.0*Parameters_.J_Hund)/(3.0*sqrt(2.0)));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (-4*JH/3*sqrt(2)) [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_m1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2","1by2_m1by2","3by2_m3by2","3by2_m3by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","3by2_m3by2","3by2_m3by2","1by2_m1by2","3by2_m1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_m1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2","1by2_m1by2","1by2_m1by2","3by2_m3by2"};
        TEMP_OP_STATE_j_ = {"3by2_m3by2","1by2_m1by2","1by2_m1by2","3by2_m3by2","3by2_m3by2","3by2_m1by2","3by2_m3by2","3by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((-4.0*Parameters_.J_Hund)/(3.0*sqrt(2.0)));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Terms (JH*sqrt(2)/3), n*a^{dagger}*a part of 2nd last term, [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","1by2_m1by2","3by2_1by2","3by2_1by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","3by2_1by2","3by2_1by2","1by2_m1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","1by2_m1by2","1by2_m1by2","3by2_1by2"};
        TEMP_OP_STATE_j_ = {"3by2_1by2","1by2_m1by2","1by2_m1by2","3by2_1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (JH*sqrt(2)/3), n*a^{dagger}*a part of 2nd last term, [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Terms (JH*sqrt(2)/3), a^{dagger}*a*a^{dagger}*a part of 2nd last term, [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2","1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2"};
        TEMP_Signs_ = {-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (JH*sqrt(2)/3), a^{dagger}*a*a^{dagger}*a part of 2nd last term, [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2","1by2_m1by2","3by2_1by2","3by2_1by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","3by2_1by2","3by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_1by2","1by2_m1by2","1by2_m1by2","3by2_1by2"};
        TEMP_OP_STATE_j_ = {"3by2_1by2","1by2_m1by2","1by2_m1by2","3by2_1by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Terms (JH*sqrt(2/3)) [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","3by2_1by2","1by2_1by2","3by2_1by2","3by2_3by2","3by2_m1by2","3by2_m1by2","3by2_3by2"};
        TEMP_STATE_j_ =    {"3by2_3by2","3by2_m1by2","3by2_m1by2","3by2_3by2","1by2_1by2","3by2_1by2","1by2_1by2","3by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","1by2_1by2","3by2_1by2","1by2_1by2","3by2_m1by2","3by2_3by2","3by2_3by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m1by2","3by2_3by2","3by2_3by2","3by2_m1by2","3by2_1by2","1by2_1by2","3by2_1by2","1by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real(0.0);temp_value_.imag((Parameters_.J_Hund)*sqrt(2.0/3.0));
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (JH*sqrt(2/3)) [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_m1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2","3by2_m3by2","3by2_1by2","3by2_1by2","3by2_m3by2"};
        TEMP_STATE_j_ =    {"3by2_m3by2","3by2_1by2","3by2_1by2","3by2_m3by2","1by2_m1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m1by2","1by2_m1by2","3by2_m1by2","1by2_m1by2","3by2_1by2","3by2_m3by2","3by2_m3by2","3by2_1by2"};
        TEMP_OP_STATE_j_ = {"3by2_1by2","3by2_m3by2","3by2_m3by2","3by2_1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2","1by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real(0.0);temp_value_.imag((Parameters_.J_Hund)*sqrt(2.0/3.0));
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        for(int term_no=0;term_no<VALUES_.size();term_no++){
            for(int term=0;term<STATE_i_[term_no].size();term++){
                state_i=convert_jm_to_int(STATE_i_[term_no][term]);
                state_j=convert_jm_to_int(STATE_j_[term_no][term]);
                op_state_i=convert_jm_to_int(OP_STATE_i_[term_no][term]);
                op_state_j=convert_jm_to_int(OP_STATE_j_[term_no][term]);
                n_i=Coordinates_.Nc_dof(site,state_i);
                n_j=Coordinates_.Nc_dof(site,state_j);
                op_i=Coordinates_.Nc_dof(site,op_state_i);
                op_j=Coordinates_.Nc_dof(site,op_state_j);
                Classical_Energy += Signs_[term_no][term]*VALUES_[term_no]*
                        (MFParams_.OParams[site][op_state_i][op_state_j])*MFParams_.OParams[site][state_i][state_j];
            }
        }


        //--------------------//

    }

    EClassical = Classical_Energy.real()*(-1.0)*(1/2.0);


    return EClassical;


} // ----------



void Hamiltonian::InteractionsCreate(){

    int space=Coordinates_.no_dof_;
    int n_i, n_j;
    int op_i, op_j;
    int op_state_i, op_state_j;
    int state_i, state_j;

    Mat_2_string STATE_i_, STATE_j_, OP_STATE_i_, OP_STATE_j_;
    Mat_2_doub Signs_;
    Mat_1_Complex_doub VALUES_;

    Mat_1_string TEMP_STATE_i_, TEMP_STATE_j_, TEMP_OP_STATE_i_, TEMP_OP_STATE_j_;
    Mat_1_doub TEMP_Signs_;
    complex<double> temp_value_;


    for(int i=0;i<space;i++) {
        for(int j=0;j<space;j++) {
            Ham_(i,j)=HTB_(i,j);
        }
    }

    for(int site=0;site<ns_;site++) {



        //-----H_{1/2}--------//
        for(state_i=4;state_i<6;state_i++){
            n_i=Coordinates_.Nc_dof(site,state_i);
            op_state_i=4 + abs((state_i - 4) - 1);
            op_i=Coordinates_.Nc_dof(site,op_state_i);

            for(state_j=4;state_j<6;state_j++){
                n_j=Coordinates_.Nc_dof(site,state_j);
                op_state_j=4 + abs((state_j - 4) - 1);
                op_j=Coordinates_.Nc_dof(site,op_state_j);

                //Value in front of (a^{dagger}_{ni}a_{nj})=OP[op_i][op_j]
                temp_value_ = one_complex*((pow(-1.0,(state_i-state_j)*1.0))*(Parameters_.U_onsite - ((4*Parameters_.J_Hund)/3.0)));
                Ham_(n_i,n_j) += temp_value_*MFParams_.OParams[site][op_state_i][op_state_j];

            }
        }

        //--------------------//



        //-----H_{3/2}--------//
        STATE_i_.clear();STATE_j_.clear();
        OP_STATE_i_.clear();OP_STATE_j_.clear();
        Signs_.clear();
        VALUES_.clear();


        //Term for (U-JH)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2"};
        TEMP_STATE_j_ =    {"3by2_m1by2","3by2_1by2","3by2_1by2","3by2_m1by2","3by2_m3by2","3by2_3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_j_ = {"3by2_1by2","3by2_m1by2","3by2_m1by2","3by2_1by2","3by2_3by2","3by2_m3by2","3by2_m3by2","3by2_3by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - (Parameters_.J_Hund));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Term for U- (7JH/3) [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_3by2","3by2_1by2","3by2_3by2","3by2_m1by2","3by2_3by2","3by2_m1by2","3by2_3by2"};
        TEMP_STATE_j_ =    {"3by2_1by2","3by2_3by2","3by2_3by2","3by2_1by2","3by2_m1by2","3by2_3by2","3by2_3by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_3by2","3by2_1by2","3by2_3by2","3by2_1by2","3by2_3by2","3by2_m1by2","3by2_3by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_3by2","3by2_1by2","3by2_1by2","3by2_3by2","3by2_3by2","3by2_m1by2","3by2_m1by2","3by2_3by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - (7.0*Parameters_.J_Hund/3.0));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Term for U- (7JH/3) [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_m3by2","3by2_1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2"};
        TEMP_STATE_j_ =    {"3by2_1by2","3by2_m3by2","3by2_m3by2","3by2_1by2","3by2_m1by2","3by2_m3by2","3by2_m3by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_1by2","3by2_m3by2","3by2_1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m3by2","3by2_1by2","3by2_1by2","3by2_m3by2","3by2_m3by2","3by2_m1by2","3by2_m1by2","3by2_m3by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - (7.0*Parameters_.J_Hund/3.0));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);




        //Terms (-4JH/3)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_1by2","3by2_m1by2","3by2_m1by2","3by2_1by2"};
        TEMP_STATE_j_ =    {"3by2_1by2","3by2_m1by2","3by2_m1by2","3by2_1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m1by2","3by2_1by2","3by2_1by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m1by2","3by2_1by2","3by2_1by2","3by2_m1by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2"};
        TEMP_Signs_ = {-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real(-4.0*Parameters_.J_Hund/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //--------------------//



        //-----H_{mix}--------//

        //Terms (U-2JH)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - (2.0*Parameters_.J_Hund));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (U-8JH/3)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","3by2_m3by2","3by2_3by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_3by2","1by2_m1by2","1by2_1by2","3by2_m3by2","3by2_3by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m3by2","3by2_3by2","1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2"};
        TEMP_Signs_ = {1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - ((8.0*Parameters_.J_Hund)/3.0));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (U-7JH/3)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","3by2_m1by2","3by2_1by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m1by2","3by2_1by2","1by2_m1by2","1by2_1by2","3by2_m1by2","3by2_1by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m1by2","3by2_1by2","1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - ((7.0*Parameters_.J_Hund)/3.0));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (U-5JH/3)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_3by2","3by2_m3by2","1by2_1by2","1by2_m1by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2"};
        TEMP_Signs_ = {1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0};
        temp_value_.real(Parameters_.U_onsite - ((5.0*Parameters_.J_Hund)/3.0));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (2JH/3)
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_m1by2","3by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2"};
        TEMP_STATE_j_ =    {"3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2","1by2_m1by2","3by2_1by2","1by2_m1by2","3by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","1by2_m1by2","3by2_1by2","1by2_m1by2","1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2","3by2_1by2","1by2_m1by2","3by2_1by2","1by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((2.0*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (JH/sqrt(3)) [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","3by2_1by2","1by2_1by2","3by2_1by2","3by2_3by2","1by2_m1by2","1by2_m1by2","3by2_3by2"};
        TEMP_STATE_j_ =    {"3by2_3by2","1by2_m1by2","1by2_m1by2","3by2_3by2","1by2_1by2","3by2_1by2","1by2_1by2","3by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","1by2_1by2","3by2_1by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_3by2","1by2_m1by2"};
        TEMP_OP_STATE_j_ = {"1by2_m1by2","3by2_3by2","3by2_3by2","1by2_m1by2","3by2_1by2","1by2_1by2","3by2_1by2","1by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real(0.0);temp_value_.imag((Parameters_.J_Hund)/sqrt(3.0));
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (JH/sqrt(3)) [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_m1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2","3by2_m3by2","1by2_1by2","1by2_1by2","3by2_m3by2"};
        TEMP_STATE_j_ =    {"3by2_m3by2","1by2_1by2","1by2_1by2","3by2_m3by2","1by2_m1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m1by2","1by2_m1by2","3by2_m1by2","1by2_m1by2","1by2_1by2","3by2_m3by2","3by2_m3by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"1by2_1by2","3by2_m3by2","3by2_m3by2","1by2_1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2","1by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real(0.0);temp_value_.imag((Parameters_.J_Hund)/sqrt(3.0));
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (-5JH/3)[m=1by2]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","1by2_1by2","1by2_m1by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"1by2_1by2","1by2_m1by2","1by2_m1by2","1by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((-5.0*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (-5JH/3)[m=3by2]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2","1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","1by2_1by2","1by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2","1by2_1by2","1by2_m1by2","1by2_m1by2","1by2_1by2"};
        TEMP_OP_STATE_j_ = {"1by2_1by2","1by2_m1by2","1by2_m1by2","1by2_1by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((-5.0*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (JH*sqrt(2)/3) [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_m1by2","3by2_3by2","3by2_m1by2","3by2_3by2","1by2_m1by2","3by2_3by2","3by2_3by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","3by2_3by2","3by2_3by2","1by2_m1by2","3by2_m1by2","3by2_3by2","3by2_m1by2","3by2_3by2"};
        TEMP_OP_STATE_i_ = {"3by2_3by2","3by2_m1by2","3by2_3by2","3by2_m1by2","3by2_3by2","1by2_m1by2","1by2_m1by2","3by2_3by2"};
        TEMP_OP_STATE_j_ = {"3by2_3by2","1by2_m1by2","1by2_m1by2","3by2_3by2","3by2_3by2","3by2_m1by2","3by2_3by2","3by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Terms (JH*sqrt(2)/3) [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_m3by2","3by2_1by2","3by2_m3by2","1by2_1by2","3by2_m3by2","3by2_m3by2","1by2_1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","3by2_m3by2","3by2_m3by2","1by2_1by2","3by2_1by2","3by2_m3by2","3by2_1by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_1by2","3by2_m3by2","3by2_1by2","3by2_m3by2","1by2_1by2","1by2_1by2","3by2_m3by2"};
        TEMP_OP_STATE_j_ = {"3by2_m3by2","1by2_1by2","1by2_1by2","3by2_m3by2","3by2_m3by2","3by2_1by2","3by2_m3by2","3by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);


        //Terms (-4*JH/3*sqrt(2)) [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_3by2","3by2_1by2","3by2_3by2","1by2_1by2","3by2_3by2","3by2_3by2","1by2_1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","3by2_3by2","3by2_3by2","1by2_1by2","3by2_1by2","3by2_3by2","3by2_1by2","3by2_3by2"};
        TEMP_OP_STATE_i_ = {"3by2_3by2","3by2_1by2","3by2_3by2","3by2_1by2","3by2_3by2","1by2_1by2","1by2_1by2","3by2_3by2"};
        TEMP_OP_STATE_j_ = {"3by2_3by2","1by2_1by2","1by2_1by2","3by2_3by2","3by2_3by2","3by2_1by2","3by2_3by2","3by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((-4.0*Parameters_.J_Hund)/(3.0*sqrt(2.0)));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (-4*JH/3*sqrt(2)) [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_m1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2","1by2_m1by2","3by2_m3by2","3by2_m3by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","3by2_m3by2","3by2_m3by2","1by2_m1by2","3by2_m1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_m1by2","3by2_m3by2","3by2_m1by2","3by2_m3by2","1by2_m1by2","1by2_m1by2","3by2_m3by2"};
        TEMP_OP_STATE_j_ = {"3by2_m3by2","1by2_m1by2","1by2_m1by2","3by2_m3by2","3by2_m3by2","3by2_m1by2","3by2_m3by2","3by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((-4.0*Parameters_.J_Hund)/(3.0*sqrt(2.0)));temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Terms (JH*sqrt(2)/3), n*a^{dagger}*a part of 2nd last term, [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","1by2_m1by2","3by2_1by2","3by2_1by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","3by2_1by2","3by2_1by2","1by2_m1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","1by2_m1by2","1by2_m1by2","3by2_1by2"};
        TEMP_OP_STATE_j_ = {"3by2_1by2","1by2_m1by2","1by2_m1by2","3by2_1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (JH*sqrt(2)/3), n*a^{dagger}*a part of 2nd last term, [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2","3by2_m1by2","3by2_1by2","3by2_m1by2","3by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Terms (JH*sqrt(2)/3), a^{dagger}*a*a^{dagger}*a part of 2nd last term, [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2","1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2"};
        TEMP_STATE_j_ =    {"1by2_1by2","3by2_m1by2","3by2_m1by2","1by2_1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m1by2","1by2_1by2","1by2_1by2","3by2_m1by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2"};
        TEMP_Signs_ = {-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (JH*sqrt(2)/3), a^{dagger}*a*a^{dagger}*a part of 2nd last term, [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2","1by2_m1by2","3by2_1by2","3by2_1by2","1by2_m1by2"};
        TEMP_STATE_j_ =    {"1by2_m1by2","3by2_1by2","3by2_1by2","1by2_m1by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_m3by2"};
        TEMP_OP_STATE_i_ = {"3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2","3by2_1by2","1by2_m1by2","1by2_m1by2","3by2_1by2"};
        TEMP_OP_STATE_j_ = {"3by2_1by2","1by2_m1by2","1by2_m1by2","3by2_1by2","3by2_m3by2","3by2_3by2","3by2_m3by2","3by2_3by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0};
        temp_value_.real((sqrt(2.0)*Parameters_.J_Hund)/3.0);temp_value_.imag(0.0);
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        //Terms (JH*sqrt(2/3)) [s=+1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_1by2","3by2_1by2","1by2_1by2","3by2_1by2","3by2_3by2","3by2_m1by2","3by2_m1by2","3by2_3by2"};
        TEMP_STATE_j_ =    {"3by2_3by2","3by2_m1by2","3by2_m1by2","3by2_3by2","1by2_1by2","3by2_1by2","1by2_1by2","3by2_1by2"};
        TEMP_OP_STATE_i_ = {"3by2_1by2","1by2_1by2","3by2_1by2","1by2_1by2","3by2_m1by2","3by2_3by2","3by2_3by2","3by2_m1by2"};
        TEMP_OP_STATE_j_ = {"3by2_m1by2","3by2_3by2","3by2_3by2","3by2_m1by2","3by2_1by2","1by2_1by2","3by2_1by2","1by2_1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real(0.0);temp_value_.imag((Parameters_.J_Hund)*sqrt(2.0/3.0));
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);

        //Terms (JH*sqrt(2/3)) [s=-1]
        TEMP_STATE_i_.clear();TEMP_STATE_j_.clear();
        TEMP_OP_STATE_i_.clear();TEMP_OP_STATE_j_.clear();TEMP_Signs_.clear();
        TEMP_STATE_i_ =    {"1by2_m1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2","3by2_m3by2","3by2_1by2","3by2_1by2","3by2_m3by2"};
        TEMP_STATE_j_ =    {"3by2_m3by2","3by2_1by2","3by2_1by2","3by2_m3by2","1by2_m1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2"};
        TEMP_OP_STATE_i_ = {"3by2_m1by2","1by2_m1by2","3by2_m1by2","1by2_m1by2","3by2_1by2","3by2_m3by2","3by2_m3by2","3by2_1by2"};
        TEMP_OP_STATE_j_ = {"3by2_1by2","3by2_m3by2","3by2_m3by2","3by2_1by2","3by2_m1by2","1by2_m1by2","3by2_m1by2","1by2_m1by2"};
        TEMP_Signs_ = {1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0};
        temp_value_.real(0.0);temp_value_.imag((Parameters_.J_Hund)*sqrt(2.0/3.0));
        STATE_i_.push_back(TEMP_STATE_i_);STATE_j_.push_back(TEMP_STATE_j_);
        OP_STATE_i_.push_back(TEMP_OP_STATE_i_);OP_STATE_j_.push_back(TEMP_OP_STATE_j_);
        Signs_.push_back(TEMP_Signs_);
        VALUES_.push_back(temp_value_);



        for(int term_no=0;term_no<VALUES_.size();term_no++){
            for(int term=0;term<STATE_i_[term_no].size();term++){
                state_i=convert_jm_to_int(STATE_i_[term_no][term]);
                state_j=convert_jm_to_int(STATE_j_[term_no][term]);
                op_state_i=convert_jm_to_int(OP_STATE_i_[term_no][term]);
                op_state_j=convert_jm_to_int(OP_STATE_j_[term_no][term]);
                n_i=Coordinates_.Nc_dof(site,state_i);
                n_j=Coordinates_.Nc_dof(site,state_j);
                op_i=Coordinates_.Nc_dof(site,op_state_i);
                op_j=Coordinates_.Nc_dof(site,op_state_j);
                Ham_(n_i,n_j) += Signs_[term_no][term]*VALUES_[term_no]*
                        (MFParams_.OParams[site][op_state_i][op_state_j]);
            }
        }


        //--------------------//


    }


} // ----------


int Hamiltonian::convert_jm_to_int(string jm_val){

    int val;
    if(jm_val=="3by2_m3by2"){val=0;}
    if(jm_val=="3by2_3by2"){val=1;}
    if(jm_val=="3by2_m1by2"){val=2;}
    if(jm_val=="3by2_1by2"){val=3;}
    if(jm_val=="1by2_m1by2"){val=4;}
    if(jm_val=="1by2_1by2"){val=5;}
    return val;
}

void Hamiltonian::Check_Hermiticity()

{
    complex<double> temp(0,0);
    complex<double>temp2;

    for(int i=0;i<Ham_.n_row();i++) {
        for(int j=0;j<Ham_.n_row();j++) {
            if(
                    abs(Ham_(i,j) - conj(Ham_(j,i)))>0.00001
                    ) {
                cout<<Ham_(i,j)<<endl;
                cout<<conj(Ham_(j,i))<<endl;

            }
            assert(
                        abs(Ham_(i,j) - conj(Ham_(j,i)))<0.00001
                        ); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}





void Hamiltonian::Diagonalize(char option){

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);


    char jobz=option;
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


}


void Hamiltonian::HTBCreate(){

    int mx=Parameters_.TBC_mx;
    int my=Parameters_.TBC_my;
    complex<double> phasex, phasey;
    int l,m,a,b;
    complex<double> Boundary_val;
    string boundary_cond="PBC";

    if(boundary_cond=="OBC"){
        Boundary_val=zero_complex;
    }
    else if(boundary_cond=="PBC"){
        Boundary_val=one_complex;
    }

    HTB_.fill(0.0);

    for(l=0;l<ns_;l++) {

        // * +x direction Neighbor
        if(Coordinates_.lx_>1){

            if(Coordinates_.indx(l)==(Coordinates_.lx_ -1)){
                phasex=Boundary_val*exp(iota_complex*2.0*(1.0*mx)*PI/(1.0*Parameters_.TBC_cellsX));
                phasey=one_complex;
            }
            else{
                phasex=one_complex;
                phasey=one_complex;
            }
            m = Coordinates_.neigh(l,0);
            for (int state=0;state<6;state++){
                b=Coordinates_.Nc_dof(l,state);
                for(int state_neigh=0;state_neigh<6;state_neigh++) {
                    a=Coordinates_.Nc_dof(m,state_neigh);
                    //value*c^{\dagger}_{a}c_{b}
                    assert (a!=b);
                    if(a!=b){
                        HTB_(a,b)=Parameters_.Hopping_jm_NN(state_neigh,state)*phasex;
                        HTB_(b,a)=conj(HTB_(a,b));
                    }
                }
            }

        }


        // * +y direction Neighbor
        if(Coordinates_.ly_>1){

            if(Coordinates_.indy(l)==(Coordinates_.ly_ -1)){
                phasex=one_complex;
                phasey=Boundary_val*exp(iota_complex*2.0*(1.0*my)*PI/(1.0*Parameters_.TBC_cellsY));
            }
            else{
                phasex=one_complex;
                phasey=one_complex;
            }
            m = Coordinates_.neigh(l,2);
            for (int state=0;state<6;state++){
                b=Coordinates_.Nc_dof(l,state);
                for(int state_neigh=0;state_neigh<6;state_neigh++) {
                    a=Coordinates_.Nc_dof(m,state_neigh);
                    //value*c^{\dagger}_{a}c_{b}
                    assert (a!=b);
                    if(a!=b){
                        HTB_(a,b)=Parameters_.Hopping_jm_NN(state_neigh,state)*phasey;
                        HTB_(b,a)=conj(HTB_(a,b));
                    }
                }
            }

        }

    }



    //Spin-orbit coupling
    for(int i=0;i<ns_;i++) {
        for(int state=0;state<6;state++) {
            a=Coordinates_.Nc_dof(i,state);
            if(state<4){
                HTB_(a,a)+=complex<double>(-0.5,0.0)*Parameters_.Lambda_SOC;
            }
            else{
                HTB_(a,a)+=complex<double>(1.0,0.0)*Parameters_.Lambda_SOC;
            }
        }
    }


    //t2g_crystal_field
    for(int i=0;i<ns_;i++) {
        for(int state=0;state<6;state++) {
            b=Coordinates_.Nc_dof(i,state);
            for(int state_p=0;state_p<6;state_p++) {
                a=Coordinates_.Nc_dof(i,state_p);
                //if(a!=b){
                HTB_(a,b)+=(Parameters_.Crystal_Field_jm(state_p,state));
                //HTB_(b,a)=conj(HTB_(a,b));
                //}
            }
        }
    }


    //HTB_.print();
} // ----------



void Hamiltonian::Hoppings(){
    //DOES SOMETHING EXACT i.e NOTHING :)

} // ----------

void Hamiltonian::copy_eigs(int i){

    int space=2*ns_;

    if (i == 0) {
        for(int j=0;j<space;j++) {
            eigs_[j] = eigs_saved_[j];
        }
    }
    else {
        for(int j=0;j<space;j++) {
            eigs_saved_[j] = eigs_[j];
        }
    }

}


#endif
