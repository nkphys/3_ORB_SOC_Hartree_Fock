#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "tensor_type.h"

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#define PI acos(-1.0)

class Observables{
public:

    Observables(Parameters& Parameters__, Coordinates& Coordinates__,
                MFParams& MFParams__, Hamiltonian& Hamiltonian__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns)
    {
        Initialize();
    }

    void Initialize();

    void OccDensity();
    void Calculate_Akw();
    void Calculate_Nw();
    void Calculate_SpinSpincorrelations();
    void Calculate_SpinSpincorrelations_Smartly();
    void Calculate_IPR();
    void Calculate_Optical_Conductivity();
    void Get_Non_Interacting_dispersion();
    double Lorentzian(double x, double brd);
    void TotalOccDensity();
    void DensityOfStates();

    void Total_Energy_Average(double Curr_QuantE, double CurrE);

    void OccDensity(int tlabel);
    void Calculate_Local_Density();
    void Calculate_Order_Params();
    void Get_OrderParameters_diffs();
    void Update_OrderParameters(int iter);
    void Update_OrderParameters_Modified_Broyden_in_ComplexSpace(int iter);

    double Omega(int i);

    double BandWidth;
    double nia_t,nib_t,nic_t,n_t;
    Matrix<complex<double>> SiSjQ_, SiSjQ_Mean_, SiSjQ_square_Mean_;
    Matrix<complex<double>> Transformation;
    Matrix<double> SiSj_Mean_, SiSj_square_Mean_;
    double Nematic_order_mean_, Nematic_order_square_mean_;
    Parameters& Parameters_;
    Coordinates& Coordinates_;
    MFParams& MFParams_;
    Hamiltonian& Hamiltonian_;
    int lx_,ly_,ns_;
    double dosincr_,tpi_;
    vector<double> nia_,nib_,nic_;
    Matrix<double> SiSj_,dos;
    vector<double> sx_,sy_,sz_;
    Mat_3_Complex_doub OParams_;

    // Declare Fields
    Matrix<double> Sz_obs, Sx_obs, Sy_obs;
    Matrix<double> Local_density_obs;
    double Error_OP_;
    double AVG_Total_Energy, AVG_Total_Energy_sqr;


    //Declare Broyden_Mixing vectors
    vector<double> F_n; //F_n=x_n_out - x_n_in []
    vector<double> F_nm1;
    vector<double> DeltaF_n; //DeltaF_n=F_n - F-nm1;
    vector<double> Delta_x_n; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_doub Jinv_n;
    Mat_2_doub Jinv_np1;


    //Declarations for Modified Broyden Method in Complex Space
    Mat_2_Complex_doub _Delta_F;
    Mat_2_Complex_doub _u;
    Mat_1_Complex_doub _Delta_n;
    Mat_2_Complex_doub _A;




};
/*
 * ***********
 *  Functions in Class Observables ------
 *  ***********
*/




void Observables::Get_OrderParameters_diffs(){

    Error_OP_=0.0;

    for(int site=0;site<ns_;site++){
        for(int state_i=0;state_i<6;state_i++){
            for(int state_j=0;state_j<6;state_j++){
                Error_OP_ += abs((OParams_[site][state_i][state_j] -
                                 MFParams_.OParams[site][state_i][state_j])*conj(OParams_[site][state_i][state_j] -
                                          MFParams_.OParams[site][state_i][state_j]));

            }}}

    Error_OP_ +=(Parameters_.mus - Parameters_.mu_old)*(Parameters_.mus - Parameters_.mu_old);

    Error_OP_=sqrt(Error_OP_);

}

void Observables::Update_OrderParameters(int iter){


    //Simple mixing
    double alpha_OP=Parameters_.alpha_OP;

    if(Parameters_.Simple_Mixing==true){

        if(iter==0){
            cout<<"Using Simple Mixing to gain Self-Consistency"<<endl;
        }

        for(int site=0;site<ns_;site++){
            for(int state_i=0;state_i<6;state_i++){
                for(int state_j=0;state_j<6;state_j++){

                    MFParams_.OParams[site][state_i][state_j] = (1-alpha_OP)*MFParams_.OParams[site][state_i][state_j]
                            + alpha_OP*OParams_[site][state_i][state_j];
                }}}

        Parameters_.mu_old = (1-alpha_OP)*Parameters_.mu_old + alpha_OP*Parameters_.mus;

    }



    /*
    //Declared in initialize function();
    vector<double> F_n; //F_n=x_n_out - x_n_in []
    vector<double> F_nm1;
    vector<double> DeltaF_n; //DeltaF_n=F_n - F-nm1;
    vector<double> Delta_x_n; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_doub Jinv_n;
    Mat_2_doub Jinv_np1;
     */

    vector<double> vec_V, vec_U;
    double Denominator_;
    vector<double> vec_L;
    int site;
    complex<double> temp_comp;


    if(Parameters_.Broyden_Mixing==true){
        //assert(false);

        if(iter==0){
            cout<<"Using Broyden Mixing to gain Self-Consistency"<<endl;


            //Get Jinv_np1
            for(int i=0;i<(2*6*6*ns_)+1;i++){
                for(int j=0;j<(2*6*6*ns_)+1;j++){
                    if(i==j){
                        Jinv_np1[i][j]=-1.0*alpha_OP;
                    }
                    else{
                        Jinv_np1[i][j]=0.0;
                    }
                }}

            //Get F_n
            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site=Coordinates_.Nc(i,j);

                    for(int state_i=0;state_i<6;state_i++){
                        for(int state_j=0;state_j<6;state_j++){
                            F_n[site*36 + state_i*6 + state_j]=OParams_[site][state_i][state_j].real()
                                    - MFParams_.OParams[site][state_i][state_j].real();
                            F_n[36*ns_ + site*36 + state_i*6 + state_j]=OParams_[site][state_i][state_j].imag()
                                    - MFParams_.OParams[site][state_i][state_j].imag();

                        }
                    }
                }
            }
            F_n[ns_*72]=(Parameters_.mus - Parameters_.mu_old)*1.0;

            for(int i=0;i<(2*6*6*ns_)+1;i++){
                Delta_x_n[i] =0.0;
                for(int j=0;j<(2*6*6*ns_)+1;j++){
                    Delta_x_n[i] +=  -1.0*Jinv_np1[i][j]*F_n[j];
                }
            }



            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site=Coordinates_.Nc(i,j);
                    for(int state_i=0;state_i<6;state_i++){
                        for(int state_j=0;state_j<6;state_j++){
                            temp_comp.real(Delta_x_n[site*36 + state_i*6 + state_j]);
                            temp_comp.imag(Delta_x_n[36*ns_ + site*36 + state_i*6 + state_j]);
                            MFParams_.OParams[site][state_i][state_j] +=(one_complex)*temp_comp;

                        }
                    }
                }
            }
            Parameters_.mu_old +=Delta_x_n[ns_*72];


            //Copy Jinv_np1 to Jinv_n
            Jinv_n = Jinv_np1;

            //Copy F_n to F_nm1
            F_nm1=F_n;

        }

        else{
            //Get F_n
            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site=Coordinates_.Nc(i,j);

                    for(int state_i=0;state_i<6;state_i++){
                        for(int state_j=0;state_j<6;state_j++){
                            F_n[site*36 + state_i*6 + state_j]=OParams_[site][state_i][state_j].real()
                                    - MFParams_.OParams[site][state_i][state_j].real();
                            F_n[36*ns_ + site*36 + state_i*6 + state_j]=OParams_[site][state_i][state_j].imag()
                                    - MFParams_.OParams[site][state_i][state_j].imag();

                        }
                    }
                }
            }
            F_n[ns_*72]=(Parameters_.mus - Parameters_.mu_old)*1.0;

            //Get DeltaF_n
            for (int i=0;i<(2*6*6*ns_)+1;i++){
                DeltaF_n[i] = F_n[i] - F_nm1[i];
            }

            //Get vec_V = Jinv_n*DeltaF_n
            vec_V.clear();
            vec_V.resize(2*6*6*ns_ + 1);

            for(int i=0;i<2*6*6*ns_ + 1;i++){
                vec_V[i] =0.0;
                for(int j=0;j<2*6*6*ns_ + 1;j++){
                    vec_V[i] += Jinv_n[i][j]*DeltaF_n[j];
                }
            }

            //Get vec_U = Delta_x_n^dagg*Jinv_n
            vec_U.clear();
            vec_U.resize(2*6*6*ns_ + 1);

            for(int i=0;i<2*6*6*ns_ + 1;i++){
                vec_U[i] =0.0;
                for(int j=0;j<2*6*6*ns_ + 1;j++){
                    vec_U[i] += Delta_x_n[j]*Jinv_n[j][i];
                }
            }

            // Get Denominator_=<Delta_x_n|vec_V>
            Denominator_=0.0;
            for(int i=0;i<2*6*6*ns_ + 1;i++){
                Denominator_ +=Delta_x_n[i]*vec_V[i];
            }


            //Get vec_L=  Delta_x_n - vec_V;
            vec_L.clear();
            vec_L.resize(2*6*6*ns_ + 1);
            for(int i=0;i<2*6*6*ns_ + 1;i++){
                vec_L[i] = Delta_x_n[i] - vec_V[i];
            }


            //Get Mat_Temp [Remember to clear later on];
            Mat_2_doub Mat_Temp;
            Mat_Temp.resize(2*6*6*ns_ + 1);
            for(int i=0;i<2*6*6*ns_ + 1;i++){
                Mat_Temp[i].resize(2*6*6*ns_ + 1);
                for(int j=0;j<2*6*6*ns_ + 1;j++){
                    Mat_Temp[i][j] = (vec_L[i]*vec_U[j])/(Denominator_);
                }
            }


            //Get Jinv_np1

            for(int i=0;i<2*6*6*ns_ + 1;i++){
                for(int j=0;j<2*6*6*ns_ + 1;j++){
                    Jinv_np1[i][j]  = Jinv_n[i][j]  + Mat_Temp[i][j];
                }
            }

            for(int i=0;i<2*6*6*ns_ + 1;i++){
                Delta_x_n[i] =0.0;
                for(int j=0;j<2*6*6*ns_ + 1;j++){
                    Delta_x_n[i] +=  -1.0*Jinv_np1[i][j]*F_n[j]; //HERE
                }
            }



            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site=Coordinates_.Nc(i,j);
                    for(int state_i=0;state_i<6;state_i++){
                        for(int state_j=0;state_j<6;state_j++){
                            temp_comp.real(Delta_x_n[site*36 + state_i*6 + state_j]);
                            temp_comp.imag(Delta_x_n[36*ns_ + site*36 + state_i*6 + state_j]);
                            MFParams_.OParams[site][state_i][state_j] +=(one_complex)*(temp_comp);

                        }
                    }
                }
            }
            Parameters_.mu_old +=Delta_x_n[ns_*72];

            //Copy Jinv_np1 to Jinv_n
            Jinv_n = Jinv_np1;

            //Copy F_n to F_nm1
            F_nm1=F_n;


            //Clear Mat_Temp
            for(int i=0;i<2*6*6*ns_ + 1;i++){
                Mat_Temp[i].clear();
            }
            Mat_Temp.clear();

        }



        //Checking MFParams_.OParams symmetries
        /*
        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);
                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                    if(state_i==state_j){
                    if(MFParams_.OParams[site][state_i][state_j].imag()!=0.0){
                        cout<<"local densities order params have problem"<<endl;
                    }
                      assert(MFParams_.OParams[site][state_i][state_j].imag()==0.0);
                    }
                    else{
                    if(
                      abs(MFParams_.OParams[site][state_i][state_j] - conj(MFParams_.OParams[site][state_j][state_i]))>0.00001
                       ){
                        cout<<"Hermiticity of order params is an issue"<<endl;

                        cout<<site<<"   "<<state_i<<"    "<<state_j<<"    "<<MFParams_.OParams[site][state_i][state_j]<<endl;
                        cout<<site<<"   "<<state_j<<"    "<<state_i<<"    "<<MFParams_.OParams[site][state_j][state_i]<<endl;
                    }
                    assert(
                          abs(MFParams_.OParams[site][state_i][state_j] - conj(MFParams_.OParams[site][state_j][state_i]))<0.00001
                          );
                    }

                    }}}}
                    */


        //Maintain hermiticity of order params, because after many iterations they tend to loose it
        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);
                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=state_i;state_j<6;state_j++){
                        if(state_j>state_i){
                MFParams_.OParams[site][state_i][state_j]=conj(MFParams_.OParams[site][state_j][state_i]);
                        }
                        else{
                            assert(state_j==state_i);
                            MFParams_.OParams[site][state_i][state_j].imag(0.0);
                        }
                    }}}}


    }


}



void Update_OrderParameters_Modified_Broyden_in_ComplexSpace(int iter){

//Modified Broyden is used from "D. D. Johnson, Phys. Rev. B 38, 12807, 1988".
//Look into this "https://arxiv.org/pdf/0805.4446.pdf" as well.
//Please note the m=iter+1, because m starts from "1".







}

void Observables::Calculate_Local_Density(){

    int c1, c2;
    int site;
    complex<double> value;
    vector<double> avg_den_jm, avg_den_t2g;
    avg_den_jm.resize(6);avg_den_t2g.resize(6);

    ofstream local_density_out_file("local_density_out.txt");

    local_density_out_file<<"#ix   iy   3by2_m3by2 3by2_3by2   3by2_m1by2   3by2_1by2    1by2_m1by2    1by2_1by2   yz_up    xz_up      xy_up    yz_dn    xz_dn      xy_dn"<<endl;



    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            local_density_out_file<<i<<setw(10)<<j;
            site=Coordinates_.Nc(i,j);

            for(int state=0;state<6;state++){
                c1=Coordinates_.Nc_dof(site,state);
                value=zero_complex;
                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    value+=( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                             (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                             );

                }
                local_density_out_file<<setw(15)<<value.real();
                avg_den_jm[state] +=value.real();
            }


            for(int t2g_type=0;t2g_type<6;t2g_type++){

                value=zero_complex;
                for(int jm_type1=0;jm_type1<6;jm_type1++){
                    c1=Coordinates_.Nc_dof_(site,jm_type1);

                    for(int jm_type2=0;jm_type2<6;jm_type2++){
                        c2=Coordinates_.Nc_dof_(site,jm_type2);

                        if( (Transformation(t2g_type,jm_type1) != zero_complex)
                                &&
                                (Transformation(t2g_type,jm_type2) != zero_complex) ){

                            for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                                value += conj(Hamiltonian_.Ham_(c2,n))*conj(Transformation(t2g_type,jm_type1))*
                                        Hamiltonian_.Ham_(c1,n)*Transformation(t2g_type,jm_type2)*
                                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                            }

                        }
                    }
                }
                local_density_out_file<<setw(15)<<value.real();
                avg_den_t2g[t2g_type] +=value.real();

            }


            local_density_out_file<<endl;

        }
        local_density_out_file<<endl;
    }

    local_density_out_file<<"#Average densities in same order"<<endl;
    local_density_out_file<<"#"<<setw(15);
    for(int state=0;state<6;state++){
        local_density_out_file<<avg_den_jm[state]/(1.0*ns_)<<setw(15);
    }
    for(int state=0;state<6;state++){
        local_density_out_file<<avg_den_t2g[state]/(1.0*ns_)<<setw(15);
    }
    local_density_out_file<<endl;


}



void Observables::Calculate_IPR(){
    /*
    double IPR;
    int c1, site;
    double eta = 0.001;
    int n_chosen=(Parameters_.ns*Parameters_.Fill*2.0) - 1;
    IPR=0.0;
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;
                //  IPR += abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                //         abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                //         Lorentzian( Parameters_.mus - Hamiltonian_.eigs_[n], eta);

                IPR += abs(Hamiltonian_.Ham_(c1,n_chosen))*abs(Hamiltonian_.Ham_(c1,n_chosen))*
                        abs(Hamiltonian_.Ham_(c1,n_chosen))*abs(Hamiltonian_.Ham_(c1,n_chosen));




            }
        }
    }

    cout<<"IPR for state no. "<<n_chosen<<" = "<<IPR<<endl;



    IPR=0.0;
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){


                    IPR += abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                            abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                            Lorentzian( Parameters_.mus - Hamiltonian_.eigs_[n], eta);

                }
            }
        }
    }

    cout<<"IPR for near mu(with eta = "<<eta<<") = "<<IPR<<endl;




    string fileout_FermiState="Fermi_state_probability.txt";
    ofstream file_FermiState_out(fileout_FermiState.c_str());
    file_FermiState_out<<"#ix   iy   site   |Psi_{Fermi,up}(ix,iy)|^2 + |Psi_{Fermi,dn}(ix,iy)|^2"<<endl;

    double value;
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);
            value = 0.0;

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;

                value += abs(Hamiltonian_.Ham_(c1,n_chosen))*abs(Hamiltonian_.Ham_(c1,n_chosen));
            }

            file_FermiState_out<<i<<"\t"<<j<<"\t"<<site<<"\t"<<value<<endl;

        }
        file_FermiState_out<<endl;
    }





    string fileout_Near_mu="Near_mu_probability.txt";
    ofstream file_Near_mu_out(fileout_Near_mu.c_str());
    file_Near_mu_out<<"#ix   iy   site   sum_{n}Lorentz(near_mu)*|Psi_{n,up}(ix,iy)|^2 + |Psi_{n,dn}(ix,iy)|^2"<<endl;


    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);
            value = 0.0;

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;
                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    value += abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                            Lorentzian(Parameters_.mus - Hamiltonian_.eigs_[n], eta);;
                }
            }

            file_Near_mu_out<<i<<"\t"<<j<<"\t"<<site<<"\t"<<value<<endl;

        }
        file_Near_mu_out<<endl;
    }



*/
}

void Observables::Calculate_Order_Params(){

    int ci,cj;
    OParams_.resize(ns_);
    for(int site=0;site<ns_;site++){
        OParams_[site].resize(6);
        for(int state=0;state<6;state++){
            OParams_[site][state].resize(6);
        }
    }

    for(int site=0;site<ns_;site++){
        for(int state_i=0;state_i<6;state_i++){
            ci=Coordinates_.Nc_dof(site,state_i);

            for(int state_j=0;state_j<6;state_j++){
                cj=Coordinates_.Nc_dof(site,state_j);

                OParams_[site][state_i][state_j]=zero_complex;
                for(int n=0;n<Coordinates_.no_dof_;n++){
                    OParams_[site][state_i][state_j] += conj(Hamiltonian_.Ham_(cj,n))*Hamiltonian_.Ham_(ci,n)*
                            (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mu_old)*Parameters_.beta ) + 1.0));
                }

            }
        }
    }


    //checking Hermiticity of order parameters
    for(int site=0;site<ns_;site++){
        for(int state_i=0;state_i<6;state_i++){
            for(int state_j=state_i;state_j<6;state_j++){
                if(state_i==state_j){
                    assert(OParams_[site][state_i][state_j].imag()==0.0);
                }
                else{
                    assert(OParams_[site][state_i][state_j].imag()==(-1.0*OParams_[site][state_j][state_i].imag()));
                    assert(OParams_[site][state_i][state_j].real()==(1.0*OParams_[site][state_j][state_i].real()));
                }
            }

        }
    }


}

void Observables::Calculate_Akw(){
    /*

    //---------Read from input file-----------------------//
    string fileout="Akw.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.001;
    omega_min=-1.6;omega_max=2.6;d_omega=0.03;
    //---------------------------------------------------//


    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Akw_out(fileout.c_str());

    int c1,c2;

    Mat_3_Complex_doub A_up_00;
    Mat_3_Complex_doub A_dn_00;
    A_up_00.resize(Parameters_.ns);
    A_dn_00.resize(Parameters_.ns);

    for (int i=0;i<Parameters_.ns;i++){
        A_up_00[i].resize(Parameters_.ns);
        A_dn_00[i].resize(Parameters_.ns);

        for(int j=0;j<Parameters_.ns;j++){
            A_up_00[i][j].resize(omega_index_max);
            A_dn_00[i][j].resize(omega_index_max);
        }
    }


    complex<double> Nup_check(0,0);
    complex<double> Ndn_check(0,0);

    for (int j=0;j<Parameters_.ns;j++){
        for (int l=0;l<Parameters_.ns;l++){
            for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
                A_up_00[j][l][omega_ind]=zero_complex;
                A_dn_00[j][l][omega_ind]=zero_complex;

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

                    //c= l + or1*ns_ + ns_*orbs_*spin;

                    //Hamiltonian_.Ham_(c2,n) is nth eigenvector and c2th component [checked];
                    c1 = l + ns_; c2 = j+ ns_;
                    A_dn_00[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = l; c2 = j;
                    A_up_00[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                }

                if(j==l){
                    Nup_check += (A_up_00[j][l][omega_ind])*d_omega;
                    Ndn_check += (A_dn_00[j][l][omega_ind])*d_omega;
                }
            }
        }
    }

    cout << "Nup_check = "<<Nup_check<<endl;
    cout << "Ndn_check = "<<Ndn_check<<endl;

    complex<double> temp_up_00;
    complex<double> temp_dn_00;
    double kx,ky;
    int kx_i,ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;

    //--------\Gamma to X-----------------
    ky_i=0;
    for(kx_i=0;kx_i<=(Parameters_.lx/2);kx_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------X to M-----------------
    kx_i=(Parameters_.lx/2);
    for(ky_i=1;ky_i<=(Parameters_.lx/2);ky_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------M to \Gamma[with one extra point,
    //                  because in gnuplor use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i=(Parameters_.lx/2) - 1;
    ky_i=(Parameters_.lx/2) - 1;
    for(kx_i=(Parameters_.lx/2) - 1;kx_i>=-1;kx_i--){
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    double k22_offset=0;
    for(int k_point=0;k_point<k_path.size();k_point++){

        kx_i=k_path[k_point].first;
        ky_i=k_path[k_point].second;
        kx=(2.0*PI*kx_i)/(1.0*Parameters_.lx);
        ky=(2.0*PI*ky_i)/(1.0*Parameters_.ly);

        for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
            temp_up_00=zero_complex;
            temp_dn_00=zero_complex;

            for(int j=0;j<ns_;j++){
                for(int l=0;l<ns_;l++){
                    temp_up_00 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_up_00[j][l][omega_ind];

                    temp_dn_00 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_dn_00[j][l][omega_ind];

                }
            }
            //Use 1:6:7----for gnuplot
            file_Akw_out<< k_point<<"   "<<kx_i<<"   "<<ky_i<<"   "<<(ky_i*Parameters_.lx) + kx_i<<"    "<<
                           omega_min + (d_omega*omega_ind)<<"   "<<omega_ind<<"    "
                        <<temp_up_00.real()<<"    "<<temp_dn_00.real()<<"    "
                       <<temp_up_00.imag()<<"     "<<temp_dn_00.imag()<<"    "<<endl;

        }
        file_Akw_out<<endl;
    }


*/
}


void Observables::Calculate_Nw(){

    //---------Read from input file-----------------------//
    string fileout="Nw_jm.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.02;
    omega_min=-20;omega_max=20.0;d_omega=0.005;
    //---------------------------------------------------//

    int c1;
    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );


    //---------------------------------------------------------------------------------//
    //************************************Nw_jm****************************************//
    //---------------------------------------------------------------------------------//

    ofstream file_Nw_out(fileout.c_str());
    double temp_val ;

    file_Nw_out<<"#(w-mu)    jm_3by2_m3by2     jm_3by2_3by2     ";
    file_Nw_out<<"jm_3by2_m1by2       jm_3by2_1by2     jm_1by2_m1by2   jm_1by2_1by2"<<endl;

    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
        file_Nw_out<<omega_min + (omega_ind*d_omega)<<"       ";


        for(int state_type=0;state_type<6;state_type++){
            temp_val=0.0;

            for(int site=0;site<Coordinates_.ns_;site++){
                c1=Coordinates_.Nc_dof_(site,state_type);

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    temp_val +=  (conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                                  Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta)).real();

                }
            }
            file_Nw_out<<temp_val<<"      ";
        }
        file_Nw_out<<endl;
    }

    file_Nw_out<<"#actual mu = "<< Parameters_.mus<<", but shifted to 0"<<endl;

    //---------------------------------------------------------------------------------//
    //*********************************************************************************//



    //---------------------------------------------------------------------------------//
    //************************************Nw_jm****************************************//
    //---------------------------------------------------------------------------------//

    int c2;

    string fileout_t2g="Nw_t2g.txt";
    ofstream file_Nw_out_t2g(fileout_t2g.c_str());

    file_Nw_out_t2g<<"#(w-mu)    yz_up    xz_up      xy_up     ";
    file_Nw_out_t2g<<"yz_dn    xz_dn      xy_dn"<<endl;

    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
        file_Nw_out_t2g<<omega_min + (omega_ind*d_omega)<<"       ";


        for(int t2g_type=0;t2g_type<6;t2g_type++){
            temp_val=0.0;

            for(int site=0;site<Coordinates_.ns_;site++){

                for(int jm_type1=0;jm_type1<6;jm_type1++){
                    c1=Coordinates_.Nc_dof_(site,jm_type1);

                    for(int jm_type2=0;jm_type2<6;jm_type2++){
                        c2=Coordinates_.Nc_dof_(site,jm_type2);

                        if( (Transformation(t2g_type,jm_type1) != zero_complex)
                                &&
                                (Transformation(t2g_type,jm_type2) != zero_complex) ){

                            for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                                temp_val += ( conj(Hamiltonian_.Ham_(c2,n))*conj(Transformation(t2g_type,jm_type1))*
                                              Hamiltonian_.Ham_(c1,n)*Transformation(t2g_type,jm_type2)*
                                              Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta)).real();

                            }
                        }

                    }
                }

            }
            file_Nw_out_t2g<<temp_val<<"      ";
        }
        file_Nw_out_t2g<<endl;
    }

    file_Nw_out_t2g<<"#actual mu = "<< Parameters_.mus<<", but shifted to 0"<<endl;

    //---------------------------------------------------------------------------------//
    //*********************************************************************************//






    int n_chosen=Parameters_.Total_Particles - 1;
    cout<<"Gap = "<<Hamiltonian_.eigs_[n_chosen+1] - Hamiltonian_.eigs_[n_chosen]<<endl;

    string fileout_Eigen="Eigen_spectrum.txt";
    ofstream file_Eigen_out(fileout_Eigen.c_str());

    for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
        file_Eigen_out<<n<<"\t"<<Hamiltonian_.eigs_[n]<<endl;
    }


}


void Observables::Calculate_Optical_Conductivity(){
    /*

    string fileout_sigma_w = "Optical_conductivity.txt";
    ofstream file_sigma_w_out(fileout_sigma_w.c_str());
    file_sigma_w_out<<"#omega   Sigma_xx  Sigma_yy"<<endl;

    //--------------------------------------------------//
    double omega_min, omega_max, d_omega;
    double eta = 0.05;
    omega_min=0.00001;omega_max=10.0;d_omega=0.001;
    //---------------------------------------------------//


    Mat_2_doub PSI_x, PSI_y;
    complex<double> value_x, value_y;
    int ipx, ipy;




    PSI_x.resize(2*ns_);
    PSI_y.resize(2*ns_);
    for(int n=0;n<2*ns_;n++){
        PSI_x[n].resize(2*ns_);
        PSI_y[n].resize(2*ns_);
    }



    for(int n=0;n<2*ns_;n++){
        for(int m=0;m<2*ns_;m++){

            value_x=zero_complex;
            value_y=zero_complex;

            for(int i=0;i<ns_;i++){
                ipx = Coordinates_.neigh(i,0);
                ipy = Coordinates_.neigh(i,2);



                for(int spin=0;spin<2;spin++){
                    value_x += ( conj(Hamiltonian_.Ham_(ipx + (ns_*spin),n))*Hamiltonian_.Ham_(i + (ns_*spin),m) )
                            -  ( conj(Hamiltonian_.Ham_(i + (ns_*spin),n))*Hamiltonian_.Ham_(ipx + (ns_*spin),m) );

                    value_y += ( conj(Hamiltonian_.Ham_(ipy + (ns_*spin),n))*Hamiltonian_.Ham_(i + (ns_*spin),m) )
                            -  ( conj(Hamiltonian_.Ham_(i + (ns_*spin),n))*Hamiltonian_.Ham_(ipy + (ns_*spin),m) );

                }


            }


            PSI_x[n][m] = abs(value_x)*abs(value_x);
            PSI_y[n][m] = abs(value_y)*abs(value_y);


        }
    }




    double sigma_x, sigma_y;
    double omega_val;


    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );


    //cout<<"omega_index_max = "<<omega_index_max<<endl;
    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
        omega_val = omega_min + (omega_ind*d_omega);

        //cout<<omega_ind<<endl;

        sigma_x=0.0; sigma_y=0.0;
        for(int n=0;n<2*ns_;n++){
            for(int m=0;m<2*ns_;m++){

                if(n!=m){
                    sigma_x += (PSI_x[n][m])
                            *((1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)))
                            *((1.0/( exp((Parameters_.mus-Hamiltonian_.eigs_[m])*Parameters_.beta ) + 1.0)))
                            *Lorentzian( omega_min + (omega_ind*d_omega) + Hamiltonian_.eigs_[n] - Hamiltonian_.eigs_[m], eta);

                    sigma_y += (PSI_y[n][m])
                            *((1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)))
                            *((1.0/( exp((Parameters_.mus-Hamiltonian_.eigs_[m])*Parameters_.beta ) + 1.0)))
                            *Lorentzian( omega_min + (omega_ind*d_omega) + Hamiltonian_.eigs_[n] - Hamiltonian_.eigs_[m], eta);


                }
            }
        }

        sigma_x = sigma_x*PI*(1.0 - exp(-1.0*Parameters_.beta*(omega_val)))*(1.0/(omega_val*ns_));
        sigma_y = sigma_y*PI*(1.0 - exp(-1.0*Parameters_.beta*(omega_val)))*(1.0/(omega_val*ns_));

        file_sigma_w_out<<omega_val<<"     "<<sigma_x<<"     "<<sigma_y<<endl;

    }



*/
}


void Observables::Calculate_SpinSpincorrelations(){
    /*

    string S2_out = "Local_S2.txt";
    ofstream file_S2_out(S2_out.c_str());
    file_S2_out<<"#site_i    ix   iy   S^2[site_i]"<<endl;


    string Sq_out = "Sq.txt";
    ofstream file_Sq_out(Sq_out.c_str());
    file_Sq_out<<"#qx  qy   qx_index    qy_index   S(qx,qy)"<<endl;

    string SSr_out = "SSr.txt";
    ofstream file_SSr_out(SSr_out.c_str());
    file_SSr_out<<"#site_i   site_i(x)    site_i(y)    site_j   site_j(x)    site_j(y)  SS[site_i][site_j]"<<endl;

    int c_iup, c_idn, c_jup, c_jdn;

    Mat_2_Complex_doub SS_nm;
    SS_nm.resize(Hamiltonian_.Ham_.n_row());
    for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
        SS_nm[n].resize(Hamiltonian_.Ham_.n_row());
    }


    Mat_2_Complex_doub SS_ri_rj;
    SS_ri_rj.resize(ns_);
    for(int site_i=0;site_i<ns_;site_i++){
        SS_ri_rj[site_i].resize(ns_);
    }


    for(int site_i=0;site_i<ns_;site_i++){
        for(int site_j=0;site_j<ns_;site_j++){

            c_iup=site_i;
            c_jup=site_j;
            c_idn=site_i + ns_;
            c_jdn=site_j + ns_;

            SS_ri_rj[site_i][site_j]=zero_complex;
            for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                for(int m=0;m<Hamiltonian_.Ham_.n_row();m++){

                    if(n==m){

                        //SzSz*f(En)*..
                        SS_ri_rj[site_i][site_j] +=(0.25*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,n))*Hamiltonian_.Ham_(c_jup,n) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,n))*Hamiltonian_.Ham_(c_jdn,n) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,n))*Hamiltonian_.Ham_(c_jdn,n) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,n))*Hamiltonian_.Ham_(c_jup,n) )

                                    );

                        //0.5*(S+S- + S-S+)*f(En)*..
                        SS_ri_rj[site_i][site_j] +=(0.5*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,n))*Hamiltonian_.Ham_(c_jdn,n) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,n))*Hamiltonian_.Ham_(c_jup,n) )
                                    );



                    }

                    else{

                        //SzSz*f(En)*f(Em)..
                        SS_ri_rj[site_i][site_j] +=(0.25*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jup,m) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jdn,m) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jdn,m) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jup,m) )

                                    );


                        //SzSz*f(En)*(1.0 - f(Em))..
                        SS_ri_rj[site_i][site_j] +=(0.25*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (1.0 - (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0)) )*
                                (
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,m)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jup,n) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,m)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jdn,n) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,m)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jdn,n) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,m)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jup,n) )

                                    );


                        //0.5*(S+S- + S-S+)*f(En)*f(Em)..
                        SS_ri_rj[site_i][site_j] +=(0.5*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jdn,m) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jup,m) )

                                    );

                        //0.5*(S+S- + S-S+)*f(En)*(1-f(Em))..
                        SS_ri_rj[site_i][site_j] +=(0.5*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0)) )*
                                (
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_iup,m)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jdn,n) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_idn,m)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jup,n) )

                                    );


                    }


                }



            }



            cout<<site_i<<"\t"<<site_j<<" done"<<endl;
        }

    }


    int site_i, site_j;

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            site_i=Coordinates_.Nc(ix,iy);

            for(int jx=0;jx<lx_;jx++){
                for(int jy=0;jy<ly_;jy++){
                    site_j=Coordinates_.Nc(jx,jy);
                    file_SSr_out<<site_i<<"\t"<<ix<<"\t"<<iy<<"\t"<<site_j<<"\t"<<
                                  jx<<"\t"<<jy<<"\t"<<
                                  real(SS_ri_rj[site_i][site_j])<<"\t"<<
                                  imag(SS_ri_rj[site_i][site_j])<<endl;
                }
            }

        }
    }


    double qx_, qy_;
    complex<double> value;



    for(int qx=0;qx<lx_;qx++){
        for(int qy=0;qy<ly_;qy++){
            value = zero_complex;
            qx_ = (2.0*qx*PI)*(1.0/(1.0*lx_));
            qy_ = (2.0*qy*PI)*(1.0/(1.0*ly_));

            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){
                    site_i=Coordinates_.Nc(ix,iy);

                    for(int jx=0;jx<lx_;jx++){
                        for(int jy=0;jy<ly_;jy++){
                            site_j=Coordinates_.Nc(jx,jy);

                            value += one_complex*(exp(iota_complex*( (qx_*(ix - jx)) + (qy_*(iy-jy)) )))*
                                    SS_ri_rj[site_i][site_j]*(1.0/(1.0*ns_));

                        }}
                }}

            file_Sq_out<<qx_<<"\t"<<qy_<<"\t"<<qx<<"\t"<<qy<<"\t"<<real(value)<<"\t"<<imag(value)<<endl;

        }
    }



    complex<double> Avg_S2=0.0;

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            site_i=Coordinates_.Nc(ix,iy);
            file_S2_out<<site_i<<"\t"<<ix<<"\t"<<iy<<"\t"<<real(SS_ri_rj[site_i][site_i])
                      <<"\t"<<imag(SS_ri_rj[site_i][site_i])<<endl;

            Avg_S2 += one_complex*(SS_ri_rj[site_i][site_i]);
        }}

    cout<<"Avg Local Moment (S^2) = "<<real(Avg_S2)/(1.0*ns_)<<"\t"<<imag(Avg_S2)/(1.0*ns_)<<endl;


*/
}


void Observables::Calculate_SpinSpincorrelations_Smartly(){
    /*

    string S2_out = "Local_S2.txt";
    ofstream file_S2_out(S2_out.c_str());
    file_S2_out<<"#site_i    ix   iy   S^2[site_i]"<<endl;


    string Sq_out = "Sq.txt";
    ofstream file_Sq_out(Sq_out.c_str());
    file_Sq_out<<"#qx  qy   qx_index    qy_index   S(qx,qy)"<<endl;

    string SSr_out = "SSr.txt";
    ofstream file_SSr_out(SSr_out.c_str());
    file_SSr_out<<"#site_i   site_i(x)    site_i(y)    site_j   site_j(x)    site_j(y)  SS[site_i][site_j]"<<endl;

    int c_iup, c_idn, c_jup, c_jdn;

    Mat_2_Complex_doub UP_UP_Fermi, DOWN_DOWN_Fermi, UP_DOWN_Fermi, DOWN_UP_Fermi ;
    Mat_2_Complex_doub UP_UP_1mFermi, DOWN_DOWN_1mFermi, UP_DOWN_1mFermi, DOWN_UP_1mFermi ;

    UP_UP_Fermi.resize(ns_);DOWN_DOWN_Fermi.resize(ns_);UP_DOWN_Fermi.resize(ns_);DOWN_UP_Fermi.resize(ns_);
    UP_UP_1mFermi.resize(ns_);DOWN_DOWN_1mFermi.resize(ns_);UP_DOWN_1mFermi.resize(ns_);DOWN_UP_1mFermi.resize(ns_);

    for(int n=0;n<ns_;n++){
        UP_UP_Fermi[n].resize(ns_);DOWN_DOWN_Fermi[n].resize(ns_);
        UP_DOWN_Fermi[n].resize(ns_);DOWN_UP_Fermi[n].resize(ns_);

        UP_UP_1mFermi[n].resize(ns_);DOWN_DOWN_1mFermi[n].resize(ns_);
        UP_DOWN_1mFermi[n].resize(ns_);DOWN_UP_1mFermi[n].resize(ns_);
    }


    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){

            UP_UP_Fermi[i][j] = zero_complex;
            DOWN_DOWN_Fermi[i][j] = zero_complex;
            UP_DOWN_Fermi[i][j] = zero_complex;
            DOWN_UP_Fermi[i][j] = zero_complex;

            UP_UP_1mFermi[i][j] = zero_complex;
            DOWN_DOWN_1mFermi[i][j] = zero_complex;
            UP_DOWN_1mFermi[i][j] = zero_complex;
            DOWN_UP_1mFermi[i][j] = zero_complex;


            for(int n=0;n<2*ns_;n++){
                UP_UP_Fermi[i][j] +=  conj(Hamiltonian_.Ham_(i,n))*Hamiltonian_.Ham_(j,n)*
                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                DOWN_DOWN_Fermi[i][j] +=  conj(Hamiltonian_.Ham_(i+ns_,n))*Hamiltonian_.Ham_(j+ns_,n)*
                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                UP_DOWN_Fermi[i][j] +=  conj(Hamiltonian_.Ham_(i,n))*Hamiltonian_.Ham_(j+ns_,n)*
                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                DOWN_UP_Fermi[i][j] +=  conj(Hamiltonian_.Ham_(i+ns_,n))*Hamiltonian_.Ham_(j,n)*
                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                UP_UP_1mFermi[i][j] +=  conj(Hamiltonian_.Ham_(i,n))*Hamiltonian_.Ham_(j,n)*
                        ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                DOWN_DOWN_1mFermi[i][j] +=  conj(Hamiltonian_.Ham_(i+ns_,n))*Hamiltonian_.Ham_(j+ns_,n)*
                        ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                UP_DOWN_1mFermi[i][j] +=  conj(Hamiltonian_.Ham_(i,n))*Hamiltonian_.Ham_(j+ns_,n)*
                        ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                DOWN_UP_1mFermi[i][j] +=  conj(Hamiltonian_.Ham_(i+ns_,n))*Hamiltonian_.Ham_(j,n)*
                        ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

            }
        }
    }

    Mat_2_Complex_doub SS_ri_rj;
    SS_ri_rj.resize(ns_);
    for(int site_i=0;site_i<ns_;site_i++){
        SS_ri_rj[site_i].resize(ns_);
    }


    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){

            SS_ri_rj[i][j]=zero_complex;

            //SzSz..
            SS_ri_rj[i][j] +=(0.25*one_complex)*(

                        (UP_UP_Fermi[i][j]*UP_UP_1mFermi[j][i])
                        - (UP_DOWN_Fermi[i][j]*DOWN_UP_1mFermi[j][i])
                        + (DOWN_DOWN_Fermi[i][j]*DOWN_DOWN_1mFermi[j][i])
                        - (DOWN_UP_Fermi[i][j]*UP_DOWN_1mFermi[j][i])
                        + (UP_UP_Fermi[i][i]*UP_UP_Fermi[j][j])
                        - (UP_UP_Fermi[i][i]*DOWN_DOWN_Fermi[j][j])
                        + (DOWN_DOWN_Fermi[i][i]*DOWN_DOWN_Fermi[j][j])
                        - (DOWN_DOWN_Fermi[i][i]*UP_UP_Fermi[j][j])

                        );


            //0.5*(S+S- + S-S+)
            SS_ri_rj[i][j] +=(0.5*one_complex)*(

                        (UP_UP_Fermi[i][j]*DOWN_DOWN_1mFermi[j][i])
                        + (DOWN_DOWN_Fermi[i][j]*UP_UP_1mFermi[j][i])
                        + (UP_DOWN_Fermi[i][i]*DOWN_UP_Fermi[j][j])
                        + (DOWN_UP_Fermi[i][i]*UP_DOWN_Fermi[j][j])

                        );

            //cout<<i<<"\t"<<j<<" done"<<endl;
        }

    }


    int site_i, site_j;

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            site_i=Coordinates_.Nc(ix,iy);

            for(int jx=0;jx<lx_;jx++){
                for(int jy=0;jy<ly_;jy++){
                    site_j=Coordinates_.Nc(jx,jy);
                    file_SSr_out<<site_i<<"\t"<<ix<<"\t"<<iy<<"\t"<<site_j<<"\t"<<
                                  jx<<"\t"<<jy<<"\t"<<
                                  real(SS_ri_rj[site_i][site_j])<<"\t"<<
                                  imag(SS_ri_rj[site_i][site_j])<<endl;


                }
            }

        }
    }


    double qx_, qy_;
    complex<double> value;



    for(int qx=0;qx<lx_;qx++){
        for(int qy=0;qy<ly_;qy++){
            value = zero_complex;
            qx_ = (2.0*qx*PI)*(1.0/(1.0*lx_));
            qy_ = (2.0*qy*PI)*(1.0/(1.0*ly_));

            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){
                    site_i=Coordinates_.Nc(ix,iy);

                    for(int jx=0;jx<lx_;jx++){
                        for(int jy=0;jy<ly_;jy++){
                            site_j=Coordinates_.Nc(jx,jy);

                            value += one_complex*(exp(iota_complex*( (qx_*(ix - jx)) + (qy_*(iy-jy)) )))*
                                    SS_ri_rj[site_i][site_j]*(1.0/(1.0*ns_));

                        }}
                }}

            file_Sq_out<<qx_<<"\t"<<qy_<<"\t"<<qx<<"\t"<<qy<<"\t"<<real(value)<<"\t"<<imag(value)<<endl;

        }
        file_Sq_out<<endl;

    }



    complex<double> Avg_S2=0.0;

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            site_i=Coordinates_.Nc(ix,iy);
            file_S2_out<<site_i<<"\t"<<ix<<"\t"<<iy<<"\t"<<real(SS_ri_rj[site_i][site_i])
                      <<"\t"<<imag(SS_ri_rj[site_i][site_i])<<endl;

            Avg_S2 += one_complex*(SS_ri_rj[site_i][site_i]);
        }
        file_S2_out<<endl;
    }

    cout<<"Avg Local Moment (S^2) = "<<real(Avg_S2)/(1.0*ns_)<<"\t"<<imag(Avg_S2)/(1.0*ns_)<<endl;


*/
}




void Observables::Get_Non_Interacting_dispersion(){

}


double Observables::Lorentzian(double x, double brd){
    double temp;

    temp = (1.0/PI)*( (brd/2.0)/ ( (x*x) + ((brd*brd)/4.0) ) );

    return temp;

}

void Observables::DensityOfStates(){
    //-----------Calculate Bandwidth------//
    BandWidth=2.0;
    //-----------------------------------//

} // ----------


void Observables::OccDensity(){

} // ----------


void Observables::TotalOccDensity(){

} // ----------



void Observables::Total_Energy_Average(double Curr_QuantE, double CurrE){

    AVG_Total_Energy += Curr_QuantE + CurrE;
    AVG_Total_Energy_sqr += (Curr_QuantE + CurrE)*(Curr_QuantE + CurrE);
}

void Observables::OccDensity(int tlabel){

} // ----------



void Observables::Initialize(){


    F_n.resize(ns_*6*6*2 + 1); //F_n=x_n_out - x_n_in []
    F_nm1.resize(ns_*6*6*2 + 1);
    DeltaF_n.resize(ns_*6*6*2 + 1); //DeltaF_n=F_n - F-nm1;
    Delta_x_n.resize(ns_*6*6*2 + 1); //Delta_x_n= x_n_in - x_nm1_in;
    Jinv_n.resize(ns_*6*6*2 + 1);
    Jinv_np1.resize(ns_*6*6*2 + 1);



    for(int i=0;i<ns_*6*6*2 + 1;i++){
        Jinv_n[i]. resize(ns_*6*6*2 + 1);
        Jinv_np1[i]. resize(ns_*6*6*2 + 1);
    }


    Transformation.resize(6,6);
    //Saves a_{jm} ---to---> c_{spin,orbital} [NOT dagger]
    Transformation(0,1)=iota_complex*(-1.0/sqrt(2));
    Transformation(0,2)=one_complex*(1.0/sqrt(6));
    Transformation(0,4)=one_complex*(-1.0/sqrt(3));

    Transformation(1,1)=one_complex*(1.0/sqrt(2));
    Transformation(1,2)=iota_complex*(-1.0/sqrt(6));
    Transformation(1,4)=iota_complex*(1.0/sqrt(3));

    Transformation(2,3)=one_complex*(2.0/sqrt(6));
    Transformation(2,5)=one_complex*(1.0/sqrt(3));

    Transformation(3,0)=iota_complex*(1.0/sqrt(2));
    Transformation(3,3)=one_complex*(-1.0/sqrt(6));
    Transformation(3,5)=one_complex*(1.0/sqrt(3));

    Transformation(4,0)=one_complex*(1.0/sqrt(2));
    Transformation(4,3)=iota_complex*(-1.0/sqrt(6));
    Transformation(4,5)=iota_complex*(1.0/sqrt(3));

    Transformation(5,2)=one_complex*(2.0/sqrt(6));
    Transformation(5,4)=one_complex*(1.0/sqrt(3));







} // ----------


double Observables::Omega(int i){
    return -20.0+double(i)*dosincr_;
} // ----------









#endif // OBSERVABLES_H
