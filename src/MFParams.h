#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class
#define MFParams_class

class MFParams{ 
public:
    /* Convention
(site-0, 3/2,-3/2)=0 (site-0, 3/2,3/2)=1 (site-0, 3/2,-1/2)=2 (site-0, 3/2,1/2)=3
(site-0, 1/2,-1/2)=4 (site-0, 1/2,1/2)=5
(site-1, 3/2,-3/2)=6 (site-1, 3/2,3/2)=7 (site-1, 3/2,-1/2)=8 (site-1, 3/2,1/2)=9
(site-1, 1/2,-1/2)=10 (site-1, 1/2,1/2)=11
(site-2)........
*/
    // Define Fields
    Mat_3_Complex_doub OParams;
    Matrix<double> Disorder;

    // Constructor
    MFParams(Parameters& Parameters__, Coordinates&  Coordinates__, mt19937_64& Generator1__ , mt19937_64& Generator2__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__), Generator1_(Generator1__), Generator2_(Generator2__)
    {
        //setupm_arr();
        initialize();
    }


    double random1();
    double random2();
    void initialize();


    Parameters &Parameters_;
    Coordinates &Coordinates_;
    mt19937_64 &Generator1_; //for random fields
    mt19937_64 &Generator2_; //for random disorder
    int lx_,ly_,ns_, no_dof_;

    uniform_real_distribution<double> dis1_;//for random fields
    uniform_real_distribution<double> dis2_;//for random disorder

    //mt19937_64 mt_rand(Parameters_.RandomSeed);


};



double MFParams::random1(){

    return dis1_(Generator1_);

}

double MFParams::random2(){

    return dis2_(Generator2_);

}


void MFParams::initialize(){

    lx_=Coordinates_.lx_;
    ly_=Coordinates_.ly_;
    no_dof_=Coordinates_.no_dof_;
    ns_=Coordinates_.ns_;

    // srand(Parameters_.RandomSeed);

    Disorder.resize(lx_,ly_);


    OParams.resize(ns_);
    for(int site=0;site<ns_;site++){
        OParams[site].resize(6);
        for(int state=0;state<6;state++){
            OParams[site][state].resize(6);
        }
    }


    ofstream Disorder_conf_file("Disorder_conf_used");
    Disorder_conf_file<<"#seed="<<Parameters_.RandomDisorderSeed<<
                        " for mt19937_64 Generator is used"<<endl;
    Disorder_conf_file<<"#ix   iy    Dis[ix,iy]"<<endl;

    ofstream Initial_OrderParams_file("Initial_OrderParams_values_generated.txt");

    if(!Parameters_.Read_OPs){
        for(int site=0;site<ns_;site++){

            for(int i=0;i<6;i++){
                for(int j=i;j<6;j++){

                    if(i!=j){
                        OParams[site][i][j].real(random1());
                        OParams[site][i][j].imag(random1());
                    }
                    else{
                        OParams[site][i][j].real(random1());
                        OParams[site][i][j].imag(0.0);
                    }
                }
            }

            for(int i=0;i<6;i++){
                for(int j=i+1;j<6;j++){
                    OParams[site][j][i]=conj(OParams[site][i][j]);
                }
            }



        }
        Initial_OrderParams_file<<"#seed="<<Parameters_.RandomSeed<<
                                  " for mt19937_64 Generator is used"<<endl;
    }
    else{
        vector<string> OPstring;
        OPstring.clear();
        OPstring.push_back("All_OP");

        for(int op_no=0;op_no<OPstring.size();op_no++){
            string fl_initial_OP_in = Parameters_.File_OPs_in;
            ifstream file_initial_OP_in(fl_initial_OP_in.c_str());
            string temp1,temp2,temp3,temp4,temp5;
            int site_temp, x,y;
            double val_real, val_imag;

            file_initial_OP_in>>temp1>>temp2>>temp3>>temp4>>temp5;
            cout<<temp1<<"  "<<temp2<<" "<<temp3<<"   "<<temp4<<"   "<<temp5<<endl;

           // file_initial_OP_in>>temp1>>temp2>>temp3;
           // cout<<temp1<<"  "<<temp2<<" "<<temp3<<endl;

            for(int site=0;site<ns_;site++){

                for(int i=0;i<6;i++){
                    for(int j=0;j<6;j++){
                        file_initial_OP_in>>site_temp>>x>>y>>val_real>>val_imag;
                        //cout<<ix<<"\t"<<iy<<"\t"<<x<<"\t"<<y<<"\t"<<val<<endl;
                        assert(site==site_temp);
                        assert(x==i);
                        assert(y==j);
                        OParams[site][i][j].real(val_real);
                        OParams[site][i][j].imag(val_imag);
                    }
                }

            }

        }

        Initial_OrderParams_file<<"#OParams are read from "<<Parameters_.File_OPs_in<<" file"<<endl;

    }



    Initial_OrderParams_file<<"#site     state_jm=i      state_jm=j     <a^{dagger}_{i}a_{j}>"<<endl;

    for(int site=0;site<ns_;site++){
        for(int i=0;i<6;i++){
            for(int j=i;j<6;j++){
                Initial_OrderParams_file<<site<<setw(15)<<i<<setw(15)<<j<<setw(15)<<
                                          OParams[site][i][j].real()<<setw(15)<<OParams[site][i][j].imag()<<endl;
            }
            Initial_OrderParams_file<<endl;
        }
    }



    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            //RANDOM Disorder
            Disorder(i,j)=Parameters_.Disorder_Strength*((2.0*random2())-1.0);
            Disorder_conf_file<<i<<"  "<<j<<"  "<<Disorder(i,j)<<endl;
        }
        Disorder_conf_file<<endl;
    }


} // ----------

#endif
