#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>

#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

const double eps=  std::numeric_limits<double>::epsilon();

typedef struct {
    bool AIDrag;
    bool TIDrag;
    bool TipLoss;
    bool HubLoss;
    bool TanInd;
    int IndToler;
    double acorr;
} ccBlade_option_t;

typedef struct {
    // inputs
    double lambda_r;
    double pitch;
    double twist;
    double sigma_p;
    int af_n;
    const double *af_alpha;
    const double *af_cl;
    const double *af_cd;
    double r;
    double r_tip;
    double r_hub;
    double B;
    ccBlade_option_t options;

    // return values
    double a;
    double ap;
    double cl;
    double cd;
} ccBlade_param_t;

typedef struct {
    int n;
    const double *alpha;
    const double *cl;
    const double *cd;
} ccAirFoil_t;

typedef struct {
    int n;
    double B;
    double rho;
    const double *R;
    const double *chord;
    const double *twist;
    int *af_idx;
    ccAirFoil_t *af;
    ccBlade_option_t options;
} ccBlade_t;

typedef struct {
    int32_t *converge;
    double *phi;
    double *a;
    double *ap;
    double *cl;
    double *cd;

    double *v_res;
    double *Fl;
    double *Fd;
    double *Max;
    double *Mtan;
    double *Fax;
    double *Ftan;
    double *Fflap;
    double *Fedge;
    double *cp_i;
    double *ct_i;
    double *cf_i;
    double *ce_i;
    double *cp;
    double *ct;
} ccResult_t;

bool tryGetData(double &value, const char *name, const mxArray *mxData);
int tryGetDataArray(const double **value, const char *name, const mxArray *mxData, int idx, int n);
int tryGetDataIntArray(int **value, const char *name, const mxArray *mxData, int n);
int tryGetAirFoil(ccAirFoil_t **value, const mxArray *mxData);

void CCBlade(double lambda, double pitch, ccBlade_t &ccBlade, ccResult_t &result);
void aeroForces(double lambda, double pitch, double v_wind, ccBlade_t &ccBlade, ccResult_t &result);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs!=4) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of arguments. Expecting (data, lambda, pitch, v_wind)"); return; }
    if(nlhs!=1) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of return values. Expecting [result]"); return; }
    
    if(!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1])!=1) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input lambda must be scalar", 1); return; }
    if(!mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input pitch must be scalar", 1); return; }
    if(!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input v_wind must be scalar", 1); return; }
    
    const mxArray *mxData= prhs[0];
    if(!mxIsStruct(mxData) || mxGetNumberOfElements(mxData)!=1) {
        mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input data must be a scalar struct.\n");
        return;
    }
    
    double lambda= mxGetPr(prhs[1])[0];
    double pitch= mxGetPr(prhs[2])[0];
    double v_wind= mxGetPr(prhs[3])[0];

    ccBlade_t ccBlade;

    double v;
    tryGetData(v, "AIDrag", mxData);
    ccBlade.options.AIDrag= v;
    tryGetData(v, "TIDrag", mxData);
    ccBlade.options.TIDrag= v;
    tryGetData(v, "TipLoss", mxData);
    ccBlade.options.TipLoss= v;
    tryGetData(v, "HubLoss", mxData);
    ccBlade.options.HubLoss= v;
    tryGetData(v, "TanInd", mxData);
    ccBlade.options.TanInd= v;
    tryGetData(v, "IndToler", mxData);
    ccBlade.options.IndToler= v;
    tryGetData(v, "acorr", mxData);
    ccBlade.options.acorr= v;
    tryGetData(v, "B", mxData);
    ccBlade.B= v;
    tryGetData(v, "rho", mxData);
    ccBlade.rho= v;

    int n= tryGetDataArray(&ccBlade.R, "R", mxData, 0, -1);
    ccBlade.n= n;
    tryGetDataArray(&ccBlade.chord, "chord", mxData, 0, n);
    tryGetDataArray(&ccBlade.twist, "twist", mxData, 0, n);
    tryGetDataIntArray(&ccBlade.af_idx, "airfoil_idx", mxData, n);
    int n_af= tryGetAirFoil(&ccBlade.af, mxData);

    for(int i= 0; i<n; ++i)
        if(ccBlade.af_idx[i]>=n_af) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "AirFoil index of node %d is greater than number of air foils.\n", i);
            return;
        }
        
    ccResult_t result;
    
    const char *field_names[]= {"converge", "phi", "a", "ap", "cl", "cd", "v_res", "Fl", "Fd", "Max", "Mtan", "Fax", "Ftan", "Fflap", "Fedge", "cp_i", "ct_i", "cp", "ct", "cf_i", "ce_i"};
    plhs[0]= mxCreateStructMatrix(1, 1, 21, field_names);
    mxArray* field_value;

    field_value= mxCreateNumericMatrix(n, 1, mxINT32_CLASS, mxREAL);
    mxSetField(plhs[0], 0, "converge", field_value);
    result.converge= (int32_T *)mxGetData(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "phi", field_value);
    result.phi= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "a", field_value);
    result.a= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "ap", field_value);
    result.ap= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "cl", field_value);
    result.cl= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "cd", field_value);
    result.cd= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "v_res", field_value);
    result.v_res= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "Fl", field_value);
    result.Fl= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "Fd", field_value);
    result.Fd= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "Max", field_value);
    result.Max= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "Mtan", field_value);
    result.Mtan= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "Fax", field_value);
    result.Fax= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "Ftan", field_value);
    result.Ftan= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "Fflap", field_value);
    result.Fflap= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "Fedge", field_value);
    result.Fedge= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "cp_i", field_value);
    result.cp_i= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "ct_i", field_value);
    result.ct_i= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "cf_i", field_value);
    result.cf_i= mxGetPr(field_value);
    
    field_value= mxCreateDoubleMatrix(n, 1, mxREAL);
    mxSetField(plhs[0], 0, "ce_i", field_value);
    result.ce_i= mxGetPr(field_value);
    
    field_value= mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "cp", field_value);
    result.cp= mxGetPr(field_value);

    field_value= mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "ct", field_value);
    result.ct= mxGetPr(field_value);

    CCBlade(lambda, pitch, ccBlade, result);
    aeroForces(lambda, pitch, v_wind, ccBlade, result);

    mxFree(ccBlade.af_idx);
    mxFree(ccBlade.af);
}

bool tryGetData(double &value, const char *name, const mxArray *mxData) {
    const mxArray *mxField;
    if((mxField= mxGetField(mxData, 0, name))!=NULL) {
        int m_= mxGetM(mxField);
        int n_= mxGetN(mxField);
        if(mxIsSparse(mxField) || !mxIsDouble(mxField) || (m_!=1 && n_!=1)) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Data field '%s' must be a scalar.\n", name);
            return false;
        }
        value= mxGetPr(mxField)[0];
        
        return true;
    } else {
        mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Data field '%s' not found.\n", name);
        return false;
    }
}

int tryGetDataArray(const double **value, const char *name, const mxArray *mxData, int idx, int n) {
    const mxArray *mxField;
    if((mxField= mxGetField(mxData, idx, name))!=NULL) {
        int m_= mxGetM(mxField);
        int n_= mxGetN(mxField);
        if(n<0)
            n= m_*n_;

        if(mxIsSparse(mxField) || !mxIsDouble(mxField) || (m_*n_!=n)) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Data field '%s' must be a vector with %d elements.\n", name, n);
            return -1;
        }
        value[0]= mxGetPr(mxField);
    } else {
        mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Data field '%s' not found.\n", name);
        return -1;
    }

    return n;
}

int tryGetAirFoil(ccAirFoil_t **value, const mxArray *mxData) {
    const mxArray *mxField;
    if((mxField= mxGetField(mxData, 0, "AirFoil"))!=NULL) {
        int n_af= mxGetNumberOfElements(mxField);
        if(mxIsSparse(mxField) || !mxIsStruct(mxField)) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Data field 'AirFoil' must be a struct.\n");
            return -1;
        }
        (*value)= (ccAirFoil_t*)mxMalloc(n_af*sizeof(ccAirFoil_t));

        for(int i= 0; i<n_af; ++i) {
            (*value)[i].n= tryGetDataArray(&((*value)[i].alpha), "alpha", mxField, i, -1);
            tryGetDataArray(&((*value)[i].cl), "cl", mxField, i, (*value)[i].n);
            tryGetDataArray(&((*value)[i].cd), "cd", mxField, i, (*value)[i].n);
        }
        
        return n_af;
    } else {
        mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Data field 'AirFoil' not found.\n");
        return -1;
    }
}

int tryGetDataIntArray(int **value, const char *name, const mxArray *mxData, int n) {
    const mxArray *mxField;
    if((mxField= mxGetField(mxData, 0, name))!=NULL) {
        int m_= mxGetM(mxField);
        int n_= mxGetN(mxField);
        if(n<0)
            n= m_*n_;

        if(mxIsSparse(mxField) || !mxIsDouble(mxField) || (m_*n_!=n)) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Data field '%s' must be a vector with %d elements.\n", name, n);
            return -1;
        }
        *value= (int *)mxMalloc(n*sizeof(int));
        for(int i= 0; i<n; ++i)
            (*value)[i]= mxGetPr(mxField)[i]-1;
    } else {
        mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Data field '%s' not found.\n", name);
        return -1;
    }

    return n;
}

int brent(double(*f)(double, void *), double &xs, double a, double b, void *data, double t= 1e-6, int maxiter= 1000) {

    double fa= f(a, data);
    double fb= f(b, data);

    if(fa*fb>0)
        return -1;

    double c= a;
    double fc= fa;

    double d= b-a;
    double e= d;

    int iter;
    for(iter= 0; iter<maxiter; ++iter) {

        if(fb*fc>0) {
            c= a;
            fc= fa;
            d= b-a;
            e= d;
        }

        if(fabs(fc)<fabs(fb)) {
            a= b;
            b= c;
            c= a;
            fa= fb;
            fb= fc;
            fc= fa;
        }

        double tol= 2*eps*fabs(b)+t;
        double m= (c-b)/2;

        if((fabs(m)>tol) && (fabs(fb)>0)) {
            double s, p, q;

            if((fabs(e)<tol) || (fabs(fa)<=fabs(fb))) {
                d= m;
                e= m;
            } else {
                s= fb/fa;
                if(a==c) {
                    p= 2*m*s;
                    q= 1-s;
                } else {
                    q= fa/fc;
                    double r= fb/fc;
                    p= s*(2*m*q*(q-r)-(b-a)*(r-1));
                    q= (q-1)*(r-1)*(s-1);
                }

                if(p>0)
                    q= -q;
                else
                    p= -p;

                s= e;
                e= d;
                if(( 2*p<3*m*q-abs(tol*q) ) && (p<abs(s*q/2))) {
                    d= p/q;
                } else {
                    d= m;
                    e= m;
                }
            }
            a= b;
            fa= fb;

            if(fabs(d)>tol) {
                b= b+d;
            } else {
                if(m>0)
                    b= b+tol;
                else
                    b= b-tol;
            }
        } else
            break;
            
        fb= f(b, data);
    }
            
    if(iter>=maxiter)
        return -2;
            
    xs= b;
    return iter;
}

int getIndex(const double *xx, int n, double x) {
    int mid, ofs, n_;
    const double *xx_;
    
    if(x<xx[0]) return -1;
    if(x>xx[n-1]) return n-1;
    
    xx_= xx;
    ofs= 0;
    n_= n+1;
    while(n>2 && n<n_) {
        n_= n;
        mid= (n-1)/2;
        if(x>xx_[mid]) {
            xx_+= mid;
            n-= mid;
            ofs+= mid;
        } else {
            n= mid+1;
        }
    }
    
    return ofs;
}

double interp1lin(const double *xx, const double *yy, int n, double x) {
    int iLow;
    double mm;
    
    iLow= getIndex(xx, n, x);
    if(iLow<0) return yy[0];
    if(iLow>=n-1) return yy[n-1];
    
    mm= (x-xx[iLow])/(xx[iLow+1]-xx[iLow]);
    return mm*yy[iLow+1] + (1.0-mm)*yy[iLow];
}

double phiResidual(double phi, void * data_) {
    ccBlade_param_t *data= (ccBlade_param_t *)data_;

    double alpha= phi - data->pitch - data->twist;

    double cl= interp1lin(data->af_alpha, data->af_cl, data->af_n, alpha);
    double cd= interp1lin(data->af_alpha, data->af_cd, data->af_n, alpha);


    double sphi= sin(phi);
    double cphi= cos(phi);

    double cn;
    double ct;

    if(data->options.AIDrag)
        cn= cl*cphi;
    else
        cn= cl*cphi + cd*sphi;

    if(data->options.TIDrag)
        ct= cl*sphi;
    else
        ct= cl*sphi - cd*cphi;
            
    double Ftip= 1.0;
    if(data->options.TipLoss) {
        double factortip= data->B/2.0*(data->r_tip - data->r)/(data->r*fabs(sphi));
        Ftip= 2.0/M_PI*acos(exp(-factortip)) + eps;
    }

    double Fhub= 1.0;
    if(data->options.HubLoss) {
        double factorhub= data->B/2.0*(data->r - data->r_hub)/(data->r*fabs(sphi));
        Fhub= 2.0/M_PI*acos(exp(-factorhub)) + eps;
    }

    double F= Ftip * Fhub;

    double k= data->sigma_p*cn/4.0/F/sphi/sphi;
    double kp= data->sigma_p*ct/4.0/F/sphi/cphi;

    double a;
    double ap;
    
    // compute axial induction factor
    if(phi > 0.0) {  // momentum/empirical
        const double acorr= data->options.acorr;
        // update axial induction factor
        if(k < -acorr/(acorr-1.0)) {  // momentum state
            a= k/(1.0+k);
        } else {
            // Glauert correction
            double k_= 1.0/k;
            double k__= pow((k_*(1.0-2.0*acorr)+2), 2) + 4.0*(k_*acorr*acorr-1.0);
            if(k__<0)
                k__= eps;

            a= 0.5 * ( 2.0 + k_ * (1.0 - 2.0*acorr) - sqrt( k__ ) );
            k= -(a-1.0)/a;
        }
    } else  { // propeller brake region (a and ap not directly used but update anyway)
        if(k > 1.0)
            a= k/(k-1.0);
        else
            a= 0;  // dummy value
    }

    // compute tangential induction factor
    ap= kp/(1.0-kp);
    if(!data->options.TanInd) {
        ap= 0;
        kp= 0;
    }


    data->a= a;
    data->ap= ap;
    data->cl= cl;
    data->cd= cd;

    // error function
    if(phi > 0.0) // momentum/empirical
        return sphi/(1.0+eps-a) - cphi/data->lambda_r*(1.0-kp);
    else  // propeller brake region
        return sphi*(1.0-k) - cphi/data->lambda_r*(1.0-kp);
}

void CCBlade(double lambda, double pitch, ccBlade_t &ccBlade, ccResult_t &result) {
    // ------ BEM solution method see (Ning, doi:10.1002/we.1636) ------
    const double epsilon= 1e-6;
    int n_nodes= ccBlade.n;

    ccBlade_param_t ccParams;
    ccParams.options= ccBlade.options;
    ccParams.B= ccBlade.B;

    ccParams.pitch= pitch;

    ccParams.r_tip= ccBlade.R[n_nodes-1];
    ccParams.r_hub= ccBlade.R[0];

    for(int node= 0; node<n_nodes; ++node) {
        ccParams.sigma_p= ccBlade.B/2.0/M_PI*ccBlade.chord[node]/ccBlade.R[node];
        ccParams.lambda_r= lambda/ccParams.r_tip*ccBlade.R[node];
        ccParams.twist= ccBlade.twist[node];
        
        int af_idx= ccBlade.af_idx[node];
        ccParams.af_n= ccBlade.af[af_idx].n;
        ccParams.af_alpha= ccBlade.af[af_idx].alpha;
        ccParams.af_cl= ccBlade.af[af_idx].cl;
        ccParams.af_cd= ccBlade.af[af_idx].cd;
        ccParams.r= ccBlade.R[node];

        double phi_lower= epsilon;
        double phi_upper= M_PI/2.0;
        

        double res_lower= phiResidual(phi_lower, &ccParams);
        double res_upper= phiResidual(phi_upper, &ccParams);

        if(res_lower*res_upper > 0) {  // an uncommon but possible case
            double res_pi4= phiResidual(-M_PI/4, &ccParams);
            double res_eps= phiResidual(epsilon, &ccParams);
            if(res_pi4 < 0 && res_eps > 0) {
                phi_lower= -M_PI/4.0;
                phi_upper= -epsilon;
            } else {
                phi_lower= M_PI/2.0;
                phi_upper= M_PI - epsilon;
            }
        }  
                
        if(phi_upper==-epsilon && phiResidual(-epsilon, &ccParams)<0) {
            result.phi[node]= -epsilon;
        } else {
            result.converge[node]= brent(phiResidual, result.phi[node], phi_lower, phi_upper, &ccParams, ccBlade.options.IndToler);
        }
//         
        result.a[node]= ccParams.a;
        result.ap[node]= ccParams.ap;
        result.cl[node]= ccParams.cl;
        result.cd[node]= ccParams.cd;
    }
}

double int_torque(double r0, double r1, double a0, double a1) {
    double dR= r1-r0;

    return ((2.0*a1+a0)*dR*dR + (3.0*r0*a1 + 3.0*r0*a0)*dR)/6.0;
}

void aeroForces(double lambda, double pitch, double v_wind, ccBlade_t &ccBlade, ccResult_t &result) {
    int n_nodes= ccBlade.n;
    double Fwind= ccBlade.rho/2.0 * M_PI*pow(ccBlade.R[n_nodes-1], 2.0) * v_wind*v_wind;
    double Pwind= Fwind*v_wind;
    double omega= lambda*v_wind/ccBlade.R[n_nodes-1];

    double fl_last, fd_last, fax_last, ftan_last;
    for(int node= 0; node<n_nodes; ++node) {
        double v_rot= lambda*v_wind * ccBlade.R[node]/ccBlade.R[n_nodes-1];

        double v_res= sqrt(pow(v_wind*(1.0-result.a[node]), 2.0) + pow(v_rot*(1.0+result.ap[node]), 2.0));
        result.v_res[node]= v_res;

        double sphi= sin(result.phi[node]);
        double cphi= cos(result.phi[node]);

        double fl= ccBlade.rho/2.0*ccBlade.chord[node] * v_res*v_res * result.cl[node];
        double fd= ccBlade.rho/2.0*ccBlade.chord[node] * v_res*v_res * result.cd[node];
        double fax= fl*cphi + fd*sphi;
        double ftan= -fd*cphi + fl*sphi;
        
        if(node>0) {
            double dR= ccBlade.R[node]-ccBlade.R[node-1];

            double Fl_sect= (fl+fl_last)*0.5*dR;
            double Fd_sect= (fd+fd_last)*0.5*dR;
            double Max_sect= int_torque(ccBlade.R[node-1], ccBlade.R[node], fax_last, fax);
            double Mtan_sect= int_torque(ccBlade.R[node-1], ccBlade.R[node], ftan_last, ftan);

            result.Fl[node]= 0.5*Fl_sect;
            result.Fd[node]= 0.5*Fd_sect;
            result.Max[node]=  0.5*Max_sect;
            result.Mtan[node]=  0.5*Mtan_sect;
            if(node>1) {
                result.Fl[node-1]+= 0.5*Fl_sect;
                result.Fd[node-1]+= 0.5*Fd_sect;
                result.Max[node-1]+= 0.5*Max_sect;
                result.Mtan[node-1]+= 0.5*Mtan_sect;
            } else {
                result.Fl[node-1]= 0.5*Fl_sect;
                result.Fd[node-1]= 0.5*Fd_sect;
                result.Max[node-1]= 0.5*Max_sect;
                result.Mtan[node-1]= 0.5*Mtan_sect;
            }
        }

        fl_last= fl;
        fd_last= fd;
        fax_last= fax;
        ftan_last= ftan;
    }
    
    double Mtan_sum= 0.0;
    double Fax_sum= 0.0;
    for(int node= 0; node<n_nodes; ++node) {
        double sphi= sin(result.phi[node]);
        double cphi= cos(result.phi[node]);

        double Fax= result.Fl[node]*cphi + result.Fd[node]*sphi;
        double Ftan= -result.Fd[node]*cphi + result.Fl[node]*sphi;

        result.Fax[node]= Fax;
        result.Ftan[node]= Ftan;

        result.Fflap[node]= Fax*cos(pitch) + Ftan*sin(pitch);
        result.Fedge[node]= -Fax*sin(pitch) + Ftan*cos(pitch);

        Mtan_sum+= result.Mtan[node];
        Fax_sum+= result.Fax[node];

        double ProtVec= 3.0*result.Mtan[node]*omega;
        result.cp_i[node]= ProtVec / Pwind;
        result.ct_i[node]= 3.0*result.Fax[node] / Fwind;
        
        result.cf_i[node]= result.Fflap[node] / Fwind;
        result.ce_i[node]= result.Fedge[node] / Fwind;
    }   
    
    double Prot= 3.0*Mtan_sum*omega;

    result.cp[0]= Prot / Pwind;
    result.ct[0]= 3.0*Fax_sum / Fwind;
}
