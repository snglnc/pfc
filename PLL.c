

/**
 * main.c
 */
#include <stdio.h>
#include <math.h>
// to generate more cycle, increase cycle number
// then duplicate the input based on cycle number
// for example : in_a[] = {Va[0], .... , Va[DATA_LENGTH - 1],
//                         Va[0], .... , Va[DATA_LENGTH - 1],
//                         . . .};
#define CYCLE 50
#define DATA_LENGTH 20

float Va[] = {
    156.63, 246.59, 294.72, 305.51, 300.66,
    268.03, 204.18, 125.41, 42.954, -48.322,
    -154.08, -243.95, -293.12, -303.09, -297.98,
    -264.13, -202.1, -122.25, -39.893, 51.818
};
float Vb[] = {
    -308.4, -280.19, -240.66, -186.6, -99.744,
    -0.54547, 92.853, 181.46, 262.05, 312.39,
    311.44, 283.76, 245.04, 188.62, 102.16,
    2.9662, -89.395, -176.17, -259.16, -309.96
};
float Vc[] = {
    156.11, 82.694, -21.783, -128.37, -213.06,
    -269.49, -309.58, -313.4, -273.73, -214.81,
    -154.29, -79.64, 24.679, 132.16, 216.63,
    274.14, 311.11, 315.76, 276.27, 216.22
};





typedef struct _DDATA{
    float *in_a;
    float *in_b;
    float *in_c;
    float *F_est;
    float *Theta_est;
    float *Harmonics;
    float Ts;
    float Kc1;    // Kc are controller gains
    float Kc2;    // choose your controller and
    float Kc3;    // gains accordingly to get satisfied result

}DDATA;

DDATA ddata = {
    .in_a = Va,
    .in_b = Vb,
    .in_c = Vc,
    .Ts = 0.001,
};



unsigned int Data_index;
float Angel_sine=0,Angel_cosine=0;
float vGrid_a_pu,vGrid_b_pu,vGrid_c_pu;
float pu_coff ;
//
// Typedefs
//

//! \brief          Defines the ABC_DQ0_POS transform structure
//!
//! \details        This block stores variables used to transforms
//!                 ABC domain variables to dq domain
//!
typedef struct{
    float alpha;   //!< Output: Alpha component (abc-> alpha beta)
    float beta;    //!< Output: Beta component (abc-> alpha beta)
    float d;       //!< Output: D axis component (alpha beta -> d,q,z)
    float q;       //!< Output: Q axis component (alpha beta -> d,q,z)
    float z;       //!< Output: Z axis component (alpha beta -> d,q,z)
}ABC_DQ0_POS;




//! \brief       Resets internal data to zero
//! \param *v    The ABC_DQ0_POS structure pointer
//! \return None
//!
static inline void ABC_DQ0_POS_reset(ABC_DQ0_POS *v)
{
    v->alpha = 0;
    v->beta = 0;
    v->d = 0;
    v->q = 0;
    v->z = 0;
}

//! \brief             Run ABC_DQ0_POS routine
//! \param *v          The ABC_DQ0_POS structure pointer
//! \param a           Phase a value
//! \param b           Phase b value
//! \param c           Phase c value
//! \param sine_val    sine value of the grid angle
//! \param cosine_val  cosine value of the grid angle
//! \return None
//!
static inline void ABC_DQ0_POS_run(ABC_DQ0_POS *v,
                                   float a, float b, float c,
                                   float sine_val, float cosine_val)
{
    v->alpha = (0.66666666677f) * (a - 0.5f * (b + c));
    v->beta  = (0.57735026913f) * (b - c);
    v->z     = (0.57735026913f) * (a + b + c);
    v->d     =  v->alpha * cosine_val + v->beta * sine_val;
    v->q     = -v->alpha * sine_val   + v->beta * cosine_val;
}

//! \brief          Defines the ABC_DQ0_NEG transform structure
//!
//! \details        This block stores variables used to transforms
//!                 ABC domain variables to dq domain
//!
typedef struct{
    float alpha;   //!< Output: Alpha component (abc-> alpha beta)
    float beta;    //!< Output: Beta component (abc-> alpha beta)
    float d;       //!< Output: D axis component (alpha beta -> d,q,z)
    float q;       //!< Output: Q axis component (alpha beta -> d,q,z)
    float z;       //!< Output: Z axis component (alpha beta -> d,q,z)
} ABC_DQ0_NEG;

//! \brief       Resets internal data to zero
//! \param *v    The ABC_DQ0_NEG structure pointer
//! \return None
//!
static inline void ABC_DQ0_NEG_reset(ABC_DQ0_NEG *v)
{
    v->alpha = 0;
    v->beta = 0;
    v->d = 0;
    v->q = 0;
    v->z = 0;
}

//! \brief             Runs ABC_DQ0_NEG routine
//! \param *v          The ABC_DQ0_NEG structure pointer
//! \param a           Phase a value
//! \param b           Phase b value
//! \param c           Phase c value
//! \param sine_val    sine value of the grid angle
//! \param cosine_val  cosine value of the grid angle
//! \return None
//!
static inline void ABC_DQ0_NEG_run(ABC_DQ0_NEG *v,
                                   float a, float b, float c,
                                   float sine_val, float cosine_val)
{
    v->alpha = (0.66666666677f) * (a - 0.5f * (b + c));
    v->beta  = (0.57735026913f) * (b - c);
    v->z     = (0.57735026913f) * (a + b + c);
    v->d     = v->alpha * cosine_val - v->beta * sine_val;
    v->q     = v->alpha * sine_val   + v->beta * cosine_val;
}

ABC_DQ0_POS vGrid_dq0_pos;
ABC_DQ0_NEG vGrid_dq0_neg;


//! \brief          Defines the coefficients for a loop filter
//!
//! \details        Loop filter coefficients
//!
typedef struct{
    float b1;
    float b0;
} SPLL_3PH_DDSRF_LPF_COEFF;

//! \brief          Defines the SPLL_3PH_DDSRF structure
//!
//! \details        This software module implements a software phase lock loop
//!                 based on decoupled double synchronous reference frame for
//!                 grid connection to three phase grid
//! \return None
//!
typedef struct{
    float d_p_decoupl;  //!< Positive Rotating reference Frame D-axis value
    float d_n_decoupl;  //!< Negative Rotating reference Frame D-axis value
    float q_p_decoupl;  //!< Positive Rotating reference Frame Q-axis value
    float q_n_decoupl;  //!< Negative Rotating reference Frame Q-axis value

    float cos_2theta;   //!< Cos of twice the grid frequency angle
    float sin_2theta;   //!< Sin of twice the grid frequency angle

    float y[2];         //!< Used to store history for filtering the decoupled D and Q axis components
    float x[2];         //!< Used to store history for filtering the decoupled D and Q axis components
    float w[2];         //!< Used to store history for filtering the decoupled D and Q axis components
    float z[2];         //!< Used to store history for filtering the decoupled D and Q axis components
    float k1;           //!< Lpf coefficient
    float k2;           //!< Lpf coefficient
    float d_p_decoupl_lpf;  //!< Decoupled positive sequence D-axis component filtered
    float d_n_decoupl_lpf;  //!< Decoupled negative sequence D-axis component filtered
    float q_p_decoupl_lpf;  //!< Decoupled positive sequence Q-axis component filtered
    float q_n_decoupl_lpf;  //!< Decoupled negative sequence Q-axis component filtered

    float v_q[2];
    float theta[2];     //!< Grid phase angle
    float ylf[2];       //!< Internal Data Buffer for Loop Filter output
    float fo;           //!< Instantaneous Grid Frequency in Hz
    float fn;           //!< Nominal Grid Frequency in Hz
    float delta_t;      //!< 1/Frequency of calling the PLL routine
    SPLL_3PH_DDSRF_LPF_COEFF lpf_coeff;
} SPLL_3PH_DDSRF;


SPLL_3PH_DDSRF spll;

//! \brief              Initialize SPLL_3PH_DDSRF module
//! \param grid_freq    The grid frequency
//! \param delta_t      1/Frequency of calling the PLL routine
//! \param k1           parameter
//! \param k2           parameter
//! \param *spll_obj    The SPLL_3PH_DDSRF structure
//! \return None
//!
static inline void SPLL_3PH_DDSRF_init(float grid_freq, float delta_t,
                                       float k1, float k2,
                                       SPLL_3PH_DDSRF *spll_obj)
{
    spll_obj->d_p_decoupl = (float)(0.0);
    spll_obj->d_n_decoupl = (float)(0.0);

    spll_obj->q_p_decoupl = (float)(0.0);
    spll_obj->q_n_decoupl = (float)(0.0);

    spll_obj->d_p_decoupl_lpf = (float)(0.0);
    spll_obj->d_n_decoupl_lpf = (float)(0.0);

    spll_obj->q_p_decoupl_lpf = (float)(0.0);
    spll_obj->q_n_decoupl_lpf = (float)(0.0);

    spll_obj->y[0] = (float)(0.0);
    spll_obj->y[1] = (float)(0.0);

    spll_obj->x[0] = (float)(0.0);
    spll_obj->x[1] = (float)(0.0);

    spll_obj->w[0] = (float)(0.0);
    spll_obj->w[1] = (float)(0.0);

    spll_obj->z[0] = (float)(0.0);
    spll_obj->z[1] = (float)(0.0);

    spll_obj->k1 = k1;
    spll_obj->k2 = k2;

    spll_obj->v_q[0] = (float)(0.0);
    spll_obj->v_q[1] = (float)(0.0);

    spll_obj->ylf[0] = (float)(0.0);
    spll_obj->ylf[1] = (float)(0.0);

    spll_obj->fo = (float)(0.0);
    spll_obj->fn = (float)(grid_freq);

    spll_obj->theta[0] = (float)(0.0);
    spll_obj->theta[1] = (float)(0.0);

    spll_obj->delta_t = delta_t;
}

//
//! \brief              Reset SPLL_3PH_DDSRF module
//! \param *spll_obj    The SPLL_3PH_DDSRF structure
//! \return None
//!
static inline void SPLL_3PH_DDSRF_reset(SPLL_3PH_DDSRF *spll_obj)
{
    spll_obj->d_p_decoupl = (float)(0.0);
    spll_obj->d_n_decoupl = (float)(0.0);

    spll_obj->q_p_decoupl = (float)(0.0);
    spll_obj->q_n_decoupl = (float)(0.0);

    spll_obj->d_p_decoupl_lpf = (float)(0.0);
    spll_obj->d_n_decoupl_lpf = (float)(0.0);

    spll_obj->q_p_decoupl_lpf = (float)(0.0);
    spll_obj->q_n_decoupl_lpf = (float)(0.0);

    spll_obj->y[0] = (float)(0.0);
    spll_obj->y[1] = (float)(0.0);

    spll_obj->x[0] = (float)(0.0);
    spll_obj->x[1] = (float)(0.0);

    spll_obj->w[0] = (float)(0.0);
    spll_obj->w[1] = (float)(0.0);

    spll_obj->z[0] = (float)(0.0);
    spll_obj->z[1] = (float)(0.0);

    spll_obj->v_q[0] = (float)(0.0);
    spll_obj->v_q[1] = (float)(0.0);

    spll_obj->ylf[0] = (float)(0.0);
    spll_obj->ylf[1] = (float)(0.0);

    spll_obj->fo = (float)(0.0);

    spll_obj->theta[0] = (float)(0.0);
    spll_obj->theta[1] = (float)(0.0);
}

//
//! \brief              Run spll_3PH_srf module
//! \param *spll_obj    The spll_3PH_ddsrf structure
//! \param d_p          D Positive seq component of the grid voltage
//! \param d_n          D Negative seq component of the grid voltage
//! \param q_p          Q Positive seq component of the grid voltage
//! \param q_n          Q Negative seq component of the grid voltage
//! \return None
//!
static inline void SPLL_3PH_DDSRF_run(SPLL_3PH_DDSRF *spll_obj,
                                      float d_p, float d_n,
                                      float q_p, float q_n)
{
    //
    // before calling this routine run the ABC_DQ0_Pos & Neg run routines
    // pass updated values for d_p,d_n,q_p,q_n
    // and update the cos_2theta and sin_2theta values with the prev angle
    //

    //
    // Decoupling Network
    //
    spll_obj->d_p_decoupl = d_p
                           - (spll_obj->d_n_decoupl_lpf * spll_obj->cos_2theta)
                           - (spll_obj->q_n_decoupl * spll_obj->sin_2theta);
    spll_obj->q_p_decoupl = q_p
                           + (spll_obj->d_n_decoupl_lpf * spll_obj->sin_2theta)
                           - (spll_obj->q_n_decoupl * spll_obj->cos_2theta);

    spll_obj->d_n_decoupl = d_n
                           - (spll_obj->d_p_decoupl_lpf * spll_obj->cos_2theta)
                           + (spll_obj->q_p_decoupl * spll_obj->sin_2theta);
    spll_obj->q_n_decoupl = q_n
                           - (spll_obj->d_p_decoupl_lpf * spll_obj->sin_2theta)
                           - (spll_obj->q_p_decoupl * spll_obj->cos_2theta);

    //
    // Low pass filter
    //

    spll_obj->y[1] = (spll_obj->d_p_decoupl * spll_obj->k1)
                   - (spll_obj->y[0] * spll_obj->k2);
    spll_obj->d_p_decoupl_lpf = spll_obj->y[1] + spll_obj->y[0];
    spll_obj->y[0] = spll_obj->y[1];

    spll_obj->x[1] = (spll_obj->q_p_decoupl * spll_obj->k1)
                  - (spll_obj->x[0] * spll_obj->k2);
    spll_obj->q_p_decoupl_lpf = spll_obj->x[1] + spll_obj->x[0];
    spll_obj->x[0] = spll_obj->x[1];

    spll_obj->w[1] = (spll_obj->d_n_decoupl * spll_obj->k1)
                  - (spll_obj->w[0] * spll_obj->k2);
    spll_obj->d_n_decoupl_lpf = spll_obj->w[1] + spll_obj->w[0];
    spll_obj->w[0] = spll_obj->w[1];

    spll_obj->z[1] = (spll_obj->q_n_decoupl * spll_obj->k1)
                  - (spll_obj->z[0] * spll_obj->k2);
    spll_obj->q_n_decoupl_lpf = spll_obj->z[1] + spll_obj->z[0];
    spll_obj->z[0] = spll_obj->z[1];

    spll_obj->v_q[0] = spll_obj->q_p_decoupl;

    //
    // Loop Filter
    //
    spll_obj->ylf[0] = spll_obj->ylf[1]
                    + (spll_obj->lpf_coeff.b0 * spll_obj->v_q[0])
                    + (spll_obj->lpf_coeff.b1 * spll_obj->v_q[1]);
    spll_obj->ylf[1] = spll_obj->ylf[0];
    spll_obj->v_q[1] = spll_obj->v_q[0];

    //
    // VCO
    //
    spll_obj->fo = spll_obj->fn + spll_obj->ylf[0];

    spll_obj->theta[0] = spll_obj->theta[1] +
             ((spll_obj->fo * spll_obj->delta_t)
             * (float)(2.0f * 3.1415926f));

    if(spll_obj->theta[0] > (float)(2.0f * 3.1415926f))
    {
        spll_obj->theta[0] = spll_obj->theta[0] -
                  (float)(2.0f * 3.1415926f);
    }


    spll_obj->theta[1] = spll_obj->theta[0];

    spll_obj->cos_2theta = cosf(spll_obj->theta[1] * 2.0f);
    spll_obj->sin_2theta = sinf(spll_obj->theta[1] * 2.0f);
}



//float inv_a[200],grid_a[200],harmonic[200];
//int t_index;

void estimateFrequencyAndTheta(DDATA *d, int dataSize)
{
    // Implementation for estimating frequency and theta

    vGrid_a_pu =  d->in_a[Data_index]*pu_coff;
    vGrid_b_pu =  d->in_b[Data_index]*pu_coff;
    vGrid_c_pu =  d->in_c[Data_index]*pu_coff;

    Angel_sine   = sinf(spll.theta[1]);
    Angel_cosine = cosf(spll.theta[1]);

    ABC_DQ0_POS_run(&vGrid_dq0_pos,
                    vGrid_a_pu,vGrid_b_pu,vGrid_c_pu,
                    Angel_sine,Angel_cosine);

    ABC_DQ0_NEG_run(&vGrid_dq0_neg,
                    vGrid_a_pu,vGrid_b_pu,vGrid_c_pu,
                    Angel_sine,Angel_cosine);

    SPLL_3PH_DDSRF_run(&spll,vGrid_dq0_pos.d,vGrid_dq0_neg.d,vGrid_dq0_pos.q,vGrid_dq0_neg.q);


    d->Theta_est = spll.theta[1];
    d->F_est = spll.fo;


//    if(++t_index>199)
//        t_index = 0;
//
//    inv_a[t_index]  = vGrid_a_pu;
//    grid_a[t_index] = Angel_cosine;
//    harmonic = Angel_cosine*vGrid_dq0_pos.d;

}


void getHarmonicAmplitudes(DDATA *d, int dataSize)
{
    // Implementation for getting harmonic amplitudes
    //. . .





}


int main()
{
  //  int i = 0;

    ABC_DQ0_POS_reset(&vGrid_dq0_pos);
    ABC_DQ0_NEG_reset(&vGrid_dq0_neg);
    SPLL_3PH_DDSRF_reset(&spll);

    SPLL_3PH_DDSRF_init(50.0f,0.001,0.00188141f,-0.99623717f,&spll);
    spll.lpf_coeff.b0 = 3.807f;
    spll.lpf_coeff.b1 = -3.674f;

    pu_coff = 1.0f/350.0f;

    Data_index = 0;

 //   for(i = 0; i < DATA_LENGTH * CYCLE; i++)
    for(;;)
    {
        estimateFrequencyAndTheta(&ddata, DATA_LENGTH * CYCLE);
        getHarmonicAmplitudes(&ddata, DATA_LENGTH * CYCLE);

        if(++Data_index >=20)
            Data_index=0;
    }
  //  return 0;
}
















