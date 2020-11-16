#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects

#define one_minus_oneOverE 0.6321205588285576784044762298 //1-(1/e)
#define euler 2.7182818284590452353602874713526
#define LP 0
#define HP 1
#define M 2
#define LOPASS 0
#define HIPASS 1

// struct to represent the object's state
typedef struct _multiComp {
	t_pxobject		ob;			// the object itself (t_pxobject in MSP instead of t_object)
    double          threshold[3];
    double          ratio[3];
    double          knee[3];
    double          attack[3];
    double          release[3];
    double          bandMakeup[3];
    double          makeup;
    long            bypass; //toggle
    double          mix;
    long            auto_makeup; //on/off
    double          samplerate;
    
    //user frequency controles for filter
    double          lp_freq;
    double          hp_freq;
    
    
    //coefficients for attack and release filter coefficients
    double          aA[3]; //attack filter coefficient
    double          aR[3];  //release filter coefficient
    
    //storage band inputs
    double          comp_in[3];
    double          factor[3]; //gain factors after compression
    
    //storage for gain smoothing
    double          y1_prev[3]; //y[n-1] stored for the smoothing filter
    double          yL_prev[3]; //level storage
    
    //coefficients for LP and HP filters
    double          lp_a0;
    double          lp_a1;
    double          lp_a2;
    double          lp_b1;
    double          lp_b2;
    
    double          lp2_a0;
    double          lp2_a1;
    double          lp2_a2;
    double          lp2_b1;
    double          lp2_b2;
    
    double          hp_a0;
    double          hp_a1;
    double          hp_a2;
    double          hp_b1;
    double          hp_b2;
    
    double          hp2_a0;
    double          hp2_a1;
    double          hp2_a2;
    double          hp2_b1;
    double          hp2_b2;
    
    //storage for LP and HP filters
    double          lo_x[3];
    double          lo_y[3];
    double          hi_x[3];
    double          hi_y[3];
    double          mid_y[3];
    
    //additional storage for 4-pole filters
    double          lo2_x[3];
    double          lo2_y[3];
    double          hi2_x[3];
    double          hi2_y[3];
    
    //bands on/off
    long lo_on;
    long mid_on;
    long hi_on;
} t_multiComp;


// method prototypes
void *multiComp_new(t_symbol *s, long argc, t_atom *argv);
void multiComp_free(t_multiComp *x);
void multiComp_dsp64(t_multiComp *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void multiComp_perform64(t_multiComp *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void multiComp_perform64s(t_multiComp *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void multiComp_anything(t_multiComp *x, t_symbol *s, long argc, t_atom *argv);
void multiComp_assist(t_multiComp *x, void *b, long io, long index, char *s);

// global class pointer variable
static t_class *multiComp_class = NULL;


//***********************************************************************************************

//helper functions.

//convert amplitude to decibel.
//how to prevent -inf without hardcoding to -96 (=16bit audio) or -144 (=24bit audio) ??
double _atodb(double a) {
    return a == 0 ? -96 : 20 * log10(a);
}

//convert decibel to amplitude
double _dbtoa(double a) {
    return pow(10, a/20);
}

//for filters: returns the center frequency represented as fraction of the sample rate.
double _freqtolin (double frequency, t_multiComp * x) {
    return frequency / x->samplerate;
}

//returns single-pole filter coefficient for given t (attack time or release time)
double _coeff(double t, t_multiComp *x) {
    double smpls = (t / 1000) * x->samplerate;
    return pow((1-one_minus_oneOverE), (1.0/smpls));
}

//calculates filter coefficients for Chebyshev LP or HP filters with up to 20 poles.
//at the end of the function, coefficients for one single filter are stored in a[] and b[] arrays.
//Note that due to round-off-noise, it might yield better results to implement filters with 4+ poles
//as a cascade of 2-pole-filters instead of one single filter. The coefficients for the corresponding
//2-pole filters can be extracted at the end of the subroutine in a[0] a[1] a[2] and b[1] b[2].
//for more information, do some research on dspguide.com (which is a good idea anyway, always)
void _coeff_chebishev(double freq, int hplp, t_multiComp *xx) {
    double a[22];
    double b[22];
    double ta[22];
    double tb[22];
    double fc = freq; //frequency center
    int type = hplp; //0: lopass, 1: hipass
    double pr = 0.5; //passband ripple in %
    int np = 2; //number of poles
    for(int i = 0; i<22; i++) {
        a[i] = 0; b[i] = 0;
    }
    a[2] = 1; b[2] = 1;
    
    //subroutine: called once per 2-pole-stage.
    for(int p = 1; p <= np/2; p++) {
        double rp = -cos(PI / (np*2) + (p-1) * (PI / np));
        double ip = sin(PI / (np*2) + (p-1) * (PI / np));
        
        if (pr != 0) {
            double es = pow(
                pow((100 / (100-pr)), 2) - 1,
                            0.5);
            double vx = (1.0/np) * log(
                (1.0/es) + pow(1.0 / pow(es, 2) + 1, 0.5)
            );
            double kx = (1.0/np) * log(
                (1.0/es) + pow(1.0 / pow(es, 2) - 1, 0.5)
            );
            kx = (pow(euler, kx) + pow(euler, -kx)) / 2;
            rp = rp * ( (pow(euler, vx) - pow(euler, -vx)) / 2) / kx;
            ip = ip * ( (pow(euler, vx) + pow(euler, -vx)) / 2) / kx;
        }
        double t = 2 * tan(0.5);
        double w = 2 * PI * fc;
        double m = pow(rp, 2) + pow(ip, 2);
        double d = 4 - 4 * rp * t + m * pow(t, 2);
        double x[3];
        double y[3];
        x[0] = pow(t, 2) / d;
        x[1] = 2 * pow(t, 2) / d;
        x[2] = pow(t, 2) / d;
        y[1] =  (8 - 2*m*pow(t, 2)) / d;
        y[2] = (-4 - 4*rp*t - m*pow(t, 2)) /d;
        double k;
        if(type == HIPASS) k = -cos(w/2 + 0.5) / cos(w/2 - 0.5);
        else k = sin(0.5 - w/2.0) / sin(0.5 + w/2.0);
        d = 1 + y[1]*k - y[2]*pow(k, 2);
        a[0] = (x[0] - x[1]*k + x[2] * pow(k, 2)) / d;
        a[1] = (-2*x[0]*k + x[1] + x[1]*pow(k, 2) - 2*x[2]*k) / d;
        a[2] = (x[0]*pow(k, 2) - x[1]*k +x[2]) / d;
        b[1] = (2*k + y[1] + y[1]*pow(k, 2) - 2*y[2]*k) / d;
        b[2] = (-pow(k, 2) - y[1] * k + y[2]) / d;
        if(type == HIPASS) {
            a[1] = -a[1];
            b[1] = -b[1];
            if(p == 1) {
                xx->hp_a0 = a[0];
                xx->hp_a1 = a[1];
                xx->hp_a2 = a[2];
                xx->hp_b1 = b[1];
                xx->hp_b2 = b[2];
            }
            else if(p ==2) {
                xx->hp2_a0 = a[0];
                xx->hp2_a1 = a[1];
                xx->hp2_a2 = a[2];
                xx->hp2_b1 = b[1];
                xx->hp2_b2 = b[2];
            }
        }
        else {
            if(p == 1) {
                xx->lp_a0 = a[0];
                xx->lp_a1 = a[1];
                xx->lp_a2 = a[2];
                xx->lp_b1 = b[1];
                xx->lp_b2 = b[2];
            }
            else if(p == 2) {
                xx->lp2_a0 = a[0];
                xx->lp2_a1 = a[1];
                xx->lp2_a2 = a[2];
                xx->lp2_b1 = b[1];
                xx->lp2_b2 = b[2];
            }
        }
    }
    //subroutine end
    
    for(int i = 0; i<22; i++) {
        ta[i] = a[i];
        tb[i] = b[i];
    }
    for(int i = 2; i<22; i++) {
        a[i] = a[0] * ta[i] + a[1] * ta[i-1] + a[2] * ta[i-2];
        b[i] = tb[i] - b[1] * tb[i-1] - b[2] * tb[i-2];
    }
    
    b[2] = 0;
    for(int i = 0; i<20; i++) {
        a[i] = a[i+2];
        b[i] = -b[i+2];
    }
    
    double sa = 0.0;
    double sb = 0.0;
    for(int i = 0; i<20; i++) {
        if(type == LOPASS) {
            sa += a[i];
            sb += b[i];
        }
        else {
            sa += a[i] * pow((-1), i);
            sb += b[i] * pow((-1), i);
        }
    }
    
    double gain = sa / (1.0 - sb);
    
    for(int i = 0; i < 20; i++) {
        a[i] = a[i] / gain;
    }
    /*lp_a0 = a[0];
    lp_a1 = a[1];
    lp_a2 = a[2];
    lp_b1 = b[1];
    lp_b2 = b[2];*/
    
}

double _max(double a, double b) {
    return (a > b ? a : b);
}

void ext_main(void *r)
{
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	t_class *c = class_new("multiComp~", (method)multiComp_new, (method)dsp_free, (long)sizeof(t_multiComp), 0L, A_GIMME, 0);

	class_addmethod(c, (method)multiComp_dsp64,		"dsp64",	A_CANT, 0);
    class_addmethod(c, (method)multiComp_anything,    "anything",    A_GIMME, 0);
    class_addmethod(c, (method)multiComp_assist, "assist", A_CANT, 0);

	class_dspinit(c);
	class_register(CLASS_BOX, c);
	multiComp_class = c;
}


void *multiComp_new(t_symbol *s, long argc, t_atom *argv)
{
	t_multiComp *x = (t_multiComp *)object_alloc(multiComp_class);

	if (x) {
		dsp_setup((t_pxobject *)x, 2);	// MSP inlets: arg is # of inlets and is REQUIRED!
		// use 0 if you don't need inlets
		outlet_new(x, "signal"); 		// signal outlet (note "signal" rather than NULL)
        for(int i = 0; i<3; i++) {
            x->threshold[i] = 1.0;
            x->ratio[i] = 1.0;
            x->knee[i] = 0.0;
            x->attack[i] = 0.0;
            x->release[i] = 0.0;
            x->bandMakeup[i] = 0.0;
            x->y1_prev[i] = 0.0;
            x->yL_prev[i] = 0.0;
            x->aA[i] = 0.01;
            x->aR[i] = 0.05;
            x->lo_x[i] = 0.0;
            x->lo_y[i] = 0.0;
            x->hi_x[i] = 0.0;
            x->hi_y[i] = 0.0;
        }

        x->makeup = 0.0;
        x->bypass = 0;
        x->auto_makeup = 0;
        x->lp_freq = -1;
        x->hp_freq = -1;
        x->lo_on = 1;
        x->mid_on = 1;
        x->hi_on = 1;
        x->mix = 1.0;
	}
	return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void multiComp_free(t_multiComp *x)
{
	;
}

void multiComp_assist(t_multiComp *x, void *b, long io, long index, char *s)
{
    switch (io) {
        case 1:
            switch (index) {
                case 0:
                    strncpy_zero(s, "Inlet 1: the audio signal to compress", 512);
                    break;
                case 1:
                    strncpy_zero(s, "Inlet 2 (optional): the audio signal to use as sidechain", 512);
                    break;
            }
            break;
        case 2:
            strncpy_zero(s, "Outlet 1: the compressed audio signal", 512);
            break;
    }
}

void multiComp_anything(t_multiComp *x, t_symbol *s, long argc, t_atom *argv) {
    char *selector = s->s_name;
    
    /* globals */
    if( strncmp(selector, "makeup", 6) == 0 ) {
        x->makeup = atom_getfloat(argv);
        post("updated makeup gain: %f", x->makeup);
    }
    else  if ( strncmp(selector, "lo_freq", 7) == 0 ) {
        x->lp_freq = atom_getfloat(argv);
        _coeff_chebishev(_freqtolin(atom_getfloat(argv), x), 0, x);
        post("updated lo pass coefficients: %f, %f, %f // %f, %f", x->lp_a0, x->lp_a1, x->lp_a2, x->lp_b1, x->lp_b2);
    }
    else  if ( strncmp(selector, "hi_freq", 7) == 0 ) {
        x->hp_freq = atom_getfloat(argv);
        _coeff_chebishev(_freqtolin(atom_getfloat(argv), x), 1, x);
        post("updated hi pass coefficients: %f, %f, %f // %f, %f", x->hp_a0, x->hp_a1, x->hp_a2, x->hp_b1, x->hp_b2);
    }
    else if ( strncmp(selector, "bypass", 6) == 0 ) {
        int flag = (int)atom_getfloat(argv);
        x->bypass = flag > 0 ? 1 : 0;
        post("updated bypass: %i", x->bypass);
    }
    
    /*low ferquency band*/
    else  if ( strncmp(selector, "lo_threshold", 12) == 0 ) {
        x->threshold[0] = atom_getfloat(argv);
        post("updated lo threshold: %f", x->threshold[0]);
    }
    else  if ( strncmp(selector, "lo_ratio", 8) == 0 ) {
        x->ratio[0] = atom_getfloat(argv);
        post("updated lo ratio: %f", x->ratio[0]);
    }
    else  if ( strncmp(selector, "lo_attack", 9) == 0 ) {
        x->attack[0] = atom_getfloat(argv);
        x->aA[0] = 0.1;
        x->aA[0] = _coeff(x->attack[0], x);
        post("updated lo attack coeff: %f", x->aA[0]);
    }
    else  if ( strncmp(selector, "lo_release", 10) == 0 ) {
        x->release[0] = atom_getfloat(argv);
        x->aR[0] = 0.5;
        x->aR[0] = _coeff(x->release[0], x);
        post("updated lo release coeff: %f", x->aR[0]);
    }
    else  if ( strncmp(selector, "lo_knee", 7) == 0 ) {
        x->knee[0] = atom_getfloat(argv);
        post("updated lo knee: %f dB", x->knee[0]);
    }
    else if ( strncmp(selector, "lo_gain", 7) == 0) {
        x->bandMakeup[0] = atom_getfloat(argv);
        post("updated lo output level: %f dB", x->bandMakeup[0]);
    }
    else if ( strncmp(selector, "lo_on", 7) == 0) {
        int flag = (int)atom_getfloat(argv);
        x->lo_on = flag > 0 ? 1 : 0;
        post("updated lo on/off: %i", x->lo_on);
    }
    
    /*mid ferquency band*/
    else  if ( strncmp(selector, "mid_threshold", 13) == 0 ) {
        x->threshold[1] = atom_getfloat(argv);
        post("updated mid threshold: %f", x->threshold[1]);
    }
    else  if ( strncmp(selector, "mid_ratio", 9) == 0 ) {
        x->ratio[1] = atom_getfloat(argv);
        post("updated mid ratio: %f", x->ratio[1]);
    }
    else  if ( strncmp(selector, "mid_attack", 10) == 0 ) {
        x->attack[1] = atom_getfloat(argv);
        x->aA[1] = 0.1;
        x->aA[1] = _coeff(x->attack[1], x);
        post("updated mid attack coeff: %f", x->aA[1]);
    }
    else  if ( strncmp(selector, "mid_release", 11) == 0 ) {
        x->release[1] = atom_getfloat(argv);
        x->aR[1] = 0.5;
        x->aR[1] = _coeff(x->release[1], x);
        post("updated mid release coeff: %f", x->aR[1]);
    }
    else  if ( strncmp(selector, "mid_knee", 8) == 0 ) {
        x->knee[1] = atom_getfloat(argv);
        post("updated mid knee: %f dB", x->knee[1]);
    }
    else if ( strncmp(selector, "mid_gain", 8) == 0) {
        x->bandMakeup[1] = atom_getfloat(argv);
        post("updated mid output level: %f dB", x->bandMakeup[1]);
    }
    else if ( strncmp(selector, "mid_on", 8) == 0) {
        int flag = (int)atom_getfloat(argv);
        x->mid_on = flag > 0 ? 1 : 0;
        post("updated mid on/off: %i", x->mid_on);
    }
    
    /* high frequency band */
    else  if ( strncmp(selector, "hi_threshold", 12) == 0 ) {
        x->threshold[2] = atom_getfloat(argv);
        post("updated hi threshold: %f", x->threshold[2]);
    }
    else  if ( strncmp(selector, "hi_ratio", 8) == 0 ) {
        x->ratio[2] = atom_getfloat(argv);
        post("updated hi ratio: %f", x->ratio[2]);
    }
    else  if ( strncmp(selector, "hi_attack", 9) == 0 ) {
        x->attack[2] = atom_getfloat(argv);
        x->aA[2] = 0.1;
        x->aA[2] = _coeff(x->attack[2], x);
        post("updated hi attack coeff: %f", x->aA[2]);
    }
    else  if ( strncmp(selector, "hi_release", 10) == 0 ) {
        x->release[2] = atom_getfloat(argv);
        x->aR[2] = 0.5;
        x->aR[2] = _coeff(x->release[2], x);
        post("updated hi release coeff: %f", x->aR[2]);
    }
    else  if ( strncmp(selector, "hi_knee", 7) == 0 ) {
        x->knee[2] = atom_getfloat(argv);
        post("updated hi knee: %f dB", x->knee[2]);
    }
    else if ( strncmp(selector, "hi_gain", 7) == 0) {
        x->bandMakeup[2] = atom_getfloat(argv);
        post("updated hi output level: %f dB", x->bandMakeup[2]);
    }
    else if ( strncmp(selector, "hi_on", 7) == 0) {
        int flag = (int)atom_getfloat(argv);
        x->hi_on = flag > 0 ? 1 : 0;
        post("updated hi on/off: %i", x->hi_on);
    }
    
    else {
        post("unknown parameter. The following parameters are valid:");
        post("GLOBAL:");
        post("lo_freq");
        post("hi_freq");
        post("makeup");
        post("LOWER BAND:");
        post("lo_threshold");
        post("lo_ratio");
        post("lo_attack");
        post("lo_release");
        post("lo_knee");
        post("lo_gain");
        post("lo_on (lo output on/off)");
        post("MID BAND:");
        post("mid_threshold");
        post("mid_ratio");
        post("mid_attack");
        post("mid_release");
        post("mid_knee");
        post("mid_gain");
        post("mid_on (mid output on/off)");
        post("HI BAND:");
        post("hi_threshold");
        post("hi_ratio");
        post("hi_attack");
        post("hi_release");
        post("hi_knee");
        post("hi_gain");
        post("hi_on (hi output on/off)");
    }
}

// registers a function for the signal chain in Max
void multiComp_dsp64(t_multiComp *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
    x->samplerate = samplerate;
	post("my sampling rate is: %f", samplerate);
    if(x->lp_freq < 0) _coeff_chebishev(0.5, 0, x);
    else _coeff_chebishev(_freqtolin(x->lp_freq, x), 0, x);
    if(x->hp_freq < 0) _coeff_chebishev(0.0, 0, x);
    else _coeff_chebishev(_freqtolin(x->hp_freq, x), 0, x);

    //if right inlet is connected, apply side chain compression.
	if(count[1]) object_method(dsp64, gensym("dsp_add64"), x, multiComp_perform64s, 0, NULL);
    else object_method(dsp64, gensym("dsp_add64"), x, multiComp_perform64, 0, NULL);
}

// this is the 64-bit perform method audio vectors
void multiComp_perform64(t_multiComp *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	t_double *inL = ins[0];		// we get audio for each inlet of the object from the **ins argument
	t_double *outL = outs[0];	// we get audio for each outlet of the object from the **outs argument
	int n = sampleframes;
    
    double xG; //current compressor input sample.
    double tmp; //intermediate value for attenuation calculation.
    double yG; //current compressor output sample.
    double xL; //current gain smoothing input sample.
    double y1; //intermediate value for smoothing filter.
    double yL; //current gain smoothing output sample.
    double cdb; //final compressor gain in dB
    
    //filters: linear-phase, not zero-phase. via FFT?
    
    
	//perform calculations
    while (n--) {
        /*band separation routine starts here*/
        /*The band separation uses two Chebyshev IIR filters*/
        //apply lo pass filter to input signal.
        x->lo_x[M] = (*inL);
        x->lo_y[M] = (x->lp_a0 * x->lo_x[M]) + (x->lp_a1 * x->lo_x[M-1]) + (x->lp_a2 * x->lo_x[M-2]) + (x->lp_b1 * x->lo_y[M-1]) + (x->lp_b2 * x->lo_y[M-2]);
        x->lo_x[M-2] = x->lo_x[M-1];
        x->lo_x[M-1] = x->lo_x[M];
        x->lo_y[M-2] = x->lo_y[M-1];
        x->lo_y[M-1] = x->lo_y[M];
        
        //apply hi pass filter to difference between input signal and low-passed signal.
        x->hi_x[M] = x->lo_x[M] - x->lo_y[M];
        x->hi_y[M] = (x->hp_a0 * x->hi_x[M]) + (x->hp_a1 * x->hi_x[M-1]) + (x->hp_a2 * x->hi_x[M-2]) + (x->hp_b1 * x->hi_y[M-1]) + (x->hp_b2 * x->hi_y[M-2]);
        x->hi_x[M-2] = x->hi_x[M-1];
        x->hi_x[M-1] = x->hi_x[M];
        x->hi_y[M-2] = x->hi_y[M-1];
        x->hi_y[M-1] = x->hi_y[M];
        
        //the mid band is just the difference between hi-pass input and hi-pass output (no further filtering needed).
        x->mid_y[M] = x->hi_x[M] - x->hi_y[M];
        
        /* compression routine starts here */
        x->comp_in[0] = x->lo_y[M];
        x->comp_in[1] = x->mid_y[M];
        x->comp_in[2] = x->hi_y[M];
        //compress each band individually
        for(int i = 0; i<3; i++) {
            //convert sample to dB so we can operate in log domain
            xG = _atodb(fabs(x->comp_in[i]));
            tmp = 2*(xG - x->threshold[i]);
            //apply threshold, ratio, knee (eq. 4)
            if(tmp < -(x->knee[i])) yG = xG;
            else if (fabs(tmp) <= x->knee[i]) yG = xG + (1/x->ratio[i]-1) * pow(xG - x->threshold[i] + x->knee[i]/2, 2) / (2*x->knee[i]);
            else if ((tmp > x->knee[i])) yG = x->threshold[i] + (xG - x->threshold[i]) / x->ratio[i];
            else error("error computing yL.");
            //gain smoothing by unipolar IIR filter (eq. 17)
            xL = xG - yG;
            y1 = _max(xL, x->aR[i] * x->y1_prev[i] + (1 - x->aR[i]) * xL);
            yL = x->aA[i] * x->yL_prev[i] + (1 - x->aA[i]) * y1;
            //add makeup gain
            if(x->auto_makeup == 1 && (x->ratio[i] > 1 && x->threshold[i] < 0)) {
                cdb = -yL + x->makeup + fabs(x->threshold[i] / x->ratio[i]);
            }
            else cdb = -yL + x->makeup + x->bandMakeup[i];
            //convert back to linear domain
            x->factor[i] = _dbtoa(cdb);
            //store data for next iteration
            x->y1_prev[i] = y1;
            x->yL_prev[i] = yL;
        }
        
        //add bands back together in output
        *outL = 0.0;
        if(x->lo_on) *outL += x->factor[0] * x->lo_y[M];
        if(x->mid_on) *outL += x->factor[1] * x->mid_y[M];
        if(x->hi_on) *outL += x->factor[2] * x->hi_y[M];
        
        //*outL = (x->mix * (*outL)) + (1 - x->mix * *inL);
        
        //move pointers to next sample
        *outL++;
        *inL++;
    }
}

// this is the 64-bit perform method audio vectors with sidechain compressing enabled
void multiComp_perform64s(t_multiComp *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    t_double *inS;
    inS = numins > 1 ? ins[1] : NULL; //set sidechain input if connected
    t_double *inL = ins[0];        // we get audio for each inlet of the object from the **ins argument
    t_double *outL = outs[0];    // we get audio for each outlet of the object from the **outs argument
    int n = sampleframes;
    
    double xG; //current compressor input sample.
    double tmp; //intermediate value for attenuation calculation.
    double yG; //current compressor output sample.
    double xL; //current gain smoothing input sample.
    double y1; //intermediate value for smoothing filter.
    double yL; //current gain smoothing output sample.
    double cdb; //final compressor gain in dB
    
    //filters: linear-phase, not zero-phase. via FFT?
    
    
    //perform calculations
    while (n--) {
        /*band separation routine starts here*/
        /*The band separation uses two Chebyshev IIR filters*/
        //apply lo pass filter to input signal.
        x->lo_x[M] = (*inL);
        x->lo_y[M] = (x->lp_a0 * x->lo_x[M]) + (x->lp_a1 * x->lo_x[M-1]) + (x->lp_a2 * x->lo_x[M-2]) + (x->lp_b1 * x->lo_y[M-1]) + (x->lp_b2 * x->lo_y[M-2]);
        x->lo_x[M-2] = x->lo_x[M-1];
        x->lo_x[M-1] = x->lo_x[M];
        x->lo_y[M-2] = x->lo_y[M-1];
        x->lo_y[M-1] = x->lo_y[M];
        
        //apply hi pass filter to difference between input signal and low-passed signal.
        x->hi_x[M] = x->lo_x[M] - x->lo_y[M];
        x->hi_y[M] = (x->hp_a0 * x->hi_x[M]) + (x->hp_a1 * x->hi_x[M-1]) + (x->hp_a2 * x->hi_x[M-2]) + (x->hp_b1 * x->hi_y[M-1]) + (x->hp_b2 * x->hi_y[M-2]);
        x->hi_x[M-2] = x->hi_x[M-1];
        x->hi_x[M-1] = x->hi_x[M];
        x->hi_y[M-2] = x->hi_y[M-1];
        x->hi_y[M-1] = x->hi_y[M];
        
        //the mid band is just the difference between hi-pass input and hi-pass output (no further filtering needed).
        x->mid_y[M] = x->hi_x[M] - x->hi_y[M];
        
        /* compression routine starts here */
        x->comp_in[0] = inS ? *inS : x->lo_y[M];
        x->comp_in[1] = inS ? *inS : x->mid_y[M];
        x->comp_in[2] = inS ? *inS : x->hi_y[M];
        //compress each band individually
        for(int i = 0; i<3; i++) {
            //convert sample to dB so we can operate in log domain
            xG = _atodb(fabs(x->comp_in[i]));
            tmp = 2*(xG - x->threshold[i]);
            //apply threshold, ratio, knee (eq. 4)
            if(tmp < -(x->knee[i])) yG = xG;
            else if (fabs(tmp) <= x->knee[i]) yG = xG + (1/x->ratio[i]-1) * pow(xG - x->threshold[i] + x->knee[i]/2, 2) / (2*x->knee[i]);
            else if ((tmp > x->knee[i])) yG = x->threshold[i] + (xG - x->threshold[i]) / x->ratio[i];
            else error("error computing yL.");
            //gain smoothing by unipolar IIR filter (eq. 17)
            xL = xG - yG;
            y1 = _max(xL, x->aR[i] * x->y1_prev[i] + (1 - x->aR[i]) * xL);
            yL = x->aA[i] * x->yL_prev[i] + (1 - x->aA[i]) * y1;
            //add makeup gain
            if(x->auto_makeup == 1 && (x->ratio[i] > 1 && x->threshold[i] < 0)) {
                cdb = -yL + x->makeup + fabs(x->threshold[i] / x->ratio[i]);
            }
            else cdb = -yL + x->makeup + x->bandMakeup[i];
            //convert back to linear domain
            x->factor[i] = _dbtoa(cdb);
            //store data for next iteration
            x->y1_prev[i] = y1;
            x->yL_prev[i] = yL;
        }
        
        //add bands back together in output
        *outL = 0.0;
        if(x->lo_on) *outL += x->factor[0] * x->lo_y[M];
        if(x->mid_on) *outL += x->factor[1] * x->mid_y[M];
        if(x->hi_on) *outL += x->factor[2] * x->hi_y[M];
        
        //*outL = (x->mix * (*outL)) + (1 - x->mix * *inL);
        
        //move pointers to next sample
        *outL++;
        *inL++;
        if(inS) *inS++;
    }
}
