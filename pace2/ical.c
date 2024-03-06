/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__ICaL
#define _nrn_initial _nrn_initial__ICaL
#define nrn_cur _nrn_cur__ICaL
#define _nrn_current _nrn_current__ICaL
#define nrn_jacob _nrn_jacob__ICaL
#define nrn_state _nrn_state__ICaL
#define _net_receive _net_receive__ICaL 
#define rates rates__ICaL 
#define states states__ICaL 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gCaL _p[0]
#define gCaL_columnindex 0
#define ica _p[1]
#define ica_columnindex 1
#define m _p[2]
#define m_columnindex 2
#define n _p[3]
#define n_columnindex 3
#define h _p[4]
#define h_columnindex 4
#define Dm _p[5]
#define Dm_columnindex 5
#define Dn _p[6]
#define Dn_columnindex 6
#define Dh _p[7]
#define Dh_columnindex 7
#define cai _p[8]
#define cai_columnindex 8
#define _g _p[9]
#define _g_columnindex 9
#define _ion_cai	*_ppvar[0]._pval
#define _ion_ica	*_ppvar[1]._pval
#define _ion_dicadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_alp(void);
 static void _hoc_bet(void);
 static void _hoc_ce(void);
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_ICaL", _hoc_setdata,
 "alp_ICaL", _hoc_alp,
 "bet_ICaL", _hoc_bet,
 "ce_ICaL", _hoc_ce,
 "rates_ICaL", _hoc_rates,
 0, 0
};
#define alp alp_ICaL
#define bet bet_ICaL
#define ce ce_ICaL
 extern double alp( double , double );
 extern double bet( double , double );
 extern double ce( double );
 /* declare global and static user variables */
#define hinf hinf_ICaL
 double hinf = 0;
#define mtau mtau_ICaL
 double mtau = 0;
#define minf minf_ICaL
 double minf = 0;
#define ntau ntau_ICaL
 double ntau = 0;
#define ninf ninf_ICaL
 double ninf = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gCaL_ICaL", 0, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "mtau_ICaL", "ms",
 "ntau_ICaL", "ms",
 "gCaL_ICaL", "S/cm2",
 "ica_ICaL", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double n0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "minf_ICaL", &minf_ICaL,
 "ninf_ICaL", &ninf_ICaL,
 "hinf_ICaL", &hinf_ICaL,
 "mtau_ICaL", &mtau_ICaL,
 "ntau_ICaL", &ntau_ICaL,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"ICaL",
 "gCaL_ICaL",
 0,
 "ica_ICaL",
 0,
 "m_ICaL",
 "n_ICaL",
 "h_ICaL",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 10, _prop);
 	/*initialize range parameters*/
 	gCaL = 0.0002476;
 	_prop->param = _p;
 	_prop->param_size = 10;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _ical_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 10, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ICaL ical.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Cardiac L-type Calcium channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double, double);
 static int _deriv1_advance = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist2[3]; static double _dlist2[3];
 static double _savstate1[3], *_temp1 = _savstate1;
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v , cai ) ;
   Dm = ( minf - m ) / mtau ;
   Dn = ( ninf - n ) / ntau ;
   Dh = ( hinf - h ) / 2.0 ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v , cai ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ntau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / 2.0 )) ;
  return 0;
}
 /*END CVODE*/
 
static int states () {_reset=0;
 { static int _recurse = 0;
 int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 {int _id; for(_id=0; _id < 3; _id++) { _savstate1[_id] = _p[_slist1[_id]];}}
 error = newton(3,_slist2, _p, states, _dlist2);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   rates ( _threadargscomma_ v , cai ) ;
   Dm = ( minf - m ) / mtau ;
   Dn = ( ninf - n ) / ntau ;
   Dh = ( hinf - h ) / 2.0 ;
   {int _id; for(_id=0; _id < 3; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _p[_dlist1[_id]] - (_p[_slist1[_id]] - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _p[_slist1[_id]] - _savstate1[_id];}}}
 } }
 return _reset;}
 
double alp (  double _lv , double _li ) {
   double _lalp;
 double _lq10 ;
 _lq10 = pow( 3.0 , ( ( celsius - 37.0 ) / 10.0 ) ) ;
   if ( _li  == 0.0 ) {
     _lalp = ( 1.0 - exp ( - ( _lv + 10.0 ) / 6.24 ) ) / ( 0.035 * ( _lv + 10.0 ) * ( 1.0 + exp ( - ( _lv + 10.0 ) / 6.24 ) ) ) / ( _lq10 * 1.0 ) ;
     }
   else if ( _li  == 1.0 ) {
     _lalp = 9.0 / ( 0.0197 * exp ( - pow( 0.0337 , 2.0 ) * pow( ( _lv + 10.0 ) , 2.0 ) ) + 0.02 ) / ( _lq10 * 1.0 ) ;
     }
   
return _lalp;
 }
 
static void _hoc_alp(void) {
  double _r;
   _r =  alp (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double bet (  double _lv , double _li ) {
   double _lbet;
 _lv = _lv ;
   if ( _li  == 0.0 ) {
     _lbet = 1.0 / ( 1.0 + exp ( - ( _lv + 10.0 ) / 8.0 ) ) ;
     }
   else if ( _li  == 1.0 ) {
     _lbet = 1.0 / ( 1.0 + exp ( ( _lv + 28.0 ) / 6.9 ) ) ;
     }
   
return _lbet;
 }
 
static void _hoc_bet(void) {
  double _r;
   _r =  bet (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double ce (  double _lcai ) {
   double _lce;
 _lce = 1.0 / ( 1.0 + ( _lcai / 0.00035 ) ) ;
   
return _lce;
 }
 
static void _hoc_ce(void) {
  double _r;
   _r =  ce (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv , double _lcai ) {
   double _la , _lb , _lc ;
 _la = alp ( _threadargscomma_ _lv , 0.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 0.0 ) ;
   mtau = _la ;
   minf = _lb ;
   _la = alp ( _threadargscomma_ _lv , 1.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 1.0 ) ;
   ntau = _la ;
   ninf = _lb ;
   _lc = ce ( _threadargscomma_ _lcai ) ;
   hinf = _lc ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
  n = n0;
 {
   rates ( _threadargscomma_ v , cai ) ;
   m = minf ;
   n = ninf ;
   h = hinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  cai = _ion_cai;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ica = gCaL * m * n * h * ( v - 65.0 ) ;
   }
 _current += ica;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  cai = _ion_cai;
 _g = _nrn_current(_v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  cai = _ion_cai;
 { error = _deriv1_advance = 1;
 derivimplicit(_ninits, 3, _slist1, _dlist1, _p, &t, dt, states, &_temp1);
_deriv1_advance = 0;
 if(error){fprintf(stderr,"at line 46 in file ical.mod:\n	SOLVE states METHOD derivimplicit\n"); nrn_complain(_p); abort_run(error);}
    if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = m_columnindex;  _dlist1[0] = Dm_columnindex;
 _slist1[1] = n_columnindex;  _dlist1[1] = Dn_columnindex;
 _slist1[2] = h_columnindex;  _dlist1[2] = Dh_columnindex;
 _slist2[0] = h_columnindex;
 _slist2[1] = m_columnindex;
 _slist2[2] = n_columnindex;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "ical.mod";
static const char* nmodl_file_text = 
  "TITLE Cardiac L-type Calcium channel\n"
  ": Hodgkin - Huxley type calcium channel from Courtemanche et al Am J Physiol 1998 275:H301b with voltage and calcium dependent inactivation\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX ICaL\n"
  "	USEION ca READ cai WRITE ica\n"
  "	RANGE gCaL, ica\n"
  "	GLOBAL minf, ninf, hinf, mtau, ntau \n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "        (mM) = (milli/liter)\n"
  "        \n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gCaL=0.0002476 (S/cm2) <0,1e9> \n"
  "               \n"
  "}\n"
  "\n"
  "STATE {\n"
  "	m n h\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	celsius (degC) : 37\n"
  "        cai (mM)\n"
  "	ica (mA/cm2)\n"
  "	minf ninf hinf\n"
  "	mtau (ms)\n"
  "	ntau (ms)\n"
  "       \n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	rates(v, cai)\n"
  "	m = minf\n"
  "	n = ninf\n"
  "    h = hinf   \n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD derivimplicit\n"
  "	ica = gCaL*m*n*h*(v - 65)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "	rates(v, cai)\n"
  "	m' = (minf - m)/mtau\n"
  "    n' = (ninf - n)/ntau\n"
  "	h' = (hinf - h)/2\n"
  "}\n"
  "\n"
  "FUNCTION alp(v(mV),i) (/ms) { LOCAL  q10\n"
  "	q10 = 3^((celsius - 37(degC))/10(degC))\n"
  "	if (i==0) {\n"
  "		alp = (1 - exp(-(v + 10)/6.24))/(0.035*(v + 10)*(1 + exp(-(v + 10)/6.24)))/(q10*1(/ms))\n"
  "	}else if (i==1){\n"
  "		alp = 9/(0.0197*exp(-0.0337^2*(v + 10)^2) + 0.02)/(q10*1(/ms))\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION bet(v(mV),i)(/ms) { \n"
  "	 v = v \n"
  "	\n"
  "	if (i==0) {\n"
  "		bet = 1/(1 + exp(-(v + 10)/8(mV)))\n"
  "	}else if (i==1){\n"
  "		bet = 1/(1 + exp((v + 28)/6.9(mV)))\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION ce(cai(mM)) { \n"
  "            \n"
  "           ce = 1/(1 + (cai/0.00035))\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "PROCEDURE rates(v(mV), cai (mM)) (/ms) {LOCAL a, b, c \n"
  "	\n"
  "	a = alp(v,0)  b=bet(v,0)\n"
  "	mtau = a\n"
  "	minf = b\n"
  "	a = alp(v,1)  b = bet(v,1)\n"
  "	ntau = a\n"
  "	ninf = b\n"
  "        c = ce(cai)\n"
  "	hinf = c\n"
  "}\n"
  ;
#endif
