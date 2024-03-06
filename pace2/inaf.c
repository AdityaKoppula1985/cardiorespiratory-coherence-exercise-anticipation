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
 
#define nrn_init _nrn_init__INaf
#define _nrn_initial _nrn_initial__INaf
#define nrn_cur _nrn_cur__INaf
#define _nrn_current _nrn_current__INaf
#define nrn_jacob _nrn_jacob__INaf
#define nrn_state _nrn_state__INaf
#define _net_receive _net_receive__INaf 
#define rates rates__INaf 
#define states states__INaf 
 
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
#define gnabar _p[0]
#define gnabar_columnindex 0
#define Tauact _p[1]
#define Tauact_columnindex 1
#define Tauinactf _p[2]
#define Tauinactf_columnindex 2
#define Tauinacts _p[3]
#define Tauinacts_columnindex 3
#define ina _p[4]
#define ina_columnindex 4
#define m _p[5]
#define m_columnindex 5
#define h _p[6]
#define h_columnindex 6
#define n _p[7]
#define n_columnindex 7
#define Dm _p[8]
#define Dm_columnindex 8
#define Dh _p[9]
#define Dh_columnindex 9
#define Dn _p[10]
#define Dn_columnindex 10
#define ena _p[11]
#define ena_columnindex 11
#define _g _p[12]
#define _g_columnindex 12
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
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
 static void _hoc_expMe(void);
 static void _hoc_expMc(void);
 static void _hoc_expMd(void);
 static void _hoc_expMb(void);
 static void _hoc_expMa(void);
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
 "setdata_INaf", _hoc_setdata,
 "alp_INaf", _hoc_alp,
 "bet_INaf", _hoc_bet,
 "expMe_INaf", _hoc_expMe,
 "expMc_INaf", _hoc_expMc,
 "expMd_INaf", _hoc_expMd,
 "expMb_INaf", _hoc_expMb,
 "expMa_INaf", _hoc_expMa,
 "rates_INaf", _hoc_rates,
 0, 0
};
#define alp alp_INaf
#define bet bet_INaf
#define expMe expMe_INaf
#define expMc expMc_INaf
#define expMd expMd_INaf
#define expMb expMb_INaf
#define expMa expMa_INaf
 extern double alp( double , double );
 extern double bet( double , double );
 extern double expMe( double , double );
 extern double expMc( double , double );
 extern double expMd( double , double );
 extern double expMb( double , double );
 extern double expMa( double , double );
 /* declare global and static user variables */
#define htau htau_INaf
 double htau = 0;
#define hinf hinf_INaf
 double hinf = 0;
#define mtau mtau_INaf
 double mtau = 0;
#define minf minf_INaf
 double minf = 0;
#define ntau ntau_INaf
 double ntau = 0;
#define ninf ninf_INaf
 double ninf = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gnabar_INaf", 0, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "mtau_INaf", "ms",
 "htau_INaf", "ms",
 "ntau_INaf", "ms",
 "gnabar_INaf", "S/cm2",
 "Tauact_INaf", "ms",
 "Tauinactf_INaf", "ms",
 "Tauinacts_INaf", "ms",
 "ina_INaf", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double n0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "minf_INaf", &minf_INaf,
 "hinf_INaf", &hinf_INaf,
 "ninf_INaf", &ninf_INaf,
 "mtau_INaf", &mtau_INaf,
 "htau_INaf", &htau_INaf,
 "ntau_INaf", &ntau_INaf,
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
"INaf",
 "gnabar_INaf",
 "Tauact_INaf",
 "Tauinactf_INaf",
 "Tauinacts_INaf",
 0,
 "ina_INaf",
 0,
 "m_INaf",
 "h_INaf",
 "n_INaf",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	gnabar = 0.0008;
 	Tauact = 0.2;
 	Tauinactf = 1;
 	Tauinacts = 1;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
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

 void _inaf_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 13, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 INaf inaf.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zk ;
static int _reset;
static char *modelname = "Cardiac fast sodium current";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double);
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
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Dn = ( ninf - n ) / ntau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ntau )) ;
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
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Dn = ( ninf - n ) / ntau ;
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
 _lv = _lv ;
   _lq10 = pow( 3.0 , ( ( celsius - 37.0 ) / 10.0 ) ) ;
   if ( _li  == 0.0 ) {
     _lalp = _lq10 * 0.32 * expMa ( _threadargscomma_ _lv * 1.0 + 47.13 , 47.13 ) ;
     }
   else if ( _li  == 1.0 ) {
     _lalp = _lq10 * 0.135 * expMb ( _threadargscomma_ _lv * 1.0 + 80.0 , 6.8 ) ;
     }
   else if ( _li  == 2.0 ) {
     _lalp = _lq10 * 1.0 * expMd ( _threadargscomma_ _lv * 1.0 , 79.23 ) ;
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
 double _lq10 ;
 _lv = _lv ;
   _lq10 = pow( 3.0 , ( ( celsius - 37.0 ) / 10.0 ) ) ;
   if ( _li  == 0.0 ) {
     _lbet = _lq10 * 0.08 * exp ( - _lv / 11.0 ) ;
     }
   else if ( _li  == 1.0 ) {
     _lbet = _lq10 * 1.0 * expMc ( _threadargscomma_ _lv * 1.0 , 11.1 ) ;
     }
   else if ( _li  == 2.0 ) {
     _lbet = _lq10 * 1.0 * expMe ( _threadargscomma_ _lv * 1.0 , 40.14 ) ;
     }
   
return _lbet;
 }
 
static void _hoc_bet(void) {
  double _r;
   _r =  bet (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double expMb (  double _lx , double _ly ) {
   double _lexpMb;
 if ( v < - 40.0 ) {
     _lexpMb = ( exp ( - _lx / _ly ) ) ;
     }
   else {
     _lexpMb = 0.0 ;
     }
   
return _lexpMb;
 }
 
static void _hoc_expMb(void) {
  double _r;
   _r =  expMb (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double expMc (  double _lx , double _ly ) {
   double _lexpMc;
 if ( v < - 40.0 ) {
     _lexpMc = 3.56 * exp ( _lx * 0.079 ) + 3.1e5 * exp ( 0.35 * _lx ) ;
     }
   else {
     _lexpMc = 1.0 / ( 0.13 * ( 1.0 + ( exp ( - ( _lx + 10.66 ) / _ly ) ) ) ) ;
     }
   
return _lexpMc;
 }
 
static void _hoc_expMc(void) {
  double _r;
   _r =  expMc (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double expMd (  double _lx , double _ly ) {
   double _lexpMd;
 if ( v < - 40.0 ) {
     _lexpMd = ( - 1.2740e5 * exp ( _lx * 0.2444 ) - 3.474e-5 * exp ( - 0.0439 * _lx ) ) * ( _lx + 37.78 ) / ( 1.0 + exp ( 0.311 * ( _lx + _ly ) ) ) ;
     }
   else {
     _lexpMd = 0.0 ;
     }
   
return _lexpMd;
 }
 
static void _hoc_expMd(void) {
  double _r;
   _r =  expMd (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double expMe (  double _lx , double _ly ) {
   double _lexpMe;
 if ( v < - 40.0 ) {
     _lexpMe = 0.1212 * exp ( - _lx * 0.01052 ) / ( 1.0 + exp ( - 0.137 * ( _lx + _ly ) ) ) ;
     }
   else {
     _lexpMe = 0.3 * exp ( - _lx * 2.535e-7 ) / ( 1.0 + exp ( - 0.1 * ( _lx + 32.0 ) ) ) ;
     }
   
return _lexpMe;
 }
 
static void _hoc_expMe(void) {
  double _r;
   _r =  expMe (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double expMa (  double _lx , double _ly ) {
   double _lexpMa;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lexpMa = 10.0 * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lexpMa = _lx / ( 1.0 - exp ( - 0.1 * _lx ) ) ;
     }
   
return _lexpMa;
 }
 
static void _hoc_expMa(void) {
  double _r;
   _r =  expMa (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv ) {
   double _la , _lb ;
 _la = alp ( _threadargscomma_ _lv , 0.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 0.0 ) ;
   mtau = 1.0 / ( _la + _lb ) * Tauact ;
   minf = _la / ( _la + _lb ) ;
   _la = alp ( _threadargscomma_ _lv , 1.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 1.0 ) ;
   htau = 1.0 / ( _la + _lb ) * Tauinactf ;
   hinf = _la / ( _la + _lb ) ;
   _la = alp ( _threadargscomma_ _lv , 2.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 2.0 ) ;
   ntau = 1.0 / ( _la + _lb ) * Tauinacts ;
   ninf = _la / ( _la + _lb ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
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
  ena = _ion_ena;
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
  ena = _ion_ena;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
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
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   n = ninf ;
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
  ena = _ion_ena;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ina = gnabar * m * m * m * h * n * ( v - ena ) ;
   }
 _current += ina;

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
  ena = _ion_ena;
 _g = _nrn_current(_v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
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
  ena = _ion_ena;
 { error = _deriv1_advance = 1;
 derivimplicit(_ninits, 3, _slist1, _dlist1, _p, &t, dt, states, &_temp1);
_deriv1_advance = 0;
 if(error){fprintf(stderr,"at line 51 in file inaf.mod:\n	SOLVE states METHOD derivimplicit\n"); nrn_complain(_p); abort_run(error);}
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
 _slist1[1] = h_columnindex;  _dlist1[1] = Dh_columnindex;
 _slist1[2] = n_columnindex;  _dlist1[2] = Dn_columnindex;
 _slist2[0] = h_columnindex;
 _slist2[1] = m_columnindex;
 _slist2[2] = n_columnindex;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "inaf.mod";
static const char* nmodl_file_text = 
  "TITLE Cardiac fast sodium current\n"
  ": Hodgkin - Huxley type sodium channel from Courtemanche et al Am J Physiol 1998 275:H301\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX INaf\n"
  "	USEION na READ ena WRITE ina\n"
  "	RANGE gnabar, ina, Tauact, Tauinactf, Tauinacts\n"
  "	GLOBAL minf, hinf, ninf, mtau, htau, ntau\n"
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
  "	gnabar=.0008(S/cm2) <0,1e9> \n"
  "	Tauact= 0.2 (ms) \n"
  "	Tauinactf=1 (ms) \n"
  "	Tauinacts=1 (ms)\n"
  "   \n"
  "               \n"
  "}\n"
  "\n"
  "STATE {\n"
  "	m h n\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	celsius (degC) : 37\n"
  "	ena (mV)\n"
  "	ina (mA/cm2)\n"
  "	minf hinf ninf\n"
  "	mtau (ms)\n"
  "	htau (ms)\n"
  "        ntau (ms)  \n"
  "}\n"
  "LOCAL k\n"
  "INITIAL {\n"
  "	rates(v)\n"
  "	m = minf\n"
  "	h = hinf\n"
  "        n = ninf   \n"
  "	\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD derivimplicit\n"
  "	ina = gnabar*m*m*m*h*n*(v - ena)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "	rates(v)\n"
  "	m' = (minf - m)/mtau\n"
  "	h' = (hinf - h)/htau\n"
  "        n' = (ninf - n)/ntau\n"
  "}\n"
  "UNITSOFF\n"
  "FUNCTION alp(v(mV),i) (/ms) { LOCAL q10 :rest = -70  order m,h,j\n"
  "	v = v \n"
  "	q10 = 3^((celsius - 37(degC))/10(degC))\n"
  "	if (i==0) {\n"
  "		alp = q10*0.32(/ms)*expMa(v *1(/mV) + 47.13, 47.13)\n"
  "	}else if (i==1){\n"
  "		alp = q10*0.135(/ms)*expMb(v *1(/mV) + 80, 6.8)\n"
  "	}else if (i==2){\n"
  "		alp = q10*1(/ms)*expMd(v *1(/mV), 79.23)\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION bet(v(mV),i)(/ms) { LOCAL q10 :rest = -70  order m,h,n\n"
  "	 v = v \n"
  "	q10 = 3^((celsius - 37(degC))/10(degC))\n"
  "	if (i==0) {\n"
  "		bet = q10*0.08(/ms)*exp(-v/11(mV))\n"
  "	}else if (i==1){\n"
  "		bet = q10*1(/ms)*expMc(v *1(/mV), 11.1)\n"
  "	}else if (i==2){\n"
  "		bet = q10*1(/ms)*expMe(v *1(/mV), 40.14)\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION expMb(x,y) {\n"
  "	if (v < -40) {\n"
  "		expMb = (exp(-x/y))\n"
  "	} else{\n"
  "	expMb = 0\n"
  "		\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION expMc(x,y) {\n"
  "	if (v < -40) {\n"
  "		expMc = 3.56*exp(x*0.079) + 3.1e5*exp(0.35*x)\n"
  "	}else{\n"
  "		expMc = 1/(0.13*(1 + (exp(-(x + 10.66)/y))))\n"
  "	}\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION expMd(x,y) {\n"
  "	if (v < -40) {\n"
  "		expMd = (-1.2740e5*exp(x*0.2444) - 3.474e-5*exp(-0.0439*x))*(x + 37.78)/(1 + exp(0.311*(x + y)))\n"
  "	}else{\n"
  "		expMd = 0\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION expMe(x,y) {\n"
  "	if (v < -40) {\n"
  "		expMe = 0.1212*exp(-x*0.01052)/(1 + exp(-0.137*(x + y)))\n"
  "	}else{\n"
  "		expMe = 0.3*exp(-x*2.535e-7)/(1 + exp(-0.1*(x + 32)))\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION expMa(x,y) {\n"
  "	if (fabs(x/y) < 1e-6) {\n"
  "		expMa = 10*(1 - x/y/2)\n"
  "	}else{\n"
  "		expMa = x/(1 - exp(-0.1*x))\n"
  "	}\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE rates(v) {LOCAL a, b\n"
  "	:TABLE minf, hinf, ninf, mtau, htau, ntau DEPEND celsius FROM -100 TO 100 WITH 200\n"
  "	a = alp(v,0)  b=bet(v,0)\n"
  "	mtau = 1/(a + b)*Tauact\n"
  "	minf = a/(a + b)\n"
  "	a = alp(v,1)  b=bet(v,1)\n"
  "	htau = 1/(a + b)*Tauinactf\n"
  "	hinf = a/(a + b)\n"
  "a = alp(v,2)  b=bet(v,2)\n"
  "	ntau = 1/(a + b)*Tauinacts\n"
  "	ninf = a/(a + b)\n"
  "}\n"
  "UNITSON\n"
  ;
#endif
