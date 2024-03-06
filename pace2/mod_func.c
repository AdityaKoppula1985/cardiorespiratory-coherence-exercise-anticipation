#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _ca_a1g_reg();
extern void _hhx_reg();
extern void _htc_reg();
extern void _ical_reg();
extern void _ik1_reg();
extern void _ikr_reg();
extern void _iks_reg();
extern void _inaf_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," ca_a1g.mod");
fprintf(stderr," hhx.mod");
fprintf(stderr," htc.mod");
fprintf(stderr," ical.mod");
fprintf(stderr," ik1.mod");
fprintf(stderr," ikr.mod");
fprintf(stderr," iks.mod");
fprintf(stderr," inaf.mod");
fprintf(stderr, "\n");
    }
_ca_a1g_reg();
_hhx_reg();
_htc_reg();
_ical_reg();
_ik1_reg();
_ikr_reg();
_iks_reg();
_inaf_reg();
}
