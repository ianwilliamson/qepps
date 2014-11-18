#ifndef QEPPS_HELP
#define QEPPS_HELP

static const char help[] = "qepps : Quadratic Eigenvalue Problem Parameter-sweeper\n\
Loads matrices representing quadratic eigenvalue problem (QEP) and solves\n\
(lambda^2*E+lambda*D+K)*U=0, where E, D, and K are matrices, U is a\n\
vector, and lambda is an eigenvalue. The matrix inputs (E, D, and K) are\n\
each be specified in components that can be combined in the form\n\
E=E0+p*E1+p^2*E2, where p is a parameter that will be swept.\n\
\n\
Parameters:\n\ 
    Problem data:\n\
    -freqs <input_file> : Specifies the ASCII text file from which the frequencies (one per line) should be read\n\
    -E0,E1,E2 <input_file> : E matricies in COMSOL QEP\n\
    -D0,D1,D2 <input_file> : D matricies in COMSOL QEP\n\
    -K0,K1,K2 <input_file> : K matricies in COMSOL QEP.\n\
\n\
    Relevant petsc and slepc options:\n\
    -qep_target_magnitude\n\
    -qep_target <complex_number>\n\
    -qep_nev <number>\n\
    -qep_ncv <number>\n\
    -qep_mpd <number>\n\
    -qep_type <type>\n\
    -st_shift <type>\n";

#endif
