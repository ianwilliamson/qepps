//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QEPPS: Quadratic eigenvalue problem parameter sweeper
//
// Copyright (C) 2014 Lab for Active Nano Devices, UT ECE 
// Developed by Ian Williamson 
// Supervised by Dr. Zheng Wang 
//
//-----------------------------------------------------------------------el-
// 
// The help message that's fed into SlepcInitialize()
// 
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef QEPPS_HELP
#define QEPPS_HELP

static const char help[] = "QEPPS: Quadratic eigenvalue problem parameter-sweeper\n\
Loads matrices representing quadratic eigenvalue problem (QEP) and solves\n\
( lambda^2*E + lambda*D + K )*U = 0, where E, D, and K are matrices, U is a\n\
vector, and lambda is an eigenvalue. The matrix inputs (E, D, and K) are\n\
each be specified in components that are each scaled by a function of the sweep\n\
parameter and then combined.";

#endif
