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
// The driver function and its auxiliaries
// 
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef QEPPS_SWEEPER
#define QEPPS_SWEEPER

/*!
 *  This assembles the system matricies for each parameter value and handles all function evaluation
 *  within the LUA state. Handles scaling/combining the component matricies.
 */
void qeppsSweeper(void);

#endif
