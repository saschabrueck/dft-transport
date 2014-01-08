/*******************************************************************************
!                              INTEL CONFIDENTIAL
!   Copyright(C) 2008 Intel Corporation. All Rights Reserved.
!   The source code contained  or  described herein and all documents related to
!   the source code ("Material") are owned by Intel Corporation or its suppliers
!   or licensors.  Title to the  Material remains with  Intel Corporation or its
!   suppliers and licensors. The Material contains trade secrets and proprietary
!   and  confidential  information of  Intel or its suppliers and licensors. The
!   Material  is  protected  by  worldwide  copyright  and trade secret laws and
!   treaty  provisions. No part of the Material may be used, copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way without Intel's prior express written permission.
!   No license  under any  patent, copyright, trade secret or other intellectual
!   property right is granted to or conferred upon you by disclosure or delivery
!   of the Materials,  either expressly, by implication, inducement, estoppel or
!   otherwise.  Any  license  under  such  intellectual property  rights must be
!   express and approved by Intel in writing.
!******************************************************************************/

#ifndef _SPIKE_C
#include "spike.h"
#endif

void setAll(void *p, void *value, int type); 
void setAll_D(void *p, void *value); 
void setAll_I(void *p, void *value); 

void setRange(array_dim1 *p, void *value, int start, int end, int type);
void setRange_D(array_dim1 *p, double value, int start, int end);
void setRange_I(array_dim1 *p, int value, int start, int end);

void setRow(array_dim2 *p, int row, void *value, int type); 
void setRow_D(array_dim2 *p, int row, double value); 
void setRow_I(array_dim2 *p, int row, int value); 

void setElem1D_D(array_dim1 *p, int element_No, double value);
void setElem1D_I(array_dim1 *p, int element_No, int value);

void setElem2D(array_dim2 *p, int row, int column, void *value, int type); 
void setElem2D_D(array_dim2 *p, int row, int column, double value); 
void setElem2D_I(array_dim2 *p, int row, int column, int value); 

double getElem1D_D(array_dim1 *p, int element_No); 
int getElem1D_I(array_dim1 *p, int element_No);

double getElem2D_D(array_dim2 *p, int row, int column);  
int getElem2D_I(array_dim2 *p, int row, int column);

void dealloc(void *p);
void alloc1D_I(array_dim1 *p, int element_No);
void alloc2D_I(array_dim2 *p, int rowNo, int columNo);
void alloc3D_I(array_dim3 *p, int firstNo, int secondNo, int thirdNo);
void alloc1D_D(array_dim1 *p, int element_No);
void alloc2D_D(array_dim2 *p, int rowNo, int columNo);
void alloc3D_D(array_dim3 *p, int firstNo, int secondNo, int thirdNo);

