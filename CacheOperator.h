/**
* @file  CacheOperator.h
* @brief Defines the data and methods of the CacheOperator class
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef CACHEOPERATOR_H
#define CACHEOPERATOR_H

#include "Operator.h"

class CacheOperator : public Operator {
public:
	CacheOperator(Cloud * const C) : Operator(C) {}
	~CacheOperator() {}
    
	void operation1(const double currentTime);
	void operation2(const double currentTime);
	void operation3(const double currentTime);
	void operation4(const double currentTime);
};

#endif // CACHEOPERATOR_H
