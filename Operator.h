/**
* @file  Operator.h
* @brief Defines the data and methods of the abstract Operator class
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef OPERATOR_H
#define OPERATOR_H

#include "Cloud.h"

class Operator {
public:
	Cloud * const cloud;	//<! The dust cloud for the simulation

	/**
	* @brief Constructor for the Operator class
	**/
	Operator(Cloud * const C) : cloud(C) {}

	/**
	* @brief Destructor for the Operator class
	**/
	virtual ~Operator() {}
    
    /**
    * @brief The first Operator operation
    **/
	virtual void operation1(const double currentTime)=0;

	/**
    * @brief The second Operator operation
    **/
	virtual void operation2(const double currentTime)=0;

	/**
    * @brief The third Operator operation
    **/
	virtual void operation3(const double currentTime)=0;

	/**
    * @brief The fourth Operator operation
    **/
	virtual void operation4(const double currentTime)=0;
};

#endif // OPERATOR_H
