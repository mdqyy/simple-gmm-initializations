#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

#include "../base/base.h"

std::string getName(AlgorithmID algorithm);
std::string getName(MeasurementError me);
void getID(std::string algorithm, AlgorithmID& target);
void getInitMethod(std::string initmethod, InitMethod& target);
void getMeasurementError(std::string measurementError, MeasurementError& target);
void parseTestLabConfiguration(int, char*[]);


#endif