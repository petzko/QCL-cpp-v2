/*
 * SimFileParser.hpp
 *
 *  Created on: Apr 8, 2016
 *      Author: petzko
 */

#ifndef SRC_PARSERS_SIMFILEPARSER_HPP_
#define SRC_PARSERS_SIMFILEPARSER_HPP_

#include "../common/utils.hpp"
#include <stdlib.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

#include "../sim/SimSettings.hpp"

int parseFile(char * filename,SimSettings* simset);
char** get_option_from_file(char* filename, char* const optionname);



#endif /* SRC_PARSERS_SIMFILEPARSER_HPP_ */
