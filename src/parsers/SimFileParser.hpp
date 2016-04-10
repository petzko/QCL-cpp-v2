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


#include "../sim/SimSettings.hpp"

int parseFile(char * filename,SimSettings* simset);

std::vector<std::string>  get_option_from_file(char* filename, char* const optionname);
std::vector<std::string> tokenize(std::string str, std::string DELIMITER);
std::string normalizeStr(std::string str);


#endif /* SRC_PARSERS_SIMFILEPARSER_HPP_ */
