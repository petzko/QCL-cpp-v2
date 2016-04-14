/*
 * SimFileParser.hpp
 *
 *  Created on: Apr 8, 2016
 *      Author: petzko
 */

#ifndef SRC_PARSERS_SIMFILEPARSER_HPP_
#define SRC_PARSERS_SIMFILEPARSER_HPP_

#include "../common/utils.hpp"
#include "../parsers/SimFileParser.hpp"


#include <stdlib.h>
#include <stdlib.h>

std::vector<std::string>  getOptionFromFile(char* filename, std::string option);
std::vector<std::string> tokenize(std::string str, std::string DELIMITER);
std::string normalizeStr(std::string str);




#endif /* SRC_PARSERS_SIMFILEPARSER_HPP_ */
