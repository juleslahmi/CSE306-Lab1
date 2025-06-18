#pragma once
#include "polygon.hpp"
#include <vector>
#include <string>
#include <iostream>

void save_svg(const std::vector<MyPolygon>&, const std::string&, const std::string& fillcol = "none");
void save_svg_animated(const std::vector<MyPolygon>&, std::string, int, int);