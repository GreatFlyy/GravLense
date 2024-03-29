#pragma once
#include<CCfits/CCfits>
#include<math.h>
#include<vector>

using namespace CCfits;

int readImageP(std::vector<std::vector<double>>& output, std::string& path);
int readImageE(std::vector<std::vector<double>>& output, std::string& path, short nExt);
int readHeader();
int writeImage(std::vector<std::vector<double>>& source, string& path);
int write2Image(std::vector<std::vector<double>>& source1, std::vector<std::vector<double>>& source2, string& path);


