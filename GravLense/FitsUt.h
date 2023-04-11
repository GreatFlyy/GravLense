#pragma once
#include<CCfits/CCfits>
#include<math.h>
#include<vector>

using namespace CCfits;

int readImage(std::vector<std::vector<double>>& output, std::auto_ptr<FITS>& pInfile, short nExt);
int readHeader();
int writeImage(std::vector<std::vector<double>>& source, string& path);


