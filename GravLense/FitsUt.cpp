#include"FitsUt.h"

int readImage(std::vector<std::vector<double>>& output, std::string& path, short nExt)
{
    std::auto_ptr<FITS> pInfile(new FITS(path, Read, true));
    ExtHDU& image = pInfile->extension(nExt);
    std::valarray<double>  contents(1, 0);
    // read all user-specifed, coordinate, and checksum keys in the image
    image.readAllKeys();
    image.read(contents);
    // this doesn't print the data, just header info.

    long ax1(image.axis(0));
    long ax2(image.axis(1));

    output.clear();
    output = std::vector<std::vector<double>>(ax2, std::vector<double>(ax1, 0.));

    for (int j = 0; j < ax2; j++)
    {
        for (int i = 0; i < ax1; i++)
        {
            output[j][i] = contents[j * ax1 + i];
        }
    }
    return 0;
}

int readHeader()
{
    const string SCI("ERR");
    // read a particular HDU within the file. This call reads just the header 
    // information from SPECTRUM
    std::auto_ptr<FITS> pInfile(new FITS("C:/Users/GreatFly/Desktop/HST/MAST_2022-11-10T0750/HST/na1l91010/na1l91010_mos.fits", Read, SCI));
    // define a reference for clarity. (std::auto_ptr<T>::get returns a pointer
    ExtHDU& table = pInfile->extension(SCI);
    // read all the keywords, excluding those associated with columns.
    table.readAllKeys();
    // print the result.
    std::cout << table << std::endl;
    return 0;
}




int writeImage(std::vector<std::vector<double>>& source, string& path)
{
    // Create a FITS primary array containing a 2-D image
    // declare axis arrays.
    long naxis = 2;
    long naxes[2] = { 0, 0 };

    naxes[0] = source.size();
    naxes[1] = source[0].size();

    // declare auto-pointer to FITS at function scope. Ensures no resources
    // leaked if something fails in dynamic allocation.
    std::auto_ptr<FITS> pFits(0);

    try
    {
        // overwrite existing file if the file already exists.

        const std::string fileName = path;

        // Create a new FITS object, specifying the data type and axes for the primary
        // image. Simultaneously create the corresponding file.

        // this image is unsigned short data, demonstrating the cfitsio extension
        // to the FITS standard.

        pFits.reset(new FITS(fileName, DOUBLE_IMG, naxis, naxes));
    }
    catch (FITS::CantCreate)
    {
        // ... or not, as the case may be.
        return -1;
    }

    // references for clarity.

    long& vectorLength = naxes[0];
    long& numberOfRows = naxes[1];
    long nelements(1);


    // Find the total size of the array. 
    // this is a little fancier than necessary ( It's only
    // calculating naxes[0]*naxes[1]) but it demonstrates  use of the 
    // C++ standard library accumulate algorithm.

    nelements = std::accumulate(&naxes[0], &naxes[naxis], 1, std::multiplies<long>());




    // create a new image extension with a 300x300 array containing float data.

    /*std::vector<long> extAx(2, 300);
    string newName("NEW-EXTENSION");
    ExtHDU* imageExt = pFits->addImage(newName, FLOAT_IMG, extAx);*/

    // create a dummy row with a ramp. Create an array and copy the row to 
    // row-sized slices. [also demonstrates the use of valarray slices].   
    // also demonstrate implicit type conversion when writing to the image:
    // input array will be of type float.

    std::valarray<int> row(vectorLength);
    for (long j = 0; j < vectorLength; ++j) row[j] = j;
    std::valarray<int> array(nelements);
    /* for (int i = 0; i < numberOfRows; ++i)
     {
         array[std::slice(vectorLength * static_cast<int>(i), vectorLength, 1)] = row + i;
     }*/



    for (int i = 0; i < vectorLength; i++)
    {
        for (int j = 0; j < numberOfRows; j++)
        {
            array[i * numberOfRows + j] = source[i][j];
        }
    }



    // create some data for the image extension.
    /*long extElements = std::accumulate(extAx.begin(), extAx.end(), 1, std::multiplies<long>());
    std::valarray<float> ranData(extElements);
    const float PIBY(3.14 / 150.);
    for (int jj = 0; jj < extElements; ++jj)
    {
        float arg = PIBY * jj;
        ranData[jj] = std::cos(arg);
    }*/

    long  fpixel(1);

    // write the image extension data: also demonstrates switching between
    // HDUs.


    //imageExt->write(fpixel, extElements, ranData);


    //add two keys to the primary header, one long, one complex.

    /*long exposure(1500);
    std::complex<float> omega(std::cos(2 * 3.14 / 3.), std::sin(2 * 3.14 / 3));
    pFits->pHDU().addKey("EXPOSURE", exposure, "Total Exposure Time");
    pFits->pHDU().addKey("OMEGA", omega, " Complex cube root of 1 ");*/





    // The function PHDU& FITS::pHDU() returns a reference to the object representing 
    // the primary HDU; PHDU::write( <args> ) is then used to write the data.

    pFits->pHDU().write(fpixel, nelements, array);


    // PHDU's friend ostream operator. Doesn't print the entire array, just the
    // required & user keywords, and is provided largely for testing purposes [see 
    // readImage() for an example of how to output the image array to a stream].

    std::cout << pFits->pHDU() << std::endl;
    return 0;
}