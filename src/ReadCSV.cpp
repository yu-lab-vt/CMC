
#include "ReadCSV.h" 
#include <string>
#include <cstring>


// UINT (1D)
uint* ReadCSV_UINT(char* fileName, uint nHeaderlines)
{

    string line;
    uint Len_line=0;
    int num = 0;
    std::string::size_type sz;
    uint offset = 0;
    uint iCol = 0, iRow = 0;

//   ifstream file("L:/CMC_C_R_code/CMC_V0/Debug/Omega2.txt");
    ifstream file(fileName);
    if (file.fail()) {
        std::cout << fileName << " does not exist!\n";
        exit(0);
    }

    
    int numMatrix[MAX_ROW_CSV];
    //char sep = ',';

    // Bypass the headerlines
    for (uint i=0; (i< nHeaderlines) & (!file.eof()); i++)
        getline(file, line, '\n');

//    std::cout << "Bypass the headerlines!\n";


    // Read data
    while (!file.eof())
    {
        getline(file, line, '\n');
        Len_line = line.length();
        if (Len_line == 0)
            continue;
        for (offset = 0, iCol=0; offset < Len_line; offset += (sz+1), iCol++ )
        {
            num = std::stoi(line.substr(offset), &sz);
            numMatrix[iRow] = num;
        }

        if (iCol > 1) {
            std::cout << "The column counts of the " << iRow + 1 + nHeaderlines <<" row is larger than 1!\n";
            exit(0);
            }

        iRow++;
        if (iRow + 1 > MAX_ROW_CSV ) {
            std::cout << "The the number of rows exceed the maximum allowed rows(" << MAX_ROW_CSV << ")\n";
            exit(0);
        }
    }

    // Restore to a malloc memory (1D)
    uint* numOut;
    numOut = CALLOC(uint,iRow + 1);
    for (uint i = 0; i <= iRow; i++)
        numOut[i] = numMatrix[i];

    std::cout << "Finish the reading!\n";

    return numOut;
}

// double (1D)
double* ReadCSV_DOUBLE(char* fileName, uint nHeaderlines, uint* p_nRow)
{

    string line;
    uint Len_line = 0;
    double num = 0;
    std::string::size_type sz;
    uint offset = 0;
    uint iCol = 0, iRow = 0;

    //   ifstream file("L:/CMC_C_R_code/CMC_V0/Debug/Omega2.txt");
    ifstream file(fileName);
    if (file.fail()) {
        std::cout << fileName << " does not exist!\n";
        exit(0);
    }


    double numMatrix[MAX_ROW_CSV] = {0};
    //char sep = ',';

    // Bypass the headerlines
    for (uint i = 0; (i < nHeaderlines) & (!file.eof()); i++)
        getline(file, line, '\n');

    //    std::cout << "Bypass the headerlines!\n";


        // Read data
    while (!file.eof())
    {
        getline(file, line, '\n');
        Len_line = line.length();
        if (Len_line == 0)
            continue;
        for (offset = 0, iCol = 0; offset < Len_line; offset += (sz + 1), iCol++)
        {
            num = std::stod(line.substr(offset), &sz);
            numMatrix[iRow] = num;
        }

        if (iCol > 1) {
            std::cout << "The column counts of the " << iRow + 1 + nHeaderlines << " row is larger than 1!\n";
            exit(0);
        }

        iRow++;
        if (iRow + 1 > MAX_ROW_CSV) {
            std::cout << "The the number of rows exceed the maximum allowed rows(" << MAX_ROW_CSV << ")\n";
            exit(0);
        }
    }

    // Restore to a malloc memory (1D)
    double* numOut;
    numOut = CALLOC(double, iRow);
    for (uint i = 0; i <= iRow; i++)
        numOut[i] = numMatrix[i];

    std::cout << "Finish the reading!\n";

    *p_nRow = iRow;
    return numOut;
}

// UINT (2D)
uint** ReadCSV_UINT_2D(char* fileName, uint nHeaderlines)
{

    string line;
    uint Len_line = 0;
    int num = 0;
    std::string::size_type sz;
    uint offset = 0;
    uint iCol = 0, iRow = 0, nCol = 0;

    //   ifstream file("L:/CMC_C_R_code/CMC_V0/Debug/Omega2.txt");
    ifstream file(fileName);
    if (file.fail()) {
        std::cout << fileName << " does not exist!\n";
        exit(0);
    }

    int numMatrix[MAX_ROW_CSV][MAX_COL_CSV];
    //char sep = ',';

    // Bypass the headerlines
    for (uint i = 0; (i < nHeaderlines) & (!file.eof()); i++)
        getline(file, line, '\n');

    // Read data
    while (!file.eof())
    {
        getline(file, line, '\n');
        Len_line = line.length();
        if (Len_line == 0)
            continue;
        for (offset = 0, iCol = 0; offset < Len_line; offset += (sz + 1), iCol++)
        {
            num = std::stoi(line.substr(offset), &sz);
            numMatrix[iRow][iCol] = num;
        }
        if (iRow == 0)
            nCol = iCol + 1;
        else {
            if (nCol != (iCol + 1)) {
                std::cout << "The column counts of the " << iRow + 1 + nHeaderlines << " row is not the same as the first row!\n";
                exit(0);
            }
        }

        iRow++;
        if (iRow + 1 > MAX_ROW_CSV) {
            std::cout << "The the number of rows exceed the maximum allowed rows(" << MAX_ROW_CSV << ")\n";
            exit(0);
        }

    }

    // Restore to a malloc memory (2D)
    uint** numOut;
    numOut = CALLOC_2D_UINT(iRow + 1, iCol + 1);
    for (uint i = 0; i <= iRow; i++)
        for (uint j = 0; j <= iCol; j++)
            numOut[i][j] = numMatrix[i][j];

    std::cout << "Finish the reading!\n";


    return numOut;
}



// readBin (unknow size; malloc memory : free before the end of main() to aviod memory leak)
void readBin(double** buff, char* fileName, uint* size) {

    FILE* fp;
    exit(0);
    //fopen_s(&fp, fileName, "rb");
    if (fp == NULL) {
        cout << "Can not open the file " << fileName << "\n";
        exit(0);
    }


    fseek(fp, 0L, SEEK_END);
    uint sz = ftell(fp) / sizeof(double);
    fseek(fp, 0L, SEEK_SET);
    *size = sz;

    // Malloc memory (1D)
    *buff = CALLOC(double, sz);

    // Read
    size_t ret_code = fread(*buff, sizeof * buff, sz, fp); // reads an array of doubles
    if (ret_code == sz) {
        //printf("Array read successfully.\n");
    }
    else { // error handling
        if (feof(fp))
            printf("Error reading test.bin: unexpected end of file\n");
        else if (ferror(fp)) {
            perror("Error reading test.bin");
        }
        exit(0);
    }
    fclose(fp);
}



// readBin for unsigned int  // limit to file with size less than 2GB
void readBin(uint** buff, char* fileName, uint* size) {

    FILE* fp;
    exit(0);
    //fopen_s(&fp, fileName, "rb");
    if (fp == NULL) {
        cout << "Can not open the file " << fileName << "\n";
        exit(0);
    }


    fseek(fp, 0L, SEEK_END);
    uint sz = ftell(fp) / sizeof(uint);
    fseek(fp, 0L, SEEK_SET);
    *size = sz;

    // Malloc memory (1D)
    *buff = CALLOC(uint, sz);

    // Read
    size_t ret_code = fread(*buff, sizeof * *buff, sz, fp); // reads an array of doubles
    if (ret_code == sz) {
        //printf("Array read successfully.\n");
    }
    else { // error handling
        if (feof(fp))
            printf("Error reading *.bin: unexpected end of file\n");
        else if (ferror(fp)) {
            perror("Error reading *.bin");
        }
        exit(0);
    }
    fclose(fp);
}


// readBin for bool (same size with char) // for file more than 2GB;
void readBin(bool** buff, char* fileName, unsigned long long int* size) {

    FILE* fp;
    exit(0);
    //fopen_s(&fp, fileName, "rb");
    if (fp == NULL) {
        cout << "Can not open the file " << fileName << "\n";
        exit(0);
    }


    fseek(fp, 0L, SEEK_END);
    //size_t sz = _ftelli64(fp) / sizeof(bool);
    size_t sz = 0;
    fseek(fp, 0L, SEEK_SET);
    *size = sz;

    // Malloc memory (1D)
    *buff = CALLOC(bool, sz);
    
    // Read
    size_t ret_code = fread(*buff, sizeof * *buff, sz, fp); // reads an array of doubles
    if (ret_code == sz) {
        //printf("Array read successfully.\n");
    }
    else { // error handling
        if (feof(fp))
            printf("Error reading *.bin: unexpected end of file\n");
        else if (ferror(fp)) {
            perror("Error reading *.bin");
        }
        exit(0);
    }
    fclose(fp);
}

