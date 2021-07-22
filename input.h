#ifndef INPUT_H_
#define INPUT_H_

#include <iostream>
#include <stdlib.h>
#include <istream>
#include <fstream>
#include <assert.h>
#include <string>
#include "htslib/kseq.h"
#include "htslib/hts.h"

using namespace std;

// declare the type of file handler and the read() function  
// KSEQ_INIT(FILE*, read);
enum InputType {FASTA};

class Input 
{
public:
    InputType inputType;	
    ifstream strm;
    htsFile *htsfp;

    Input() 
    {
        htsfp = NULL;
    }
  
    ~Input() {};

    bool Initialize(char* filename) 
    {
        // Check if input is valid
        strm.open(filename);
        if (strm.good() == false or strm.eof()) {
        	cerr << "Cannot open target " << filename << endl;
        	exit(1);	
        }

        // Check if input is FASTA
        htsfp = hts_open(filename, "r");
        const htsFormat *fmt = hts_get_format(htsfp);
        string format = hts_format_file_extension(fmt);
        	
        if (format == "fa") {
        	inputType = FASTA;
        	return true;
        }
        else {
        	cerr << "Cannot determine type of input " << endl;
        	exit(1);
        }

        return false;
    }
};

#endif

