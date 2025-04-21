/******************************************
Copyright (c) 2016, Mate Soos

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***********************************************/

#pragma once

#include <string>
#include <vector>
#include <fstream>

#include "main_common.h"
#include "solverconf.h"
#include "cmsgen.h"

using std::string;
using std::vector;

using namespace CMSGen;

class Main: public MainCommon
{
    public:
        Main(int argc, char** argv);
        ~Main() {
            delete solver;
        }

        void parseCommandLine();
        virtual int solve();

    private:
        //arguments
        int argc;
        char** argv;
        string var_elim_strategy;
        void manually_parse_some_options();
        void parse_restart_type();
        void parse_polarity_type();


    protected:
        //Options
        virtual void add_supported_options();
        virtual void call_after_parse() {}
        SATSolver* solver = NULL;

        //File reading
        void readInAFile(SATSolver* solver2, const string& filename);
        void readInStandardInput(SATSolver* solver2);
        void parseInAllFiles(SATSolver* solver2);

        //Helper functions
        void printResultFunc(
            std::ostream* os
            , const bool toFile
            , const lbool ret
        );
        void printVersionInfo();
        int correctReturnValue(const lbool ret) const;
        lbool multi_solutions();

        //Config
        std::string resultFilename = "samples.out";
        std::string debugLib;
        int printResult = true;
        string command_line;
        uint32_t max_nr_of_solutions = 100;
        int sql = 0;
        string decisions_for_model_fname;

        //Sampling vars
        vector<uint32_t> sampling_vars;
        std::string sampling_vars_str = "";
        bool only_sampling_solution = false;
        std::string assump_filename;
        vector<Lit> assumps;


        //Files to read & write
        bool fileNamePresent;
        string fileToRead;
        std::ofstream* resultfile = NULL;
        string dump_red_fname;
        uint32_t dump_red_max_len = 10000;
        uint32_t dump_red_max_glue = 1000;

        //Drat checker
        bool clause_ID_needed = false;
};
