/*
Copyright (c) 2010-2015 Mate Soos

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
*/

#include <ctime>
#include <cstring>
#include <errno.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <string.h>

#include "main.h"
#include "time_mem.h"
#include "dimacsparser.h"
#include "cmsgen.h"
#include "signalcode.h"
#include "argparse.hpp"

using namespace CMSGen;
argparse::ArgumentParser program = argparse::ArgumentParser("cmsgen");

using std::cout;
using std::cerr;
using std::endl;

Main::Main(int _argc, char** _argv) :
    argc(_argc)
    , argv(_argv)
    , fileNamePresent (false)
{
}

void Main::readInAFile(SATSolver* solver2, const string& filename)
{
    if (conf.verbosity) {
        cout << "c Reading file '" << filename << "'" << endl;
    }
    #ifndef USE_ZLIB
    FILE * in = fopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<FILE*, FN> > parser(solver2, &debugLib, conf.verbosity);
    #else
    gzFile in = gzopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<gzFile, GZ> > parser(solver2, &debugLib, conf.verbosity);
    #endif

    if (in == NULL) {
        std::cerr
        << "ERROR! Could not open file '"
        << filename
        << "' for reading: " << strerror(errno) << endl;

        std::exit(1);
    }

    bool strict_header = true;
    if (!parser.parse_DIMACS(in, strict_header)) {
        exit(-1);
    }

    sampling_vars.swap(parser.sampling_vars);
    if (sampling_vars.empty()) {
        if (only_sampling_solution) {
            cout << "ERROR: only sampling vars are requested in the solution, but no sampling vars have been set!" << endl;
            exit(-1);
        }
    } else {
        solver2->set_sampling_vars(&sampling_vars);
        if (sampling_vars.size() > 100) {
            cout
            << "c Sampling var set contains over 100 variables, not displaying"
            << endl;
        } else {
            cout << "c Sampling vars set (total num: "
            << sampling_vars.size() << " ) : ";
            for(size_t i = 0; i < sampling_vars.size(); i++) {
                const uint32_t v = sampling_vars[i];
                cout << v+1;
                if (i+1 != sampling_vars.size())
                    cout << ",";
            }
            cout << endl;
        }
    }
    #ifndef USE_ZLIB
        fclose(in);
    #else
        gzclose(in);
    #endif
}

void Main::readInStandardInput(SATSolver* solver2)
{
    if (conf.verbosity) {
        cout
        << "c Reading from standard input... Use '-h' or '--help' for help."
        << endl;
    }

    #ifndef USE_ZLIB
    FILE * in = stdin;
    #else
    gzFile in = gzdopen(0, "rb"); //opens stdin, which is 0
    #endif

    if (in == NULL) {
        std::cerr << "ERROR! Could not open standard input for reading" << endl;
        std::exit(1);
    }

    #ifndef USE_ZLIB
    DimacsParser<StreamBuffer<FILE*, FN> > parser(solver2, &debugLib, conf.verbosity);
    #else
    DimacsParser<StreamBuffer<gzFile, GZ> > parser(solver2, &debugLib, conf.verbosity);
    #endif

    if (!parser.parse_DIMACS(in, false)) {
        exit(-1);
    }

    #ifdef USE_ZLIB
        gzclose(in);
    #endif
}

void Main::parseInAllFiles(SATSolver* solver2)
{
    const double myTimeTotal = cpuTimeTotal();
    const double myTime = cpuTime();
    if (!fileNamePresent) readInStandardInput(solver2);
    else readInAFile(solver2, fileToRead.c_str());

    if (conf.verbosity) {
        if (num_threads > 1) {
            cout
            << "c Sum parsing time among all threads (wall time will differ): "
            << std::fixed << std::setprecision(2)
            << (cpuTimeTotal() - myTimeTotal)
            << " s" << endl;
        } else {
            cout
            << "c Parsing time: "
            << std::fixed << std::setprecision(2)
            << (cpuTime() - myTime)
            << " s" << endl;
        }
    }
}

void Main::printResultFunc(
    std::ostream* os
    , const bool toFile
    , const lbool ret
) {
    assert(toFile);
    assert(ret == l_True || ret == l_False);
    if (ret == l_True) {
        for (uint32_t var = 0; var < solver->nVars(); var++) {
            if (solver->get_model()[var] != l_Undef) {
                *os << ((solver->get_model()[var] == l_True)? "" : "-") << var+1 << " ";
            }
        }
        *os << "0" << endl;
    } else if (ret == l_False) {
        cout << "WARNING: No samples generated, CNF is unsatisfiable" << endl;
    }
}

/* clang-format off */
void Main::add_supported_options() {
    // Declare the supported options.
    program.add_argument("-v", "--verb")
        .action([&](const auto& a) {conf.verbosity = std::atoi(a.c_str());})
        .default_value(conf.verbosity)
        .help("verbosity");
    program.add_argument("-s", "--seed")
        .action([&](const auto& a) {conf.origSeed = std::atoi(a.c_str());})
        .default_value(conf.origSeed)
        .help("Seed");
    program.add_argument("--samples")
        .action([&](const auto& a) {max_nr_of_solutions= std::atoi(a.c_str());})
        .default_value(conf.origSeed)
        .help("Number of samples needed");
    program.add_argument("--fixedconfl")
        .action([&](const auto& a) {conf.fixed_restart_num_confl = std::atoi(a.c_str());})
        .default_value(conf.fixed_restart_num_confl)
        .help("In case fixed restart strategy is used, how many conflicts should elapse between restarts");
    program.add_argument("--samplefile")
        .action([&](const auto& a) {resultFilename = a;})
        .help("Write sample(s) to this file");
    program.add_argument("file").remaining().help("input CNF file");
}
/* clang-format on */

void Main::manually_parse_some_options() {
    resultfile = new std::ofstream;
    resultfile->open(resultFilename.c_str());
    if (!(*resultfile)) {
        cout
        << "ERROR: Couldn't open file '"
        << resultFilename
        << "' for writing samples!"
        << endl;
        std::exit(-1);
    }
    conf.polarity_mode = PolarityMode::polarmode_weighted;
    conf.restartType = Restart::fixed;

    vector<std::string> file;
    try {
        file = program.get<std::vector<std::string>>("file");
        if (file.size() > 1) {
            cout << "ERROR: you can only pass at most one positional parameters: an INPUT file" << endl;
            exit(-1);
        }
        fileToRead = file[0];
        fileNamePresent = true;
    } catch (std::logic_error& e) {
        fileNamePresent = false;
    }

    if (conf.verbosity >= 3) cout << "c Outputting solution to console" << endl;
}

void Main::parseCommandLine()
{
    need_clean_exit = 0;

    //Reconstruct the command line so we can emit it later if needed
    for(int i = 0; i < argc; i++) {
        command_line += string(argv[i]);
        if (i+1 < argc) command_line += " ";
    }

    add_supported_options();
    try {
        program.parse_args(argc, argv);
        if (program.is_used("--help")) {

            cout
            << "A sampler that tries to output uniform samples" << endl
            << "Input "
            #ifndef USE_ZLIB
            << "must be plain"
            #else
            << "can be either plain or gzipped"
            #endif
            << " DIMACS with XOR extension" << endl << endl;

            cout
            << "USAGE: " << argv[0] << " --samples [NUM]  -s [SEED] --samplefile [OUTFILE] [CNFFILE]" << endl;
            cout << program << endl;
            std::exit(0);
        }
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        exit(-1);
    }

    printVersionInfo();
    if (program["version"] == true) std::exit(0);
    cout << "c executed with command line: " << command_line << endl;

    manually_parse_some_options();
}

int Main::solve() {
    double myTime = cpuTime();
    solver = new SATSolver((void*)&conf);
    solverToInterrupt = solver;

    //Print command line used to execute the solver: for options and inputs
    if (conf.verbosity) {
        printVersionInfo();
        cout << "c Executed with command line: " << command_line << endl;
    }

    //Parse in DIMACS (maybe gzipped) files
    parseInAllFiles(solver);

    lbool ret = multi_solutions();
    //printResultFunc(&cout, false, ret);
    assert(resultfile);
    printResultFunc(resultfile, true, ret);
    if (ret == l_True) {
        cout << "c Finished generating all " << max_nr_of_solutions << " samples" << endl;
    }
    cout << "c Total time: " << std::setprecision(2) << (cpuTime()-myTime) << " s " << endl;

    return correctReturnValue(ret);
}

lbool Main::multi_solutions() {
    cout << "c Writing samples to file: " << resultFilename << endl;
    unsigned long current_nr_of_solutions = 0;
    lbool ret = l_True;
    while(current_nr_of_solutions < max_nr_of_solutions && ret == l_True) {
        ret = solver->solve(NULL, only_sampling_solution);
        current_nr_of_solutions++;
        if (ret == l_True && current_nr_of_solutions < max_nr_of_solutions) {
            //printResultFunc(&cout, false, ret);
            if (resultfile) printResultFunc(resultfile, true, ret);
            if (current_nr_of_solutions % 10 == 0) {
                cout
                << "c Number of samples found until now: "
                << std::setw(6) << current_nr_of_solutions
                << endl;
            }
        }
    }
    return ret;
}

void Main::printVersionInfo() {
    cout << solver->get_text_version_info();
}

int Main::correctReturnValue(const lbool ret) const {
    int retval = -1;
    if (ret == l_True) {
        retval = 10;
    } else if (ret == l_False) {
        retval = 20;
    } else if (ret == l_Undef) {
        retval = 15;
    } else {
        std::cerr << "Something is very wrong, output is neither l_Undef, nor l_False, nor l_True" << endl;
        exit(-1);
    }

    if (zero_exit_status) {
        return 0;
    } else {
        return retval;
    }
}
