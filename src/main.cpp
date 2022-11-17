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

#if defined(__GNUC__) && defined(__linux__)

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <fenv.h>
#endif

#include <ctime>
#include <cstring>
#include <errno.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <fstream>
#include <sys/stat.h>
#include <string.h>
#include <list>
#include <array>
#include <thread>

#include "main.h"
#include "main_common.h"
#include "time_mem.h"
#include "dimacsparser.h"
#include "cmsgen/cryptominisat.h"
#include "signalcode.h"

using namespace CMSat;

using std::cout;
using std::cerr;
using std::endl;

struct WrongParam
{
    WrongParam(string _param, string _msg) :
        param(_param)
        , msg(_msg)
    {}

    const string& getMsg() const
    {
        return msg;
    }

    const string& getParam() const
    {
        return param;
    }

    string param;
    string msg;
};

bool fileExists(const string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}


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

    if (!sampling_vars_str.empty() && !parser.sampling_vars.empty()) {
        cerr << "ERROR! Sampling vars set in console but also in CNF." << endl;
        exit(-1);
    }

    if (!sampling_vars_str.empty()) {
        assert(sampling_vars.empty());

        std::stringstream ss(sampling_vars_str);
        uint32_t i;
        while (ss >> i)
        {
            const uint32_t var = i-1;
            sampling_vars.push_back(var);

            if (ss.peek() == ',' || ss.peek() == ' ')
                ss.ignore();
        }
    } else {
        sampling_vars.swap(parser.sampling_vars);
    }

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
    call_after_parse();

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

    //First read normal extra files
    if (!debugLib.empty() && filesToRead.size() > 1) {
        cout
        << "ERROR: debugLib must be OFF"
        << " to parse in more than one file"
        << endl;

        std::exit(-1);
    }

    for (const string& fname: filesToRead) {
        readInAFile(solver2, fname.c_str());
    }

    if (!fileNamePresent) {
        readInStandardInput(solver2);
    }

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
void Main::add_supported_options()
{
    // Declare the supported options.
    generalOptions.add_options()
    ("help,h", "Print simple help")
    ("version,v", "Print version info")
    ("verb", po::value(&conf.verbosity)->default_value(conf.verbosity)
        , "[0-10] Verbosity")
    ("seed,s", po::value(&conf.origSeed)->default_value(conf.origSeed)
        , "[0..] Random seed")
    ("samples", po::value(&max_nr_of_solutions)->default_value(max_nr_of_solutions)
        , "Number of samples needed")
    ("fixedconfl", po::value(&conf.fixed_restart_num_confl)->default_value(conf.fixed_restart_num_confl)
        , "In case fixed restart strategy is used, how many conflicts should elapse between restarts")
    ("samplefile", po::value(&resultFilename)->default_value(resultFilename)
        , "Write sample(s) to this file")
    ;

    hiddenOptions.add_options()
    ("input", po::value< vector<string> >(), "file(s) to read")
    ;

    help_options_complicated
    .add(generalOptions)
    ;
}
/* clang-format on */

string remove_last_comma_if_exists(std::string s)
{
    std::string s2 = s;
    if (s[s.length()-1] == ',')
        s2.resize(s2.length()-1);
    return s2;
}

void Main::check_options_correctness()
{
    try {
        po::store(po::command_line_parser(argc, argv).options(all_options).positional(p).run(), vm);
        if (vm.count("help"))
        {
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
            cout << help_options_simple << endl;
            std::exit(0);
        }

        po::notify(vm);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::unknown_option> >& c
    ) {
        cerr
        << "ERROR: Some option you gave was wrong. Please give '--help' to get help" << endl
        << "       Unknown option: " << c.what() << endl;
        std::exit(-1);
    } catch (boost::bad_any_cast &e) {
        std::cerr
        << "ERROR! You probably gave a wrong argument type" << endl
        << "       Bad cast: " << e.what()
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_option_value> >& what
    ) {
        cerr
        << "ERROR: Invalid value '" << what.what() << "'" << endl
        << "       given to option '" << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::multiple_occurrences> >& what
    ) {
        cerr
        << "ERROR: " << what.what() << " of option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::required_option> >& what
    ) {
        cerr
        << "ERROR: You forgot to give a required option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::too_many_positional_options_error> >& what
    ) {
        cerr
        << "ERROR: You gave too many positional arguments. Only at most two can be given:" << endl
        << "       the 1st the CNF file input, and optionally, the 2nd the DRAT file output" << endl
        << "    OR (pre-processing)  1st for the input CNF, 2nd for the simplified CNF" << endl
        << "    OR (post-processing) 1st for the solution file" << endl
        ;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::ambiguous_option> >& what
    ) {
        cerr
        << "ERROR: The option you gave was not fully written and matches" << endl
        << "       more than one option. Please give the full option name." << endl
        << "       The option you gave: '" << what.get_option_name() << "'" <<endl
        << "       The alternatives are: ";
        for(size_t i = 0; i < what.alternatives().size(); i++) {
            cout << what.alternatives()[i];
            if (i+1 < what.alternatives().size()) {
                cout << ", ";
            }
        }
        cout << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_command_line_syntax> >& what
    ) {
        cerr
        << "ERROR: The option you gave is missing the argument or the" << endl
        << "       argument is given with space between the equal sign." << endl
        << "       detailed error message: " << what.what() << endl
        ;
        std::exit(-1);
    }
}

void Main::manually_parse_some_options()
{
    if (conf.yalsat_max_mems < 1) {
        cout << "ERROR: '--walkmems' must be at least 1" << endl;
        exit(-1);
    }

    if (conf.sls_every_n < 1) {
        cout << "ERROR: '--walkeveryn' must be at least 1" << endl;
        exit(-1);
    }

    if (conf.maxXorToFind > MAX_XOR_RECOVER_SIZE) {
        cout << "ERROR: The '--maxxorsize' parameter cannot be larger than " << MAX_XOR_RECOVER_SIZE << endl;
        exit(-1);
    }

    if (conf.shortTermHistorySize <= 0) {
        cout
        << "You MUST give a short term history size (\"--gluehist\")" << endl
        << "  greater than 0!"
        << endl;

        std::exit(-1);
    }

    if (!decisions_for_model_fname.empty()) {
        conf.need_decisions_reaching = true;
    }

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

    if (vm.count("input")) {
        filesToRead = vm["input"].as<vector<string> >();
        fileNamePresent = true;
    } else {
        fileNamePresent = false;
    }

    if (conf.verbosity >= 3) {
        cout << "c Outputting solution to console" << endl;
    }
}

void Main::parseCommandLine()
{
    need_clean_exit = 0;

    //Reconstruct the command line so we can emit it later if needed
    for(int i = 0; i < argc; i++) {
        commandLine += string(argv[i]);
        if (i+1 < argc) {
            commandLine += " ";
        }
    }

    add_supported_options();
    p.add("input", 1);
    p.add("drat", 1);
    all_options.add(help_options_complicated);
    all_options.add(hiddenOptions);

    help_options_simple
    .add(generalOptions)
    ;

    check_options_correctness();
    if (vm.count("version")) {
        printVersionInfo();
        std::exit(0);
    }

    try {
        manually_parse_some_options();
    } catch(WrongParam& wp) {
        cerr << "ERROR: " << wp.getMsg() << endl;
        exit(-1);
    }
}

int Main::solve()
{
    double myTime = cpuTime();
    solver = new SATSolver((void*)&conf);
    solverToInterrupt = solver;
    if (dratf) {
        solver->set_drat(dratf, clause_ID_needed);
    }
    solver->set_num_threads(num_threads);

    //Print command line used to execute the solver: for options and inputs
    if (conf.verbosity) {
        printVersionInfo();
        cout
        << "c Executed with command line: "
        << commandLine
        << endl;
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
    cout << "c Total time: " << (cpuTime()-myTime) << endl;

    return correctReturnValue(ret);
}

lbool Main::multi_solutions()
{
    cout << "c Writing samples to file: " << resultFilename << endl;
    unsigned long current_nr_of_solutions = 0;
    lbool ret = l_True;
    while(current_nr_of_solutions < max_nr_of_solutions && ret == l_True) {
        ret = solver->solve(&assumps, only_sampling_solution);
        current_nr_of_solutions++;
        if (ret == l_True && !decisions_for_model_fname.empty()) {
            solver->add_empty_cl_to_drat();
            assert(max_nr_of_solutions == 1);
        }

        if (ret == l_True && current_nr_of_solutions < max_nr_of_solutions) {
            //printResultFunc(&cout, false, ret);
            if (resultfile) {
                printResultFunc(resultfile, true, ret);
            }

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

///////////
// Useful helper functions
///////////

void Main::printVersionInfo()
{
    cout << solver->get_text_version_info();
}

int Main::correctReturnValue(const lbool ret) const
{
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
