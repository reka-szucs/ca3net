

#include "objects.h"
#include "synapses_classes.h"
#include "brianlib/clocks.h"
#include "brianlib/dynamic_array.h"
#include "brianlib/stdint_compat.h"
#include "network.h"
#include<random>
#include<vector>
#include<iostream>
#include<fstream>
#include<map>
#include<tuple>
#include<cstdlib>
#include<string>

namespace brian {

std::string results_dir = "results/";  // can be overwritten by --results_dir command line arg

// For multhreading, we need one generator for each thread. We also create a distribution for
// each thread, even though this is not strictly necessary for the uniform distribution, as
// the distribution is stateless.
std::vector< RandomGenerator > _random_generators;

//////////////// networks /////////////////
Network magicnetwork;

void set_variable_from_value(std::string varname, char* var_pointer, size_t size, char value) {
    #ifdef DEBUG
    std::cout << "Setting '" << varname << "' to " << (value == 1 ? "True" : "False") << std::endl;
    #endif
    std::fill(var_pointer, var_pointer+size, value);
}

template<class T> void set_variable_from_value(std::string varname, T* var_pointer, size_t size, T value) {
    #ifdef DEBUG
    std::cout << "Setting '" << varname << "' to " << value << std::endl;
    #endif
    std::fill(var_pointer, var_pointer+size, value);
}

template<class T> void set_variable_from_file(std::string varname, T* var_pointer, size_t data_size, std::string filename) {
    ifstream f;
    streampos size;
    #ifdef DEBUG
    std::cout << "Setting '" << varname << "' from file '" << filename << "'" << std::endl;
    #endif
    f.open(filename, ios::in | ios::binary | ios::ate);
    size = f.tellg();
    if (size != data_size) {
        std::cerr << "Error reading '" << filename << "': file size " << size << " does not match expected size " << data_size << std::endl;
        return;
    }
    f.seekg(0, ios::beg);
    if (f.is_open())
        f.read(reinterpret_cast<char *>(var_pointer), data_size);
    else
        std::cerr << "Could not read '" << filename << "'" << std::endl;
    if (f.fail())
        std::cerr << "Error reading '" << filename << "'" << std::endl;
}

//////////////// set arrays by name ///////
void set_variable_by_name(std::string name, std::string s_value) {
    size_t var_size;
    size_t data_size;
    // C-style or Python-style capitalization is allowed for boolean values
    if (s_value == "true" || s_value == "True")
        s_value = "1";
    else if (s_value == "false" || s_value == "False")
        s_value = "0";
    // non-dynamic arrays
    if (name == "spikegeneratorgroup._spikespace") {
        var_size = 8001;
        data_size = 8001*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_spikegeneratorgroup__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_spikegeneratorgroup__spikespace, data_size, s_value);
        }
        return;
    }
    // dynamic arrays (1d)
    if (name == "synapses.A_postsyn") {
        var_size = _dynamic_array_synapses_A_postsyn.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_synapses_A_postsyn[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_synapses_A_postsyn[0], data_size, s_value);
        }
        return;
    }
    if (name == "synapses.A_presyn") {
        var_size = _dynamic_array_synapses_A_presyn.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_synapses_A_presyn[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_synapses_A_presyn[0], data_size, s_value);
        }
        return;
    }
    if (name == "synapses.delay") {
        var_size = _dynamic_array_synapses_delay.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_synapses_delay[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_synapses_delay[0], data_size, s_value);
        }
        return;
    }
    if (name == "synapses.delay") {
        var_size = _dynamic_array_synapses_delay_1.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_synapses_delay_1[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_synapses_delay_1[0], data_size, s_value);
        }
        return;
    }
    if (name == "synapses.lastupdate") {
        var_size = _dynamic_array_synapses_lastupdate.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_synapses_lastupdate[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_synapses_lastupdate[0], data_size, s_value);
        }
        return;
    }
    if (name == "synapses.w") {
        var_size = _dynamic_array_synapses_w.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_synapses_w[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_synapses_w[0], data_size, s_value);
        }
        return;
    }
    std::cerr << "Cannot set unknown variable '" << name << "'." << std::endl;
    exit(1);
}
//////////////// arrays ///////////////////
double * _array_defaultclock_dt;
const int _num__array_defaultclock_dt = 1;
double * _array_defaultclock_t;
const int _num__array_defaultclock_t = 1;
int64_t * _array_defaultclock_timestep;
const int _num__array_defaultclock_timestep = 1;
int32_t * _array_spikegeneratorgroup__lastindex;
const int _num__array_spikegeneratorgroup__lastindex = 1;
int32_t * _array_spikegeneratorgroup__period_bins;
const int _num__array_spikegeneratorgroup__period_bins = 1;
int32_t * _array_spikegeneratorgroup__spikespace;
const int _num__array_spikegeneratorgroup__spikespace = 8001;
int32_t * _array_spikegeneratorgroup_i;
const int _num__array_spikegeneratorgroup_i = 8000;
double * _array_spikegeneratorgroup_period;
const int _num__array_spikegeneratorgroup_period = 1;
int32_t * _array_synapses_N;
const int _num__array_synapses_N = 1;

//////////////// dynamic arrays 1d /////////
std::vector<int32_t> _dynamic_array_spikegeneratorgroup__timebins;
std::vector<int32_t> _dynamic_array_spikegeneratorgroup_neuron_index;
std::vector<int32_t> _dynamic_array_spikegeneratorgroup_spike_number;
std::vector<double> _dynamic_array_spikegeneratorgroup_spike_time;
std::vector<int32_t> _dynamic_array_synapses__synaptic_post;
std::vector<int32_t> _dynamic_array_synapses__synaptic_pre;
std::vector<double> _dynamic_array_synapses_A_postsyn;
std::vector<double> _dynamic_array_synapses_A_presyn;
std::vector<double> _dynamic_array_synapses_delay;
std::vector<double> _dynamic_array_synapses_delay_1;
std::vector<double> _dynamic_array_synapses_lastupdate;
std::vector<int32_t> _dynamic_array_synapses_N_incoming;
std::vector<int32_t> _dynamic_array_synapses_N_outgoing;
std::vector<double> _dynamic_array_synapses_w;

//////////////// dynamic arrays 2d /////////

/////////////// static arrays /////////////
int32_t * _static_array__dynamic_array_spikegeneratorgroup__timebins;
const int _num__static_array__dynamic_array_spikegeneratorgroup__timebins = 155437;
double * _static_array__dynamic_array_spikegeneratorgroup_neuron_index;
const int _num__static_array__dynamic_array_spikegeneratorgroup_neuron_index = 155437;
int64_t * _static_array__dynamic_array_spikegeneratorgroup_spike_number;
const int _num__static_array__dynamic_array_spikegeneratorgroup_spike_number = 155437;
double * _static_array__dynamic_array_spikegeneratorgroup_spike_time;
const int _num__static_array__dynamic_array_spikegeneratorgroup_spike_time = 155437;

//////////////// synapses /////////////////
// synapses
SynapticPathway synapses_post(
    _dynamic_array_synapses__synaptic_post,
    0, 8000);
SynapticPathway synapses_pre(
    _dynamic_array_synapses__synaptic_pre,
    0, 8000);

//////////////// clocks ///////////////////
Clock defaultclock;  // attributes will be set in run.cpp

// Profiling information for each code object
}

void _init_arrays()
{
    using namespace brian;

    // Arrays initialized to 0
    _array_defaultclock_dt = new double[1];
    
    for(int i=0; i<1; i++) _array_defaultclock_dt[i] = 0;

    _array_defaultclock_t = new double[1];
    
    for(int i=0; i<1; i++) _array_defaultclock_t[i] = 0;

    _array_defaultclock_timestep = new int64_t[1];
    
    for(int i=0; i<1; i++) _array_defaultclock_timestep[i] = 0;

    _array_spikegeneratorgroup__lastindex = new int32_t[1];
    
    for(int i=0; i<1; i++) _array_spikegeneratorgroup__lastindex[i] = 0;

    _array_spikegeneratorgroup__period_bins = new int32_t[1];
    
    for(int i=0; i<1; i++) _array_spikegeneratorgroup__period_bins[i] = 0;

    _array_spikegeneratorgroup__spikespace = new int32_t[8001];
    
    for(int i=0; i<8001; i++) _array_spikegeneratorgroup__spikespace[i] = 0;

    _array_spikegeneratorgroup_i = new int32_t[8000];
    
    for(int i=0; i<8000; i++) _array_spikegeneratorgroup_i[i] = 0;

    _array_spikegeneratorgroup_period = new double[1];
    
    for(int i=0; i<1; i++) _array_spikegeneratorgroup_period[i] = 0;

    _array_synapses_N = new int32_t[1];
    
    for(int i=0; i<1; i++) _array_synapses_N[i] = 0;

    _dynamic_array_spikegeneratorgroup__timebins.resize(155437);
    
    for(int i=0; i<155437; i++) _dynamic_array_spikegeneratorgroup__timebins[i] = 0;


    // Arrays initialized to an "arange"
    _array_spikegeneratorgroup_i = new int32_t[8000];
    
    for(int i=0; i<8000; i++) _array_spikegeneratorgroup_i[i] = 0 + i;


    // static arrays
    _static_array__dynamic_array_spikegeneratorgroup__timebins = new int32_t[155437];
    _static_array__dynamic_array_spikegeneratorgroup_neuron_index = new double[155437];
    _static_array__dynamic_array_spikegeneratorgroup_spike_number = new int64_t[155437];
    _static_array__dynamic_array_spikegeneratorgroup_spike_time = new double[155437];

    // Random number generator states
    std::random_device rd;
    for (int i=0; i<1; i++)
        _random_generators.push_back(RandomGenerator());
}

void _load_arrays()
{
    using namespace brian;

    ifstream f_static_array__dynamic_array_spikegeneratorgroup__timebins;
    f_static_array__dynamic_array_spikegeneratorgroup__timebins.open("static_arrays/_static_array__dynamic_array_spikegeneratorgroup__timebins", ios::in | ios::binary);
    if(f_static_array__dynamic_array_spikegeneratorgroup__timebins.is_open())
    {
        f_static_array__dynamic_array_spikegeneratorgroup__timebins.read(reinterpret_cast<char*>(_static_array__dynamic_array_spikegeneratorgroup__timebins), 155437*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_spikegeneratorgroup__timebins." << endl;
    }
    ifstream f_static_array__dynamic_array_spikegeneratorgroup_neuron_index;
    f_static_array__dynamic_array_spikegeneratorgroup_neuron_index.open("static_arrays/_static_array__dynamic_array_spikegeneratorgroup_neuron_index", ios::in | ios::binary);
    if(f_static_array__dynamic_array_spikegeneratorgroup_neuron_index.is_open())
    {
        f_static_array__dynamic_array_spikegeneratorgroup_neuron_index.read(reinterpret_cast<char*>(_static_array__dynamic_array_spikegeneratorgroup_neuron_index), 155437*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_spikegeneratorgroup_neuron_index." << endl;
    }
    ifstream f_static_array__dynamic_array_spikegeneratorgroup_spike_number;
    f_static_array__dynamic_array_spikegeneratorgroup_spike_number.open("static_arrays/_static_array__dynamic_array_spikegeneratorgroup_spike_number", ios::in | ios::binary);
    if(f_static_array__dynamic_array_spikegeneratorgroup_spike_number.is_open())
    {
        f_static_array__dynamic_array_spikegeneratorgroup_spike_number.read(reinterpret_cast<char*>(_static_array__dynamic_array_spikegeneratorgroup_spike_number), 155437*sizeof(int64_t));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_spikegeneratorgroup_spike_number." << endl;
    }
    ifstream f_static_array__dynamic_array_spikegeneratorgroup_spike_time;
    f_static_array__dynamic_array_spikegeneratorgroup_spike_time.open("static_arrays/_static_array__dynamic_array_spikegeneratorgroup_spike_time", ios::in | ios::binary);
    if(f_static_array__dynamic_array_spikegeneratorgroup_spike_time.is_open())
    {
        f_static_array__dynamic_array_spikegeneratorgroup_spike_time.read(reinterpret_cast<char*>(_static_array__dynamic_array_spikegeneratorgroup_spike_time), 155437*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_spikegeneratorgroup_spike_time." << endl;
    }
}

void _write_arrays()
{
    using namespace brian;

    ofstream outfile__array_defaultclock_dt;
    outfile__array_defaultclock_dt.open(results_dir + "_array_defaultclock_dt_1978099143", ios::binary | ios::out);
    if(outfile__array_defaultclock_dt.is_open())
    {
        outfile__array_defaultclock_dt.write(reinterpret_cast<char*>(_array_defaultclock_dt), 1*sizeof(_array_defaultclock_dt[0]));
        outfile__array_defaultclock_dt.close();
    } else
    {
        std::cout << "Error writing output file for _array_defaultclock_dt." << endl;
    }
    ofstream outfile__array_defaultclock_t;
    outfile__array_defaultclock_t.open(results_dir + "_array_defaultclock_t_2669362164", ios::binary | ios::out);
    if(outfile__array_defaultclock_t.is_open())
    {
        outfile__array_defaultclock_t.write(reinterpret_cast<char*>(_array_defaultclock_t), 1*sizeof(_array_defaultclock_t[0]));
        outfile__array_defaultclock_t.close();
    } else
    {
        std::cout << "Error writing output file for _array_defaultclock_t." << endl;
    }
    ofstream outfile__array_defaultclock_timestep;
    outfile__array_defaultclock_timestep.open(results_dir + "_array_defaultclock_timestep_144223508", ios::binary | ios::out);
    if(outfile__array_defaultclock_timestep.is_open())
    {
        outfile__array_defaultclock_timestep.write(reinterpret_cast<char*>(_array_defaultclock_timestep), 1*sizeof(_array_defaultclock_timestep[0]));
        outfile__array_defaultclock_timestep.close();
    } else
    {
        std::cout << "Error writing output file for _array_defaultclock_timestep." << endl;
    }
    ofstream outfile__array_spikegeneratorgroup__lastindex;
    outfile__array_spikegeneratorgroup__lastindex.open(results_dir + "_array_spikegeneratorgroup__lastindex_987837788", ios::binary | ios::out);
    if(outfile__array_spikegeneratorgroup__lastindex.is_open())
    {
        outfile__array_spikegeneratorgroup__lastindex.write(reinterpret_cast<char*>(_array_spikegeneratorgroup__lastindex), 1*sizeof(_array_spikegeneratorgroup__lastindex[0]));
        outfile__array_spikegeneratorgroup__lastindex.close();
    } else
    {
        std::cout << "Error writing output file for _array_spikegeneratorgroup__lastindex." << endl;
    }
    ofstream outfile__array_spikegeneratorgroup__period_bins;
    outfile__array_spikegeneratorgroup__period_bins.open(results_dir + "_array_spikegeneratorgroup__period_bins_4200411184", ios::binary | ios::out);
    if(outfile__array_spikegeneratorgroup__period_bins.is_open())
    {
        outfile__array_spikegeneratorgroup__period_bins.write(reinterpret_cast<char*>(_array_spikegeneratorgroup__period_bins), 1*sizeof(_array_spikegeneratorgroup__period_bins[0]));
        outfile__array_spikegeneratorgroup__period_bins.close();
    } else
    {
        std::cout << "Error writing output file for _array_spikegeneratorgroup__period_bins." << endl;
    }
    ofstream outfile__array_spikegeneratorgroup__spikespace;
    outfile__array_spikegeneratorgroup__spikespace.open(results_dir + "_array_spikegeneratorgroup__spikespace_37288933", ios::binary | ios::out);
    if(outfile__array_spikegeneratorgroup__spikespace.is_open())
    {
        outfile__array_spikegeneratorgroup__spikespace.write(reinterpret_cast<char*>(_array_spikegeneratorgroup__spikespace), 8001*sizeof(_array_spikegeneratorgroup__spikespace[0]));
        outfile__array_spikegeneratorgroup__spikespace.close();
    } else
    {
        std::cout << "Error writing output file for _array_spikegeneratorgroup__spikespace." << endl;
    }
    ofstream outfile__array_spikegeneratorgroup_i;
    outfile__array_spikegeneratorgroup_i.open(results_dir + "_array_spikegeneratorgroup_i_1329498599", ios::binary | ios::out);
    if(outfile__array_spikegeneratorgroup_i.is_open())
    {
        outfile__array_spikegeneratorgroup_i.write(reinterpret_cast<char*>(_array_spikegeneratorgroup_i), 8000*sizeof(_array_spikegeneratorgroup_i[0]));
        outfile__array_spikegeneratorgroup_i.close();
    } else
    {
        std::cout << "Error writing output file for _array_spikegeneratorgroup_i." << endl;
    }
    ofstream outfile__array_spikegeneratorgroup_period;
    outfile__array_spikegeneratorgroup_period.open(results_dir + "_array_spikegeneratorgroup_period_3457314764", ios::binary | ios::out);
    if(outfile__array_spikegeneratorgroup_period.is_open())
    {
        outfile__array_spikegeneratorgroup_period.write(reinterpret_cast<char*>(_array_spikegeneratorgroup_period), 1*sizeof(_array_spikegeneratorgroup_period[0]));
        outfile__array_spikegeneratorgroup_period.close();
    } else
    {
        std::cout << "Error writing output file for _array_spikegeneratorgroup_period." << endl;
    }
    ofstream outfile__array_synapses_N;
    outfile__array_synapses_N.open(results_dir + "_array_synapses_N_483293785", ios::binary | ios::out);
    if(outfile__array_synapses_N.is_open())
    {
        outfile__array_synapses_N.write(reinterpret_cast<char*>(_array_synapses_N), 1*sizeof(_array_synapses_N[0]));
        outfile__array_synapses_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_synapses_N." << endl;
    }

    ofstream outfile__dynamic_array_spikegeneratorgroup__timebins;
    outfile__dynamic_array_spikegeneratorgroup__timebins.open(results_dir + "_dynamic_array_spikegeneratorgroup__timebins_4247931087", ios::binary | ios::out);
    if(outfile__dynamic_array_spikegeneratorgroup__timebins.is_open())
    {
        if (! _dynamic_array_spikegeneratorgroup__timebins.empty() )
        {
            outfile__dynamic_array_spikegeneratorgroup__timebins.write(reinterpret_cast<char*>(&_dynamic_array_spikegeneratorgroup__timebins[0]), _dynamic_array_spikegeneratorgroup__timebins.size()*sizeof(_dynamic_array_spikegeneratorgroup__timebins[0]));
            outfile__dynamic_array_spikegeneratorgroup__timebins.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_spikegeneratorgroup__timebins." << endl;
    }
    ofstream outfile__dynamic_array_spikegeneratorgroup_neuron_index;
    outfile__dynamic_array_spikegeneratorgroup_neuron_index.open(results_dir + "_dynamic_array_spikegeneratorgroup_neuron_index_2789266935", ios::binary | ios::out);
    if(outfile__dynamic_array_spikegeneratorgroup_neuron_index.is_open())
    {
        if (! _dynamic_array_spikegeneratorgroup_neuron_index.empty() )
        {
            outfile__dynamic_array_spikegeneratorgroup_neuron_index.write(reinterpret_cast<char*>(&_dynamic_array_spikegeneratorgroup_neuron_index[0]), _dynamic_array_spikegeneratorgroup_neuron_index.size()*sizeof(_dynamic_array_spikegeneratorgroup_neuron_index[0]));
            outfile__dynamic_array_spikegeneratorgroup_neuron_index.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_spikegeneratorgroup_neuron_index." << endl;
    }
    ofstream outfile__dynamic_array_spikegeneratorgroup_spike_number;
    outfile__dynamic_array_spikegeneratorgroup_spike_number.open(results_dir + "_dynamic_array_spikegeneratorgroup_spike_number_3584615111", ios::binary | ios::out);
    if(outfile__dynamic_array_spikegeneratorgroup_spike_number.is_open())
    {
        if (! _dynamic_array_spikegeneratorgroup_spike_number.empty() )
        {
            outfile__dynamic_array_spikegeneratorgroup_spike_number.write(reinterpret_cast<char*>(&_dynamic_array_spikegeneratorgroup_spike_number[0]), _dynamic_array_spikegeneratorgroup_spike_number.size()*sizeof(_dynamic_array_spikegeneratorgroup_spike_number[0]));
            outfile__dynamic_array_spikegeneratorgroup_spike_number.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_spikegeneratorgroup_spike_number." << endl;
    }
    ofstream outfile__dynamic_array_spikegeneratorgroup_spike_time;
    outfile__dynamic_array_spikegeneratorgroup_spike_time.open(results_dir + "_dynamic_array_spikegeneratorgroup_spike_time_2775435892", ios::binary | ios::out);
    if(outfile__dynamic_array_spikegeneratorgroup_spike_time.is_open())
    {
        if (! _dynamic_array_spikegeneratorgroup_spike_time.empty() )
        {
            outfile__dynamic_array_spikegeneratorgroup_spike_time.write(reinterpret_cast<char*>(&_dynamic_array_spikegeneratorgroup_spike_time[0]), _dynamic_array_spikegeneratorgroup_spike_time.size()*sizeof(_dynamic_array_spikegeneratorgroup_spike_time[0]));
            outfile__dynamic_array_spikegeneratorgroup_spike_time.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_spikegeneratorgroup_spike_time." << endl;
    }
    ofstream outfile__dynamic_array_synapses__synaptic_post;
    outfile__dynamic_array_synapses__synaptic_post.open(results_dir + "_dynamic_array_synapses__synaptic_post_1801389495", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses__synaptic_post.is_open())
    {
        if (! _dynamic_array_synapses__synaptic_post.empty() )
        {
            outfile__dynamic_array_synapses__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_synapses__synaptic_post[0]), _dynamic_array_synapses__synaptic_post.size()*sizeof(_dynamic_array_synapses__synaptic_post[0]));
            outfile__dynamic_array_synapses__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_synapses__synaptic_pre;
    outfile__dynamic_array_synapses__synaptic_pre.open(results_dir + "_dynamic_array_synapses__synaptic_pre_814148175", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses__synaptic_pre.is_open())
    {
        if (! _dynamic_array_synapses__synaptic_pre.empty() )
        {
            outfile__dynamic_array_synapses__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_synapses__synaptic_pre[0]), _dynamic_array_synapses__synaptic_pre.size()*sizeof(_dynamic_array_synapses__synaptic_pre[0]));
            outfile__dynamic_array_synapses__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_synapses_A_postsyn;
    outfile__dynamic_array_synapses_A_postsyn.open(results_dir + "_dynamic_array_synapses_A_postsyn_4029683827", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_A_postsyn.is_open())
    {
        if (! _dynamic_array_synapses_A_postsyn.empty() )
        {
            outfile__dynamic_array_synapses_A_postsyn.write(reinterpret_cast<char*>(&_dynamic_array_synapses_A_postsyn[0]), _dynamic_array_synapses_A_postsyn.size()*sizeof(_dynamic_array_synapses_A_postsyn[0]));
            outfile__dynamic_array_synapses_A_postsyn.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_A_postsyn." << endl;
    }
    ofstream outfile__dynamic_array_synapses_A_presyn;
    outfile__dynamic_array_synapses_A_presyn.open(results_dir + "_dynamic_array_synapses_A_presyn_2710595424", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_A_presyn.is_open())
    {
        if (! _dynamic_array_synapses_A_presyn.empty() )
        {
            outfile__dynamic_array_synapses_A_presyn.write(reinterpret_cast<char*>(&_dynamic_array_synapses_A_presyn[0]), _dynamic_array_synapses_A_presyn.size()*sizeof(_dynamic_array_synapses_A_presyn[0]));
            outfile__dynamic_array_synapses_A_presyn.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_A_presyn." << endl;
    }
    ofstream outfile__dynamic_array_synapses_delay;
    outfile__dynamic_array_synapses_delay.open(results_dir + "_dynamic_array_synapses_delay_3246960869", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_delay.is_open())
    {
        if (! _dynamic_array_synapses_delay.empty() )
        {
            outfile__dynamic_array_synapses_delay.write(reinterpret_cast<char*>(&_dynamic_array_synapses_delay[0]), _dynamic_array_synapses_delay.size()*sizeof(_dynamic_array_synapses_delay[0]));
            outfile__dynamic_array_synapses_delay.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_delay." << endl;
    }
    ofstream outfile__dynamic_array_synapses_delay_1;
    outfile__dynamic_array_synapses_delay_1.open(results_dir + "_dynamic_array_synapses_delay_1_3310102259", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_delay_1.is_open())
    {
        if (! _dynamic_array_synapses_delay_1.empty() )
        {
            outfile__dynamic_array_synapses_delay_1.write(reinterpret_cast<char*>(&_dynamic_array_synapses_delay_1[0]), _dynamic_array_synapses_delay_1.size()*sizeof(_dynamic_array_synapses_delay_1[0]));
            outfile__dynamic_array_synapses_delay_1.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_delay_1." << endl;
    }
    ofstream outfile__dynamic_array_synapses_lastupdate;
    outfile__dynamic_array_synapses_lastupdate.open(results_dir + "_dynamic_array_synapses_lastupdate_3710850267", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_lastupdate.is_open())
    {
        if (! _dynamic_array_synapses_lastupdate.empty() )
        {
            outfile__dynamic_array_synapses_lastupdate.write(reinterpret_cast<char*>(&_dynamic_array_synapses_lastupdate[0]), _dynamic_array_synapses_lastupdate.size()*sizeof(_dynamic_array_synapses_lastupdate[0]));
            outfile__dynamic_array_synapses_lastupdate.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_lastupdate." << endl;
    }
    ofstream outfile__dynamic_array_synapses_N_incoming;
    outfile__dynamic_array_synapses_N_incoming.open(results_dir + "_dynamic_array_synapses_N_incoming_1151751685", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_N_incoming.is_open())
    {
        if (! _dynamic_array_synapses_N_incoming.empty() )
        {
            outfile__dynamic_array_synapses_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_synapses_N_incoming[0]), _dynamic_array_synapses_N_incoming.size()*sizeof(_dynamic_array_synapses_N_incoming[0]));
            outfile__dynamic_array_synapses_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_synapses_N_outgoing;
    outfile__dynamic_array_synapses_N_outgoing.open(results_dir + "_dynamic_array_synapses_N_outgoing_1673144031", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_N_outgoing.is_open())
    {
        if (! _dynamic_array_synapses_N_outgoing.empty() )
        {
            outfile__dynamic_array_synapses_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_synapses_N_outgoing[0]), _dynamic_array_synapses_N_outgoing.size()*sizeof(_dynamic_array_synapses_N_outgoing[0]));
            outfile__dynamic_array_synapses_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_synapses_w;
    outfile__dynamic_array_synapses_w.open(results_dir + "_dynamic_array_synapses_w_441891901", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_w.is_open())
    {
        if (! _dynamic_array_synapses_w.empty() )
        {
            outfile__dynamic_array_synapses_w.write(reinterpret_cast<char*>(&_dynamic_array_synapses_w[0]), _dynamic_array_synapses_w.size()*sizeof(_dynamic_array_synapses_w[0]));
            outfile__dynamic_array_synapses_w.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_w." << endl;
    }

    // Write last run info to disk
    ofstream outfile_last_run_info;
    outfile_last_run_info.open(results_dir + "last_run_info.txt", ios::out);
    if(outfile_last_run_info.is_open())
    {
        outfile_last_run_info << (Network::_last_run_time) << " " << (Network::_last_run_completed_fraction) << std::endl;
        outfile_last_run_info.close();
    } else
    {
        std::cout << "Error writing last run info to file." << std::endl;
    }
}

void _dealloc_arrays()
{
    using namespace brian;


    // static arrays
    if(_static_array__dynamic_array_spikegeneratorgroup__timebins!=0)
    {
        delete [] _static_array__dynamic_array_spikegeneratorgroup__timebins;
        _static_array__dynamic_array_spikegeneratorgroup__timebins = 0;
    }
    if(_static_array__dynamic_array_spikegeneratorgroup_neuron_index!=0)
    {
        delete [] _static_array__dynamic_array_spikegeneratorgroup_neuron_index;
        _static_array__dynamic_array_spikegeneratorgroup_neuron_index = 0;
    }
    if(_static_array__dynamic_array_spikegeneratorgroup_spike_number!=0)
    {
        delete [] _static_array__dynamic_array_spikegeneratorgroup_spike_number;
        _static_array__dynamic_array_spikegeneratorgroup_spike_number = 0;
    }
    if(_static_array__dynamic_array_spikegeneratorgroup_spike_time!=0)
    {
        delete [] _static_array__dynamic_array_spikegeneratorgroup_spike_time;
        _static_array__dynamic_array_spikegeneratorgroup_spike_time = 0;
    }
}

