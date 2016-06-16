#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<sys/time.h>
#include<algorithm>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "../src/Stopwatch.h"
#include "../src/parameters.h"
#include "../src/readindata.h"
#include "../src/fitCF.h"
#include "../src/generate_processing_record.h"
#include "../src/plumberglib.h"
#include "../src/sorter.h"
#include "fit_correlation_function.h"

using namespace std;

int main(int argc, char *argv[])
{
	Stopwatch sw;
	Stopwatch sw_total;
	sw_total.tic();
	sw.tic();

	//check if proper number of command-line arguments were passed
	if (argc != 1 && argc != 3)
	{
		cerr << "You passed the wrong number of parameters to fit_correlation_function()!  Exiting..." << endl;
		exit(1);
	}

	bool generatedcorrfuncs = false;
	string currentworkingdirectory = get_selfpath();

	initialize_PRfile(currentworkingdirectory);

	ostringstream filename_stream;
	filename_stream << currentworkingdirectory << "/Processing_record.txt";
	ofstream output(filename_stream.str().c_str(), ios::app);

	output << "/**********Processing output**********/" << endl;
	output << "entering folder: " << currentworkingdirectory << endl;

	//load freeze out and particle information
	int particle_idx = 1;  //for pion+

	//read particle resonance decay table
	particle_info * particle = new particle_info [Maxparticle];
	int Nparticle = read_resonance(particle);
	output << "Read in total " << Nparticle << " particles!" << endl;
	output << "Calculating "<< particle[particle_idx].name << endl;

	//HBT calculations begin ...
	double localy = 0.0e0;
	sw.tic();
	if(fabs(localy) > 1e-16)
	{
		output << "Case of y != 0 not yet supported.  Exiting..." << endl;
		return 0;
	}

	//defines events to include in ensemble average
	//handles single event automatically
	vector<int> chosen_events;
	int lower_limit = -1;
	int upper_limit = -1;
	if (argc == 1)
	{
		int folderindex = get_folder_index(currentworkingdirectory);
		lower_limit = folderindex;
		upper_limit = folderindex;
	}
	else if (argc == 3)
	{
		lower_limit = atoi(argv[1]);
		upper_limit = atoi(argv[2]);
	}
	for (int i = lower_limit; i <= upper_limit; ++i)	//set chosen_events from command line
		chosen_events.push_back(i);

	//for (int i = 0; i < (int)chosen_events.size(); ++i)
	//	cout << "chosen_events[" << i << "] = " << chosen_events[i] << endl;

	//cout << "chosen_events.size() = " << chosen_events.size() << endl;

	FitCF fit_CF(particle, Nparticle, particle_idx, chosen_events, output);

	//fit_CF.read_in_all_dN_dypTdpTdphi = false;
	//fit_CF.output_all_dN_dypTdpTdphi = !(fit_CF.read_in_all_dN_dypTdpTdphi);
	fit_CF.Set_path(currentworkingdirectory);
	fit_CF.Set_use_delta_f(true);
	//fit_CF.Set_ofstream(output);

	if (argc == 1)
		fit_CF.currentfolderindex = get_folder_index(currentworkingdirectory);
		//if no CMD args passed, assume we're actually in the folder

	////////////////////////////////////////////
	// Actual calculations start here...
	////////////////////////////////////////////

	output << "Calculating correlation function with all resonance decays..." << endl;
	//do calculations
	fit_CF.Cal_correlationfunction();	//if we didn't compute resonance decays, must read them in from files

	//output results
	//fit_CF.Output_correlationfunction();

	output << "Finished calculating correlation function with all resonance decays..." << endl;

	//if there's a full 3D grid to fit over, do the Gaussian fit and get the HBT radii too
	if (qxnpts > 1 && qynpts > 1 && qznpts > 1)
	{
		output << "Calculating HBT radii via Gaussian fit method..." << endl;
		fit_CF.Get_GF_HBTradii();
		fit_CF.Output_results(0);	//0 means do GF R2ij
		output << "Finished calculating HBT radii via Gaussian fit method" << endl;
	}

	sw.toc();
	output << "Finished in " << sw.takeTime() << " sec." << endl;
	sw_total.toc();
	output << "Program totally finished in " << sw_total.takeTime() << " sec." << endl;

	output << "/**********End of processing output**********/" << endl;

	output.close();

	//checkforfiles_PRfile(currentworkingdirectory, folderindex, generatedcorrfuncs);

	finalize_PRfile(currentworkingdirectory);

	return 0;
}

//End of file
