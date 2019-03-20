#pragma once
//#include "model.h"
//#include "parametrs.h"

//#include <highfive/H5Attribute.hpp>
//#include <highfive/H5File.hpp>
//#include <highfive/H5DataSet.hpp>
//#include <highfive/H5DataSpace.hpp>

#include "functions-grammar.h"
//#include "larsw.h" -- This is not needed here as we optimize globaly

using namespace std;
using namespace HighFive;

class Nlreg
{
    public:
	Nlreg(QString func_file_name, std::vector<std::string> meas, int nF, int wl, int nC, int nG, int print_trace = 0) : funcs_file_name(func_file_name), measurements(meas) {
/*		std::cout << "Got " << h5_file_name.toStdString() << " file" << std::endl;
		std::cout << "Got " << funcs_file_name.toStdString() << " funcs" << std::endl;*/
		nFunctions = nF;
		wordlength = wl;
		num_of_climate_vars = nC;
		num_of_gt_vars = nG;
		PRINT_TRACE = print_trace;
	}
    private:
		int num_of_climate_vars; // Number of columns in the weather table == number of consts to read after func's and betas
		int num_of_gt_vars; // Number of genotype data columns (snp or qtl or locations). Number of betas=number of funcs + number of func * number of gt vars
		int nFunctions;
		int wordLength;
		int PRINT_TRACE = 1;
		std::vector<std::string> measurements;
		GrammarContainer* grc;
		vector<double> climate_var;
		vector<double> beta;
		QString funcs_file_name;
		double *phenotype;
		int *phenomask;
		vector<int> read_genotype(vector<double>& climate_var, vector<double>& beta) {
			QFile file(funcs_file_name);
			vector<int> gt;
			if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
				cout << "Can't open " << funcs_file_name.toStdString() << endl;
				return gt;
			}
			QTextStream in(&file);
/* Alway true - remove option
 * 			if (_param.nl_param.force) {*/
				QString nnm;
				in >> nnm;
				/*			cout << nnm.toStdString () << endl;*/
/*			}
*/
			for (size_t i = 0; i < _param.nl_param.nFunctions * _param.nl_param.wordLength; ++i) {
				int arg;
				in >> arg;
				gt.push_back(arg);
				//			cout << gt[i] << endl;
			}
			vector<double> be;
			for (size_t i = 0; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
				double arg;
				in >> arg;
				be.push_back(arg);
				//			cout << gt[i] << endl;
			}
			beta = concs;
			vector<double> concs;
			for (size_t i = 0; i < num_of_climate_vars; ++i) {
				double arg;
				in >> arg;
				concs.push_back(arg);
				//			cout << gt[i] << endl;
			}
			climate_var = concs;
			//		for (size_t i = 0; i < nFunctions * wordlength; ++i) { cout << gt[i] << endl; }
			return gt;
		};
		double get_func_value(vector<double> clim_arg)
		{
			double val = 0;
			GrammarNode *retFtn;
			for (size_t i = 0; i < nFunctions; ++i) {
				retFtn = grc->get_nth_tree(i);
				double fval = retFtn->eval(clim_arg);
				val += beta[i] * fval
				for (size_t j = 0; j < num_of_gt_vars; ++j) {
					val += beta[j + nFunctions] * fval;
				}
			}
			return val;
		}
		void nlreg_build()
		{
			vector<int> genotype = read_genotype(climate_var, beta);
			grc = new GrammarContainer(measurements, nFunctions);
			phenotype = new double[nFunctions * (wordlength + 1) + nFunctions * num_of_gt_vars];
			phenomask = new int[nFunctions * (wordlength + 1) + nFunctions * num_of_gt_vars];
			int first = 0, last = wordLength;
			for (int i = 0; i < nFunctions; ++i) {
				std::vector<int> gt = std::vector<int>(genotype.begin() + first, genotype.begin() + last);
				first += wordLength;
				last += wordLength;
				grc->build_nth_tree(gt, climate_var, i, &phenotype[i * wordLength], &phenomask[i * wordlength]);
				if (PRINT_TRACE > 0) {
					for (size_t i = 0; i < nFunctions; ++i) {
						retFtn = grc->get_nth_tree(i);
						cout << i << " ";
						retFtn->coprint();
					}
				}
			}
		}
};


