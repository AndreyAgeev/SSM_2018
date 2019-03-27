#include "functions-grammar.h"

using namespace std;
using namespace HighFive;

class Nlreg
{
public:
	Nlreg(QString func_file_name, std::vector<std::string> meas, std::vector<std::string> grn, int nF, int wl, int rF, int print_trace = 0) : funcs_file_name(func_file_name), measurements(meas), gr_names(grn) {
		/*		std::cout << "Got " << h5_file_name.toStdString() << " file" << std::endl;
				std::cout << "Got " << funcs_file_name.toStdString() << " funcs" << std::endl;*/
		nFunctions = nF;
		wordLength = wl;
		num_of_climate_vars = measurements.size();
		num_of_gt_vars = gr_names.size();
		PRINT_TRACE = print_trace;
		read_flag = rF;
		read_genotype();
	}
	double get_func_value(vector<double> clim_arg, vector<double> gt_vars)
	{
		double val = 0;
		GrammarNode *retFtn;
		for (size_t i = 0; i < nFunctions; ++i) {
			retFtn = grc->get_nth_tree(i);
			double fval = retFtn->eval(clim_arg);
			val += beta[i] * fval;
			for (size_t j = 0; j < num_of_gt_vars; ++j) {
				val += beta[j + nFunctions] * gt_vars[j] * fval;
			}
		}
		return val;
	}
	double get_cbd()
	{
		return CBD;
	}
	double get_l1_pen()
	{
		double val = 0;
		for (size_t i = 0; i < nFunctions; ++i) {
			val += (beta[i] > 0) ? beta[i] : -beta[i];
			for (size_t j = 0; j < num_of_gt_vars; ++j) {
				val += (beta[j + nFunctions] > 0) ? beta[j + nFunctions] : -beta[j + nFunctions];
			}
		}
		return val;
	}
	double get_l2_pen()
	{
		double val = 0;
		for (size_t i = 0; i < nFunctions; ++i) {
			val += beta[i] * beta[i];
			for (size_t j = 0; j < num_of_gt_vars; ++j) {
				val += beta[j + nFunctions] * beta[j + nFunctions];
			}
		}
		return val;
	}
	void nlreg_build()
	{
		GrammarNode *retFtn;
		grc = new GrammarContainer(measurements, nFunctions);
		phenotype = new double[nFunctions * wordLength];
		phenomask = new int[nFunctions * wordLength];
		int first = 0, last = wordLength;
		for (int i = 0; i < nFunctions; ++i) {
			std::vector<int> gt = std::vector<int>(genotype.begin() + first, genotype.begin() + last);
			first += wordLength;
			last += wordLength;
			grc->build_nth_tree(gt, climate_var, i, &phenotype[i * wordLength], &phenomask[i * wordLength]);
		}
	}
	void print_trace(int arg) {
		GrammarNode *retFtn;
		size_t i, j;
		if (arg == 0) {
			for (i = 0; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
				cout << beta[i] << " ";
			}
			cout << endl;
			cout << "F= ";
			int flag = 0;
			for (i = 0; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
				if (fabs(beta[i]) > 0) {
					if (flag == 1 && beta[i] > 0) cout << "+";
					cout << beta[i];
					j = i;
					if (j < nFunctions) {
						cout << "*";
						retFtn = grc->get_nth_tree(j);
						retFtn->coprint();
					}
					else if (j >= nFunctions) {
						int jj = j - nFunctions;
						int ii = jj / num_of_gt_vars;
						int kk = jj % num_of_gt_vars;
						cout << "*" << gr_names[kk] << "*" << "f[" << ii << "]";
					}
					flag = 1;
				}
			}
			cout << endl;
			for (size_t i = 0; i < nFunctions * wordLength; ++i) { // + nEcovar
				cout << phenotype[i] << " ";
			}
			cout << endl;
			for (size_t i = 0; i < nFunctions * wordLength; ++i) { // + nEcovar
				cout << phenomask[i] << " ";
			}
			cout << endl;
		}
		else {
			for (size_t i = 0; i < nFunctions; ++i) {
				retFtn = grc->get_nth_tree(i);
				cout << i << " ";
				retFtn->coprint();
				cout << endl;
			}
		}
	}
private:
	int num_of_climate_vars; // Number of columns in the weather table == number of consts to read after func's and betas
	int num_of_gt_vars; // Number of genotype data columns (snp or qtl or locations). Number of betas=number of funcs + number of func * number of gt vars
	int nFunctions;
	int wordLength;
	int PRINT_TRACE = 1;
	int read_flag = -1;
	std::vector<std::string> measurements;
	std::vector<std::string> gr_names;
	GrammarContainer* grc;
	vector<double> climate_var;
	vector<double> beta;
	double CBD;
	double MB;
	QString funcs_file_name;
	double *phenotype;
	int *phenomask;
	vector<int> genotype;
	void read_genotype() {
		QFile file(funcs_file_name);
		vector<int> gt;
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
			cout << "Can't open " << funcs_file_name.toStdString() << endl;
			return;
		}
		QTextStream in(&file);
		QString nnm;
		in >> nnm;
		for (size_t i = 0; i < nFunctions * wordLength; ++i) {
			int arg;
			in >> arg;
			gt.push_back(arg);
			//							cout << gt[i] << " " << arg << endl;
		}
		if (read_flag > -1) {
			vector<double> be;
			for (size_t i = 0; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
				double arg;
				in >> arg;
				be.push_back(arg);
				//			cout << gt[i] << endl;
			}
			beta = be;
		}
		vector<double> concs;
		for (size_t i = 0; i < num_of_climate_vars; ++i) {
			double arg;
			in >> arg;
			concs.push_back(arg);
			//			cout << gt[i] << endl;
		}
		climate_var = concs;
		if (read_flag > 0) {
			in >> CBD;
			in >> MB;
			for (size_t i = nFunctions; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
				double be = (beta[i] > 0) ? beta[i] : -beta[i];
				beta[i] = (be < MB) ? 0.0 : beta[i];
			}
		}
		//					for (size_t i = 0; i < nFunctions * wordLength; ++i) { cout << gt[i] << endl; }
		genotype = gt;
	};
};



