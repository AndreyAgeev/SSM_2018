#include "functions-grammar.h"
#include "data.h"
using namespace std;
using namespace HighFive;
#include <QSettings>

class Nlreg
{
public:
	Nlreg(QString func_file_name, std::vector<std::string> meas, std::vector<std::string> grn, int nF, int wl, int rF, int print_trace = 0, int crops = 0) : funcs_file_name(func_file_name), measurements(meas), gr_names(grn) {
		/*		std::cout << "Got " << h5_file_name.toStdString() << " file" << std::endl;
				std::cout << "Got " << funcs_file_name.toStdString() << " funcs" << std::endl;*/
		nFunctions = nF;
		wordLength = wl;
		num_of_climate_vars = measurements.size();
		num_of_gt_vars = gr_names.size();
		PRINT_TRACE = print_trace;
		read_flag = rF;
		read_genotype(crops);
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
	void delete_all()
	{
		delete[] phenotype;
		delete[] phenomask;
	}
/*	double get_func_value(vector<double> clim_arg)
	{
		double val = 0;
		GrammarNode *retFtn;
		for (size_t i = 0; i < nFunctions; ++i) {
			retFtn = grc->get_nth_tree(i);
			double fval = retFtn->eval(clim_arg);
			val += beta[i] * fval;
		}
		return val;
	}*/
	double get_cbd()
	{
		return nlCBD;
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
	//	cout << "measurements "  << measurements.size()  << "nFunction " << nFunctions << "wordLength " << wordLength  <<  "climate_var " << climate_var.size() << endl;
		grc = new GrammarContainer(measurements, nFunctions);
		phenotype = new double[nFunctions * wordLength];
		phenomask = new int[nFunctions * wordLength];
	/*	for (int i = 0; i < nFunctions * wordLength; i++)
		{
			phenotype[i] = 0.0;
			phenomask[i] = 0;
		}*/
		int first = 0, last = wordLength;
		for (int i = 0; i < nFunctions; ++i) {
			std::vector<int> gt = std::vector<int>(genotype.begin() + first, genotype.begin() + last);
			first += wordLength;
			last += wordLength;
		//	cout << "first " << first << " last " << last << endl;
		//	cout << "first " << first << " last " << last << endl;
			grc->build_nth_tree(gt, climate_var, i, &phenotype[i * wordLength], &phenomask[i * wordLength]);
		//	cout << " phenotype["<<i << " * wordLength] = " << phenotype[i * wordLength] << " phenomask[" << i << " * wordLength] " << phenomask[i * wordLength] << endl;
		}
	/*	cout << " BEFORE" << endl;

		for (size_t i = 0; i < 5 * 11; ++i) { // + nEcovar
			cout << " phenotype[" << i << " ] = " << phenotype[i] << endl;
		}
		for (size_t i = 0; i < 5 * 11; ++i) { // + nEcovar
			cout << " phenomask[" << i << " ] " << phenomask[i] << endl;
		}*/
	}
	void print_trace(std::string func_name, int arg) {
		std::ofstream out_func;
		out_func.open(func_name + "_" + "out_func.txt", std::ios::app);
		GrammarNode *retFtn;
		size_t i, j;
		if (arg == 0) {
			for (i = 0; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
				out_func << beta[i] << " ";
			}
			out_func << endl;
			out_func << "F= ";
			int flag = 0;
			for (i = 0; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
				if (fabs(beta[i]) > 0.0) {
					if (flag == 1 && beta[i] > 0) out_func << "+";
					out_func << beta[i];
					j = i;
					if (j < nFunctions) {
						out_func << "*";
						retFtn = grc->get_nth_tree(j);
						retFtn->coprint(out_func);
					}
					else if (j >= nFunctions) {
						int jj = j - nFunctions;
						int ii = jj / num_of_gt_vars;
						int kk = jj % num_of_gt_vars;
						out_func << "*" << gr_names[kk] << "*" << "f[" << ii << "]";
					}
					flag = 1;
				}
			}
			out_func << endl; 
			for (size_t i = 0; i < nFunctions * wordLength; ++i) { // + nEcovar
				out_func << phenotype[i] << " ";
		//		cout << " phenotype[" << i << " ] = " << phenotype[i] << endl;
			}
			out_func << endl;
			for (size_t i = 0; i < nFunctions * wordLength; ++i) { // + nEcovar
				out_func << phenomask[i] << " ";
		//		cout << " phenomask[" << i << " ] " << phenomask[i] << endl;
			}
			out_func << endl;
		}
		else {
			for (size_t i = 0; i < nFunctions; ++i) {
				retFtn = grc->get_nth_tree(i);
				out_func << i << " ";
				retFtn->coprint(out_func);
				out_func << endl;
			}
		}
		out_func << endl;
		out_func.close();
	}
private:
	int num_of_climate_vars; // Number of columns in the weather table == number of consts to read after func's and betas
	int num_of_gt_vars; // Number of genotype data columns (snp or qtl or locations). Number of betas=number of funcs + number of func * number of gt vars
	int nFunctions;
	int wordLength;
	int PRINT_TRACE = 1;
	int read_flag = 2;//-1
	std::vector<std::string> measurements;
	std::vector<std::string> gr_names;
	GrammarContainer* grc;
	vector<double> climate_var;
	vector<double> beta;
	double nlCBD;
	double MB;
	QString funcs_file_name;
	double *phenotype;
    int *phenomask;
	vector<int> genotype;
	void read_genotype(const int crops) {
	//	if (crops == 0)
	//	{
			QFile file(funcs_file_name);
			vector<int> gt;
			if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
				cout << "Can't open " << funcs_file_name.toStdString() << endl;
				return;
			}
			QTextStream in(&file);
			QString nnm;
			in >> nnm;
	//		cout << "NEED =  " << nFunctions * wordLength + nFunctions + nFunctions * num_of_gt_vars + num_of_climate_vars + 2 << endl;
	//		cout << "ge  " << nFunctions * wordLength << endl;
	//		cout << "n+nf*num_of_gt  " << nFunctions  << endl;
	//		cout << "nf*num_of_gt  " << num_of_gt_vars << endl;

			
			for (size_t i = 0; i < nFunctions * wordLength; ++i) {
				int arg1;
				in >> arg1;
				gt.push_back(arg1);
			//									cout << gt[i] << endl;
			}
			if (read_flag > -1) {
				vector<double> be;
				for (size_t i = 0; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
					double arg2;
					in >> arg2;
					be.push_back(arg2);
		//							cout << be[i] << endl;
				}
				beta = be;
			}
			vector<double> concs;
		//	cout << "num_of_climate_vars " << num_of_climate_vars << endl;

			for (size_t i = 0; i < num_of_climate_vars; ++i) {
				double arg3;
				in >> arg3;
					//cout << arg << endl;
				concs.push_back(arg3);
			}
			climate_var = concs;
			if (read_flag > 0) {
				in >> nlCBD;
				in >> MB;
				for (size_t i = nFunctions; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
					double be = (beta[i] > 0.0) ? beta[i] : -beta[i];
					beta[i] = (be < MB) ? 0.0 : beta[i];
				}
			}
		//	for (size_t i = 0; i < nFunctions * wordLength; ++i) { cout << gt[i] << " "; }
			genotype = gt;
			file.close();

	//	}
		/*else
		{
			vector<int> gt;
			QSettings sett(funcs_file_name, QSettings::IniFormat);
			sett.beginGroup("Funcs");
			for (size_t i = 0; i < nFunctions * wordLength; ++i) {
				int arg;
				string codon = "codon";
				string str = std::to_string(i);
				string value = codon + str;
				QString val = QString::fromUtf8(value.c_str());
				arg = sett.value(val, 1).toInt();
				gt.push_back(arg);
			}
			if (read_flag > -1) {
				vector<double> be;
				for (size_t i = 0; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
					int arg;
					string codon = "beta";
					string str = std::to_string(i);
					string value = codon + str;
					QString val = QString::fromUtf8(value.c_str());
					arg = sett.value(val, 46).toDouble();
					be.push_back(arg);
					//					cout << be[i] << endl;
				}
				beta = be;
			}
			vector<double> concs;
			for (size_t i = 0; i < num_of_climate_vars; ++i) {
				int arg;
				string codon = "climate_vars";
				string str = std::to_string(i);
				string value = codon + str;
				QString val = QString::fromUtf8(value.c_str());
				arg = sett.value(val, 46).toDouble();
				//		cout << arg << endl;
				concs.push_back(arg);
			}
			climate_var = concs;
			if (read_flag > 0) {
				CBD = sett.value("CBD", 46).toDouble();
				MB = sett.value("MB", 46).toDouble();
				for (size_t i = nFunctions; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
					double be = (beta[i] > 0) ? beta[i] : -beta[i];
					beta[i] = (be < MB) ? 0.0 : beta[i];
				}
			}
			//					for (size_t i = 0; i < nFunctions * wordLength; ++i) { cout << gt[i] << endl; }
			genotype = gt;
			sett.endGroup();
		}*/

	}
};
