#include "functions-grammar.h"
#include "data.h"
using namespace std;
using namespace HighFive;

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
	Data::Data_phen data_p;
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
	double get_func_value(vector<double> clim_arg)
	{
		double val = 0;
		GrammarNode *retFtn;
		for (size_t i = 0; i < nFunctions; ++i) {
			retFtn = grc->get_nth_tree(i);
			double fval = retFtn->eval(clim_arg);
			val += beta[i] * fval;
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
				if (fabs(beta[i]) > 0) {
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
			}
			out_func << endl;
			for (size_t i = 0; i < nFunctions * wordLength; ++i) { // + nEcovar
				out_func << phenomask[i] << " ";
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
	void read_genotype(const int crops) {
		QFile file(funcs_file_name);
		vector<int> gt;
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
			cout << "Can't open " << funcs_file_name.toStdString() << endl;
			return;
		}
		QTextStream in(&file);
		QString nnm;
		in >> nnm;
		//cout << "NEED =  " << nFunctions * wordLength + nFunctions + nFunctions * num_of_gt_vars + num_of_climate_vars + 2 << endl;
		for (size_t i = 0; i < nFunctions * wordLength; ++i) {
			int arg;
			in >> arg;
			gt.push_back(arg);
		//								cout << gt[i] << endl;
		}
		if (read_flag > -1) {
		//	cout << "NUM_OF_GT" << num_of_gt_vars << endl;
			vector<double> be;
			for (size_t i = 0; i < nFunctions + nFunctions * num_of_gt_vars; ++i) {
				double arg;
				in >> arg;
				be.push_back(arg);
		//					cout << be[i] << endl;
			}
			beta = be;
		}
		vector<double> concs;
		for (size_t i = 0; i < num_of_climate_vars; ++i) {
			double arg;
			in >> arg;
	//		cout << arg << endl;
			concs.push_back(arg);
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
		if (crops == 1)
		{
			double arg1;
			in >> arg1;
			data_p.phyl = arg1;
			double arg2;
			in >> arg2;
			data_p.PLACON = arg2;
			double arg3;
			in >> arg3;
			data_p.PLAPOW30 = arg3;
			double arg4;
			in >> arg4;
			data_p.SLA = arg4;
			double arg5;
			in >> arg5;
			data_p.TBRUE = arg5;
			double arg6;
			in >> arg6;
			data_p.TP1RUE = arg6;
			double arg7;
			in >> arg7;
			data_p.TP2RUE = arg7;
			double arg8;
			in >> arg8;
			data_p.TCRUE = arg8;
			double arg9;
			in >> arg9;
			data_p.KPAR = arg9;
			double arg10;
			in >> arg10;
			data_p.IRUE1 = arg10;
			double arg11;
			in >> arg11;
			data_p.IRUE2 = arg11;
			double arg12;
			in >> arg12;
			data_p.FLF1A = arg12;
			double arg13;
			in >> arg13;
			data_p.FLF1B = arg13;
			double arg14;
			in >> arg14;
			data_p.WTOPL = arg14;
			double arg15;
			in >> arg15;
			data_p.FLF2 = arg15;
			double arg16;
			in >> arg16;
			data_p.FRTRL = arg16;
			double arg17;
			in >> arg17;
			data_p.GCF = arg17;
			double arg18;
			in >> arg18;
			data_p.PDHI = arg18;
			double arg19;
			in >> arg19;
			data_p.WDHI1 = arg19;
			double arg20;
			in >> arg20;
			data_p.WDHI2 = arg20;
			double arg21;
			in >> arg21;
			data_p.WDHI3 = arg21;
			double arg22;
			in >> arg22;
			data_p.WDHI4 = arg22;
			double arg23;
			in >> arg23;
			data_p.DEPORT = arg23;
			double arg24;
			in >> arg24;
			data_p.EED = arg24;
			double arg25;
			in >> arg25;
			data_p.GRTDP = arg25;
			double arg26;
			in >> arg26;
			data_p.TEC = arg26;
			double arg27;
			in >> arg27;
			data_p.WSSG = arg27;
			double arg28;
			in >> arg28;
			data_p.WSSL = arg28;
			double arg29;
			in >> arg29;
			data_p.WSSD = arg29;
			double arg30;
			in >> arg30;
			data_p.SLNG = arg30;
			double arg31;
			in >> arg31;
			data_p.SLNS = arg31;
			double arg32;
			in >> arg32;
			data_p.SNCG = arg32;
			double arg33;
			in >> arg33;
			data_p.SNCS = arg33;
			double arg34;
			in >> arg34;
			data_p.GNC = arg34;
			double arg35;
			in >> arg35;
			data_p.MXNUP = arg35;
			double arg36;
			in >> arg36;
			data_p.WSSN = arg36;	
			double arg37;
			in >> arg37;
			data_p.TBD = arg37;
			double arg38;
			in >> arg38;
			data_p.TP1D = arg38;
			double arg39;
			in >> arg39;
			data_p.TP2D = arg39;
			double arg40;
			in >> arg40;
			data_p.TCD = arg40;
			double arg41;
			in >> arg41;
			data_p.CPP = arg41;
			double arg42;
			in >> arg42;
			data_p.ppsen = arg42;
			double arg43;
			in >> arg43;
			data_p.ttSWEM = arg43;
			double arg44;
			in >> arg44;
			data_p.ttEMR1 = arg44;
			double arg45;
			in >> arg45;
			data_p.ttR1R3 = arg45;
			double arg46;
			in >> arg46;
			data_p.ttR3R5 = arg46;
			double arg47;
			in >> arg47;
			data_p.ttR5R7 = arg47;
			double arg48;
			in >> arg48;
			data_p.ttR7R8 = arg48;
			double arg49;
			in >> arg49;
			data_p.ttBRP = arg49;
			double arg50;
			in >> arg50;
			data_p.ttTRP = arg50;
			double arg51;
			in >> arg51;
			data_p.ttWSD = arg51;
			double arg52;
			in >> arg52;
			data_p.ttR1TLM = arg52;
			double arg53;
			in >> arg53;
			data_p.ttR1TLP = arg53;
			double arg54;
			in >> arg54;
			data_p.ttRUE = arg54;
			double arg55;
			in >> arg55;
			data_p.ttBSG = arg55;
			double arg56;
			in >> arg56;
			data_p.ttTSG = arg56;
			double arg57;
			in >> arg57;
			data_p.ttBRG = arg57;
			double arg58;
			in >> arg58;
			data_p.ttTRG = arg58;
			double arg59;
			in >> arg59;
			data_p.ttBNF = arg59;
			double arg60;
			in >> arg60;
			data_p.TRESH = arg60;
			double arg61;
			in >> arg61;
			data_p.ttDKill = arg61;
			double arg62;
			in >> arg62;
			data_p.LtFtsw = arg62;
			double arg63;
			in >> arg63;
			data_p.LtWdDur = arg63;
			double arg64;
			in >> arg64;
			data_p.vpd_resp = arg64;
			double arg65;
			in >> arg65;
			data_p.vpd_cr = arg65;
		}
	}
};

