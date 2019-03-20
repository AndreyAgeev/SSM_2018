#pragma once
#include "model.h"
#include "parametrs.h"

#include <highfive/H5Attribute.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

#include "functions-grammar.h"
#include <mlpack/methods/lars/lars.hpp>
#include <mlpack/core/data/load.hpp>
#include <cmath>
using namespace std;
using namespace HighFive;

class Nlreg: public QObject
{	
  private:
	QString h5_file_name;
	QString funcs_file_name;
	double * phenotype;
	int nFunctions;
	int wordlength;
	double epsilon;
	double alpha;
	int PRINT_TRACE;
	bool force;
	bool scale;
	bool intersept;
	bool positive;
	bool zeromean;
	bool extra_covar;
	int n_cycles;
	int n_weights;
	int n_summands;
	bool make_corrections;
	bool training_error_flag;
	

	int ROW;
	int *phenomask;
	int phenolength;
	GrammarContainer* grc;

	int nSamples;
	int ns;
	int nPredictors;
	int nn;

	std::vector<std::string> measurements;
	std::vector<std::string> species;
	std::vector<std::vector<double> > vec, resp;
	std::vector<std::vector<double> > gr_covar;
	std::vector<std::string> gr_names;
	int nGrCovar;

	std::vector<std::vector<double> > ecovar;
	std::vector<std::string> ec_names;
	int nEcovar;

	vector<double> conc;
	vector<int> genotype;
	vector<double> b_opt; ///////////////my_code
	arma::vec response;

	bool useCholesky;

	int iter_max;
	double err_tol;
	double par_alpha;

	arma::mat Fpoints;

	vector< vector<double> > eval_temp_photo;

	vector<int> read_genotype(vector<double>& conc) {
		funcs_file_name = "C:\\project\\SSM\\SSM_improved\\SSM_improved\\deep.txt";
		QFile file(funcs_file_name);
		vector<int> gt;
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
			cout << "Can't open " << funcs_file_name.toStdString() << endl;
			return gt;
		}
		QTextStream in(&file);
		if (force) {
			QString nnm;
			in >> nnm;
		//	cout << nnm.toStdString () << endl;
		}
		for (size_t i = 0; i < nFunctions * wordlength; ++i) {
			int arg;
			in >> arg;
		//	cout << arg << endl;
			gt.push_back(arg);
		//	cout << "genotype " << endl;
		}
		////////////////////////my_code
		for (size_t i = 0; i < nFunctions; ++i) {
			double arg;
			in >> arg;
			b_opt.push_back(arg);
		}





		vector<double> concs;
		if (force) {
			QString nnm;
			in >> nnm;
		}
		for (size_t i = 0; i < nFunctions * wordlength; ++i) {
			double arg;
			in >> arg;
			concs.push_back(arg);
		}
		conc = concs;
	//			for (size_t i = 0; i < nFunctions * wordlength; ++i) { cout << gt[i] << endl; }
		return gt;
	};
	void print_trace(arma::vec betaOpt, arma::uvec i_sub, arma::vec opt_weights, arma::rowvec col_min, arma::rowvec col_max) {
		GrammarNode *retFtn;
		double er;
		arma::vec predictions;
		for (size_t i = 0; i < nFunctions; ++i) {
			cout << i + 1 << " = ";
			retFtn = grc->get_nth_tree(i);
			retFtn->coprint();
			cout << " = " << col_max(i + 1) << " = " << col_min(i + 1);
			cout << endl;
		}
		for (size_t i = 0; i < nGrCovar; ++i) {
			for (size_t k = 0; k < nFunctions; ++k) {
				cout << i * nFunctions + k + nFunctions + 1 << " = ";
				cout << gr_names[i] << "*" << "f[" << k + 1 << "]";
				cout << " = " << col_max(i * nFunctions + k + nFunctions + 1) << " = " << col_min(i * nFunctions + k + nFunctions + 1);
				cout << endl;
			}
		}
		/*
				for (size_t i = 0; i < nEcovar; ++i) {
					cout << i + nFunctions + 1 << " ";// + nFunctions * nGrCovar
					cout << ec_names[i];
					cout << " " << col_max(i + nFunctions + 1) << " " << col_min(i + nFunctions + 1);// + nFunctions * nGrCovar
					cout << endl;
				}
		*/
		i_sub.t().print();
		size_t i, j;
		cout << "nel " << betaOpt.n_elem << endl;
		for (i = 0; i < betaOpt.n_elem; ++i) {
			cout << betaOpt[i] << " ";
		}
		cout << endl;
		cout << "F= ";
		int flag = 0;
		for (i = 0; i < betaOpt.n_elem; ++i) {
			if (fabs(betaOpt[i]) > 0) {
				if (flag == 1 && betaOpt[i] > 0) cout << "+";
				cout << betaOpt[i];
				j = i_sub(i);
				if (0 < j && j <= nFunctions) {
					cout << "*";
					retFtn = grc->get_nth_tree(j - 1);
					retFtn->coprint();
				}
				else if (j > nFunctions) {
					int jj = j - nFunctions - 1;
					int ii = jj / nFunctions;
					int kk = jj % nFunctions;
					cout << "*" << gr_names[ii] << "*" << "f[" << kk + 1 << "]";
				}
				/*
								} else if (j > nFunctions) { // + nFunctions * nGrCovar
									cout << "*" << ec_names[j - nFunctions - 1];
								}
				*/
				flag = 1;
			}
		}
		cout << endl;
		arma::mat func_points = Fpoints.cols(i_sub);
		predictions = func_points * betaOpt;
		er = norm(predictions - response);
		cout << "error=" << er << endl;
		arma::vec d = response - mean(response);
		arma::vec d1 = response - predictions;
		double r2 = dot(d1, d1) / dot(d, d);
		cout << "R2 = " << 1.0 - r2 << endl;
		cout << "R2adj =" << 1.0 - r2 * ((double)(nSamples - 1)) / ((double)(nSamples - i_sub.n_elem - 1)) << endl;
		for (size_t i = 0; i < nFunctions * (wordlength + 1) + 1 + nFunctions * nGrCovar; ++i) { // + nEcovar
			cout << phenotype[i] << " ";
		}
		cout << endl;
		opt_weights.t().print();
	}
	double run_training(arma::vec& betaOpt, arma::uvec& i_sub, arma::vec& opt_weights, arma::rowvec& col_min, arma::rowvec& col_max) {
		GrammarNode *retFtn;
		double training_error;
		int offst = 0;
		for (size_t j = 0; j < nSamples; j++) {
			Fpoints(j, 0) = 1;
		}
		eval_temp_photo.resize(nFunctions);
		std::ofstream in("eval.txt");
		for (size_t i = 0; i < nFunctions; i++) {
			eval_temp_photo[i].resize(nSamples);
			retFtn = grc->get_nth_tree(i);
			for (size_t j = 0; j < nSamples; j++) {
				Fpoints(j, i + 1) = retFtn->eval(vec[j]);///////////////////////нужно брать погоду,то есть семлов у нас намного больше или как?
			//	cout << "EVAL = " << retFtn->eval(vec[j]) << endl;
				eval_temp_photo[i][j] = retFtn->eval(vec[j]);
				in << "EVAL = " << eval_temp_photo[i][j] << endl;
			}
		}
		in.close();
	/*	for (size_t i = 0; i < nGrCovar; i++) {
			for (size_t k = 0; k < nFunctions; k++) {
				for (size_t j = 0; j < nSamples; j++) {
					Fpoints(j, nFunctions + 1 + i * nFunctions + k) = gr_covar[j][i] * Fpoints(j, k + 1);
				}
			}
		}*/
		/*
				for (size_t i = 0; i < nEcovar; ++i) {
					for (size_t j = 0; j < nSamples; ++j) {
						Fpoints(j, i + nFunctions + 1) = ecovar[j][i]; // + nFunctions * nGrCovar
					}
				}
		*/
		col_min = min(Fpoints, 0);
		if (intersept) {
			col_min(0) = 0.5;
			offst = 1;
		}
		col_max = max(Fpoints, 0);
		arma::rowvec col_nan = mean(Fpoints, 0);

		/*
				if (positive) {
					i_sub = arma::find(col_min >= 0 && col_max > col_min && col_max != arma::datum::inf && col_min != arma::datum::inf);
				} else {
					i_sub = arma::find(col_max > col_min && col_max != arma::datum::inf && col_min != arma::datum::inf);
				}
		*/
		if (positive) {
			//			i_sub = arma::unique(arma::join_vert(arma::find(col_min >= 0 && col_max > col_min), arma::find_finite(col_nan)));
			//			i_sub = arma::intersect(arma::find(col_min >= 0 && col_max > col_min), arma::find_finite(col_nan));
			i_sub = arma::find(col_min >= 0 && col_max > col_min);
			i_sub = i_sub.elem(arma::find_finite(col_nan.elem(i_sub)));
		}
		else {
			//			i_sub = arma::unique(arma::join_vert(arma::find(col_max > col_min), arma::find_finite(col_nan)));
			//			i_sub = arma::intersect(arma::find(col_max > col_min), arma::find_finite(col_nan));
			i_sub = arma::find(col_max > col_min);
			i_sub = i_sub.elem(arma::find_finite(col_nan.elem(i_sub)));
		}
		if (i_sub.n_elem < 1) {
			return 2.76e16;
			//			cout << "error=" << 2.76e16 << endl;
			//			exit(-200);
		}
		/*
		arma::mat func_points = Fpoints.cols(i_sub);

	//	LARSW larsw(useCholesky, alpha, func_points, response);
		larsw.setPrintTrace(PRINT_TRACE);
		//training_error = larsw.TrainMaxWeight (betaOpt, epsilon);
		if (n_weights > 0 && nSamples > n_weights) {
			training_error = larsw.TrainRngWeight(betaOpt, epsilon, n_weights);
		}
		else {
			training_error = larsw.TrainMaxWeight(betaOpt, epsilon);
		}*/
		for (size_t i = 0; i < betaOpt.n_elem; ++i) {
			phenotype[i_sub(i) + nFunctions * wordlength] = betaOpt[i];
		}
		//opt_weights = larsw.get_opt_weights();
		return training_error;
		//*/
	}
	/*
	 * Nelder-Mead
	 * modified from
	 * https://github.com/kthohr/optim/blob/master/src/unconstrained/nm.cpp
	 */
	 //
	double nm_obj_fun(arma::vec tmp, arma::vec& betaOpt, arma::uvec& i_sub, arma::vec& opt_weights, arma::rowvec& col_min, arma::rowvec& col_max) {
		double training_error;
		size_t j = 0;
		for (size_t i = 0; i < phenolength; i++) {
			if (phenomask[i] == 1) {
				phenotype[i] = tmp(j++);
			}
		}
		training_error = run_training(betaOpt, i_sub, opt_weights, col_min, col_max);
		return training_error;
	}
	double run_nm_training(arma::vec& betaOpt, arma::uvec& i_sub, arma::vec& opt_weights, arma::rowvec& col_min, arma::rowvec& col_max) {
		double training_error;
		arma::vec init_out_vals(phenolength);
		size_t j = 0;
		for (size_t i = 0; i < phenolength; i++) {
			if (phenomask[i] == 1) {
				init_out_vals(j) = phenotype[i];
				j++;
			}
		}
		size_t n_vals = j;
		init_out_vals.resize(n_vals);
		// expansion / contraction parameters
		double par_par = 0.1;
		double par_par0 = 0.01 * par_par * par_par;
		double par_beta = 0.75 - 1.0 / (2.0 * n_vals);
		double par_gamma = 1.0 + 2.0 / n_vals;
		double par_delta = 1.0 - 1.0 / n_vals;
		if (PRINT_TRACE > 1) {
			std::cout << "Nelder-Mead: beginning search..." << endl;
			std::cout << "n_vals = " << n_vals << " par_beta = " << par_beta << " par_gamma = " << par_gamma << " par_delta = " << par_delta << endl;
			arma::cout << "Initial guess:\n" << init_out_vals.t() << "\n";
		}
		//
		// setup
		arma::vec simplex_fn_vals(n_vals + 1);
		arma::mat simplex_points(n_vals + 1, n_vals);

		simplex_fn_vals(0) = run_training(betaOpt, i_sub, opt_weights, col_min, col_max);
		simplex_points.row(0) = init_out_vals.t();

		// for (size_t i=1; i < n_vals + 1; i++) {
		//     simplex_points.row(i) = init_out_vals.t() + 0.05*arma::trans(unit_vec(i-1,n_vals));
		//     simplex_fn_vals(i) = opt_objfn(simplex_points.row(i).t(),nullptr,opt_data);
		// }

		for (size_t i = 1; i < n_vals + 1; i++) {
			simplex_points.row(i) = init_out_vals.t();
			if (init_out_vals(i - 1) != 0.0) {
				simplex_points(i, i - 1) += par_par * init_out_vals(i - 1);
			}
			else {
				simplex_points(i, i - 1) += par_par0;
			}
			simplex_fn_vals(i) = nm_obj_fun(simplex_points.row(i).t(), betaOpt, i_sub, opt_weights, col_min, col_max);
		}
		double min_val = simplex_fn_vals.min();
		int ind_min_val = index_min(simplex_fn_vals);
		//
		// begin loop

		if (PRINT_TRACE > 1) {
			std::cout << "Initialization Phase:";
			arma::cout << " Objective function value at each vertex:" << simplex_fn_vals.t() << endl;
			if (PRINT_TRACE > 2) {
				printf("\n");
				arma::cout << "Simplex matrix:\n" << simplex_points << "\n";
			}
		}
		int iter = 0;
		double err = 2 * err_tol;
		while (err > err_tol && iter < iter_max) {
			iter++;
			bool next_iter = false;
			// step 1
			arma::uvec sort_vec = arma::sort_index(simplex_fn_vals); // sort from low (best) to high (worst) values
			simplex_fn_vals = simplex_fn_vals(sort_vec);
			simplex_points = simplex_points.rows(sort_vec);
			// step 2
			arma::vec centroid = arma::trans(arma::sum(simplex_points.rows(0, n_vals - 1), 0)) / static_cast<double>(n_vals);
			arma::vec x_r = centroid + par_alpha * (centroid - simplex_points.row(n_vals).t());
			double f_r = nm_obj_fun(x_r, betaOpt, i_sub, opt_weights, col_min, col_max);
			if (f_r >= simplex_fn_vals(0) && f_r < simplex_fn_vals(n_vals - 1)) {   // reflected point is neither best nor worst in the new simplex
				simplex_points.row(n_vals) = x_r.t();
				next_iter = true;
			}
			// step 3
			if (!next_iter && f_r < simplex_fn_vals(0)) {   // reflected point is better than the current best; try to go farther along this direction
				arma::vec x_e = centroid + par_gamma * (x_r - centroid);
				double f_e = nm_obj_fun(x_e, betaOpt, i_sub, opt_weights, col_min, col_max);
				if (f_e < f_r) {
					simplex_points.row(n_vals) = x_e.t();
				}
				else {
					simplex_points.row(n_vals) = x_r.t();
				}
				next_iter = true;
			}
			// steps 4, 5, 6
			if (!next_iter && f_r >= simplex_fn_vals(n_vals - 1)) {   // reflected point is still worse than x_n; contract
				// steps 4 and 5
				if (f_r < simplex_fn_vals(n_vals)) {   // outside contraction
					arma::vec x_oc = centroid + par_beta * (x_r - centroid);
					double f_oc = nm_obj_fun(x_oc, betaOpt, i_sub, opt_weights, col_min, col_max);
					if (f_oc <= f_r) {
						simplex_points.row(n_vals) = x_oc.t();
						next_iter = true;
					}
				}
				else {   // inside contraction: f_r >= simplex_fn_vals(n_vals)
				 // x_ic = centroid - par_beta*(x_r - centroid);
					arma::vec x_ic = centroid + par_beta * (simplex_points.row(n_vals).t() - centroid);
					double f_ic = nm_obj_fun(x_ic, betaOpt, i_sub, opt_weights, col_min, col_max);
					if (f_ic < simplex_fn_vals(n_vals)) {
						simplex_points.row(n_vals) = x_ic.t();
						next_iter = true;
					}
				}
			}
			// step 6
			if (!next_iter) {   // neither outside nor inside contraction was acceptable; shrink the simplex toward x(0)
				for (size_t i = 1; i < n_vals + 1; i++) {
					simplex_points.row(i) = simplex_points.row(0) + par_delta * (simplex_points.row(i) - simplex_points.row(0));
				}
			}
			// check change in fn_val
			for (size_t i = 0; i < n_vals + 1; i++) {
				simplex_fn_vals(i) = nm_obj_fun(simplex_points.row(i).t(), betaOpt, i_sub, opt_weights, col_min, col_max);
			}
			//
			err = std::abs(min_val - simplex_fn_vals.max());
			min_val = simplex_fn_vals.min();
			ind_min_val = index_min(simplex_fn_vals);
			// printing
			if (PRINT_TRACE > 1) {
				std::cout << "Iteration = " << iter;
				std::cout << " min_val = " << min_val;
				std::cout << " err = " << err;
				std::cout << " ind = " << ind_min_val;
				arma::cout << " Current optimal input values: ";
				arma::cout << simplex_points.row(ind_min_val);
				if (PRINT_TRACE > 2) {
					printf("\n");
					arma::cout << "    Objective function value at each vertex:\n" << simplex_fn_vals.t() << "\n";
					arma::cout << "    Simplex matrix:\n" << simplex_points << "\n";
				}
			}
		}
		if (PRINT_TRACE > 1) {
			std::cout << "Nelder-Mead: search completed.\n";
		}
		//
		training_error = nm_obj_fun(simplex_points.row(index_min(simplex_fn_vals)).t(), betaOpt, i_sub, opt_weights, col_min, col_max);
		return training_error;
	}
public:
	Nlreg(QString file_name, QString func_file_name, int nF, int wl, double alfa = 0, double epsi = 0.003, int print_trace = 0, QObject *parent = 0) : QObject(parent), h5_file_name(file_name), funcs_file_name(func_file_name) {
				std::cout << "Got " << h5_file_name.toStdString() << " file" << std::endl;
				std::cout << "Got " << funcs_file_name.toStdString() << " funcs" << std::endl;
		nFunctions = nF;
		wordlength = wl;
		epsilon = epsi;
		alpha = alfa;
		PRINT_TRACE = print_trace;
		force = false;
		scale = false;
		intersept = false;
		positive = true;////////было false
		zeromean = false;
		extra_covar = false;
		make_corrections = false;
		n_cycles = 3;
		n_summands = 0;
		n_weights = 0;
		training_error_flag = false;
		useCholesky = true;
		iter_max = 20;
		err_tol = 1e-3;
		par_alpha = 0.5;
	}
	void setForce(bool arg) { force = arg; }
	void setScale(bool arg) { scale = arg; }
	void setIntersept(bool arg) { intersept = arg; }
	void setPositive(bool arg) { positive = arg; }
	void setZeroMean(bool arg) { zeromean = arg; }
	void setExtraCovar(bool arg) { extra_covar = arg; }
	void setMakeCorrections(bool arg) { make_corrections = arg; }
	void setNCycles(int arg) { n_cycles = arg; }
	void setNWeights(int arg) { n_weights = arg; }
	void setNSummands(int arg) { n_summands = arg; }
	void setTrainingError(bool arg) { training_error_flag = arg; }
public:
	void create_data(const Data & data)
	{
		int Height = data.data_h5.tmax.size();
		int Width = 4;
		vec.resize(Height);
		for (int i = 0; i < Height; i++)
			vec[i].resize(Width);
		for (int i = 0; i < Height; i++)
		{
			ROW += 1;
			for (int j = 0; j < Width; j++)
			{
				if (j == 0)
					vec[i][j] = (data.data_h5.tmax[ROW]);//TMAX
				if (j == 1)
					vec[i][j] = (data.data_h5.tmin[ROW]);//TMIN
				if (j == 2)
					vec[i][j] = (data.data_h5.srad[ROW]);//SRAD
				if (j == 3)
					vec[i][j] = (data.data_h5.rain[ROW]);//RAIN
			}
		}

	}
	void run_nlreg(const Data & data)
	{
		// Load data
		nSamples = 100;
		ns = -1;
		nn = -1;
		nGrCovar = 0;
		nEcovar = 0;
		genotype = read_genotype(conc);
	/*	if (zeromean) {
			for (size_t i = 0; i < nFunctions * wordlength; i++) {
				conc[i] -= (double)genotype[i];
				//				cout << i << " " << conc[i] << " " << genotype[i] << endl;
			}
		}
	/*	try {
			File file(h5_file_name.toStdString(), File::ReadOnly);
			DataSet a0_read = file.getDataSet("data");
			DataSpace space = a0_read.getSpace();
			nSamples = space.getDimensions()[0]; // number of samples
			nPredictors = space.getDimensions()[1]; // number of measurements
			a0_read.read(vec);
			DataSet b_read = file.getDataSet("response");
			space = b_read.getSpace();
			ns = space.getDimensions()[0]; // number of samples
			nn = space.getDimensions()[1]; // number of measurements
			b_read.read(resp);
			assert(ns == nSamples);
			assert(nn == 1);
			DataSet a1_read = file.getDataSet("species");
			a1_read.read(species);
			DataSet a2_read = file.getDataSet("measurements");
			a2_read.read(measurements);
			if (extra_covar) {
				DataSet c0_read = file.getDataSet("gr_covar");
				DataSpace space = c0_read.getSpace();
				nGrCovar = space.getDimensions()[1]; // number of measurements
				c0_read.read(gr_covar);
				DataSet c1_read = file.getDataSet("gr_names");
				c1_read.read(gr_names);
			}
			if (false /*extra_covar) {
				DataSet c0_read = file.getDataSet("ecovar");
				DataSpace space = c0_read.getSpace();
				nEcovar = space.getDimensions()[1]; // number of measurements
				c0_read.read(ecovar);
				DataSet c1_read = file.getDataSet("ec_names");
				c1_read.read(ec_names);
			}
		}
		catch (Exception& err) {
			std::cerr << err.what() << std::endl;
		}*/
		/////////////////////////////////my_code
		measurements = { "TMAX", "TMIN", "RAIN", "SRAD" };
		nn = nPredictors = measurements.size();
		nSamples = 1;///?
		nGrCovar = 0;///?
	//	GiveSize(data);// прописать measurments - тк это предикторы,а species - места посадки
		//response = std2arvec(resp, ns, 0);////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Load and decode functions
		ofstream f("outdataread.txt");
		int begin_count_samples = ROW;
		create_data(data);
		/*for (int i = 0; i < vec.size(); i++)
		{
			
			f << endl;
			for (int j = 0; j < vec[i].size(); j++)
			{
				cout << "VEC[I] SIZE == " << vec[i].size(); 
				f << vec[i][j] << " ";
			}
		}
		f.close();*/
		nSamples = ROW - begin_count_samples;
		phenotype = new double[nFunctions * (wordlength + 1) + 1 + nFunctions * nGrCovar]; /* +nEcovar */
		phenomask = new int[nFunctions * (wordlength + 1) + 1 + nFunctions * nGrCovar]; /* +nEcovar */
		phenolength = nFunctions * (wordlength + 1) + 1 + nFunctions * nGrCovar;
		grc = new GrammarContainer(measurements, nFunctions, n_cycles, make_corrections);
		int first = 0, last = wordlength;
		for (size_t i = 0; i < nFunctions; ++i) {
			std::vector<int> gt = std::vector<int>(genotype.begin() + first, genotype.begin() + last);
			std::vector<double> cn = std::vector<double>(conc.begin() + first, conc.begin() + last);
			first += wordlength;
			last += wordlength;
			grc->build_nth_tree(gt, cn, i, &phenotype[i * wordlength], &phenomask[i * wordlength]);
		}
		// Train coefficients
	    arma::vec betaOpt;
		arma::uvec i_sub;
		Fpoints = arma::mat(nSamples, nFunctions + 1 + nFunctions * nGrCovar); /* +nEcovar ;*/
		arma::vec opt_weights;
		arma::rowvec col_min;
		arma::rowvec col_max;
		double training_error;
		if (positive) {
			training_error = run_training(betaOpt, i_sub, opt_weights, col_min, col_max);
		}
	/*	else {
			training_error = run_nm_training(betaOpt, i_sub, opt_weights, col_min, col_max);
		}
		// Print results
		cout << "terror=" << training_error << endl;
		if (PRINT_TRACE > 0) {
			print_trace(betaOpt, i_sub, opt_weights, col_min, col_max);
		}*/
	}
   	vector<vector<double > > createFunction(const Data& data, int _ROW)
    {
		ROW = _ROW;
	    run_nlreg(data);
		//cout << "begin print nlreg inside" << endl;
		for (int i = 0; i < eval_temp_photo.size(); i++)
		{
		//	cout << endl;
			for (int j = 0; j < eval_temp_photo[i].size(); j++)
			{
				if (std::fpclassify(eval_temp_photo[i][j]) == FP_INFINITE || std::fpclassify(eval_temp_photo[i][j]) == FP_NAN)
					eval_temp_photo[i][j] = 0.0;
			   // cout << eval_temp_photo[i][j] << endl;

			}
		}
	//	cout << "end print nlreg inside" << endl;
	//	cout << "end inside nlreg" << endl;
		return eval_temp_photo;
	}
};



