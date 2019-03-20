#pragma once
#include "model.h"
#include "parametrs.h"

#include <highfive/H5Attribute.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

#include "functions-grammar.h"
#include "larsw.h"

using namespace std;
using namespace HighFive;

class Nlreg
{
	Parametrs _param;
	int PRINT_TRACE = 1;
	int training_error_flag;
	double *phenotype;
    public:
   	    void createFunction(const QString& file_name, Parametrs &param)
	    {
			_param = param;
			read_ini_nlreg_param(file_name);
			run_nlreg();
			param = _param;
	    }
    private:
		void read_ini_nlreg_param(const QString& file_name)
		{
			QSettings sett(file_name, QSettings::IniFormat);
			sett.beginGroup("Params");

			_param.nl_param.nFunctions = sett.value("n", 6).toInt();
			_param.nl_param.epsilon = sett.value("e", 10).toDouble();
			_param.nl_param.l = sett.value("l", 0.003).toDouble();
			_param.nl_param.alpha = sett.value("a", 0.05).toDouble();
			_param.nl_param.t = sett.value("t", 1).toDouble();
			_param.nl_param.force = sett.value("f", 1).toBool();
			_param.nl_param.intersept = sett.value("i", 1).toBool();
			_param.nl_param.c = sett.value("c", 1).toDouble();
			_param.nl_param.m = sett.value("m", 1).toDouble();
			_param.nl_param.d = sett.value("d", 1).toDouble();
			_param.nl_param.r = sett.value("r", 1).toDouble();
			_param.nl_param.q = sett.value("q", 1).toDouble();
			_param.nl_param.g = sett.value("g", 1).toDouble();
			_param.nl_param.b = sett.value("b", 20).toDouble();
			_param.nl_param.scale = sett.value("s", 0).toBool();
			_param.nl_param.wordLength = sett.value("w", 20).toInt();
			_param.nl_param.h5_file_name = sett.value("file_name", QString("chickpea-csminds.h5")).toString();
			_param.nl_param.funcs_file_name = sett.value("func_file_name", QString("chickpea-csminds.h5")).toString();//////////
			_param.nl_param.n_cycles = 3;
			_param.nl_param.positive = false;
			_param.nl_param.make_corrections = false;
			sett.endGroup();
		}
		vector<int> read_genotype(vector<double>& conc) {
			QFile file(_param.nl_param.funcs_file_name);
			vector<int> gt;
			if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
				cout << "Can't open " << _param.nl_param.funcs_file_name.toStdString() << endl;
				return gt;
			}
			QTextStream in(&file);
			if (_param.nl_param.force) {
				QString nnm;
				in >> nnm;
				/*			cout << nnm.toStdString () << endl;*/
			}
			for (size_t i = 0; i < _param.nl_param.nFunctions * _param.nl_param.wordLength; ++i) {
				int arg;
				in >> arg;
				gt.push_back(arg);
				//			cout << gt[i] << endl;
			}
			vector<double> concs;
			if (_param.nl_param.force) {
				QString nnm;
				in >> nnm;
				/*			cout << nnm.toStdString () << endl;*/
			}
			for (size_t i = 0; i < _param.nl_param.nFunctions * _param.nl_param.wordLength; ++i) {
				double arg;
				in >> arg;
				concs.push_back(arg);
				//			cout << gt[i] << endl;
			}
			conc = concs;
			//		for (size_t i = 0; i < nFunctions * wordlength; ++i) { cout << gt[i] << endl; }
			return gt;
		};
		void run_nlreg()
		{
			std::vector<std::string> measurements;
			std::vector<std::string> species;
			std::vector<std::vector<double> > vec, resp;

			int nSamples = 100, ns = -1, nPredictors;
			int nn = -1;

			//		vector<int> genotype = {0, 1, 0, 0, 1, 1, 0, 0};
			//		vector<int> genotype = {0, 4, 0, 4, 1};
			vector<double> conc;
			vector<int> genotype = read_genotype(conc);

			try {
				// Open a new file using the default property lists.
				File file(_param.nl_param.h5_file_name.toStdString(), File::ReadOnly);

				DataSet a0_read = file.getDataSet("data");
				DataSpace space = a0_read.getSpace();

				//			std::cout << "as " << space.getDimensions().size() << std::endl;
				//			std::cout << "ax " << space.getDimensions()[0] << std::endl;
				//			std::cout << "ay " << space.getDimensions()[1] << std::endl;

				nSamples = space.getDimensions()[0]; // number of samples
				nPredictors = space.getDimensions()[1]; // number of measurements

				a0_read.read(vec);

				DataSet b_read = file.getDataSet("response");
				space = b_read.getSpace();

				//			std::cout << "bs " << space.getDimensions().size() << std::endl;
				//			std::cout << "bx " << space.getDimensions()[0] << std::endl;
				//			std::cout << "by " << space.getDimensions()[1] << std::endl;

				ns = space.getDimensions()[0]; // number of samples
				nn = space.getDimensions()[1]; // number of measurements

				b_read.read(resp);

				assert(ns == nSamples);
				assert(nn == 1);
				/*
							for (size_t i = 0; i < 10; ++i) {
								for (size_t j = 0; j < 10; ++j) {
									std::cout << vec[i][j] << " " ;
								}
								std::cout << std::endl;
							}
				*/
				DataSet a1_read = file.getDataSet("species");
				a1_read.read(species);

				//		for (size_t i = 0; i < species.size(); ++i) {
				/*			for (size_t i = 0; i < 10; ++i) {
								std::cout << species[i] << " " << std::endl;
							}
				*/
				DataSet a2_read = file.getDataSet("measurements");
				a2_read.read(measurements);

				//		for (size_t i = 0; i < measurements.size(); ++i) {
				/*			for (size_t i = 0; i < 10; ++i) {
								std::cout << measurements[i] << " " << std::endl;
							}
				*/
			}
			catch (Exception& err) {
				// catch and print any HDF5 error
				std::cerr << err.what() << std::endl;
			}

			std::vector<double> col_min, col_max, col_max_min;
			double resp_min, resp_max, resp_max_min;
			if (_param.nl_param.scale) {
				for (size_t i = 0; i < nPredictors; ++i) {
					double c_m = DBL_MAX;
					double c_x = -DBL_MAX;
					for (size_t j = 0; j < nSamples; ++j) {
						c_m = (c_m < vec[j][i]) ? c_m : vec[j][i];
						c_x = (c_x > vec[j][i]) ? c_x : vec[j][i];
					}
					col_min.push_back(c_m);
					col_max.push_back(c_x);
				}
				resp_min = DBL_MAX;
				resp_max = -DBL_MAX;
				for (size_t j = 0; j < nSamples; ++j) {
					resp_min = (resp_min < resp[j][0]) ? resp_min : resp[j][0];
					resp_max = (resp_max > resp[j][0]) ? resp_max : resp[j][0];
				}
				resp_max_min = resp_max - resp_min;
				resp_max_min = (resp_max_min == 0) ? 1 : resp_max_min;
				for (size_t i = 0; i < nPredictors; ++i) {
					double res = col_max[i] - col_min[i];
					res = (res == 0) ? 1 : res;
					col_max_min.push_back(res);
				}

				for (size_t i = 0; i < nSamples; ++i) {
					resp[i][0] = (resp[i][0] - resp_min) / resp_max_min;
					for (size_t j = 0; j < nPredictors; ++j) {
						vec[i][j] = (vec[i][j] - col_min[j]) / col_max_min[j];
					}
				}
			}

			GrammarNode *retFtn;
			/*
					retFtn = new ConstNode();
					std::cout << *retFtn << endl;

					retFtn = new InputNode(3, "name");
					std::cout << *retFtn << endl;

					retFtn = new Add();
					std::cout << *retFtn << endl;

					retFtn = new Subtract();
					std::cout << *retFtn << endl;

					retFtn = new Multiply();
					std::cout << *retFtn << endl;

					retFtn = new Divide();
					std::cout << *retFtn << endl;
			*/
			bool useCholesky = true;
			double lambda = 0.8;
			arma::vec betaOpt;
			arma::vec predictions;
			double er;
			double training_error;
			if (_param.nl_param.intersept) {
				phenotype = new double[_param.nl_param.nFunctions * (_param.nl_param.wordLength + 1) + 1];
				GrammarContainer* grc = new GrammarContainer(measurements, _param.nl_param.nFunctions, _param.nl_param.n_cycles, _param.nl_param.make_corrections);
				arma::mat Fpoints(nSamples, _param.nl_param.nFunctions + 1);
				int first = 0, last = _param.nl_param.wordLength;
				for (size_t j = 0; j < nSamples; ++j) {
					Fpoints(j, 0) = 1;
				}
				for (size_t i = 0; i < _param.nl_param.nFunctions; ++i) {
					std::vector<int> gt = std::vector<int>(genotype.begin() + first, genotype.begin() + last);
					std::vector<double> cn = std::vector<double>(conc.begin() + first, conc.begin() + last);
					first += _param.nl_param.wordLength;
					last += _param.nl_param.wordLength;
					grc->build_nth_tree(gt, cn, i, &phenotype[i * _param.nl_param.wordLength]);
					retFtn = grc->get_nth_tree(i);
					if (_param.nl_param.scale) retFtn->setScale(col_max_min, col_min);
					for (size_t j = 0; j < nSamples; ++j) {
						Fpoints(j, i + 1) = retFtn->eval(vec[j]);
					}
				}
				arma::rowvec col_min = min(Fpoints, 0);
				col_min(0) = 0.5;
				arma::rowvec col_max = max(Fpoints, 0);
				arma::uvec i_sub;
				if (_param.nl_param.positive) {
					i_sub = arma::find(col_min > 0 && col_max > col_min && col_max != arma::datum::inf && col_min != arma::datum::inf);
				}
				else {
					i_sub = arma::find(col_max > col_min && col_max != arma::datum::inf && col_min != arma::datum::inf);
				}
				if (PRINT_TRACE > 0) {
					for (size_t i = 0; i < _param.nl_param.nFunctions; ++i) {
						retFtn = grc->get_nth_tree(i);
						cout << i + 1 << " ";
						retFtn->coprint();
						cout << " " << col_max(i + 1) << " " << col_min(i + 1);
						cout << endl;
					}
					i_sub.t().print();
				}
				if (i_sub.n_elem < 1) {
					cout << "error=" << 2.76e16 << endl;
					exit(-200);
				}
				Fpoints = Fpoints.cols(i_sub);
				arma::vec response = std2arvec(resp, ns, 0);
				/*
							if (!Fpoints.is_finite()) {
					Fpoints.each_col( [](arma::vec& a){ cout << "o =" << a.is_finite() << endl; } );
					cout << "error=" << 2.76e16 << endl;
					exit(-200);
				}
				*/
				LARSW larsw(useCholesky, _param.nl_param.alpha, Fpoints, response);
				larsw.setPrintTrace(PRINT_TRACE);
				training_error = larsw.TrainMaxWeight(betaOpt, _param.nl_param.epsilon);
				if (PRINT_TRACE > 0) {
					size_t i;
					cout << "nel " << betaOpt.n_elem << endl;
					for (i = 0; i < betaOpt.n_elem; ++i) {
						cout << betaOpt[i] << " ";
						phenotype[i_sub(i) + _param.nl_param.nFunctions * _param.nl_param.wordLength] = betaOpt[i];
					}
					cout << endl;
					cout << "F= ";
					int flag = 0;
					if (fabs(betaOpt[0]) != 0) {
						cout << betaOpt[0];
						flag = 1;
					}
					for (i = 1; i < betaOpt.n_elem - 1; ++i) {
						retFtn = grc->get_nth_tree(i_sub(i) - 1);
						if (fabs(betaOpt[i]) > 0) {
							if (flag == 1 && betaOpt[i] > 0) cout << "+";
							cout << betaOpt[i] << "*";
							retFtn->coprint();
							flag = 1;
						}
					}
					i = betaOpt.n_elem - 1;
					retFtn = grc->get_nth_tree(i_sub(i) - 1);
					if (fabs(betaOpt[i]) > 0) {
						if (flag == 1 && betaOpt[i] > 0) cout << "+";
						cout << betaOpt[i] << "*";
						retFtn->coprint();
					}
					cout << endl;
				}
				cout << "terror=" << training_error << endl;
				if (!training_error_flag || PRINT_TRACE > 0) {
					predictions = Fpoints * betaOpt;
					er = norm(predictions - response);
					cout << "error=" << er << endl;
				}
				if (PRINT_TRACE > 0) {
					arma::vec d = response - mean(response);
					arma::vec d1 = response - predictions;
					double r2 = dot(d1, d1) / dot(d, d);
					cout << "R2 = " << 1.0 - r2 << endl;
					cout << "R2adj =" << 1.0 - r2 * ((double)(nSamples - 1)) / ((double)(nSamples - _param.nl_param.nFunctions - 1)) << endl;
					for (size_t i = 0; i < _param.nl_param.nFunctions * (_param.nl_param.wordLength + 1) + 1; ++i) {
						cout << phenotype[i] << " ";
					}
					cout << endl;
				}
			}
			else {
				phenotype = new double[_param.nl_param.nFunctions * (_param.nl_param.wordLength + 1)];
				GrammarContainer* grc = new GrammarContainer(measurements, _param.nl_param.nFunctions, _param.nl_param.n_cycles, _param.nl_param.make_corrections);
				arma::mat Fpoints(nSamples, _param.nl_param.nFunctions);
				int first = 0, last = _param.nl_param.wordLength;
				for (int i = 0; i < _param.nl_param.nFunctions; ++i) {
					std::vector<int> gt = std::vector<int>(genotype.begin() + first, genotype.begin() + last);
					std::vector<double> cn = std::vector<double>(conc.begin() + first, conc.begin() + last);
					first += _param.nl_param.wordLength;
					last += _param.nl_param.wordLength;
					grc->build_nth_tree(gt, cn, i, &phenotype[i * _param.nl_param.wordLength]);
					retFtn = grc->get_nth_tree(i);
					if (_param.nl_param.scale) retFtn->setScale(col_max_min, col_min);
					for (size_t j = 0; j < nSamples; ++j) {
						Fpoints(j, i) = retFtn->eval(vec[j]);
					}
				}
				arma::rowvec col_min = min(Fpoints, 0);
				arma::rowvec col_max = max(Fpoints, 0);
				arma::uvec i_sub;
				if (_param.nl_param.positive) {
					i_sub = arma::find(col_min > 0 && col_max > col_min && col_max != arma::datum::inf && col_min != arma::datum::inf);
				}
				else {
					i_sub = arma::find(col_max > col_min && col_max != arma::datum::inf && col_min != arma::datum::inf);
				}
				if (PRINT_TRACE > 0) {
					for (size_t i = 0; i < _param.nl_param.nFunctions; ++i) {
						retFtn = grc->get_nth_tree(i);
						cout << i << " ";
						retFtn->coprint();
						cout << " " << col_max(i) << " " << col_min(i);
						cout << endl;
					}
					i_sub.t().print();
				}
				if (i_sub.n_elem < 1) {
					cout << "error=" << 2.76e16 << endl;
					exit(-200);
				}
				Fpoints = Fpoints.cols(i_sub);
				arma::vec response = std2arvec(resp, ns, 0);

				LARSW larsw(useCholesky, _param.nl_param.alpha, Fpoints, response);
				larsw.setPrintTrace(PRINT_TRACE);
				/*
						 * Currently we are interested only in the last form of training.
						 *
						larsw.Train(lambda, betaOpt);
						if (PRINT_TRACE > 0) {
							cout << "nel " << betaOpt.n_elem << endl;
							for (size_t i = 0; i < betaOpt.n_elem; ++i)
								cout << i << " " << betaOpt[i] << endl;
						}

						arma::vec predictions = Fpoints * betaOpt;
						double er = norm(predictions - response);
						cout << "error " << er << endl;
						arma::vec weights(nSamples, 1);
						for (size_t i = 0; i < nSamples; ++i) {
							weights (i) = 1/((double)nSamples);
						}
						larsw.TrainWeighted (betaOpt, weights);
						if (PRINT_TRACE > 0) {
							cout << "nel " << betaOpt.n_elem << endl;
							for (size_t i = 0; i < betaOpt.n_elem; ++i)
								cout << i << " " << betaOpt[i] << endl;
						}

						predictions = Fpoints * betaOpt;
						er = norm(predictions - response);
						cout << "error " << er << endl;
				*/
				training_error = larsw.TrainMaxWeight(betaOpt, _param.nl_param.epsilon);
				if (PRINT_TRACE > 0) {
					size_t i;
					cout << "nel " << betaOpt.n_elem << endl;
					for (i = 0; i < betaOpt.n_elem; ++i) {
						cout << betaOpt[i] << " ";
						phenotype[i_sub(i) + _param.nl_param.nFunctions * _param.nl_param.wordLength] = betaOpt[i];
					}
					cout << endl;
					cout << "F= ";
					int flag = 0;
					for (i = 0; i < betaOpt.n_elem - 1; ++i) {
						retFtn = grc->get_nth_tree(i_sub(i));
						if (fabs(betaOpt[i]) > 0) {
							if (flag == 1 && betaOpt[i] > 0) cout << "+";
							cout << betaOpt[i] << "*";
							retFtn->coprint();
							flag = 1;
						}
					}
					i = betaOpt.n_elem - 1;
					retFtn = grc->get_nth_tree(i_sub(i));
					if (fabs(betaOpt[i]) > 0) {
						if (flag == 1 && betaOpt[i] > 0) cout << "+";
						cout << betaOpt[i] << "*";
						retFtn->coprint();
					}
					cout << endl;
				}
				cout << "terror=" << training_error << endl;
				if (!training_error_flag || PRINT_TRACE > 0) {
					predictions = Fpoints * betaOpt;
					er = norm(predictions - response);
					cout << "error=" << er << endl;
				}
				if (PRINT_TRACE > 0) {
					arma::vec d = response - mean(response);
					double r2 = er * er / dot(d, d);
					cout << "R2 = " << 1.0 - r2 << endl;
					cout << "R2adj =" << 1.0 - r2 * ((double)(nSamples - 1)) / ((double)(nSamples - _param.nl_param.nFunctions - 1)) << endl;
					for (size_t i = 0; i < _param.nl_param.nFunctions * (_param.nl_param.wordLength + 1); ++i) {
						cout << phenotype[i] << " ";
					}
					cout << endl;
				}
			}
		}
};


