/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * larsw.h
 * Copyright (C) 2018 Konstantin Kozlov <mackoel@gmail.com>
 *
 * nlreg is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * nlreg is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _LARSW_H_
#define _LARSW_H_

#include <mlpack/methods/lars/lars.hpp>
#include <mlpack/core/data/load.hpp>

arma::mat std2armat(std::vector<std::vector<double> > &vec, int n_rows, int n_cols, int offset);

arma::mat std2arvec(std::vector<std::vector<double> > &vec, int n_rows, int offset);

arma::mat std2armat_full(std::vector<std::vector<double> > &vec, int n_rows, int n_cols);

class LARSW
{
public:
	LARSW(bool useCholesky_, double alpha_, std::vector<std::vector<double> >& vec, int nSamples, int nColumns, int nFunctions, int offset);
	LARSW(bool useCholesky_, double alpha_, arma::mat& points, arma::vec& response, int nSamples, int nFunctions);
	LARSW(bool useCholesky_, double alpha_, arma::mat& points, arma::vec& response);
	void Train(double lambda, arma::vec& beta);
	void TrainWeighted(double lambda, arma::vec& beta, arma::vec& weights);
	void TrainWeighted(arma::vec& beta, arma::vec& weights, int nfolds = 4, int nlambdas = 100, double lambda_max_init = -1, double lambda_max_ratio = 2e3);
	double TrainMaxWeight(arma::vec& beta, double epsilon, int nfolds = 4, int nlambdas = 100, double lambda_max_init = -1, double lambda_max_ratio = 2e3);
	double TrainRngWeight(arma::vec& beta, double epsilon, int nweights = 5, int nfolds = 4, int nlambdas = 100, double lambda_max_init = -1, double lambda_max_ratio = 2e3);
	double getError();
	void setPrintTrace(int arg) { PRINT_TRACE = arg; }
	arma::mat getX() { return X; }
	arma::vec get_opt_weights() { return opt_weights; };
private:
	double alpha;
	bool useCholesky;
	arma::vec y;
	arma::mat X;
	bool transposeData;
	bool rowMajor;
	int PRINT_TRACE;
	arma::vec opt_weights;
};

#endif // _LARSW_H_

