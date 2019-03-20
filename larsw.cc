/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * larsw.cc
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
#include <iostream>
#include "larsw.h"

#include <mlpack/methods/lars/lars.hpp>
#include <mlpack/core/data/load.hpp>

using namespace mlpack;
using namespace mlpack::regression;

LARSW::LARSW(bool useCholesky_, double alpha_, std::vector<std::vector<double> >& vec, int nSamples, int nColumns, int nFunctions, int offset) {
	useCholesky = useCholesky_;
	alpha = alpha_;
	y = std2arvec(vec, nSamples, offset);
	X = std2armat(vec, nSamples, nColumns, offset);
	transposeData = false;
	rowMajor = !transposeData;
	PRINT_TRACE = 0;
}

LARSW::LARSW(bool useCholesky_, double alpha_, arma::mat& points, arma::vec& response, int nSamples, int nFunctions) {
	useCholesky = useCholesky_;
	alpha = alpha_;
	y = response;
	X = points;
	transposeData = false;
	rowMajor = !transposeData;
	PRINT_TRACE = 0;
}

LARSW::LARSW(bool useCholesky_, double alpha_, arma::mat& points, arma::vec& response) {
	useCholesky = useCholesky_;
	alpha = alpha_;
	y = response;
	X = points;
	transposeData = false;
	rowMajor = !transposeData;
	PRINT_TRACE = 0;
}

void LARSW::Train(double lambda, arma::vec& beta) {
	double lambda1 = lambda * (1 - alpha), lambda2 = 0.5 * lambda * alpha;
	arma::vec betaOpt;
	LARS lars(useCholesky, lambda1, lambda2);
	lars.Train(X, y, betaOpt, transposeData);
/*
		std::cout << "nel " << betaOpt.n_elem << std::endl;
		for (size_t i = 0; i < betaOpt.n_elem; ++i)
			std::cout << i << " " << betaOpt[i] << std::endl;

	arma::vec predictions, p = X * betaOpt;
	lars.Predict(X, predictions, rowMajor);
	for (size_t i = 0; i < y.n_elem; ++i)
		std::cout << i << " " << y[i] << " " << p[i] << " " << predictions[i] << std::endl;
*/
	beta = betaOpt;
}

void LARSW::TrainWeighted(double lambda, arma::vec& beta, arma::vec& weights) {
	double lambda1 = lambda * (1 - alpha), lambda2 = 0.5 * lambda * alpha;
	arma::vec betaOpt;
	LARS lars(useCholesky, lambda1, lambda2);
	arma::vec wy = sqrt(weights) % y;
	arma::mat wX = diagmat(sqrt(weights)) * X;
	lars.Train(wX, wy, betaOpt, transposeData);
	beta = betaOpt;
}

/*
# Unless a lambda sequence is provided by the user, generate it
n = nrow(X)
 lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
*/

void LARSW::TrainWeighted(arma::vec& beta, arma::vec& weights, int nfolds, int nlambdas, double lambda_max_init, double lambda_max_ratio) {
	arma::ivec folds_ind = arma::randi(y.n_elem, arma::distr_param(0, nfolds - 1));
	arma::uvec i_sub;
	arma::vec wyi;
	arma::mat wXi;
	arma::vec predictions;
	arma::vec betaOpt;
	arma::vec wy = sqrt(weights) % y;
	arma::mat wX = diagmat(sqrt(weights)) * X;
	int n = y.n_elem;
	double n1 = 1/((double)n);
	double lambda_max = (lambda_max_init > 0) ? lambda_max_init : n1 * max(abs(trans(X)*y));
	std::cout << "n1 = " << n1 << " lambda max = " << lambda_max << std::endl;
	double lambda_min = lambda_max / lambda_max_ratio;
	std::cout << "ratio = " << lambda_max_ratio << " lambda min = " << lambda_min << std::endl;
	double lambda_inc = pow(lambda_min / lambda_max, 1.0 / ((double)nlambdas));
	double k = 0, k_inc = 1 / ((double)nlambdas);
	double lambda = lambda_max;
	double lambda1, lambda2;
	double er;
	for (size_t i = 0; i < nlambdas; ++i) {
		std::cout << "i = " << i << " lambda = " << lambda;
		lambda1 = lambda * (1 - alpha);
		lambda2 = 0.5 * lambda * alpha;
		LARS lars(useCholesky, lambda1, lambda2);
		for (size_t j = 0; j < nfolds; ++j) {
			i_sub = arma::find(folds_ind != j);
			wyi = wy(i_sub);
			wXi = wX.rows(i_sub);
			lars.Train(wXi, wyi, betaOpt, transposeData);
			i_sub = arma::find(folds_ind == j);
			wyi = y(i_sub);
			wXi = X.rows(i_sub);
			predictions = wXi * betaOpt;
			er = norm(predictions - wyi);
			std::cout << " er[" << j << "] = " << er;
		}
		std::cout << std::endl;
		k += k_inc;
		lambda *= lambda_inc;
	}
	beta = betaOpt;
}

/*
    n <- length(y)
    epln <- epsilon / n
    n.inv <- 1 / n
    cvfit.beta <- lapply(1:n, FUN=function(i) {
	k <- i:n
	w <- epln * (k - 1) / (n - k + 1) + n.inv
	if (i > 1) {
	    w <- c(rep((n.inv - epln), times = (i - 1)), w)
	}
 *
 *
 */

double LARSW::TrainMaxWeight(arma::vec& beta, double epsilon, int nfolds, int nlambdas, double lambda_max_init, double lambda_max_ratio) {
	arma::ivec folds_ind = arma::randi(y.n_elem, arma::distr_param(0, nfolds - 1));
	arma::uvec i_sub;
	arma::vec wyi;
	arma::mat wXi;
	arma::vec predictions;
	arma::vec betaOpt, betaLambda;
	arma::vec wy;
	arma::mat wX;
	int n = y.n_elem, skip;
	double n1 = 1/((double)n);
	double epln = epsilon * n1;
	double lambda_max = (lambda_max_init > 0) ? lambda_max_init : n1 * max(abs(trans(X)*y));
	if (PRINT_TRACE > 1) {
		std::cout << "n1 = " << n1 << " lambda max = " << lambda_max << std::endl;
	}
	double lambda_min = lambda_max / lambda_max_ratio;
	if (PRINT_TRACE > 1) {
		std::cout << "ratio = " << lambda_max_ratio << " lambda min = " << lambda_min << std::endl;
	}
	double lambda_inc = pow(lambda_min / lambda_max, 1.0 / ((double)nlambdas));
	double k = 0, k_inc = 1 / ((double)nlambdas);
	double lambda;
	double lambda1, lambda2;
	double er;
	arma::vec weights(n, 1);
	double result_min, result_max, lambda_result_min, result_mean, lambda_result;
	result_max = -DBL_MAX;
	for (size_t I = 0; I < n; ++I) {
		for (size_t J = 0; J < I; ++J) {
			weights(J) = n1 - epln;
		}
		for (size_t J = I; J < n; ++J) {
			weights(J) = epln * J / ((double)(n - J)) + n1;
		}
		wy = sqrt(weights) % y;
		wX = diagmat(sqrt(weights)) * X;
		result_min = DBL_MAX;
		lambda = lambda_max;
		for (size_t i = 0; i < nlambdas; ++i) {
			if (PRINT_TRACE > 1) {
				std::cout << "I = " << I << " i = " << i << " lambda = " << lambda;
			}
			lambda1 = lambda * (1 - alpha);
			lambda2 = 0.5 * lambda * alpha;
			LARS lars(useCholesky, lambda1, lambda2);
			result_mean = 0;
			skip = 0;
			for (size_t j = 0; j < nfolds; ++j) {
				i_sub = arma::find(folds_ind != j);
				wyi = wy(i_sub);
				wXi = wX.rows(i_sub);
				lars.Train(wXi, wyi, betaOpt, transposeData);
				if (betaOpt.n_elem != X.n_cols) {
					skip++;
				} else {
					i_sub = arma::find(folds_ind == j);
					wyi = y(i_sub);
					wXi = X.rows(i_sub);
					predictions = wXi * betaOpt;
					er = norm(predictions - wyi);
					result_mean += er;
				}
				if (PRINT_TRACE > 1) {
					std::cout << " er[" << j << "] = " << er;
				}
			}
			if (nfolds - skip) {
				result_mean /= (nfolds - skip);
			} else {
				result_mean = DBL_MAX;
			}
			if (result_min > result_mean) {
				result_min = result_mean;
				lambda_result_min = lambda;
				betaLambda = betaOpt;
			}
			if (PRINT_TRACE > 1) {
				std::cout << std::endl;
			}
			k += k_inc;
			lambda *= lambda_inc;
		}
		if (result_max < result_min) {
			result_max = result_min;
			lambda_result = lambda_result_min;
			beta = betaLambda;
			opt_weights = weights;
		}
	}
	return result_max;
}

double LARSW::TrainRngWeight(arma::vec& beta, double epsilon, int nweights, int nfolds, int nlambdas, double lambda_max_init, double lambda_max_ratio) {
	arma::ivec folds_ind = arma::randi(y.n_elem, arma::distr_param(0, nfolds - 1));
	arma::uvec i_sub;
	arma::vec wyi;
	arma::mat wXi;
	arma::vec predictions;
	arma::vec betaOpt, betaLambda;
	arma::vec wy;
	arma::mat wX;
	int n = y.n_elem, skip;
	double n1 = 1/((double)n);
	double epln = epsilon * n1;
	double lambda_max = (lambda_max_init > 0) ? lambda_max_init : n1 * max(abs(trans(X)*y));
	if (PRINT_TRACE > 1) {
		std::cout << "n1 = " << n1 << " lambda max = " << lambda_max << std::endl;
	}
	double lambda_min = lambda_max / lambda_max_ratio;
	if (PRINT_TRACE > 1) {
		std::cout << "ratio = " << lambda_max_ratio << " lambda min = " << lambda_min << std::endl;
	}
	double lambda_inc = pow(lambda_min / lambda_max, 1.0 / ((double)nlambdas));
	double k = 0, k_inc = 1 / ((double)nlambdas);
	double lambda;
	double lambda1, lambda2;
	double er;
	arma::vec weights(n, 1);
	arma::uvec uw = arma::randi<arma::uvec> (nweights, arma::distr_param(0, n));
	double result_min, result_max, lambda_result_min, result_mean, lambda_result;
	result_max = -DBL_MAX;
	for (size_t II = 0; II < nweights; ++II) {
		size_t I = uw(II);
		for (size_t J = 0; J < I; ++J) {
			weights(J) = n1 - epln;
		}
		for (size_t J = I; J < n; ++J) {
			weights(J) = epln * J / ((double)(n - J)) + n1;
		}
		wy = sqrt(weights) % y;
		wX = diagmat(sqrt(weights)) * X;
		result_min = DBL_MAX;
		lambda = lambda_max;
		for (size_t i = 0; i < nlambdas; ++i) {
			if (PRINT_TRACE > 1) {
				std::cout << "I = " << I << " i = " << i << " lambda = " << lambda;
			}
			lambda1 = lambda * (1 - alpha);
			lambda2 = 0.5 * lambda * alpha;
			LARS lars(useCholesky, lambda1, lambda2);
			result_mean = 0;
			skip = 0;
			for (size_t j = 0; j < nfolds; ++j) {
				i_sub = arma::find(folds_ind != j);
				wyi = wy(i_sub);
				wXi = wX.rows(i_sub);
				lars.Train(wXi, wyi, betaOpt, transposeData);
				if (betaOpt.n_elem != X.n_cols) {
					skip++;
				} else {
					i_sub = arma::find(folds_ind == j);
					wyi = y(i_sub);
					wXi = X.rows(i_sub);
					predictions = wXi * betaOpt;
					er = norm(predictions - wyi);
					result_mean += er;
				}
				if (PRINT_TRACE > 1) {
					std::cout << " er[" << j << "] = " << er;
				}
			}
			if (nfolds - skip) {
				result_mean /= (nfolds - skip);
			} else {
				result_mean = DBL_MAX;
			}
			if (result_min > result_mean) {
				result_min = result_mean;
				lambda_result_min = lambda;
				betaLambda = betaOpt;
			}
			if (PRINT_TRACE > 1) {
				std::cout << std::endl;
			}
			k += k_inc;
			lambda *= lambda_inc;
		}
		if (result_max < result_min) {
			result_max = result_min;
			lambda_result = lambda_result_min;
			beta = betaLambda;
			opt_weights = weights;
		}
	}
	return result_max;
}

double LARSW::getError() {
	double err = 0;

	return err;
}

arma::mat std2armat(std::vector<std::vector<double> > &vec, int n_rows, int n_cols, int offset) {
	arma::mat X(n_rows, n_cols);
	for (size_t i = 0; i < n_rows; ++i) {
			for (size_t j = 0; j < n_cols; ++j) {
				if (j < offset) {
        			X(i, j) = vec[i][j];
				} else if (j > offset) {
					X(i, j) = vec[i][j + 1];
				}
			}
	}
	return X;
}

arma::mat std2arvec(std::vector<std::vector<double> > &vec, int n_rows, int offset) {
	arma::vec Y(n_rows, 1);
	for (size_t i = 0; i < n_rows; ++i) {
       	Y(i) = vec[i][offset];
	}
	return Y;
}

arma::mat std2armat_full(std::vector<std::vector<double> > &vec, int n_rows, int n_cols) {
		arma::mat X(n_rows, n_cols);
		for (size_t i = 0; i < n_rows; ++i) {
				for (size_t j = 0; j < n_cols; ++j) {
	        		X(i, j) = vec[i][j];
				}
		}
		return X;
}
