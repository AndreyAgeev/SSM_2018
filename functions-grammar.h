/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * functions-grammar.h
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

#ifndef _FUNCTIONS_GRAMMAR_H_
#define _FUNCTIONS_GRAMMAR_H_

// based on github.com:pbharrin/Genetic-Prog.git

#include <string>
#include <vector>
using namespace std;

class GrammarNode
{
public:
	string label;
	int precendence;
	int numChildren;
	GrammarNode* children[2];
	int child_type[2];
	GrammarNode() { } // to be used later
	virtual ~GrammarNode() { } // to be used later
	virtual double eval(vector<double>& inVal) = 0;  //setting the 0 makes it a PURE
	virtual void coprint() = 0;
	virtual GrammarNode* clone() = 0; //make a deep copy of the current tree
	virtual void setScale(vector<double>& a, vector<double>& b) = 0;
	virtual string getLabel() {	return label; }
// Print function. It's virtual, so only one operator << is needed.
	virtual void printOn (ostream& os) { os << label; }
};

// Print operator for all classes that inherit from GrammarNode
inline ostream& operator << (ostream& os, GrammarNode& gnd)
{
  gnd.printOn (os);
  return os;
}

//class for storing constant values
class ConstNode : public GrammarNode {
	double constVal;
	int constInd;
public:
	ConstNode();
	ConstNode(int preSetVal);
	ConstNode(int ind, double preSetVal);
	virtual double eval(vector<double>& inVal);
	ConstNode* clone();
	virtual void coprint();
	virtual void setScale(vector<double>& a, vector<double>& b) { }
//	virtual string getLabel();
//	virtual void printOn (ostream& os);
};

//class for using inputs
class InputNode : public GrammarNode {
	int inputIndex;
	string inpname;
	bool scalePresent;
	double scale;
	double offset;
public:
	InputNode ();
	InputNode(int inputId, string pname);
	virtual double eval(vector<double>& inVal);
	InputNode* clone();
	virtual void coprint();
	virtual void setScale(vector<double>& a, vector<double>& b) { 
		if (inputIndex >= 0) {
			scale = a[inputIndex]; offset = b[inputIndex]; scalePresent = true; 
		} else if (children[0]) {
			children[0]->setScale(a, b);
		}
	}
//	void setValues(int inIndex);
//	virtual string getLabel();
//	virtual void printOn (ostream& os);
};

//addition
class Add : public GrammarNode {
public:
	Add();
	virtual double eval(vector<double>& inVal);
	Add* clone();
	virtual void coprint();
	virtual void setScale(vector<double>& a, vector<double>& b) { 
		children[0]->setScale(a, b);
		children[1]->setScale(a, b);
	}
//	virtual string getLabel();
//	virtual void printOn (ostream& os);
};

//subtraction
class Subtract : public GrammarNode {
public:
	Subtract();
	virtual double eval(vector<double>& inVal);
	Subtract* clone();
	virtual void coprint();
	virtual void setScale(vector<double>& a, vector<double>& b) { 
		children[0]->setScale(a, b);
		children[1]->setScale(a, b);
	}
//	virtual string getLabel();
//	virtual void printOn (ostream& os);
};

//multiplication
class Multiply : public GrammarNode {
public:
	Multiply();
	virtual double eval(vector<double>& inVal);
	Multiply* clone();
	virtual void coprint();
	virtual void setScale(vector<double>& a, vector<double>& b) { 
		children[0]->setScale(a, b);
		children[1]->setScale(a, b);
	}
//	virtual string getLabel();
//	virtual void printOn (ostream& os);
};

//division
class Divide : public GrammarNode {
public:
	Divide();
	virtual double eval(vector<double>& inVal);
	Divide* clone();
	virtual void coprint();
	virtual void setScale(vector<double>& a, vector<double>& b) { 
		children[0]->setScale(a, b);
		children[1]->setScale(a, b);
	}
//	virtual string getLabel();
//	virtual void printOn (ostream& os);
};

//subtraction of const
class InputMinusConst : public GrammarNode {
public:
	InputMinusConst();
	virtual double eval(vector<double>& inVal);
	InputMinusConst* clone();
	virtual void coprint();
	virtual void setScale(vector<double>& a, vector<double>& b) { 
		children[0]->setScale(a, b);
		children[1]->setScale(a, b);
	}
//	virtual string getLabel();
//	virtual void printOn (ostream& os);
};

//rec subtraction of const
class RecInputMinusConst : public GrammarNode {
public:
	RecInputMinusConst();
	virtual double eval(vector<double>& inVal);
	RecInputMinusConst* clone();
	virtual void coprint();
	virtual void setScale(vector<double>& a, vector<double>& b) { 
		children[0]->setScale(a, b);
		children[1]->setScale(a, b);
	}
//	virtual string getLabel();
//	virtual void printOn (ostream& os);
};
	
class GrammarContainer {
public:	
	GrammarContainer (vector<string>& measurements, int n_t, int n_c = 3, bool m_c = false) : predictors(measurements), tree(n_t) {
		n_trees = n_t;
		n_predictors = predictors.size();
		n_nodes_type_0 = 7;
		n_nodes_type_1 = 1;
		n_nodes_type_2 = 1;
		n_cycles = n_c;
		make_corrections = m_c;
	}
	void build_nth_tree(vector<int>& genotype, vector<double>& conc, int n, double *phenotype);
	GrammarNode*get_nth_tree(int n) { return tree[n]; }
private:
	vector<string> predictors;
	vector<GrammarNode*> tree;
	int n_trees;
	int n_predictors;
	int n_nodes_type_0;
	int n_nodes_type_1;
	int n_nodes_type_2;
	int n_cycles;
	bool make_corrections;
	GrammarNode* build_tree(vector<int>& genotype, vector<double>& conc, double *phenotype);
	GrammarNode* find_node(int type, int gen, double conc, double& phenotype);
	GrammarNode* find_node_type_0(int gen, double& phenotype);
	GrammarNode* find_node_type_1(int gen, double conc, double& phenotype);
	GrammarNode* find_node_type_2(int gen, double& phenotype);
};

#endif // _FUNCTIONS_GRAMMAR_H_
