/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * functions-grammar.cc
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

// based on github.com:pbharrin/Genetic-Prog.git

#include <iostream>
#include <math.h>
#include <QStack>
#include "functions-grammar.h"

#ifndef ARGVAL
#define ARGVAL 0.00390625 /* 1/256 */
#endif

using namespace std;

double int_to_double(int arg) {
/*	int i_arg = ( arg < 22 ) ? (int)floor(exp(arg)) : (int)floor(exp(arg % 21));
    int left;
	unsigned short left16;
	unsigned short right16;
	right16 = (unsigned short)i_arg;
	left = i_arg >> (8 * sizeof(unsigned char));
	left16 = (unsigned short)left;
	double merged = (double)left16/(double)right16;*/
	double argval = ARGVAL;
	double merged = argval * arg;
	return merged;
}

//
//		ConstNode
//
ConstNode::ConstNode() {
	numChildren = 0;
	constVal = rand()/(double)RAND_MAX;
	char str[20] = "";
	sprintf(str, "Const: %f", constVal);
	label = str;
	precendence = -1;
}

ConstNode::ConstNode(int preSetVal) {
	numChildren = 0;
	constVal = int_to_double(preSetVal);
	char str[20] = "";
	sprintf(str, "Const: %f", constVal);
	label = str;
	precendence = -1;
}

ConstNode::ConstNode(int ind, double preSetVal) {
	numChildren = 0;
	constVal = preSetVal;
	constInd = ind;
	char str[20] = "";
	sprintf(str, "Const: %f", constVal);
	label = str;
	precendence = -1;
}

double ConstNode::eval(vector<double>& inVal) {
	return constVal;
}

void ConstNode::coprint() {
	cout << constVal;
}

ConstNode* ConstNode::clone() {
	ConstNode* retTree = new ConstNode(constVal);
	return retTree; 
}

//
//		InputNode
//
InputNode::InputNode(int inputId, string pname) {
	numChildren = 0;
	scalePresent = false;
	inputIndex = inputId;
	inpname = pname;
	char str[256] = "";
	sprintf(str, "%s[%d]", pname.c_str(), inputIndex);
//	setValues(inputIndex);
	label = str;
	precendence = -1;
}

InputNode::InputNode() {
	numChildren = 1;
	scalePresent = false;
	inputIndex = -1;
	child_type[0] = 2;
	precendence = -1;
}

double InputNode::eval(vector<double>& inVal){
	if (inputIndex >= 0) {
		return inVal[inputIndex];
	} else if (children[0]) {
		return children[0]->eval(inVal);
	} else {
		cerr << "not defined in input" << endl;
	}
}

void InputNode::coprint(){
	if (inputIndex >= 0) {
		if (scalePresent) {
			cout << "((" << inpname << "-" << offset << ")/" << scale << ")";
		} else {
			cout << inpname;
		}
	} else if (children[0]) {
		children[0]->coprint();
	} else {
		cerr << "not defined in input coprint" << endl;
	}
}

//void InputNode::setValues(int inIndex){
//	char str[20] = "";
//	sprintf(str, "InputVal: %d", inIndex);
//	label = str;
//}

InputNode* InputNode::clone(){
	InputNode* retTree = new InputNode(inputIndex, inpname);
	return retTree; 
}

//
//		Add
//
Add::Add(){
	numChildren = 2;
	child_type[0] = 0;
	child_type[1] = 0;
	label = "Add";
	precendence = 6;
}

double Add::eval(vector<double>& inVal){
	if (children[0] && children[1]){
		return children[0]->eval(inVal) + children[1]->eval(inVal);
	}
	else {
		cerr << "left and right not defined in add"<<endl;
		return -1.0;
	}
}

void Add::coprint(){
	if (children[0] && children[1]){
		cout << " (";
		children[0]->coprint();
		cout << " + ";
		children[1]->coprint();
		cout << ") ";
	}
	else {
		cerr << "left and right not defined in add"<<endl;
	}
}

Add* Add::clone(){
	Add* retNode = new Add();
	for (int i=0; i<numChildren; i++) {
		retNode->children[i] = children[i]->clone();
	}
	return retNode;
}

//
//		Subtract
//
Subtract::Subtract(){
	numChildren = 2;
	child_type[0] = 0;
	child_type[1] = 0;
	label = "Subtract";
	precendence = 6;
}

double Subtract::eval(vector<double>& inVal){
	if (children[0] && children[1]){
		return children[0]->eval(inVal) - children[1]->eval(inVal);
	}
	else {
		cerr << "left and right not defined in subtract"<<endl;
		return -1.0;
	}
}

void Subtract::coprint(){
	if (children[0] && children[1]){
		cout << " (";
		children[0]->coprint();
		cout << " - ";
		children[1]->coprint();
		cout << ") ";
	}
	else {
		cerr << "left and right not defined in subtract"<<endl;
	}
}

Subtract* Subtract::clone(){
	Subtract* retNode = new Subtract();
	for (int i=0; i<numChildren; i++) {
		retNode->children[i] = children[i]->clone();
	}
	return retNode;
}

//
//		Multiply
//
Multiply::Multiply(){
	numChildren = 2;
	child_type[0] = 0;
	child_type[1] = 0;
	label = "Multiply";
	precendence = 5;
}

double Multiply::eval(vector<double>& inVal){
	if (children[0] && children[1]){
		return children[0]->eval(inVal) * children[1]->eval(inVal);
	}
	else {
		cerr << "left and right not defined in multiply"<<endl;
		return -1.0;
	}
}

void Multiply::coprint(){
	if (children[0] && children[1]){
		cout << " (";
		children[0]->coprint();
		cout << " * ";
		children[1]->coprint();
		cout << ") ";
	}
	else {
		cerr << "left and right not defined in multiply"<<endl;
	}
}

Multiply* Multiply::clone(){
	Multiply* retNode = new Multiply();
	for (int i=0; i<numChildren; i++) {
		retNode->children[i] = children[i]->clone();
	}
	return retNode;
}

//
//		Divide
//
Divide::Divide() {
	numChildren = 2;
	child_type[0] = 0;
	child_type[1] = 0;
	label = "Divide";
	precendence = 5;
}

double Divide::eval(vector<double>& inVal){
	if (children[0] && children[1]){
		return children[0]->eval(inVal) / children[1]->eval(inVal);
	}
	else {
		cerr << "left and right not defined in divide"<<endl;
		return -1.0;
	}
}

void Divide::coprint(){
	if (children[0] && children[1]){
		cout << " (";
		children[0]->coprint();
		cout << " / ";
		children[1]->coprint();
		cout << ") ";
	}
	else {
		cerr << "left and right not defined in divide"<<endl;
	}
}

Divide* Divide::clone(){
	Divide* retNode = new Divide;
	for (int i=0; i<numChildren; i++) {
		retNode->children[i] = children[i]->clone();
	}
	return retNode;
}

//
//		InputMinusConst
//
InputMinusConst::InputMinusConst(){
	numChildren = 2;
	child_type[0] = 2;
	child_type[1] = 1;
	label = "InputMinusConst";
	precendence = 6;
}

double InputMinusConst::eval(vector<double>& inVal){
	if (children[0] && children[1]){
		return children[0]->eval(inVal) - children[1]->eval(inVal);
	}
	else {
		cerr << "left and right not defined in InputMinusConst" << endl;
		return -1.0;
	}
}

void InputMinusConst::coprint(){
	if (children[0] && children[1]){
		cout << " (";
		children[0]->coprint();
		cout << " - ";
		children[1]->coprint();
		cout << ") ";
	}
	else {
		cerr << "left and right not defined in InputMinusConst" << endl;
	}
}

InputMinusConst* InputMinusConst::clone(){
	InputMinusConst* retNode = new InputMinusConst();
	for (int i=0; i<numChildren; i++) {
		retNode->children[i] = children[i]->clone();
	}
	return retNode;
}

//
//		RecInputMinusConst
//
RecInputMinusConst::RecInputMinusConst(){
	numChildren = 2;
	child_type[0] = 2;
	child_type[1] = 1;
	label = "RecInputMinusConst";
	precendence = 6;
}

double RecInputMinusConst::eval(vector<double>& inVal){
	if (children[0] && children[1]){
		return 1/(children[0]->eval(inVal) - children[1]->eval(inVal));
	}
	else {
		cerr << "left and right not defined in RecInputMinusConst" << endl;
		return -1.0;
	}
}

void RecInputMinusConst::coprint(){
	if (children[0] && children[1]){
		cout << " (1/(";
		children[0]->coprint();
		cout << " - ";
		children[1]->coprint();
		cout << ")) ";
	}
	else {
		cerr << "left and right not defined in InputMinusConst" << endl;
	}
}

RecInputMinusConst* RecInputMinusConst::clone(){
	RecInputMinusConst* retNode = new RecInputMinusConst();
	for (int i=0; i<numChildren; i++) {
		retNode->children[i] = children[i]->clone();
	}
	return retNode;
}

/*
 * expr = grule(op(expr, expr), X, (X - Y), 1/(X - Y)),
 * Y = grule(0.1, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0),
 * op = grule(`+`, `-`, `*`, `/`),
 * X = grule(@PREDICTORS@)

 * expr = grule((expr + expr), (expr - expr), (expr * expr), (expr / expr), X, (X - Y), 1/(X - Y)),
 * Y = grule(0.1, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0),
 * X = grule(@PREDICTORS@)

 
 */


GrammarNode* GrammarContainer::find_node_type_0(int gen, double& phenotype)
{
	GrammarNode*n;
	int ge = gen % n_nodes_type_0;
	switch (ge) {
		case 0:
			n = new Add();
			break;
		case 1:
			n = new Subtract();
			break;
		case 2:
			n = new Multiply();
			break;
		case 3:
			n = new Divide();
			break;
		case 4:
			n = new InputNode();
			break;
		case 5:
			n = new InputMinusConst();
			break;
		case 6:
			n = new RecInputMinusConst();
			break;
		default:
			cout << "Error find node type 0" << gen << endl;
	}
	phenotype = (double)ge;
	return n;
}

GrammarNode* GrammarContainer::find_node_type_1(int gen, double conc, double& phenotype)
{
	GrammarNode*n;
	n = new ConstNode(gen, conc);
	phenotype = conc;
	return n;
}

GrammarNode* GrammarContainer::find_node_type_2(int gen, double& phenotype)
{
	GrammarNode*n;
	int ge = gen % n_predictors;
	n = new InputNode(ge, predictors[ge]);
	phenotype = (double)ge;
	return n;
}

GrammarNode* GrammarContainer::find_node(int type, int gen, double conc, double& phenotype)
{
	GrammarNode*n;
	switch (type) {
		case 0:
			n = find_node_type_0(gen, phenotype);
			break;
		case 1:
			n = find_node_type_1(gen, conc, phenotype);
			break;
		case 2:
			n = find_node_type_2(gen, phenotype);
			break;
		default:
			cout << "Error find node " << type << endl;
	}
	return n;
}

/* type:
 * 0 - expr
 * 1 - const
 * 2 - predictor
 * (3 - operation?)
 */

GrammarNode* GrammarContainer::build_tree(vector<int>& genotype, vector<double>& conc, double *phenotype)
{
	struct StackNode {
		GrammarNode*kernel;
		int side;
		int child_type;
		StackNode(GrammarNode*orig, int p, int ct) {
			kernel = orig;
			side = p;
			child_type = ct;
		}
	};
	GrammarNode*root_expr, *child;
	StackNode *curr;
	QStack<StackNode*> stack;

	int i = 0;
	int curr_type = 0;
	int cycle = 0;
	
	root_expr = find_node(curr_type, genotype[i], conc[i], phenotype[i]);
	i++;
//	cout << "g " << genotype[i - 1] << ' ' << i - 1 << ' ' << *root_expr << endl;

	if (root_expr->numChildren == 2) {
		stack.push(new StackNode(root_expr, 1, root_expr->child_type[1]));
	}
	if (root_expr->numChildren > 0) {
		stack.push(new StackNode(root_expr, 0, root_expr->child_type[0]));
	}
	i = (i < genotype.size()) ? i : 0 ; 
	while(!stack.isEmpty() && cycle < n_cycles) {
		curr = stack.pop();
		curr_type = curr->child_type;
		child = find_node(curr_type, genotype[i], conc[i], phenotype[i]);
		i++;
//		cout << "g " << genotype[i - 1] << ' ' << i - 1 << ' ' << *child << endl;
		curr->kernel->children[curr->side] = child;
		if (child->numChildren == 2) {
			stack.push(new StackNode(child, 1, child->child_type[1]));
		}
		if (child->numChildren > 0) {
			stack.push(new StackNode(child, 0, child->child_type[0]));
		}
		if (i >= genotype.size()) { 
			i = 0;
			cycle++;
		}
	}
/**/
	while(make_corrections == true && !stack.isEmpty()) {
		curr = stack.pop();
		curr->kernel->children[curr->side] = new ConstNode(-1, curr->side);
	}/**/
/* This should be an exception probably */
	if (!stack.isEmpty()) {
		cout << "error=" << 1.76e16 << endl;
		exit(-100);
	}
	return root_expr;
}

void GrammarContainer::build_nth_tree(vector<int>& genotype, vector<double>& conc, int n, double *phenotype)
{
	tree[n] = build_tree(genotype, conc, phenotype);
}
