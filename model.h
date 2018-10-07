#pragma once
#include <string>
#include <iostream>
#include <QtCore>
#include <QDebug>

#include <highfive/H5Attribute.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include "data.h"
#include "parametrs.h"
using namespace std;
using namespace HighFive;
class Model: public QObject
{
	Q_OBJECT
public:
	Data data;
	Parametrs param;
	explicit Model(QString file_name, Parametrs::TypeFile type, QObject *parent = 0) : QObject(parent), param(file_name, type)
	{
		if (type == Parametrs::TypeFile::H5)
			run_h5();
		else if (type == Parametrs::TypeFile::SQL)
			run_sql();
		else
			return;
	}
public slots:
	void run_sql()
	{
		data.read_xls(param.sql_file_name);
		//TO DO/////////////////////////////////////////////////////////////////////
		for (int x = 0; x < 4; x++)//to do: x - params ;
		{
			for (int y = 0; y < 4; y++)
			{

			}
		}
	}
	void run_h5()
	{
	}
};