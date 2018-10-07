#pragma once
#include <QStandardItemModel>
#include <QtSql/qsqldatabase.h>
#include <QtSql/qsqlquery.h>
using namespace std;
class Data
{
public:
	struct Data_h5 
	{
			
	}; 
	struct Data_xls 
	{
            vector<int> years;
			vector<int> doy;
			vector<double> srad;
			vector<double> tmax;
			vector<double> tmin;
			vector<double> rain;
	};
	Data_h5 data_h5;
	Data_xls data_xls;
	Data() {}
	void parse_param(QString line)
	{
		QStringList list;
    	list = line.split(QRegExp("\\W+"));//WRONG separation!!!
		for (int i = 0; i < list.size(); i++)
		{
			cout << list.at(i).toLocal8Bit().constData() << endl;
			string elem = list.at(i).toLocal8Bit().constData();
			if(i == 0)
			    data_xls.years.push_back(atoi(elem.c_str()));
			else if(i == 1)
				data_xls.doy.push_back(atoi(elem.c_str()));
			else if (i == 2)
				data_xls.srad.push_back(atof(elem.c_str()));
			else if (i == 3)
				data_xls.tmax.push_back(atof(elem.c_str()));
			else if (i == 4)
				data_xls.tmin.push_back(atof(elem.c_str()));
			else if (i == 5)
				data_xls.rain.push_back(atof(elem.c_str()));
		}
	} 
	int read_xls(QString file_name)
	{
		int index = 0;
		QStandardItemModel * csvModel = new QStandardItemModel();
		QFile file(file_name);
		if (!file.open(QFile::ReadOnly | QFile::Text)) {
			qDebug() << "File not exists";
		}
		else {
			// Create a thread to retrieve data from a file
			QTextStream in(&file);
			//Reads the data up to the end of file
			while (!in.atEnd())
			{
				QString line = in.readLine();
				parse_param(line);
			}
			file.close();
		}
		return 1;
	}
};