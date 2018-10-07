#pragma once

class Parametrs
{
	
public:
	enum TypeFile
	{
		SQL = 0,
		H5 = 1
	};
	QString h5_file_name;
	QString sql_file_name;
	int type;
	Parametrs(QString file_name, TypeFile type)
	{
		if (type == SQL)
			sql_file_name = file_name;
		else
			h5_file_name = file_name;
	}
};