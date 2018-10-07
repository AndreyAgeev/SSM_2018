#include <QtCore/QCoreApplication>
#include "model.h"
int main(int argc, char *argv[])
{
	QCoreApplication a(argc, argv);
	QCommandLineParser parser;
	parser.addHelpOption();
	parser.addVersionOption();
	parser.addOptions({ { { "f", "parametr file"} ,
		QCoreApplication::translate("main", "Type of file(sql or h5)."),
		QCoreApplication::translate("main", "F") } });
	parser.addPositionalArgument("file", "The file to open.");
	parser.process(a);
	//TO DO - Add all parameters
	QStringList arqs = parser.positionalArguments();
	Model * model;
	if (arqs.size())
	{
		QString file_name = arqs.at(0);
		if (parser.isSet("f"))
		{
			Parametrs::TypeFile type;
			const QString type_str = parser.value("parametr file");
			const int type_int = type_str.toInt();
			cout << type_int;
			if (type_int == 0)
				type = Parametrs::TypeFile::SQL;
			else if (type_int == 1)
				type = Parametrs::TypeFile::H5;
			else
				return 0;//give an error or make a default
			model = new Model(file_name, type);
		}
		
	}
	return a.exec();
}
