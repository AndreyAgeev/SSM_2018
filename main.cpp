#include <QtCore/QCoreApplication>
#include "model.h"
#define add(name) param.## name
int main(int argc, char *argv[])
{

	QCoreApplication a(argc, argv);
	QCommandLineParser parser;
	parser.addHelpOption();
	parser.addVersionOption();
	parser.addOptions({
	{{"l", "Latitude"},
		QCoreApplication::translate("main", "Latitude."),
		QCoreApplication::translate("main", "l")},
	{{"s", "VPDF"},
		QCoreApplication::translate("main", "VPDF."),
		QCoreApplication::translate("main", "s")},
	{{"y", "yno"},
		QCoreApplication::translate("main", "Number of years."),
		QCoreApplication::translate("main", "Y")},
	{{"R", "R"},
		QCoreApplication::translate("main", "function mode"),
		QCoreApplication::translate("main", "R")},
	{{"P", "P"},
		QCoreApplication::translate("main", "threshold"),
		QCoreApplication::translate("main", "P")},
	{{"f", "FirstYear"},
		QCoreApplication::translate("main", "FirstYear."),
		QCoreApplication::translate("main", "F")},
	{{"i", "FixFind"},
		QCoreApplication::translate("main", "FixFind")},
	{{"p", "Pdoy"},
		QCoreApplication::translate("main", "Pdoy."),
		QCoreApplication::translate("main", "P")},
	{{"T", "print-trace"},
		QCoreApplication::translate("main", "Print trace."),
		QCoreApplication::translate("main", "T")},
	{{"D", "num-days"},
		QCoreApplication::translate("main", "Number of days."),
		QCoreApplication::translate("main", "D")},
	{{"N", "number-of-funcs"},
		QCoreApplication::translate("main", "Number of functions."),
		QCoreApplication::translate("main", "N")},
	{{"L", "length-of-word"},
		QCoreApplication::translate("main", "Length of a word to represent one function."),
		QCoreApplication::translate("main", "L")},
	{{"Q", "extra_covar"},
		QCoreApplication::translate("main", "Read binary covariates.")},
	{{"d", "SearchDur"},
		QCoreApplication::translate("main", "SearchDur.")},
	{{"w", "SowWat"},
		QCoreApplication::translate("main", "SowWat."),
		QCoreApplication::translate("main", "M")},
	{{"n", "PDEN"},////////////////////////////
		QCoreApplication::translate("main", "Number of functions."),
		QCoreApplication::translate("main", "N")},
	{{"m", "water"},
		QCoreApplication::translate("main", "Make corrections.")},
	{{"r", "IRGLVL"},
		QCoreApplication::translate("main", "IRGLVL."),
		QCoreApplication::translate("main", "L")},
	{{"t", "SOLDEP"},
		QCoreApplication::translate("main", "SOLDEP."),
		QCoreApplication::translate("main", "P")},
	{{"q", "DEP1"},
		QCoreApplication::translate("main", "DEP1")},
	{{"z", "SALB"},
		QCoreApplication::translate("main", "SALB.")},
	{{"x", "CN2"},
		QCoreApplication::translate("main", "CN2.")},
	{{"c", "DRAINF"},
		QCoreApplication::translate("main", "DRAINF.")},
	{{"b", "SAT"},
		QCoreApplication::translate("main", "SAT.")},
	{{"k", "DUL"},
		QCoreApplication::translate("main", "DUL.")},
	{{"E", "EXTR"},
		QCoreApplication::translate("main", "EXTR.")},
	{{"j", "MAI1"},
		QCoreApplication::translate("main", "MAI1.")},
	{{"a", "MAI"},
		QCoreApplication::translate("main", "MAI.")},
	{{"o", "INSOL"},
		QCoreApplication::translate("main", "INSOL.")},
	{{"U", "U"},
		QCoreApplication::translate("main", "U.")},
	});
//	parser.addPositionalArgument("cropsoinifile", "The file to open.");
	parser.addPositionalArgument("samplesfile", "The file to open.");
	parser.addPositionalArgument("weatherfile", "The file to open.");
//	parser.addPositionalArgument("funcsfile", "The file to read funcs.");
	parser.addPositionalArgument("cropsoinifile", "The file to open.");
	parser.process(a);
	const QStringList args = parser.positionalArguments();
	Model * model;
	Parametrs param;
	if (args.size())
	{
//		param.func_file_name = args.at(3); // funcs
		param.crops_ini_file = args.at(2); // crops
		param.h5_file_name = args.at(1); // samples
		param.h5_table_name = args.at(0); // weather
//		param.crops_ini_file = args.at(0); // crops
//		cout << "name file: " << " " << param.h5_file_name.toStdString() << " " << param.h5_table_name.toStdString() << " " << param.crops_ini_file.toStdString() << endl;
		const QString TParameter = parser.value("print-trace");
		const int T = TParameter.toInt();
		if (T < 0) {
			std::cout << "Bad p: " + T;
		}

		const QString RParameter = parser.value("R");
		const int R = RParameter.toInt();
		if (R != 0 && R != 1) {
			std::cout << "Bad p: " + R;
		}

		const QString PParameter = parser.value("P");
		const int P = PParameter.toInt();
		if (P != 0 && P != 1 && P != 3 && P != 5 && P != 7 && P != 8 && P != -1) {
			std::cout << "Bad p: " + P;
		}

		const QString NParameter = parser.value("number-of-funcs");
		const int N = NParameter.toInt();
		if (N < 0) {
			std::cout << "Bad nf: " + N;
		}
		const QString LParameter = parser.value("length-of-word");
		const int L = LParameter.toInt();
		if (L < 0) {
			std::cout << "Bad wl: " + L;
		}
		const QString DParameter = parser.value("num-days");
		const int D = DParameter.toInt();
		if (D <= 0) {
			std::cout << "Bad nd: " + D;
		}
		const QString nfParameter2 = parser.value("Latitude");
		param.Latitude = nfParameter2.toDouble();
		const QString nfParameter3 = parser.value("VPDF");
		param.VPDF = nfParameter3.toDouble();
		const QString nfParameter4 = parser.value("yno");
		param.yno = nfParameter4.toDouble();
		const QString nfParameter5 = parser.value("FirstYear");
		param.FirstYear = nfParameter5.toInt();
		const QString nfParameter6 = parser.value("FixFind");
		param.FixFind = nfParameter6.toInt();
		const QString nfParameter7 = parser.value("Pdoy");
		param.Pdoy = nfParameter7.toInt();
		/*		const QString nfParameter8 = parser.value("SearchDur");
				param.SearchDur = nfParameter8.toInt();
				const QString nfParameter9 = parser.value("SowWat");
				param.SowWat = nfParameter9.toInt();
				const QString nfParameter10 = parser.value("PDEN");
				param.PDEN = nfParameter10.toInt();
				const QString nfParameter11 = parser.value("water");
				param.water = nfParameter11.toInt();
				const QString nfParameter12 = parser.value("IRGLVL");
				param.IRGLVL = nfParameter12.toInt();
				const QString nfParameter13 = parser.value("SOLDEP");
				param.SOLDEP = nfParameter13.toDouble();
				const QString nfParameter14 = parser.value("DEP1");
				param.DEP1 = nfParameter14.toInt();
				const QString nfParameter15 = parser.value("SALB");
				param.SALB = nfParameter15.toDouble();
				const QString nfParameter16 = parser.value("CN2");
				param.CN2 = nfParameter16.toInt();
				const QString nfParameter17 = parser.value("DRAINF");
				param.DRAINF = nfParameter17.toDouble();
				const QString nfParameter18 = parser.value("SAT");
				param.SAT = nfParameter18.toDouble();
				const QString nfParameter = parser.value("DUL");
				param.DUL = nfParameter.toDouble();
				const QString nfParameter19 = parser.value("EXTR");
				param.EXTR = nfParameter19.toDouble();
				const QString nfParameter20 = parser.value("MAI1");
				param.MAI1 = nfParameter20.toInt();
				const QString nfParameter21 = parser.value("MAI");
				param.MAI = nfParameter21.toInt();
				const QString nfParameter22 = parser.value("INSOL");
				param.INSOL = nfParameter22.toInt();
				const QString nfParameter23 = parser.value("U");
				param.U = nfParameter23.toInt();
				const QString nfParameter24 = parser.value("CropColNo");
				param.CropColNo = nfParameter24.toInt();*/
				//		param.Print();
						////////////////////////////////////////////////////////
//		param.Latitude = 36.41;
//		param.VPDF = 0.75;
		param.yno = 10;
		param.FirstYear = 2014;
		param.FixFind = 1;
		param.Pdoy = 90;
		param.SearchDur = 60;
		param.SowWat = 15;
		param.PDEN = 15;
		param.water = 2;
		param.IRGLVL = 0.5;
		param.SOLDEP = 900;
		param.DEP1 = 200;
		param.SALB = 0.13;
		param.CN2 = 79;
		param.DRAINF = 0.5;
		param.SAT = 0.36;
		param.DUL = 0.264;
		param.EXTR = 0.13;
		param.MAI1 = 1;
		param.MAI = 1;
		param.INSOL = 3;
		param.U = 6;
		param.ROW = 7;
		param.file_mode = false;
		param.nF = N;
		param.wL = L;
		param.nD = D;
		param.rT = 1;
		param.ecovar = parser.isSet("extra_covar");
		param.print_trace = T;
		param.function_mode = R;
		param.threshold = P;
	}
	if (param.print_trace > 0) {
		param.Print();
	}
	model = new Model(param, &a);
	QObject::connect(model, SIGNAL(finished()), &a, SLOT(quit()));
	QTimer::singleShot(0, model, SLOT(run_h5()));
	return a.exec();
}
