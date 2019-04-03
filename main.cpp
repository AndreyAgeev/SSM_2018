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
					QCoreApplication::translate("main", "Latitude.")},
					                {{"3", "VPDF"},
										QCoreApplication::translate("main", "VPDF.")},
										{{"R", "R"},
										QCoreApplication::translate("main", "function mode"),
										QCoreApplication::translate("main", "R")},
										{{"P", "P"},
										QCoreApplication::translate("main", "threshold"),
										QCoreApplication::translate("main", "P")},
									{{"i", "FixFind"},
										QCoreApplication::translate("main", "FixFind")},
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
									{{"f", "PDEN"},////////////////////////////
														QCoreApplication::translate("main", "Number of functions."),
														QCoreApplication::translate("main", "f")},
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
													{{"4", "EXTR"},
														QCoreApplication::translate("main", "EXTR.")},
													{{"j", "MAI1"},
														QCoreApplication::translate("main", "MAI1.")},
													{{"a", "MAI"},
														QCoreApplication::translate("main", "MAI.")},
													{{"o", "INSOL"},
														QCoreApplication::translate("main", "INSOL.")},
													{{"1", "U"},
														QCoreApplication::translate("main", "U.")},
		});
	parser.addPositionalArgument("cropsoinifile", "The file to open.");
	parser.addPositionalArgument("samplesfile", "The file to open.");
	parser.addPositionalArgument("weatherfile", "The file to open.");
	parser.addPositionalArgument("funcsfile", "The file to read funcs.");
	parser.process(a);
	const QStringList args = parser.positionalArguments();
	Model * model;
	Parametrs param;
	if (args.size())
	{
		param.func_file_name = args.at(3); // funcs
		param.h5_file_name = args.at(2); // samples
		param.h5_table_name = args.at(1); // weather
		param.crops_ini_file = args.at(0); // weather
	//	cout << "name file: " << param.func_file_name.toStdString() << " " << param.h5_file_name.toStdString() << " " << param.h5_table_name.toStdString() << param.crops_ini_file.toStdString() << param.crops_ini_file.toStdString() << endl;
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
			std::cout << "Bad nd: " + D << endl;
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/*	const QString LatParameter = parser.value("Latitude");
			const double Lat = LatParameter.toDouble();
			param.Latitude = Lat;
		
			const QString VPDFParameter = parser.value("VPDF");
			const int VPDF = VPDFParameter.toDouble();	
			param.VPDF = VPDF;

		const QString FixFindParameter = parser.value("FixFind");
		const int FixFind = FixFindParameter.toInt();
		param.FixFind = FixFind;


		const QString PdoyParameter = parser.value("Pdoy");
		const int Pdoy = PdoyParameter.toInt();
		param.Pdoy = Pdoy;

		const QString SearchDurParameter = parser.value("SearchDur");
		const double SearchDur = SearchDurParameter.toDouble();
		param.SearchDur = SearchDur;


		const QString SowWatParameter = parser.value("SowWat");
		const double SowWat = SowWatParameter.toDouble();
		param.SowWat = SowWat;


				const QString PDENParameter = parser.value("PDEN");
				const double PDEN = PDENParameter.toDouble();
				param.PDEN = PDEN;


				const QString waterParameter = parser.value("water");
				const int water = waterParameter.toInt();
				param.water = water;

				const QString IRGLVLParameter = parser.value("IRGLVL");
				const double IRGLVL = IRGLVLParameter.toDouble();
				param.IRGLVL = IRGLVL;


				const QString SOLDEPParameter = parser.value("SOLDEP");
				const double SOLDEP = SOLDEPParameter.toDouble();
				param.SOLDEP = SOLDEP;

				const QString DEP1Parameter = parser.value("DEP1");
				const double DEP1 = DEP1Parameter.toDouble();
				param.DEP1 = DEP1;

				const QString SALBParameter = parser.value("SALB");
				const double SALB = SALBParameter.toDouble();
				param.SALB = SALB;

				const QString CN2Parameter = parser.value("CN2");
				const double CN2 = CN2Parameter.toDouble();
				param.CN2 = CN2;

				const QString DRAINFParameter = parser.value("DRAINF");
				const double DRAINF = DRAINFParameter.toDouble();
				param.DRAINF = DRAINF;

				const QString SATParameter = parser.value("SAT");
				const double SAT = SATParameter.toDouble();
				param.SAT = SAT;


				const QString DULParameter = parser.value("DUL");
				const double DUL = DULParameter.toDouble();
				param.DUL = DUL;

				const QString EXTRParameter = parser.value("EXTR");
				const double EXTR = EXTRParameter.toDouble();
				param.EXTR = EXTR;

				const QString MAI1Parameter = parser.value("MAI1");
				const int MAI1 = MAI1Parameter.toInt();
				param.MAI1 = MAI1;


				const QString MAIParameter = parser.value("MAI");
				const int MAI = MAIParameter.toInt();
				param.MAI = MAI;

				const QString INSOLParameter = parser.value("INSOL");
				const int INSOL = INSOLParameter.toInt();
				param.INSOL = INSOL;


				const QString UParameter = parser.value("U");
				const double U = UParameter.toDouble();
				param.U = U;

				*/
						
				//////////////////////////////////////////////////////////////////
		param.Latitude = 36.41;
		param.VPDF = 0.75;
		param.FixFind = 1;
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
