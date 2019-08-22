#include <QtCore/QCoreApplication>
#include "model.h"
#include <WinSDKVer.h>
#define _WIN32_WINNT 0x0601
#include <SDKDDKVer.h>

#define add(name) param.## name
int main(int argc, char *argv[])
{

	QCoreApplication a(argc, argv);
	QCommandLineParser parser;
	parser.addHelpOption();
	parser.addVersionOption();
	parser.addOptions({
	{{"l", "Latitude"},
	    QCoreApplication::translate("main", "Latitude"),
	    QCoreApplication::translate("main", "l") },
	{{"s", "VPDF"},
		QCoreApplication::translate("main", "VPDF."),
		QCoreApplication::translate("main", "s") },
	{{"R", "R"},
		QCoreApplication::translate("main", "function mode"),
		QCoreApplication::translate("main", "R")},
	{{"P", "P"},
		QCoreApplication::translate("main", "threshold"),
		QCoreApplication::translate("main", "P")},
	{{"i", "FixFind"},
		QCoreApplication::translate("main", "FixFind"),
	    QCoreApplication::translate("main", "i")},
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
		QCoreApplication::translate("main", "SearchDur."),
		QCoreApplication::translate("main", "d")},
	{{"w", "SowWat"},
		QCoreApplication::translate("main", "SowWat."),
		QCoreApplication::translate("main", "M")},
	{{"f", "PDEN"},
		QCoreApplication::translate("main", "Number of functions."),
		QCoreApplication::translate("main", "f")},
	{{"m", "water"},
		QCoreApplication::translate("main", "Make corrections."),
		QCoreApplication::translate("main", "m")},
	{{"r", "IRGLVL"},
		QCoreApplication::translate("main", "IRGLVL."),
		QCoreApplication::translate("main", "r")},
	{{"t", "SOLDEP"},
		QCoreApplication::translate("main", "SOLDEP."),
		QCoreApplication::translate("main", "t")},
	{{"q", "DEP1"},
		QCoreApplication::translate("main", "DEP1"),
		QCoreApplication::translate("main", "q")},
	{{"z", "SALB"},
		QCoreApplication::translate("main", "SALB."),
		QCoreApplication::translate("main", "z")},
	{{"x", "CN2"},
		QCoreApplication::translate("main", "CN2."),
		QCoreApplication::translate("main", "x")},
	{{"c", "DRAINF"},
		QCoreApplication::translate("main", "DRAINF."),
		QCoreApplication::translate("main", "c")},
	{{"b", "SAT"},
		QCoreApplication::translate("main", "SAT."),
		QCoreApplication::translate("main", "b")},
	{{"k", "DUL"},
		QCoreApplication::translate("main", "DUL."),
		QCoreApplication::translate("main", "k")},
	{{"E", "EXTR"},
		QCoreApplication::translate("main", "EXTR."),
		QCoreApplication::translate("main", "E")},
    {{"j", "MAI1"},
		QCoreApplication::translate("main", "MAI1."),
		QCoreApplication::translate("main", "j")},
	{{"a", "MAI"},
		QCoreApplication::translate("main", "MAI."), 
		QCoreApplication::translate("main", "a")},
	{{"o", "INSOL"},
		QCoreApplication::translate("main", "INSOL."),
		QCoreApplication::translate("main", "o")},
	{{"y", "CROPS"},//////////////////////////////////////////////////////////////////////new
		QCoreApplication::translate("main", "CROPS."),
		QCoreApplication::translate("main", "y")},
	{{"U", "U"},
		QCoreApplication::translate("main", "U."),
		QCoreApplication::translate("main", "U")},
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
		const QString CROPSParameter = parser.value("CROPS");///////////////////////////////////////////////////new
		const int CROPS = CROPSParameter.toInt();
		if (CROPS != 0 && CROPS != 1) {////0 - outside file, 1 - inside
			std::cout << "Bad p: " + CROPS;
		}
		if (CROPS == 0)
		{
			param.func_file_name = args.at(3); // funcs
			param.h5_file_name = args.at(2); // weather
			param.h5_table_name = args.at(1); // samples
			param.crops_ini_file = args.at(0);
		}
		else
		{
			param.func_file_name = args.at(2); // funcs
			param.h5_file_name = args.at(1); // weather
			param.h5_table_name = args.at(0); // samples
		}
		param.crops = CROPS;


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
		//if (P != 0 && P != 1 && P != 3 && P != 5 && P != 7 && P != 8 && P != 2) {
		//	std::cout << "Bad p: " + P;
	//	}
		
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

		const QString LatParameter = parser.value("Latitude");
		const double Lat = LatParameter.toDouble();
		param.Latitude = Lat;

		const QString VPDFParameter = parser.value("VPDF");
		const double VPDF = VPDFParameter.toDouble();
		param.VPDF = VPDF;

		const QString FixFindParameter = parser.value("FixFind");
		const int FixFind = FixFindParameter.toInt();
		param.FixFind = FixFind;

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
	//if (param.print_trace > 0) {
	//	param.Print();
	//}
	model = new Model(param, &a);
	QObject::connect(model, SIGNAL(finished()), &a, SLOT(quit()));
	QTimer::singleShot(0, model, SLOT(run_h5()));
	return a.exec();
}
