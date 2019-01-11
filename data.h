#pragma once
#include <QtSql/qsqldatabase.h>
#include <QtSql/qsqlquery.h>

#include <highfive/H5Attribute.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <QSettings>
#include <QVariant>
using namespace std;
using namespace HighFive;
class Data
{
public:

	struct Data_f
	{
            vector<int> years;
			vector<int> doy;
			vector<double> srad;
			vector<double> tmax;
			vector<double> tmin;
			vector<double> rain;
	};
	struct Data_phen
	{
		double phyl;
		double PLACON;
		double PLAPOW30;
		double SLA;
		/////
		double TBRUE;
		double TP1RUE;
		double TP2RUE;
		double TCRUE;
		double KPAR;
		double IRUE1;
		double IRUE2;
		/////
		double FLF1A;
		double FLF1B;
		double WTOPL;
		double FLF2;
		//////
		double FRTRL;
		double GCF;
		double PDHI;
		double WDHI1;
		double WDHI2;
		double WDHI3;
		double WDHI4;
		//////
		double DEPORT;
		double EED;
		double GRTDP;
		double TEC;
		double WSSG;
		double WSSL;
		double WSSD;
		////////
		double SLNG;
		double SLNS;
		double SNCG;
		double SNCS;
		double GNC;
		double MXNUP;
		double WSSN;
		////////////////////////////
		double TBD;
		double TP1D;
		double TP2D;
		double TCD;

		double CPP;
		double ppsen;

		double ttSWEM;
		double  ttEMR1;
		double  ttR1R3;
		double  ttR3R5;
		double  ttR5R7;
		double  ttR7R8;

		double  ttBRP;
		double  ttTRP;

		double ttWSD;
		double  ttR1TLM;
		double  ttR1TLP;
		double  ttRUE;
		double  ttBSG;
		double  ttTSG;
		double  ttBRG;
		double  ttTRG;
		double  ttBNF;
		double  TRESH;
		double  ttDKill;
		double  LtFtsw;
		double  LtWdDur;
		double  vpd_resp;
		double  vpd_cr;
	};
	Data_f data_h5;
	Data_phen data_p;
	Data() {}
	void read_ini(void)
	{
		QSettings sett("crops.ini", QSettings::IniFormat);
		sett.beginGroup("Jam");
		data_p.phyl = sett.value("phyl", 46).toDouble();
		data_p.PLACON = sett.value("PLACON", 1).toDouble();
		data_p.PLAPOW30 = sett.value("PLAPOW30", 2.158).toDouble();
		data_p.SLA = sett.value("SLA", 0.021).toDouble();
		data_p.TBRUE = sett.value("TBRUE", 2).toDouble();
		data_p.TP1RUE = sett.value("TP1RUE", 14).toDouble();
		data_p.TP2RUE = sett.value("TP2RUE", 30).toDouble();
		data_p.TCRUE = sett.value("TCRUE", 38).toDouble();
		data_p.KPAR = sett.value("KPAR", 0.65).toDouble();
		data_p.IRUE1 = sett.value("IRUE1", 1.8).toDouble();
		data_p.IRUE2 = sett.value("IRUE2", 1.8).toDouble();
		data_p.FLF1A = sett.value("FLF1A", 0.53).toDouble();
		data_p.FLF1B = sett.value("FLF1B", 0.3).toDouble();
		data_p.WTOPL = sett.value("WTOPL", 180).toDouble();
		data_p.FLF2 = sett.value("FLF2", 0.13).toDouble();
		data_p.FRTRL = sett.value("FRTRL", 0.22).toDouble();
		data_p.GCF = sett.value("GCF", 1).toDouble();
		data_p.PDHI = sett.value("PDHI", 0.02).toDouble();
		data_p.WDHI1 = sett.value("WDHI1", 0).toDouble();
		data_p.WDHI2 = sett.value("WDHI2", 0).toDouble();
		data_p.WDHI3 = sett.value("WDHI3", 450).toDouble();
		data_p.WDHI4 = sett.value("WDHI4", 2000).toDouble();
		data_p.DEPORT = sett.value("DEPORT", 200).toDouble();
		data_p.EED = sett.value("EED", 1000).toDouble();
		data_p.GRTDP = sett.value("GRTDP", 17).toDouble();
		data_p.TEC = sett.value("TEC", 5).toDouble();
		data_p.WSSG = sett.value("WSSG", 0.3).toDouble();
		data_p.WSSL = sett.value("WSSL", 0.4).toDouble();

		data_p.WSSD = sett.value("WSSD", 0.4).toDouble();
		data_p.SLNG = sett.value("SLNG", 2.3).toDouble();
		data_p.SLNS = sett.value("SLNS", 0.78).toDouble();
		data_p.SNCG = sett.value("SNCG", 0.025).toDouble();
		data_p.SNCS = sett.value("SNCS", 0.0078).toDouble();
		data_p.GNC = sett.value("GNC", 0.043).toDouble();
		data_p.MXNUP = sett.value("MXNUP", 0.45).toDouble();
		data_p.WSSN = sett.value("WSSN", 0.5).toDouble();
		data_p.TBD = sett.value("TBD", 2).toDouble();
		data_p.TP1D = sett.value("TP1D", 21).toDouble();
		data_p.TP2D = sett.value("TP2D", 30).toDouble();

		data_p.TCD = sett.value("TCD", 40).toDouble();
		data_p.CPP = sett.value("CPP", 18).toDouble();
		data_p.ppsen = sett.value("ppsen", 0.1).toDouble();
		data_p.ttSWEM = sett.value("ttSWEM", 5.0).toDouble();
		data_p.ttEMR1 = sett.value("ttEMR1", 23.0).toDouble();
		data_p.ttR1R3 = sett.value("ttR1R3", 9.0).toDouble();
		data_p.ttR3R5 = sett.value("ttR3R5 ", 3.0).toDouble();
		data_p.ttR5R7 = sett.value("ttR5R7", 18.0).toDouble();
		data_p.ttR7R8 = sett.value("ttR7R8", 11.0).toDouble();
		data_p.ttBRP = sett.value("ttBRP", 5.0).toDouble();
		data_p.ttTRP = sett.value("ttTRP", 28.0).toDouble();

		data_p.ttWSD = sett.value("ttWSD", 5.0).toDouble();
		data_p.ttR1TLM = sett.value("ttR1TLM", 9.0).toDouble();
		data_p.ttR1TLP = sett.value("ttR1TLP", 12.0).toDouble();
		data_p.ttRUE = sett.value("ttRUE", 23.0).toDouble();
		data_p.ttBSG = sett.value("ttBSG", 40.0).toDouble();
		data_p.ttTSG = sett.value("ttTSG", 64.0).toDouble();
		data_p.ttBRG = sett.value("ttBRG", 5.0).toDouble();
		data_p.ttTRG = sett.value("ttTRG", 40.0).toDouble();
		data_p.ttBNF = sett.value("ttBNF", 14.7).toDouble();
		data_p.TRESH = sett.value("TRESH", 0.75).toDouble();
		data_p.ttDKill = sett.value("ttDKill", 40.0).toDouble();

		data_p.LtFtsw = sett.value("LtFtsw", 0).toDouble();
		data_p.LtWdDur = sett.value("LtWdDur ", 2).toDouble();
		data_p.vpd_resp = sett.value("vpd_resp ", 1).toDouble();
		data_p.vpd_cr = sett.value("vpd_cr", 20.0).toDouble();
		sett.endGroup();
	}
	void read_h5(QString file_name)
	{
		try {
			cout << "BEGIN READ" << endl;
            File file(file_name.toStdString(), File::ReadOnly);
			DataSet doy_read = file.getDataSet("doy");
			doy_read.read(data_h5.doy);
			cout << "END READ" << endl;
			DataSet rain_read = file.getDataSet("rain");
			rain_read.read(data_h5.rain);
			DataSet srad_read = file.getDataSet("srad");
			srad_read.read(data_h5.srad);
			DataSet tmax_read = file.getDataSet("tmax");
			tmax_read.read(data_h5.tmax);
			DataSet tmin_read = file.getDataSet("tmin");
			tmin_read.read(data_h5.tmin);
			DataSet year_read = file.getDataSet("year");
			year_read.read(data_h5.years);
		}
		catch (Exception &err){
			std::cerr << err.what() << std::endl;
		}

	}
};
