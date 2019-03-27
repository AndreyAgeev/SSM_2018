#pragma once
//#include <QStandardItemModel>
//#include <QtSql/qsqldatabase.h>
//#include <QtSql/qsqlquery.h>

#include <highfive/H5Attribute.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <QSettings>
#include <QVariant>
#include <armadillo>
using namespace std;
using namespace HighFive;
class Data
{
public:
	struct Data_f
	{
            vector<int> years;
			vector<int> doy;
			vector<int> geo_id;
			vector<double> srad;
			vector<double> tmax;
			vector<double> tmin;
			vector<double> rain;
			vector<int> month;
			vector<double> dl;
			vector<string> clim_names = { "tmax", "tmin", "rain", "dl", "srad" };
			int nWeather;
			int nClCovar = 5;

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
		int  vpd_resp;
		double  vpd_cr;
	};
	struct Data_a
	{
		vector<int> years;
		vector<int> doy;
		vector<int> geo_id;
		arma::vec response;
		vector<vector<double>> resp;
		vector<string> species;
		vector<int> month;
		vector<vector<double>> gr_covar;
		vector<string> gr_names;
		int nSamples;
		int nGrCovar;
	};


	Data_a data_a5;

	Data_f data_h5;
	Data_phen data_p;
	//Nlreg_param param_nlreg;
	Data() {}
	void read_ini(void)
	{
		/*QSettings sett("C:\project\SSM\SSM_improved\SSM_improved\crops.ini.src", QSettings::IniFormat);
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
		data_p.CPP = sett.value("CPP", 11).toDouble();
		data_p.ppsen = sett.value("ppsen", 0.143).toDouble();
		data_p.ttSWEM = sett.value("ttSWEM", 3.0).toDouble();
		data_p.ttEMR1 = sett.value("ttEMR1", 30.3).toDouble();
		data_p.ttR1R3 = sett.value("ttR1R3", 8.3).toDouble();
		data_p.ttR3R5 = sett.value("ttR3R5 ", 7.3).toDouble();
		data_p.ttR5R7 = sett.value("ttR5R7", 36.9).toDouble();
		data_p.ttR7R8 = sett.value("ttR7R8", 6.0).toDouble();
		data_p.ttBRP = sett.value("ttBRP", 3.0).toDouble();
		data_p.ttTRP = sett.value("ttTRP", 33.3).toDouble();

		data_p.ttWSD = sett.value("ttWSD", 3.0).toDouble();
		data_p.ttR1TLM = sett.value("ttR1TLM", 8.3).toDouble();
		data_p.ttR1TLP = sett.value("ttR1TLP", 15.6).toDouble();
		data_p.ttRUE = sett.value("ttRUE", 30.3).toDouble();
		data_p.ttBSG = sett.value("ttBSG", 48.9).toDouble();
		data_p.ttTSG = sett.value("ttTSG", 85.8).toDouble();
		data_p.ttBRG = sett.value("ttBRG", 3.0).toDouble();
		data_p.ttTRG = sett.value("ttTRG", 48.9).toDouble();
		data_p.ttBNF = sett.value("ttBNF", 14.7).toDouble();
		data_p.TRESH = sett.value("TRESH", 0.75).toDouble();
		data_p.ttDKill = sett.value("ttDKill", 48.9).toDouble();

		data_p.LtFtsw = sett.value("LtFtsw", 0).toDouble();
		data_p.LtWdDur = sett.value("LtWdDur ", 2).toDouble();
		data_p.vpd_resp = sett.value("vpd_resp ", 1).toDouble();
		data_p.vpd_cr = sett.value("vpd_cr", 20.0).toDouble();
		sett.endGroup();*/
		QSettings sett("crops.ini.src", QSettings::IniFormat);
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
		data_p.ppsen = sett.value("ppsen", 0.11).toDouble();
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
			File file(file_name.toStdString(), File::ReadOnly);
			DataSet doy_read = file.getDataSet("doy");
			DataSpace space = doy_read.getSpace();
			data_h5.nWeather = space.getDimensions()[0];
			doy_read.read(data_h5.doy);
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
			DataSet dl_read = file.getDataSet("dl");
			dl_read.read(data_h5.dl);
			DataSet month_read = file.getDataSet("month");
			month_read.read(data_h5.month);
			DataSet geo_id_read = file.getDataSet("geo_id");
			geo_id_read.read(data_h5.geo_id);
		}
		catch (Exception &err) {
			std::cerr << err.what() << std::endl;
		}

	}
	void read_spieces(QString file_name, bool extra_covar)
	{
		try {
			File file(file_name.toStdString(), File::ReadOnly);

			DataSet doy_read = file.getDataSet("doy");
			DataSpace space = doy_read.getSpace();
			data_a5.nSamples = space.getDimensions()[0];
			doy_read.read(data_a5.doy);
			DataSet year_read = file.getDataSet("year");
			year_read.read(data_a5.years);
			DataSet month_read = file.getDataSet("month");
			month_read.read(data_a5.month);
			DataSet geo_id_read = file.getDataSet("geo_id");
			geo_id_read.read(data_a5.geo_id);
			DataSet b_read = file.getDataSet("response");
			space = b_read.getSpace();
			int ns = space.getDimensions()[0]; // number of samples
			int nn = space.getDimensions()[1]; // number of measurements
			b_read.read(data_a5.resp);
			assert(ns == data_a5.nSamples);
			assert(nn == 1);
			DataSet a1_read = file.getDataSet("species");
			a1_read.read(data_a5.species);
			data_a5.nGrCovar = 0;
			if (extra_covar) {
				DataSet c0_read = file.getDataSet("gr_covar");
				DataSpace space = c0_read.getSpace();
				data_a5.nGrCovar = space.getDimensions()[1]; // number of measurements
				c0_read.read(data_a5.gr_covar);
				DataSet c1_read = file.getDataSet("gr_names");
				c1_read.read(data_a5.gr_names);
			}
		}
		catch (Exception& err) {
			std::cerr << err.what() << std::endl;
		}
		data_a5.response = Data::std2arvec(data_a5.resp, data_a5.nSamples, 0);
	}
private:
	arma::mat std2arvec(std::vector<std::vector<double> > &vec, int n_rows, int offset) {
		arma::vec Y(n_rows, 1);
		for (size_t i = 0; i < n_rows; ++i) {
			Y(i) = vec[i][offset];
		}
		return Y;
	}


};

