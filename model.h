#ifndef _MODEL_H_
#define _MODEL_H_
#include <string>
#include <iostream>
#include <QtCore>
#include <QDebug>
#include <armadillo>


#include <highfive/H5Attribute.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include "data.h"
#include "parametrs.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <nlreg.h>
using namespace std;
using namespace HighFive;
#define SOLTANI_FUNC 0
class Model : public QObject
{
	Q_OBJECT
public:
	int ROW = 0;
	int MAT;
	int index_lai = -1;
	int iniPheno;//INT
	int iniLai;
	int iniDMP;
	int iniDMD;
	int iniSW;
	int iniPNB;
	int CumFind;

	//SoilWater
	double CLL;
	double ISATSW;
	double ATSW;
	double TTSW;
	double FTSW;
	double WSTORG;
	double ATSW1;
	double TTSW1;
	double FTSW1;
	double WLL1;
	double WAT1;
	double WSAT1;
	double WLL;
	double WATRT;
	double WSAT;
	double EOSMIN;
	double WETWAT;
	double KET;
	double CALB;
	double CTR;
	double CE;
	double CRAIN;
	double CRUNOF;
	double CIRGW;
	double IRGNO;
	double DDMP;

	double LAI;
	double LtDrCntr;
	double SE2C;
	double SE1MX;
	double DSR;
	double SSE1;
	double SSE;

	double IRGW;

	double DRAIN1;
	double DRAIN;
	double GRTD;

	double CBD;

	double EWAT;
	double RUNOF;
	double s;
	double SWER;
	double ETLAI;
	double BSGLAI;
	double TD;
	double ALBEDO;
	double EEQ;
	double PET;

	double EOS;
	double SEVP;
	double DYSE;

	int semethod = 1;
	double FLUX1;
	int vpdtp;
	double VPTMIN;
	double VPTMAX;
	double VPD;
	double TR;


	double TR1;
	double RT1;
	double WSFN;
	double WSFG;
	double WSFL;
	double WSFD;

	double LtDrCnt;

	double TMP;

	double DOY;
	int dtEM;
	int dtR1;
	int dtR3;
	int dtR5;
	int dtR7;
	int  dtR8;
	double  MXLAI;
	double BSGDM;
	double  WTOP;
	double WGRN;
	double  NLF;
	double  NST;
	double  NVEG;
	double  NGRN;
	double  CNUP;

	double bdEM;

	double bdR1;
	double bdR3;

	double bdR5;
	double bdR7;
	double bdR8;

	double  bdBLG;

	double  bdTLM;
	double bdTLP;

	double cbdEM;

	double cbdR1;
	double cbdR3;

	double cbdR5;
	double cbdR7;
	double cbdR8;


	int DAP;
	double tempfun;
	double DTT;

	double SABH;
	double Pi;
	double RDN;
	double ALPHA;
	double SMA3;
	double LANDA;
	double DEC;
	double TALSOC;

	double CEDSOC;
	double SOCRA;
	vector <double> DL;
	double pp;
	double ppfun;
	double bd;


	double PLAPOW;
	double MSNN;
	double PLA2;
	double PLA1;

	double GLAI;

	double DLAI;

	double INODE;
	double GLF;

	double TCFRUE;
	double RUE;
	double VPDcr;
	double LAT;
	double DECL;
	double SINLD;
	double COSLD;
	double AOB;
	double AOB2;
	double DAYL;
	double DSINB;
	double DSINBE;
	double SC;
	double DSO;
	double Bnoon;
	double DTR;
	double TC;
	double P;
	double SUNRIS;
	double SUNSET;
	std::ofstream out;
	std::ofstream out_s;
	std::ofstream out_LAI;
	std::ofstream out_DAP;
	std::ofstream out_FTSW;
	std::ofstream out_CE;
	std::ofstream out_CTR;
	std::ofstream out_MSNN;
	std::ofstream out_CNUP;
	std::ofstream out_NVEG;
	std::ofstream out_NGRN;
	std::ofstream out_error;
	std::ofstream out_cbd;
	std::ofstream out_dtR1;
	std::ofstream out_bd;


	double TSUNST;
	double NIGHTL;
	double TEMP1;
	double TMINA;

	double TMAXB;
	double SINB;
	double BET;
	double SRAD1;
	double BETA;
	double FINT;
	double DDMP1;

	double VPTEMP1;
	double VPD1;



	double NUP;
	double XNLF;
	double XNST;
	double INLF;
	double INST;
	double INGRN;

	double NFC;
	double PDNF;
	double DNF;
	double TRLN;
	double FXLF;
	double NUP2;

	double WLF;
	double WST;
	double WVEG;
	double TRANSL;
	double SGR;
	double DHI;
	double DHIF;
	double TRLDM;
	double DDMP2;
	double FLF1;
	double GST;
	double HI;



	Nlreg *nl;
	Data data;
	Parametrs param;


	bool write_check = false;


	int NSAM = 0;
	bool check_last_cbd = false;

	Model(Parametrs new_param, QObject *parent = 0) : QObject(parent), param(new_param)
	{
		if (param.print_trace > 1) cout << "begin read" << endl;
		data.read_h5(param.h5_file_name);
		data.read_spieces(param.h5_table_name, param.ecovar);
		data.read_ini(param.crops_ini_file);
	//	if (param.crops == 0)
		nl = new Nlreg(param.func_file_name, data.data_h5.clim_names, data.data_a5.gr_names, param.nF, param.wL, param.rT, param.print_trace, param.crops);
//		else
	//		nl = new Nlreg(param.crops_ini_file, data.data_h5.clim_names, data.data_a5.gr_names, param.nF, param.wL, param.rT, param.print_trace, param.crops);
		nl->nlreg_build();
		if (param.print_trace > 1) cout << "end read" << endl;
	}
	void SoilWater()
	{
		if (iniSW == 0)
		{
			CLL = param.DUL - param.EXTR;
			//Parameters and Initials
			ISATSW = param.SOLDEP * param.EXTR * param.MAI;
			ATSW = data.data_p.DEPORT * param.EXTR * param.MAI1;
			TTSW = data.data_p.DEPORT * param.EXTR;
			FTSW = ATSW / TTSW;
			WSTORG = ISATSW - ATSW;

			ATSW1 = param.DEP1 *  param.EXTR * param.MAI1;
			TTSW1 = param.DEP1 *  param.EXTR;
			FTSW1 = ATSW1 / TTSW1;

			WLL1 = param.DEP1 * CLL;
			WAT1 = WLL1 + ATSW1;
			WSAT1 = param.DEP1  * param.SAT;

			WLL = data.data_p.DEPORT * CLL;
			WATRT = WLL + ATSW;
			WSAT = data.data_p.DEPORT *  param.SAT;

			EOSMIN = 1.5;
			WETWAT = 10;
			KET = 0.5;
			CALB = 0.23;

			DYSE = 1.0;
			CTR = 0.0;
			CE = 0.0;
			CRAIN = 0.0;
			CRUNOF = 0.0;
			CIRGW = 0.0;
			IRGNO = 0.0;
			DDMP = 0.0;
			LAI = 0.0;
			LtDrCntr = 0.0;

			SE2C = 3.5;
			SE1MX = param.U;
			DSR = 1.0;
			SSE1 = param.U;
			SSE = param.U + SE2C;
			iniSW = 1;
		}
		//Irrigation
		if (param.water == 1 && FTSW <= param.IRGLVL && CBD < data.data_p.ttTSG)
		{
			IRGW = (TTSW - ATSW);
			IRGNO = IRGNO + 1.0;
		}
		else
			IRGW = 0.0;

		CIRGW = CIRGW + IRGW;

		// Drainage
		if (ATSW1 <= TTSW1)
			DRAIN1 = 0.0;
		else if (ATSW1 > TTSW1)
			DRAIN1 = (ATSW1 - TTSW1) * param.DRAINF;

		if (ATSW <= TTSW)
			DRAIN = 0.0;
		else if (ATSW > TTSW)
			DRAIN = (ATSW - TTSW) * param.DRAINF;
		WSTORG = WSTORG + DRAIN - EWAT;
		if (WSTORG < 0.0)
			WSTORG = 0.0;

		//Water exploitation by root growth
		GRTD = data.data_p.GRTDP;//
		if (CBD < data.data_p.ttBRG)
			GRTD = 0.0;
		if (CBD > data.data_p.ttTRG)
			GRTD = 0.0;
		if (DDMP == 0.0)
			GRTD = 0.0;
		if (data.data_p.DEPORT >= param.SOLDEP)
			GRTD = 0.0;
		if (data.data_p.DEPORT >= data.data_p.EED)
			GRTD = 0.0;
		if (WSTORG == 0.0)
			GRTD = 0.0;
		data.data_p.DEPORT = data.data_p.DEPORT + GRTD;

		EWAT = GRTD * param.EXTR;
		if (EWAT > WSTORG)
			EWAT = WSTORG;


		RUNOF = 0;
		if (param.water == 2 && data.data_h5.rain[ROW] > 0.01)
		{
			s = 254.0 * (100.0 / param.CN2 - 1.0);
			SWER = 0.15 * ((WSAT1 - WAT1) / (WSAT1 - WLL1));
			if (SWER < 0.0)
				SWER = 0.0;
			if ((data.data_h5.rain[ROW] - SWER * s) > 0.0)
				RUNOF = ((data.data_h5.rain[ROW] - SWER * s) * (data.data_h5.rain[ROW] - SWER * s)) / (data.data_h5.rain[ROW] + (1.0 - SWER) * s);
			else
				RUNOF = 0.0;
		}

		if ((WAT1 - DRAIN1) > WSAT1)
			RUNOF = RUNOF + (WAT1 - DRAIN1 - WSAT1);

		CRAIN = CRAIN + data.data_h5.rain[ROW];
		CRUNOF = CRUNOF + RUNOF;

		//LAI for soil evaporation
		if (CBD <= data.data_p.ttBSG)//bdBSG)
			ETLAI = LAI;
		else
			ETLAI = BSGLAI;

		// Potential ET
		TD = 0.6 * data.data_h5.tmax[ROW] + 0.4 * data.data_h5.tmin[ROW];
		ALBEDO = CALB * (1.0 - exp(-KET * ETLAI)) + param.SALB * exp((-KET) * ETLAI);
		EEQ = data.data_h5.srad[ROW] * (0.004876 - 0.004374 * ALBEDO) * (TD + 29.0);
		PET = EEQ * 1.1;
		if (data.data_h5.tmax[ROW] > 34.0)
			PET = EEQ * ((data.data_h5.tmax[ROW] - 34.0) * 0.05 + 1.1);
		if (data.data_h5.tmax[ROW] < 5.0)
			PET = EEQ * 0.01 * exp(0.18 * (data.data_h5.tmax[ROW] + 20.0));

		// Soil evaporation
		EOS = PET * exp(((-KET) * ETLAI));
		if (PET > EOSMIN && EOS < EOSMIN)
			EOS = EOSMIN;

		SEVP = EOS;
		if ((data.data_h5.rain[ROW] + IRGW) > WETWAT)
			DYSE = 1.0;
		if (DYSE > 1.0 || FTSW < 0.5 || ATSW1 <= 2.0)
		{
			SEVP = EOS * (pow((DYSE + 1.0), 0.5) - pow(DYSE, 0.5));
			DYSE = DYSE + 1.0;
		}
		if (semethod == 2)
		{
			if (ATSW1 < 0.0)
				SEVP = 0.0;
			else
			{

				if (SSE1 < SE1MX)
				{
					//	'Stage I evaporation
					SEVP = EOS;
					if (SEVP > (SE1MX - SSE1))
						SEVP = SE1MX - SSE1;
					SSE1 = SSE1 + SEVP;
					SSE = SSE + SEVP;
					if (SSE1 >= SE1MX)
					{
						//Transition from Stage I to Stage II
						SEVP = SEVP + SE2C * (pow(DSR, 0.5) - pow((DSR - 1.0), 0.5)) * (1.0 - (SEVP / EOS));
						DSR = DSR + 1.0 - SEVP / EOS;
					}
				}
				else
				{
					//Stage II evaporation
					SEVP = SE2C * (pow(DSR, 0.5) - pow((DSR - 1.0), 0.5));
					if (SEVP > EOS)
						SEVP = EOS;
					DSR = DSR + 1.0;
					SSE = SSE + SEVP;
				}

				FLUX1 = data.data_h5.rain[ROW] + IRGW - RUNOF;
				if (FLUX1 >= SSE1)
				{
					SSE = SSE - FLUX1;
					if (SSE < 0)
					{
						SSE = 0.0;
						SSE1 = 0.0;
						DSR = 1.0 + pow((SSE / SE2C), 0.5);
					}
				}
				else
				{
					SSE = SSE - FLUX1;
					SSE1 = SSE1 - FLUX1;
					DSR = 1.0 + pow((SSE / SE2C), 0.5);
				}
			}
		}
		CE = CE + SEVP;
		vpdtp = data.data_p.vpd_resp;
		VPDcr = data.data_p.vpd_cr;
		//Plant transpiration
		if (vpdtp == 1 || vpdtp == 3)
		{
			VPTMIN = 0.6108 * exp(17.27 * data.data_h5.tmin[ROW] / (data.data_h5.tmin[ROW] + 237.3));
			VPTMAX = 0.6108 * exp(17.27 * data.data_h5.tmax[ROW] / (data.data_h5.tmax[ROW] + 237.3));
			VPD = param.VPDF * (VPTMAX - VPTMIN);
			TR = DDMP * VPD / data.data_p.TEC;    //    'VPD in kPa, TEC in Pa
		}
		else if (vpdtp == 2)
		{
		}

		if (TR < 0.0)
			TR = 0.0;
		CTR = CTR + TR;

		if (data.data_p.DEPORT <= param.DEP1)
			TR1 = TR;
		else if (data.data_p.DEPORT > param.DEP1)
		{
			if (FTSW1 > data.data_p.WSSG)
				RT1 = 1.0;
			else
				RT1 = FTSW1 / data.data_p.WSSG;
			TR1 = TR * RT1;
		}
		//  Updating
		ATSW1 = ATSW1 + data.data_h5.rain[ROW] + IRGW - DRAIN1 - RUNOF - TR1 - SEVP;
		if (ATSW1 < 0.0)
			ATSW1 = 0.0;
		FTSW1 = ATSW1 / TTSW1;
		WAT1 = WLL1 + ATSW1;

		ATSW = ATSW + data.data_h5.rain[ROW] + IRGW + EWAT - DRAIN - RUNOF - TR - SEVP;
		if (ATSW < 0.0)
			ATSW = 0.0;
		TTSW = data.data_p.DEPORT  * param.EXTR;
		FTSW = ATSW / TTSW;

		WLL = data.data_p.DEPORT  * CLL;
		WATRT = WLL + ATSW;
		WSAT = data.data_p.DEPORT  * param.SAT;

		// Water-stress-factors
		if (FTSW > data.data_p.WSSN)
			WSFN = 1.0;
		else
			WSFN = FTSW / data.data_p.WSSN;
		if (FTSW > data.data_p.WSSL)
			WSFL = 1.0;
		else
			WSFL = FTSW / data.data_p.WSSL;
		if (FTSW > data.data_p.WSSG)
			WSFG = 1.0;
		else
			WSFG = FTSW / data.data_p.WSSG;
		WSFD = (1.0 - WSFG) * data.data_p.WSSD + 1.0;

		if (WATRT > (0.95 * WSAT))
		{
			WSFN = 0.0;
			WSFG = 0.0;
			WSFL = 0.0;
			WSFD = 0.0;
		}

		// Crop termination by water stress
		if (CBD > data.data_p.ttDKill/*bdDKill*/ && (CBD < data.data_p.ttTSG/*bdTSG*/))
		{
			//If LAI < 0.05 Then CBD = bdTSG
			//If FTSW < 0.02 And VPD > 2.2 Then CBD = bdTSG
			//If FTSW <= 0 And VPD > 1.8 Then CBD = bdTSG
			if (FTSW <= data.data_p.LtFtsw)
			{
				LtDrCntr = LtDrCntr + 1.0;
			}
			else
				LtDrCnt = 0.0;

		}
		if (LtDrCntr >= data.data_p.LtWdDur && CBD < data.data_p.ttTSG)
		{
			cout << "change cbd not opt" << endl;
			CBD = data.data_p.ttTSG;
		}
	}

	void FindSowingData(int geo_id, int start_day, int start_year)
	{
		/*if (param.FixFind == 2)
		{
			CumFind = 0.0;
			for (int row = ROW; row < data.data_h5.years.size(); row++)
			{
				ROW = row;
				DOY = data.data_h5.doy[ROW];;
				CBD = 0.0;
				SoilWater();
				CumFind = CumFind + 1;
				if (CumFind > param.SearchDur)
					MAT = 1;
				if (MAT == 1)
				{
					DOY = -1;
					dtEM = 0.0;
					dtR1 = 0.0;
					dtR3 = 0.0;
					dtR5 = 0.0;
					dtR7 = 0.0;
					dtR8 = 0.0;
					MXLAI = 0.0;
					BSGLAI = 0.0;
					BSGDM = 0.0;
					WTOP = 0.0;
					WGRN = 0.0;
					NLF = 0.0;
					NST = 0.0;
					NVEG = 0.0;
					NGRN = 0.0;
					CNUP = 0.0;
					param.INSOL = 0.0;
				}
				if (ATSW1 + WSTORG >= param.SowWat || MAT == 1)
					break;
			}
			Pdoy = data.data_h5.doy[ROW];
		}*/
	}
	void Weather(void)
	{
		//	ROW += 1;
		TMP = (data.data_h5.tmax[ROW] + data.data_h5.tmin[ROW]) / 2.0;
	}

	void Phenology()
	{
		int threshold_index;
		if (iniPheno == 0)
		{
			bdEM = data.data_p.ttSWEM;
			bdR1 = bdEM + data.data_p.ttEMR1;//Thermal time from sowing to emergence
			bdR3 = bdR1 + data.data_p.ttR1R3;
			bdR5 = bdR3 + data.data_p.ttR3R5;
			bdR7 = bdR5 + data.data_p.ttR5R7;
			bdR8 = bdR7 + data.data_p.ttR7R8;
			bdBLG = bdEM;
			bdTLM = bdR1 + data.data_p.ttR1TLM;
			bdTLP = bdR1 + data.data_p.ttR1TLP;
			DAP = 0;
			CBD = 0.0;
			WSFD = 1.0;
			iniPheno = 1;
		}
		switch (param.threshold)
		{
		case 0:
			bdEM = nl->get_cbd();
			break;
		case 1:
			bdR1 = nl->get_cbd();
			break;
		case 3:
			bdR3 = nl->get_cbd();
			break;
		case 5:
			bdR5 = nl->get_cbd();
			break;
		case 7:
			bdR7 = nl->get_cbd();
			break;
		case 8:
			bdR8 = nl->get_cbd();
			break;
		}
		if (param.function_mode == SOLTANI_FUNC)
		{
			// Thermal time calculation
			if (TMP <= data.data_p.TBD || TMP >= data.data_p.TCD)
				tempfun = 0.0;
			else if (TMP > data.data_p.TBD && TMP < data.data_p.TP1D)
				tempfun = (TMP - data.data_p.TBD) / (data.data_p.TP1D - data.data_p.TBD);
			else if (TMP > data.data_p.TP2D && TMP < data.data_p.TCD)
				tempfun = (data.data_p.TCD - TMP) / (data.data_p.TCD - data.data_p.TP2D);
			else if (TMP >= data.data_p.TP1D && TMP <= data.data_p.TP2D)
				tempfun = 1.0;

			DTT = (data.data_p.TP1D - data.data_p.TBD) * tempfun;
			if (CBD > data.data_p.ttWSD)
				tempfun = tempfun * WSFD;
			if (CBD > data.data_p.ttWSD)
				DTT = DTT * WSFD;
			// Photoperiod function
			SABH = 6.0;
			Pi = 3.141592654;
			RDN = Pi / 180.0;
			ALPHA = 90.0 + SABH;
			SMA3 = 0.9856 * DOY - 3.251;
			LANDA = SMA3 + 1.916 * sin(SMA3 * RDN) + 0.02 * sin(2.0 * SMA3 * RDN) + 282.565;
			DEC = 0.39779 * sin(LANDA * RDN);
			DEC = atan(DEC / sqrt(1.0 - (DEC * DEC)));
			DEC = DEC / RDN;
			TALSOC = 1.0 / cos(LAI * RDN);
			CEDSOC = 1.0 / cos(DEC * RDN);
			SOCRA = (cos(ALPHA * RDN) * TALSOC * CEDSOC) - (tan(LAI * RDN) * tan(DEC * RDN));
			//	DL = Pi / 2.0 - (atan(SOCRA / (((1.0 - (SOCRA *  SOCRA))) * (1.0 - (SOCRA *  SOCRA)))));
			//	DL = Pi / 2.0 - (atan(SOCRA / sqrt(1.0 - (SOCRA * SOCRA))));
			//	DL = DL / RDN;
			//	pp = 2.0/ 15.0 * DL;
			pp = 2.0 / 15.0 * data.data_h5.dl[ROW];

			if (data.data_p.ppsen >= 0.0)
				ppfun = 1.0 - data.data_p.ppsen * (data.data_p.CPP - pp);
			else if (data.data_p.ppsen < 0.0)
				ppfun = 1.0 - (-data.data_p.ppsen) * (pp - data.data_p.CPP);
			if (ppfun > 1.0)
				ppfun = 1.0;
			if (ppfun < 0.0)
				ppfun = 0.0;

			if (CBD < data.data_p.ttBRP)
				ppfun = 1.0;
			if (CBD > data.data_p.ttTRP)
				ppfun = 1.0;

			bd = tempfun * ppfun;
			CBD = CBD + bd;
		}
		else
		{

			vector<double> clim_covar = { data.data_h5.tmax[ROW], data.data_h5.tmin[ROW], data.data_h5.rain[ROW], data.data_h5.dl[ROW], data.data_h5.srad[ROW] };
			//		if (param.ecovar == 1)// 0 - без снипов, 1 со снипами
				//	cout << "gr_covar[nsam]= " << data.data_a5.gr_covar[NSAM].size() << endl;
	//		cout << "    TEMP DATA: " << endl;		
	//		cout << "_________________________ " << endl;
	//		cout << "TMAX  TMIN  RAIN  DL  SRAD " << endl;

	//		cout << data.data_h5.tmax[ROW] << " " << data.data_h5.tmin[ROW] << " " <<  data.data_h5.rain[ROW] << " " <<  data.data_h5.dl[ROW] << " " <<  data.data_h5.srad[ROW] << endl;

	//		cout << "_________________________ " << endl;
		//	   bd = nl->get_func_value(clim_covar);
		if(param.snp_mode == 1)
	       bd = nl->get_func_value(clim_covar, data.data_a5.gr_covar[NSAM]);
		else
			bd = nl->get_func_value(clim_covar);
		//	out_bd.open(param.h5_table_name.toStdString() + "_" + "bd.txt", std::ios::app);
		//	out_bd << "bd = " << bd << endl;
	//		cout << "THIS bd = " << bd;
		//		else
		//		{
		//			
		//		}
			CBD += Heaviside(bd) * bd;
		//	out_bd << "CBD = " << CBD << endl;

	//		cout << " THIS CBD = " << CBD;
	//		out_bd.close();

		}

		DAP = DAP + 1;
		//	cout << "THIS DAP = " << DAP << endl;
		if (CBD < bdEM)
		{
			dtEM = DAP + 1.0;  // 'Saving days to EMR
			cbdEM = CBD;
		}

		if (CBD < bdR1)
		{
		//	cout << "THIS bdR1 = " << bdR1 << endl;

			dtR1 = DAP + 1; // 'Saving days to R1
			cbdR1 = CBD;
		//	cout << "THIS cbdR1 = " << cbdR1 << endl;
		}
	//	if (CBD >= bdR1 && check_last_cbd == false)
	//	{
	//		cbdR1 = CBD;
	//		check_last_cbd = true;
	//	}

		if (CBD < bdR3)
		{
			dtR3 = DAP + 1;  // 'Saving days to R3
			cbdR3 = CBD;
		}
		if (CBD < bdR5)
		{
			dtR5 = DAP + 1;//  'Saving days to R5
			cbdR5 = CBD;
		}
		if (CBD < bdR7)
		{
			dtR7 = DAP + 1; //  'Saving days to R7
			cbdR7 = CBD;
		}
		if (CBD < bdR8)
		{
			dtR8 = DAP + 1; // 'Saving days to R8
			cbdR8 = CBD;
		}
	//	cout << "BDR8" << bdR8 << endl;
		// Maturity ?
		if (CBD > bdR8)
			MAT = 1;
	}
	int get_curr_day(void)
	{
		if (param.threshold == 0)
			return dtEM;
		if (param.threshold == 1)
			return dtR1;
		if (param.threshold == 3)
			return dtR3;
		if (param.threshold == 5)
			return dtR5;
		if (param.threshold == 7)
			return dtR7;
		if (param.threshold == 8)
			return dtR8;
		return (double)(-param.nD);
	}
/*	double get_phase_change(void)
	{
		if (param.threshold == 0)
			return bdEM;
		if (param.threshold == 1)
			return bdR1;
		if (param.threshold == 3)
			return bdR3;
		if (param.threshold == 5)
			return bdR5;
		if (param.threshold == 7)
			return bdR7;
		if (param.threshold == 8)
			return bdR8;
	}*/
	double get_cbd(void)
	{
		if (param.threshold == 0)
			return cbdEM;
		if (param.threshold == 1)
			return cbdR1;
		if (param.threshold == 3)
			return cbdR3;
		if (param.threshold == 5)
			return cbdR5;
		if (param.threshold == 7)
			return cbdR7;
		if (param.threshold == 8)
			return cbdR8;
	return CBD; 
	}
	void CropLAIN(void)
	{
		if (iniLai == 0)
		{
			PLAPOW = data.data_p.PLAPOW30 * (-0.002 * param.PDEN + 1.0612);//       'Valid for (10<pden<70)

			MSNN = 1.0;
			PLA2 = 0.0;
			PLA1 = 0.0;
			LAI = 0.0;
			MXLAI = 0.0;
			WSFL = 1.0;
			data.data_p.SLNG = 2.0;
			iniLai = 1;
		}
		// Yesterday LAI to intercept PAR today
		if (GLAI > (INLF / data.data_p.SLNG))
			GLAI = (INLF / data.data_p.SLNG);
		LAI = LAI + GLAI - DLAI;
		if (LAI < 0.0)
			LAI = 0.0;
		if (LAI > MXLAI)
			MXLAI = LAI;  //'Saving maximum LAI

		// Daily increase and decrease in LAI
		if (CBD <= bdBLG)
			GLAI = 0.0;
		else if (CBD > bdBLG && CBD <= bdTLM)
		{
			INODE = DTT / data.data_p.phyl;
			MSNN = MSNN + INODE * WSFL;
			PLA2 = data.data_p.PLACON * pow(MSNN, PLAPOW);
			GLAI = (((PLA2 - PLA1) * (double)param.PDEN) / 10000.0);
			PLA1 = PLA2;
		}
		else if (CBD > bdTLM && CBD <= bdTLP)
		{
			GLAI = GLF * data.data_p.SLA;
			BSGLAI = LAI;
		}
		//Saving LAI at BSG
		else if (CBD > bdTLP)
			GLAI = 0.0;
		DLAI = XNLF / (data.data_p.SLNG - data.data_p.SLNS);
	}
	//dmproduction

	void DMProduction(void)
	{
		//'------------------------------- Parameters and Initials
		if (iniDMP == 0)
		{
			WSFG = 1.0;
			iniDMP = 1;
		}

		//'------------------------------- Adjustment of RUE
		if (TMP <= data.data_p.TBRUE || TMP >= data.data_p.TCRUE)
		{
			TCFRUE = 0.0;
		}
		else if (TMP > data.data_p.TBRUE && TMP < data.data_p.TP1RUE)
		{
			TCFRUE = (TMP - data.data_p.TBRUE) / (data.data_p.TP1RUE - data.data_p.TBRUE);
		}
		else if (TMP > data.data_p.TP2RUE && TMP < data.data_p.TCRUE)
		{
			TCFRUE = (data.data_p.TCRUE - TMP) / (data.data_p.TCRUE - data.data_p.TP2RUE);
		}
		else if (TMP >= data.data_p.TP1RUE && TMP <= data.data_p.TP2RUE)
		{
			TCFRUE = 1.0;
		}

		if (CBD <= data.data_p.ttRUE)
			RUE = data.data_p.IRUE1 * TCFRUE * WSFG;
		else if (CBD > data.data_p.ttRUE)
			RUE = data.data_p.IRUE2 * TCFRUE * WSFG;


		//'------------------------------ CODES FOR RESPONSE TO VPD
		vpdtp = data.data_p.vpd_resp;
		VPDcr = data.data_p.vpd_cr;

		if (vpdtp == 2 || vpdtp == 3) {
			//'__________________________ Hourly calcs ___________________________
			Pi = 3.141592654;
			RDN = Pi / 180.0;
			DEC = sin(23.45 * RDN) * cos(2.0 * Pi * (DOY + 10.0) / 365.0);
			DEC = atan(DEC / sqrt(1.0 - (DEC *  DEC))) * (-1.0);
			DECL = DEC * 57.29578;
			SINLD = sin(RDN * LAT) * sin(DEC);
			COSLD = cos(RDN * LAT) * cos(DEC);
			AOB = SINLD / COSLD;
			AOB2 = atan(AOB / sqrt(1.0 - (AOB *  AOB)));
			DAYL = 12.0 * (1.0 + 2.0 * AOB2 / Pi);
			DSINB = 3600.0 * (DAYL * SINLD + 24.0 * COSLD * sqrt(1.0 - AOB * AOB) / Pi);
			DSINBE = 3600.0 * (DAYL * (SINLD + 0.4 * (SINLD * SINLD + COSLD * COSLD * 0.5)) + 12.0 * COSLD * (2.0 + 3.0 * 0.4 * SINLD) * sqrt(1.0 - AOB * AOB) / Pi);
			SC = 1370.0 * (1.0 + 0.033 * cos(2.0 * Pi * DOY / 365.0));
			DSO = SC * DSINB;
			Bnoon = 90.0 - (LAT - DECL);
			DTR = data.data_h5.srad[ROW] * 1000000.0; //   'from MJ m-2 d-1 to J m-2 d-1
			TC = 4.0;
			P = 1.5;
			SUNRIS = 12.0 - 0.5 * DAYL;
			SUNSET = 12.0 + 0.5 * DAYL;

			if (vpdtp == 2)
			{
				VPTMIN = 0.6108 * exp(17.27 * data.data_h5.tmin[ROW] / (237.3 + data.data_h5.tmin[ROW]));
				TR = 0.0;
			}

			DDMP = 0.0;

			for (int H = 1; H <= 24; H++)
			{
				if ((double)H > SUNRIS && (double)H < SUNSET)//TMINA?????
				{
					if ((double)H < SUNRIS)
					{
						TSUNST = data.data_h5.tmin[ROW] + (TMAXB - data.data_h5.tmin[ROW]) * sin(Pi * (DAYL / (DAYL + 2.0 * P)));
						NIGHTL = 24.0 - DAYL;
						TEMP1 = (data.data_h5.tmin[ROW] - TSUNST * exp(-NIGHTL / TC) + (TSUNST - data.data_h5.tmin[ROW]) * exp(-((double)H + 24.0 - SUNSET) / TC)) / (1.0 - exp(-NIGHTL / TC));
					}
					else if ((double)H < 13.5)
					{
						TEMP1 = data.data_h5.tmin[ROW] + (data.data_h5.tmax[ROW] - data.data_h5.tmin[ROW]) * sin(Pi * ((double)H - SUNRIS) / (DAYL + 2.0 * P));
					}
					else if ((double)H < SUNSET)
					{
						TEMP1 = TMINA + (data.data_h5.tmax[ROW] - TMINA) * sin(Pi * ((double)H - SUNRIS) / (DAYL + 2.0 * P));
					}
					else
					{
						TSUNST = TMINA + (data.data_h5.tmax[ROW] - TMINA) * sin(Pi * (DAYL / (DAYL + 2.0 * P)));
						NIGHTL = 24.0 - DAYL;
						TEMP1 = (TMINA - TSUNST * exp(-NIGHTL / TC) + (TSUNST - TMINA) * exp(-((double)H - SUNSET) / TC)) / (1.0 - exp(-NIGHTL / TC));
					}

					SINB = SINLD + COSLD * cos(2.0 * Pi * ((double)H + 12.0) / 24.0);
					if (SINB < 0.0)
						SINB = 0.0;

					BET = atan(SINB / sqrt(1.0 - (SINB *  SINB)));     //      'BETA IN RADIAN
					BETA = BET * 57.29578;         //       'BETA IN DEGREE
					SRAD1 = DTR * SINB * (1.0 + 0.4 * SINB) / DSINBE; // 'J m-2 s-1
					SRAD1 = SRAD1 * 3600.0 / 1000000.0;      //   'to J m-2 h-1; then MJ m-2 h-1
					if (SRAD1 < 0.0)
						SRAD1 = 0.0;

					FINT = 1.0 - exp(-data.data_p.KPAR * LAI);
					DDMP1 = SRAD1 * 0.48 * FINT * RUE;
				}
				if (vpdtp == 2)
				{
					VPTEMP1 = 0.6108 * exp(17.27 * TEMP1 / (237.3 + TEMP1));
					VPD1 = (VPTEMP1 - VPTMIN) * (param.VPDF / 0.75);//       '(VPDF/0.75) is a correction factor when VPDF is different from 0.75
					TR1 = DDMP1 * VPD1 / data.data_p.TEC;   //'VPD in kPa, TEC in Pa

					if (VPD1 > VPDcr)	//'correction of TR and DBP for VDP>VPDcr
					{
						TR1 = DDMP1 * VPDcr / data.data_p.TEC;
						DDMP1 = TR1 * data.data_p.TEC / VPDcr;
						TR1 = DDMP1 * VPDcr / data.data_p.TEC;
						DDMP1 = TR1 * data.data_p.TEC / VPDcr;
					}
					TR = TR + TR1;
				}
				DDMP = DDMP + DDMP1;
			}
		}
		else if (vpdtp == 1)
		{

			FINT = 1.0 - exp(-data.data_p.KPAR * LAI);
			DDMP = data.data_h5.srad[ROW] * 0.48 * FINT * RUE;
		}

		if (CBD < bdEM || CBD > data.data_p.ttTSG)
			DDMP = 0.0;

	}

	void DMDistribution(void)
	{

		//	'------------------------------- Parameters and Initials
		if (iniDMD == 0)
		{
			WLF = 0.5;
			WST = 0.5;
			WVEG = WLF + WST;
			WGRN = 0.0;
			iniDMD = 1;
		}

		//	'------------------------------- Seed dry matter growth
		if (CBD <= data.data_p.ttBSG) {
			TRANSL = 0.0;
			SGR = 0.0;
			BSGDM = WTOP;   //             'Saving WTOP at BSG
			if (BSGDM <= data.data_p.WDHI1 || BSGDM >= data.data_p.WDHI4)
				DHIF = 0.0;
			else if (BSGDM > data.data_p.WDHI1 && BSGDM < data.data_p.WDHI2)
				DHIF = (BSGDM - data.data_p.WDHI1) / (data.data_p.WDHI2 - data.data_p.WDHI1);
			else if (BSGDM > data.data_p.WDHI3 && BSGDM < data.data_p.WDHI4)
				DHIF = (data.data_p.WDHI4 - BSGDM) / (data.data_p.WDHI4 - data.data_p.WDHI3);
			else if (BSGDM >= data.data_p.WDHI2 && BSGDM <= data.data_p.WDHI3)
				DHIF = 1.0;

			DHI = data.data_p.PDHI * DHIF;
			TRLDM = BSGDM * data.data_p.FRTRL;
		}
		else if (CBD > data.data_p.ttBSG && CBD <= data.data_p.ttTSG)
		{
			SGR = DHI * (WTOP + DDMP) + DDMP * HI;
			if (LAI == 0.0 && NST <= (WST * data.data_p.SNCS))
				SGR = 0.0; // 'There is no N for seed filling
			if (SGR < 0.0)
				SGR = 0.0;

			if ((SGR * data.data_p.GCF) > DDMP)
			{
				TRANSL = (SGR * data.data_p.GCF) - DDMP;
				if (TRANSL > TRLDM)
					TRANSL = TRLDM;
			}
			else if ((SGR * data.data_p.GCF) <= DDMP)
				TRANSL = 0.0;

			TRLDM = TRLDM - TRANSL;
			if (SGR > ((DDMP + TRANSL) / data.data_p.GCF))
				SGR = (DDMP + TRANSL) / data.data_p.GCF;
		}
		else if (CBD > data.data_p.ttTSG)
		{
			TRANSL = 0.0;
			SGR = 0.0;
		}
		//'------------------------------- DM avail. for Leaf & stem
		DDMP2 = DDMP - SGR * data.data_p.GCF;
		if (DDMP2 < 0.0)
			DDMP2 = 0.0;

		//	'------------------------------- Leaf dry matter growth
		if (CBD <= bdBLG || CBD > bdTLP)
			GLF = 0.0;
		else if (CBD > bdBLG && CBD <= bdTLM)
		{
			if (WTOP < data.data_p.WTOPL)
				FLF1 = data.data_p.FLF1A;
			else
				FLF1 = data.data_p.FLF1B;
			GLF = FLF1 * DDMP2;
		}
		else if (CBD > bdTLM && CBD <= bdTLP)
			GLF = data.data_p.FLF2 * DDMP2;

		//	'------------------------------- Stem dry matter growth
		GST = DDMP2 - GLF;

		//'------------------------------- Organs accumulated mass
		WLF = WLF + GLF;
		WST = WST + GST;
		WGRN = WGRN + SGR;
		WVEG = WVEG + DDMP - (SGR * data.data_p.GCF);
		WTOP = WVEG + WGRN;
		HI = WGRN / WTOP;
	}

	void LegumPlant(void)
	{
		if (iniPNB == 0)
		{
			NST = WST * data.data_p.SNCG;
			NLF = LAI * data.data_p.SLNG;
			WSFN = 1.0;
			CNUP = NST + NLF;
			NGRN = 0.0;
			iniPNB = 1;

		}
		if (CBD <= bdEM || CBD > data.data_p.ttTSG)
		{
			NUP = 0.0;
			XNLF = 0.0;
			XNST = 0.0;
			INLF = 0.0;
			INST = 0.0;
			INGRN = 0.0;
		}

		else if (CBD > bdEM && CBD < data.data_p.ttBSG)
		{
			INGRN = 0.0;
			NUP = (GST * data.data_p.SNCG) + (GLAI * data.data_p.SLNG); // '+ NSTDF    '<---- -
			if (CBD < data.data_p.ttBNF && CNUP > param.INSOL)
				NUP = 0.0;
			if (NUP > data.data_p.MXNUP)
				NUP = data.data_p.MXNUP;

			NFC = NFC * 3.0 / 4.0 + NUP / WVEG * (1.0 / 4.0);//   'from Sinclair et al. 2003
			NUP = NUP * WSFN;
			if (NUP < 0.0)
				NUP = 0.0;
			// 'If FTSW > 1 Then NUP = 0  'Inactivated
			if (DDMP == 0.0)
				NUP = 0.0;

			if (NST <= (WST * data.data_p.SNCS))
			{
				INST = WST * data.data_p.SNCS - NST;
				XNST = 0.0;
				if (INST >= NUP)
				{
					INLF = 0.0;
					XNLF = INST - NUP;
				}
				else if (INST < NUP)
					INLF = GLAI * data.data_p.SLNG;
				if (INLF > (NUP - INST))
					INLF = NUP - INST;
				INST = NUP - INLF;
				XNLF = 0.0;
			}
			else if (NST > (WST * data.data_p.SNCS))
			{

				INLF = GLAI * data.data_p.SLNG;
				XNLF = 0.0;
			}
			if (INLF >= NUP)
			{
				INST = 0.0;
				XNST = INLF - NUP;
				if (XNST > (NST - WST * data.data_p.SNCS))
					XNST = NST - WST * data.data_p.SNCS;
				INLF = NUP + XNST;
			}
			else if (INLF < NUP)
			{
				INST = NUP - INLF;
				XNST = 0.0;
			}

		}
		//'       TRLN = LAI * (SLNG - SLNS) + (NST + INST - WST * SNCS)
		//'       FXLF = LAI * (SLNG - SLNS) / TRLN

		else if (CBD >= data.data_p.ttBSG && CBD <= data.data_p.ttTSG)
		{
			INGRN = SGR * data.data_p.GNC;
			NUP = INGRN + (GST * data.data_p.SNCG) + (GLAI * data.data_p.SLNG);
			if (NUP > data.data_p.MXNUP)
				NUP = data.data_p.MXNUP;
			PDNF = NFC * WVEG;
			if (PDNF > NUP)
				PDNF = NUP;
			DNF = PDNF * WSFN;
			//   'If FTSW > 1 Then DNF = 0  'Inactivated
			if (DNF < 0.0)
				DNF = 0.0;
			if (DDMP <= (SGR * data.data_p.GCF))
				DNF = 0.0;
			if (DDMP == 0.0)
				DNF = 0.0;
			NUP = DNF;

			if (NUP > (SGR * data.data_p.GNC))
			{
				//   'N is excess of seed needs
				NUP2 = NUP - SGR * data.data_p.GNC;
				if (NST <= (WST * data.data_p.SNCS))
				{
					INST = WST * data.data_p.SNCS - NST;
					XNST = 0.0;

					if (INST >= NUP2)
					{
						INLF = 0.0;
						XNLF = INST - NUP2;
					}
					else if (INST < NUP2)
					{
						INLF = GLAI * data.data_p.SLNG;
						if (INLF > (NUP2 - INST))
						{
							INLF = NUP2 - INST;
							INST = NUP2 - INLF;
							XNLF = 0.0;
						}
					}
				}
				else if (NST > (WST * data.data_p.SNCS))
				{
					INLF = GLAI * data.data_p.SLNG;
					XNLF = 0.0;
					if (INLF >= NUP2)
					{
						INST = 0.0;
						XNST = INLF - NUP2;

						if (XNST > (NST - WST * data.data_p.SNCS))
							XNST = NST - WST * data.data_p.SNCS;
						INLF = NUP2 + XNST;
					}
					else if (INLF < NUP2)
					{
						INST = NUP2 - INLF;
						XNST = 0.0;
					}

				}

				TRLN = LAI * (data.data_p.SLNG - data.data_p.SLNS) + (NST + INST - WST * data.data_p.SNCS);
				FXLF = LAI * (data.data_p.SLNG - data.data_p.SLNS) / TRLN;
			}
			else if (NUP <= (SGR * data.data_p.GNC))
			{
				/// 'Need to transfer N from vegetative tissue
				INLF = 0.0;
				INST = 0.0;
				XNLF = (SGR * data.data_p.GNC - NUP) * FXLF;
				XNST = (SGR * data.data_p.GNC - NUP) * (1.0 - FXLF);
			}

		}
		NST = NST + INST - XNST;
		NLF = NLF + INLF - XNLF;
		NVEG = NLF + NST;
		NGRN = NGRN + INGRN;
		CNUP = CNUP + NUP;

		TRLN = LAI * (data.data_p.SLNG - data.data_p.SLNS) + (NST - WST * data.data_p.SNCS);
		FXLF = LAI * (data.data_p.SLNG - data.data_p.SLNS) / (TRLN + 0.000000000001);///////////
		if (FXLF > 1.0)
			FXLF = 1.0;
		if (FXLF < 0.0)
			FXLF = 0.0;
	}
	void DailyPrintOut(void)
	{
		
		//out_LAI.open(param.func_file_name.toStdString() + "_" + "output_LAI" + param.h5_file_name.toStdString() + ".txt", std::ios::app);
		out_LAI.open(param.h5_table_name.toStdString() + "_" + "output_LAI.txt", std::ios::app);
		out_DAP.open(param.h5_table_name.toStdString() + "_" + "DAP.txt" , std::ios::app);
		out_FTSW.open(param.h5_table_name.toStdString() + "_" + "output_ftsw.txt" , std::ios::app);
		out_CE.open(param.h5_table_name.toStdString() + "_" + "output_CE.txt", std::ios::app);
		out_CTR.open(param.h5_table_name.toStdString() + "_" + "output_CTR.txt", std::ios::app);
		out_MSNN.open(param.h5_table_name.toStdString() + "_" + "output_MSNN.txt", std::ios::app);
		out_CNUP.open(param.h5_table_name.toStdString() + "_" + "output_CNUP.txt", std::ios::app);
		out_NVEG.open(param.h5_table_name.toStdString() + "_" + "output_NVEG.txt", std::ios::app);
		out_NGRN.open(param.h5_table_name.toStdString() + "_" + "output_NGRN.txt", std::ios::app);

		out.open(param.h5_table_name.toStdString() + "_" + "output.txt", std::ios::app);
		if (write_check == false) {
			out <<"GEO_ID" <<"NSAM= " << "ROW= " << " YEARS= " << " DOY= " << "DAP= " << " TMP= " << " DTT= " << "CBD= " << "MSNN= " << "GLAI= " << " DLAI=" << "LAI = " << "TCFRUE = " << "FINT= " << " DDMP= " << "GLF= " << "GST= " << " SGR= " << " WLF= " << "WST= " << "WVEG=  " << " WGRN= " << "WTOP= " << "DEPORT= " << "RAIN= " << "IRGW= " << "RUNOF= " << "PET= " << "SEVP= " << "TR= " << "ATSW= " << " FTSW=  " << " CRAIN= " << " CIRGW= " << "IRGNO= " << " CRUNOF= " << "CE= " << " CTR= " << "WSTORG= " << "NUP= " << " NLF= " << "NST= " << "NVEG= " << " NGRN= " << "CNUP= " << "MAT= " << endl;
			write_check = true;
		}
		out_LAI << LAI << endl;
		out_DAP << DAP << endl;
		out_FTSW << FTSW << endl;
		out_CE << CE << endl;
		out_CTR << CTR << endl;
		out_MSNN << MSNN << endl;
		out_CNUP << CNUP << endl;
		out_NVEG << NVEG << endl;
		out_NGRN << NGRN << endl;

		out << data.data_a5.geo_id[NSAM] << " ";
		out << NSAM << "  ";
		out << ROW << "  ";
		out << data.data_h5.years[ROW] << "  ";
		out << data.data_h5.doy[ROW] << "  ";
		out << DAP << "  ";
		out << TMP << "  ";
		out << DTT << "  ";
		out << CBD << "  ";
		out << MSNN << "  ";
		out << GLAI << "  ";
		out << DLAI << "  ";
		out << LAI << "  ";
		out << TCFRUE << "  ";
		out << FINT << "  ";
		out << DDMP << "  ";
		out << GLF << "  ";
		out << GST << "  ";
		out << SGR << "  ";
		out << WLF << "  ";
		out << WST << "  ";
		out << WVEG << "  ";
		out << WGRN << "  ";
		out << WTOP << "  ";
		out << data.data_p.DEPORT << "  ";
		out << data.data_h5.rain[ROW] << "  ";
		out << IRGW << "  ";
		out << RUNOF << "  ";
		out << PET << "  ";
		out << SEVP << "  ";
		out << TR << "  ";
		out << ATSW << "  ";
		out << FTSW << "  ";
		out << CRAIN << "  ";
		out << CIRGW << "  ";
		out << IRGNO << "  ";
		out << CRUNOF << "  ";
		out << CE << "  ";
		out << CTR << "  ";
		out << WSTORG << "  ";
		out << NUP << "  ";
		out << NLF << "  ";
		out << NST << "  ";
		out << NVEG << "  ";
		out << NGRN << "  ";
		out << CNUP << "  ";
		out << MAT << endl;
		out.close();
		out_LAI.close();
		out_DAP.close();
		out_FTSW.close();
		out_CE.close();
		out_CTR.close();
		out_MSNN.close();
		out_CNUP.close();
		out_NVEG.close();
		out_NGRN.close();
	}
	void SummaryPrintOut()
	{
		out_s.open(param.func_file_name.toStdString() + "_" + "output_summary.txt", std::ios::app);
		//out_s.open(param.func_file_name.toStdString() + "_" + "output_summary.txt", std::ios::app);
		out_s << "[" << NSAM << "]" << endl;
		out_s <<"geo_id = " << data.data_a5.geo_id[NSAM] << endl;
		out_s << "CBD = " << CBD << endl;
		out_s << "year = " << data.data_h5.years[ROW] << endl;
		out_s << "dtEM = " << dtEM << endl;
		out_s << "dtR1 = " << dtR1 << endl;
		out_s << "dtR3 = " << dtR3 << endl;
		out_s << "dtR5 = " << dtR5 << endl;
		out_s << "dtR7 = " << dtR7 << endl;
		out_s << "dtR8 = " << dtR8 << endl;
		out_s << "cbdR1 = " << cbdR1 << endl;
		out_s << "MSNN = " << MSNN << endl;
		out_s << "MXLAI = " << MXLAI << endl;
		out_s << "BSGLAI = " << BSGLAI << endl;
		out_s << "BSGDM = " << BSGDM << endl;
		out_s << "WTOP = " << WTOP << endl;
		out_s << "WGRN = " << WGRN << endl;
		out_s << " WGRN / WTOP * 100 = " << WGRN / WTOP * 100 << endl;
		out_s << "WGRN / data.data_p.TRESH = " << WGRN / data.data_p.TRESH << endl;
		out_s << "ISATSW = " << ISATSW << endl;
		out_s << "CRAIN = " << CRAIN << endl;
		out_s << "CIRGW = " << CIRGW << endl;
		out_s << "IRGNO = " << IRGNO << endl;
		out_s << "ATSW = " << ATSW << endl;
		out_s << "CRUNOF = " << CRUNOF << endl;
		out_s << "CE = " << CE << endl;
		out_s << "CTR" << CTR << endl;
		out_s << "WSTORG = " << WSTORG << endl;
		out_s << "CE + CTR = " << CE + CTR << endl;
		out_s << "CE / (CE + CTR + 0.00001) = " << CE / (CE + CTR + 0.00001) << endl;
		out_s << "NLF = " << NLF << endl;
		out_s << "(NLF + NST) = " << (NLF + NST) << endl;
		out_s << " NGRN = " << NGRN << endl;
		out_s << "CNUP = " << CNUP << endl;
		out_s << " param.INSOL = " << param.INSOL << endl;
		out_s << "CIRGW = " << CIRGW << endl;
		out_s.close();
	}
	double Heaviside(double arg) { return (arg > 0.0) ? 1.0 : 0.0; }
	void calculation(void)
	{
	//	arma::mat Fpoints;
	//	return;
		double phase_change;
		double cbd;
		double training_error = 0;
		double curr_error = 0;
		int nDays = param.nD;
		int START_YEAR_TO_PRINT;
	/*	string last_name_file = param.h5_table_name.toStdString();
		string res;
		if (last_name_file.length() == 108)
		{
			char o = last_name_file[last_name_file.length() - 4];
			string one(1, o);
			res = one;
		}
		else if (last_name_file.length() == 109)
		{
			char o = last_name_file[last_name_file.length() - 5];
			char t = last_name_file[last_name_file.length() - 4];
			string one(1, o);
			string two(1, t);
			res = one + two;
		}*/
	//	cout << last_name_file << endl;
	//	cout << last_name_file.length() << endl;
	//	cout << res << endl;
	//	cout << nDays << endl;
	//	nl->print_trace(param.func_file_name.toStdString(), 0);
	//	double * interpol_error = new double[data.data_a5.nSamples];
	//	if (param.print_trace > 0 && param.function_mode != SOLTANI_FUNC)
	//	{
	//		nl->print_trace(param.func_file_name.toStdString(), 0);
	//	}			
	//	cout << data.data_a5.nSamples << endl;
	//	cout << data.data_a5.nSamples << endl;
		 for (size_t nsam = 0; nsam < data.data_a5.nSamples; ++nsam)
		 { 
		//	 cout << endl;
		//	cout << "NEW CHICKPEA - NSAM: " <<  nsam << endl;
			NSAM = nsam;
			ROW = -1;
			MAT = 0;
			CBD = 0.0;
			iniPheno = 0;
			iniLai = 0;
			iniDMP = 0;
			iniDMD = 0;
			iniSW = 0;
			iniPNB = 0;
			//DAP = 0.0;
			
			int start_day = data.data_a5.doy[nsam];
			int start_year = data.data_a5.years[nsam];
			int geo_id = data.data_a5.geo_id[nsam];
			double event_day = data.data_a5.response[nsam];
		//	out_bd.open(param.h5_table_name.toStdString() + "_" + "bd.txt", std::ios::app);
		//	out_bd << "chickpea = " << nsam << endl;
		//	for (int k = 0; k < 18; k++)
		//		out_bd << data.data_a5.gr_covar[NSAM][k] << " ";
		//	out_bd << endl;
		//	out_bd << "start_day = " << start_day << endl;
		//	out_bd << "start_year = " << start_year << endl;
		//	out_bd << "geo_id = " << geo_id << endl;

		//	out_bd.close();
		//	cout << "start_day = " << start_day << endl;
		//	cout << "start_year = " << start_year << endl;
		//	cout << "geo_id = " << geo_id << endl;

			//	FindSowingData(geo_id, start_day, start_year);
			bool check = false;
			size_t j = 0;
			START_YEAR_TO_PRINT = start_year;
			for (size_t nd = 0; nd < data.data_h5.nWeather; ++nd) {
				if (geo_id == data.data_h5.geo_id[nd] && data.data_h5.doy[nd] == start_day && data.data_h5.years[nd] == start_year) {
					j = nd;
					check = true;
					break;
				}
			}
			ROW = j;
		//	cout << "ROW = " << ROW << endl;
			if (check == false)
			{
				continue;
			}
 		//	if (param.print_trace > 0) cout << "nsam = " << nsam << " geo id = " << geo_id << " BEGIN ROW = " << ROW;
			int curr_day = -nDays;
			for (size_t nd = 0; nd < nDays; nd++)
			{
				ROW = j + nd;
				if (ROW >= data.data_h5.nWeather)
					break;

				Weather();
				Phenology();
				CropLAIN();
				DMProduction();
				DMDistribution();
				LegumPlant();
				SoilWater();
			//	if (param.print_trace > 0)
			//	{
			//		DailyPrintOut();
			//	}
				if (MAT == 1)
				{
					break;
				}
		//		double rhs = 0.1 * ((double)nd + 6.0 - event_day);
			//	cout << "CBD = " << CBD << endl;
		//		rhs = (rhs > 0.9) ? 0.9 : rhs;
		//		rhs = (rhs < 0.1) ? 0.1 : rhs;
		//		cout << "RHS = " << rhs << endl;
		//		interpol_error[nsam] += (CBD - rhs) * (CBD - rhs);
			}
			phase_change = nl->get_cbd();//get_phase_change();
			cbd = get_cbd();
			curr_day = get_curr_day();
			//	last_name_file = last_name_file[last_name_file.length() - 4];
	  //	last_name_file = last_name_file[last_name_file.length() - 4];
		//	out_cbd.open(param.h5_table_name.toStdString() + "_"+ "cbd.txt", std::ios::app);
		//	out_cbd << cbd << endl;
		//	out_cbd.close();
			check_last_cbd = false;
			//		cout << "    TEMP DATA: " << endl;		
	
			//cout << last_name_file << endl;
		  //  out_dtR1.open( param.h5_file_name.toStdString() +  "_" + "dtR1.txt", std::ios::app);
			out_dtR1.open(param.func_file_name.toStdString()+param.h5_file_name.toStdString() + param.h5_table_name.toStdString()+ "_" + "dtR1.txt", std::ios::app);
		//	cout << param.h5_file_name.toStdString() + "_" + to_string(start_year) + "_" + "dtR1.txt" << endl;
			out_dtR1 << curr_day << endl;
			out_dtR1.close();
			
		//	SummaryPrintOut();

		//	if (param.print_trace > 0)
		//	{
			//	SummaryPrintOut();
			//	cout << " END ROW = " << ROW;
			//	cout << " MAT = " << MAT;
			//	cout << " curr_day = " << curr_day;
			//	cout << " event_dat = " << event_day;
			//	cout << " CBD = " << cbd;
			//	cout << " phase_change = " << phase_change << endl;
		//	}
		//	cout << "answer:" << endl;
		//	cout << curr_day << endl;
		//	cout << event_day << endl;

			training_error += (curr_day - event_day) * (curr_day - event_day);
			curr_error += (cbd - phase_change) * (cbd - phase_change);
		}
	//	cout << endl << "final print:" << endl;
		//out_error.open(param.h5_file_name.toStdString() + "_" + "ERROR.txt", std::ios::app);
		out_error.open(param.func_file_name.toStdString() + param.h5_file_name.toStdString() + param.h5_table_name.toStdString()  + "_" + "ERROR.txt", std::ios::app);

		out_error << training_error << endl;
		out_error << curr_error << endl;
		out_error << nl->get_l1_pen() << endl;
      	out_error << nl->get_l2_pen() << endl;
		out_error.close();
		cout << training_error << endl;////////////////////////////////////////////////
		cout << curr_error << endl;
		cout << nl->get_l1_pen() << endl;//////////////////////
		cout << nl->get_l2_pen() << endl;/////////////////////////////////////
		nl->print_trace(param.func_file_name.toStdString(), 0);
	//	if (param.print_trace > 0 && param.function_mode != SOLTANI_FUNC)
	//	{
	//		nl->print_trace(param.func_file_name.toStdString(), 0);
	//	}
	//	delete [] interpol_error;
	//	delete nl;
		nl->delete_all();
	}
public slots:
	void run_h5()
	{

		//if (param.print_trace > 0) cout << "BEGIN CALC" << endl;//////////////////////////////
		calculation();
		//if (param.print_trace > 0) cout << "END CALC" << endl;///////////////////////////////
		emit finished();
	}
signals:
	void finished();
};

#endif // _MODEL_H_
