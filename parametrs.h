#pragma once

class Parametrs
{
	
public:

	QString file_name;
	QString sName;//Name of site or scenario. This name will appear in output
	double Latitude;//Latitude of location. NOTE: south latitude should be NEGETIVE
	double VPDF;//Vapor pressure deficit factor.
	double yno; // Number of years of simulation for each site/scenario
	int FirstYear; // First year of weather data.
	int FixFind; // 1=fixed sowing date, 2=first possible sowing date based on pre-defined conditions.
	int Pdoy; // Sowing date as Day of Year if Fix/Find=1 and day of year that model must start searching for sowing date from if Fix/Find=2.
	double SearchDur; // Specifies the end of sowing window
	double SowWat; // Amount of water in the soil that makes sowing possible (mm) 
	double PDEN; // Plant density (# per m2)
	int water; //  1= irrigated 2 = rainfed
	double IRGLVL; // Irrigation level based on FTSW if the farming is irrigated. 
	double SOLDEP; // Soil depth (mm)
	double DEP1; // Depth of top layer (mm)
	double SALB; // Soil albedo
	double CN2; // Soil curve number. 
	double DRAINF; // Drainage factor
	double SAT; // Soil saturation limit (m3 m-3 or mm mm-1).
	double DUL;// Soil drained upper limit (m3 m-3 or mm mm-1).
	double EXTR; // Soil extractable moisture (m3 m-3 or mm mm-1).
	int MAI1; // Soil moisture availabilty index for top layer. 1=DUL, 0=LL. 
	int MAI;// Soil moisture availability index (1=DUL, 0=LL). 
	int INSOL; // Initial soil nitrogen that can be uptaken by the crop before BNF is activated (g N m-2)
	double U;// Cumulative evaporation that stage I terminates and stage II starts. Not used by the model (inactive)
	int CropColNo; // Column number in "Crops" sheeet that contains cultivar parameters estimates
	int ROW;
	bool file_mode;// true - with dl or false - without dl
	Parametrs() {}
	void Print()
	{
		std::cout << Latitude << std::endl;
		std::cout << VPDF << std::endl;
		std::cout << yno << std::endl;
		std::cout << FirstYear << std::endl;
		std::cout << FixFind << std::endl;
		std::cout << Pdoy << std::endl;
		std::cout << SearchDur << std::endl;
		std::cout << SowWat << std::endl;
		std::cout << PDEN << std::endl;
		std::cout << water << std::endl;
		std::cout << IRGLVL << std::endl;
		std::cout << SOLDEP << std::endl;
		std::cout << DEP1 << std::endl;
		std::cout << SALB << std::endl;
		std::cout << CN2 << std::endl;
		std::cout << DRAINF << std::endl;
		std::cout << SAT << std::endl;
		std::cout << DUL << std::endl;
		std::cout << EXTR << std::endl;
		std::cout << MAI1 << std::endl;
		std::cout << MAI << std::endl;
		std::cout << INSOL << std::endl;
		std::cout << U << std::endl;
		std::cout << CropColNo << std::endl;
		std::cout << file_mode << endl;
	}
};
