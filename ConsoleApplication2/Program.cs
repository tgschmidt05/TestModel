using System;
using System.Windows;
using System.IO;
using System.Collections.Generic;
using System.Windows.Markup;
using OxyPlot;
using OxyPlot.Series;


//////////////////////////////
///////////////////////////////
//PLEASE NOTE: This model can only operate when no precipitation occurs during running period and atmosphere remains statically neutral.
//Does not account for changing surface soil temperatures or solar radiation.
///////////////////////////////
///////////////////////////////
namespace ConsoleApplication2
{
    public class GridPoint
    {
        private Boolean maxIBoundaryPoint;
        private Boolean minIBoundaryPoint;
        private Boolean maxJBoundaryPoint;
        private Boolean minJBoundaryPoint;
        private Boolean maxKBoundaryPoint;
        private Boolean minKBoundaryPoint;
        private Boolean ghostCell;
        private double longitude; //originally separated by 0.25 degrees
        private double latitude; //originally separated by 0.25 degrees
        private int iCoord;
        private int jCoord;
        private int kCoord;
        private double deltaX; // m (Multiplier of 8)
        private double deltaY; // m (Multiplier of 16)
        private double deltaZ; // m (needs conversion from pa)
        private double geopotentialHeight; // gpm
        private double netRadiation; // w/m^2
        private double sfcTemp; // k (soil temperature) <CONSTANT>
        private double nearSfcPressure; // pa (Consider how this changes)
        private double turbTemp;
        private double turbMixingRatio;
        private double turbUWind;
        private double turbVWind;
        private double relativeHumidity; // 0 to 100
        private double verticalVelocity; // Pa/s
        private double windU; // m/s
        private double windV; // m/s
        private double windW; // m/s (needs conversion from pa/s)
        private double airDensity; // kg/m^3
        private double temperature; // k
        private double airPressure; // pa
        private double mixingRatio; // kg/kg
        //private double planetaryBoundaryLayer; //height in m

        public GridPoint(int xLen, int yLen, int zLen, int iCoord, int jCoord, int kCoord, double windU, double windV, 
            double windW, double airDensity, double temperature, double airPressure, double mixingRatio,double deltaX, double deltaY,
            double deltaZ, double latitude, double longitude, double geopotentialHeight, double netRadiation, double sfcTemp,
            double nearSfcPressure, double relativeHumidity, double verticalVelocity)
        {
            this.iCoord = iCoord;
            this.jCoord = jCoord;
            this.kCoord = kCoord;
            this.windU = windU;
            this.windV = windV;
            this.windW = windW;
            this.airDensity = airDensity;
            this.temperature = temperature;
            this.airPressure = airPressure;
            this.mixingRatio = mixingRatio;
            maxIBoundaryPoint = (iCoord == (xLen - 2)); //Takes into account ghost cells
            minIBoundaryPoint = (iCoord == 1); //Takes into account ghost cells
            maxJBoundaryPoint = (jCoord == (yLen - 2)); //Takes into account ghost cells
            minJBoundaryPoint = (jCoord == 1); //Takes into account ghost cells
            maxKBoundaryPoint = (kCoord == (zLen - 2)); //Takes into account ghost cells
            minKBoundaryPoint = (kCoord == 0); //Takes into account ghost cells
            ghostCell = ((iCoord == 0) || (iCoord == (xLen - 1)) || (jCoord == 0) || (jCoord == (yLen - 1)) || (kCoord == (zLen - 1)));
            this.deltaX = deltaX;
            this.deltaY = deltaY;
            this.deltaZ = deltaZ;
            this.latitude = latitude;
            this.longitude = longitude;
            this.geopotentialHeight = geopotentialHeight;
            this.netRadiation = netRadiation;
            this.sfcTemp = sfcTemp;
            this.nearSfcPressure = nearSfcPressure;
            this.relativeHumidity = relativeHumidity;
            this.verticalVelocity = verticalVelocity;
        }

        //This constructor allows for cases where blank gridpoints need to be made
        public GridPoint(int xLen, int yLen, int zLen, int iCoord, int jCoord, int kCoord)
        {
            this.iCoord = iCoord;
            this.jCoord = jCoord;
            this.kCoord = kCoord;
            maxIBoundaryPoint = (iCoord == (xLen - 2)); //Takes into account ghost cells
            minIBoundaryPoint = (iCoord == 1); //Takes into account ghost cells
            maxJBoundaryPoint = (jCoord == (yLen - 2)); //Takes into account ghost cells
            minJBoundaryPoint = (jCoord == 1); //Takes into account ghost cells
            maxKBoundaryPoint = (kCoord == (zLen - 2)); //Takes into account ghost cells
            minKBoundaryPoint = (kCoord == 0); //Takes into account ghost cells
            ghostCell = ((iCoord == 0) || (iCoord == (xLen - 1)) || (jCoord == 0) || (jCoord == (yLen - 1)) || (kCoord == (zLen - 1)));
        }

        public void inputUpdatedData(int xLen, int yLen, int zLen, int iCoord, int jCoord, int kCoord, double windU,
            double windV,
            double windW, double airDensity, double temperature, double airPressure, double mixingRatio, double deltaX,
            double deltaY,
            double deltaZ, double latitude, double longitude, double geopotentialHeight, double netRadiation,
            double sfcTemp,
            double nearSfcPressure, double relativeHumidity, double verticalVelocity)
        {
            this.iCoord = iCoord;
            this.jCoord = jCoord;
            this.kCoord = kCoord;
            this.windU = windU;
            this.windV = windV;
            this.windW = windW;
            this.airDensity = airDensity;
            this.temperature = temperature;
            this.airPressure = airPressure;
            this.mixingRatio = mixingRatio;
            maxIBoundaryPoint = (iCoord == (xLen - 2)); //Takes into account ghost cells
            minIBoundaryPoint = (iCoord == 1); //Takes into account ghost cells
            maxJBoundaryPoint = (jCoord == (yLen - 2)); //Takes into account ghost cells
            minJBoundaryPoint = (jCoord == 1); //Takes into account ghost cells
            maxKBoundaryPoint = (kCoord == (zLen - 2)); //Takes into account ghost cells
            minKBoundaryPoint = (kCoord == 0); //Takes into account ghost cells
            ghostCell = ((iCoord == 0) || (iCoord == (xLen - 1)) || (jCoord == 0) || (jCoord == (yLen - 1)) || (kCoord == (zLen - 1)));
            this.deltaX = deltaX;
            this.deltaY = deltaY;
            this.deltaZ = deltaZ;
            this.latitude = latitude;
            this.longitude = longitude;
            this.geopotentialHeight = geopotentialHeight;
            this.netRadiation = netRadiation;
            this.sfcTemp = sfcTemp;
            this.nearSfcPressure = nearSfcPressure;
            this.relativeHumidity = relativeHumidity;
            this.verticalVelocity = verticalVelocity;
        }

        public Boolean isTarget(int inputI, int inputJ, int inputK)
        {
            return (inputI == iCoord && inputJ == jCoord && inputK == kCoord);
        }
        public double getNetRadiation()
        {
            return netRadiation;
        }
        public double getDeltaX()
        {
            return deltaX;
        }
        public double getDeltaY()
        {
            return deltaY;
        }
        public double getDeltaZ()
        {
            return deltaZ;
        }
        public double getTemp()
        {
            return temperature;
        }
        public double getPressure()
        {
            return airPressure;
        }
        public double getMixingRatio()
        {
            return mixingRatio;
        }
        public double getAirDensity()
        {
            return airDensity;
        }
        public double getWindU()
        {
            return windU;
        }
        public double getWindV()
        {
            return windV;
        }
        public double getWindW()
        {
            return windW;
        }
        public double getWindSpeed()
        {
            double uSqaured = Math.Pow(windU, 2);
            double vSqaured = Math.Pow(windV, 2);
            double wSqaured = Math.Pow(windW, 2);
            return Math.Pow((uSqaured+vSqaured+wSqaured),0.5);
        }
        public double getSfcTemp()
        {
            return sfcTemp;
        }
        public double getTurbTemp()
        {
            return turbTemp;
        }
        public Boolean isOnBottomBoundary()
        {
            return minKBoundaryPoint;
        }
        public Boolean isOnTopBoundary()
        {
            return maxKBoundaryPoint;
        }
        public Boolean isOnLeftBoundary()
        {
            return maxJBoundaryPoint;
        }
        public Boolean isOnRightBoundary()
        {
            return minJBoundaryPoint;
        }
        public Boolean isOnFrontBoundary()
        {
            return maxIBoundaryPoint;
        }
        public Boolean isOnBackBoundary()
        {
            return minIBoundaryPoint;
        }
        public Boolean isGhostCell()
        {
            return ghostCell;
        }
        public double getRelativeHumidity()
        {
            return relativeHumidity;
        }
        public double getNearSfcPressure()
        {
            return nearSfcPressure;
        }
        public void setTurbTemp(double input)
        {
            turbTemp = input;
        }
        public void setTemp(double temp)
        {
            this.temperature = temp;
        }
        public void setMixingRatio(double mixingRatio)
        {
            this.mixingRatio = mixingRatio;
        }
        public void setAirDensity(double airDensity)
        {
            this.airDensity = airDensity;
        }
        public void setAirPressure(double pressure)
        {
            this.airPressure = pressure;
        }
        public void setWindU(double windU)
        {
            this.windU = windU;
        }
        public void setWindV(double windV)
        {
            this.windV = windV;
        }
        public void setWindW(double windW)
        {
            this.windW = windW;
        }

        public int getICoord()
        {
            return iCoord;
        }
        public int getJCoord()
        {
            return jCoord;
        }
        public int getKCoord()
        {
            return kCoord;
        }

        public void setTurbMixingRatio(double turbMixingRatio)
        {
            this.turbMixingRatio = turbMixingRatio;
        }

        public void setTurbUWind(double turbUWind)
        {
            this.turbUWind = turbUWind;
        }
        
        public void setTurbVWind(double turbVWind)
        {
            this.turbVWind = turbVWind;
        }

        public void setDeltaXYZ(double x, double y, double z)
        {
            this.deltaX = x;
            this.deltaY = y;
            this.deltaZ = z;
        }

        public void setLatLon(double lat, double lon)
        {
            this.latitude = lat;
            this.longitude = lon;
        }
        public void setGeopotentialHeight(double geopotential)
        {
            this.geopotentialHeight = geopotential;
        }
        public void setNetRadiation(double netRadiation)
        {
            this.netRadiation = netRadiation;
        }
        public void setSoilTemperature(double soilTemp)
        {
            this.sfcTemp = soilTemp;
        }
        public void setNearSfcPressure(double sfcPressure)
        {
            this.nearSfcPressure = sfcPressure;
        }
        public void setRelativeHumidity(double RH)
        {
            this.relativeHumidity = RH;
        }
        public void setVerticalVelocity(double verticalVelocity)
        {
            this.verticalVelocity = verticalVelocity;
        }
        public double getTurbMixingRatio()
        {
            return turbMixingRatio;
        }

        public double getTurbUWind()
        {
            return turbUWind;
        }

        public double getTurbVWind()
        {
            return turbVWind;
        }

        public double getLatitude()
        {
            return latitude;
        }

        public double getLongitude()
        {
            return longitude;
        }

        public double getGeopotentialHeight()
        {
            return geopotentialHeight;
        }

        public double getVerticalVelocity()
        {
            return verticalVelocity;
        }
    }

    public class DataReceiver
    {
        public const double lengthQuarterDegreeLon = 22000; //For Kansas
        public const double lengthQaurterDegreeLat = 28000; //For Kansas
        public const double a = 0.61; //g/g
        public const double Rd = 287.053; //dry-air gas constant
        private int resolutionMultiplierX; //Determines how many gridpoints will be placed between the 0.25 lon points
        private int resolutionMultiplierY; //Determines how many gridpoints will be placed between the 0.25 lat points
        private int resolutionMultiplierZ; //Determines how many gridpoints will be placed between the 25 mb z points
        private int xLen;
        private int yLen;
        private int zLen;
        private double Zi; //ABL depth
        private double bottomLeftLat;
        private double bottomLeftLon;
        private double[,,] deltaX;
        private double[,,] deltaY;
        private double[,,] deltaZ;
        private double[,,] longitudeData;
        private double[,,] latitudeData;
        private double[,,] geopotentialHeightData;
        private double[,,] netRadiationData;
        private double[,,] soilTemperatureData;
        private double[,,] groundAirPressData;
        private double[,,] relativeHumidityData;
        private double[,,] verticalVelocityData;
        private double[,,] tempData; // k
        private double[,,] pressureData; // pa
        private double[,,] mixingRatioData; // g/kg
        private double[,,] densityData; // kg/m^3
        private double[,,] windUData; // m/s
        private double[,,] windVData; // m/s
        private double[,,] windWData; // m/s

        public DataReceiver(double Zi, int resolutionMultiplierX, int resolutionMultiplierY, int resolutionMultiplierZ, double bottomLeftLat, double bottomLeftLon)
        {
            this.Zi = Zi;
            this.bottomLeftLat = bottomLeftLat;
            this.bottomLeftLon = bottomLeftLon;
            this.resolutionMultiplierX = resolutionMultiplierX;
            this.resolutionMultiplierY = resolutionMultiplierY;
            this.resolutionMultiplierZ = resolutionMultiplierZ;
        }
        
        private List<string> loadCsvFile(string filePath)
        {
            var reader = new StreamReader(File.OpenRead(filePath));
            List<string> searchList = new List<string>();
            while(!reader.EndOfStream)
            {
                var line = reader.ReadLine();
                searchList.Add(line);
            }
            return searchList;
        }

        public void inputDataSet(Boolean isAnalysisData)
        {
            string appPath = System.IO.Path.GetDirectoryName(
                System.Reflection.Assembly.GetExecutingAssembly().Location);
            string appPath1;
            string appPath2;
            if (isAnalysisData)
            {
                appPath1 = String.Concat(appPath, @"\out1.csv"); //CSV file1 name here
                appPath2 = String.Concat(appPath, @"\out2.csv"); //CSV file2 name here
            }
            else
            {
                appPath1 = String.Concat(appPath, @"\out1future.csv"); //CSV file1 name here
                appPath2 = String.Concat(appPath, @"\out2future.csv"); //CSV file2 name here
            }
            string appPathFuture = String.Concat(appPath, @"\out1future.csv"); //CSV file1 name here
            List<string> fullDataSetList = loadCsvFile(appPath1);
            List<string> moreDataList = loadCsvFile(appPath2);
            List<string> radiationDataSet = loadCsvFile(appPathFuture);
            fullDataSetList.AddRange(moreDataList);
            List<double> allLons = new List<double>();
            List<double> allLats = new List<double>();
            
            //Generate 2 arrays for lats and lons
            foreach (string i in fullDataSetList)
            {
                string[] oneVariable = i.Split(',');
                double[] currentLonLat = new double[2]; //First index is longitude, second is latitude
                currentLonLat[0] = Convert.ToDouble(oneVariable[4]);
                currentLonLat[1] = Convert.ToDouble(oneVariable[5]);
                Boolean addCoordLon = true;
                Boolean addCoordLat = true;
                ///////Probably unnecessary <------------------ (Just use List.contains)///////////
                foreach (double j in allLons)
                {
                    if (j == currentLonLat[0])
                    {
                        addCoordLon = false;
                    }
                }
                if (addCoordLon)
                {
                    allLons.Add(currentLonLat[0]);
                }
                foreach (double j in allLats)
                {
                    if (j == currentLonLat[1])
                    {
                        addCoordLat = false;
                    }
                }
                if (addCoordLat)
                {
                    allLats.Add(currentLonLat[1]);
                }
                ///////////////////////////////////////////////////////////////////////////////////
            }

            //This determines the size of grid needed to accommodate the given latlons
            Boolean needSfcPressure = true;
            Boolean needSfcGeopotentialHeight = true;
            double sfcPressureCorn = -1;
            double sfcHGTCorn = -1;
            //This loop may need rethinking <------
            foreach (string i in fullDataSetList)
            {
                string[] oneVariable = i.Split(',');
                if (needSfcPressure && oneVariable[2] == "\"PRES\"" && oneVariable[3] == "\"surface\"")
                {
                    sfcPressureCorn = Convert.ToDouble(oneVariable[6]);
                    sfcPressureCorn = sfcPressureCorn / 100; //Convert to mbar from pa
                    needSfcPressure = false;
                }
                if (needSfcGeopotentialHeight && oneVariable[2] == "\"HGT\"" && oneVariable[3] == "\"surface\"")
                {
                    sfcHGTCorn = Convert.ToDouble(oneVariable[6]); //ASSUMING SAME COORDS AS ABOVE PRESSURE
                    needSfcGeopotentialHeight = false;
                }
            }
            double ABLHGTCorn = sfcHGTCorn + Zi;
            double differenceInHGTMag = 10000;
            String ABLPressureCornString = "-1";
            foreach (string i in fullDataSetList)
            {
                string[] oneVariable = i.Split(',');
                if (oneVariable[2] == "\"HGT\"" && oneVariable[3].Contains("mb"))
                {
                    double currentDifference = Math.Abs(Convert.ToDouble(oneVariable[6]) - ABLHGTCorn);
                    if (currentDifference <= differenceInHGTMag)
                    {
                        ABLPressureCornString = oneVariable[3];
                        differenceInHGTMag = currentDifference;
                    }
                }
            }
            String[] ABLPressureSplit = ABLPressureCornString.Split(' ');
            double ABLPressureCorn = Convert.ToDouble(ABLPressureSplit[0].TrimStart('\"'));
            double pressureCornDifference = sfcPressureCorn - ABLPressureCorn; //<----- Need for later
            List<double> allPossiblePressures = new List<double>();
            for (int i = 25; i < 1300; i = i + 25) //<--------- This System may need alot of work (1500 mbar may be out of range but unlikely)
            {
                allPossiblePressures.Add(i);
            }
            double differenceInSfcPressMag = 10000;
            double differenceInABLPressMag = 10000;
            double roundedSfcPressureCorn = -1;
            double roundedABLPressureCorn = -1;
            foreach (double i in allPossiblePressures)
            {
                double currentDifferenceSfc = Math.Abs(sfcPressureCorn - i);
                double currentDifferenceABL = Math.Abs(ABLPressureCorn - i);
                if (currentDifferenceSfc <= differenceInSfcPressMag)
                {
                    roundedSfcPressureCorn = i;
                    differenceInSfcPressMag = currentDifferenceSfc;
                }
                if (currentDifferenceABL <= differenceInABLPressMag)
                {
                    roundedABLPressureCorn = i;
                    differenceInABLPressMag = currentDifferenceABL;
                }
            }
            double zCountDouble = (roundedSfcPressureCorn - roundedABLPressureCorn) / 25;
            int zCount = Convert.ToInt32(zCountDouble);
            xLen = allLons.Count * resolutionMultiplierX - (resolutionMultiplierX - 1) + 2; //2 is included for ghost points
            yLen = allLats.Count * resolutionMultiplierY - (resolutionMultiplierY - 1) + 2; //2 is included for ghost points
            zLen = zCount * resolutionMultiplierZ - (resolutionMultiplierZ - 1) + 1; //1 is included for ghost points

            deltaX = new double[xLen, yLen, zLen];
            deltaY = new double[xLen, yLen, zLen];
            deltaZ = new double[xLen, yLen, zLen];
            latitudeData = new double[xLen,yLen,zLen];
            longitudeData = new double[xLen,yLen,zLen];
            pressureData = new double[xLen, yLen, zLen];
            geopotentialHeightData = new double[xLen,yLen,zLen];
            netRadiationData = new double[xLen,yLen,zLen];
            soilTemperatureData = new double[xLen,yLen,zLen];
            groundAirPressData = new double[xLen,yLen,zLen];
            relativeHumidityData = new double[xLen,yLen,zLen];
            verticalVelocityData = new double[xLen,yLen,zLen];
            windUData = new double[xLen,yLen,zLen];
            windVData = new double[xLen,yLen,zLen];
            windWData = new double[xLen,yLen,zLen];
            densityData = new double[xLen,yLen,zLen];
            tempData = new double[xLen,yLen,zLen];
            mixingRatioData = new double[xLen,yLen,zLen];
            for (int i = 0; i < xLen; i++)
            {
                for (int j = 0; j < yLen; j++)
                {
                    for (int k = 0; k < zLen; k++)
                    {
                        deltaX[i,j,k] = Double.NaN;
                        deltaY[i,j,k] = Double.NaN;
                        deltaZ[i,j,k] = Double.NaN;
                        latitudeData[i,j,k] = Double.NaN;
                        longitudeData[i,j,k] = Double.NaN;
                        pressureData[i,j,k] = Double.NaN;
                        geopotentialHeightData[i,j,k] = Double.NaN;
                        netRadiationData[i,j,k] = Double.NaN;
                        soilTemperatureData[i,j,k] = Double.NaN;
                        groundAirPressData[i,j,k] = Double.NaN;
                        relativeHumidityData[i,j,k] = Double.NaN;
                        verticalVelocityData[i,j,k] = Double.NaN;
                        windUData[i,j,k] = Double.NaN;
                        windVData[i,j,k] = Double.NaN;
                        windWData[i,j,k] = Double.NaN;
                        densityData[i,j,k] = Double.NaN;
                        tempData[i,j,k] = Double.NaN;
                        mixingRatioData[i,j,k] = Double.NaN;
                    }
                }
            }
            double currentLat = bottomLeftLat;
            double currentLon = bottomLeftLon;
            int lastCoordI = 1;
            int lastCoordJ = 1;
            //size of the gridpoints//
            for (int i = 0; i < xLen; i++)
            {
                for (int j = 0; j < yLen; j++)
                {
                    for (int k = 0; k < zLen; k++)
                    {
                        //These grids are totally complete
                        deltaX[i,j,k] = lengthQuarterDegreeLon / resolutionMultiplierX;
                        deltaY[i,j,k] = lengthQaurterDegreeLat / resolutionMultiplierY;
                        //TODO: Consider how this changes with time
                        deltaZ[i,j,k] = Zi / zLen;
                    }
                }
            }
            //Latlons//
            for (int i = 1; i < xLen; i = i + resolutionMultiplierX)
            {
                for (int j = 1; j < yLen; j = j + resolutionMultiplierY)
                {
                    for (int k = 0; k < zLen; k = k + resolutionMultiplierZ)
                    {
                        longitudeData[i, j, k] = currentLon;
                    }
                }
                currentLon = currentLon + 0.25;
            }
            for (int j = 1; j < yLen; j = j + resolutionMultiplierY)
            {
                for (int i = 1; i < xLen; i = i + resolutionMultiplierX)
                {
                    for (int k = 0; k < zLen; k = k + resolutionMultiplierZ)
                    {
                        latitudeData[i, j, k] = currentLat;
                    }
                }
                currentLat = currentLat + 0.25;
            }
            List<double> longitudeList = new List<double>(xLen);
            List<double> latitudeList = new List<double>(yLen);
            List<double>[,] pressureList = new List<double>[xLen,yLen];
            double[,] surfacePressureList = new double[xLen,yLen];
            for (int i = 0; i < xLen; i++)
            {
                longitudeList.Add(0);
            }
            for (int i = 0; i < yLen; i++)
            {
                latitudeList.Add(0);
            }
            for (int i = 0; i < xLen; i++)
            {
                for (int j = 0; j < yLen; j++)
                {
                    pressureList[i, j] = null;
                }
            }
            currentLat = bottomLeftLat;
            currentLon = bottomLeftLon;
            for (int i = 1; i < xLen; i = i + resolutionMultiplierX)
            {
            longitudeList[i] = currentLon;
            currentLon = currentLon + 0.25;
            }
            for (int i = 1; i < yLen; i = i + resolutionMultiplierY)
            {
            latitudeList[i] = currentLat;
            currentLat = currentLat + 0.25;
            }
            foreach (string i in fullDataSetList)
            {
                double roundedPressure = -1;
                double differenceInPressMag = 10000;
                string[] oneVariable = i.Split(',');
                //Press//
                if (oneVariable[2] == "\"PRES\"" && oneVariable[3] == "\"surface\"")
                {
                    foreach (double j in allPossiblePressures)
                    {
                        double currentDifference = Math.Abs(Convert.ToDouble(oneVariable[6])/100 - j); //Converted to mbar
                        if (currentDifference <= differenceInPressMag)
                        {
                            roundedPressure = j;
                            differenceInPressMag = currentDifference;
                        }
                    }
                    List<double> onePressureList = new List<double>(zLen);
                    for (int j = 0; j < zLen; j++)
                    {
                        onePressureList.Add(0);
                    }
                    surfacePressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                        latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))] = roundedPressure;
                    roundedPressure = roundedPressure - 25;
                    for (int j = resolutionMultiplierZ; j < zLen; j = j + resolutionMultiplierZ)
                    {
                        pressureData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), j] =
                            roundedPressure * 100; //Converted to pa
                        onePressureList[j] = roundedPressure; //Kept as mbar
                        roundedPressure = roundedPressure - 25;
                    }
                    pressureData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                           latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), 0] =
                        Convert.ToDouble(oneVariable[6]);
                    pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                        latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))] = onePressureList;
                    //Ground Pressure//
                    for (int j = 0; j < zLen; j++)
                    {
                        groundAirPressData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), j] =
                            Convert.ToDouble(oneVariable[6]);
                    }
                }
            }
            
            foreach (string i in fullDataSetList)
            {
                string[] oneVariable = i.Split(',');
                //HGT//
                if (oneVariable[2] == "\"HGT\"")
                {
                    if (oneVariable[3] == "\"surface\"")
                    {
                        geopotentialHeightData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), 0] =
                            Convert.ToDouble(oneVariable[6]);
                    }
                    else if (oneVariable[3].Contains("mb") && !oneVariable[3].Contains("above"))
                    {
                        String[] currentPressureString = oneVariable[3].Split(' ');
                        double currentPressure = Convert.ToDouble(currentPressureString[0].TrimStart('\"'));
                        if (pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                            latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].Contains(currentPressure))
                        {
                            int pressureIndex = pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].IndexOf(currentPressure);
                            geopotentialHeightData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                    latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), pressureIndex] =
                                Convert.ToDouble(oneVariable[6]);
                        } 
                    }
                }
                //TSOIL//
                else if (oneVariable[2] == "\"TSOIL\"" && oneVariable[3] == "\"0-0.1 m below ground\"")
                {
                    for (int j = 0; j < zLen; j++)
                    {
                        soilTemperatureData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), j] =
                            Convert.ToDouble(oneVariable[6]);
                    }
                }
                //RH//
                else if (oneVariable[2] == "\"RH\"")
                {
                    if (oneVariable[3] == "\"2 m above ground\"")
                    {
                        relativeHumidityData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), 0] =
                            Convert.ToDouble(oneVariable[6]);
                    }
                    else if (oneVariable[3].Contains("mb") && !oneVariable[3].Contains("above"))
                    {
                        String[] currentPressureString = oneVariable[3].Split(' ');
                        double currentPressure = Convert.ToDouble(currentPressureString[0].TrimStart('\"'));
                        if (pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                            latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].Contains(currentPressure))
                        {
                            int pressureIndex = pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].IndexOf(currentPressure);
                            relativeHumidityData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                    latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), pressureIndex] =
                                Convert.ToDouble(oneVariable[6]);
                        } 
                    }
                }
                else if (oneVariable[2] == "\"VVEL\"")
                {
                    if (oneVariable[3].Contains("mb") && !oneVariable[3].Contains("above"))
                    {
                        String[] currentPressureString = oneVariable[3].Split(' ');
                        double currentPressure = Convert.ToDouble(currentPressureString[0].TrimStart('\"'));
                        if (pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                            latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].Contains(currentPressure))
                        {
                            int pressureIndex = pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].IndexOf(currentPressure);
                            verticalVelocityData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                    latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), pressureIndex] =
                                Convert.ToDouble(oneVariable[6]);
                        } 
                        else if (currentPressure == surfacePressureList[
                                     longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                     latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))])
                        {
                            verticalVelocityData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                    latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), 0] =
                                Convert.ToDouble(oneVariable[6]);
                        }
                    }
                }
                //TEMP//
                else if (oneVariable[2] == "\"TMP\"")
                {
                    if (oneVariable[3] == "\"surface\"")
                    {
                        tempData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), 0] =
                            Convert.ToDouble(oneVariable[6]);
                    }
                    else if (oneVariable[3].Contains("mb") && !oneVariable[3].Contains("above"))
                    {
                        String[] currentPressureString = oneVariable[3].Split(' ');
                        double currentPressure = Convert.ToDouble(currentPressureString[0].TrimStart('\"'));
                        if (pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                            latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].Contains(currentPressure))
                        {
                            int pressureIndex = pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].IndexOf(currentPressure);
                            tempData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                    latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), pressureIndex] =
                                Convert.ToDouble(oneVariable[6]);
                        } 
                    }
                }
                //Mixing Ratio//
                else if (oneVariable[2] == "\"SPFH\"")
                {
                    if (oneVariable[3] == "\"2 m above ground\"")
                    {
                        mixingRatioData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), 0] =
                            Math.Pow(1 / (Convert.ToDouble(oneVariable[6])) - 1,-1);
                    }
                    else if (oneVariable[3].Contains("mb") && !oneVariable[3].Contains("above"))
                    {
                        String[] currentPressureString = oneVariable[3].Split(' ');
                        double currentPressure = Convert.ToDouble(currentPressureString[0].TrimStart('\"'));
                        if (pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                            latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].Contains(currentPressure))
                        {
                            int pressureIndex = pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].IndexOf(currentPressure);
                            mixingRatioData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                    latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), pressureIndex] =
                                Math.Pow(1 / (Convert.ToDouble(oneVariable[6])) - 1,-1);
                        } 
                    }
                }
                //U Wind//
                else if (oneVariable[2] == "\"UGRD\"")
                {
                    if (oneVariable[3] == "\"10 m above ground\"")
                    {
                        windUData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), 0] =
                            Convert.ToDouble(oneVariable[6]);
                    }
                    else if (oneVariable[3].Contains("mb") && !oneVariable[3].Contains("above"))
                    {
                        String[] currentPressureString = oneVariable[3].Split(' ');
                        double currentPressure = Convert.ToDouble(currentPressureString[0].TrimStart('\"'));
                        if (pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                            latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].Contains(currentPressure))
                        {
                            int pressureIndex = pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].IndexOf(currentPressure);
                            windUData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                    latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), pressureIndex] =
                                Convert.ToDouble(oneVariable[6]);
                        } 
                    }
                }
                //V Wind//
                else if (oneVariable[2] == "\"VGRD\"")
                {
                    if (oneVariable[3] == "\"10 m above ground\"")
                    {
                        windVData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), 0] =
                            Convert.ToDouble(oneVariable[6]);
                    }
                    else if (oneVariable[3].Contains("mb") && !oneVariable[3].Contains("above"))
                    {
                        String[] currentPressureString = oneVariable[3].Split(' ');
                        double currentPressure = Convert.ToDouble(currentPressureString[0].TrimStart('\"'));
                        if (pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                            latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].Contains(currentPressure))
                        {
                            int pressureIndex = pressureList[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))].IndexOf(currentPressure);
                            windVData[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                                    latitudeList.IndexOf(Convert.ToDouble(oneVariable[5])), pressureIndex] =
                                Convert.ToDouble(oneVariable[6]);
                        } 
                    }
                }
            }
            //Net Radiation//
            //NOTE: This is done with an assumption about how the radiation changes with altitude. Revision may be required.//
            double[,] downwardShortWaveFlux = new double[xLen, yLen];
            double[,] downwardLongWaveFlux = new double[xLen, yLen];
            double[,] upwardShortWaveFlux = new double[xLen, yLen];
            double[,] upwardLongWaveFlux = new double[xLen, yLen];
            foreach (string i in radiationDataSet)
            {
                string[] oneVariable = i.Split(',');
                if (oneVariable[2] == "\"DSWRF\"" && oneVariable[3] == "\"surface\"")
                {
                    downwardShortWaveFlux[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                            latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))] =
                        Convert.ToDouble(oneVariable[6]);
                }
                else if (oneVariable[2] == "\"DLWRF\"" && oneVariable[3] == "\"surface\"")
                {
                    downwardLongWaveFlux[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                            latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))] =
                        Convert.ToDouble(oneVariable[6]);
                }
                else if (oneVariable[2] == "\"USWRF\"" && oneVariable[3] == "\"surface\"")
                {
                    upwardShortWaveFlux[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                            latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))] =
                        Convert.ToDouble(oneVariable[6]);
                }
                else if (oneVariable[2] == "\"ULWRF\"" && oneVariable[3] == "\"surface\"")
                {
                    upwardLongWaveFlux[longitudeList.IndexOf(Convert.ToDouble(oneVariable[4])),
                            latitudeList.IndexOf(Convert.ToDouble(oneVariable[5]))] =
                        Convert.ToDouble(oneVariable[6]);
                }
            }
            for (int i = 1; i < xLen; i = i + resolutionMultiplierX)
            {
                for (int j = 1; j < yLen; j = j + resolutionMultiplierY)
                {
                    for (int k = 0; k < zLen; k++)
                    {
                        double downShort = downwardShortWaveFlux[i, j];
                        double downLong = downwardLongWaveFlux[i, j];
                        double upShort = upwardShortWaveFlux[i, j];
                        double upLong = upwardLongWaveFlux[i, j];
                        netRadiationData[i, j, k] = (downShort + downLong - upLong - upShort) + 2 * k; //This "2" is what I estimate. May change.
                    }
                }
            }
            //W Wind//
            for (int i = 1; i < xLen; i = i + resolutionMultiplierX)
            {
                for (int j = 1; j < yLen; j = j + resolutionMultiplierY)
                {
                    for (int k = 0; k < zLen; k = k +resolutionMultiplierZ)
                    {
                        double conversionRatio = deltaZ[i,j,k] / 2500; //<----------- MAJOR CHANGE///////
                        windWData[i,j,k] = -verticalVelocityData[i,j,k] * conversionRatio;
                    }
                }
            }
            //Air Density//
            for (int i = 1; i < xLen; i = i + resolutionMultiplierX)
            {
                for (int j = 1; j < yLen; j = j + resolutionMultiplierY)
                {
                    for (int k = 0; k < zLen; k = k + resolutionMultiplierZ)
                    {
                        double virtualTemp = tempData[i,j,k] * (1 + (a * mixingRatioData[i,j,k])); //***1.21***
                        densityData[i, j, k] = pressureData[i, j, k] / (Rd * virtualTemp); //***1.23***
                    }
                }
            }
        }

        public double[,,] getTempData()
        {
            return tempData;
        }
        public double[,,] getPressureData()
        {
            return pressureData;
        }
        public double[,,] getMixingRatioData()
        {
            return mixingRatioData;
        }
        public double[,,] getDensityData()
        {
            return densityData;
        }
        public double[,,] getWindUData()
        {
            return windUData;
        }
        public double[,,] getWindVData()
        {
            return windVData;
        }
        public double[,,] getWindWData()
        {
            return windWData;
        }
        public double[,,] getDeltaXData()
        {
            return deltaX;
        }
        public double[,,] getDeltaYData()
        {
            return deltaY;
        }
        public double[,,] getDeltaZData()
        {
            return deltaZ;
        }
        public double[,,] getLatitudeData()
        {
            return latitudeData;
        }
        public double[,,] getLongitudeData()
        {
            return longitudeData;
        }
        public double[,,] getGeopotentialHeightData()
        {
            return geopotentialHeightData;
        }
        public double[,,] getNetRadiationData()
        {
            return netRadiationData;
        }
        public double[,,] getSoilTemperatureData()
        {
            return soilTemperatureData;
        }
        public double[,,] getGroundPressureData()
        {
            return groundAirPressData;
        }
        public double[,,] getRelativeHumidityData()
        {
            return relativeHumidityData;
        }
        public double[,,] getVerticalVelocityData()
        {
            return verticalVelocityData;
        }
        public int getXLen()
        {
            return xLen;
        }
        public int getYLen()
        {
            return yLen;
        }
        public int getZLen()
        {
            return zLen;
        }
    }
    
    public class Calculator
    {
        public const double g = 9.81; //gravity
        public const double R = 8.3144598; //gas constant <NEEDS A REVIEW>
        public const double Rd = 287.053; //dry-air gas constant
        public const double Cpd = 1004; //dry-air specific heat
        public const double bH = 5e-4; //convective transport coefficient
        public const double RdCp = 0.28571; //Dimensionless
        public const double a = 0.61; //g/g
        public const double CH = 1.1e-2; //Assuming an average number here for a mostly smooth terrain. <-----
        public const double CD = 2e-3; //Assuming a smooth surface.
        public const double twoOmega = 1.458423e-4; //s^-1
        private double Zi; //ABL depth (MUST ALIGN WITH deltaZ and Zlen!!!!)
        private double wB;
        private Grid presentGrid;
        private Grid pastGrid;
        private Grid futureGrid;
        private GridPoint targetPoint;
        private GridPoint pastTargetPoint;
        private GridPoint topPoint;
        private GridPoint infrontPoint;
        private GridPoint behindPoint;
        private GridPoint leftPoint;
        private GridPoint rightPoint;
        private GridPoint bottomPoint;
        private GridPoint futureTargetPoint;
        
        
        public Calculator(double Zi)
        {
            this.Zi = Zi;
        }

        public void update(Grid presentGrid, Grid pastGrid, Grid futureGrid)
        {
            this.presentGrid = presentGrid; //full present grid
            this.pastGrid = pastGrid;
            this.futureGrid = futureGrid;
            targetPoint = this.presentGrid.getTargetPoint();
            pastTargetPoint = this.pastGrid.getTargetPoint();
            topPoint = this.presentGrid.getTopPoint();
            infrontPoint = this.presentGrid.getInfrontPoint();
            behindPoint = this.presentGrid.getBehindPoint();
            leftPoint = this.presentGrid.getLeftPoint();
            rightPoint = this.presentGrid.getRightPoint();
            bottomPoint = this.presentGrid.getBottomPoint();
            futureTargetPoint = this.futureGrid.getTargetPoint();
        }

        public void temp()
        {
            targetPoint.setTurbTemp(0);
            targetPoint.setTurbUWind(0);
            targetPoint.setTurbVWind(0);
            targetPoint.setTurbMixingRatio(0);
            
        }
        
        //This method may not sample air from high enough in atmosphere. Revision may be required. <-------- BIG TIME
        //For revisions refer to RT turb for 3.38
        //Use before calculations
        public void initTurbTemp()
        {
            if ((targetPoint.isOnBottomBoundary()) && (targetPoint.getWindSpeed() > 1.543)) //threshold for calm winds
            {
                double FH = CH * targetPoint.getWindSpeed() * (targetPoint.getSfcTemp() - targetPoint.getTemp()); //***3.35***
                double turbTemp = -1.2 * FH/Zi; //***3.41***
                turbTemp = turbTemp / presentGrid.getZLen(); //averaged over each grid point
                targetPoint.setTurbTemp(turbTemp);
            }
            else if ((targetPoint.isOnBottomBoundary()) && (targetPoint.getWindSpeed() <= 1.543))
            {
                GridPoint MLPoint = presentGrid.getPoint(targetPoint.getICoord(), targetPoint.getJCoord(), 5); ///<----------SET ALL TOPPOINTS TO THIS
                double potentialTempAir = MLPoint.getTemp() *
                                          Math.Pow((targetPoint.getNearSfcPressure() / MLPoint.getPressure()), RdCp); //***3.12***
                double potentialTempSfc = targetPoint.getSfcTemp() *
                                          Math.Pow((targetPoint.getNearSfcPressure() / targetPoint.getPressure()),
                                              RdCp); //***3.12*** <------------------------------------------------------- Unused because of below
                double potentialTempSfcBot = targetPoint.getTemp() *
                                             Math.Pow((targetPoint.getNearSfcPressure() / targetPoint.getPressure()), RdCp); //***3.12*** (TEMP)
                double virtualPotentialTempAir = potentialTempAir * (1 + (a * MLPoint.getMixingRatio())); //***3.13***
                double virtualTemperatureAir = MLPoint.getTemp() * (1 + (a * MLPoint.getMixingRatio())); //***1.21***
                double virtualPotentialTemperatureSurface =
                    potentialTempSfcBot * (1 + (a * targetPoint.getMixingRatio())); //***3.13*** <------------ CHANGED TO AIR TEMP FOR TESTING INSTEAD OF SOIL TEMP
                wB = Math.Pow(
                    g * Zi * (virtualPotentialTemperatureSurface - virtualPotentialTempAir) / virtualTemperatureAir,
                    0.5); //***3.38***
                double FH = bH * wB * (targetPoint.getSfcTemp() - potentialTempAir); //***3.36***
                double turbTemp = -1.2 * FH / Zi; //***3.41***
                turbTemp = turbTemp / presentGrid.getZLen(); //averaged over each grid point
                //////////////////////THIS SYSTEM NEEDS WORK///////////////////////////////////////
                if (!double.IsNaN(turbTemp))
                {
                    targetPoint.setTurbTemp(turbTemp);
                }
                else
                {
                    targetPoint.setTurbTemp(0);
                    wB = 0;
                }
                //////////////////////////////////////////////////////////////////////////////////////
            }
            else if (!targetPoint.isOnBottomBoundary())
            {
                targetPoint.setTurbTemp(bottomPoint.getTurbTemp());
            }
            else
            {
                Console.WriteLine("Major Error Encountered in temp turbulence initialization!");
            }
        }

        //Use before calculations
        //Use print to verify wB when testing. Should be between 10 to 50 m/s.
        //May need custom input for ABL height later.
        public void initTurbMixingRatioandUVWind()
        {
            GridPoint ABLPointTop = presentGrid.getPoint(targetPoint.getICoord(), targetPoint.getJCoord(), presentGrid.getZLen()-1);
            GridPoint ABLPointBot = presentGrid.getPoint(targetPoint.getICoord(), targetPoint.getJCoord(), presentGrid.getZLen()-3);
            GridPoint sfcPointBot = presentGrid.getPoint(targetPoint.getICoord(), targetPoint.getJCoord(), 0); //MAYBE CHANGE?
            GridPoint sfcPointTop = presentGrid.getPoint(targetPoint.getICoord(), targetPoint.getJCoord(), 2); //MAYBE CHANGE?
            double potentialTempSfcTop = sfcPointTop.getTemp() *
                                      Math.Pow((sfcPointTop.getNearSfcPressure() / sfcPointTop.getPressure()), RdCp); //***3.12***
            double potentialTempSfcBot = sfcPointBot.getTemp() *
                                         Math.Pow((sfcPointBot.getNearSfcPressure() / sfcPointBot.getPressure()), RdCp); //***3.12***
            double potentialTempABLBot = ABLPointBot.getTemp() *
                                         Math.Pow((ABLPointBot.getNearSfcPressure() / ABLPointBot.getPressure()), RdCp); //***3.12***
            double potentialTempABLTop = ABLPointTop.getTemp() *
                                         Math.Pow((ABLPointTop.getNearSfcPressure() / ABLPointTop.getPressure()), RdCp); //***3.12***
            double potentialTempSfcSoil = sfcPointBot.getSfcTemp() * Math.Pow((sfcPointBot.getNearSfcPressure() / sfcPointBot.getPressure()),
                                          RdCp); //***3.12*** <----------------------------------------------- Unused because of below
            double deltaPotentialABL = potentialTempABLTop - potentialTempABLBot;
            double deltaPotentialSfc = potentialTempSfcTop - potentialTempSfcBot;
            double deltaRatioABL = ABLPointTop.getMixingRatio() - ABLPointBot.getMixingRatio();
            double deltaRatioSfc = sfcPointTop.getMixingRatio() - sfcPointBot.getMixingRatio();
            double virtualPotentialTempABL = potentialTempABLBot * (1 + (a * ABLPointBot.getMixingRatio())); //***3.13***
            double virtualPotentialTempSfc = potentialTempSfcBot * (1 + (a * sfcPointBot.getMixingRatio())); //***3.13*** <------------ CHANGED TO AIR TEMP FOR TESTING INSTEAD OF SOIL TEMP
            double virtualTempABL = ABLPointBot.getTemp() * (1 + (a * ABLPointBot.getMixingRatio())); //***1.21***
            wB = Math.Pow(
                (g * Zi * (virtualPotentialTempSfc - virtualPotentialTempABL) / virtualTempABL), 0.5); //***3.38***
            double turbMixingRatio = -bH * wB * (0.2 * deltaRatioABL * Math.Abs(deltaPotentialSfc / deltaPotentialABL) +
                                                 Math.Abs(deltaRatioSfc)) / Zi; //***4.56***
            double wT = CD * targetPoint.getWindSpeed(); //***10.21*** (Assuming neutral conditions)
            double turbUWind = -wT * targetPoint.getWindU() / Zi; //***10.19a***
            double turbVWind = -wT * targetPoint.getWindV() / Zi; //***10.19b***
            /////////////////////////////////THIS SYSTEM NEEDS WORK///////////////////////////////
            if (!double.IsNaN(turbMixingRatio)) 
            {
                targetPoint.setTurbMixingRatio(turbMixingRatio);
            }
            else
            {
                targetPoint.setTurbMixingRatio(0);
                wB = 0;
            }
            //////////////////////////////////////////////////////////////////////////////////////
            targetPoint.setTurbUWind(turbUWind);
            targetPoint.setTurbVWind(turbVWind);
        }
        
        public void calcTemp(int deltaTime)
        {
            double Cp = Cpd * (1 + 1.84 * targetPoint.getMixingRatio()); //***3.3***
            double output;
            if (!targetPoint.isOnBottomBoundary()) {
                output = pastTargetPoint.getTemp() + 2 * deltaTime * (
                                    -targetPoint.getWindU() * (infrontPoint.getTemp() - behindPoint.getTemp()) /
                                    (2 * targetPoint.getDeltaX())
                                    - targetPoint.getWindV() * (leftPoint.getTemp() - rightPoint.getTemp()) /
                                    (2 * targetPoint.getDeltaY())
                                    - targetPoint.getWindW() * (topPoint.getTemp() - bottomPoint.getTemp()) /
                                    (2 * targetPoint.getDeltaZ())
                                    - R * targetPoint.getTemp() / (targetPoint.getPressure() * Cp) *
                                    (topPoint.getNetRadiation() - bottomPoint.getNetRadiation()) /
                                    (2 * targetPoint.getDeltaZ())
                                    - targetPoint.getTurbTemp()); //***20.17*** (No Rain)
            }
            else
            {
                output = pastTargetPoint.getTemp() + 2 * deltaTime * (
                                    -targetPoint.getWindU() * (infrontPoint.getTemp() - behindPoint.getTemp()) /
                                    (2 * targetPoint.getDeltaX())
                                    - targetPoint.getWindV() * (leftPoint.getTemp() - rightPoint.getTemp()) /
                                    (2 * targetPoint.getDeltaY())
                                    - targetPoint.getWindW() * (topPoint.getTemp() - targetPoint.getTemp()) /
                                    (targetPoint.getDeltaZ())
                                    - R * targetPoint.getTemp() / (targetPoint.getPressure() * Cp) *
                                    (topPoint.getNetRadiation() - targetPoint.getNetRadiation()) /
                                    (targetPoint.getDeltaZ())
                                    - targetPoint.getTurbTemp()); //***20.17*** (No Rain)
            }
            futureTargetPoint.setTemp(output);
        }

        public void calcMixingRatio(int deltaTime)
        {
            double output;
            if (!targetPoint.isOnBottomBoundary())
            {
                output = pastTargetPoint.getMixingRatio() + 2 * deltaTime * (
                                    -targetPoint.getWindU() *
                                    (infrontPoint.getMixingRatio() - behindPoint.getMixingRatio()) /
                                    (2 * targetPoint.getDeltaX())
                                    - targetPoint.getWindV() *
                                    (leftPoint.getMixingRatio() - rightPoint.getMixingRatio()) /
                                    (2 * targetPoint.getDeltaY())
                                    - targetPoint.getWindW() *
                                    (topPoint.getMixingRatio() - bottomPoint.getMixingRatio()) /
                                    (2 * targetPoint.getDeltaZ()) -
                                    targetPoint.getTurbMixingRatio()); //***20.5*** (No rain, neutral)
            }
            else
            {
                output = pastTargetPoint.getMixingRatio() + 2 * deltaTime * (
                                    -targetPoint.getWindU() *
                                    (infrontPoint.getMixingRatio() - behindPoint.getMixingRatio()) /
                                    (2 * targetPoint.getDeltaX())
                                    - targetPoint.getWindV() *
                                    (leftPoint.getMixingRatio() - rightPoint.getMixingRatio()) /
                                    (2 * targetPoint.getDeltaY())
                                    - targetPoint.getWindW() *
                                    (topPoint.getMixingRatio() - targetPoint.getMixingRatio()) /
                                    (targetPoint.getDeltaZ()) -
                                    targetPoint.getTurbMixingRatio()); //***20.5*** (No rain, neutral)
            }
            futureTargetPoint.setMixingRatio(output);
            
        }

        public void calcWindU(int deltaTime)
        {
            double output;
            if (!targetPoint.isOnBottomBoundary())
            {
                output = pastTargetPoint.getWindU() + 2 * deltaTime * (
                                    -targetPoint.getWindU() * (infrontPoint.getWindU() - behindPoint.getWindU()) /
                                    (2 * targetPoint.getDeltaX())
                                    - targetPoint.getWindV() * (leftPoint.getWindU() - rightPoint.getWindU()) /
                                    (2 * targetPoint.getDeltaY())
                                    - targetPoint.getWindW() * (topPoint.getWindU() - bottomPoint.getWindU()) /
                                    (2 * targetPoint.getDeltaZ()) -
                                    (infrontPoint.getPressure() - behindPoint.getPressure()) /
                                    (2 * targetPoint.getDeltaX() * targetPoint.getAirDensity()) +
                                    twoOmega * Math.Sin(targetPoint.getLatitude() * Math.PI / 180) *
                                    targetPoint.getWindV() -
                                    targetPoint.getTurbUWind()); //***20.1 and 10.6***
            }
            else
            {
                output = pastTargetPoint.getWindU() + 2 * deltaTime * (
                                    -targetPoint.getWindU() * (infrontPoint.getWindU() - behindPoint.getWindU()) /
                                    (2 * targetPoint.getDeltaX())
                                    - targetPoint.getWindV() * (leftPoint.getWindU() - rightPoint.getWindU()) /
                                    (2 * targetPoint.getDeltaY())
                                    - targetPoint.getWindW() * (topPoint.getWindU() - targetPoint.getWindU()) /
                                    (targetPoint.getDeltaZ()) -
                                    (infrontPoint.getPressure() - behindPoint.getPressure()) /
                                    (2 * targetPoint.getDeltaX() * targetPoint.getAirDensity()) +
                                    twoOmega * Math.Sin(targetPoint.getLatitude() * Math.PI / 180) *
                                    targetPoint.getWindV() -
                                    targetPoint.getTurbUWind()); //***20.1 and 10.6***
            }

            futureTargetPoint.setWindU(output);
        }
        
        public void calcWindV(int deltaTime)
        {
            double output;
            if (!targetPoint.isOnBottomBoundary())
            {
                output = pastTargetPoint.getWindV() + 2 * deltaTime * (
                                    -targetPoint.getWindU() * (infrontPoint.getWindV() - behindPoint.getWindV()) /
                                    (2 * targetPoint.getDeltaX())
                                    - targetPoint.getWindV() * (leftPoint.getWindV() - rightPoint.getWindV()) /
                                    (2 * targetPoint.getDeltaY())
                                    - targetPoint.getWindW() * (topPoint.getWindV() - bottomPoint.getWindV()) /
                                    (2 * targetPoint.getDeltaZ()) -
                                    (leftPoint.getPressure() - rightPoint.getPressure()) /
                                    (2 * targetPoint.getDeltaY() * targetPoint.getAirDensity()) -
                                    twoOmega * Math.Sin(targetPoint.getLatitude() * Math.PI / 180) *
                                    targetPoint.getWindU() -
                                    targetPoint.getTurbVWind()); //***20.2 and 10.6***
            }
            else
            {
                output = pastTargetPoint.getWindV() + 2 * deltaTime * (
                             -targetPoint.getWindU() * (infrontPoint.getWindV() - behindPoint.getWindV()) /
                             (2 * targetPoint.getDeltaX())
                             - targetPoint.getWindV() * (leftPoint.getWindV() - rightPoint.getWindV()) /
                             (2 * targetPoint.getDeltaY())
                             - targetPoint.getWindW() * (topPoint.getWindV() - targetPoint.getWindV()) /
                             (targetPoint.getDeltaZ()) -
                             (leftPoint.getPressure() - rightPoint.getPressure()) /
                             (2 * targetPoint.getDeltaY() * targetPoint.getAirDensity()) -
                             twoOmega * Math.Sin(targetPoint.getLatitude() * Math.PI / 180) *
                             targetPoint.getWindU() -
                             targetPoint.getTurbVWind()); //***20.2 and 10.6***
            }

            futureTargetPoint.setWindV(output);
        }
        
        public void calcWindW(int deltaTime)
        {
            double output;
            if (!targetPoint.isOnBottomBoundary())
            {
                double tempBarTarget = 288.15 - 6.5 * targetPoint.getGeopotentialHeight() / 1000; //Convert to gpkm
                double pressureBarTarget = 101325 * Math.Pow((288.15 / tempBarTarget), -5.255877);
                double virtualTempTarget =
                    targetPoint.getTemp() * (1 + (a * targetPoint.getMixingRatio())); //***1.21***
                double densityBarTarget = pressureBarTarget / (Rd * virtualTempTarget); //***1.23***
                double densityPrimeTarget = targetPoint.getAirDensity() - densityBarTarget;
                double tempBarTop = 288.15 - 6.5 * topPoint.getGeopotentialHeight() / 1000; //Convert to gpkm
                double pressureBarTop = 101325 * Math.Pow((288.15 / tempBarTop), -5.255877);
                double pressurePrimeTop = topPoint.getPressure() - pressureBarTop;
                double tempBarBot = 288.15 - 6.5 * bottomPoint.getGeopotentialHeight() / 1000; //Convert to gpkm
                double pressureBarBot = 101325 * Math.Pow((288.15 / tempBarBot), -5.255877);
                double pressurePrimeBot = bottomPoint.getPressure() - pressureBarBot;
                //Console.WriteLine(pressureBarTarget);
                //Console.WriteLine(wB);
                output = pastTargetPoint.getWindW() + 2 * deltaTime * (
                                    -targetPoint.getWindU() * (infrontPoint.getWindW() - behindPoint.getWindW()) /
                                    (2 * targetPoint.getDeltaX())
                                    - targetPoint.getWindV() * (leftPoint.getWindW() - rightPoint.getWindW()) /
                                    (2 * targetPoint.getDeltaY())
                                    - targetPoint.getWindW() * (topPoint.getWindW() - bottomPoint.getWindW()) /
                                    (2 * targetPoint.getDeltaZ()) -
                                    (pressurePrimeTop - pressurePrimeBot) /
                                    (2 * targetPoint.getDeltaZ() * densityBarTarget) -
                                    densityPrimeTarget * g /
                                    densityBarTarget - 0.08* wB); //***20.3*** + ***19.22*** (sort of)
            }
            else
            {
                double tempBarTarget = 288.15 - 6.5 * targetPoint.getGeopotentialHeight() / 1000; //Convert to gpkm
                double pressureBarTarget = 101325 * Math.Pow((288.15 / tempBarTarget), -5.255877);
                double virtualTempTarget =
                    targetPoint.getTemp() * (1 + (a * targetPoint.getMixingRatio())); //***1.21***
                double densityBarTarget = pressureBarTarget / (Rd * virtualTempTarget); //***1.23***
                double densityPrimeTarget = targetPoint.getAirDensity() - densityBarTarget;
                double tempBarTop = 288.15 - 6.5 * topPoint.getGeopotentialHeight() / 1000; //Convert to gpkm
                double pressureBarTop = 101325 * Math.Pow((288.15 / tempBarTop), -5.255877);
                double pressurePrimeTop = topPoint.getPressure() - pressureBarTop;
                double pressurePrimeTarget = targetPoint.getPressure() - pressureBarTarget;
                output = pastTargetPoint.getWindW() + 2 * deltaTime * (
                             -targetPoint.getWindU() * (infrontPoint.getWindW() - behindPoint.getWindW()) /
                             (2 * targetPoint.getDeltaX())
                             - targetPoint.getWindV() * (leftPoint.getWindW() - rightPoint.getWindW()) /
                             (2 * targetPoint.getDeltaY())
                             - targetPoint.getWindW() * (topPoint.getWindW() - targetPoint.getWindW()) /
                             (targetPoint.getDeltaZ()) -
                             (pressurePrimeTop - pressurePrimeTarget) /
                             (targetPoint.getDeltaZ() * densityBarTarget) -
                             densityPrimeTarget * g /
                             densityBarTarget - 0.08 * wB); //***20.3*** + ***19.22*** (sort of)
            }
            futureTargetPoint.setWindW(output);
        }

        public void calcAirDensity(int deltaTime)
        {
            double output;
            if (!targetPoint.isOnBottomBoundary())
            {
                output = pastTargetPoint.getAirDensity() + 2 * deltaTime * (
                                    -targetPoint.getWindU() *
                                    (infrontPoint.getAirDensity() - behindPoint.getAirDensity()) /
                                    (2 * targetPoint.getDeltaX())
                                    - targetPoint.getWindV() *
                                    (leftPoint.getAirDensity() - rightPoint.getAirDensity()) /
                                    (2 * targetPoint.getDeltaY())
                                    - targetPoint.getWindW() *
                                    (topPoint.getAirDensity() - bottomPoint.getAirDensity()) /
                                    (2 * targetPoint.getDeltaZ()) - targetPoint.getAirDensity() *
                                    ((infrontPoint.getWindU() - behindPoint.getWindU()) /
                                     (2 * targetPoint.getDeltaX()) +
                                     (leftPoint.getWindV() - rightPoint.getWindV()) / (2 * targetPoint.getDeltaY()) +
                                     (topPoint.getWindW() - bottomPoint.getWindW()) /
                                     (2 * targetPoint.getDeltaZ()))); //***20.6***
            }
            else
            {
                output = pastTargetPoint.getAirDensity() + 2 * deltaTime * (
                             -targetPoint.getWindU() *
                             (infrontPoint.getAirDensity() - behindPoint.getAirDensity()) /
                             (2 * targetPoint.getDeltaX())
                             - targetPoint.getWindV() *
                             (leftPoint.getAirDensity() - rightPoint.getAirDensity()) /
                             (2 * targetPoint.getDeltaY())
                             - targetPoint.getWindW() *
                             (topPoint.getAirDensity() - targetPoint.getAirDensity()) /
                             (targetPoint.getDeltaZ()) - targetPoint.getAirDensity() *
                             ((infrontPoint.getWindU() - behindPoint.getWindU()) /
                              (2 * targetPoint.getDeltaX()) +
                              (leftPoint.getWindV() - rightPoint.getWindV()) / (2 * targetPoint.getDeltaY()) +
                              (topPoint.getWindW() - targetPoint.getWindW()) /
                              (targetPoint.getDeltaZ()))); //***20.6***
            }
            futureTargetPoint.setAirDensity(output);
            
        }

        //To be executed last in loop
        public void calcAirPressure(int deltaTime)
        {
            double futureVirtualTemp = futureTargetPoint.getTemp() * (1 + (a * futureTargetPoint.getMixingRatio())); //***1.21***
            double output = futureTargetPoint.getAirDensity() * Rd * futureVirtualTemp; //***1.23***
            futureTargetPoint.setAirPressure(output);
        }

        public void calcRest()
        {
            futureTargetPoint.setDeltaXYZ(targetPoint.getDeltaX(),targetPoint.getDeltaY(),targetPoint.getDeltaZ());
            futureTargetPoint.setLatLon(targetPoint.getLatitude(),targetPoint.getLongitude());
            futureTargetPoint.setGeopotentialHeight(targetPoint.getGeopotentialHeight());
            futureTargetPoint.setNetRadiation(targetPoint.getNetRadiation());
            futureTargetPoint.setSoilTemperature(targetPoint.getSfcTemp());
            futureTargetPoint.setNearSfcPressure(targetPoint.getNearSfcPressure());
            futureTargetPoint.setRelativeHumidity(targetPoint.getRelativeHumidity());
            futureTargetPoint.setVerticalVelocity(targetPoint.getVerticalVelocity());
            futureTargetPoint.setWindW(targetPoint.getWindW());
        }

        public Grid getPresentGrid()
        {
            return presentGrid;
        }

        public Grid getFutureGrid()
        {
            return futureGrid;
        }
    }
    
    //TODO: Restructure to remove 3D grids that are not gridpoints
    public class Grid
    {
        public const int timeBetweenGivenData = 3600; //PASS THIS IN MAYBE?//
        private GridPoint[,,] grid; //3d grid
        private int xLen; //number of gridpoints
        private int yLen; //number of gridpoints
        private int zLen; //number of gridpoints
        private int deltaTime;
        private GridPoint targetPoint;
        private GridPoint infrontPoint; //i + 1
        private GridPoint behindPoint; //i - 1
        private GridPoint leftPoint; //j + 1
        private GridPoint rightPoint; //j - 1
        private GridPoint bottomPoint; //k - 1
        private GridPoint topPoint; //k + 1
        private double Zi;
        private int resolutionMultiplierX;
        private int resolutionMultiplierY;
        private int resolutionMultiplierZ;
        private double botLeftLat; //pass to DataReceiver {LATER} <----------
        private double botLeftLon; //pass to DataReceiver {LATER} <----------
        private double[,,] tempData;
        private double[,,] pressureData;
        private double[,,] mixingRatioData;
        private double[,,] densityData;
        private double[,,] windUData;
        private double[,,] windVData;
        private double[,,] windWData;
        private double[,,] deltaXData;
        private double[,,] deltaYData;
        private double[,,] deltaZData;
        private double[,,] latitudeData;
        private double[,,] longitudeData;
        private double[,,] geopotentialData;
        private double[,,] netRadiationData;
        private double[,,] soilTemperatureData;
        private double[,,] groundPressureData;
        private double[,,] relativeHumidityData;
        private double[,,] verticalVelocityData;

        public Grid(double inputBotLeftLat, double inputBotLeftLon, int deltaTime, double Zi, int resolutionMultiplierX, int resolutionMultiplierY, int resolutionMultiplierZ)
        {
            botLeftLat = inputBotLeftLat;
            botLeftLon = inputBotLeftLon;
            this.Zi = Zi;
            this.resolutionMultiplierX = resolutionMultiplierX;
            this.resolutionMultiplierY = resolutionMultiplierY;
            this.resolutionMultiplierZ = resolutionMultiplierZ;
            this.deltaTime = deltaTime;
        }

        public void deepCopy(Grid inputGrid)
        {
            for (int i = 0; i < xLen; i++)
            {
                for (int j = 0; j < yLen; j++)
                {
                    for (int k = 0; k < zLen; k++)
                    {
                        grid[i,j,k].inputUpdatedData(xLen,yLen,zLen, i,j,k,(inputGrid.getGrid()[i, j, k].getWindU()),
                            (inputGrid.getGrid()[i, j, k].getWindV()),(inputGrid.getGrid()[i, j, k].getWindW()),
                            (inputGrid.getGrid()[i, j, k].getAirDensity()),(inputGrid.getGrid()[i, j, k].getTemp()),
                            (inputGrid.getGrid()[i, j, k].getPressure()),
                            (inputGrid.getGrid()[i, j, k].getMixingRatio()),inputGrid.getGrid()[i,j,k].getDeltaX(),
                            inputGrid.getGrid()[i,j,k].getDeltaY(),inputGrid.getGrid()[i,j,k].getDeltaZ(),(inputGrid.getGrid()[i, j, k].getLatitude()),
                            (inputGrid.getGrid()[i, j, k].getLongitude()),
                            (inputGrid.getGrid()[i, j, k].getGeopotentialHeight()),
                            (inputGrid.getGrid()[i, j, k].getNetRadiation()),
                            (inputGrid.getGrid()[i, j, k].getSfcTemp()),
                            (inputGrid.getGrid()[i, j, k].getNearSfcPressure()),
                            (inputGrid.getGrid()[i, j, k].getRelativeHumidity()),
                            (inputGrid.getGrid()[i, j, k].getVerticalVelocity()));
                    }
                }
            }
        }
        
        private double[,,] intragridAverage(double[,,] inputGrid)
        {
            Boolean moreAvergingNeeded = true;
            double lastValue;
            int lastValueIndex;
            while (moreAvergingNeeded)
            {
                moreAvergingNeeded = false;
                for (int k = 0; k < zLen - 1; k = k + resolutionMultiplierZ)
                {
                    for (int i = 1; i < xLen - 1; i = i + resolutionMultiplierX)
                    {
                        lastValue = inputGrid[i, 1, k];
                        lastValueIndex = 1;
                        for (int j = 1; j < yLen - 1; j++)
                        {
                            if (Double.IsNaN(inputGrid[i,j,k]))
                            {
                                moreAvergingNeeded = true;
                            }
                            else if (j != 1 && moreAvergingNeeded && !Double.IsNaN(inputGrid[i,j,k]))
                            {
                                int avgIndex = ((j - lastValueIndex) / 2) + lastValueIndex;
                                inputGrid[i,avgIndex,k] = (lastValue + inputGrid[i,j,k]) / 2;
                                lastValue = inputGrid[i, j, k];
                                lastValueIndex = j;
                            }
                        }
                    }
                }
            }
            moreAvergingNeeded = true;
            while (moreAvergingNeeded)
            {
                moreAvergingNeeded = false;
                for (int k = 0; k < zLen - 1; k = k + resolutionMultiplierZ)
                {
                    for (int j = 1; j < yLen - 1; j++)
                    {
                        lastValue = inputGrid[1, j, k];
                        lastValueIndex = 1;
                        for (int i = 1; i < xLen - 1; i++)
                        {
                            if (Double.IsNaN(inputGrid[i, j, k]))
                            {
                                moreAvergingNeeded = true;
                            }
                            else if (i != 1 && moreAvergingNeeded && !Double.IsNaN(inputGrid[i,j,k]))
                            {
                                int avgIndex = ((i - lastValueIndex) / 2) + lastValueIndex;
                                inputGrid[avgIndex, j, k] = (lastValue + inputGrid[i, j, k]) / 2;
                                lastValue = inputGrid[i, j, k];
                                lastValueIndex = i;
                            }
                        }
                    }
                }
            }
            moreAvergingNeeded = true;
                while (moreAvergingNeeded)
                {
                    moreAvergingNeeded = false;
                    for (int i = 1; i < xLen - 1; i++)
                    {
                        for (int j = 1; j < yLen - 1; j++)
                        {
                            lastValue = inputGrid[i, j, 0];
                            lastValueIndex = 0;
                            for (int k = 0; k < zLen - 1; k++)
                            {
                                if (Double.IsNaN(inputGrid[i,j,k]))
                                {
                                    moreAvergingNeeded = true;
                                }
                                else if (k != 0 && moreAvergingNeeded && !Double.IsNaN(inputGrid[i,j,k]))
                                {
                                    int avgIndex = ((k - lastValueIndex) / 2) + lastValueIndex;
                                    inputGrid[i,j,avgIndex] = (lastValue + inputGrid[i,j,k]) / 2;
                                    lastValue = inputGrid[i, j, k];
                                    lastValueIndex = k;
                                }
                            }
                        }
                    }
                }
            return inputGrid;
        }

        private double[,,] intergridAverage(double[,,] presentGrid, double[,,] futureGrid)
        {
            int currentTime = timeBetweenGivenData;
            double[,,] avgGrid = new double[xLen,yLen,zLen];
            while (currentTime >= deltaTime)
            {
                currentTime = currentTime / 2;
                for (int i = 1; i < xLen - 1; i++)
                {
                    for (int j = 1; j < yLen - 1; j++)
                    {
                        for (int k = 0; k < zLen - 1; k++)
                        {
                            avgGrid[i, j, k] = (presentGrid[i, j, k] + futureGrid[i, j, k]) / 2;
                        }
                    }
                }
                futureGrid = avgGrid;
            }
            return avgGrid;
        }
        
        public void createEmptyGrid(Grid otherGrid)
        {
            xLen = otherGrid.getXLen();
            yLen = otherGrid.getYLen();
            zLen = otherGrid.getZLen();
            grid = new GridPoint[xLen,yLen,zLen];
            for (int i = 0; i < xLen; i++)
            {
                for (int j = 0; j < yLen; j++)
                {
                    for (int k = 0; k < zLen; k++)
                    {
                        grid[i, j, k] = new GridPoint(xLen, yLen, zLen, i, j, k);
                    }
                }
            }
        }
        
        //TODO: Remove temp comments
        //Data must be initialized from the ground up! (Sfc to ABL)
        public void populate(Boolean isAnalysis)
        {
            DataReceiver inputData = new DataReceiver(Zi,resolutionMultiplierX,resolutionMultiplierY,resolutionMultiplierZ,botLeftLat,botLeftLon);
            
            inputData.inputDataSet(isAnalysis);

            xLen = inputData.getXLen();
            yLen = inputData.getYLen();
            zLen = inputData.getZLen();
            tempData = inputData.getTempData();
            pressureData = inputData.getPressureData();
            mixingRatioData = inputData.getMixingRatioData();
            densityData = inputData.getDensityData();
            windUData = inputData.getWindUData();
            windVData = inputData.getWindVData();
            windWData = inputData.getWindWData();
            deltaXData = inputData.getDeltaXData();
            deltaYData = inputData.getDeltaYData();
            deltaZData = inputData.getDeltaZData();
            latitudeData = inputData.getLatitudeData();
            longitudeData = inputData.getLongitudeData();
            geopotentialData = inputData.getGeopotentialHeightData();
            netRadiationData = inputData.getNetRadiationData();
            soilTemperatureData = inputData.getSoilTemperatureData();
            groundPressureData = inputData.getGroundPressureData();
            relativeHumidityData = inputData.getRelativeHumidityData();
            verticalVelocityData = inputData.getVerticalVelocityData();

            tempData = intragridAverage(tempData);
            pressureData = intragridAverage(pressureData);
            mixingRatioData = intragridAverage(mixingRatioData);
            densityData = intragridAverage(densityData);
            windUData = intragridAverage(windUData);
            windVData = intragridAverage(windVData);
            windWData = intragridAverage(windWData);
            latitudeData = intragridAverage(latitudeData);
            longitudeData = intragridAverage(longitudeData);
            geopotentialData = intragridAverage(geopotentialData);
            netRadiationData = intragridAverage(netRadiationData);
            soilTemperatureData = intragridAverage(soilTemperatureData);
            groundPressureData = intragridAverage(groundPressureData);
            relativeHumidityData = intragridAverage(relativeHumidityData);
            verticalVelocityData = intragridAverage(verticalVelocityData);

            grid = new GridPoint[xLen,yLen,zLen];
            for (int i = 0; i < xLen; i++)
            {
                for (int j = 0; j < yLen; j++)
                {
                    for (int k = 0; k < zLen; k++)
                    {
                        grid[i,j,k] = new GridPoint(xLen,yLen,zLen,i,j,k,windUData[i,j,k],windVData[i,j,k],
                            windWData[i,j,k],densityData[i,j,k],tempData[i,j,k],
                            pressureData[i,j,k],mixingRatioData[i,j,k],deltaXData[i,j,k],deltaYData[i,j,k],
                            deltaZData[i,j,k],latitudeData[i,j,k],longitudeData[i,j,k],geopotentialData[i,j,k],
                            netRadiationData[i,j,k],soilTemperatureData[i,j,k],groundPressureData[i,j,k],
                            relativeHumidityData[i,j,k],verticalVelocityData[i,j,k]);
                    }
                }
            }
        }

        //Only used after populate() if this is future grid. Assumes present grid already populated
        public void futuregridTimeDividing(Grid presentGrid)
        {
            double[,,] tempDataPresent = presentGrid.getTempData();
            double[,,] pressureDataPresent = presentGrid.getPressureData();
            double[,,] mixingRatioDataPresent = presentGrid.getMixingRatioData();
            double[,,] densityDataPresent = presentGrid.getDensityData();
            double[,,] windUDataPresent = presentGrid.getWindUData();
            double[,,] windVDataPresent = presentGrid.getWindVData();
            double[,,] windWDataPresent = presentGrid.getWindWData();
            double[,,] geopotentialDataPresent = presentGrid.getGeopotentialHeightData();
            double[,,] netRadiationDataPresent = presentGrid.getNetRadiationData();
            double[,,] soilTemperatureDataPresent = presentGrid.getSoilTemperatureData();
            double[,,] groundPressureDataPresent = presentGrid.getGroundPressureData();
            double[,,] relativeHumidityDataPresent = presentGrid.getRelativeHumidityData();
            double[,,] verticalVelocityDataPresent = presentGrid.getVerticalVelocityData();
            
            tempData = intergridAverage(tempDataPresent,tempData);
            pressureData = intergridAverage(pressureDataPresent,pressureData);
            mixingRatioData = intergridAverage(mixingRatioDataPresent,mixingRatioData);
            densityData = intergridAverage(densityDataPresent,densityData);
            windUData = intergridAverage(windUDataPresent,windUData);
            windVData = intergridAverage(windVDataPresent, windVData);
            windWData = intergridAverage(windWDataPresent,windWData);
            geopotentialData = intergridAverage(geopotentialDataPresent,geopotentialData);
            netRadiationData = intergridAverage(netRadiationDataPresent, netRadiationData);
            soilTemperatureData = intergridAverage(soilTemperatureDataPresent, soilTemperatureData);
            groundPressureData = intergridAverage(groundPressureDataPresent, groundPressureData);
            relativeHumidityData = intergridAverage(relativeHumidityDataPresent, relativeHumidityData);
            verticalVelocityData = intergridAverage(verticalVelocityDataPresent, verticalVelocityData);
            
            for (int i = 0; i < xLen; i++)
            {
                for (int j = 0; j < yLen; j++)
                {
                    for (int k = 0; k < zLen; k++)
                    {
                        grid[i,j,k] = new GridPoint(xLen,yLen,zLen,i,j,k,windUData[i,j,k],windVData[i,j,k],
                            windWData[i,j,k],densityData[i,j,k],tempData[i,j,k],
                            pressureData[i,j,k],mixingRatioData[i,j,k],deltaXData[i,j,k],deltaYData[i,j,k],
                            deltaZData[i,j,k],latitudeData[i,j,k],longitudeData[i,j,k],geopotentialData[i,j,k],
                            netRadiationData[i,j,k],soilTemperatureData[i,j,k],groundPressureData[i,j,k],
                            relativeHumidityData[i,j,k],verticalVelocityData[i,j,k]);
                    }
                }
            }
        }

        public void createGhosts()
        {
            for (int j = 1; j < yLen - 1; j++)
            {
                for (int k = 1; k < zLen - 1; k++)
                {
                    grid[0,j,k].inputUpdatedData(xLen,yLen,zLen, 0,j,k,(2 * grid[1, j, k].getWindU() - grid[2, j, k].getWindU()),
                        (2 * grid[1, j, k].getWindV() - grid[2, j, k].getWindV()),(2 * grid[1, j, k].getWindW() - grid[2, j, k].getWindW()),
                        (2 * grid[1, j, k].getAirDensity() - grid[2, j, k].getAirDensity()),(2 * grid[1, j, k].getTemp() - grid[2, j, k].getTemp()),
                        (2 * grid[1, j, k].getPressure() - grid[2, j, k].getPressure()),
                        (2 * grid[1, j, k].getMixingRatio() - grid[2, j, k].getMixingRatio()),grid[1,j,k].getDeltaX(),
                        grid[1,j,k].getDeltaY(),grid[1,j,k].getDeltaZ(),(2 * grid[1, j, k].getLatitude() - grid[2, j, k].getLatitude()),
                        (2 * grid[1, j, k].getLongitude() - grid[2, j, k].getLongitude()),
                        (2 * grid[1, j, k].getGeopotentialHeight() - grid[2, j, k].getGeopotentialHeight()),
                        (2 * grid[1, j, k].getNetRadiation() - grid[2, j, k].getNetRadiation()),
                        (2 * grid[1, j, k].getSfcTemp() - grid[2, j, k].getSfcTemp()),
                        (2 * grid[1, j, k].getNearSfcPressure() - grid[2, j, k].getNearSfcPressure()),
                        (2 * grid[1, j, k].getRelativeHumidity() - grid[2, j, k].getRelativeHumidity()),
                        (2 * grid[1, j, k].getVerticalVelocity() - grid[2, j, k].getVerticalVelocity()));
                    
                    grid[xLen - 1,j,k].inputUpdatedData(xLen,yLen,zLen, xLen-1,j,k,(2 * grid[xLen - 2, j, k].getWindU() - grid[xLen - 3, j, k].getWindU()),
                        (2 * grid[xLen - 2, j, k].getWindV() - grid[xLen - 3, j, k].getWindV()),(2 * grid[xLen - 2, j, k].getWindW() - grid[xLen - 3, j, k].getWindW()),
                        (2 * grid[xLen - 2, j, k].getAirDensity() - grid[xLen - 3, j, k].getAirDensity()),(2 * grid[xLen - 2, j, k].getTemp() - grid[xLen - 3, j, k].getTemp()),
                        (2 * grid[xLen - 2, j, k].getPressure() - grid[xLen - 3, j, k].getPressure()),
                        (2 * grid[xLen - 2, j, k].getMixingRatio() - grid[xLen - 3, j, k].getMixingRatio()),grid[xLen - 2,j,k].getDeltaX(),
                        grid[xLen - 2,j,k].getDeltaY(),grid[xLen - 2,j,k].getDeltaZ(),(2 * grid[xLen - 2, j, k].getLatitude() - grid[xLen - 3, j, k].getLatitude()),
                        (2 * grid[xLen - 2, j, k].getLongitude() - grid[xLen - 3, j, k].getLongitude()),
                        (2 * grid[xLen - 2, j, k].getGeopotentialHeight() - grid[xLen - 3, j, k].getGeopotentialHeight()),
                        (2 * grid[xLen - 2, j, k].getNetRadiation() - grid[xLen - 3, j, k].getNetRadiation()),
                        (2 * grid[xLen - 2, j, k].getSfcTemp() - grid[xLen - 3, j, k].getSfcTemp()),
                        (2 * grid[xLen - 2, j, k].getNearSfcPressure() - grid[xLen - 3, j, k].getNearSfcPressure()),
                        (2 * grid[xLen - 2, j, k].getRelativeHumidity() - grid[xLen - 3, j, k].getRelativeHumidity()),
                        (2 * grid[xLen - 2, j, k].getVerticalVelocity() - grid[xLen - 3, j, k].getVerticalVelocity()));
                }
            }
            for (int i = 1; i < xLen - 1; i++)
            {
                for (int k = 1; k < zLen - 1; k++)
                {
                    grid[i,0,k].inputUpdatedData(xLen,yLen,zLen, i,0,k,(2 * grid[i, 1, k].getWindU() - grid[i, 2, k].getWindU()),
                        (2 * grid[i, 1, k].getWindV() - grid[i, 2, k].getWindV()),(2 * grid[i, 1, k].getWindW() - grid[i, 2, k].getWindW()),
                        (2 * grid[i, 1, k].getAirDensity() - grid[i, 2, k].getAirDensity()),(2 * grid[i, 1, k].getTemp() - grid[i, 2, k].getTemp()),
                        (2 * grid[i, 1, k].getPressure() - grid[i, 2, k].getPressure()),
                        (2 * grid[i, 1, k].getMixingRatio() - grid[i, 2, k].getMixingRatio()),grid[i,1,k].getDeltaX(),
                        grid[i,1,k].getDeltaY(),grid[i,1,k].getDeltaZ(),(2 * grid[i, 1, k].getLatitude() - grid[i, 2, k].getLatitude()),
                        (2 * grid[i, 1, k].getLongitude() - grid[i, 2, k].getLongitude()),
                        (2 * grid[i, 1, k].getGeopotentialHeight() - grid[i, 2, k].getGeopotentialHeight()),
                        (2 * grid[i, 1, k].getNetRadiation() - grid[i, 2, k].getNetRadiation()),
                        (2 * grid[i, 1, k].getSfcTemp() - grid[i, 2, k].getSfcTemp()),
                        (2 * grid[i, 1, k].getNearSfcPressure() - grid[i, 2, k].getNearSfcPressure()),
                        (2 * grid[i, 1, k].getRelativeHumidity() - grid[i, 2, k].getRelativeHumidity()),
                        (2 * grid[i, 1, k].getVerticalVelocity() - grid[i, 2, k].getVerticalVelocity()));
                    
                    grid[i,yLen-1,k].inputUpdatedData(xLen,yLen,zLen, i,yLen-1,k,(2 * grid[i, yLen-2, k].getWindU() - grid[i, yLen-3, k].getWindU()),
                        (2 * grid[i, yLen-2, k].getWindV() - grid[i,yLen-3, k].getWindV()),(2 * grid[i, yLen-2, k].getWindW() - grid[i,yLen - 3, k].getWindW()),
                        (2 * grid[i, yLen-2, k].getAirDensity() - grid[i,yLen - 3, k].getAirDensity()),(2 * grid[i, yLen-2, k].getTemp() - grid[i,yLen-3, k].getTemp()),
                        (2 * grid[i, yLen-2, k].getPressure() - grid[i,yLen-3, k].getPressure()),
                        (2 * grid[i, yLen-2, k].getMixingRatio() - grid[i,yLen-3, k].getMixingRatio()),grid[i, yLen-2, k].getDeltaX(),
                        grid[i, yLen-2, k].getDeltaY(),grid[i, yLen-2, k].getDeltaZ(),(2 * grid[i, yLen-2, k].getLatitude() - grid[i,yLen-3, k].getLatitude()),
                        (2 * grid[i, yLen-2, k].getLongitude() - grid[i,yLen-3, k].getLongitude()),
                        (2 * grid[i, yLen-2, k].getGeopotentialHeight() - grid[i,yLen-3, k].getGeopotentialHeight()),
                        (2 * grid[i, yLen-2, k].getNetRadiation() - grid[i,yLen-3, k].getNetRadiation()),
                        (2 * grid[i, yLen-2, k].getSfcTemp() - grid[i,yLen-3, k].getSfcTemp()),
                        (2 * grid[i, yLen-2, k].getNearSfcPressure() - grid[i,yLen-3, k].getNearSfcPressure()),
                        (2 * grid[i, yLen-2, k].getRelativeHumidity() - grid[i,yLen-3, k].getRelativeHumidity()),
                        (2 * grid[i, yLen-2, k].getVerticalVelocity() - grid[i,yLen-3, k].getVerticalVelocity()));
                }
            }
             for (int i = 1; i < xLen - 1; i++)
            {
                for (int j = 1; j < yLen - 1; j++)
                {
                    grid[i,j,0].inputUpdatedData(xLen,yLen,zLen, i,j,0,(2 * grid[i, j, 1].getWindU() - grid[i, j, 2].getWindU()),
                        (2 * grid[i, j, 1].getWindV() - grid[i, j, 2].getWindV()),(2 * grid[i, j, 1].getWindW() - grid[i, j, 2].getWindW()),
                        (2 * grid[i, j, 1].getAirDensity() - grid[i, j, 2].getAirDensity()),(2 * grid[i, j, 1].getTemp() - grid[i, j, 2].getTemp()),
                        (2 * grid[i, j, 1].getPressure() - grid[i, j, 2].getPressure()),
                        (2 * grid[i, j, 1].getMixingRatio() - grid[i, j, 2].getMixingRatio()),grid[i, j, 1].getDeltaX(),
                        grid[i, j, 1].getDeltaY(),grid[i, j, 1].getDeltaZ(),(2 * grid[i, j, 1].getLatitude() - grid[i, j, 2].getLatitude()),
                        (2 * grid[i, j, 1].getLongitude() - grid[i, j, 2].getLongitude()),
                        (2 * grid[i, j, 2].getGeopotentialHeight() - grid[i, j, 2].getGeopotentialHeight()),
                        (2 * grid[i, j, 1].getNetRadiation() - grid[i, j, 2].getNetRadiation()),
                        (2 * grid[i, j, 1].getSfcTemp() - grid[i, j, 2].getSfcTemp()),
                        (2 * grid[i, j, 1].getNearSfcPressure() - grid[i, j, 2].getNearSfcPressure()),
                        (2 * grid[i, j, 1].getRelativeHumidity() - grid[i, j, 2].getRelativeHumidity()),
                        (2 * grid[i, j, 1].getVerticalVelocity() - grid[i, j, 2].getVerticalVelocity()));
                    
                    grid[i,j,zLen-1].inputUpdatedData(xLen,yLen,zLen, i,j,zLen-1,(2 * grid[i, j, zLen-2].getWindU() - grid[i, j, zLen-3].getWindU()),
                        (2 * grid[i, j, zLen-2].getWindV() - grid[i,j, zLen-3].getWindV()),(2 * grid[i, j, zLen-2].getWindW() - grid[i,j, zLen-3].getWindW()),
                        (2 * grid[i, j, zLen-2].getAirDensity() - grid[i,j,zLen-3].getAirDensity()),(2 * grid[i, j, zLen-2].getTemp() - grid[i,j,zLen-3].getTemp()),
                        (2 * grid[i, j, zLen-2].getPressure() - grid[i,j,zLen-3].getPressure()),
                        (2 * grid[i, j, zLen-2].getMixingRatio() - grid[i,j,zLen-3].getMixingRatio()),grid[i, j, zLen-2].getDeltaX(),
                        grid[i, j, zLen-2].getDeltaY(),grid[i, j, zLen-2].getDeltaZ(),(2 * grid[i, j, zLen-2].getLatitude() - grid[i,j,zLen-3].getLatitude()),
                        (2 * grid[i, j, zLen-2].getLongitude() - grid[i,j,zLen-3].getLongitude()),
                        (2 * grid[i, j, zLen-2].getGeopotentialHeight() - grid[i,j,zLen-3].getGeopotentialHeight()),
                        (2 * grid[i, j, zLen-2].getNetRadiation() - grid[i,j,zLen-3].getNetRadiation()),
                        (2 * grid[i, j, zLen-2].getSfcTemp() - grid[i,j,zLen-3].getSfcTemp()),
                        (2 * grid[i, j, zLen-2].getNearSfcPressure() - grid[i,j,zLen-3].getNearSfcPressure()),
                        (2 * grid[i, j, zLen-2].getRelativeHumidity() - grid[i,j,zLen-3].getRelativeHumidity()),
                        (2 * grid[i, j, zLen-2].getVerticalVelocity() - grid[i,j,zLen-3].getVerticalVelocity()));
                }
            }
        }

        public void initGhosts()
        {
            for (int j = 1; j < yLen - 1; j++)
            {
                for (int k = 1; k < zLen - 1; k++)
                {
                    grid[0,j,k].inputUpdatedData(xLen,yLen,zLen, 0,j,k,(2 * grid[1, j, k].getWindU() - grid[2, j, k].getWindU()),
                        (2 * grid[1, j, k].getWindV() - grid[2, j, k].getWindV()),(2 * grid[1, j, k].getWindW() - grid[2, j, k].getWindW()),
                        (2 * grid[1, j, k].getAirDensity() - grid[2, j, k].getAirDensity()),(2 * grid[1, j, k].getTemp() - grid[2, j, k].getTemp()),
                        (2 * grid[1, j, k].getPressure() - grid[2, j, k].getPressure()),
                        (2 * grid[1, j, k].getMixingRatio() - grid[2, j, k].getMixingRatio()),grid[1,j,k].getDeltaX(),
                        grid[1,j,k].getDeltaY(),grid[1,j,k].getDeltaZ(),(2 * grid[1, j, k].getLatitude() - grid[2, j, k].getLatitude()),
                        (2 * grid[1, j, k].getLongitude() - grid[2, j, k].getLongitude()),
                        (2 * grid[1, j, k].getGeopotentialHeight() - grid[2, j, k].getGeopotentialHeight()),
                        (2 * grid[1, j, k].getNetRadiation() - grid[2, j, k].getNetRadiation()),
                        (2 * grid[1, j, k].getSfcTemp() - grid[2, j, k].getSfcTemp()),
                        (2 * grid[1, j, k].getNearSfcPressure() - grid[2, j, k].getNearSfcPressure()),
                        (2 * grid[1, j, k].getRelativeHumidity() - grid[2, j, k].getRelativeHumidity()),
                        (2 * grid[1, j, k].getVerticalVelocity() - grid[2, j, k].getVerticalVelocity()));
                    
                    grid[xLen - 1,j,k].inputUpdatedData(xLen,yLen,zLen, xLen-1,j,k,(2 * grid[xLen - 2, j, k].getWindU() - grid[xLen - 3, j, k].getWindU()),
                        (2 * grid[xLen - 2, j, k].getWindV() - grid[xLen - 3, j, k].getWindV()),(2 * grid[xLen - 2, j, k].getWindW() - grid[xLen - 3, j, k].getWindW()),
                        (2 * grid[xLen - 2, j, k].getAirDensity() - grid[xLen - 3, j, k].getAirDensity()),(2 * grid[xLen - 2, j, k].getTemp() - grid[xLen - 3, j, k].getTemp()),
                        (2 * grid[xLen - 2, j, k].getPressure() - grid[xLen - 3, j, k].getPressure()),
                        (2 * grid[xLen - 2, j, k].getMixingRatio() - grid[xLen - 3, j, k].getMixingRatio()),grid[xLen - 2,j,k].getDeltaX(),
                        grid[xLen - 2,j,k].getDeltaY(),grid[xLen - 2,j,k].getDeltaZ(),(2 * grid[xLen - 2, j, k].getLatitude() - grid[xLen - 3, j, k].getLatitude()),
                        (2 * grid[xLen - 2, j, k].getLongitude() - grid[xLen - 3, j, k].getLongitude()),
                        (2 * grid[xLen - 2, j, k].getGeopotentialHeight() - grid[xLen - 3, j, k].getGeopotentialHeight()),
                        (2 * grid[xLen - 2, j, k].getNetRadiation() - grid[xLen - 3, j, k].getNetRadiation()),
                        (2 * grid[xLen - 2, j, k].getSfcTemp() - grid[xLen - 3, j, k].getSfcTemp()),
                        (2 * grid[xLen - 2, j, k].getNearSfcPressure() - grid[xLen - 3, j, k].getNearSfcPressure()),
                        (2 * grid[xLen - 2, j, k].getRelativeHumidity() - grid[xLen - 3, j, k].getRelativeHumidity()),
                        (2 * grid[xLen - 2, j, k].getVerticalVelocity() - grid[xLen - 3, j, k].getVerticalVelocity()));
                }
            }
            for (int i = 1; i < xLen - 1; i++)
            {
                for (int k = 1; k < zLen - 1; k++)
                {
                    grid[i,0,k].inputUpdatedData(xLen,yLen,zLen, i,0,k,(2 * grid[i, 1, k].getWindU() - grid[i, 2, k].getWindU()),
                        (2 * grid[i, 1, k].getWindV() - grid[i, 2, k].getWindV()),(2 * grid[i, 1, k].getWindW() - grid[i, 2, k].getWindW()),
                        (2 * grid[i, 1, k].getAirDensity() - grid[i, 2, k].getAirDensity()),(2 * grid[i, 1, k].getTemp() - grid[i, 2, k].getTemp()),
                        (2 * grid[i, 1, k].getPressure() - grid[i, 2, k].getPressure()),
                        (2 * grid[i, 1, k].getMixingRatio() - grid[i, 2, k].getMixingRatio()),grid[i,1,k].getDeltaX(),
                        grid[i,1,k].getDeltaY(),grid[i,1,k].getDeltaZ(),(2 * grid[i, 1, k].getLatitude() - grid[i, 2, k].getLatitude()),
                        (2 * grid[i, 1, k].getLongitude() - grid[i, 2, k].getLongitude()),
                        (2 * grid[i, 1, k].getGeopotentialHeight() - grid[i, 2, k].getGeopotentialHeight()),
                        (2 * grid[i, 1, k].getNetRadiation() - grid[i, 2, k].getNetRadiation()),
                        (2 * grid[i, 1, k].getSfcTemp() - grid[i, 2, k].getSfcTemp()),
                        (2 * grid[i, 1, k].getNearSfcPressure() - grid[i, 2, k].getNearSfcPressure()),
                        (2 * grid[i, 1, k].getRelativeHumidity() - grid[i, 2, k].getRelativeHumidity()),
                        (2 * grid[i, 1, k].getVerticalVelocity() - grid[i, 2, k].getVerticalVelocity()));
                    
                    grid[i,yLen-1,k].inputUpdatedData(xLen,yLen,zLen, i,yLen-1,k,(2 * grid[i, yLen-2, k].getWindU() - grid[i, yLen-3, k].getWindU()),
                        (2 * grid[i, yLen-2, k].getWindV() - grid[i,yLen-3, k].getWindV()),(2 * grid[i, yLen-2, k].getWindW() - grid[i,yLen - 3, k].getWindW()),
                        (2 * grid[i, yLen-2, k].getAirDensity() - grid[i,yLen - 3, k].getAirDensity()),(2 * grid[i, yLen-2, k].getTemp() - grid[i,yLen-3, k].getTemp()),
                        (2 * grid[i, yLen-2, k].getPressure() - grid[i,yLen-3, k].getPressure()),
                        (2 * grid[i, yLen-2, k].getMixingRatio() - grid[i,yLen-3, k].getMixingRatio()),grid[i, yLen-2, k].getDeltaX(),
                        grid[i, yLen-2, k].getDeltaY(),grid[i, yLen-2, k].getDeltaZ(),(2 * grid[i, yLen-2, k].getLatitude() - grid[i,yLen-3, k].getLatitude()),
                        (2 * grid[i, yLen-2, k].getLongitude() - grid[i,yLen-3, k].getLongitude()),
                        (2 * grid[i, yLen-2, k].getGeopotentialHeight() - grid[i,yLen-3, k].getGeopotentialHeight()),
                        (2 * grid[i, yLen-2, k].getNetRadiation() - grid[i,yLen-3, k].getNetRadiation()),
                        (2 * grid[i, yLen-2, k].getSfcTemp() - grid[i,yLen-3, k].getSfcTemp()),
                        (2 * grid[i, yLen-2, k].getNearSfcPressure() - grid[i,yLen-3, k].getNearSfcPressure()),
                        (2 * grid[i, yLen-2, k].getRelativeHumidity() - grid[i,yLen-3, k].getRelativeHumidity()),
                        (2 * grid[i, yLen-2, k].getVerticalVelocity() - grid[i,yLen-3, k].getVerticalVelocity()));
                }
            }
             for (int i = 1; i < xLen - 1; i++)
            {
                for (int j = 1; j < yLen - 1; j++)
                {
                    grid[i,j,zLen-1].inputUpdatedData(xLen,yLen,zLen, i,j,zLen-1,(2 * grid[i, j, zLen-2].getWindU() - grid[i, j, zLen-3].getWindU()),
                        (2 * grid[i, j, zLen-2].getWindV() - grid[i,j, zLen-3].getWindV()),(2 * grid[i, j, zLen-2].getWindW() - grid[i,j, zLen-3].getWindW()),
                        (2 * grid[i, j, zLen-2].getAirDensity() - grid[i,j,zLen-3].getAirDensity()),(2 * grid[i, j, zLen-2].getTemp() - grid[i,j,zLen-3].getTemp()),
                        (2 * grid[i, j, zLen-2].getPressure() - grid[i,j,zLen-3].getPressure()),
                        (2 * grid[i, j, zLen-2].getMixingRatio() - grid[i,j,zLen-3].getMixingRatio()),grid[i, j, zLen-2].getDeltaX(),
                        grid[i, j, zLen-2].getDeltaY(),grid[i, j, zLen-2].getDeltaZ(),(2 * grid[i, j, zLen-2].getLatitude() - grid[i,j,zLen-3].getLatitude()),
                        (2 * grid[i, j, zLen-2].getLongitude() - grid[i,j,zLen-3].getLongitude()),
                        (2 * grid[i, j, zLen-2].getGeopotentialHeight() - grid[i,j,zLen-3].getGeopotentialHeight()),
                        (2 * grid[i, j, zLen-2].getNetRadiation() - grid[i,j,zLen-3].getNetRadiation()),
                        (2 * grid[i, j, zLen-2].getSfcTemp() - grid[i,j,zLen-3].getSfcTemp()),
                        (2 * grid[i, j, zLen-2].getNearSfcPressure() - grid[i,j,zLen-3].getNearSfcPressure()),
                        (2 * grid[i, j, zLen-2].getRelativeHumidity() - grid[i,j,zLen-3].getRelativeHumidity()),
                        (2 * grid[i, j, zLen-2].getVerticalVelocity() - grid[i,j,zLen-3].getVerticalVelocity()));
                }
            }
        }

        public void tempCreateGhosts()
        {
            for (int j = 1; j < yLen - 1; j++)
            {
                for (int k = 0; k < zLen - 1; k++)
                {
                    grid[0,j,k].inputUpdatedData(xLen,yLen,zLen, 0,j,k,(grid[1, j, k].getWindU()),
                        (grid[1, j, k].getWindV()),(grid[1, j, k].getWindW()),
                        (grid[1, j, k].getAirDensity()),(grid[1, j, k].getTemp()),
                        (grid[1, j, k].getPressure()),
                        (grid[1, j, k].getMixingRatio()),grid[1,j,k].getDeltaX(),
                        grid[1,j,k].getDeltaY(),grid[1,j,k].getDeltaZ(),(grid[1, j, k].getLatitude()),
                        (grid[1, j, k].getLongitude()),
                        (grid[1, j, k].getGeopotentialHeight()),
                        (grid[1, j, k].getNetRadiation()),
                        (grid[1, j, k].getSfcTemp()),
                        (grid[1, j, k].getNearSfcPressure()),
                        (grid[1, j, k].getRelativeHumidity()),
                        (grid[1, j, k].getVerticalVelocity()));
                    
                    grid[xLen-1,j,k].inputUpdatedData(xLen,yLen,zLen, xLen-1,j,k,(grid[xLen-2, j, k].getWindU()),
                        (grid[xLen-2, j, k].getWindV()),(grid[xLen-2, j, k].getWindW()),
                        (grid[xLen-2, j, k].getAirDensity()),(grid[xLen-2, j, k].getTemp()),
                        (grid[xLen-2, j, k].getPressure()),
                        (grid[xLen-2, j, k].getMixingRatio()),grid[xLen-2, j, k].getDeltaX(),
                        grid[xLen-2, j, k].getDeltaY(),grid[xLen-2, j, k].getDeltaZ(),(grid[xLen-2, j, k].getLatitude()),
                        (grid[xLen-2, j, k].getLongitude()),
                        (grid[xLen-2, j, k].getGeopotentialHeight()),
                        (grid[xLen-2, j, k].getNetRadiation()),
                        (grid[xLen-2, j, k].getSfcTemp()),
                        (grid[xLen-2, j, k].getNearSfcPressure()),
                        (grid[xLen-2, j, k].getRelativeHumidity()),
                        (grid[xLen-2, j, k].getVerticalVelocity()));
                }
            }
            for (int i = 1; i < xLen - 1; i++)
            {
                for (int k = 0; k < zLen - 1; k++)
                {
                    grid[i,0,k].inputUpdatedData(xLen,yLen,zLen, i,0,k,(grid[i, 1, k].getWindU()),
                        (grid[i, 1, k].getWindV()),(grid[i, 1, k].getWindW()),
                        (grid[i, 1, k].getAirDensity()),(grid[i, 1, k].getTemp()),
                        (grid[i, 1, k].getPressure()),
                        (grid[i, 1, k].getMixingRatio()),grid[i, 1, k].getDeltaX(),
                        grid[i, 1, k].getDeltaY(),grid[i, 1, k].getDeltaZ(),(grid[i, 1, k].getLatitude()),
                        (grid[i, 1, k].getLongitude()),
                        (grid[i, 1, k].getGeopotentialHeight()),
                        (grid[i, 1, k].getNetRadiation()),
                        (grid[i, 1, k].getSfcTemp()),
                        (grid[i, 1, k].getNearSfcPressure()),
                        (grid[i, 1, k].getRelativeHumidity()),
                        (grid[i, 1, k].getVerticalVelocity()));
                    
                    grid[i,yLen-1,k].inputUpdatedData(xLen,yLen,zLen, i,yLen-1,k,(grid[i, yLen-2, k].getWindU()),
                        (grid[i, yLen-2, k].getWindV()),(grid[i, yLen-2, k].getWindW()),
                        (grid[i, yLen-2, k].getAirDensity()),(grid[i, yLen-2, k].getTemp()),
                        (grid[i, yLen-2, k].getPressure()),
                        (grid[i, yLen-2, k].getMixingRatio()),grid[i, yLen-2, k].getDeltaX(),
                        grid[i, yLen-2, k].getDeltaY(),grid[i, yLen-2, k].getDeltaZ(),(grid[i, yLen-2, k].getLatitude()),
                        (grid[i, yLen-2, k].getLongitude()),
                        (grid[i, yLen-2, k].getGeopotentialHeight()),
                        (grid[i, yLen-2, k].getNetRadiation()),
                        (grid[i, yLen-2, k].getSfcTemp()),
                        (grid[i, yLen-2, k].getNearSfcPressure()),
                        (grid[i, yLen-2, k].getRelativeHumidity()),
                        (grid[i, yLen-2, k].getVerticalVelocity()));
                }
            }
             for (int i = 1; i < xLen - 1; i++)
            {
                for (int j = 1; j < yLen - 1; j++)
                {
                    grid[i,j,zLen-1].inputUpdatedData(xLen,yLen,zLen, i,j,zLen-1,(grid[i, j, zLen-2].getWindU()),
                        (grid[i, j, zLen-2].getWindV()),(grid[i, j, zLen-2].getWindW()),
                        (grid[i, j, zLen-2].getAirDensity()),(grid[i, j, zLen-2].getTemp()),
                        (grid[i, j, zLen-2].getPressure()),
                        (grid[i, j, zLen-2].getMixingRatio()),grid[i, j, zLen-2].getDeltaX(),
                        grid[i, j, zLen-2].getDeltaY(),grid[i, j, zLen-2].getDeltaZ(),(grid[i, j, zLen-2].getLatitude()),
                        (grid[i, j, zLen-2].getLongitude()),
                        (grid[i, j, zLen-2].getGeopotentialHeight()),
                        (grid[i, j, zLen-2].getNetRadiation()),
                        (grid[i, j, zLen-2].getSfcTemp()),
                        (grid[i, j, zLen-2].getNearSfcPressure()),
                        (grid[i, j, zLen-2].getRelativeHumidity()),
                        (grid[i, j, zLen-2].getVerticalVelocity()));
                }
            }
        }
        
        public void setTarget(int targetI, int targetJ, int targetK)
        {
            targetPoint = grid[targetI, targetJ, targetK];
            behindPoint = grid[targetI - 1, targetJ, targetK];
            infrontPoint = grid[targetI + 1, targetJ, targetK];
            rightPoint = grid[targetI, targetJ - 1, targetK];
            leftPoint = grid[targetI, targetJ + 1, targetK];
            topPoint = grid[targetI, targetJ, targetK + 1];
            if (targetK != 0)
            {
                bottomPoint = grid[targetI, targetJ, targetK - 1];
            }
            if (targetPoint.getRelativeHumidity() >= 95)
            {
                //Console.WriteLine("Warning: Model humidity exceeding expected parameters. Intervention required."); //TEMP COMMENTED OUT///
            }
        }

        public void setXLen(int xLen)
        {
            this.xLen = xLen;
        }

        public void setYLen(int yLen)
        {
            this.yLen = yLen;
        }

        public void setZLen(int zLen)
        {
            this.zLen = zLen;
        }
        
        public GridPoint getTargetPoint()
        {
            return targetPoint;
        }
        public GridPoint getInfrontPoint()
        {
            return infrontPoint;
        }
        public GridPoint getBehindPoint()
        {
            return behindPoint;
        }
        public GridPoint getLeftPoint()
        {
            return leftPoint;
        }
        public GridPoint getRightPoint()
        {
            return rightPoint;
        }
        public GridPoint getTopPoint()
        {
            return topPoint;
        }
        public GridPoint getBottomPoint()
        {
            return bottomPoint;
        }
        public GridPoint getPoint(int i, int j, int k)
        {
            return grid[i, j, k];
        }
        public double[,,] getTempData()
        {
            return tempData;
        }
        public double[,,] getPressureData()
        {
            return pressureData;
        }
        public double[,,] getMixingRatioData()
        {
            return mixingRatioData;
        }
        public double[,,] getDensityData()
        {
            return densityData;
        }
        public double[,,] getWindUData()
        {
            return windUData;
        }
        public double[,,] getWindVData()
        {
            return windVData;
        }
        public double[,,] getWindWData()
        {
            return windWData;
        }
        public double[,,] getDeltaXData()
        {
            return deltaXData;
        }
        public double[,,] getDeltaYData()
        {
            return deltaYData;
        }
        public double[,,] getDeltaZData()
        {
            return deltaZData;
        }
        public double[,,] getLatitudeData()
        {
            return latitudeData;
        }
        public double[,,] getLongitudeData()
        {
            return longitudeData;
        }
        public double[,,] getGeopotentialHeightData()
        {
            return geopotentialData;
        }
        public double[,,] getNetRadiationData()
        {
            return netRadiationData;
        }
        public double[,,] getSoilTemperatureData()
        {
            return soilTemperatureData;
        }
        public double[,,] getGroundPressureData()
        {
            return groundPressureData;
        }
        public double[,,] getRelativeHumidityData()
        {
            return relativeHumidityData;
        }
        public double[,,] getVerticalVelocityData()
        {
            return verticalVelocityData;
        }
        public int getXLen()
        {
            return xLen;
        }
        public int getYLen()
        {
            return yLen;
        }
        public int getZLen()
        {
            return zLen;
        }
        public GridPoint[,,] getGrid()
        {
            return grid;
        }
    }
    
    public partial class MainWindow : Window
    {
        public MainWindow(List<int> timeData, List<Double> inputData)
        {
            
            DataContext = this;
            this.Title = "test";
            var tmp = new PlotModel {Title = "Test", Subtitle = "axis"};
            var s1 = new LineSeries
            {
                StrokeThickness = 0,
                MarkerSize = 3,
                //MarkerFill = OxyColors.Blue,
                MarkerStroke=OxyColors.ForestGreen,
                MarkerType = MarkerType.Plus
            };

            for (int i = 0; i < inputData.Count; i++)
            {
                s1.Points.Add(new DataPoint(timeData[i],inputData[i]));
            }

            //s1.Show();
            tmp.Series.Add(s1);
            //tmp.Show();
            this.ScatterModel = tmp;
        }
        public PlotModel ScatterModel { get; set; }
    }
    
    //TODO: SOLVE NUMERICAL INSTABILITY
    public class MainClass
    {
        public const double Zi = 2000; //ABL depth (MUST ALIGN WITH deltaZ and Zlen!!!!) (WAS 2000)
        public const double botLeftLat = 36.5; //I am using kansas
        public const double botLeftLon = -102.5; //I am using kansas
        public const int deltaTime = 1;
        public const int forcastEndTime = 60;
        public const int resolutionMultiplierX = 8; //Determines how many gridpoints will be placed between the given 0.25 lon points (WAS 8) <---------
        public const int resolutionMultiplierY = 16; //Determines how many gridpoints will be placed between the given 0.25 lat points (WAS 16) <---------
        public const int resolutionMultiplierZ = 2; //Determines how many gridpoints will be placed between the given 25 mb z points (WAS 2) <-------------
        
        [STAThread]
        static void Main(string[] args)
        {
            int currentTime = 0;
            Console.WriteLine("Loading all variables...");
            Grid pastGrid = new Grid(botLeftLat,botLeftLon,deltaTime,Zi,resolutionMultiplierX,resolutionMultiplierY,resolutionMultiplierZ);
            Grid presentGrid = new Grid(botLeftLat,botLeftLon,deltaTime,Zi,resolutionMultiplierX,resolutionMultiplierY,resolutionMultiplierZ);
            Grid futureGrid = new Grid(botLeftLat, botLeftLon, deltaTime, Zi, resolutionMultiplierX,resolutionMultiplierY, resolutionMultiplierZ);
            List<double> data = new List<double>(); //TEMP///
            List<int> time = new List<int>(); //TEMP////
            pastGrid.populate(true); //For creating past grid
            presentGrid.populate(false); //For future grid
            presentGrid.futuregridTimeDividing(pastGrid); //Brings the grid to a couple seconds into future
            double xLen = presentGrid.getXLen();
            double yLen = presentGrid.getYLen();
            double zLen = presentGrid.getZLen();
            //pastGrid.createGhosts();
            //presentGrid.createGhosts();
            pastGrid.tempCreateGhosts();
            presentGrid.tempCreateGhosts();
            //pastGrid.initGhosts();
            //presentGrid.initGhosts();
            futureGrid.createEmptyGrid(pastGrid); //Gets the dimensions of grid from other grids. (Maybe reconsider using "presentGrid")
            Calculator calc = new Calculator(Zi);
            Console.WriteLine("Load complete. Beginning forecast...");
            
            presentGrid.setTarget(40,40,5);
            GridPoint inputPoint = presentGrid.getTargetPoint();
            //Console.WriteLine(inputPoint.getLatitude());
            //Console.WriteLine(inputPoint.getLongitude());
            //Console.WriteLine(inputPoint.getMixingRatio());
            double later = inputPoint.getTemp(); //TEMP///////
            double later2 = inputPoint.getWindU();
            
            while (currentTime != forcastEndTime)
            {
                for (int i = 1; i < xLen - 1; i++)
                {
                    for (int j = 1; j < yLen - 1; j++)
                    {
                        for (int k = 0; k < zLen - 1; k++)
                        {
                            pastGrid.setTarget(i, j, k);
                            presentGrid.setTarget(i, j, k);
                            futureGrid.setTarget(i, j, k);
                            calc.update(presentGrid,pastGrid,futureGrid); //Inputs the new grids into the calculator
                            calc.initTurbTemp();
                            calc.initTurbMixingRatioandUVWind();
                            //calc.temp();
                            //if (k != 0)
                            //{
                                calc.calcTemp(deltaTime);
                                calc.calcMixingRatio(deltaTime);
                                calc.calcWindU(deltaTime);
                                calc.calcWindV(deltaTime);
                                //calc.calcWindW(deltaTime);
                                calc.calcAirDensity(deltaTime);
                                calc.calcAirPressure(deltaTime);
                                calc.calcRest();
                            //}
                            futureGrid = calc.getFutureGrid();
                            presentGrid = calc.getPresentGrid();
                            GridPoint target = presentGrid.getTargetPoint();
                            //////////////////////////////////////TEMP/////////////////////////////////////
                            if (Double.IsNaN(target.getWindW()) || Double.IsNaN(target.getWindU()) ||
                                Double.IsNaN(target.getWindV()) || Double.IsNaN(target.getAirDensity()) ||
                                Double.IsNaN(target.getTemp()) || Double.IsNaN(target.getPressure()) ||
                                Double.IsNaN(target.getMixingRatio()) || Double.IsNaN(target.getTurbMixingRatio()) ||
                                Double.IsNaN(target.getTurbTemp()) || Double.IsNaN(target.getTurbUWind()) ||
                                Double.IsNaN(target.getTurbVWind()))
                            {
                                Console.WriteLine("Problem");
                                Console.WriteLine(target.getICoord());
                                Console.WriteLine(target.getJCoord());
                                Console.WriteLine(target.getKCoord());
                                Console.WriteLine(target.getLatitude());
                                Console.WriteLine(target.getLongitude());
                                Console.WriteLine(target.getWindU());
                                Console.WriteLine(target.getWindV());
                                Console.WriteLine(target.getWindW());
                                Console.WriteLine(target.getAirDensity());
                                Console.WriteLine(target.getTemp());
                                Console.WriteLine(target.getPressure());
                                Console.WriteLine(target.getMixingRatio());
                                Console.WriteLine(target.getTurbMixingRatio());
                                Environment.Exit(0);
                            }
                            if (deltaTime > target.getDeltaX() / Math.Abs(target.getWindU()))
                            {
                                Console.WriteLine("MASSIVE ERROR IN U");
                                Console.WriteLine(target.getICoord());
                                Console.WriteLine(target.getJCoord());
                                Console.WriteLine(target.getKCoord());
                                Console.WriteLine(target.getLatitude());
                                Console.WriteLine(target.getLongitude());
                                Console.WriteLine(target.getWindU());
                                Console.WriteLine(target.getWindV());
                                Console.WriteLine(target.getWindW());
                                Console.WriteLine(target.getAirDensity());
                                Console.WriteLine(target.getTemp());
                                Console.WriteLine(target.getPressure());
                                Console.WriteLine(target.getMixingRatio());
                                Console.WriteLine(target.getTurbMixingRatio());
                                Environment.Exit(0);
                            }
                            else if (deltaTime > target.getDeltaY() / Math.Abs(target.getWindV()))
                            {
                                Console.WriteLine("MASSIVE ERROR IN V");
                                Console.WriteLine(target.getICoord());
                                Console.WriteLine(target.getJCoord());
                                Console.WriteLine(target.getKCoord());
                                Console.WriteLine(target.getLatitude());
                                Console.WriteLine(target.getLongitude());
                                Console.WriteLine(target.getWindU());
                                Console.WriteLine(target.getWindV());
                                Console.WriteLine(target.getWindW());
                                Console.WriteLine(target.getAirDensity());
                                Console.WriteLine(target.getTemp());
                                Console.WriteLine(target.getPressure());
                                Console.WriteLine(target.getMixingRatio());
                                Console.WriteLine(target.getTurbMixingRatio());
                                Environment.Exit(0);
                            }
                            else if (deltaTime > target.getDeltaZ() / Math.Abs(target.getWindW()))
                            {
                                Console.WriteLine("MASSIVE ERROR IN W");
                                Console.WriteLine(target.getICoord());
                                Console.WriteLine(target.getJCoord());
                                Console.WriteLine(target.getKCoord());
                                Console.WriteLine(target.getLatitude());
                                Console.WriteLine(target.getLongitude());
                                Console.WriteLine(target.getWindU());
                                Console.WriteLine(target.getWindV());
                                Console.WriteLine(target.getWindW());
                                Console.WriteLine(target.getAirDensity());
                                Console.WriteLine(target.getTemp());
                                Console.WriteLine(target.getPressure());
                                Console.WriteLine(target.getMixingRatio());
                                Console.WriteLine(target.getTurbMixingRatio());
                                Environment.Exit(0);
                            }
                            ///////////////////////////////////////////////////////////////////////////
                        }
                    }
                }
                //futureGrid.createGhosts();
                futureGrid.tempCreateGhosts();
                pastGrid.deepCopy(presentGrid);
                presentGrid.deepCopy(futureGrid);
                pastGrid.setTarget(40,40,5); //TEMP/////
                Console.WriteLine(pastGrid.getTargetPoint().getWindU());
                Console.WriteLine(pastGrid.getTargetPoint().getTemp());
                currentTime = currentTime + deltaTime;
                Console.WriteLine(currentTime); //TEMP//////////////
                if (currentTime % 5 == 0)
                {
                    data.Add(pastGrid.getTargetPoint().getWindU());
                    time.Add(currentTime);
                }
            }
            
            pastGrid.setTarget(40,40,5);
            GridPoint outputPoint = pastGrid.getTargetPoint();
            Console.WriteLine("Input Temp:");
            Console.WriteLine(later); //TEMP///////////
            Console.WriteLine(later2);
            Console.WriteLine("------------");
            Console.WriteLine("Output Temp:");
            Console.WriteLine(outputPoint.getTemp());
            Console.WriteLine(outputPoint.getWindU());
            MainWindow derp = new MainWindow(time, data);
            derp.ShowDialog();
        }
    }
}