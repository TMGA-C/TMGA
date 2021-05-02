#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "NGA.h"
#include "TMGA.h"
#include "CGA.h"
#include "LWSGA.h"
#include "HGA.h"
#include "HEFT.h"
#include "IHEFT3.h"

using namespace std;

int main() {
    srand((int) time(0));
    //{set the runtime (termination) time of the algorithm} -xy4
    map<string, double> SchTime;
    SchTime["Montage25_0.4"] =  1.128 ;
    SchTime["Montage25_0.7"] =  1.430 ;
    SchTime["Montage25_1.0"] =  1.915 ;
    SchTime["Montage100_0.4"] = 8.657 ;
    SchTime["Montage100_0.7"] = 17.988 ;
    SchTime["Montage100_1.0"] = 37.370 ;
    SchTime["Montage1000_0.4"] =567.317 ;
    SchTime["Montage1000_0.7"] =944.537 ;
    SchTime["Montage1000_1.0"] =2231.823 ;

    SchTime["Epigenomics24_0.4"] = 1.087 ;
    SchTime["Epigenomics24_0.7"] = 1.322 ;
    SchTime["Epigenomics24_1.0"] = 2.064 ;
    SchTime["Epigenomics100_0.4"] =6.920 ;
    SchTime["Epigenomics100_0.7"] =19.656 ;
    SchTime["Epigenomics100_1.0"] =55.043 ;
    SchTime["Epigenomics997_0.4"] =635.331 ;
    SchTime["Epigenomics997_0.7"] =1284.426 ;
    SchTime["Epigenomics997_1.0"] =1722.925 ;

    SchTime["Ligo30_0.4"] = 1.601 ;
    SchTime["Ligo30_0.7"] = 1.922 ;
    SchTime["Ligo30_1.0"] = 2.926 ;
    SchTime["Ligo100_0.4"] =7.144 ;
    SchTime["Ligo100_0.7"] =11.420 ;
    SchTime["Ligo100_1.0"] =21.433 ;
    SchTime["Ligo1000_0.4"]=349.066 ;
    SchTime["Ligo1000_0.7"]=674.353 ;
    SchTime["Ligo1000_1.0"]=1133.507 ;
    //{clear "result"}
    ofstream outfile("../result.txt", ios::out);
    outfile.close();

    string Model, NumOfTask, RscAvlRatio;
    do {
        string StrLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(0);
        }
        getline(iFile, StrLine);
        if (StrLine.size() < 1) {
            cout << "Empty input file" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(StrLine);
        is >> Model >> NumOfTask >> RscAvlRatio;  //NumOfTask, RscAvlRatio
        XmlFile = Model + "_" + NumOfTask + "_0.xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";

        cout <<endl<< Model << " " << NumOfTask << " " << RscAvlRatio << " ";
//        double HEFT_SchTime  = 0;
//        double HEFT_Result = runHEFT(XmlFile, RscAlcFile,HEFT_SchTime);
//
//        double IHEFT3_SchTime  = 0;
//        double IHEFT3_Result = runIHEFT3(XmlFile, RscAlcFile,IHEFT3_SchTime);

//        double HGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int HGA_Iteration = 0;
//        double HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration);
//
//        double NGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int NGA_Iteration = 0;
//        double NGA_Result = runNGA(XmlFile, RscAlcFile, NGA_SchTime, NGA_Iteration);
//
//        double LWSGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int LWSGA_Iteration = 0;
//        double LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration);
//
//        double CGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int CGA_Iteration = 0;
//        double CGA_Result = runCGA(XmlFile, RscAlcFile, CGA_SchTime, CGA_Iteration);
//
        double TMGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int TMGA_Iteration = 0;
        double TMGA_Result = runTMGA(XmlFile, RscAlcFile, TMGA_SchTime, TMGA_Iteration);

        //results are written into the file
        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(5);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << HEFT_Result << " " << HEFT_SchTime << " "
//                << IHEFT3_Result << " " << IHEFT3_SchTime << " "
//                << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " "
//                << NGA_Result << " " << NGA_SchTime << " " << NGA_Iteration << " "
//                << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " "
//                << CGA_Result << " " << CGA_SchTime << " " << CGA_Iteration << " "
                << TMGA_Result  << " " << TMGA_SchTime  << " " << TMGA_Iteration  << " "
                << endl;
        outfile.close();
        //delete the first line in the file
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
    return 0;
}
