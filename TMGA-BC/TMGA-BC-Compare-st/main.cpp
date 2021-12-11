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

using namespace std;

int main() {
    srand((int) time(0));
    //{set the runtime (termination) time of the algorithm}
    map<string, double> SchTime;
    SchTime["Montage25_0.4"] =  1.026  ;
    SchTime["Montage25_0.7"] =  1.673  ;
    SchTime["Montage25_1.0"] =  2.091  ;
    SchTime["Montage100_0.4"] = 9.477  ;
    SchTime["Montage100_0.7"] = 24.268 ;
    SchTime["Montage100_1.0"] = 33.797  ;
    SchTime["Montage1000_0.4"] =656.780 ;
    SchTime["Montage1000_0.7"] =983.592 ;
    SchTime["Montage1000_1.0"] =1548.255 ;

    SchTime["Epigenomics24_0.4"] = 1.163  ;
    SchTime["Epigenomics24_0.7"] = 1.356  ;
    SchTime["Epigenomics24_1.0"] = 2.248  ;
    SchTime["Epigenomics100_0.4"] =8.055  ;
    SchTime["Epigenomics100_0.7"] =16.558  ;
    SchTime["Epigenomics100_1.0"] =26.689  ;
    SchTime["Epigenomics997_0.4"] =592.233  ;
    SchTime["Epigenomics997_0.7"] =857.129 ;
    SchTime["Epigenomics997_1.0"] =1008.571 ;

    SchTime["Ligo30_0.4"] = 1.486 ;
    SchTime["Ligo30_0.7"] = 2.038 ;
    SchTime["Ligo30_1.0"] = 3.229 ;
    SchTime["Ligo100_0.4"] =7.826 ;
    SchTime["Ligo100_0.7"] =11.632  ;
    SchTime["Ligo100_1.0"] =19.225  ;
    SchTime["Ligo1000_0.4"]=379.681 ;
    SchTime["Ligo1000_0.7"]=602.163 ;
    SchTime["Ligo1000_1.0"]=951.419 ;


    map<string, double> EndFit;
    EndFit["Montage25_0.4"] =  42.42 ;
    EndFit["Montage25_0.7"] =  29.92 ;
    EndFit["Montage25_1.0"] =  29.72 ;
    EndFit["Montage100_0.4"] = 149.25  ;
    EndFit["Montage100_0.7"] = 142.84  ;
    EndFit["Montage100_1.0"] = 104.49  ;
    EndFit["Montage1000_0.4"] =1164.06  ;
    EndFit["Montage1000_0.7"] =1061.26  ;
    EndFit["Montage1000_1.0"] =1068.74  ;

    EndFit["Epigenomics24_0.4"] = 3179.97 ;
    EndFit["Epigenomics24_0.7"] = 2321.04 ;
    EndFit["Epigenomics24_1.0"] = 2320.91 ;
    EndFit["Epigenomics100_0.4"] =50457.03  ;
    EndFit["Epigenomics100_0.7"] =38885.37  ;
    EndFit["Epigenomics100_1.0"] =36933.52  ;
    EndFit["Epigenomics997_0.4"] =328901.02  ;
    EndFit["Epigenomics997_0.7"] =326692.65  ;
    EndFit["Epigenomics997_1.0"] =325704.84  ;

    EndFit["Ligo30_0.4"] =  995.50  ;
    EndFit["Ligo30_0.7"] =  692.69  ;
    EndFit["Ligo30_1.0"] =  662.39  ;
    EndFit["Ligo100_0.4"] = 2054.44  ;
    EndFit["Ligo100_0.7"] = 1863.73  ;
    EndFit["Ligo100_1.0"] = 1801.10  ;
    EndFit["Ligo1000_0.4"] =19731.92  ;
    EndFit["Ligo1000_0.7"] =19114.54  ;
    EndFit["Ligo1000_1.0"] =19040.83 ;

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

        //results are written into the file
        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(5);

//        Due to the use of multi-threaded parallel computing, each algorithm should be run separately in order to get accurate results

        double HEFT_SchTime  = 0;
        double HEFT_Result = runHEFT(XmlFile, RscAlcFile,HEFT_SchTime);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << HEFT_Result << " " << HEFT_SchTime << " "  << endl;

        double EndFitness = EndFit[Model + NumOfTask + "_" + RscAvlRatio];

        double HGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int HGA_Iteration = 0;
        double HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration , EndFitness);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " " << endl;

        double NGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int NGA_Iteration = 0;
        double NGA_Result = runNGA(XmlFile, RscAlcFile, NGA_SchTime, NGA_Iteration , EndFitness);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << NGA_Result << " " << NGA_SchTime << " " << NGA_Iteration << " " << endl;

        double LWSGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int LWSGA_Iteration = 0;
        double LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration , EndFitness);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " " << endl;;

        double CGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int CGA_Iteration = 0;
        double CGA_Result = runCGA(XmlFile, RscAlcFile, CGA_SchTime, CGA_Iteration , EndFitness);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << CGA_Result << " " << CGA_SchTime << " " << CGA_Iteration << " " << endl;

        double TMGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int TMGA_Iteration = 0;
        double TMGA_Result = runTMGA(XmlFile, RscAlcFile, TMGA_SchTime, TMGA_Iteration , EndFitness);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << TMGA_Result  << " " << TMGA_SchTime  << " " << TMGA_Iteration  << " " << endl;

        outfile.close();
        //delete the first line in the file
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
    return 0;
}
