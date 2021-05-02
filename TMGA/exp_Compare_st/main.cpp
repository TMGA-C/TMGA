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

    map<string, double> EndFit;
    EndFit["Montage25_0.4"] =  42.42 ;
    EndFit["Montage25_0.7"] =  29.88 ;
    EndFit["Montage25_1.0"] =  28.92 ;
    EndFit["Montage100_0.4"] = 149.32 ;
    EndFit["Montage100_0.7"] = 142.85 ;
    EndFit["Montage100_1.0"] = 104.35 ;
    EndFit["Montage1000_0.4"] =1075.49 ;
    EndFit["Montage1000_0.7"] =1037.27 ;
    EndFit["Montage1000_1.0"] =1036.25 ;

    EndFit["Epigenomics24_0.4"] = 3179.97 ;
    EndFit["Epigenomics24_0.7"] = 2321.04 ;
    EndFit["Epigenomics24_1.0"] = 2320.91 ;
    EndFit["Epigenomics100_0.4"] =47104.73 ;
    EndFit["Epigenomics100_0.7"] =36627.68 ;
    EndFit["Epigenomics100_1.0"] =36108.50 ;
    EndFit["Epigenomics997_0.4"] =330492.48 ;
    EndFit["Epigenomics997_0.7"] =325458.56 ;
    EndFit["Epigenomics997_1.0"] =325041.26 ;

    EndFit["Ligo30_0.4"] =  995.50 ;
    EndFit["Ligo30_0.7"] =  692.69 ;
    EndFit["Ligo30_1.0"] =  664.28 ;
    EndFit["Ligo100_0.4"] = 2024.36 ;
    EndFit["Ligo100_0.7"] = 1863.73 ;
    EndFit["Ligo100_1.0"] = 1774.48 ;
    EndFit["Ligo1000_0.4"] =19207.31 ;
    EndFit["Ligo1000_0.7"] =19053.96 ;
    EndFit["Ligo1000_1.0"] =19039.35 ;


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

//        double HEFT_SchTime  = 0;
//        double HEFT_Result = runHEFT(XmlFile, RscAlcFile,HEFT_SchTime);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << HEFT_Result << " " << HEFT_SchTime << " "  << endl;

        double EndFitness = EndFit[Model + NumOfTask + "_" + RscAvlRatio];

//        double HGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int HGA_Iteration = 0;
//        double HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration , EndFitness);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " " << endl;

//        double NGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int NGA_Iteration = 0;
//        double NGA_Result = runNGA(XmlFile, RscAlcFile, NGA_SchTime, NGA_Iteration , EndFitness);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << NGA_Result << " " << NGA_SchTime << " " << NGA_Iteration << " " << endl;

//        double LWSGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int LWSGA_Iteration = 0;
//        double LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration , EndFitness);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " " << endl;;

//        double CGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int CGA_Iteration = 0;
//        double CGA_Result = runCGA(XmlFile, RscAlcFile, CGA_SchTime, CGA_Iteration , EndFitness);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << CGA_Result << " " << CGA_SchTime << " " << CGA_Iteration << " " << endl;

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
