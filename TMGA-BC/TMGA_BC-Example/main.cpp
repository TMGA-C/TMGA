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
#include "Exhaustive.h"
#include "IHEFT3.h"

using namespace std;

int main() {
    srand((int) time(0));
    //{set the runtime (termination) time of the algorithm}

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
//
//        double IHEFT3_SchTime  = 0;
//        double IHEFT3_Result = runIHEFT3(XmlFile, RscAlcFile,IHEFT3_SchTime);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << IHEFT3_Result << " " << IHEFT3_SchTime << " "  << endl;

//        double Ex_SchTime  = 0;
//        chromosome Ex_Chrom;
//        Ex_Chrom = runExhaustive(XmlFile, RscAlcFile,Ex_SchTime);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << Ex_Chrom.FitnessValue << " " << Ex_SchTime << " "  << endl;
//        for (int i = 0; i < comConst.NumOfTsk; ++i) {
//            outfile<<Ex_Chrom.TskSchLst[i]<<" ";
//        }
//        outfile<<";";
//        for (int i = 0; i < comConst.NumOfTsk; ++i) {
//            outfile<<Ex_Chrom.RscAlcLst[i]<<" ";
//        }

//        double EndFitness = 77;
//        double HGA_SchTime  = 0;
//        int HGA_Iteration = 0;
//        double HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration , EndFitness);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " " << endl;

//        double EndFitness = 82;
//        double NGA_SchTime  = 0;
//        int NGA_Iteration = 0;
//        double NGA_Result = runNGA(XmlFile, RscAlcFile, NGA_SchTime, NGA_Iteration , EndFitness);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << NGA_Result << " " << NGA_SchTime << " " << NGA_Iteration << " " << endl;
//
//        double EndFitness = 77;
//        double LWSGA_SchTime  = 0;
//        int LWSGA_Iteration = 0;
//        double LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration , EndFitness);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " " << endl;;
//
//
//        double EndFitness = 73;
//        double CGA_SchTime  = 0;
//        int CGA_Iteration = 0;
//        double CGA_Result = runCGA(XmlFile, RscAlcFile, CGA_SchTime, CGA_Iteration , EndFitness);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << CGA_Result << " " << CGA_SchTime << " " << CGA_Iteration << " " << endl;
//
//
        double EndFitness = 73;
        double TMGA_SchTime  = 0;
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
