#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "NGA.h"
#include "CGA.h"
#include "LWSGA.h"
#include "HGA.h"

using namespace std;

int main() {
    srand((int) time(0));
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

        double HGA_SchTime  = 0;
        int HGA_Iteration = 0;
        double HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration);

        double NGA_SchTime  = 0;
        int NGA_Iteration = 0;
        double NGA_Result = runNGA(XmlFile, RscAlcFile, NGA_SchTime, NGA_Iteration);

        double LWSGA_SchTime  = 0;
        int LWSGA_Iteration = 0;
        double LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration);

        double CGA_SchTime  = 0;
        int CGA_Iteration = 0;
        double CGA_Result = runCGA(XmlFile, RscAlcFile, CGA_SchTime, CGA_Iteration);

        //results are written into the file
        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(5);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " "
                << NGA_Result << " " << NGA_SchTime << " " << NGA_Iteration << " "
                << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " "
                << CGA_Result << " " << CGA_SchTime << " " << CGA_Iteration << " "
                << endl;
        outfile.close();
        //delete the first line in the file
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
    return 0;
}
