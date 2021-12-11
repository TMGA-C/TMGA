#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "TMGA.h"


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
    SchTime["Montage1000_0.7"] =1003.126 ;
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

    string strLine;
    ifstream iFile("../exp.txt");
    if (!iFile) {
        cout << "filelist open failed!\n";
        exit(1);
    }
    while(getline(iFile,strLine)){
        istringstream is(strLine);
        Orthogonal TemOrthogonal;
        is >> TemOrthogonal.SubPopSizeFac >> TemOrthogonal.NumOfSubPopM >> TemOrthogonal.ExchangeFactor
           >> TemOrthogonal.ElistRate >> TemOrthogonal.MutationRate >> TemOrthogonal.TerFacForSta1;
        orthogonal.push_back(TemOrthogonal);
    }
    iFile.close();

    string Model, NumOfTask, RscAvlRatio;
    do {
        string strLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(1);
        }
        getline(iFile, strLine);
        if (strLine.size() < 1) {
            cout << "Empty input file(fileList)" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(strLine);
        is >> Model >> NumOfTask >> RscAvlRatio;
        XmlFile = Model + "_" + NumOfTask + "_0.xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        int index =0 ;
        for(Orthogonal TemOrthogonal: orthogonal){
            ++index;
            ofstream outfile("../print/result"+to_string(index)+".txt", ios::app);
            if (!outfile) {
                cout << "Open the file failure...\n";
                exit(0);
            }
            outfile.setf(ios::fixed, ios::floatfield);
            outfile.precision(5);
            for (int times = 0; times < 5; ++times) {
                double TMGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
                int TMGA_Iteration = 0;
                double TMGA_Result = runTMGA(XmlFile, RscAlcFile, TemOrthogonal, TMGA_SchTime, TMGA_Iteration );
                outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                        << TMGA_Result << " " << TMGA_SchTime << " " << TMGA_Iteration
                        << endl;
            }
            outfile.close();
        }
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
}
