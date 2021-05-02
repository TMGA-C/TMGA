#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "TMGA.h"


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
