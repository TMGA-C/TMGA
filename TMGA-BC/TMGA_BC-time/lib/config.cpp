#include <fstream>
#include <sstream>
#include "config.h"
#include "pugixml.hpp"

void DeleteFirstLineInFile(string fileName) {
    vector<string> VecContent;
    string StrLine;
    ifstream iFile(fileName);
    if (!iFile) {
        cout << fileName << "filelist open failed!\n";
        exit(0);
    }
    //Read all content in the document into "VecContent"
    while (iFile) {
        getline(iFile, StrLine);
        VecContent.push_back(StrLine);
    }
    iFile.close();
    VecContent.erase(VecContent.begin()); // delete the first line
    ofstream oFile(fileName);
    if (!oFile) {
        cout << fileName << "filelist open failed!\n";
        exit(0);
    }
    vector<string>::const_iterator iter = VecContent.begin();
    //{Rewrite the contents of "vecContent" into the file}.
    for (; VecContent.end() != iter; ++iter) {
        oFile.write((*iter).c_str(), (*iter).size());
        oFile << '\n';
    }
    oFile.close();
}

int ReadID(string id) {
    int ret = 0;
    for (int i = 0; i < id.length(); ++i) {
        if (id[i] <= '9' && id[i] >= '0')
            ret = ret * 10 + id[i] - '0';
    }
    return ret;
}

//｛Read the model information from the file and store them in the data structures}
// In "Epigenomics_997" there are following unreasonable cases:
// 1) the runtime of some tasks is less than 0 or equal to 0;
// 2) the size of some files is less than 0.
// So we let the length of task be equal to max{fabs(runtime),0.001} and let the size of file be equal to fabs(size)
void ReadFile(string XmlFile, string RscAlcFile) {
    string FilePath = "../data";
    string XmlPath = FilePath + "/" + XmlFile;
    pugi::xml_document doc;
    int w = doc.load_file((XmlPath).c_str());
    pugi::xml_node root = doc.child("adag");
    for (pugi::xml_node job = root.child("job"); job; job = job.next_sibling("job")) {
        Task task = Task();
        task.length = XY_MAX(fabs(atof(job.attribute("runtime").value())),0.001);// read the length of task
        //{read file}
        for (pugi::xml_node uses = job.child("uses"); uses; uses = uses.next_sibling("uses")) {
            vfile file = vfile();
            file.source = -1;
            if (!strcmp(uses.attribute("link").value(), "input")) {    //read input file
                file.FileName = uses.attribute("file").value();
                file.size = fabs(atof(uses.attribute("size").value()));
                task.IFile.push_back(file);
            } else { //read output file
                file.FileName = uses.attribute("file").value();
                file.size = fabs(atof(uses.attribute("size").value()));
                task.OFile.push_back(file);
            }
        }
        Tasks.push_back(task);
    }
    comConst.NumOfTsk = Tasks.size();

    // read the info of relation between tasks
    for (pugi::xml_node child = root.child("child"); child; child = child.next_sibling("child")) {
        int ChildIndex = ReadID(child.attribute("ref").value());
        for (pugi::xml_node parent = child.child("parent"); parent; parent = parent.next_sibling("parent")) {
            int ParentIndex = ReadID(parent.attribute("ref").value());
            Tasks[ChildIndex].parents.push_back(ParentIndex);
            Tasks[ParentIndex].children.push_back(ChildIndex);
        }
    }

    //{calculate the transfer data size among tasks}
    ParChildTranFileSizeSum.resize(comConst.NumOfTsk);
    for (int k = 0; k < comConst.NumOfTsk; ++k) {
        ParChildTranFileSizeSum[k].resize(comConst.NumOfTsk,0);
    }
    for (int i = 0; i < Tasks.size(); ++i) {
        if (Tasks[i].parents.size() == 0) continue;
        for (int p = 0; p < Tasks[i].IFile.size(); ++p) {               //two loop (for p and j) can be switched in order
            string IName = Tasks[i].IFile[p].FileName;
            for (int j = 0; j < Tasks[i].parents.size(); ++j) {         //Traverse the parent task
                int Parent = Tasks[i].parents[j];
                for (int q = 0; q < Tasks[Parent].OFile.size(); ++q) {  //Traverse the output files of the parent task
                    string OName = Tasks[Parent].OFile[q].FileName;
                    if (IName.compare(OName) == 0) {                    // judge whether two file names are the same; 0: same; -1: not same
                        ParChildTranFileSizeSum[Parent][i] += Tasks[i].IFile[p].size;
                        //If multiple identical files from different parent tasks are transferred to the same child task, the "source" records the last parent task
                        Tasks[i].IFile[p].source = Parent;
                        break;
                    }
                }
            }
        }
    }

    //{Rsc can be added here}
    Resource Rsc_0 = Resource(0, 1, 20);//-w
    Resource Rsc_1 = Resource(0, 1, 20);
    Resource Rsc_2 = Resource(0, 2, 30);
    Resource Rsc_3 = Resource(1, 2, 30);
    Resource Rsc_4 = Resource(1, 3, 40);
    Resource Rsc_5 = Resource(1, 3, 40);
    Rscs.push_back(Rsc_0);
    Rscs.push_back(Rsc_1);
    Rscs.push_back(Rsc_2);
    Rscs.push_back(Rsc_3);
    Rscs.push_back(Rsc_4);
    Rscs.push_back(Rsc_5);
    comConst.NumOfRsc = Rscs.size();

    //read the RscAlc file to task data structure
    //in the RscAlc file, each resource can perform at least one task and each task can be performed by at least one resource
    char line[4096] = {0};
    int TskIndex = 0;
    int RscIndex = 0;
    string RscAlcPath = FilePath + "/" + RscAlcFile;  //RscAlcPath xy4
    ifstream fin(RscAlcPath, ios::in);
    if (!fin) {
        cout << "Error at open Rsc file" << endl;
        exit(0);
    } else {
        while (fin.getline(line, sizeof(line))) {
            stringstream Word(line);
            while (1) {
                Word >> TskIndex;
                if (Word.fail()) break;
                Tasks[TskIndex].ElgRsc.push_back(RscIndex);
                Rscs[RscIndex].ElgTsk.push_back(TskIndex);
            }
            ++RscIndex;
        }
    }
    ModelScale = 0;  //-xy4
    for (int i = 0; i < comConst.NumOfTsk; ++i ){
        ModelScale += Tasks[i].ElgRsc.size();
    }
    fin.close();
}

void ClearALL() {
    Tasks.clear();
    TskLstInLvl.clear();
    LevelIdOfTask.clear();
    Rscs.clear();
    population.clear();
    populations.clear();
    ParChildTranFileSizeSum.clear();
}

void ConfigParameter_CGA() {
    Parameter_CGA.NumOfChormPerPop = 2 * Tasks.size();
    Parameter_CGA.MutationRate = 0.2;
    Parameter_CGA.CrossoverRate = 0.8;
}

void ConfigParameter_HGA() {
    Parameter_HGA.NumOfChormPerPop = 100;
    Parameter_HGA.MutationRate = 0.02;
    Parameter_HGA.EliteRate = 0.2;
}

void ConfigParameter_LWSGA() {
    Parameter_LWSGA.NumOfChormPerPop = 70;
    Parameter_LWSGA.CrossoverRate = 0.7;
}

void ConfigParameter_NGA() {
    Parameter_NGA.NumOfChormPerPop = 4 * Tasks.size();
    Parameter_NGA.MutationRate = 0.2;
    Parameter_NGA.CrossoverRate = 0.8;
    Parameter_NGA.EliteRate = 0.2;
}

void ConfigParameter_TMGA() {
    Parameter_TMGA.NumOfChrInSubPop = Tasks.size()*2;                                   //set the subpopulation size
    if (Parameter_TMGA.NumOfChrInSubPop % 2 == 1) {
        ++Parameter_TMGA.NumOfChrInSubPop;
    }
    Parameter_TMGA.NumOfSubPop = 3;                                                     //set the number of subpopulations
    Parameter_TMGA.interval = 8;                                                       // set the exchange interval/0.8
    Parameter_TMGA.NumOfEliteOfPop = ceil(Parameter_TMGA.NumOfChrInSubPop * 0.5);   //set the number of elitist for exchange -xy
    Parameter_TMGA.MutationRate_TMGA = 0.125;                                            //set the mutation rate
    Parameter_TMGA.TrmThresholdOfStg1 = ceil(10 * ModelScale / Parameter_TMGA.NumOfChrInSubPop / Parameter_TMGA.NumOfSubPop);  //set the termination condition of stage 1
}




