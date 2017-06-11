#ifndef PR_H
#define PR_H
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <vector>
#include "date.h"
#include "pair.h"
#include "part.h"
#include "subtree.h"

using namespace std;

typedef struct Pr
{
    string  inFile;         //Nom du fichier d'arbres
    string  inDateFile;     //Nom du fichier contenant les dates
    string partitionFile;   //Nom du fichier contenant la partition des arretes pour different taux
    string  outFile;        //Nom du fichier de resultats.
    string treeFile1;       //Nom du fichier d'abres sorties Nexus
    string treeFile2;       //Nom du fichier d'arbres sorties Newick avec des longueurs de branches mesures par substiution par site
    string treeFile3;       //Nom du fichier d'arbres sorties Newick avec des longueurs de branches mesures par annees
    bool relative;         //=true if all the leaves have the same date, the program estimate the relative dates
    double mrca;
    double leaves;
    int    seqLength;      //Longueur des sequences dans l'alignement
    int    nbData;         //Nombre de cas a  traiter (dans le cas de bootstrap)
    string fnOutgroup;
    string rate;           //le fichier contient les taux en entree
    string estimate_root;    //Method to estimate root
    bool rooted;
    bool constraint;       //Impose the constraints or not
    int variance;         //Use the variances or not
    bool ci;         //Compute confidence interval or not
    int  c;                //var = b+c/s;
    double rho_min;
    int nbINodes;
    int nbBranches;
    double rho;
    double* multiplierRate;
    bool givenRate;
    double objective;
    int nbSampling;
    bool keepOutgroup;
    vector<Part* > ratePartition;
    vector<Date*> internalConstraints;
    Pr(int n,int m){
        nbINodes=n;
        nbBranches=m;
    }
    void initConstraints(){
        internalConstraints.clear();
    }
    void copy(Pr* pr){
        relative=pr->relative;         //=true if all the leaves have the same date, the program estimate the relative dates
        mrca=pr->mrca;
        leaves=pr->leaves;
        seqLength=pr->seqLength;      //Longueur des sequences dans l'alignement
        nbData=pr->nbData;         //Nombre de cas a  traiter (dans le cas de bootstrap)
        rate=pr->rate;           //le fichier contient les taux en entree
        rooted=pr->rooted;
        constraint=pr->constraint;       //Impose the constraints or not
        variance=pr->variance;         //Use the variances or not
        ci=pr->ci;         //Compute confidence interval or not
        c=pr->c;                //var = b+c/s;
        rho=pr->rho;
        rho_min=pr->rho_min;
        givenRate=pr->givenRate;
        nbSampling=pr->nbSampling;
        ratePartition = pr->ratePartition;
        multiplierRate = new double[pr->ratePartition.size()+1];
        if (pr->ratePartition.size() > 0) {
            multiplierRate[0] = pr->multiplierRate[0];
        }
    }
    Pr()
    {
        inFile = "";
        inDateFile = "";
        partitionFile = "";
        outFile = "";
        treeFile1 = "";
        treeFile2 = "";
        treeFile3 = "";
        fnOutgroup = "";
        seqLength = 1000;
        nbData = 1;
        rate = "";
        relative = false;
        mrca=0;
        leaves=1;
        estimate_root = "";
        constraint = false;
        variance = 0;
        c = 10;
        rho_min = 1e-10;
        ci = false;
        nbSampling=100;
        rooted=true;
        keepOutgroup=false;
        givenRate=false;
        ratePartition = vector<Part* >();
    }
}Pr;
#endif


