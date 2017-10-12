//    LSD - Least Square Dating for etimating substitution rate and divergence dates
//    Copyright (C) <2015> <Thu-Hien To>

//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include "readData.h"

int tree2data(FILE * tree,Pr* pr,Node** nodes,bool& constraintConsistent){
    int inode = 0;//number of internal nodes;
    stack<int> pileNode;
    char c = readBracket(tree,"input tree");
    int a=1;
    int countleaf=0;//n;
    int s=0;
    int nbChild=0;
    do{
        c = readChar(tree,"input tree");
        if (c==')'){
            a--;inode++;
            s=0;
            nbChild=0;
            nodes[pr->nbINodes-inode]=new Node();
            nodes[pr->nbINodes-inode]->L=readSupport(tree,"input tree");
            stack<int> listSuc;
            while (!pileNode.empty() && s!=-1) {
                s=pileNode.top();pileNode.pop();
                if (s!=-1){
                    nodes[s]->P=pr->nbINodes-inode;
                    listSuc.push(s);
                    nbChild++;
                }
            }
            while (!listSuc.empty()) {
                s=listSuc.top();
                listSuc.pop();
                nodes[pr->nbINodes-inode]->suc.push_back(s);
            }
            if (a>0) nodes[pr->nbINodes-inode]->B=readdouble(tree,"input tree");
            pileNode.push(pr->nbINodes-inode);
        }
        else if (c!='(' && c!=',' && c!=-1 && c!=';' && c!='\n'){
            string lb=readLabel(c,tree,a);
            pileNode.push(pr->nbBranches-countleaf);
            nodes[pr->nbBranches-countleaf]=new Node();
            nodes[pr->nbBranches-countleaf]->L=lb;
            nodes[pr->nbBranches-countleaf]->B=readdouble(tree,"input tree");
            countleaf++;
        }
        else if (c=='(') {a++;pileNode.push(-1);}
        else if (c=='\n') {
            c=readChar(tree,"input tree");
        }
    }  while (a>0);
    if (nbChild==2) {
        pr->rooted=true;
        nodes[0]->P=-1;
        nodes[0]->B=-1;
    }
    else{
        pr->rooted=false;
        nodes[1]->P=-1;
        nodes[1]->B=-1;
        nodes[0] = new Node();
    }
    int lineNb=getLineNumber(pr->inDateFile);
    FILE * dateFile = fopen(pr->inDateFile.c_str(),"rt");
    int ino=readInt(dateFile,"Error in the date file, the file should begin with an integer (the number of temporal constrains)");
    if (lineNb-1<ino) {
        cout<<"The number of given constraints is small than the number of constraints to read. Please change the number of constraints to read at the first line of the input date file."<<endl;
        exit(EXIT_FAILURE);
    }
    for (int i=0;i<ino;i++){
        string s=readWord(dateFile,pr->inDateFile);
        int type='n';
        double v1=0;
        double v2=0;
        int k = getPosition(nodes,s,0,pr->nbBranches+1);
        vector<int> mr;
        if (k==-1 && (s.compare("mrca")==0)){
            char c='(';
            while (c!=')'){
                string s="";
                c=readCommaBracket(dateFile,pr->inDateFile,s);
                int k1=getPosition(nodes,s,0,pr->nbBranches+1);
                if (k1!=-1){
                    mr.push_back(k1);
                }
            }
            if (mr.size()>0){ k=mrca(nodes,mr);}
        }
        if (k!=-1){
            if (nodes[k]->type!='n'){
                cout<<"Warning: There are nodes that have more than one temporal constraint"<<endl;
            }
            char c = readChar(dateFile,pr->inDateFile);
            while (c<33 || c>126) c=readChar(dateFile,pr->inDateFile);
            if (c=='l' || c=='L' || c=='u' || c=='U' || c=='b' || c=='B'){//interval value
                char p = readChar(dateFile,pr->inDateFile);
                if (p=='('){
                    if (c=='l' || c=='L'){
                        type='l';
                        v1=readdouble(dateFile,pr->inDateFile);
                    }
                    else if (c=='u' || c=='U'){
                        type='u';
                        v1=readdouble(dateFile,pr->inDateFile);
                    }
                    else if (c=='b' || c=='B'){
                        type='b';
                        v1=readdouble(dateFile,pr->inDateFile);
                        if (readChar(dateFile,pr->inDateFile)==','){
                            v2=readdouble(dateFile,pr->inDateFile);
                        }
                        else{
                            cout<<"date constraint of type 'b' must have two values"<<endl;
                            exit(EXIT_FAILURE);
                        }
                        if (v1>v2) {
                            double t=v1;
                            v1=v2;
                            v2=t;
                        }
                        if (v1==v2) {
                            type='p';
                        }
                        else {
                            type='b';
                        }
                    }
                    c=readChar(dateFile,pr->inDateFile);
                    while (c<33 || c>126) c=readChar(dateFile,pr->inDateFile);
                }
                else{
                    cout<<"Error reading "<<pr->inDateFile<<" file: calibration point must be defined as either 'l(lower_bound)' or 'u(upper_bound)' or 'b(lower_bound,upper_bound)'"<<endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {//precise value
                string wd="";
                wd+=c;
                while (fscanf(dateFile,"%c",&c)==1 && c>=33 && c<=126) {
                    wd+=c;
                }
                v1=atof(wd.c_str());
                type='p';
            }
            Date* newdate;
            if (mr.size()>0){
                newdate = new Date(type,v1,v2,mr);
            }
            else{
                newdate = new Date(type,v1,v2,k);
            }
            if (k<pr->nbINodes){
                pr->internalConstraints.push_back(newdate);
            }
            else{
                bool bl = nodes[k]->addConstraint(newdate);
                delete newdate;
                constraintConsistent=constraintConsistent && bl;
            }
        }
        else{
            cout<<"Warning: There are constraints on some taxa that are not in the input tree"<<endl;
            if (i<ino-1){
                char c=readChar(dateFile,pr->inDateFile);
                while (c!='\n') c=readChar(dateFile,pr->inDateFile);
                
            }
        }
    }
    fclose(dateFile);
    return s;
}


void readPartitionFile(Pr* pr){
    //Read partition file
    ifstream partFile(pr->partitionFile.c_str());
    string line;
    while (getline(partFile,line)) {
        int pos=0;
        string groupName = readWord(line,pos);
        int s = line.find_first_of("{",0);
        int e = line.find_first_of ("}",s+1);
        Part* part = new Part(groupName);
        while (s<e && s!=-1) {
            string term = line.substr(s+1,(e-s-1));
            int ss = 0;
            int ee = ss;
            Subtree* subtree;
            bool first = true;
            while (ee!=term.length()) {
                Pair* node;
                ss = firstCharacter(term,ss);
                ee = lastCharacter(term,ss);
                if (term.at(ss)=='[' && term.at(ee-1)==']') {
                    node = new Pair(true,term.substr(ss+1,ee-ss-2));
                }
                else {
                    node = new Pair(false,term.substr(ss,ee-ss));
                }
                if (first) {
                    subtree = new Subtree(node);
                    first=false;
                }
                else{
                    subtree->tips.push_back(node);
                }
                ss=ee+1;
            }
            part->subtrees.push_back(subtree);
            s = line.find_first_of("{",e+1);
            e = line.find_first_of("}",s+1);
        }
        pr->ratePartition.push_back(part);
    }
}

int tree2dataS(FILE * tree,Pr* pr,Node** nodes){
    int inode = 0;//number of internal nodes;
    stack<int> pileNode;
    char c = readBracket(tree,"input tree");
    int a=1;
    int countleaf=0;//n;
    int s=0;
    int nbChild=0;
    do{
        c = readChar(tree,"input tree");
        if (c==')'){
            a--;inode++;
            s=0;
            nbChild=0;
            nodes[pr->nbINodes-inode]=new Node();
            nodes[pr->nbINodes-inode]->L=readSupport(tree,"input tree");
            stack<int> listSuc;
            while (!pileNode.empty() && s!=-1) {
                s=pileNode.top();pileNode.pop();
                if (s!=-1){
                    nodes[s]->P=pr->nbINodes-inode;
                    listSuc.push(s);
                    nbChild++;
                }
            }
            while (!listSuc.empty()) {
                s=listSuc.top();
                listSuc.pop();
                nodes[pr->nbINodes-inode]->suc.push_back(s);
            }
            if (a>0) nodes[pr->nbINodes-inode]->B=readdouble(tree,"input tree");
            pileNode.push(pr->nbINodes-inode);
        }
        else if (c!='(' && c!=',' && c!=-1 && c!=';' && c!='\n'){
            string lb=readLabel(c,tree,a);
            pileNode.push(pr->nbBranches-countleaf);
            nodes[pr->nbBranches-countleaf]=new Node();
            nodes[pr->nbBranches-countleaf]->L=lb;
            nodes[pr->nbBranches-countleaf]->B=readdouble(tree,"input tree");
            countleaf++;
        }
        else if (c=='(') {a++;pileNode.push(-1);}
        else if (c=='\n') {
            c=readChar(tree,"input tree");
        }
    }  while (a>0);
    nodes[0]->P=-1;
    nodes[0]->B=-1;
    while (c!='\n') {
        c=readChar(tree,"input tree");
    }
    if (nbChild==2) {
        pr->rooted=true;
    }
    else{
        pr->rooted=false;
    }
    return s;
}

list<int> path(Pr* pr, Node** nodes,int s,int t){
    int cs=s;
    int ct=t;
    list<int> froms;
    stack<int> fromt;
    bool stops=false;
    bool stopt=false;
    while (!stops || !stopt){
        if (!stops) {
            froms.push_back(cs);
            cs=nodes[cs]->P;
            if (isAncestor(nodes, cs, t)) {
                stops=true;
            }
        }
        if (!stopt) {
            fromt.push(ct);
            ct=nodes[ct]->P;
            if (isAncestor(nodes, ct, s)) {
                stopt=true;
            }
        }
    }
    concat(froms,fromt);
    return froms;
}


int getBranchOut(Pr* pr,Node** nodes,list<string> &outgroups,bool &keepBelow){
    list<string>::iterator iter=outgroups.begin();
    if (outgroups.size()==1) {
        keepBelow=false;
        return getPosition(nodes,*iter,0,pr->nbBranches+1);
    }
    else{
        list<int> out;
        while (iter!=outgroups.end()) {
            int t=getPosition(nodes,*iter,0,pr->nbBranches+1);
            if (t!=-1) {
                out.push_back(t);
            }
            iter++;
        }
        if (out.size()>0) {
            int s=pr->nbINodes;
            while (s<=pr->nbBranches && contain(s , out)) {
                s++;
            }
            list<int>::iterator t=out.begin();
            list<int> common = path(pr,nodes,s,*t);
            t++;
            while (t!=out.end()) {
                list<int> p = path(pr,nodes,s,*t);
                common = intersect(common,p);
                t++;
            }
            int nbPass0=0;
            list<int>::iterator iter=common.begin();
            int r=*iter;
            while (iter!=common.end()){
                r=*iter;
                if (nodes[r]->P==0){
                    nbPass0++;
                }
                iter++;
            }
            if (nbPass0==2) {
                keepBelow=false;
            }
            else{
                keepBelow=true;
            }
            return r;
        }
        else{
            return -1;
        }
    }
}

void extrait_outgroup(Pr* pr,list<string> &outgroups){
    FILE * tree = fopen(pr->inFile.c_str(),"rt");
    if (tree==NULL) cout<<"Can not open the tree file"<<endl;
    else{
        Node** nodes = new Node*[pr->nbBranches+1];
        string newFile = pr->inFile;
        if (pr->keepOutgroup) newFile+=".reroot";
        else newFile+=".ingroup";
        FILE* w=fopen(newFile.c_str(),"wt");
        for (int y=0;y<pr->nbData;y++){
            printf("Analyzing outgroups of tree %d ...\n",y+1);
            for (int i=0; i<=pr->nbBranches; i++) {
                nodes[i]=new Node();
            }
            int s=tree2dataS(tree,pr,nodes);
            if (!pr->rooted) {
                nodes=unrooted2rootedS(pr, nodes, s);
            }
            bool keepBelow;
            int r=getBranchOut(pr, nodes, outgroups,keepBelow);
            if (r!=-1) {
                Node** nodes_new = cloneLeaves(pr,nodes,0);
                int p_r=reroot_rootedtree(r, pr, nodes, nodes_new);
                computeSuc_polytomy(pr, nodes_new);
                if (pr->keepOutgroup) {                    
                    fprintf(w,"%s",newick(0,0,pr,nodes_new).c_str());
                }
                else{
                    if (keepBelow) {
                        fprintf(w,"%s",newick(r, r, pr, nodes_new).c_str());
                    }
                    else{
                        fprintf(w,"%s",newick(p_r, p_r,pr, nodes_new).c_str());
                    }
                }
                delete[] nodes_new;
            }
            else{
                cout<<"The outgroups are not in the tree "<<y+1<<endl;
                fprintf(w,"%s",newick(0,0,pr,nodes).c_str());
            }
        }
        delete[] nodes;
        cout<<"The new input trees are writen to the file "<<newFile<<endl;
        fclose(w);
        pr->inFile=newFile;
    }
    fclose(tree);
}

