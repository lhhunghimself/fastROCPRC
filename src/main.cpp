#include <string>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <tclap/CmdLine.h>
#include "my_sort.hpp"
#include "timers.hpp"

using namespace std;

void readGoldEdges(string goldListFile,std::unordered_set<string> *allRegulators, std::unordered_set<string> *allTargets,std::unordered_set<string> *goldNodePairs,char separator, bool threeColumns,bool noSelf,bool filtered);
void readSampleEdges(string edgeListFile,std::unordered_set<string> *allRegulators, std::unordered_set<string> *allTargets, vector<std::string> *sampleNodePairs, vector<double> *sampleProbs, char separator, bool noSelf, bool filter);
template <class T> int 	contabs(int nSamples,int nPositives,int totalEdges,vector<T> sampleProbs,std::vector<std::string> sampleSet,std::unordered_set<string> goldSet,int *TPs,int *FPs,T *TPR, T *FPR,T* confidences,int *order,vector<T> threshholds, vector<int> &outputIndices);
template <class T> void interpolatePR (int n,int nPositives,int totalPoints, int *TP, int *FP,T *precision,T *recall, bool *inter);
template <class T> T AUROC(int n,T *x,T *y);
template <class T> T AUPRC(int n,T *x,T *y);
template <class T> T AUC(int n,T *x, T *y);

template <class T> T AUCblock(int n,T *x, T *y);
 
int main(int argc, char *argv[]){    
  string edgeListFile="";
  string goldListFile="";
  string contabsFile;
  std::unordered_set<string>  allRegulators,allTargets,goldNodePairs;
  std::vector<std::string> sampleNodePairs;
  std::vector<double> sampleProbs;
  std::vector<double>threshholds;
  
  int totalEdges=0;
  bool timer=0,verbose=0,threeColumns=0,noSelf=0,filter=0,outputPR=0,AUCOnly=0;
  char separator=',';
  float increment=0.0000001;
  
		struct timespec start_time,read_time,end_time;
		
		try {  
   using namespace TCLAP;
	  CmdLine cmd("fastROCPRC flags", ' ', "1.00");
	  ValueArg<string> edgeListFileArg("e","edgeList","edgeList file to be evaluated",true,"","string");
	  ValueArg<string> goldListFileArg ("g","goldList","two or three column list of edges",true,"","string");
	  ValueArg<string> contabsFileArg ("c","contabs","file to store the contingency tables",false,"","string");
	  ValueArg<char> separatorArg("s","separator","internal separator - must be a character not found in names","false",',',"char");
	  ValueArg<int> totalEdgesArg ("t","totalEdges","Total number of edges in network",false,0,"int");
	  MultiArg <double> threshholdsArg("","threshholds","cutoff values to evaluate Precision Recall FPR\n",false,"float");
	  SwitchArg filterArg ("f","filter","filter out sample edge list regulators and targets not found as a regulator in the gold set",cmd,false);
	  SwitchArg timerArg ("","time","shows execution time in seconds not including file i/o ",cmd,false);
	  SwitchArg verboseArg("","verbose","show more parameter information",cmd,false);
	  SwitchArg threeColumnsArg("","threeColumns","gold standard has 3 columns instead of two",cmd,false);
	  SwitchArg noSelfArg("","noSelf","ignore self targets",cmd,false);
	  SwitchArg AUCOnlyArg("a","AUCOnly","only output AUROC and AUPRC in that order",cmd,false);

	  cmd.add(edgeListFileArg);
	  cmd.add(goldListFileArg);
	  cmd.add(totalEdgesArg);
	  cmd.add(contabsFileArg);
	  cmd.add(threshholdsArg);
	  cmd.parse( argc, argv );
	  // Get the value parsed by each arg
	  edgeListFile= edgeListFileArg.getValue();
	  goldListFile= goldListFileArg.getValue();
	  contabsFile= contabsFileArg.getValue();
   totalEdges= totalEdgesArg.getValue();
   verbose=verboseArg.getValue();
   separator=separatorArg.getValue();
   timer=timerArg.getValue();
			threeColumns=threeColumnsArg.getValue();
			noSelf=noSelfArg.getValue();
			threshholds=threshholdsArg.getValue();
			AUCOnly=AUCOnlyArg.getValue();
			filter=filterArg.getValue();
			if(!threshholds.size()){
				threshholds.push_back(.5);
				threshholds.push_back(.95);
			}

	 // Do what you intend. 
	  if(verbose){
	   cerr << "Input Parameters: "<< endl;
		  cerr << "edgeListFile: " << edgeListFile << endl;		 
	   cerr << "goldListFile: " << goldListFile << endl;
	   cerr << "Total edges in network: " << totalEdges << endl;
    cerr << "contabs file: " << contabsFile << endl;
    cerr << "separator: " << separator <<endl;
    cerr << "totalEdges: " << totalEdges <<endl;
    cerr << "threshholds: ";
    for (int i=0;i<threshholds.size();i++) cerr<< threshholds[i] << " ";
    cerr << endl;
    cerr << "timer: " << timer << endl;
    cerr << "threeColoumns: " << threeColumns <<endl;
    cerr << "noSelf: " << noSelf <<endl;
    cerr << "filter: " << filter <<endl;
    cerr << "AUCOnly: " << AUCOnly << endl;
			}		 	
	} 
	catch (TCLAP::ArgException &e){ cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
	
	current_utc_time(&start_time);
  
	readGoldEdges(goldListFile,&allRegulators,&allTargets,&goldNodePairs,separator,threeColumns,noSelf,filter);
	readSampleEdges(edgeListFile,&allRegulators,&allTargets,&sampleNodePairs,&sampleProbs,separator,noSelf,filter);
	current_utc_time(&read_time);
	if(timer) cerr << "read time: "<< get_elapsed_time(&start_time, &read_time) << " seconds"<<endl;
	
	int nPositives=goldNodePairs.size(); 
	int nSampleEdges=sampleNodePairs.size();
	int *order= new int [nSampleEdges];
	double pr[nPositives+1],recall[nPositives+1];
	double *TPR = new double [nSampleEdges];
	double *FPR = new double [nSampleEdges];
	double *confidences= new double [nSampleEdges];
	bool inter[nPositives+1]; //indicates whether the point is interpolated or not
	memset(inter,0,(nPositives+1)*sizeof(bool));
	sort_by_scores (nSampleEdges,&sampleProbs[0],order,0);
	if(!totalEdges){
		if(noSelf)
		 totalEdges=allRegulators.size()*allRegulators.size()-allRegulators.size();
		else totalEdges=allRegulators.size()*allRegulators.size();
	}
	fprintf(stderr,"%d sampleEdges read %d totalEdges\n",nSampleEdges,totalEdges);	
 int *TPs = new int [nSampleEdges];
 int *FPs = new int [nSampleEdges];

 //make contingency tables

	vector<int>outputIndices(threshholds.size());

	int nPoints=contabs(nSampleEdges,nPositives,totalEdges,sampleProbs,sampleNodePairs,goldNodePairs,TPs,FPs,TPR,FPR,confidences,order,threshholds,outputIndices);
	
	//print out the contingency table precision recall FPR for the threshholds
	if(!AUCOnly){
	 fprintf(stdout,"Threshhold\tTP\tFP\tTN\tFN\tPrecision\tRecall\tFPR\n");
	 for(int i=0;i<threshholds.size();i++){
	 	int k=outputIndices[i];
	 	if(k<0){
	 		//threshhold not met	 
	   int TP=0,FP=0,TN=totalEdges-nPositives,FN=nPositives;
    double precision=0,recall=0,FPRate=0;  
	 		fprintf(stdout,"%f\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n",threshholds[i],TP,FP,TN,FN,precision,recall,FPRate);
	 	}
	 	else{	
	 	 int TN=totalEdges-nPositives-FPs[k];
	 	 int FN=(nPositives-TPs[k]);
		  double precision=TPs[k]/(double)(TPs[k]+FPs[k]);
		  double recall=TPs[k]/(double)nPositives;
		  fprintf(stdout,"%f\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n",threshholds[i],TPs[k],FPs[k],TN,FN,precision,recall,FPR[k]);
		 }
	 }
	}	 
	interpolatePR (nPoints,nPositives,totalEdges,TPs,FPs,pr,recall,inter);
	fprintf(stdout,"\nAUROC\tAUPRC\n%f\t%f\n",AUROC(nPoints,FPR,TPR),AUCblock(nPositives+1,recall,pr));
	current_utc_time(&end_time);
 if(!AUCOnly){
	 FILE *outfp=0;
	 if(contabsFile== ""){
	  outfp=stdout;
	  fprintf(stdout,"\n");	
	 }
	 else{
	 		//print out contabs
	  fprintf(stderr,"opening %s to write contingency tables\n",contabsFile.c_str());
	 	outfp=fopen(contabsFile.c_str(),"w");
	 }
	 if(outfp){
	  fprintf(outfp,"Threshhold\tTP\tFP\tTN\tFN\tPrecision\tRecall\tFPR\tInter\n");
	  int k=0;
	  double lastConf=1.0;
	  for (int i=0;i<nPositives+1;i++){
	  if(inter[i]){
			 double TP=(double)i;
			 double FP= (i)? (TP-pr[i]*TP)/pr[i] : 0;
			 double TN=((double)(totalEdges-nPositives))-FP;
	   double FN=((double)nPositives)-TP;
	   fprintf(outfp,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",lastConf,TP,FP,TN,FN,pr[i],TP/(double)nPositives,FP/(double)(FP+TN),(int)inter[i]);
			}
			else{
			 int TP=i,FP=FPs[k];
	   int TN=totalEdges-nPositives-FP;
	   int FN=(nPositives-TP);
	   fprintf(outfp,"%f\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\n",confidences[k],TP,FP,TN,FN,pr[i],TP/(double)nPositives,FP/(double)(FP+TN),(int)inter[i]);
	   lastConf=confidences[k];
	   k++;
		 	}
	  }
	  if(outfp != stdout)fclose(outfp);
	 }
	}	
	if(timer) cerr << "execution time: "<< get_elapsed_time(&read_time, &end_time) << " seconds"<<endl;
	if(timer) cerr << "elapsed time: "<< get_elapsed_time(&start_time, &end_time) << " seconds"<<endl;	 
	delete [] TPs;
	delete [] FPs;
	delete [] order;
	delete []TPR;
	delete []FPR;
	delete []confidences;
	return(1);
}
void readGoldEdges(string goldListFile,std::unordered_set<string> *allRegulators, std::unordered_set<string> *allTargets,std::unordered_set<string> *goldNodePairs,char separator, bool threeColumns, bool noSelf, bool filtered){
	char line[2048];
	fprintf(stderr,"opening %s\n",goldListFile.c_str());
	vector<pair<string,string>> unfilteredList;
	FILE *fp=fopen(goldListFile.c_str(),"r");
	while(fgets(line,2048,fp)){
		char g1[1024],g2[1024];
		float value=-1;
		if(threeColumns){
		 sscanf(line,"%s %s %f",g1,g2,&value);
		}
		else{
			sscanf(line,"%s %s",g1,g2);
		}
		if(filtered){
			if(!noSelf || strcmp(g1,g2)){
				string str1(g1),str2(g2);
			 allRegulators->insert(string(g1));
		  allTargets->insert(string(g2));
		  unfilteredList.push_back(make_pair(str1,str2));
			}
		}
		else if(!noSelf || strcmp(g1,g2)){	
		 allRegulators->insert(string(g1));
		 allTargets->insert(string(g2));
		 if(!threeColumns || value > 0){
		  int len=strlen(g1);
		  g1[len]=separator;
		  g1[len+1]='\0';
		  goldNodePairs->insert(string(strcat(g1,g2)));
		 }
		}
	}
	fclose(fp);
	for(int i=0;i<unfilteredList.size();i++){
		string str1=unfilteredList[i].first;
		string str2=unfilteredList[i].second;
		if(allRegulators->count(str1) && allRegulators->count(str2)){
			str1+=separator;
		 goldNodePairs->insert(str1+str2); 	
		}
	}	
}		
void readSampleEdges(string edgeListFile,std::unordered_set<string> *allRegulators, std::unordered_set<string> *allTargets, vector<std::string> *sampleNodePairs, vector<double> *sampleProbs, char separator, bool noSelf, bool filter){
	char line[2048];
	fprintf(stderr,"opening %s\n",edgeListFile.c_str());
	FILE *fp=fopen(edgeListFile.c_str(),"r");
	while(fgets(line,2048,fp)){
		float value=0;
		char g1[1024],g2[1024];
		sscanf(line,"%s %s %f",g1,g2,&value);
		if(!noSelf || strcmp(g1,g2)){	
			if(filter){
				if(allRegulators->count(string(g1)) && allRegulators->count(string(g2))){
				//add to sample if the Regulator and  are in the gold set
				 int len=strlen(g1);
		   g1[len]=separator;
		   g1[len+1]='\0';
		   sampleNodePairs->push_back(string(strcat(g1,g2)));
		   sampleProbs->push_back(value);
				}
			}
			else{
				//assume that the prediction universe is a union of the gold and sample set
		  allRegulators->insert(string(g1));
		  allTargets->insert(string(g2));
		  int len=strlen(g1);
		  g1[len]=separator;
		  g1[len+1]='\0';
		  sampleNodePairs->push_back(string(strcat(g1,g2)));
		  sampleProbs->push_back(value);
			}
		}
	}
	fclose(fp);
}
template <class T> T AUROC(int n,T *x,T *y){
	//starts at 0 goes to 1
 T auc=AUCblock(n,x,y);
	//check edges
	if (x[0] > 0) auc+=x[0]*y[0]*0.5; 
 if (x[n-1] < 1) auc+=(1.0-x[n-1])*(1.0+y[n-1])*0.5;
 return(auc);
}	
template <class T> void interpolatePR (int n,int nPositives,int totalPoints, int *TP, int *FP,T *precision,T *recall, bool *inter){
	//write out interpolated FPs 
	int nExt=n;
	int TPext[n+2];
	int FPext[n+2];

	if(FP[0] !=0){
		precision[0]=0;
		recall[0]=0;
		TPext[0]=0;FPext[0]=0;
		memmove(TPext+1,TP,n*sizeof(int));
		memmove(FPext+1,FP,n*sizeof(int));
		inter[0]=1;
		nExt++;
	}
	else{
		memmove(TPext,TP,n*sizeof(int));
		memmove(FPext,FP,n*sizeof(int));		
	}	
	if(TP[nExt-1] !=nPositives){
		TPext[nExt]=nPositives;
		FPext[nExt]=totalPoints-nPositives;
		nExt++;
		inter[nPositives+1]=1;
	}
	//check for every true positive
	int lastTP=TPext[0];
	int lastFP=FPext[0];
	for(int i=1;i<nExt;i++){
		if(TPext[i] > lastTP){
			const T skew=(FPext[i]-lastFP)/(T)(TPext[i]-lastTP);
			for(int j=lastTP+1;j<TPext[i];j++){
				recall[j]=j/(T)nPositives;
				precision[j]=j/(T)(j+lastFP+skew*(j-lastTP));
				inter[j]=1;
			}
			precision[TPext[i]]=	(T )TPext[i]/(T)(TPext[i]+FPext[i]);
			recall[TPext[i]]=TPext[i]/(T)nPositives;
		 //fprintf(stderr,"*%d %d %f %f %f\n",TPext[i],FPext[i],skew,recall[TPext[i]],precision[TPext[i]]);
		}
		lastTP=TPext[i];
		lastFP=FPext[i];
	}
	recall[nPositives]=1;
 precision[nPositives]=nPositives/(T)totalPoints;
}
  

template <class T> T AUCtrap(int n,T *x, T *y){
	//trapezoidal area
	double auc=0;
	T xdiff,ymean;
	for(int i=1;i<n;i++){
		xdiff=x[i]-x[i-1];
		if(xdiff){
		 ymean=(y[i]+y[i-1])*0.5;
		 auc+=ymean*xdiff;
		}			 
	}
	return((T) auc);
}
template <class T> T AUCblock(int n,T *x, T *y){
	//min rectangular area
	double auc=0;
	T xdiff,ymean;
	for(int i=1;i<n;i++){
		xdiff=x[i]-x[i-1];
		if(xdiff){
			const T ymin=(y[i]<y[i-1])? y[i]:y[i-1];
		 auc+=ymin*xdiff;
		}			 
	}
	return((T) auc);
}	
template <class T> int 	contabs(int nSamples,int nPositives,int totalEdges,vector<T> sampleProbs,std::vector<std::string> sampleSet,std::unordered_set<string> goldSet,int *TPs,int *FPs,T *TPR, T *FPR, T *confidences,int *order,vector<T> threshholds, vector<int> &outputIndices){
	//order by largest threshhold first
	//returns contingency tables
	//outputIndices are where the value of  
	T oldProbs=sampleProbs[order[0]],probs=0;
	int nMatches=0; //number of contabs found matching threshold values;
	int tvectorSize=threshholds.size(); //size of the threshhold vector
	for (int i=0;i<tvectorSize;i++) outputIndices[i]=-1;
	int TP=0,FP=0,TN=totalEdges-nPositives,FN=nPositives;
	int k=0;
	for(int i=0;i<nSamples;i++){
		const int m=order[i];
  probs=sampleProbs[m];
		if(probs != oldProbs){
			double precision=TP/(double)(TP+FP);
			TPR[k]=(double)TP/(TP+FN);
			FPR[k]=FP/(double)(TN+FP);
			TPs[k]=TP;
			FPs[k]=FP;
			confidences[k]=oldProbs;
			oldProbs=probs;
			//check output indices
			if(nMatches < tvectorSize){
			 for(int j=0;j<tvectorSize;j++){
			  if(outputIndices[j] <0 && probs < threshholds[j]){
						outputIndices[j]=k;
						nMatches++;
				 }
			 }	 
		 }
			k++;
		}
		if(goldSet.count(sampleSet[m])){
			//positive positive prediction
			TP++;
			FN--;
		}
		else{
			//negative positive prediction
			FP++;
			TN--;
		}		
	}
	double precision=TP/(double)(TP+FP);
	TPR[k]=(double)TP/(TP+FN);
	FPR[k]=FP/(double)(TN+FP);
	TPs[k]=TP;
	FPs[k]=FP;
	if(nMatches < tvectorSize){
		for(int j=0;j<tvectorSize;j++){
			 if(outputIndices[j] <0 && probs < threshholds[j]){
					outputIndices[j]=k;
					nMatches++;
				}
			}	 
	 }
	k++;
	return(k);
}
