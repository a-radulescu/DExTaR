#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/config.hpp>
#include <cassert>
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/graph/howard_cycle_ratio.hpp>
#include <boost/graph/directed_graph.hpp>
#include <boost/graph/tiernan_all_cycles.hpp>


using namespace boost;
using namespace std;
const char* tab = "\t";
/*Graph structure*/
typedef directed_graph<> Graph;
typedef graph_traits<Graph>::vertex_descriptor VertexDir;
typedef graph_traits<Graph>::edge_descriptor EdgeDir;
const	int maxLength = 40;
const int readLength = 100;
typedef std::pair<int, int> Edge;

/*Vertices printer*/
template <typename OutputStream>
struct cycle_printer
{

    cycle_printer(OutputStream& stream)
        : os(stream)
    { }

    template <typename Path, typename Graph>
    void cycle(const Path& p, const Graph& g)
    {
		if (p.empty()||(p.size()>11))
			return;
        // Get the property map containing the vertex indices
        // so we can print them.
        typedef typename property_map<Graph, vertex_index_t>::const_type IndexMap;
        IndexMap indices = get(vertex_index, g);

        // Iterate over path printing each vertex that forms the cycle.
        typename Path::const_iterator i, before_end = boost::prior(p.end());
        for(i = p.begin(); i !=  before_end; ++i) {
            os << get(indices, *i)+1 << " ";
        }
        os << get(indices, *i)+1<<'\n';
    }

    OutputStream& os;
};

std::string readFile;
std::string patternFile;
typedef std::pair<string, string> ESPatterns;
typedef std::vector<ESPatterns> ESSet;
std::map <string, ESSet> patternSet;
int totalRepFound = 0;
int totalRepToFind = 0;


/*build the complement of a contig*/
std::string reverseCont(std::string contig){

	std::string revContig = "";
	for(int i = 1;i <= contig.length();i++){
		char nucleotide = contig[contig.length()-i];
		switch (nucleotide){
		case 'A': revContig = revContig+"T"; break;
		case 'C': revContig = revContig+"G"; break;
		case 'G': revContig = revContig+"C"; break;
		case 'T': revContig = revContig+"A"; break;
		}
	}
	return revContig;

}

/*Retrive the repetitions from the genome*/
void patternToFind(){

	ifstream patternF(patternFile.c_str());
	if(patternF.is_open()){
	
		std::string motif = "";
		std::string motifAncien = "";
		std::string line;
		ESSet vecES;
		
		while( getline(patternF,line)){
		
				std:string debut,fin;
				int pos;
				pos = line.find(tab);
				motif = line.substr(0,pos);
				std::transform(motif.begin(), motif.end(),motif.begin(), ::toupper);
				
				if(motif.compare(motifAncien) !=  0){
				
						if(motifAncien.length() == 0) motifAncien = motif;
						
						else{
						
								patternSet.insert(std::pair<string,ESSet>(motifAncien,vecES));
								motifAncien = motif;
								vecES.clear();
								
						}
						
				}
				
				for(int j = 0;j<3;j++){
				
					line = line.substr(pos+1);
					pos = line.find(tab);
					
				}	
				
				debut = line.substr(0,pos);				
				line = line.substr(pos+1);
				pos =  line.find(tab);
				fin = line.substr(pos+1);
				vecES.push_back(std::pair<string,string> (debut,fin));
			}
			
	}
	
	cout<<"Fin chargmenet motifs Genome, nb motifs: "<<patternSet.size()<<endl;
	
}



/* Retrive the cycles from file*/
std::vector<int> detCycle(string line){

	std::vector<int> cycle;
	int pos = line.find(" ");
	while(pos>0){
	
		int node = atoi(line.substr(0,pos).c_str());
		cycle.push_back(node);
		line = line.substr(pos+1);
		pos = line.find(" ");
		
	}
	
	return cycle;

}


int identifESTandem(set<string> & readsTandem, string pattern, vector<int> & nbCop,vector<string> & entries,vector<string> & exits){

	// cout <<"in identifESTandem "<< pattern <<endl;
	int cont = 2;
	vector<float> nbCopReads;
	
	for(int i = 0;i<nbCop.size();i++) nbCopReads.push_back(0);

	/*cout <<"in identifESTandem entries and exits"<<endl;
	for(int i = 0;i<entries.size();i++)
	cout <<entries[i]<<" "<<exits[i]<<endl;*/

	std::set<string>::iterator it = readsTandem.begin();
	
	while( it != readsTandem.end() ){
	
		bool readErase = false;
		string::size_type start = 0;
		while ((start = it->find(pattern, start)) !=  string::npos){
		
			float count = 1;
			int posFirst = start;
      		start = start+pattern.length();
      		
			while (start<it->length() &&  it->find(pattern, start) !=  string::npos && it->find(pattern, start) ==  start){
			
				count++;
				start += pattern.length();
				
			}
			
			int longLastCopy = 0;
			bool eq = true;
			
			while (start<it->length()&& eq && longLastCopy<pattern.length()-1){
			
				if(it->at(start) ==  pattern.at(longLastCopy)){
				
						longLastCopy++;
						start++;
						
				}else{ 
					eq = false;
				}
					
			}
			count = count+(longLastCopy/(float) pattern.length());
			
      		//cout<<"count "<<count<<endl;
    		if(count >= 2){
    		
				string before = it->substr(0,posFirst);
				string after = it->substr(min(start,it->length()));
				//cout<<"before "<<before<<" after "<<after<<endl;
				//cout<<"read "<<*it<<endl;
				
				if(before.length()>0 && after.length()>0){

						readErase = true;
						bool found = false;
						int j = 0;
						while(!found && j<entries.size()){
								bool foundE = false;
								if(entries[j].length()<before.length()){ 
								
									  //cout<<"pos before 1"<<before.substr(before.length()-entries[j].length()).compare(entries[j])<<endl;
								
									  string s = before.substr(before.length()-entries[j].length());
									  int comp = s.compare(entries[j]);
									  if(s.length()>0 && comp == 0) foundE = true;
								
								}else {
								
									//cout<<"pos before 2"<< entries[j].substr(entries[j].length()-before.length()).compare(before)<<endl;
									string s = entries[j].substr(entries[j].length()-before.length());
									int comp = s.compare(before);
									if(s.length()>0 && comp == 0) foundE = true;
								
								}
								bool foundS = false;
								if(exits[j].length()<after.length()){
									// cout<<"pos after 1"<<after.find(exits[j])<<endl;
									if(after.find(exits[j]) == 0) foundS = true;
									
								}else {//cout << "pos after 2"<<exits[j].find(after)<<endl;
								
									if(exits[j].find(after) == 0) foundS = true;
									
								}
									//cout<<foundE <<" "<<foundS<<endl;
								if( foundE && foundS ){
								
									found = true;
									nbCopReads[j] = max(nbCopReads[j],count);
								
								}
												
												
								j++;
						}
				}
		}
    }
    
	if (readErase) readsTandem.erase(it++);
	else ++it;

  }
	 //cout <<"in identifESTandem fin recherche"<<endl;
	bool multiplicite = true;
    /*for(int i = 0;i<nbCopReads.size();i++){
		 //cout <<nbCopReads[i]<<" ";
	}
	 //cout <<endl;*/
	int i = 0;
	
	while(i<nbCopReads.size()){
	
		if(nbCopReads[i]>0){
		
				cout <<"REPETITION InOneREAD"<<" "<<pattern<<" "<<nbCopReads[i]<<" "<<entries[i]<<" "<<exits[i]<< endl;
				totalRepFound = totalRepFound++;
				
				if(nbCopReads[i]>= nbCop[i]) multiplicite = false;
				
				nbCopReads.erase(nbCopReads.begin());
				nbCop.erase(nbCop.begin()+i);
				entries.erase(entries.begin()+i);
				exits.erase(exits.begin()+i);
				
		}else i++;
	}

	if(multiplicite == false) cont = 0;
	else if(entries.size() == 0 || exits.size() == 0) cont = 1;
	//cout<<"in identifESTandem fin "<<cont<<endl;
	return cont;
	
}


int identifESDiff(set<string> readsTandem,set<string> readsBegin,set<string> readsEnd, string pattern,
								  vector<int>  nbCop, vector<string>  entries, vector<string> exits){
								  
	 //cout<<"in identifESDiff "<<pattern<<endl;
	int cont = 1;

	vector<bool> foundEntries;
	for(int i = 0;i<entries.size();i++) foundEntries.push_back(false);
	vector<bool> foundExits;
	for(int i = 0;i<entries.size();i++) foundExits.push_back(false);

	std::set<string>::iterator it;

	//in tandem reads
	 //cout<<"in identifESDiff tandem reads"<<endl;
	if(pattern.length()<maxLength*2){
		for(it = readsTandem.begin(); it != readsTandem.end() ; ++it ){
			string::size_type start = 0;
			 //cout<<"READ "<<*it<<endl;
			while ((start = it->find(pattern, start)) !=  string::npos){
				int posFirst = start;
        float count = 1;
				start = start+pattern.length();
				 //cout<<"debut rep tandem "<<posFirst<<"start "<< start<<" next "<< it->find(pattern, start) <<endl;
				while (start<it->length() && it->find(pattern, start) !=  string::npos && it->find(pattern, start) ==  start){
          count++;
					start +=  pattern.length();
				}
				 //cout<<"fin rep tandem "<<start<<endl;
				if(count == 1) continue;
				int longLastCopy = 0;
				bool eq = true;
				while ( start<it->length() && eq && longLastCopy<pattern.length()-1){
					if(it->at(start) ==  pattern.at(longLastCopy)){
							longLastCopy++;start++;
					}else eq = false;
				}
				string before = it->substr(0,posFirst);
				string after = it->substr(min(start,it->length()));
				 //cout<<"before "<<before<<endl;

				 //cout<<"after "<<after<<endl;

				if(before.length()>0){
						int j = 0;
						while(j<entries.size()){
								if(entries[j].length()<before.length()){
                    string s = before.substr(before.length()-entries[j].length());
                    int comp = s.compare(entries[j]);
										if(s.length()>0 && comp == 0) { foundEntries[j] = true;}
                }
								else {string s = entries[j].substr(entries[j].length()-before.length());
                  int comp = s.compare(before);
                  if(s.length()>0 && comp == 0) { foundEntries[j] = true;}}
								j++;
						}
				}else if (after.length()>0){
							int j = 0;
							while(j<entries.size()){
								if(exits[j].length()<after.length())
										if(after.find(exits[j]) == 0){foundExits[j] = true;}
								else if(exits[j].find(after) == 0){foundExits[j] = true;}
								j++;
							}
				}
			}
		}
	}

	 //cout<<"in identifESDiff readsBegin"<<endl;
	//in readsBegin
	string patt = pattern;
	if(pattern.length()>= maxLength*2) patt = pattern.substr(0,maxLength*2);
	for(it = readsBegin.begin(); it != readsBegin.end() ; ++it ){

			 //cout<<"READ "<<*it<<endl;
		string::size_type start = it->find(patt);
		string before = it->substr(0,min(start,it->length()));
		 //cout<<"before "<<before<<endl;

		if(before.length()>0){
					int j = 0;
					while(j<entries.size()){
            //cout<<"entries "<<entries[j]<<endl;
							if(entries[j].length()<before.length()){//cout<<"pos before 1 "<<before.substr(before.length()-entries[j].length())<<endl;
                string s = before.substr(before.length()-entries[j].length());
                int comp = s.compare(entries[j]);
									if(s.length()>0 && comp == 0) {foundEntries[j] = true;}}
							else {//cout<<"pos before 2 "<<entries[j].substr(entries[j].length()-before.length())<<endl;
                string s = entries[j].substr(entries[j].length()-before.length());
                int comp = s.compare(before);
                if(s.length()>0 && comp == 0) {foundEntries[j] = true;}
              }
							j++;
					}
			}
	}

	 //cout<<"in identifESDiff readsEnd"<<endl;
	//in readsEnd
	patt = pattern;
	if(pattern.length()>= maxLength*2) patt = pattern.substr(pattern.length()-maxLength*2);

	 //cout<<pattern<< " "<<patt<<endl;
	for(it = readsEnd.begin(); it != readsEnd.end() ; ++it ){
			 //cout<<"READ "<<*it<<endl;
			string::size_type start = it->find_last_of(patt);
			int longLastCopy = 0;
			bool eq = true;
			while ( eq && longLastCopy<pattern.length()-1 && start<it->length()){
				if(it->at(start) ==  pattern.at(longLastCopy)){
						longLastCopy++;start++;
				}else eq = false;
			}
			 //cout<<start<<endl;
			string after = it->substr(min(start,it->length()));
			 //cout<<"after "<<after<<endl;
			if (after.length()>0){
						int j = 0;
						while(j<entries.size()){
							if(exits[j].length()<after.length())
									if(after.find(exits[j]) == 0){foundExits[j] = true;}
							else if(exits[j].find(after) == 0){foundExits[j] = true;}
							j++;
						}
			}

	}
	 //cout<<"in identifESDiff fin recherche"<<endl;
	bool trouveAuMoinsUnePaireES = false;
	for(int i = 0;i<entries.size();i++){
		 //cout<<foundEntries[i]<<" "<<foundExits[i]<<endl;
		if(foundEntries[i]&&foundExits[i]){
				trouveAuMoinsUnePaireES = true;
				cout<<"REPETITION InTwoREADS"<<" "<<pattern<<" "<<nbCop[i]<<" "<<entries[i]<<" "<<exits[i]<< endl;
				totalRepFound = totalRepFound++;
		}
	}
	if(trouveAuMoinsUnePaireES) cont = 0;

	return cont;
}
//identification in Reads
void identificationReads(std::vector<string> possiblePatterns, std::vector<int> possibleNbCopies,std::vector<string> possibleEntries,std::vector<string> possibleExits){
	 //cout<<"in The Reads patterns:"<<endl;
	/*for(int i = 0;i<possiblePatterns.size();i++)
		 //cout<<possiblePatterns[i]<<" ";
	 //cout<<"entries"<<endl;
	for(int i = 0;i<possibleEntries.size();i++)
		cout<<possibleEntries[i]<<" ";
cout<<endl;*/
	ifstream readFileCheck(readFile.c_str());
	string line;
	if(readFileCheck.is_open()){
		getline(readFileCheck,line);
		getline(readFileCheck,line);
	}else cout << "Unable to open read file"<<endl;
  vector<string>::iterator pPat;

	for(pPat = possiblePatterns.begin();pPat<possiblePatterns.end();pPat++){//for each possible pattern
			cout<<"pattern to search in the reads "<<*pPat<<endl;
				string pattern = *pPat;
				int patternLength = pattern.length();


				vector<int> nbCop;
				vector<string> entries;
				vector<string> exits;

        vector<string>::iterator pPat2;
        int j = 0;
        int m = 0;
      	for(pPat2 = possiblePatterns.begin();pPat2<possiblePatterns.end();){//for each possible pattern
	       string pattern2 = *pPat2;
					if(pattern2.compare(pattern) == 0){
              m++;
							nbCop.push_back(possibleNbCopies[j]);
							entries.push_back(possibleEntries[j]);
							exits.push_back(possibleExits[j]);
              if(m>1) 
              { 
                pPat2 = possiblePatterns.erase(pPat2);
                possibleNbCopies.erase(possibleNbCopies.begin()+j);
                possibleEntries.erase(possibleEntries.begin()+j);
                possibleExits.erase(possibleExits.begin()+j);
              }else {++pPat2;j++;}  
					}else {++pPat2;j++;}
				}
				//cout<<nbCop.size()<<" "<<entries.size()<<" "<<exits.size()<<endl;
				set<string> readsTandem;//reads qui attestent la presence en tandem
				set<string> readsEnd;//reads avec des seq avant
				set<string> readsBegin;//reads avec des seq apres

				string::size_type start = 0; //index de recherche dans les reads
		if(pattern.length()<maxLength){//case 1
				cout<<"case 1"<<endl;
			//cout<<"start "<<start<<" line length "<<line.length()<<endl;
				  pattern = pattern+pattern;	
          patternLength = pattern.length();
        //search the reads
					while ((start = line.find(pattern, start)) !=  string::npos) {

						//cout<<"start "<<start<<endl;
						unsigned int posEndRead = line.find("$",start);
						//cout<<"posEndRead "<<posEndRead<<endl;
            unsigned int posBeginRead = line.substr(0,start).find_last_of("$");
            std::string read = line.substr(posBeginRead+1, posEndRead-posBeginRead-1);
  					//std::string read = line.substr(posEndRead-readLength,readLength);
            int readLength2 = read.length();
						//cout<<"read "<<read<<endl;
						int posStartInRead  = 	readLength2-(posEndRead-start);
						//cout<<"posStartInRead "<<posStartInRead<<endl;
						std::string before = read.substr(max(0,posStartInRead-patternLength),min(posStartInRead,patternLength));
						//cout<<"before "<<before<<endl;
						std::string after;
						if(posStartInRead+patternLength>= readLength2) after = "";
						else
								after = read.substr(posStartInRead+patternLength,min(patternLength,readLength2-posStartInRead-patternLength));
						//cout<<"after "<<after<<endl;

						if(pattern.compare(before) == 0 || pattern.compare(after) == 0){
								readsTandem.insert(read);
								//cout<<"read goes in tandem"<<endl;
								start = posEndRead;
						}else if (before.length()<pattern.length() && pattern.find(before) == (pattern.length()-before.length())){
								readsEnd.insert(read);
								//cout<<"read goes in readsEnd"<<endl;
								start +=  pattern.length();
						}else if (after.length()<pattern.length() && pattern.find(after) == 0){
								readsBegin.insert(read);
								//cout<<"read goes in readsBegin"<<endl;
								start = posEndRead;
						}else start +=  pattern.length();

					}

					//cout<<"test1 "<< readsTandem.size()<<endl;
					if(readsTandem.size() == 0) continue; //false positive

					//identification entree+motif*nbCopies+sortie
					int cont1 = identifESTandem(readsTandem,pattern,nbCop,entries,exits);
					if (cont1 == 0) break;//toute la multiplicité du cyle à été utilisé
					if (cont1 == 1) continue; //il reste de la multiplicité et des entrées et de sorties

					//identification des couples des entrees/sorties
					int cont2 = identifESDiff(readsTandem,readsBegin,readsEnd,pattern,nbCop,entries,exits);
					if (cont2 == 0) break;
		}else if(pattern.length()<maxLength*2){
							cout<<"case 2"<<endl;
							while ((start = line.find(pattern, start)) !=  string::npos) {

								unsigned int posEndRead = line.find("$",start);
                unsigned int posBeginRead = line.substr(0,start).find_last_of("$");
                std::string read = line.substr(posBeginRead+1, posEndRead-posBeginRead-1);
                int readLength2 = read.length();
						
								//std::string read = line.substr(posEndRead-readLength,readLength);
								int posStartInRead  = 	readLength2-(posEndRead-start);
								std::string before = read.substr(max(0,posStartInRead-patternLength),min(posStartInRead,patternLength));
								std::string after;
								if(posStartInRead+patternLength>= readLength2) after = "";
								else
									after = read.substr(posStartInRead+patternLength,min(patternLength,readLength2-posStartInRead-patternLength));

								if(pattern.compare(before) == 0 || pattern.compare(after) == 0){
										readsTandem.insert(read);
										start = posEndRead;
								}else if (before.length()<pattern.length()  && pattern.find(before) == (pattern.length()-before.length())){
										readsEnd.insert(read);
										if( before.length()>10)
											readsTandem.insert(read);
										start +=  pattern.length();
								}else if (after.length()<pattern.length() && pattern.find(after) == 0){
                    readsBegin.insert(read);
										if( after.length()>10)
											readsTandem.insert(read);
										start = posEndRead;
								}else start +=  pattern.length();

							}
							if(readsTandem.size() == 0) continue; //false positive

					//identification des couples des entrees/sorties
					int cont2 = identifESDiff(readsTandem,readsBegin,readsEnd,pattern,nbCop,entries,exits);
					if (cont2 == 0) break;

		}else {//motif>= 80
						cout<<"case 3"<<endl;
							std::string patternBegining = pattern.substr(0,80);
							int fragLength = 80;
							while ((start = line.find(patternBegining, start)) !=  string::npos) {

								unsigned int posEndRead = line.find("$",start);
                unsigned int posBeginRead = line.substr(0,start).find_last_of("$");
                std::string read = line.substr(posBeginRead+1, posEndRead-posBeginRead-1);
                int readLength2 = read.length();
						

								//std::string read = line.substr(posEndRead-readLength,readLength);
								int posStartInRead  = 	readLength2-(posEndRead-start);
								std::string before = read.substr(max(0,posStartInRead-fragLength),min(posStartInRead,fragLength));
								std::string after;
								if(posStartInRead+fragLength>= readLength2) after = "";
								else
									after = read.substr(posStartInRead+fragLength,min(fragLength,readLength2-posStartInRead-fragLength));

								if(pattern.find(before) == pattern.length()-before.length() && (after.length() == 0 || pattern.find(after) == 80)){
										readsTandem.insert(read);
										start = posEndRead;

								}else if (after.length() == 0 || pattern.find(after) == 80){
										readsEnd.insert(read);
										start = posEndRead;
								}else start +=  patternBegining.length();

							}

							std::string patternEnding = pattern.substr(patternLength-80);
							while ((start = line.find(patternEnding, start)) !=  string::npos) {

								unsigned int posEndRead = line.find("$",start);
                unsigned int posBeginRead = line.substr(0,start).find_last_of("$");
                std::string read = line.substr(posBeginRead+1, posEndRead-posBeginRead-1);
                int readLength2 = read.length();
						

								//std::string read = line.substr(posEndRead-readLength,readLength);
								int posStartInRead  = 	readLength2-(posEndRead-start);
								std::string before = read.substr(max(0,posStartInRead-fragLength),min(posStartInRead,fragLength));
								std::string after;
								if(posStartInRead+fragLength>= readLength2) after = "";
								else
									after = read.substr(posStartInRead+fragLength,min(fragLength,readLength2-posStartInRead-fragLength));

								if(pattern.find(after) == 0 && (before.length() == 0 || pattern.find(before) == patternLength-80-before.length())){
										readsTandem.insert(read);
										start = posEndRead;

								}else if (before.length() == 0 || pattern.find(before) == (patternLength-80-before.length())){
										readsBegin.insert(read);
										start = posEndRead;
								}else start +=  patternBegining.length();

							}
							if(readsTandem.size() == 0) continue; //false positive

						//identification des couples des entrees/sorties
					int cont2 = identifESDiff(readsTandem,readsBegin,readsEnd,pattern,nbCop,entries,exits);
					if (cont2 == 0) break;
		}

	}
}
bool compEI( ESSet & vectIO, std::string startSeq, std::string endSeq){
							int l = 0;
							bool found = false;
							while(l<vectIO.size()&&!found){
								if(startSeq.length()<vectIO[l].first.length()){
		              string s = vectIO[l].first.substr(vectIO[l].first.length()-startSeq.length());
                  int comp = s.compare(startSeq);
									if(comp == 0) found = true;
								}else {		             
									string s = startSeq.substr(startSeq.length()-vectIO[l].first.length());
                  int comp = s.compare(vectIO[l].first);
									if(comp == 0) found = true;
								}
								if(found){
									if(endSeq.length()<vectIO[l].first.length()){
										string s = vectIO[l].second.substr(0,endSeq.length());
                  	int comp = s.compare(endSeq);
										if(comp == 0) found = false;
									}else {
										string s = endSeq.substr(0,vectIO[l].second.length());
                  	int comp = s.compare(vectIO[l].second);
										if(comp == 0) found = false;
									}
								}
								if(found) vectIO.erase(vectIO.begin()+l);
								else l++;
							}
	return found;
}
void retrivePatterns(std::vector<int> cycle, std::vector<Edge> edge_in, std::vector<Edge> edge_out,std::vector<string> & patternsInCycle, std::vector< std::pair<Edge,Edge> > & ioArcs,int kmerLength, std::vector<string> contigList) {
			for(int i = 0;i<edge_in.size();i++){
				for(int j = 0;j<edge_out.size();j++)	{
						std::string pattern;
						int startNode = 0;
						int endNode = 0;
						for(int l = 0;l<cycle.size();l++){
							if(cycle[l] == edge_in[i].second) startNode = l;
							if(cycle[l] == edge_out[j].first) endNode = l;
						}
						pattern = contigList[cycle[startNode]-1];
						if(startNode <= endNode){
						
							if(startNode<endNode){
								for(int l = startNode+1;l <= endNode;l++){
									pattern = pattern+contigList[cycle[l]-1].substr(kmerLength-1);
								}
							}else{
		
								for(int l = startNode+1;l<cycle.size();l++){
									pattern = pattern+contigList[cycle[l]-1].substr(kmerLength-1);
								}
								for(int l = 0;l <= endNode;l++){
									pattern = pattern+contigList[cycle[l]-1].substr(kmerLength-1);
								}
							}
						}else{

							if(startNode<cycle.size()-1)
								for(int l = startNode+1;l<cycle.size();l++){
									pattern = pattern+contigList[cycle[l]-1].substr(kmerLength-1);
								}
							for(int l = 0;l <= endNode;l++){
								pattern = pattern+contigList[cycle[l]-1].substr(kmerLength-1);
							}

					}
					int length = 2;
					string atomPattern = "";
					string smallPattern = "";
					string restPattern = "";
					int nbCopiesAtomic;
					bool repAtom = false;
					while(!repAtom&&length <= pattern.length()/2){
							smallPattern = pattern.substr(0,length);
							nbCopiesAtomic = 1;
							bool repInt = true;
							restPattern = pattern.substr(length*nbCopiesAtomic);
							while(repInt&&restPattern.length()>= length){;
									size_t pos = restPattern.find(smallPattern);
									if(restPattern.find(smallPattern) != 0) repInt = false;
									nbCopiesAtomic++;
									//cout<<smallPattern<<" "<<pos<<" "<<restPattern<<endl;
									restPattern = pattern.substr(length*nbCopiesAtomic);
							}
							if (repInt == true) {
										repAtom = true;
										atomPattern = smallPattern;
							}
							length++;
					}
					if(!repAtom){
						length = ceil(pattern.length()/2);
						nbCopiesAtomic = 1;
						while(!repAtom&&length<pattern.length()){
							smallPattern = pattern.substr(0,length);
							restPattern = pattern.substr(length);
							//cout<<smallPattern<<" "<<pos<<" "<<restPattern<<endl;
							if(smallPattern.find(restPattern) == 0) {
		 							 repAtom = true;
									 nbCopiesAtomic++;
				           atomPattern = smallPattern;
							}
							length++;
						}
					}
					if(atomPattern.length()>0) pattern = atomPattern;
					string startSeq = contigList[edge_in[i].first-1].substr(0,contigList[edge_in[i].first-1].length()-kmerLength+1);
					string endSeq = contigList[edge_out[j].second-1].substr(kmerLength-1);
					std::map <string, ESSet> ::iterator it;
					it = patternSet.find(pattern);
					if(it != patternSet.end()){
							bool found = compEI(it->second,startSeq,endSeq);
							if (found) {
									patternsInCycle.push_back(pattern);
									ioArcs.push_back(std::pair<Edge,Edge>(edge_in[i],edge_out[j]));
							}
					}else{ 
							pattern = reverseCont(pattern);
							startSeq = reverseCont(startSeq);
							endSeq = reverseCont(endSeq);
							it = patternSet.find(pattern);
							if(it !=  patternSet.end()){
									bool found = compEI(it->second,startSeq,endSeq);
									if (found) {
											patternsInCycle.push_back(pattern);
											ioArcs.push_back(std::pair<Edge,Edge>(edge_in[i],edge_out[j]));
									}
							}
					}
				}
			}
	

}

/*Main part: remove supp freq, test remaining freq and print if ok*/
void removeSuppFreq(std::vector<int> cycle, std::vector<int> freqCycle, std::vector<Edge> edge_in, std::vector<Edge> edge_out,
 std::vector<int> freq_in, std::vector<int> freq_out, float cover, int kmerLength, std::vector<string> contigList){
  //std::vector<string> & patternsInCycle, std::vector< std::pair<Edge,Edge> > & ioArcs){

	std::vector<int> newFreq;
	std::vector<int> suppFreqPos(cycle.size());
	std::vector<int> suppFreqNeg(cycle.size());
	std::vector<int> negFreqNodes;

	//compute for each node the in freq (suppFreqPos) and the out freq (suppFreqNeg) to discover the negative freq nodes
	for(int i = 0;i<cycle.size();i++){
		//cout<<cycle[i]<<" ";
		for(int j = 0;j<edge_in.size();j++){
			if(edge_in[j].second == cycle[i]) suppFreqPos[i] = suppFreqPos[i]+freq_in[j];
		}
		for(int j = 0;j<edge_out.size();j++){
			if(edge_out[j].first == cycle[i]) suppFreqNeg[i] = suppFreqNeg[i]-freq_out[j];
		}
		if(suppFreqPos[i]+suppFreqNeg[i]<0) negFreqNodes.push_back(i);
	}
	if(negFreqNodes.size()>0){
		negFreqNodes.push_back(negFreqNodes[0]);
	}else{
		negFreqNodes.push_back(0);
		negFreqNodes.push_back(0);
	}

	std::vector<int> freqReal(cycle.size());
	int freqCumulate;
	/*****if indetifReads**/
	std::vector <string> possibleEntries;
	std::vector <string> possibleExits;
	std::vector <string> possiblePatterns;
	std::vector <int> possibleNbCopies;
	/*****end if indetifReads**/
	for(int i = 0;i<edge_in.size();i++){
		for(int j = 0;j<edge_out.size();j++)	{//for each pair of possible entries and exits

			freqCumulate = 0;
			int startNode = 0;
			int endNode = 0;
			for(int l = 0;l<cycle.size();l++){
				if(cycle[l] == edge_in[i].second) startNode = l;
				if(cycle[l] == edge_out[j].first) endNode = l;
			}

			//remouve the min freq between the entry and exit = the number of times we go throught the repetiton
			if(freq_in[i]>= freq_out[j]) {
				suppFreqPos[startNode] = suppFreqPos[startNode]-freq_out[j];
				suppFreqNeg[endNode] = suppFreqNeg[endNode]+freq_out[j];

			}else{
				suppFreqPos[startNode] = suppFreqPos[startNode]-freq_in[i];
				suppFreqNeg[endNode] = suppFreqNeg[endNode]+freq_in[i];
			}

			//remouve supp feq by the negative freq nodes
			for(int l = 0; l<negFreqNodes.size()-1;l++){

				if(negFreqNodes[l]<negFreqNodes[l+1]){

					for(int m = negFreqNodes[l]+1;m <= negFreqNodes[l+1];m++){
						freqReal[m] = freqCycle[m]-suppFreqPos[m]-freqCumulate;
						freqCumulate = freqCumulate+suppFreqPos[m]+suppFreqNeg[m];
						if (freqCumulate < 0) freqCumulate = 0;
					}
				}
				else{
					for(int m = negFreqNodes[l]+1;m<cycle.size();m++){
						freqReal[m] = freqCycle[m]-suppFreqPos[m]-freqCumulate;
						freqCumulate = freqCumulate+suppFreqPos[m]+suppFreqNeg[m];
						if(freqCumulate<0) freqCumulate = 0;
					}
					for(int m = 0;m <= negFreqNodes[l+1];m++){
						freqReal[m] = freqCycle[m]-suppFreqPos[m]-freqCumulate;
						freqCumulate = freqCumulate+suppFreqPos[m]+suppFreqNeg[m];
						if(freqCumulate<0) freqCumulate = 0;
					}
				}
			}
			//replace the freq of the entry and the exit on the nodes for the next possible cupple
			int freqIn;
			if(freq_in[i]>= freq_out[j]) {
				suppFreqPos[startNode] = suppFreqPos[startNode]+freq_out[j];
				suppFreqNeg[endNode] = suppFreqNeg[endNode]-freq_out[j];
				freqIn = freq_out[j];
			}else{
				suppFreqPos[startNode] = suppFreqPos[startNode]+freq_in[i];
				suppFreqNeg[endNode] = suppFreqNeg[endNode]-freq_in[i];
				freqIn = freq_in[i];
			}

			//compute approximation of the freq by dividing it by the cover
			bool repetition = true;
			for(int l = 0;l<freqReal.size();l++){
				if(freqReal[l]<0) {
						freqReal[l] = 0;
						repetition = false; //if there is a node with freq 0 then there is no repetiton in the cycle
				}
				else {
					 float fValeur = freqReal[l]/cover;
				   float fDecimal = fValeur-floor(fValeur);
          // if (fDecimal< 0.5)
          //   freqReal[l] = floor(fValeur);
          // else
             freqReal[l] = ceil(fValeur);
				}
			}

			//test for passing at least two times in the first node and equal times in the start and the end node
			if((freqReal[startNode] != freqReal[endNode])||(freqReal[startNode] <= 1)) repetition = false;


			//test if the remaining freq can give a tandem repetiton and compute the pattern
			int nbCopies;
			std::string pattern = contigList[cycle[startNode]-1];
			if(repetition){
				nbCopies = freqReal[startNode];
				if(startNode <= endNode){
						//cout<<"if\n";
						for(int l = startNode;l <= endNode;l++){
							if(freqReal[l] != nbCopies) repetition = false;
						}
						//cout<<"Second condition:"<<repetition<<"\n";
						if(startNode<endNode){
							for(int l = startNode+1;l <= endNode;l++){
								pattern = pattern+contigList[cycle[l]-1].substr(kmerLength-1);
							}
						}else{
							if(freqReal[startNode]<2) repetition = false;
							for(int l = startNode+1;l<cycle.size();l++){
								pattern = pattern+contigList[cycle[l]-1].substr(kmerLength-1);
							}
							for(int l = 0;l <= endNode;l++){
								pattern = pattern+contigList[cycle[l]-1].substr(kmerLength-1);
							}
						}
						if(endNode<cycle.size()-1){
							for(int l = endNode+1;l<cycle.size();l++){
								if(freqReal[l] != nbCopies-1) repetition = false;
							}
						}
					  //cout<<"Third condition:"<<repetition<<"\n";
						if(startNode>0){
							for(int l = 0;l <= startNode-1;l++){
								if(freqReal[l] != nbCopies-1) repetition = false;
							}
						}
						// cout<<"Fourth condition:"<<repetition<<"\n";
				}else{
						// cout<<"else\n";
						for(int l = startNode;l<cycle.size();l++){
							if(freqReal[l] != nbCopies) repetition = false;
						}
						// cout<<"Second condition:"<<repetition<<"\n";
						if(startNode<cycle.size()-1)
							for(int l = startNode+1;l<cycle.size();l++){
								pattern = pattern+contigList[cycle[l]-1].substr(kmerLength-1);
							}
						for(int l = 0;l <= endNode;l++){
							if(freqReal[l] != nbCopies) repetition = false;
							pattern = pattern+contigList[cycle[l]-1].substr(kmerLength-1);
						}
					 	//cout<<"Third condition:"<<repetition<<"\n";
						for(int l = endNode+1;l <= startNode-1;l++){
							if(freqReal[l] != nbCopies-1) repetition = false;
						}
						// cout<<"Fourth condition:"<<repetition<<"\n";
				}

			}

			//test the length of the pattern matches the length of the cycle
			bool LengthPattern = true;
			int lengthCycle = 0;
			for (int l = 0;l<cycle.size();l++){
					lengthCycle = lengthCycle+contigList[cycle[l]-1].length();
			}
			lengthCycle = lengthCycle-(kmerLength-1)*cycle.size();
			if (lengthCycle>pattern.length()) repetition = false;


			bool repetitionOneTime = false;
			size_t testOk = 2;
			if((repetition == true)&&(startNode ==  endNode && freqReal[startNode] ==  2)){
						repetition = false;
					if( cycle.size() <= 2){
						int l = floor(kmerLength/2);
						//cout<<"l "<<l<<endl;
						string smallPattern = "";
						string restPattern = "";
						while(!repetitionOneTime&&l<floor(pattern.length()/2)){
						  smallPattern = pattern.substr(0,l);
							//cout<<"smallPattern"<<smallPattern<<endl;
							restPattern = pattern.substr(l);
							//cout<<"restPattern"<<restPattern<<endl;
							if(restPattern.find(smallPattern) == testOk) repetitionOneTime = true;
							l++;
						}
					//if (repetitionOneTime) cout<<"OneTime l final"<<l<<" "<<smallPattern<<" "<<restPattern<< " "<<restPattern.find(smallPattern) <<endl;
				}
			}


			if(repetition || repetitionOneTime){
					//compute the number of occurences
					 float fValeur = freqIn/cover;
				   float fDecimal = fValeur-floor(fValeur);
           if (fDecimal< 0.5)
             freqIn = floor(fValeur);
           else
             freqIn = ceil(fValeur);
					//compute the atomic pattern
					int nbCopiesAtomic = 0;
					int length = 2;
					string atomPattern = "";
					string smallPattern = "";
					string restPattern = "";
					bool repAtom = false;
					while(!repAtom&&length <= pattern.length()/2){
							smallPattern = pattern.substr(0,length);
							nbCopiesAtomic = 1;
							bool repInt = true;
							restPattern = pattern.substr(length*nbCopiesAtomic);
							while(repInt&&restPattern.length()>= length){;
									size_t pos = restPattern.find(smallPattern);
									if(restPattern.find(smallPattern) != 0) repInt = false;
									nbCopiesAtomic++;
									//cout<<smallPattern<<" "<<pos<<" "<<restPattern<<endl;
									restPattern = pattern.substr(length*nbCopiesAtomic);
							}
							if (repInt == true) {
										repAtom = true;
										atomPattern = smallPattern;
							}
							length++;
					}
					if(!repAtom){
						length = ceil(pattern.length()/2);
						nbCopiesAtomic = 1;
						while(!repAtom&&length<pattern.length()){
							smallPattern = pattern.substr(0,length);
							restPattern = pattern.substr(length);
							//cout<<smallPattern<<" "<<pos<<" "<<restPattern<<endl;
							if(smallPattern.find(restPattern) == 0) {
		 							 repAtom = true;
									 nbCopiesAtomic++;
				           atomPattern = smallPattern;
							}
							length++;
						}
					}
				//cout<<atomPattern<<endl;
				string startSeq = contigList[edge_in[i].first-1].substr(0,contigList[edge_in[i].first-1].length()-kmerLength+1);
				string endSeq = contigList[edge_out[j].second-1].substr(kmerLength-1);

				//cout<<"begins "<<contigList[edge_in[i].first-1]<<" "<<contigList[edge_in[i].second-1]<<" "<<kmerLength<<" "<<startSeq<<endl;

				//cout<<"ends "<<contigList[edge_out[j].first-1]<<" "<<contigList[edge_out[j].second-1]<<" "<<kmerLength<<" "<<endSeq<<endl;
	/*****if not indentifReads
				if(repAtom){
						cout<<"REPETITION "<<" "<<atomPattern<<" "<<nbCopies*nbCopiesAtomic<< " "<< freqIn <<" "<<startSeq<<" "<<endSeq<< endl;
				}else{
						cout<<"REPETITION "<<" "<<pattern<<" "<<nbCopies<< " "<< freqIn <<" "<<startSeq<<" "<<endSeq<< endl;
				}
**/
				//print cycle information
				/*cout<<"CYCLE: "<< cycle.size() << endl;
				for(int l = 0; l<cycle.size();l++){
						cout<<cycle[l]<<" "<<contigList[cycle[l]-1]<<" "<< freqCycle[l]<<"\n";
				}
				cout<<"FREQUENCY: ";
				for(int l = 0;l<freqReal.size();l++){
					cout<<freqReal[l]<<" ";
				}
				cout<<"\n";
				cout<<"edge in" <<edge_in[i].first<<"->"<<edge_in[i].second<<"; ";

				cout<<"edge in" <<contigList[edge_in[i].first-1]<<"->"<<contigList[edge_in[i].second-1]<<"; ";
				cout<<"edge out" <<edge_out[j].first<<"->"<<edge_out[j].second<<"\n";

				cout<<"edge out" <<contigList[edge_out[j].first-1]<<"->"<<contigList[edge_out[j].second-1]<<"\n";	*/
	/****/

	/*****if indetifReads*/
				possibleEntries.push_back(startSeq);
				possibleExits.push_back(endSeq);
				if(repAtom){
						possiblePatterns.push_back(atomPattern);
						possibleNbCopies.push_back(nbCopies*nbCopiesAtomic);
				}else{
						possiblePatterns.push_back(pattern);
						possibleNbCopies.push_back(nbCopies);
				}
/**/
				//initialization of frequency
				for(int l = 0;l<freqReal.size();l++){
					freqReal[l] = 0;
				}

			}/*else{
						int p = 0;
						while(p<patternsInCycle.size()){
							if(	edge_in[i].first == ioArcs[p].first.first && 
									edge_in[i].second == ioArcs[p].first.second &&
									edge_out[j].first == ioArcs[p].second.first && 
									edge_out[j].second == ioArcs[p].second.second){
						cout<<"REPETITION LOST FREQUENCY " << patternsInCycle[p]
								<<" "<<ioArcs[p].first.first<<"->"<<ioArcs[p].first.second
 								<<" "<<ioArcs[p].second.first<<"->"<<ioArcs[p].second.second
								<<endl;
						patternsInCycle.erase(patternsInCycle.begin()+p);
						ioArcs.erase(ioArcs.begin()+p);
						}else p++;
						}

			}*/
		}
	}
  cout<<"Identification Reads"<<endl;
	if(possiblePatterns.size()>0)
		identificationReads(possiblePatterns,possibleNbCopies,possibleEntries,possibleExits);

}







/*retrive the contigs smaller than 2*k for Velvet*/
void retriveVelvetShortContigs(std::vector<Edge> edge_array, std::vector<string> & contigList, int kmerLength){

	std::set<int> IDsmallContigs;
	std::set<int>::iterator it;
	//count the number of small contigs
	int smallContigs = 0;
	for(int i = 0;i<contigList.size()/2;i++){
		if((contigList[i].length() <= 2*kmerLength)&&(contigList[i].length()>0)){
   		smallContigs++;
			IDsmallContigs.insert(i);
		}
	}
	cout<<"small Contigs"<< smallContigs<<endl;
	int contigsSolved = 0;
	int steps = 0;
  int contigsNotSolvedBefore = IDsmallContigs.size();
  int contigsNotSolvedAfter = IDsmallContigs.size()-1;
	while ((contigsSolved<smallContigs)&&(steps<smallContigs+2)&&(contigsNotSolvedBefore>contigsNotSolvedAfter)){
    contigsNotSolvedBefore = IDsmallContigs.size();
		it = IDsmallContigs.begin();
	  while( it != IDsmallContigs.end()){
				int id = (int)*it;
				//cout<<"contig id"<<id+1<<endl;
				bool canBeSolved = false;
				int j = 0;
				while((j<edge_array.size())&&(!canBeSolved)){
					if((edge_array[j].second == (id+1))&&(IDsmallContigs.count(edge_array[j].first-1) == 0)) {
							canBeSolved = true;
							std::set<int>::iterator current = it++;
							IDsmallContigs.erase(current);
							contigsSolved++;
							//cout<<"canBeSolved "<<id+1 <<" with "<<edge_array[j].first<<endl;
							std::string kmerMissing = contigList[edge_array[j].first-1].substr(contigList[edge_array[j].first-1].length()-kmerLength+1);
							std::string rest = reverseCont(contigList[id+contigList.size()/2].substr(0,contigList[id].length()-kmerLength+1));
							//cout<<"contig avant "<<contigList[id]<<"contig apres"<<kmerMissing+rest<<endl;
							contigList[id] = kmerMissing+rest;
							contigList[id+contigList.size()/2] = reverseCont(contigList[id]);
					}
					j++;
				}
			if(!canBeSolved) it++;
		}
    contigsNotSolvedAfter = IDsmallContigs.size();
		steps++;
	}
	cout<<"nb of contigs not solved: "<<IDsmallContigs.size()<<endl;
	if(IDsmallContigs.size()>0){
		cout<<"contigs not solved: ";
		for(it = IDsmallContigs.begin();it != IDsmallContigs.end();it++){
			cout<<*it<<" ";
		}
		cout<<"\n";
	}

}

int main (int argc, char* argv[]) {

	//Graph gDir;//the graph for boost
	std::string graph;//the file with the edges
	std::string contig;//the file with the nodes
	std::string k1mers;//the file with the k+1 mers and theirs frequency
	std::string cyclesFile;
	float cover;//the cover of the sequencing
	int kmerLength;
	bool velvetFiles = false;
	std::string arg;
	if (argc < 6) {
		std::cout << "Usage is -g <graphFile> -c <contigFile> -k <k+1 mers file> -l <cover> -f <cycleFile> -r <readFile> -p <patternFile> -v \n";
		std::cin.get();
		exit(0);
	}else {

        for (int i = 1; i < argc; i++) {
						arg = std::string(argv[i]);
            if (arg ==  "-g") {
                graph = argv[i + 1];
            } else if (arg ==  "-c") {
                contig = argv[i + 1];
						} else if (arg ==  "-k") {
                k1mers = argv[i + 1];
            } else if (arg ==  "-l") {
                cover = atof(argv[i + 1]);
            } else if (arg ==  "-f") {
                cyclesFile = argv[i + 1];
						} else if(arg ==  "-r") {
								readFile = argv[i + 1];
						//} else if(arg ==  "-p") {
							//	patternFile = argv[i + 1];
						} else if(arg ==  "-v") {
								velvetFiles = true;
						}

        }
	}
	//patternToFind();
	ifstream graphFile(graph.c_str());
	ifstream contigFile(contig.c_str());
	ifstream k1mersFile(k1mers.c_str());
	std::vector<string> contigList;//contig sequences
	std::vector<int> contigFrequency;//contig multiplicity


	string line;
	cout<<"creation des noeuds et de leur frequence \n";

	if (contigFile.is_open()){
		while ( getline(contigFile,line) ){
		  contigList.push_back(line);
		  getline(contigFile,line);
		  int pos =  line.find_last_of(" ");
		  int freq = atoi(line.substr(pos+1).c_str());
		  contigFrequency.push_back(freq);
		}
		contigFile.close();
		cout<<"nb nodes"<<contigList.size()<<endl;
	}else cout << "Unable to open contig file \n";

	/* //construction of a graph for boost library
		const int num_nodes = contigList.size();
		vector<VertexDir> v(num_nodes*10);
		for(size_t i = 1; i  <=  num_nodes; ++i) {
		      v[i] = add_vertex(gDir);
		  }
		cout<<"nb nodes"<<num_nodes<< " " <<gDir.num_vertices()<<"\n";
	*/

	cout<<"creation des acrs \n";

	//	std::vector<int> weights; //boost
	std::vector<Edge> edge_array;
	if(graphFile.is_open()){
		getline(graphFile,line);
		while ( getline(graphFile,line) ){
			int pos =  line.find(" -> ");
			int node1 = atoi(line.substr(0,pos).c_str());
			int node2 = atoi(line.substr(pos+3).c_str());
			//cout<<node1 << " " << node2 <<" " << contigList[node2-1].length() <<"\n";
			//add_edge(v[node1],v[node2], gDir);//boost
			edge_array.push_back(Edge(node1,node2));
			//weights.push_back(contigList[node2-1].length());//boost
		}
		graphFile.close();
		cout<<"nb arcs"<<edge_array.size()<<endl;
	}else cout << "Unable to open graph file \n";
	//cout<<"nb arcs"<<edge_array.size() << " " << weights.size() << " "<<gDir.num_edges()<<"\n";

	cout<<"recuperation de frequences pour k+1 mers \n";

	std::map<string,int> k1mersFreq;
	string k1mer;
	if(k1mersFile.is_open()){
		while ( getline(k1mersFile,line) ){
			if(line != ""){
				int pos =  line.find(" ");
				k1mer = line.substr(0,pos);
				int freqk1 = atoi(line.substr(pos).c_str());
				k1mersFreq.insert(std::pair<string,int>(k1mer,freqk1));
			}
		}
		k1mersFile.close();
		kmerLength = k1mer.length()-1;
		cout<<"nb k+1mers "<<k1mersFreq.size()<< "kmer length" <<kmerLength<<"\n";
	}else cout << "Unable to open k+1mers file \n";

	if(velvetFiles) {
		cout<<"recuperation de small contigs for Velvet \n";
		retriveVelvetShortContigs(edge_array,contigList,kmerLength);
	}


/***************boost
	// Instantiate the visitor for printing cycles
    // créer un flux de sortie
    std::ostringstream oss;
    // écrire un nombre dans le flux
    oss << kmerLength ;
    // récupérer une chaîne de caractères
    std::string result = oss.str();
string fileC = "cycles"+result;
string cyclesFile1 = fileC+".txt";
ofstream outCycle;
outCycle.open(cyclesFile1.c_str(), ios::out );
   cycle_printer<ostream> vis(outCycle);
    // Use the Tiernan algorithm to visit all cycles, printing them
    // as they are found.
    tiernan_all_cycles(gDir, vis);
outCycle.close();
******/

	cout<<"Debut pb"<<endl;
  ifstream cycles(cyclesFile.c_str());
	std::vector<int> cycle;
	std::vector<string> patternsInCycle;//les motif a retrouver pour un cycle
	std::vector< std::pair<Edge ,Edge> > ioArcs;//les in et out arcs pour les motifs patternsInGenome;
	if (cycles.is_open()){
    getline(cycles,line);
    getline(cycles,line);
		while ( getline(cycles,line) ){
			cycle.clear();
			//1. retrive a cycle
      line = line.substr(line.find(" ")+1);
      cout<<"line:"<<line<<endl;
			cycle = detCycle(line);
			patternsInCycle.clear();
			ioArcs.clear();
			//2. determine the frequency of each node
			std::vector<int> freqCycle;
			for(int i = 0; i<cycle.size();i++){
				freqCycle.push_back(contigFrequency[cycle[i]-1]);
			}
      //cout<<"cycle size:"<<cycle.size()<<endl;
			//3. test if each arc is verified by a k+1 mer
			bool coherentCycle = true;
			int j = 0;
			while (coherentCycle && j<cycle.size()-1 ){
				k1mer = contigList[cycle[j]-1].substr(contigList[cycle[j]-1].length()-kmerLength)+contigList[cycle[j+1]-1][kmerLength-1];
				 if(k1mersFreq.find(k1mer) == k1mersFreq.end()) {
						coherentCycle = false;
						//cout<<"cycle noncouvert par k+1: "<<k1mer<<" "<<cycle[j]<<"->"<<cycle[j+1];
				}
				j++;
			}
			k1mer = contigList[cycle[cycle.size()-1]-1].substr(contigList[cycle[cycle.size()-1]-1].length()-kmerLength)+contigList[cycle[0]-1][kmerLength-1];
			if(k1mersFreq.find(k1mer) == k1mersFreq.end()) {
				coherentCycle = false;
				//cout<<"cycle noncouvert par k+1: "<<k1mer<<" "<<cycle[cycle.size()-1]<<"->"<<cycle[0];
			}


      cout<<"Fin test k+1mers cycle"<<endl;
			//4. retrive the in and out edges that do not belong to the cycle
			std::vector<Edge> edge_in;
			std::vector<Edge> edge_out;
			std::vector<int> freq_in;
			std::vector<int> freq_out;
			if(cycle.size() == 1){//one node in the cycle
					for(int j = 0; j<edge_array.size();j++){
						if((edge_array[j].first == cycle[0])&&(edge_array[j].second != cycle[0])) edge_out.push_back(edge_array[j]);
						else if ((edge_array[j].second == cycle[0])&&(edge_array[j].first != cycle[cycle.size()-1])) edge_in.push_back(edge_array[j]);
					}
			}
			if(cycle.size()>1){//for more than one node
				for(int j = 0; j<edge_array.size();j++){//the edges for the first node
						if((edge_array[j].first == cycle[0])&&(edge_array[j].second != cycle[1])) edge_out.push_back(edge_array[j]);
						else if ((edge_array[j].second == cycle[0])&&(edge_array[j].first != cycle[cycle.size()-1])) edge_in.push_back(edge_array[j]);
				}
				for(int j = 0; j<edge_array.size();j++){//the edges for the last node
						if((edge_array[j].first == cycle[cycle.size()-1])&&(edge_array[j].second != cycle[0])) edge_out.push_back(edge_array[j]);
						else if ((edge_array[j].second == cycle[cycle.size()-1])&&(edge_array[j].first != cycle[cycle.size()-2])) edge_in.push_back(edge_array[j]);
				}

			}
			if(cycle.size()>2){//more than two nodes
				for(int i = 1; i<cycle.size()-1; i++){//the edges for the nodes between the first one and the last one
					for(int j = 0; j<edge_array.size();j++){
						if((edge_array[j].first == cycle[i])&&(edge_array[j].second != cycle[i+1])) edge_out.push_back(edge_array[j]);
						else if ((edge_array[j].second == cycle[i])&&(edge_array[j].first != cycle[i-1])) edge_in.push_back(edge_array[j]);
					}
				}
			}

/*			retrivePatterns(cycle, edge_in, edge_out,patternsInCycle,ioArcs,kmerLength,contigList);
			totalRepToFind = totalRepToFind+patternsInCycle.size();
			if(patternsInCycle.size()>0){
				cout<<"Cycle "<<line<<endl;
				cout<<"REP_TO_FIND: "<<patternsInCycle.size()<<endl;
			}
			if(!coherentCycle){
				for(int p = 0;p<patternsInCycle.size();p++){
					cout<<"REPETITION LOST K+1MER_CYCLE " << patternsInCycle[p]
		<<" "<<ioArcs[p].first.first<<"->"<<ioArcs[p].first.second
 <<" "<<ioArcs[p].second.first<<"->"<<ioArcs[p].second.second
<<endl;
				}
			}*/
		
			if(!coherentCycle) continue;
			//cout<<"Fin recuperation arcs externs"<<endl;
			//5. test if each in and out edge is verified by a k+1 mer and remove it if not
      int i = 0;
			//for each in edge
			while(i<edge_in.size()){
				k1mer = contigList[edge_in[i].first-1].substr(contigList[edge_in[i].first-1].length()-kmerLength)+contigList[edge_in[i].second-1][kmerLength-1];
          if(k1mersFreq.find(k1mer) !=  k1mersFreq.end() && edge_in[i].first !=  edge_in[i].second){
					  float fValeur = k1mersFreq.at(k1mer)/cover;
				    float fDecimal = fValeur-floor(fValeur);
            if (fDecimal< 0.5)
              freq_in.push_back(int(floor(fValeur)*cover));
            else
            	freq_in.push_back(int(ceil(fValeur)*cover));
            i++;
          } else{
				/*		int p = 0;
						while(p<patternsInCycle.size()){
							if(edge_in[i].first == ioArcs[p].first.first && edge_in[i].second == ioArcs[p].first.second){
						cout<<"REPETITION LOST K+1MER_IN_ARC " << patternsInCycle[p]
								<<" "<<ioArcs[p].first.first<<"->"<<ioArcs[p].first.second
 								<<" "<<ioArcs[p].second.first<<"->"<<ioArcs[p].second.second
								<<endl;
						patternsInCycle.erase(patternsInCycle.begin()+p);
						ioArcs.erase(ioArcs.begin()+p);
						}else p++;
						}*/
            edge_in.erase(edge_in.begin()+i);

						
          }
			}
			i = 0;
			//for each out edge
      while( i<edge_out.size()){
				k1mer = contigList[edge_out[i].first-1].substr(contigList[edge_out[i].first-1].length()-kmerLength)+contigList[edge_out[i].second-1][kmerLength-1];
        if((k1mersFreq.find(k1mer) != k1mersFreq.end())&&(edge_out[i].first != edge_out[i].second)){
						float fValeur = k1mersFreq.at(k1mer)/cover;
				    float fDecimal = fValeur-floor(fValeur);
            if (fDecimal< 0.5)
              freq_out.push_back(int(floor(fValeur)*cover));
            else
            	freq_out.push_back(int(ceil(fValeur)*cover));
             i++;
        } else {
					/*	int p = 0;
						while(p<patternsInCycle.size()){
							if(edge_out[i].first == ioArcs[p].second.first && edge_out[i].second == ioArcs[p].second.second){
						cout<<"REPETITION LOST K+1MER_OUT_ARC " << patternsInCycle[p]
								<<" "<<ioArcs[p].first.first<<"->"<<ioArcs[p].first.second
 								<<" "<<ioArcs[p].second.first<<"->"<<ioArcs[p].second.second
								<<endl;
						patternsInCycle.erase(patternsInCycle.begin()+p);
						ioArcs.erase(ioArcs.begin()+p);
						}else p++;
						}*/
         		edge_out.erase(edge_out.begin()+i);
        }
			}

			bool coherentIOEdges = true;
			if(edge_in.size() == 0 || edge_out.size() == 0) coherentIOEdges = false;
      cout<<"Fin verification k+1 mers pour arcs externs"<<endl;
			if(coherentCycle && coherentIOEdges){
					/*retrivePatterns(cycle, edge_in, edge_out,patternsInCycle,ioArcs,kmerLength,contigList);
					totalRepToFind = totalRepToFind+patternsInCycle.size();
					if(patternsInCycle.size()>0){
						cout<<"Cycle "<<line<<endl;
						cout<<"REP_TO_FIND: "<<patternsInCycle.size()<<endl;
					}*/
          cout<<"Remouve sup freq"<<endl;
					removeSuppFreq(cycle, freqCycle, edge_in, edge_out, freq_in, freq_out,cover,kmerLength,contigList);
          //,patternsInCycle,ioArcs);
			}


    cout<<"Fin cycle"<<endl;
		}//end while for cycles
	}//end if file cycle is open
	//cout<<"Total repetitions à trouver: "<<totalRepToFind<<endl;
	//cout<<"Total repetitions trouvées: "<<totalRepFound<<endl;

	return 0;
}
