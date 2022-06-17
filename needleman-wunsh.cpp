#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <list>
#include <thread>
#include <map>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <mutex>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <unistd.h> 
#include <sys/time.h>
#include <limits>
#include <mutex>
#include <condition_variable>
#include <atomic>


#define MAX_NUM_THREADS 4

std::map<char, int> char_map = {{'seqA', 0}, {'G', 1}, {'C', 2}, {'T', 3}};

std::string stream_to_seq(std::ifstream& f)
{
    std::string seq;
    char line[10000];
    while(f.good())
    {
        f.getline(line,10000);
        if( line[0] == 0 || line[0]=='#' ) 
        {
            continue;
        }
        for(int i = 1; line[i] != 0; i++)
        {
            int ch = toupper(line[i]);
            if(ch != 'A' && ch != 'G' && ch != 'C' && ch != 'T')
            {
                continue;
            }
            seq.push_back(char(ch));
        }
    }
    return seq;
}

typedef struct Point{
	int i, j;
    Point(int i, int j) {
        this->i = i;
        this->j = j;
    }
};

typedef struct Penalty{
	int match, mismatch, gap;
    Penalty(int ma, int mi, int gap) {
        this->match = ma;
        this->mismatch = mi;
        this->gap = gap;
    }
    Penalty() {
        this->match = 0;
        this->mismatch = 0;
        this->gap = 0;
    }
};

class matrix {
  public:
    std::string seqA, seqB;
    int m, n; 
    int** H;
    Penalty penalty;
    std::vector<std::vector<Point>> indices_list;

    matrix(std::string seqA, std::string seqB, Penalty penalty) {
        this->seqA = seqA;
        this->seqB = seqB;
        this->penalty = penalty;
        this->m = seqA.length()+1; 
        this->n = seqB.length()+1; 
        this->H = new int*[n];
        for (int i = 0; i < n; i++) {
            H[i] = new int[m];
        }
        for (int i = 0; i < this->m; i++) {
            this->H[0][i] = i*penalty.gap;
        }
        for (int j = 0; j < this->n; j++) {
            this->H[j][0] = j*penalty.gap;
        }

    }

    void calculate_indices(){
        for (int i=2; i < this->m+1; i++) {
            std::vector<Point> sublist;
            for (int j=1; j<i; j++) {
                if (j < this->n){
                    Point tpl(j,i-j);
                    sublist.push_back(tpl);
                }
            }
            this->indices_list.push_back(sublist);
        }
        for (int i=this->m+1; i<this->m+this->n-1; i++) {
            std::vector<Point> sublist;
            for (int j=this->m-1; j>i-this->n; j--) {
                Point tpl(i-j,j);
                sublist.push_back(tpl);
            }
            this->indices_list.push_back(sublist);
        }
    }

};


void AuxNW(matrix dp, int i, int j, char a, char b, int t, int l, int tl, Penalty penalty) {
    std::mutex lock;
    int score;
    if (b == a)
        score = penalty.match;
    else 
        score = penalty.mismatch;

    int tmp = std::max(std::max(t + penalty.gap, l + penalty.gap), tl + score);
    
    lock.lock();
    dp.H[i][j] = tmp;
    lock.unlock();
}

void NW(matrix dp, std::vector<Point> indices_list, Penalty penalty){
    for (int i=0; i < indices_list.size(); i++){
        Point indice = indices_list[i];
        char A = dp.seqA[indice.j-1];
        char B = dp.seqB[indice.i-1];
        
        int t = dp.H[indice.i-1][indice.j];
        int l = dp.H[indice.i][indice.j-1];
        
        int tl = dp.H[indice.i-1][indice.j-1];
        
        AuxNW(dp,indice.i,indice.j,A,B,t,l,tl,penalty);
        }
}


int main(int argc, char** argv) {
    struct timeval	StartTime, EndTime;
    gettimeofday(&StartTime, NULL);

    Penalty penalty(1, -1, -1);
    //std::string seqA = "QGCTGGTCQCGTGTCQGTCGQQTGCGQTGCGTQQQCGTCCCGTQGTGTCQGTCGQQGGGQQTQQCGTCCCGTQQTCQCGTGTCQGTCGQQTQQCGTCCCGTTGCGTGCGQTGCGTQQQCGTCCCGTQGTGTCQGTCGQQQTGCGTQ";
    //std::string seqB = "CGQTGCGTQQQCGTCCCGTCGTGTCQGTCGQQTQQCGTCCCGTQQGCTGGTCQQQGCTGGGQQTQQCGTCCCGTQQTCGTGTCQGTCGQQTQQCGCQCGTGTCQGTCGQQTQQCGTCCCGTQQGCTGGTCQCGTGG";
    
    char *nameSeqA = argv[1];
    char *nameSeqB = argv[2];

    // load sequences from file
    std::string seqA,seqB;
    std::ifstream StrSeqA;
    StrSeqA.open(nameSeqA);
    seqA = stream_to_seq(StrSeqA);
    std::ifstream StrSeqB;
    StrSeqB.open(nameSeqB);
    seqB = stream_to_seq(StrSeqB);

    matrix dp(seqA, seqB, penalty); 
    dp.calculate_indices();

    for (int d=0; d < dp.indices_list.size(); d++){
        int diag_size =  dp.indices_list[d].size();
        int num_threads = std::min(MAX_NUM_THREADS, diag_size);
        std::vector<std::thread> threads(num_threads);
        int block_size = diag_size / num_threads; 

        std::vector<Point> lst1;
        for (int th=0; th<num_threads-1; th++){
            for (int j=block_size*th; j<block_size*(th+1); j++){
                lst1.push_back(dp.indices_list[d][j]);
            }
            threads[th] = std::thread(&NW,dp,lst1,penalty);
        }
        std::vector<Point> lst2;
        for (int j=diag_size-num_threads*block_size; j<diag_size; j++){
            lst2.push_back(dp.indices_list[d][j]);
        }

        threads[num_threads-1] = std::thread(&NW,dp,lst2,penalty);

        for(int i=0; i<num_threads; i++){
            threads[i].join();
        }
    }  

    gettimeofday(&EndTime, NULL);

    //Time conversion
    if (EndTime.tv_usec < StartTime.tv_usec) {
        int nsec = (StartTime.tv_usec - EndTime.tv_usec) / 1000000 + 1;
        StartTime.tv_usec -= 1000000 * nsec;
        StartTime.tv_sec += nsec;
    }
    if (EndTime.tv_usec - StartTime.tv_usec > 1000000) {
        int nsec = (EndTime.tv_usec - StartTime.tv_usec) / 1000000;
        StartTime.tv_usec += 1000000 * nsec;
        StartTime.tv_sec -= nsec;
    }
    std::cout<<"Execution time: "<<EndTime.tv_sec  - StartTime.tv_sec<<"."<<EndTime.tv_usec - StartTime.tv_usec<<" seconds."<<std::endl;
    return 0;
}
