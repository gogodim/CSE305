#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <unistd.h> 
#include <sys/time.h>
#include <iostream>
#include <omp.h>
#include <x86intrin.h>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>
#include <vector>
#include <limits>

#define MU 0.33
#define DELTA  1.33
#define BLOCK_W 4
#define BLOCK_H 64

//Global variables
std::string seqA,seqB;
float **H;
int  Na, Nb;

// Block sizes
int Bw, Bh;
const int minBw = 64;
const int minBh = 4;


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


float score(char a,char b) 
{
    float res;
    if(a==b) {
        res=1.;
    }
    else {
        res=-MU;
    }
    return res;
}


float max_array_zero(float array[],int length) {
    float maxx = 0;
    for(int i = 0; i<length; i++) 
    {
        if(array[i] > maxx) 
        {
            maxx = array[i];
        }
    }
    return maxx;           
}

typedef struct {
	int i, j;
} Point;

Point next_point(Point pnt) {
	
	int row = pnt.i;
	int	col = pnt.j;
	
	float Hij = H[row][col];
	
	if (Hij == H[row-1][col-1] + score(seqA[row-1], seqB[col-1])) {
		--pnt.i;
		--pnt.j;
	} 
	else if ( Hij == H[row][col-1] - DELTA) {
		--pnt.j;
	}
	else if ( Hij == H[row-1][col] - DELTA ) {
		--pnt.i;
	}
	return pnt;
}

void block(Point block) {
	int i, j;
	int row = block.i * Bh;
	int col = block.j * Bw; 
	int rows_len = std::min(Na, row + Bh);
	int cols_len = std::min(Nb, col + Bw);	
    
    float tmp[3];
	for (i = row + 1; i <= rows_len; ++i) 
	{	
		float temp = H[i][col];

		for (j = col + 1; j <= cols_len; ++j)
		{		
			tmp[0] = H[i-1][j-1] + score(seqA[i-1], seqB[j-1]);
			tmp[1] = H[i-1][ j ] - DELTA;
			tmp[2] = temp - DELTA;

			H[i][j] = temp = max_array_zero(tmp,3);
		}
	}
}

int main (int argc, char** argv) {
//int main () {
    char *nameSeqA = argv[1];
    char *nameSeqB = argv[2];

    //load sequences from file
    std::ifstream StrSeqA;
    StrSeqA.open(nameSeqA);
    seqA = stream_to_seq(StrSeqA);
    std::ifstream StrSeqB;
    StrSeqB.open(nameSeqB);
    seqB = stream_to_seq(StrSeqB);

    // seqA = "QGCTGGTCQCGTGTCQGTCGQQTGCGQTGCGTQQQCGTCCCGTQGTGTCQGTCGQQGGGQQTQQCGTCCCGTQQTCQCGTGTCQGTCGQQTQQCGTCCCGTTGCGTGCGQTGCGTQQQCGTCCCGTQGTGTCQGTCGQQQTGCGTQ";
    // seqB = "CGQTGCGTQQQCGTCCCGTCGTGTCQGTCGQQTQQCGTCCCGTQQGCTGGTCQQQGCTGGGQQTQQCGTCCCGTQQTCGTGTCQGTCGQQTQQCGCQCGTGTCQGTCGQQTQQCGTCCCGTQQGCTGGTCQCGTGG";
    int Na = seqA.length();
    int Nb = seqB.length();

    float** H = new float*[Na + 1];
    for(int i = 0; i < Na + 1; ++i)
        H[i] = new float[Nb + 1];

			
	Point upper = { 0, 0 }, 
		  lower = upper, tmp;	
	
	int newRows, newCols, numPhases;
	
			
    //Start timer
    struct timeval	StartTime, EndTime;
    gettimeofday(&StartTime, NULL);
    
    
    Point max_pos;
    float max_value = -std::numeric_limits<float>::infinity();
    
    #pragma omp parallel
    {
        #pragma omp single
        {
            char thr = omp_get_num_threads();
            Bh = (int)ceil(((double)Na)/(thr*BLOCK_H));
            if (Bh < minBh) Bh = minBh;
            newRows = (int)ceil(((double)Na) / Bh);

            Bw = (int)((double)Nb)/(thr*BLOCK_W);
            if (Bw < minBw) Bw = minBw;
            newCols = (int)ceil(((double)Nb) / Bw);
            
            numPhases = newRows + newCols - 1;
            
            // Create all blocks
            for (int i = 0; i < numPhases; i++) {
				tmp = lower;
				while (tmp.i >= upper.i && tmp.j <= upper.j) {
                    // create tasks for each block
					#pragma omp task firstprivate(tmp)
					block(tmp);
					--tmp.i; 
					++tmp.j;
				}
				#pragma omp taskwait
				// Wait for the whole phase (diagonal) to be done
				if (upper.j < newCols - 1) ++upper.j; else ++upper.i;
				if (lower.i < newRows - 1) ++lower.i; else ++lower.j;
			}
        }

        for (int i = 1; i <= Na; ++i) {
            for (int j = 1; j <= Nb; ++j) {
                if( H[i][j] > max_value ) {
                    max_value = H[i][j];
                    max_pos.i = i;
                    max_pos.j = j;
                }
            }
        }
    }

    // Backtracking
    Point curr = { max_pos.i, max_pos.j };
	Point next = next_point(curr);
    
    char matchA[Na+Nb+2];
    char matchB[Na+Nb+2];
    
    int iter = 0;
	while (((curr.i!=next.i) || (curr.j!=next.j)) && (next.i!=0) && (next.j!=0)) {
				
		if ( next.i == curr.i ) {
			matchA[iter] = '-';                  					
		} else {
			matchA[iter] = seqA[curr.i-1];   					
		}
		
		if ( next.j == curr.j ) { 
			matchB[iter] = '-';                  					
        } else {
			matchB[iter] = seqB[curr.j-1];  
		}
	
		++iter;
		curr = next;
		next = next_point(next);    
    }

    //End timer
    gettimeofday(&EndTime, NULL);

    std::cout<<"SeqA: ";
    for(int i=0; i<Na; i++) {
        std::cout<<seqA[i];
    };
    std::cout<<std::endl<<"SeqB: ";
    for(int i=0; i<Nb; i++) {
        std::cout<<seqB[i];
    };
    std::cout<<std::endl<<"Allignment: ";
    
    for(int i=iter-1; i>=0; i--){
        std::cout<<matchA[i];
    } 
    std::cout<<std::endl;
        for(int j=iter-1; j>=0; j--){
        std::cout<<matchB[j];
    } 
    std::cout<<std::endl;

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