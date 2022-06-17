#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>

#define MU 0.33
#define DELTA  1.33

float  find_array_max(float array[], int length);

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

float max_array(float array[],int length) {
    float maxx = array[0];
    int index = 0;
    for(int i = 1; i<length; i++) 
    {
        if(array[i] > maxx) 
        {
            maxx = array[i];
            index = i;
        }
    }
    return index;           
}


int main(int argc, char** argv) {
    
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


    int Na = seqA.length();
    int Nb = seqB.length();

    float H[Na + 1][Nb + 1] = {0.};
    unsigned short I_i[Na + 1][Nb + 1];
    unsigned short I_j[Na + 1][Nb + 1];


    //Start timer
    struct timeval	StartTime, EndTime;
    gettimeofday(&StartTime, NULL);

    //Algorithm
    float temp[4];
    int index;
    for(int i=1; i<=Na; i++) {
        for(int j=1; j<=Nb; j++) {
            temp[0] = H[i-1][j-1]+score(seqA[i-1],seqB[j-1]);
            temp[1] = H[i-1][j]-DELTA;
            temp[2] = H[i][j-1]-DELTA;
            temp[3] = 0.;
            index = max_array(temp, 4);
            H[i][j] = temp[index];

            if (index==0)
            {
                I_i[i][j] = i-1;
                I_j[i][j] = j-1;
            }
            else if (index==1)
            {
                I_i[i][j] = i-1;
                I_j[i][j] = j;
            }
            else if (index==2)
            {
                I_i[i][j] = i;
                I_j[i][j] = j-1;
            }
            else
            {
                I_i[i][j] = i;
                I_j[i][j] = j;
            }
        }
    }


    float H_max = 0.;
    int i_max=0,j_max=0;
    for(int i=1; i<=Na; i++) {
        for(int j=1; j<=Nb; j++) {
            if(H[i][j]>H_max) {
                H_max = H[i][j];
                i_max = i;
                j_max = j;
            }
        }
    }

    //Backtracking
    int i_temp = i_max;
    int j_temp = j_max;
    int i_next = I_i[i_temp][j_temp];
    int j_next = I_j[i_temp][j_temp];
    char matchA[Na+Nb+2], matchB[Na+Nb+2];

    int iter = 0;
    while(((i_temp!=i_next) || (j_temp!=j_next)) && (j_next!=0) && (i_next!=0)) 
    {

        if(i_next==i_temp)
        {  
            matchA[iter] = '-';
        }
        else    
        {
            matchA[iter] = seqA[i_temp-1];
        }

        if(j_next==j_temp)
        {  
            matchB[iter] = '-';
        }
        else    
        {
            matchB[iter] = seqB[j_temp-1];
        }

        i_temp = i_next;
        j_temp = j_next;
        i_next = I_i[i_temp][j_temp];
        j_next = I_j[i_temp][j_temp];
        iter++;
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
}