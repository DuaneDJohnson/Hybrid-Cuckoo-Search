//SCRAPs: A Multicomponent Alloy Structure Design Tool
//A novel Hybrid CUCKOO SEARCH determines combinatorial global optimum
//SuperCell Random APproximates(SCRAPs)for specified point and pair
//correlations in high-entropy alloys having proper distributions.
//Copyright (c) 2020. Duane D. Johnson. All rights reserved.
//Developed by: Duane D. Johnson, Rahul Singh, Prashant Singh, Aayush Sharma, Ganesh Balasubramanian.
//Iowa State University
//www.iastate.edu

//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the &quot;Software&quot;), to deal
//with the Software without restriction, including without limitation the
//rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
//sell copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//* Redistributions of source code must retain the above copyright notice,
//this list of conditions and the following disclaimers.
//* Redistributions in binary form must reproduce the above copyright notice,
//this list of conditions and the following disclaimers in the documentation
//and/or other materials provided with the distribution.
//* Neither the names of Duane D. Johnson, Iowa State University, nor the names
//of its contributors may be used to endorse or promote products derived from
//this Software without specific prior written permission.
//* Citation of the original archival journal should be made as is proper, i.e.
//R. Singh, et al., Nature Computational Science VOL, PAGE (2021) DOI:

//THE SOFTWARE IS PROVIDED &quot;AS IS&quot;, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH
//THE SOFTWARE.


#ifndef readInput_hpp
#define readInput_hpp

#include <stdio.h>
#include <vector>
#include <ostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "Lattice.hpp"

using namespace std;
class readInput{
public:
 // Input Variables
    int typeofatoms;
    int globalIter;
    vector<int> maxAtomNum;
    vector<int> atomNum;
    vector<int> eachElement;
    vector<int> cell_dim;
    vector<double> weights;
    int num_of_nests;
    int MCiter;
    int iterTotal;
    int num_of_shells;
// Lattice Variables
    int totalatoms;
    point3d a, b, c;
    vector<point3d> atoms;

// Functions
    void input();
    readInput();
    void readPOSCAR();
    
    
};


#endif /* readInput_hpp */
