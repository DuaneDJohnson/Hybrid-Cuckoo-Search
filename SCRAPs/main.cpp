/*
SCRAPs: A Multicomponent Alloy Structure Design Tool
A novel Hybrid CUCKOO SEARCH determines combinatorial global optimum
SuperCell Random APproximates(SCRAPs)for specified point and pair
correlations in high-entropy alloys having proper distributions.
Copyright (c) 2020. Duane D. Johnson. All rights reserved.
Developed by: Duane D. Johnson, Rahul Singh, Prashant Singh, Aayush Sharma, Ganesh Balasubramanian.
Iowa State University
www.iastate.edu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the &quot;Software&quot;), to deal
with the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimers.
* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimers in the documentation
and/or other materials provided with the distribution.
* Neither the names of Duane D. Johnson, Iowa State University, nor the names
of its contributors may be used to endorse or promote products derived from
this Software without specific prior written permission.
* Citation of the original archival journal should be made as is proper, i.e.
R. Singh, et al., Nature Computational Science VOL, PAGE (2021) DOI:

THE SOFTWARE IS PROVIDED &quot;AS IS&quot;, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH
THE SOFTWARE.
*/

#include <iostream>
#include <fstream>
#include "Lattice.hpp"
#include <unordered_map>
#include <chrono>
#include <random>
#include "readInput.hpp"
#include <mpi.h>
#include <math.h>


void sortthenests(vector<vector<int>> &_nests, vector<double> &_fitness, int &max_num_neighs)
{
    int numofnests = _nests.size();
    vector<int> temp;
    double tempfitn;
    for(int i = 0; i < numofnests; i++)
    {
        for (int j = i + 1; j < numofnests; j++) {
            if(_fitness[j] < _fitness[i])
            {
                swap(_nests[i], _nests[j]);
                swap(_fitness[i], _fitness[j]);
            }
        }
        
    }
}

void printtoscreen(vector<int> nest, const Lattice &lat, int _typesatoms, vector<int> eachElement)
{
    int numofmaxneighs = lat.maxneigh;
    int numofatoms = nest.size();
    double conc;
    vector<vector<vector<double> > > pairGr(_typesatoms, std::vector<vector <double>> (_typesatoms, std::vector<double> (numofmaxneighs, 0.0)));
    vector<vector<vector<double> > > pairSRO(_typesatoms, std::vector<vector <double>> (_typesatoms, std::vector<double> (numofmaxneighs, 0.0)));
    vector<vector<vector<double> > > countSRO(_typesatoms, std::vector<vector <double>> (_typesatoms, std::vector<double> (numofmaxneighs, 0.0)));
    for(int i = 0; i < numofatoms; i++)
    {
        for(int j = 0; j < numofatoms; j++)
        {
            //cout << nest[i] << " Next " << nest[j] << " NEXT " << lat.neighbours[i][j] << endl;
            if (lat.neighbours[i][j] > 0)
                pairGr[nest[i]][nest[j]][lat.neighbours[i][j] - 1] += 1;
        }
    }
    for(int i = 0; i < nest.size(); i++)
        {
            for(int j = 0; j < nest.size(); j++)
            {
                if (lat.neighbours[i][j] > 0){
                    //if( nest[j] != nest[i])
                        //cout << "The pair Gr[" <<  nest[i] << "][" << nest[j] << "] - > is " << pairGr[nest[i]][nest[j]][lat.neighbours[i][j] - 1] << endl;
                    conc = (double)eachElement[nest[j]]/(double)numofatoms;
                    if( nest[i] == nest[j])
                    {
                        pairSRO[nest[i]][nest[j]][lat.neighbours[i][j] - 1] =	(((double)pairGr[nest[i]][nest[j]][lat.neighbours[i][j] - 1]/lat.maxneighbourspershell[nest[i]][lat.neighbours[i][j] - 1]) - conc*eachElement[nest[j]])/(1 - conc);
                        pairSRO[nest[i]][nest[j]][lat.neighbours[i][j] - 1] += pairSRO[nest[i]][nest[j]][lat.neighbours[i][j] - 1]/eachElement[nest[i]];
                        countSRO[nest[i]][nest[j]][lat.neighbours[i][j] - 1] += 1;
                    }
                    else
                    {
                        pairSRO[nest[i]][nest[j]][lat.neighbours[i][j] - 1] = (double)eachElement[nest[i]] - (pairGr[nest[i]][nest[j]][lat.neighbours[i][j] - 1]/lat.maxneighbourspershell[nest[i]][lat.neighbours[i][j] - 1])/((double)eachElement[nest[j]]/(double)numofatoms);
                        pairSRO[nest[i]][nest[j]][lat.neighbours[i][j] - 1] += pairSRO[nest[i]][nest[j]][lat.neighbours[i][j] - 1]/eachElement[nest[i]];
                        countSRO[nest[i]][nest[j]][lat.neighbours[i][j] - 1] += 1;
                    }
                }
            }
        }
    for (int k = 0; k < numofmaxneighs; k++) {
        cout << " SRO FOR THE SHELL NO - > " << k << endl;
        for(int i = 0; i <  _typesatoms; i++){
            for(int j = 0; j < _typesatoms; j++)
            {
                if(i != j)
                {
                   pairSRO[i][j][k] /= countSRO[i][j][k];
                    cout << "The pair SRO[" << i<<"]["<< j<< "] - > is " << pairSRO[i][j][k] << endl;
                //				cout << pairSRO[i][j] << endl;
                }
            }
        }
    }
    return;
}

//Test for the goodness of fit
double fitnessvalue(vector<int> nest, const Lattice &lat, int _typesatoms, vector<int> eachElement, int &max_num_neighs, vector<double> &weights)
{
    int doublepairs = (int)_typesatoms + (int)_typesatoms*((int)_typesatoms - 1)/2;
    double targetErr = 0.0, Error = 0.0, sqrError = 0.0;
    double totalGr = 0.0;
    double targetSRO[5] = {0.0, 0.0, 0.0, 0.0,0.0};
    int numofmaxneighs = lat.maxneigh;
    //cout << "Max Neighbours are " << numofmaxneighs << endl;
    int numofatoms = nest.size();
    double errorArray[5] = {0.0};
    vector<vector<vector<double> > > pairGr(_typesatoms, std::vector<vector <double>> (_typesatoms, std::vector<double> (numofmaxneighs, 0.0)));
    vector<vector<vector<double> > > avgGr(_typesatoms, std::vector<vector <double>> (_typesatoms, std::vector<double> (numofmaxneighs, 0.0)));
    vector<vector<vector<double> > > countgr(_typesatoms, std::vector<vector <double>> (_typesatoms, std::vector<double> (numofmaxneighs, 0.0)));
    //@vector<double> weights(5, 0.0);
    int center, surr, shellnum;
    for(int i = 0; i < nest.size(); i++)
    {
        int size = lat.atoms_shell_atom.size();
        for(int j = 0; j < max_num_neighs; j++)
        {
           // cout << " The shell - - >" << j << endl;
            shellnum = lat.atoms_shell_atom[i][j].size();
            for(int k = 0; k < shellnum; k++)
            {
                center = nest[i];
                surr = nest[lat.atoms_shell_atom[i][j][k]];
                pairGr[center][surr][j] += 1;
                if(nest[center] != nest[surr] )
                {
                    avgGr[center][surr][j] += (double)eachElement[center]*(1 - targetSRO[j])*shellnum*((double)eachElement[surr]/(double)numofatoms);
                    countgr[center][surr][j] += 1;
            
                }
            }
            
        }
    }
    
   // cout << "   pair GER stsrts" << endl;
    for(int k = 0; k < max_num_neighs; k++)
    {
        for(int i = 0; i < _typesatoms; i++)
        {
            for(int j = 0; j < _typesatoms; j++)
            {
                if(i != j )
                {
                    if(countgr[i][j][k] == 0)
                        avgGr[i][j][k] = 0;
                    else
                        avgGr[i][j][k] /= countgr[i][j][k];
                    errorArray[k] += pow(pairGr[i][j][k] - avgGr[i][j][k], 2);
                }
            }
        }
        Error += weights[k]*errorArray[k];
    }
    sqrError = sqrt(Error);
//    cout << Error << endl;
    return sqrError;
}

/*
 
 The function based on fitness randomize and sorts the best nests for global optimization.
 
 */
int get_best_nest(vector<vector<int>> &_nests, vector<vector<int>> &_newnests, vector<double> &_fitness, const Lattice &lat, int _typeatoms, vector<int> eachElement, int &max_num_neighs, vector<double> &weights)
{
    int n, indexOfBest;
    n = _nests.size();
    double keepCuckoos = 0.4;
    sortthenests(_nests, _fitness, max_num_neighs);
    if(n == 1)
    {
        int fitness_new = fitnessvalue(_newnests[0], lat, _typeatoms, eachElement, max_num_neighs, weights);
        if (fitness_new < _fitness[0])
        {
            _nests = _newnests;
            _fitness[0] = fitness_new;
        }
    }
    else
    {
        for(int i = (int)(keepCuckoos*(double)n); i < n; i++)
        {
            _nests[i] = _newnests[i - (int)(keepCuckoos*(double)n)];
            _fitness[i] = fitnessvalue(_nests[i], lat, _typeatoms, eachElement, max_num_neighs, weights);
        }
        sortthenests(_nests, _fitness, max_num_neighs);
    }
    cout << "===================  BEST NEST  =========================" << endl;
    cout << "            	Fitness ::" << _fitness[0] << endl;
    cout << endl;
    indexOfBest = 0;
    return indexOfBest;
}
/*
 
 Based on goodness of "Fitness", function below discards the old and add new ones
during global optimization.

 */

int empty_nests(vector<vector<int>> &_nests, vector<vector<int>> &_newnests, int iter, vector<double> &origfitness , const Lattice &lat, const int _typesatoms, const vector<int> eachElement, int &max_num_neighs, vector<double> &weights)
{
    double len = _nests[0].size();
    double eps = 0.001;
    std::random_device rd, rd1;
    std::mt19937 gen(rd());
    std::mt19937 gen1(rd1());
    std::uniform_int_distribution<> dis(0, len - 1);
    std::uniform_int_distribution<> dis1(1, 100);
    double oldFitness, newFitness;
    
    int atom1, atom2, temphold, prob;
    vector<int> prevnest, currnest;
    for (int i = 0; i < _nests.size(); i++) {
        prevnest = _nests[i];
        currnest = _nests[i];
        oldFitness = fitnessvalue(_nests[i], lat, _typesatoms, eachElement, max_num_neighs, weights);
        origfitness[i] = oldFitness;
        for(int j = 0; j < iter; j++)
        {
            atom1 = (int) dis(gen);
            atom2 = (int) dis(gen);
            if(currnest[atom1] == currnest[atom2])
                continue;
            temphold = currnest[atom1];
            currnest[atom1] = currnest[atom2];
            currnest[atom2] = temphold;
            newFitness = fitnessvalue(currnest, lat, _typesatoms, eachElement, max_num_neighs, weights);
            if(newFitness < oldFitness)
            {
                prevnest = currnest;
                origfitness[i] = newFitness;
                oldFitness = newFitness;
            }
            else
            {
                if (abs(newFitness - oldFitness) < eps) {
                    prob = dis1(gen1);
                    if(prob > 50)
                    {
                        prevnest = currnest;
                        origfitness[i] = newFitness;
                        oldFitness = newFitness;
                    }
                }
                else
                    currnest = prevnest;
            }
        }
        _nests[i] = currnest;
    }
    return 0;
}
vector<double> fitnessofnests(vector<vector<int>> nests, const Lattice &lat, int _typesatoms, vector<int> eachElement, int &max_num_neighs, vector<double> &weights)
{
    vector<double> fitness(nests.size(), INT_MAX);
    //cout << "SIZE OF NESTS " << nests.size() << endl;
    //cout << "SIZE OF NESTS {0} - > " << nests[0].size() << endl;
    for(std::size_t i = 0; i < nests.size(); i++)
    {
        fitness[i] = fitnessvalue(nests[i], lat, _typesatoms, eachElement, max_num_neighs, weights);
        //cout << fitness[i] << endl;
    }
    return fitness;
}

vector<vector<int>> randomizeLattice(int numofnests, int _typeofatoms, int _size, vector<int> _maxAtomNum, int numofatoms)
{
// Random Number Generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, _typeofatoms - 1);
    vector<vector<int>> nests(numofnests, vector<int>(numofatoms, 0));
    vector<int> _numofatoms(_size, 0);
    vector<int> atomNum(_typeofatoms, 0);
    int i = 0, atom;
    for (int j = 0; j < numofnests; j++) {
        while( i < _size)
        {
            atom = (int)dis(gen);
            if(atomNum[atom] < _maxAtomNum[atom])
            {
                _numofatoms[i] = atom;
                i++;
                atomNum[atom] += 1;
            }
            else
                continue;
        
        }
        nests[j] = _numofatoms;
        i = 0;
        std::fill(atomNum.begin(), atomNum.end(), 0);
    }
    return nests;
}

/*
 
 This par below oldbased on fitness randomize and sorts the best nests for global optimization.
.
 
 */

int get_cuckoos(vector<vector<int>> &_nests, vector<vector<int>> &_newnests, int iter, int MCiter, const Lattice &lat, int _typesatoms, const vector<int> &eachElement, const vector<int> &maxAtomNum, int &max_num_neighs, vector<double> &weights, int nproc)
{
    int numofatoms = _nests[0].size();
    int numofnests = _nests.size();
    int localiter = (36 + numofnests)/(int)nproc; // This decides the number of MC steps for local optimization
    vector<vector<int> > localnests(randomizeLattice(localiter, _typesatoms, numofatoms, maxAtomNum, numofatoms));
    vector<double> localfitness(localiter, INT_MAX);
    vector<double> fitness(numofnests, INT_MAX);
    empty_nests(localnests, localnests, 8*MCiter, localfitness, lat, _typesatoms, eachElement, max_num_neighs, weights);
    sortthenests(localnests, localfitness, max_num_neighs);
    std::copy ( localnests.begin(), localnests.begin() + numofnests, _newnests.begin());
    empty_nests(_newnests, _newnests, iter, fitness, lat, _typesatoms, eachElement, max_num_neighs, weights);
    sortthenests(_newnests, fitness, max_num_neighs);
    cout << " NEW NESTS FITNESS =<----1----> " << fitness[0] << " <----2----> " << fitness[1] << " <----3----> " << fitness[2] << endl;
    return 0;
}
void printFinalCoords2File(vector<int> nest, const Lattice &lat, int _typesofatoms, vector<int> _eachElements)
{
    ofstream finalcoords;
    vector<string> coords(_typesofatoms, "");
    int numofatoms = nest.size();
    finalcoords.open("./SCRAPS.vasp", ios::out); //Calculates and writes the final SCRAPS supercell using generalized Hybrid Cuckoo program
    if(finalcoords.is_open()){
        finalcoords << "Generated-Hybrid-Cuckoo-Search" << "\n" << "1.0" << "\n" ;
        finalcoords << "\t" << lat.celldim[0]*lat.basisset.a.x << "\t" << lat.celldim[0]*lat.basisset.a.y << "\t" << lat.celldim[0]*lat.basisset.a.z << "\n" ;
        finalcoords << "\t" << lat.celldim[1]*lat.basisset.b.x << "\t" << lat.celldim[1]*lat.basisset.b.y << "\t" << lat.celldim[1]*lat.basisset.b.z << "\n" ;
        finalcoords << "\t" << lat.celldim[2]*lat.basisset.c.x << "\t" << lat.celldim[2]*lat.basisset.c.y << "\t" << lat.celldim[2]*lat.basisset.c.z << "\n" ;
        for(int i = 0; i < _typesofatoms; i++)
        {
            finalcoords << _eachElements[i] << "   ";
        }
        finalcoords << "\n" << "Cartesian" << "\n";
        for(int i = 0; i < numofatoms; i++){
            coords[nest[i]] += to_string(lat.supercell[i].x) + " \t"  +  to_string(lat.supercell[i].y) + "\t" + to_string(lat.supercell[i].z) + "\n";
            cout << "   " << lat.supercell[i].x << "  " << lat.supercell[i].y << "   " << lat.supercell[i].z << endl;
        }
        for(int i = 0; i < _typesofatoms; i++)
        {
            finalcoords << coords[i];
        }
        
    }
    else
        cout << "Unable to open file";
    
}
/*
Parallaization is done over numbers of nests during global optimization to 
speed up the SCRAPs generation
*/
int main(int argc, char ** argv)
{
    auto started = std::chrono::high_resolution_clock::now();
    readInput inputfile;
    int id, nproc;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);  // get number of total nodes
    MPI_Comm_rank(MPI_COMM_WORLD, &id); // get id of mynode


//  *****READ INPUT FROM "INPUT" and "POSCAR" FILES******
    
    //  ###### INITIALIZE THE SCRAPS GENERATION######
    int iters = 0;
    double error = INT_MAX, epsi = 0.001;
    int globalIter = inputfile.globalIter;
    int num_of_nests = inputfile.num_of_nests;
    int MCiter = inputfile.MCiter;
    int iterTotal = inputfile.iterTotal;

    //####    LATTICE INPUT #######
    
    point3d a = inputfile.a;
    point3d b = inputfile.b;
    point3d c = inputfile.c;
    latstruct basis(a, b, c);
    vector<int>dim = inputfile.cell_dim;
    vector<point3d> atoms = inputfile.atoms;
    int typeofatoms = inputfile.typeofatoms;
    vector<int> eachElement = inputfile.eachElement;
    vector<double> weights = inputfile.weights; // Weights can be defined later on:
    int max_num_neighs = inputfile.num_of_shells; // Maximum number of neighbours can be changed here
    
    //point3d atom1(-0.000001613, 1.838152961, 1.312359929), atom2(1.591887967 , 0.919075084, 3.937079787);
    //vector<point3d> atoms  = inputfile.atoms;
    

    cout << "Creating Lattice....." << endl;
    Lattice templat(basis, dim, atoms);
    cout << "Lattice Built....." << endl;
    
// Switch to maximum shell number if exceeds MAXIMUM allowed in given lattice
    if(templat.maxneigh < max_num_neighs)
        max_num_neighs = templat.maxneigh;
    if(weights.size() < max_num_neighs)
    {
        cout << "WEIGHTS SIZE NOT EQUAL TO SHELL NUMBER\n STOPPING......." << endl;
        
    }
    
    vector<vector<double>> distances = templat.distanceMatrix();
    vector<vector<int>> neighs = templat.neighbours;
    unordered_map<int, int> numneighs;
    int total_atoms = atoms.size()*dim[0]*dim[1]*dim[2];

    vector<int> maxAtomNum(typeofatoms, 0);
    vector<int> atomNum(typeofatoms, 0);
    
    std::uniform_int_distribution<> dis(0, typeofatoms - 1);
    for(int i = 0; i < typeofatoms; i++)
    {
        atomNum[i] = 0;
        maxAtomNum[i] = eachElement[i];
    }
    
// Global Optimization - The Cuckoo search algorithm starts here
// Set the number of nests
     cout << "Starting Optimization ....." << endl;
    double fmin, fmax, fnew, pa;
    int best_max, best_new;
    int startval = num_of_nests*id/nproc + 1;
    int endval = num_of_nests*(id + 1)/nproc;
    int nestsperproc = endval - startval + 1;
    vector<vector<int>> nests(randomizeLattice(nestsperproc, typeofatoms, total_atoms, maxAtomNum, total_atoms));
    vector<vector<int >> new_nests(randomizeLattice(nestsperproc, typeofatoms, total_atoms, maxAtomNum, total_atoms));
    vector<double> fitness(nestsperproc, INT_MAX);
    fitness = fitnessofnests(nests, templat, typeofatoms, eachElement, max_num_neighs, weights);
    vector<int> best(total_atoms, 0), bestnest(total_atoms, 0);
    best_max = get_best_nest(nests, new_nests, fitness, templat, typeofatoms, eachElement, max_num_neighs, weights);

//    empty_nests(nests, new_nests, MCiter, fitness, templat, typeofatoms, eachElement, max_num_neighs, weights);
    while (iters < iterTotal && error > epsi)
    {
// The Cuckoo search algorithm: set up cuckoo search, finds best nests, do fitness test until fins the best nest
        get_cuckoos(nests, new_nests, globalIter, MCiter, templat, typeofatoms, eachElement, maxAtomNum, max_num_neighs, weights, nproc);
        fitness = fitnessofnests(nests, templat, typeofatoms, eachElement, max_num_neighs, weights);
        best_new = get_best_nest(nests, new_nests, fitness, templat, typeofatoms, eachElement, max_num_neighs, weights);
        fitness = fitnessofnests(nests, templat, typeofatoms, eachElement, max_num_neighs, weights);
        empty_nests(nests, new_nests, MCiter, fitness, templat, typeofatoms, eachElement, max_num_neighs, weights);
        fitness = fitnessofnests(nests, templat, typeofatoms, eachElement, max_num_neighs, weights);
        globalIter = globalIter + (iters+1)*MCiter;
        sortthenests(nests, fitness, max_num_neighs);
        error = fitness[0];
        cout << " FITNESS - > " << fitness[0] << "<-- iter -->" << iters << endl;
        iters += 1;

    }
    vector<int> nest1D;
    if (id != 0) {
        for (int i = 0; i < nestsperproc; i++)
        {
            for(int j = 0; j < total_atoms; j++)
                nest1D.push_back(nests[i][j]);
        }
        MPI_Send(&nest1D[0], total_atoms*nestsperproc, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
        for (int j = 1; j < nproc; j++) {
            startval = num_of_nests*j/nproc + 1;
            endval = num_of_nests*(j + 1)/nproc;
            nestsperproc = endval - startval + 1;
            nest1D.resize(total_atoms*nestsperproc);
            MPI_Recv(&nest1D[0], total_atoms*nestsperproc, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
            for (int i = 0; i < nestsperproc; i++) {
                vector<int> temp(nest1D.begin() + i*total_atoms, nest1D.begin() + (i + 1)*total_atoms);
                nests.push_back(temp);
            }
        }
        fitness = fitnessofnests(nests, templat, typeofatoms, eachElement, max_num_neighs, weights);
        sortthenests(nests, fitness, max_num_neighs);
        printtoscreen(nests[0], templat, typeofatoms, eachElement);
        printFinalCoords2File(nests[0], templat, typeofatoms, eachElement);
    //cout << "HELOO" << endl;
        auto done = std::chrono::high_resolution_clock::now();
        std::cout << "TIME TAKEN ->  ";
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count();
        std::cout << "ms \n" << endl;
    }
    MPI_Finalize();
}
/**
 Generate new random structures
 **/
