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

#include "Lattice.hpp"
#include <Eigen/Dense>
#include <queue>
#include <set>
#include <unordered_map>

using namespace std;
point3d& point3d::operator=(const point3d &vec2)
{
    x  = vec2.x;
    y = vec2.y;
    z = vec2.z;
    return *this;
}
point3d operator+(const point3d &vec1,const point3d &vec2)
{
    double _x = vec1.x + vec2.x;
    double _y = vec1.y + vec2.y;
    double _z = vec1.z + vec2.z;
    return point3d(_x, _y, _z);
}
point3d operator*(const point3d &vec1,const point3d &vec2)
{
    double _x = vec1.x * vec2.x;
    double _y = vec1.y * vec2.y;
    double _z = vec1.z * vec2.z;
    return point3d(_x, _y, _z);
}
point3d operator*(const point3d &vec1, int m)
{
    double _x = m*vec1.x;
    double _y = m*vec1.y;
    double _z = m*vec1.z;
    return point3d(_x, _y, _z);
}
point3d operator*(int m,const point3d &vec2)
{
    return vec2*m;
}
point3d operator-(const point3d &vec1,const point3d &vec2){
    double _x, _y, _z;
    _x = vec1.x - vec2.x;
    _y = vec1.y - vec2.y;
    _z = vec1.z - vec2.z;
    return point3d(_x, _y, _z);
}
point3d cross(const point3d &vec1,const point3d &vec2)
{
    double _x, _y, _z;
    _x = vec1.y*vec2.z - vec1.z*vec2.y;
    _y = vec1.x*vec2.z - vec1.z*vec2.x;
    _z = vec1.x*vec2.y - vec1.y*vec2.x;
    return point3d(_x, _y, _z);
}
double mod(const point3d vec1)
{
    double result;
    result = sqrt(vec1.x*vec1.x + vec1.y*vec1.y + vec1.z*vec1.z);
    return result;

}

double dot(const point3d &vec1, const point3d &vec2){
    double _x, _y, _z, result;
    _x = vec1.x*vec2.x;
    _y = vec1.y*vec2.y;
    _z = vec1.z*vec2.z;
    result = _x + _y + _z;
    return result;
}
double distance(const point3d &vec1,const point3d &vec2){
    double _x, _y, _z, dist;
    _x = vec1.x - vec2.x;
    _y = vec1.y - vec2.y;
    _z = vec1.z - vec2.z;
    dist = sqrt(pow(_x, 2) + pow(_y, 2) + pow(_z, 2));
    return dist;
}


Lattice::Lattice(latstruct _basisset, vector<int> dim, vector<struct point3d> _atoms)
{
    basisset = _basisset;
    celldim = dim;
    atoms = _atoms;
    initialize();
    
}

void Lattice::initialize(){
        int count = 1;
        int x_dir = celldim[0];
        int y_dir = celldim[1];
        int z_dir = celldim[2];
        vector<point3d> comp_strcut;
        
        point3d temppoint;
        for(point3d atom:atoms)
        {
            for(int k = 0; k < z_dir; k++)
            {
                for (int j = 0; j < y_dir; j++)
                {
                    for(int i = 0; i < x_dir; i++)
                    {
                        temppoint = atom + i*basisset.a + j*basisset.b + k*basisset.c;
                        //cout << count++ << "  " << temppoint.x << " " << temppoint.y << "  " << temppoint.z << endl;
                        comp_strcut.push_back(temppoint);
                    }
                }
                
            }
        }
        
        supercell = comp_strcut;
        int latticeSites = comp_strcut.size();
        point3d latvec1 = (x_dir)*basisset.a;
        point3d latvec2 = (y_dir)*basisset.b;
        point3d latvec3 = (z_dir)*basisset.c;
        real << latvec1.x, latvec2.x, latvec3.x,
                latvec1.y, latvec2.y, latvec3.y,
                latvec1.z, latvec2.z, latvec3.z;
        frac = real.inverse();
        distances = distanceMatrix();
        neighbours = neighbourMatrix();
        maxneighbourspershell = numneighofeachshell();
        maxneighshell();
    }
Lattice operator+(const Lattice &lat1, const Lattice &lat2)
{
    vector<point3d> atoms1 = lat1.atoms;
    vector<point3d> atoms2 = lat2.atoms;
    atoms1.insert(atoms1.end(), atoms2.begin(), atoms2.end());
    Lattice newlat(lat1.basisset, lat1.celldim, atoms1);
    return newlat;
}

double Lattice::imageDistance(point3d latpoint1, point3d latpoint2)
{
    double realdis;
    Eigen::Vector3d point1, point2;
    point1 << latpoint1.x, latpoint1.y, latpoint1.z;
    point2 << latpoint2.x, latpoint2.y, latpoint2.z;
    Eigen::Vector3d fracpoint1, fracpoint2, deltafrac12, realdisp;
    fracpoint1 = frac*point1;
    fracpoint2 = frac*point2;
    deltafrac12 = fracpoint1 - fracpoint2;
    deltafrac12 = deltafrac12 - Eigen::Vector3d(round(deltafrac12(0)), round(deltafrac12(1)), round(deltafrac12(2)));
    realdisp = real*deltafrac12;
    realdis = realdisp(0)*realdisp(0) + realdisp(1)*realdisp(1) + realdisp(2)*realdisp(2);
    return sqrt(realdis);
}
unordered_map<int, int> Lattice::removedublicates(vector<double> temp)
{
    std::vector<double>::iterator it;
    unordered_map<int, int> disttoshell;
    int shell = 0;
    vector<double> new_temp;
    double prev, curr;
    int len = temp.size();
    prev = temp[0];
    disttoshell[prev] = 0;
    if(len < 2)
        return disttoshell;
    new_temp.push_back(prev);
    for (it = temp.begin() + 1; it != temp.end(); it++) {
        curr = *it;
        if(fabs(curr - prev) < 0.01)
        {
            disttoshell[int(100*curr)] = shell;
            continue;
        }
        else
        {
            shell += 1;
            prev = curr;
        }
          //  cout << curr << "  -- shell -- " << shell << endl;
    }
    return disttoshell;
}

vector<vector<double>> Lattice::distanceMatrix(){
        int latticeSites = supercell.size();
        double deltaR;
        //string strdistanc;
        std::vector<std::vector<double>> distances (latticeSites, std::vector<double>(latticeSites));
        for(int i = 0; i < latticeSites; i++){
            for(int j = 0; j < latticeSites; j++)
            {
                if(i == j)
                    deltaR = 0.0;
                else
                    deltaR = imageDistance(supercell[i],supercell[j]);
                distances[i][j] = deltaR;
          //      strdistanc += "  " + to_string(deltaR);
            }
            //cout << strdistanc << endl;
        }
        return distances;
}
vector<vector <int>> Lattice::neighbourMatrix(){
    int latticeSites = supercell.size();
    std::vector<double> neighdistances;
    unordered_map<int, int> disttoshell;
    std::vector<std::vector<int>> neighs (latticeSites, std::vector<int>(latticeSites));
    double maxneighdist = mindist();
    maxneigh = INT_MIN;
//    cout << "MAX NEIGH DIST IS - > " << maxneighdist << endl;
    vector<double> localdists;
    int shellnum;
    for(int i = 0; i < latticeSites; i++){
        localdists = distances[i];
        std::sort(localdists.begin(), localdists.end());
        disttoshell = removedublicates(localdists);
        localdists.clear();
        for(int j = 0; j < latticeSites; j++)
        {
            if(i == j || distances[i][j] >= maxneighdist - 0.0001)
                shellnum = 0;
            else
                shellnum = disttoshell[100*distances[i][j]];
            neighs[i][j] = shellnum;
//	    cout << "Atom - " << i << " Atom - " << j << " shellnum - " << shellnum << endl;
            maxneigh = max(shellnum, maxneigh); // Set the number of neighbour shells allowed in a conventional unit cell.
        }
        disttoshell.clear();
    }
    return neighs;
}

double Lattice::mindist()
{
    point3d a1 = celldim[0]*basisset.a;
    point3d b1 = celldim[1]*basisset.b;
    point3d c1 = celldim[2]*basisset.c;
    point3d ab = cross(a1, b1);
    point3d bc = cross(b1, c1);
    point3d ca = cross(a1, c1);
    double modab = sqrt(ab.x*ab.x + ab.y*ab.y + ab.z*ab.z);
    double modbc = sqrt(bc.x*bc.x + bc.y*bc.y + bc.z*bc.z);
    double modca = sqrt(ca.x*ca.x + ca.y*ca.y + ca.z*ca.z);
    double w1 = fabs(dot(cross(a1, b1), c1))/modab;
    double w2 = fabs(dot(cross(b1, c1), a1))/modbc;
    double w3 = fabs(dot(cross(a1, c1), b1))/modca;
    return 0.5*min(w1, min(w2, w3));
}

vector<double> sortneighlist(vector<double> temp, double dis)
{
    int len = temp.size();
    std::vector<double>::iterator it;
    it = temp.begin();
    double var;
    for (int i = 0; i < len; i++) {
        if(fabs(temp[i] - dis) < 0.003)
            continue;
        else
        {
            if(dis < temp[i])
            {
                temp.insert(it + i, dis);
                return temp;
            }
                
        
        }
    }
    temp.push_back(dis);
    return temp;

}
vector<vector<int>> Lattice::numneighofeachshell()
{
    int latticeSites = supercell.size();
    vector<vector<int>> numofneighs; // Number of neighbours in each shell
    vector<int> tempneighs;
    int shell = 0, countperatom;
    //string atomsinshell;
    for (int i = 0; i < latticeSites; i++) {
        shell +=1;
        countperatom = std::count(neighbours[i].begin(), neighbours[i].end(), shell);
        tempneighs.push_back(countperatom);
        while (countperatom > 0) {
            shell +=1;
            countperatom = std::count(neighbours[i].begin(), neighbours[i].end(), shell);
         //   atomsinshell = atomsinshell + " - Shell -  " + to_string(countperatom);
            tempneighs.push_back(countperatom);
        }
        //cout << atomsinshell << endl;
        numofneighs.push_back(tempneighs);
        tempneighs.clear();
        shell = 0;
    }
 //   std::cout << "maximum num of neighbours are  " << numofneighs << endl; 
    return numofneighs;
}
void Lattice::maxneighshell()
{
    int latticeSites = supercell.size();
    vector<vector<vector<int>>> temp_atoms (latticeSites, vector<vector<int>>(maxneigh, vector<int>(0)));
    //vector<int, int> temp_neigh(maxneigh);
    for(int i = 0; i < latticeSites; i++)
    {
        for(int j = 0; j < latticeSites; j++)
        {
            if(neighbours[i][j] > 0)
            {
                temp_atoms[i][neighbours[i][j] - 1].push_back(j);
                
            }
        }
    }
    atoms_shell_atom = temp_atoms;
}

